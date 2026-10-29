"""Cinética interfacial e reações balanceadas. Parâmetros empíricos não calibrados.

Inclinações beta são coeficientes efetivos (alpha*n), não alpha isolado.
Uma via HER efetiva substitui a soma duplicada de H+/água do original.
"""
from dataclasses import dataclass
import numpy as np
from scipy.optimize import brentq
from termodinamica import MotorTermodinamico


@dataclass
class MotorCinetico:
    F: float = 96485.33212
    R: float = 8.314462618
    T: float = 298.15
    C_dl: float = .2
    E_eq_Fe_padrao: float = -.44
    E_eq_O2_puro: float = 1.229
    E_eq_H2O_puro: float = 0.
    E_eq_NO3_padrao: float = .01
    i_0_Fe: float = .01
    i_0_O2_ref: float = 5e-4
    i_0_H2O: float = 1e-9
    i_0_NO3: float = 1e-4
    beta_a_Fe: float = .226
    beta_c_Fe: float = .5
    beta_c_O2: float = .604
    beta_c_H2O: float = .5
    beta_c_NO3: float = .5
    P_O2: float = .21
    P_H2: float = 1.
    K_H_O2: float = .0013
    kappa_filme: float = .1
    porosidade: float = .4
    # molar volumes in m³/mol, indicative Fe(OH)2 and Fe(OH)3 solids.
    V_FeOH2: float = 2.65e-5
    V_FeOH3: float = 3.4e-5
    k_oxidacao: float = 1.33e12  # M^-3 s^-1; Fe(II) loss rate coefficient
    k_precip: float = 1.       # s^-1
    k_diss: float = .01       # s^-1

    @property
    def C_O2_sat(self):
        return 1000*self.K_H_O2*self.P_O2

    @property
    def rt_f(self):
        return self.R*self.T/self.F

    def nernst(self, a):
        if min(a['Fe'], a['H'], a['OH'], a['NO3'], a['NO2']) <= 0:
            raise ValueError('Nernst requer atividades positivas; estado fora do domínio')
        if self.P_H2 <= 0 or self.P_O2 <= 0:
            raise ValueError('Pressões dos gases devem ser positivas, em atm relativas a 1 atm')
        return np.array([
            self.E_eq_Fe_padrao + self.rt_f/2*np.log(a['Fe']),
            self.E_eq_O2_puro + self.rt_f*np.log(a['H']) + self.rt_f/4*np.log(a['O2']/self.K_H_O2),
            self.E_eq_H2O_puro + self.rt_f*np.log(a['H']) - self.rt_f/2*np.log(self.P_H2),
            self.E_eq_NO3_padrao - self.rt_f/2*np.log(a['NO2']*a['OH']**2/a['NO3']),
        ])

    @staticmethod
    def _exp(x):
        # Não modificar a lei de Tafel: rejeitar iterados fora de representação.
        if np.any(np.abs(x) > 600):
            raise FloatingPointError('Sobretensão fora do domínio numérico exponencial')
        return np.exp(x)

    def correntes(self, delta, a, co2, h=0.):
        """Uma queda ôhmica comum a todas as reações: j=sum(jk(delta-Rf*j)).

        ORR, HER e nitrato usam aproximação catódica irreversível; o Nernst
        define a referência, mas não se afirma reversibilidade desses ramos.
        """
        eq = self.nernst(a)
        rf = h/self.kappa_filme
        def intrinsic(j):
            eta = (delta-rf*j-eq)/self.rt_f
            # i0 de ferro depende da atividade sem duplicá-la em um só ramo.
            i0 = self.i_0_Fe*(a['Fe']/1e-3)**.5
            fe = i0*(self._exp(self.beta_a_Fe*eta[0])-self._exp(-self.beta_c_Fe*eta[0]))
            oxy = -self.i_0_O2_ref*(co2/self.C_O2_sat)*self._exp(-self.beta_c_O2*eta[1])
            if h > 0:
                lim = 4*self.F*2.1e-9*self.porosidade**1.5*co2/h
                oxy = oxy/(1+abs(oxy)/lim) if lim > 0 else 0.
            # Uma única HER efetiva; atividade/pH já entram no Nernst.
            hyd = -self.i_0_H2O*self._exp(-self.beta_c_H2O*eta[2])
            nit = -self.i_0_NO3*(a['NO3']/.1)*self._exp(-self.beta_c_NO3*eta[3])
            return np.array([fe, oxy, hyd, nit])
        if rf == 0:
            cur = intrinsic(0.)
            return cur, cur.sum()
        # Construir o bracket por queda de tensão, sem avaliar um iterado
        # Rf*j bruto que poderia exceder a faixa exponencial representável.
        raw = intrinsic(0.).sum()
        if abs(raw) < 1e-25:
            return intrinsic(0.), raw
        direction = np.sign(raw)
        edge = direction*.025/rf
        for _ in range(12):
            if (edge-intrinsic(edge).sum())*direction >= 0:
                break
            edge *= 2
        else:
            raise RuntimeError('Não foi possível delimitar a queda ôhmica do filme')
        lo, hi = sorted((0.,edge))
        j = brentq(lambda x: x-intrinsic(x).sum(),lo,hi,xtol=1e-12)
        return intrinsic(j), j

    def fontes_superficie(self, cur):
        """Mol/(m² s), ordem Fe,O2,DIC,Na,NO3,NO2,AT. Positivo entra no líquido."""
        fe, oxy, hyd, nit = cur
        return np.array([fe/(2*self.F), oxy/(4*self.F), 0., 0.,
                         nit/(2*self.F), -nit/(2*self.F), -(oxy+hyd+1.5*nit)/self.F])

    def taxas_solidos(self, esp, s2, s3):
        """Extensões Fe(OH)2(s) e oxidação a Fe(OH)3(s), mol/(m³ s).

        4Fe²+ + O2 + 10H2O -> 4Fe(OH)3(s) + 8H+.
        Não há Fe³+ aquoso nem conversão fictícia goethita/magnetita.
        Dissolução Fe(OH)2 é limitada pelo inventário no solve implícito.
        """
        ratio = esp.atividade('Fe')*esp.atividade('OH')**2/MotorTermodinamico.Ksp
        precip = self.k_precip*esp.c['Fe']*np.maximum(ratio-1,0)/(1+np.maximum(ratio-1,0))
        dissol = self.k_diss*s2*np.maximum(1-ratio,0)
        oxid = .25*self.k_oxidacao*(esp.c['Fe']/1000)*(esp.c['O2']/1000)*(esp.c['OH']/1000)**2*1000
        return precip-dissol, oxid
