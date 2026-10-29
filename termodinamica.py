"""Especiação aquosa a 298.15 K; concentrações externas em mol/m³.

Davies é limitado explicitamente a I <= 0.5 mol/L. Atividades são
adimensionais, referidas a 1 mol/L. Fe designa Fe(II) total dissolvido.
"""
from dataclasses import dataclass
import numpy as np


@dataclass
class Especiacao:
    c: dict
    gamma: dict
    I: np.ndarray
    AT: np.ndarray
    pH: np.ndarray

    def atividade(self, nome):
        return self.gamma[nome] * self.c[nome] / 1000.0


# espécie: (carga, D em m²/s)
ESPECIES = {
    'H': (1, 9.31e-9), 'OH': (-1, 5.27e-9),
    'CO2': (0, 1.92e-9), 'HCO3': (-1, 1.18e-9), 'CO3': (-2, .92e-9),
    'Fe': (2, .71e-9), 'FeOH': (1, .71e-9),
    'O2': (0, 2.10e-9), 'Na': (1, 1.33e-9),
    'NO3': (-1, 1.90e-9), 'NO2': (-1, 1.90e-9), 'HNO2': (0, 1.90e-9),
}
# Projeção estequiométrica das espécies para os componentes conservados.
COMPONENTES = ('Fe', 'O2', 'DIC', 'Na', 'NO3', 'NO2', 'AT')
PESOS = {
    'Fe': {'Fe': 1, 'FeOH': 1}, 'O2': {'O2': 1},
    'DIC': {'CO2': 1, 'HCO3': 1, 'CO3': 1},
    'Na': {'Na': 1}, 'NO3': {'NO3': 1}, 'NO2': {'NO2': 1, 'HNO2': 1},
    'AT': {'OH': 1, 'H': -1, 'HCO3': 1, 'CO3': 2, 'FeOH': 1, 'NO2': 1},
}


class MotorTermodinamico:
    Kw = 1e-14
    Ka1 = 10**-6.35
    Ka2 = 10**-10.33
    Kh1 = 10**-9.5
    Ka_NO2 = 10**-3.3
    Ksp = 8e-16
    I_max = .5

    def gammas(self, I):
        I = np.asarray(I, dtype=float)
        if np.any(~np.isfinite(I)) or np.any(I < 0) or np.any(I > self.I_max):
            raise ValueError('Davies requer 0 <= I <= 0.5 mol/L; use outro modelo fora desta faixa')
        b = -.509 * (np.sqrt(I)/(1+np.sqrt(I)) - .3*I)
        return {s: 10**(z*z*b) for s, (z, _) in ESPECIES.items()}

    def de_atividade_h(self, aH, Fe, DIC, O2, Na, NO3, NO2, I):
        """Fe total = Fe²+ + FeOH+; carbono = CO2* + HCO3- + CO3²-."""
        vals = np.broadcast_arrays(*[np.asarray(x, float) for x in (aH, Fe, DIC, O2, Na, NO3, NO2, I)])
        aH, Fe, DIC, O2, Na, NO3, NO2, I = vals
        if any(np.any(~np.isfinite(v)) for v in vals) or np.any(aH <= 0):
            raise ValueError('Estado não finito ou atividade H não positiva')
        if any(np.any(v < 0) for v in (Fe, DIC, O2, Na, NO3, NO2)):
            raise ValueError('Concentração negativa')
        g = self.gammas(I)
        r1 = self.Ka1 / (aH*g['HCO3'])
        r2 = self.Ka2*g['HCO3']/(aH*g['CO3'])
        co2 = DIC/(1+r1+r1*r2)
        rf = self.Kh1*g['Fe']/(aH*g['FeOH'])
        free = Fe/(1+rf)
        nitrite_free = NO2/(1+aH*g['NO2']/self.Ka_NO2)
        c = dict(H=1000*aH/g['H'], OH=1000*self.Kw/(aH*g['OH']),
                 CO2=co2, HCO3=co2*r1, CO3=co2*r1*r2,
                 Fe=free, FeOH=free*rf, O2=O2, Na=Na, NO3=NO3, NO2=nitrite_free, HNO2=NO2-nitrite_free)
        at = sum(w*c[s] for s, w in PESOS['AT'].items())
        return Especiacao(c, g, I, at, -np.log10(aH))

    @staticmethod
    def forca_ionica(esp):
        return .5*sum(z*z*esp.c[s] for s, (z, _) in ESPECIES.items())/1000

    @staticmethod
    def carga(esp):
        return sum(z*esp.c[s] for s, (z, _) in ESPECIES.items())

    def especiar(self, AT, Fe, DIC, O2=0., Na=100., NO3=100., NO2=0., I=None):
        """Bisseção com bracket e teste do balanço; I é autoconsistente por padrão."""
        vals = np.broadcast_arrays(*[np.asarray(x, float) for x in (AT, Fe, DIC, O2, Na, NO3, NO2)])
        AT, Fe, DIC, O2, Na, NO3, NO2 = vals
        if I is None:
            ionic = np.full_like(AT, .1)
        else:
            ionic = np.broadcast_to(np.asarray(I, float), AT.shape).copy()
        for _ in range(100):
            lo, hi = np.full_like(AT, -4.), np.full_like(AT, 20.)
            def evaluate(ph):
                return self.de_atividade_h(10.**(-ph), Fe, DIC, O2, Na, NO3, NO2, ionic)
            if np.any(evaluate(lo).AT > AT) or np.any(evaluate(hi).AT < AT):
                raise ValueError('Não há raiz de alcalinidade no intervalo de atividade H suportado')
            for __ in range(82):
                mid = (lo+hi)/2
                below = evaluate(mid).AT < AT
                lo, hi = np.where(below, mid, lo), np.where(below, hi, mid)
            esp = evaluate((lo+hi)/2)
            if np.any(np.abs(esp.AT-AT) > 1e-9*(1+np.abs(AT))):
                raise ValueError('Resíduo de alcalinidade excede a tolerância')
            if I is not None:
                return esp
            new = self.forca_ionica(esp)
            self.gammas(new)  # não limitar artificialmente a força iônica
            if np.max(np.abs(new-ionic)) < 1e-12:
                return esp
            ionic = (ionic+new)/2
        raise RuntimeError('Força iônica/especiação não convergiu')

    def calcular_pH_alcalinidade(self, AT_m3, DIC_m3, I_m3, Fe2_total_m3):
        # Compatibilidade de API: esta entrada especifica I em mol/m³.
        e = self.especiar(AT_m3, Fe2_total_m3, DIC_m3, I=np.asarray(I_m3)/1000)
        return e.c['H'], e.c['OH']

    def indice_saturacao(self, esp):
        q = esp.atividade('Fe')*esp.atividade('OH')**2
        with np.errstate(divide='ignore'):
            return np.log10(q/self.Ksp)

    def calcular_SI_FeOH2(self, Fe2_m3, H_m3, OH_m3, I_m3=100.):
        """API de compatibilidade: Fe2_m3 deve ser Fe²+ LIVRE, não total."""
        g = self.gammas(np.asarray(I_m3)/1000)
        q = (g['Fe']*np.asarray(Fe2_m3)/1000)*(g['OH']*np.asarray(OH_m3)/1000)**2
        with np.errstate(divide='ignore'):
            return np.log10(q/self.Ksp)
