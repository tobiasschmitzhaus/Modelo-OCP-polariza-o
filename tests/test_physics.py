import unittest
import numpy as np
from numpy.testing import assert_allclose
from termodinamica import MotorTermodinamico, ESPECIES
from cinetica import MotorCinetico
from transporte import Malha, fluxo_sg
from modelo import Modelo


class Thermodynamics(unittest.TestCase):
    def setUp(self): self.t=MotorTermodinamico()

    def test_water_activity_and_mass_action(self):
        e=self.t.especiar(0.,0.,0.,I=.1)
        self.assertAlmostEqual(float(e.pH),7.,places=10)
        self.assertAlmostEqual(float(e.atividade('H')*e.atividade('OH')),1e-14,places=25)

    def test_component_mass_and_equilibria(self):
        e=self.t.especiar(np.array([1e-4,1.,-1.]),np.array([1e-6,.1,.01]),.014)
        assert_allclose(e.c['Fe']+e.c['FeOH'],[1e-6,.1,.01],rtol=1e-13)
        assert_allclose(e.c['CO2']+e.c['HCO3']+e.c['CO3'],.014,rtol=1e-13)
        assert_allclose(e.atividade('FeOH')*e.atividade('H')/e.atividade('Fe'),self.t.Kh1,rtol=1e-12)
        assert_allclose(self.t.forca_ionica(e),e.I,atol=2e-12)
        assert_allclose(e.AT,[1e-4,1.,-1.],atol=1e-10)

    def test_invalid_salinity_and_acidity_not_silent(self):
        with self.assertRaises(ValueError): self.t.gammas(.6)
        with self.assertRaises(ValueError): self.t.especiar(-1e12,0.,0.,I=.1)
        e=self.t.especiar(-1000.,0.,0.,I=.1)
        self.assertAlmostEqual(float(e.AT),-1000.,places=7)

    def test_saturation_same_activity_convention(self):
        e=self.t.especiar(1.,.01,.014)
        si=self.t.calcular_SI_FeOH2(e.c['Fe'],e.c['H'],e.c['OH'],e.I*1000)
        assert_allclose(si,self.t.indice_saturacao(e),atol=1e-13)


class Kinetics(unittest.TestCase):
    def setUp(self):
        self.m=MotorCinetico()
        self.e=MotorTermodinamico().especiar(1e-4,1e-6,.014,self.m.C_O2_sat,100.,100.,1e-6)
        self.a={s:float(self.e.atividade(s)) for s in self.e.c}

    def test_iron_zero_at_nernst(self):
        eq=self.m.nernst(self.a)[0]
        cur,_=self.m.correntes(eq,self.a,float(self.e.c['O2']))
        self.assertAlmostEqual(cur[0],0.,places=15)

    def test_film_common_drop(self):
        cur,j=self.m.correntes(-.2,self.a,float(self.e.c['O2']),1e-6)
        raw,_=self.m.correntes(-.2-1e-6/self.m.kappa_filme*j,self.a,float(self.e.c['O2']),0.)
        # The ORR additionally has an oxygen transport resistance in the film.
        self.assertAlmostEqual(cur.sum(),j,places=7)
        assert_allclose(cur[[0,2,3]],raw[[0,2,3]],rtol=1e-8)

    def test_nitrogen_and_charge_stoichiometry(self):
        for cur in (np.array([2.,-3.,-1.,-4.]),np.array([-1.,0.,0.,0.])):
            s=self.m.fontes_superficie(cur)
            self.assertAlmostEqual(s[4]+s[5],0.)
            # q/F = 2 Fe + Na - NO3 - AT (AT includes nitrite base).
            charge=2*s[0]+s[3]-s[4]-s[6]
            self.assertAlmostEqual(charge*self.m.F,cur.sum(),places=12)

    def test_oxidation_and_precipitation_conserve(self):
        p,ox=self.m.taxas_solidos(self.e,.001,.002)
        # Fe total including BOTH minerals and oxygen stoichiometry.
        delta=np.array([-p-4*ox,-ox,-2*p-8*ox,p,4*ox])
        self.assertAlmostEqual(delta[0]+delta[3]+delta[4],0.)
        self.assertAlmostEqual(2*delta[0]-delta[2],0.)


class Transport(unittest.TestCase):
    def test_boltzmann_zero_flux(self):
        rt=.0257; z=2; dv=.1; gl=.6; gr=.8; cl=1.
        cr=cl*gl/gr*np.exp(-z*dv/rt)
        j=fluxo_sg(cl,cr,gl,gr,z,1e-9,dv,1e-5,0.,rt)
        self.assertLess(abs(float(j)),1e-18)

    def test_linear_diffusion_flux(self):
        j=fluxo_sg(1.,3.,1.,1.,0,2e-9,0.,.01,0.,.0257)
        self.assertAlmostEqual(float(j),-4e-7,places=18)

    def test_internal_flux_cancellation(self):
        g=Malha(4,3,grading=2)
        t=MotorTermodinamico()
        e=t.especiar(np.linspace(.01,.02,g.n),np.linspace(.001,.002,g.n),.014)
        bath=t.especiar(.01,.001,.014)
        div,boundary,_=g.fluxos(e,np.zeros(g.shape),bath)
        assert_allclose(div@g.vol.ravel(),boundary,atol=1e-20)

    def test_uniform_closed_equilibrium(self):
        model=Modelo(Malha(2,2),reactions=False,closed=True)
        state=model.initial()
        next_state=model.solve(state,.1)
        assert_allclose(model.inventory(next_state),model.inventory(state),rtol=1e-10,atol=1e-15)
        self.assertLess(np.max(np.abs(model.thermo.carga(model.details(next_state.x)[0]))),1e-7)


class Capacitor(unittest.TestCase):
    def test_linear_rc_electrolyte_network(self):
        g=Malha(2,2,grading=1.)
        m=Modelo(g,reactions=False,tolerance=1e-10)
        old=m.initial(); dt=.01; imposed=1e-4
        new=m.solve(old,dt,target=imposed)
        esp=m.bath
        kappa=m.m.F/m.m.rt_f*sum(z*z*D*float(esp.c[s]) for s,(z,D) in ESPECIES.items())
        G=np.zeros((g.n,g.n))
        for cell in range(g.n):
            y,x=divmod(cell,g.nx)
            for yy,xx,area,dist in ((y,x-1,g.dy,g.dx[0]),(y,x+1,g.dy,g.dx[0]),
                                    (y-1,x,g.dx[x],g.dy),(y+1,x,g.dx[x],g.dy)):
                if 0<=xx<g.nx and 0<=yy<g.ny:
                    j=yy*g.nx+xx; conduct=kappa*area/dist
                    G[cell,cell]+=conduct; G[cell,j]-=conduct
                elif xx != -1:  # reservoir faces, not the metal wall
                    G[cell,cell]+=kappa*area/(dist/2)
        rhs=np.zeros(g.n)
        for cell in m.wall:
            conduct=m.m.C_dl*g.dy/dt
            G[cell,cell]+=conduct; rhs[cell]=conduct*imposed
        exact=np.linalg.solve(G,rhs)
        phi=m.details(new.x)[2].ravel()
        np.testing.assert_allclose(phi,exact,rtol=2e-3,atol=1e-8)
        _,_,cap,_,_=m.surface(new.x,old,dt)
        expected=m.m.C_dl*(imposed-exact[m.wall])/dt
        np.testing.assert_allclose(cap,expected,rtol=2e-3,atol=1e-8)


if __name__=='__main__': unittest.main()
