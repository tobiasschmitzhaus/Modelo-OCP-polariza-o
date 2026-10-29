import unittest
import numpy as np
from transporte import Malha
from escoamento import Escoamento


class Fluid(unittest.TestCase):
    def test_zero_flow_and_boussinesq_guard(self):
        g=Malha(3,3,grading=1.)
        f=Escoamento(g,{0:.05})
        reference=np.zeros(7); c=np.zeros((7,g.n))
        state,vx,vy,di=f.passo(None,c,reference,.01)
        self.assertLess(np.max(np.abs(vx)),1e-12)
        self.assertLess(np.max(np.abs(vy)),1e-12)
        c[0]=1000.
        with self.assertRaises(ValueError): f.passo(state,c,reference,.01)

    def test_manufactured_incompressible_polynomial(self):
        import dolfin as fe
        g=Malha(4,4,L=1.,H=1.,grading=1.)
        f=Escoamento(g)
        dt=.1
        exact=fe.Expression(('x[0]*x[0]','-2*x[0]*x[1]'),degree=2)
        force=fe.Expression(('rho_dt*x[0]*x[0]+1-0.002','-2*rho_dt*x[0]*x[1]+1'),rho_dt=1000/dt,degree=2)
        out,_,_,_=f.passo(None,np.zeros((7,g.n)),np.zeros(7),dt,force,exact)
        u=out.split(deepcopy=True)[0]
        self.assertLess(fe.errornorm(exact,u,'L2',degree_rise=1),1e-9)


if __name__=='__main__': unittest.main()
