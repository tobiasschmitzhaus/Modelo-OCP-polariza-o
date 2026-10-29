"""Convergência espacial/temporal da evolução de uma espécie neutra.

A comparação temporal usa a solução exata do operador espacial discreto;
a espacial usa a solução contínua de UM passo de Euler implícito. Assim
cada teste isola o erro que pretende medir, sem confundir h e dt.
"""
import unittest
import numpy as np
from modelo import Modelo
from transporte import Malha


def evolution(nx,steps,T=.01):
    g=Malha(nx,2,L=1e-4,H=1e-4,grading=1.)
    m=Modelo(g,reactions=False,closed=True,tolerance=1e-8)
    s=m.initial()
    k=np.pi/g.L
    avg=(np.sin(k*g.xf[1:])-np.sin(k*g.xf[:-1]))/(k*g.dx)
    profile=np.tile(avg,(g.ny,1)).ravel()
    s.x[:-1].reshape(g.n,m.nvar)[:,1]=np.log((m.reference[1]+.05*profile)/m.scales[1])
    initial=m.inventory(s)
    for _ in range(steps): s=m.solve(s,T/steps)
    actual=m.details(s.x)[1][1]
    return m,g,profile,actual,initial,s


class Convergence(unittest.TestCase):
    def test_spatial_order(self):
        errors=[]
        for nx in (4,8):
            m,g,p,c,initial,s=evolution(nx,1)
            exact=m.reference[1]+.05*p/(1+2.1e-9*(np.pi/g.L)**2*.01)
            errors.append(np.sqrt(np.mean((c-exact)**2)))
            np.testing.assert_allclose(m.inventory(s),initial,rtol=1e-9,atol=1e-15)
        self.assertGreater(errors[0]/errors[1],3.5)

    def test_temporal_order(self):
        errors=[]
        for steps in (1,2,4):
            m,g,p,c,_,_=evolution(4,steps)
            lam=4*2.1e-9*np.sin(np.pi/(2*g.nx))**2/g.dx[0]**2
            exact=m.reference[1]+.05*p*np.exp(-lam*.01)
            errors.append(np.sqrt(np.mean((c-exact)**2)))
        self.assertGreater(errors[0]/errors[1],1.8)
        self.assertGreater(errors[1]/errors[2],1.8)


if __name__=='__main__': unittest.main()
