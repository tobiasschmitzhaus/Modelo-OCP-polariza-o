"""Volumes finitos: fluxo Scharfetter–Gummel conservativo por espécie.

A neutralidade é uma restrição do solve acoplado. Phi é determinada pelas
continuidades de componentes e carga; não se omite a corrente difusiva.
"""
from dataclasses import dataclass
import numpy as np
from termodinamica import ESPECIES, PESOS, COMPONENTES


def bernoulli(x):
    x = np.asarray(x, float)
    out = np.empty_like(x)
    small = np.abs(x) < 1e-4
    out[small] = 1-x[small]/2+x[small]**2/12-x[small]**4/720
    pos, neg = (~small)&(x > 50), (~small)&(x < -50)
    out[pos] = x[pos]*np.exp(-x[pos])
    out[neg] = -x[neg]
    mid = ~(small|pos|neg)
    out[mid] = x[mid]/np.expm1(x[mid])
    return out


def fluxo_sg(cl, cr, gl, gr, z, D, dphi, distance, velocity, rt_f):
    drift = z*dphi/rt_f + np.log(gr/gl)
    return D/distance*(bernoulli(drift)*cl-bernoulli(-drift)*cr) + velocity*np.where(velocity >= 0, cl, cr)


@dataclass
class Malha:
    nx: int = 16
    ny: int = 16
    L: float = .01
    H: float = .02
    grading: float = 3.

    def __post_init__(self):
        if min(self.nx,self.ny) < 2 or min(self.L,self.H,self.grading) <= 0:
            raise ValueError('Malha inválida')
        self.xf = self.L*np.linspace(0,1,self.nx+1)**self.grading
        self.yf = np.linspace(0,self.H,self.ny+1)
        self.x = (self.xf[1:]+self.xf[:-1])/2
        self.y = (self.yf[1:]+self.yf[:-1])/2
        self.dx, self.dy = np.diff(self.xf), self.H/self.ny
        self.vol = np.broadcast_to(self.dx*self.dy,(self.ny,self.nx)).copy()
        self.shape = (self.ny,self.nx)
        self.n = self.nx*self.ny

    def fluxos(self, esp, phi, bath, vx=None, vy=None, rt_f=.02569258, fechado=False):
        # Face flux positive right/up. The left metal face is handled separately.
        vx = np.zeros((self.ny,self.nx+1)) if vx is None else vx
        vy = np.zeros((self.ny+1,self.nx)) if vy is None else vy
        components = np.zeros((7,*self.shape))
        boundary = np.zeros(7)
        species_flux = {}
        for s,(z,D) in ESPECIES.items():
            c,g = esp.c[s].reshape(self.shape), np.broadcast_to(esp.gamma[s],(self.n,)).reshape(self.shape)
            bc,bg = float(bath.c[s]), float(bath.gamma[s])
            fx,fy = np.zeros_like(vx),np.zeros_like(vy)
            fx[:,1:-1] = fluxo_sg(c[:,:-1],c[:,1:],g[:,:-1],g[:,1:],z,D,
                                      phi[:,1:]-phi[:,:-1],np.diff(self.x)[None,:],vx[:,1:-1],rt_f)
            fy[1:-1,:] = fluxo_sg(c[:-1,:],c[1:,:],g[:-1,:],g[1:,:],z,D,
                                      phi[1:,:]-phi[:-1,:],self.dy,vy[1:-1,:],rt_f)
            if not fechado:
                fx[:,-1] = fluxo_sg(c[:,-1],bc,g[:,-1],bg,z,D,-phi[:,-1],self.xf[-1]-self.x[-1],vx[:,-1],rt_f)
                fy[0,:] = fluxo_sg(bc,c[0,:],bg,g[0,:],z,D,phi[0,:],self.dy/2,vy[0,:],rt_f)
                fy[-1,:] = fluxo_sg(c[-1,:],bc,g[-1,:],bg,z,D,-phi[-1,:],self.dy/2,vy[-1,:],rt_f)
            div = (fx[:,1:]-fx[:,:-1])/self.dx[None,:] + (fy[1:,:]-fy[:-1,:])/self.dy
            net = np.sum(fx[:,-1])*self.dy + np.sum((fy[-1,:]-fy[0,:])*self.dx)
            species_flux[s] = (fx,fy)
            for k,comp in enumerate(COMPONENTES):
                weight = PESOS[comp].get(s,0)
                components[k] += weight*div
                boundary[k] += weight*net
        return components.reshape(7,-1), boundary, species_flux
