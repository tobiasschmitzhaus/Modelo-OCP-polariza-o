"""Modelo eletroneutro de componentes; integração implícita conservativa.

Os fluxos de carga são exatamente F vezes a combinação dos mesmos fluxos
usados nos balanços. Dupla camada: adsorção reversível de NO3- referida ao
estado inicial, q=Cdl(E-Phi). O inventário de superfície entra no balanço N.
"""
from dataclasses import dataclass
from types import SimpleNamespace
import numpy as np
from scipy.optimize import least_squares, brentq
from scipy.sparse import lil_matrix
from termodinamica import MotorTermodinamico, COMPONENTES
from cinetica import MotorCinetico
from transporte import Malha


@dataclass
class Estado:
    x: np.ndarray
    time: float = 0.
    fluid: object = None
    exchange: object = None
    info: object = None


class Modelo:
    nvar = 11
    scales = np.array([1e-3,.273,.014,100.,100.,1e-5,.01])
    solid_scale = 1e-3

    def __init__(self, grid=None, motor=None, bath_at=1e-4, bath_fe=1e-12,
                 bath_dic=.014, bath_no3=100., bath_no2=1e-12, flow=None,
                 reactions=True, closed=False, tolerance=1e-7):
        self.g = grid or Malha()
        self.m = motor or MotorCinetico()
        if abs(self.m.T-298.15) > 1e-9:
            raise ValueError("Constantes de equilíbrio disponíveis apenas em 298.15 K")
        self.thermo = MotorTermodinamico()
        self.flow, self.reactions, self.closed, self.tolerance = flow,reactions,closed,tolerance
        na = bath_no3+bath_at-2*bath_fe
        self.bath = self.thermo.especiar(bath_at,bath_fe,bath_dic,self.m.C_O2_sat,na,bath_no3,bath_no2)
        self.reference = np.array([bath_fe,self.m.C_O2_sat,bath_dic,na,bath_no3,bath_no2,bath_at])
        self.wall = np.arange(self.g.ny)*self.g.nx
        self.pattern = self.sparsity()
        self.initial_delta = None
        self.diagnostics = {}

    def initial(self):
        row = np.zeros(self.nvar)
        row[:6] = np.log(self.reference[:6]/self.scales[:6])
        row[6],row[7] = self.bath.pH, np.log(self.bath.I/.1)
        a = {s:float(self.bath.atividade(s)) for s in self.bath.c}
        if self.reactions:
            e = brentq(lambda E:self.m.correntes(E,a,float(self.bath.c['O2']))[1],-2.,2.)
        else:
            e = 0.
        x = np.r_[np.tile(row,(self.g.n,1)).ravel(),e/self.m.rt_f]
        self.initial_delta = np.full(self.g.ny,e)
        return Estado(x,exchange=np.zeros(9))

    def decode(self, x):
        f = x[:-1].reshape(self.g.n,self.nvar)
        c = np.exp(f[:,:6])*self.scales[:6]
        esp = self.thermo.de_atividade_h(10.**(-f[:,6]),c[:,0],c[:,2],c[:,1],c[:,3],c[:,4],c[:,5],.1*np.exp(f[:,7]))
        components = np.vstack((c.T,esp.AT))
        return esp,components,f[:,8].reshape(self.g.shape)*self.m.rt_f,

    def details(self, x):
        esp,c,phi = self.decode(x)
        f = x[:-1].reshape(self.g.n,self.nvar)
        s2,s3 = f[:,9]*self.solid_scale, f[:,10]*self.solid_scale
        E = x[-1]*self.m.rt_f
        h = (s2[self.wall]*self.m.V_FeOH2+s3[self.wall]*self.m.V_FeOH3)*self.g.dx[0]/(1-self.m.porosidade)
        return esp,c,phi,s2,s3,E,h

    def surface(self, x, old, dt):
        esp,c,phi,s2,s3,E,h = self.details(x)
        _,_,oldphi,_,_,oldE,_ = self.details(old.x)
        delta = E-phi[:,0]
        olddelta = oldE-oldphi[:,0]
        cap = self.m.C_dl*(delta-olddelta)/dt
        farad,cur,src = np.zeros(self.g.ny),np.zeros((self.g.ny,4)),np.zeros((7,self.g.ny))
        if self.reactions:
            for k,cell in enumerate(self.wall):
                a = {s:float(esp.atividade(s)[cell]) for s in esp.c}
                cur[k],farad[k] = self.m.correntes(delta[k],a,float(c[1,cell]),h[k])
                src[:,k] = self.m.fontes_superficie(cur[k])
        # E positive attracts nitrate into the compact layer; no fictitious OH flux.
        src[4] -= cap/self.m.F
        gamma = 2*self.m.C_dl/self.m.F + self.m.C_dl*(delta-self.initial_delta)/self.m.F
        return src,farad,cap,gamma,cur

    def residual(self, x, old, dt, target=None, velocities=None, return_info=False):
        esp,c,phi,s2,s3,E,h = self.details(x)
        _,cn,_,s2n,s3n,_,_ = self.details(old.x)
        vx,vy = (None,None) if velocities is None else velocities
        div,boundary,_ = self.g.fluxos(esp,phi,self.bath,vx,vy,self.m.rt_f,self.closed)
        sources,jf,jc,gamma,cur = self.surface(x,old,dt)
        source = np.zeros_like(c)
        source[:,self.wall] = sources/self.g.dx[0]
        p,ox = self.m.taxas_solidos(esp,s2,s3) if self.reactions else (np.zeros_like(s2),np.zeros_like(s2))
        source[0] -= p+4*ox
        source[1] -= ox
        source[6] -= 2*p+8*ox
        balance = c-cn+dt*(div-source)
        r = np.zeros((self.g.n,self.nvar))
        r[:,:7] = (balance/self.scales[:,None]).T
        r[:,7] = np.log(self.thermo.forca_ionica(esp)/esp.I)
        r[:,8] = self.thermo.carga(esp)/.01
        r[:,9] = (s2-s2n-dt*p)/self.solid_scale
        r[:,10] = (s3-s3n-dt*4*ox)/self.solid_scale
        # Floating electrode, or imposed potential vs SHE (same solution reference).
        last = np.mean(jf+jc) if target is None else (E-target)/self.m.rt_f
        if self.closed:
            # With reactions disabled this gauge replaces a redundant charge row.
            if self.reactions:
                raise ValueError('Célula fechada com eletrodo requer contraeletrodo explícito')
            r[0,8] = phi.ravel()[0]/self.m.rt_f
        result = np.r_[r.ravel(),last]
        if return_info:
            return result,dict(boundary=boundary,surface=sources,faradaic=jf,capacitive=jc,
                               gamma=gamma,currents=cur,balance=balance,components=c,
                               s2=s2,s3=s3,E=E,h=h,esp=esp,phi=phi)
        return result

    def sparsity(self):
        n,nv,g = self.g.n,self.nvar,self.g
        a = lil_matrix((n*nv+1,n*nv+1),dtype=int)
        for cell in range(n):
            y,x = divmod(cell,g.nx)
            near = [cell]
            for yy,xx in ((y-1,x),(y+1,x),(y,x-1),(y,x+1)):
                if 0 <= yy < g.ny and 0 <= xx < g.nx:
                    near.append(yy*g.nx+xx)
            rows = list(range(cell*nv,(cell+1)*nv))
            for neighbor in near:
                a[np.ix_(rows,list(range(neighbor*nv,(neighbor+1)*nv)))] = 1
            if x == 0:
                a[rows,-1] = 1
        for cell in self.wall:
            a[-1,cell*nv:(cell+1)*nv] = 1
        a[-1,-1] = 1
        return a.tocsr()

    def solve(self, old, dt, target=None):
        if dt <= 0 or not np.isfinite(dt):
            raise ValueError('dt deve ser positivo e finito')
        n = self.g.n
        lower = np.tile(np.array([-55]*6+[-2.,-30.,-100.,0.,0.]),(n,1)).ravel()
        upper = np.tile(np.array([25]*6+[16.,np.log(5.),100.,1e8,1e8]),(n,1)).ravel()
        lower,upper = np.r_[lower,-100.],np.r_[upper,100.]
        guess = old.x.copy()
        # Predictor only for the initial iterate: trace Fe can increase by many
        # decades in the first step. Acceptance still uses the full residual.
        if self.reactions:
            _,oldc,_=self.decode(old.x)
            rounds=6 if np.min(oldc[[0,5]]) < 1e-8 else 1
            for _ in range(rounds):
                wallsource=self.surface(guess,old,dt)[0]
                predicted=oldc[:,self.wall]+dt*wallsource/self.g.dx[0]
                # This is an iterate, not a positivity repair of an accepted state.
                # Negative explicit guesses revert to the previous positive value.
                predicted[:6]=np.where(predicted[:6]>0,predicted[:6],oldc[:6,self.wall])
                gf=guess[:-1].reshape(n,self.nvar)
                gf[self.wall,:6]=np.log(predicted[:6].T/self.scales[:6])
                try:
                    spec=self.thermo.especiar(predicted[6],predicted[0],predicted[2],predicted[1],
                                              predicted[3],predicted[4],predicted[5])
                    gf[self.wall,6]=spec.pH
                    gf[self.wall,7]=np.log(spec.I/.1)
                except (ValueError,RuntimeError):
                    pass  # retain the previous valid thermodynamic iterate
        if target is not None:
            guess[-1]=target/self.m.rt_f
        guess = np.maximum(guess,lower+1e-14)
        velocities,fluid,flowinfo = None,old.fluid,dict(Re=0.,density_ratio=0.)
        # Strong Picard coupling; every attempt starts from the SAME accepted fluid.
        previous_velocity = None
        for coupling in range(12):
            if self.flow is not None:
                _,comp,_ = self.decode(guess)
                fluid,vx,vy,flowinfo = self.flow.passo(old.fluid,comp,self.reference,dt)
                velocities = (vx,vy)
            # A small physical residual is the convergence criterion. Waiting
            # for a gradient norm can waste iterations on zero mineral bounds.
            class ResidualConverged(Exception):
                pass
            evaluations=0
            def fun(x):
                nonlocal evaluations
                evaluations += 1
                residual=self.residual(x,old,dt,target,velocities)
                if np.max(np.abs(residual)) < self.tolerance*.1:
                    raise ResidualConverged(x.copy())
                return residual
            linear = ({'tr_solver':'exact'} if n <= 16 else
                      {'jac_sparsity':self.pattern,'tr_solver':'lsmr',
                       'tr_options':{'atol':1e-13,'btol':1e-13,'maxiter':2000}})
            try:
                fit = least_squares(fun,guess,bounds=(lower,upper),
                                    xtol=1e-11,ftol=1e-11,gtol=1e-11,max_nfev=150,x_scale='jac',**linear)
            except ResidualConverged as done:
                fit=SimpleNamespace(x=done.args[0],nfev=evaluations,message='Resíduo físico convergiu')
            r,info = self.residual(fit.x,old,dt,target,velocities,True)
            if np.max(np.abs(r)) > self.tolerance or not np.all(np.isfinite(fit.x)):
                idx=int(np.argmax(np.abs(r)))
                label='controle E' if idx == n*self.nvar else (*COMPONENTES,'I','carga','FeOH2','FeOH3')[idx % self.nvar]
                raise RuntimeError(f'Solve não convergiu: {label}, resíduo={np.max(np.abs(r)):.3e}; {fit.message}')
            guess = fit.x
            if self.flow is None:
                break
            if previous_velocity is not None:
                error = max(np.max(np.abs(v-ov)) for v,ov in zip(velocities,previous_velocity))
                if error < 1e-10:
                    break
            previous_velocity = tuple(v.copy() for v in velocities)
        else:
            raise RuntimeError('Acoplamento química/escoamento não convergiu')
        # The integral of a conservative FV residual is a mass-balance error.
        integrated = np.abs(info['balance'] @ self.g.vol.ravel())
        allowed = self.tolerance*self.scales*np.sum(self.g.vol)
        if np.any(integrated > allowed):
            raise RuntimeError('Balanço integral de componentes excede a tolerância')
        solid_fraction = info['s2']*self.m.V_FeOH2+info['s3']*self.m.V_FeOH3
        if np.max(solid_fraction) > .05:
            raise RuntimeError('Fração volumétrica sólida excede 5%; geometria/porosidade fixa deixa de ser válida')
        if np.min(info['gamma']) < 0:

            raise RuntimeError('Capacidade do modelo de adsorção de nitrato excedida')
        # Exchange entries: seven bath components, Fe from metal, total N from bath.
        exchange = np.zeros(9) if old.exchange is None else old.exchange.copy()
        exchange[:7] -= dt*info['boundary']
        exchange[7] += dt*np.sum(info['currents'][:,0])/(2*self.m.F)*self.g.dy
        exchange[8] -= dt*(info['boundary'][4]+info['boundary'][5])
        result = Estado(fit.x,old.time+dt,fluid,exchange,info)
        self.diagnostics = dict(max_residual=float(np.max(np.abs(r))), nfev=fit.nfev, **flowinfo)
        return result

    def inventory(self, state):
        esp,c,phi,s2,s3,E,h = self.details(state.x)
        v = self.g.vol.ravel()
        delta = E-phi[:,0]
        gamma = 2*self.m.C_dl/self.m.F+self.m.C_dl*(delta-self.initial_delta)/self.m.F
        return np.array([np.dot(c[0]+s2+s3,v),
                         np.dot(c[4]+c[5],v)+np.sum(gamma)*self.g.dy,
                         np.dot(c[2],v),np.dot(c[3],v)])

    def error_norm(self, a, b, atol=1e-7, rtol=2e-3):
        _,ca,pa,sa,ta,Ea,_ = self.details(a.x)
        _,cb,pb,sb,tb,Eb,_ = self.details(b.x)
        # Concentrations use component-specific units; no voltage/concentration mix.
        ec = np.max(np.abs(ca-cb)/(atol*self.scales[:,None]+rtol*np.maximum(np.abs(ca),np.abs(cb))))
        es = max(np.max(np.abs(sa-sb)/(atol*self.solid_scale+rtol*np.maximum(sa,sb))),
                 np.max(np.abs(ta-tb)/(atol*self.solid_scale+rtol*np.maximum(ta,tb))))
        ep = max(np.max(np.abs(pa-pb)),abs(Ea-Eb))/(1e-6+rtol*self.m.rt_f)
        ef = 0.
        if self.flow is not None:
            fe=self.flow.fe
            va=a.fluid.split(deepcopy=True)[0]; vb=b.fluid.split(deepcopy=True)[0]
            ef=fe.errornorm(va,vb,'L2',degree_rise=0)/(1e-10*np.sqrt(self.g.L*self.g.H)+rtol*max(fe.norm(va,'L2'),fe.norm(vb,'L2')))
        return max(ec,es,ep,ef)
