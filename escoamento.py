"""Oseen transiente P2/P1: paredes no metal, tração zero nos lados abertos.

A advecção usa velocidade anterior, a tensão é 2 mu sym(grad u), e a
tração aberta fixa a referência de pressão sem âncora pontual adicional.
O solve é COMM_SELF: somente o coordenador MPI executa o modelo.
"""
import numpy as np


class Escoamento:
    def __init__(self, malha, density_coefficients=None):
        import dolfin as fe
        self.fe, self.grid = fe, malha
        self.coeff = density_coefficients or {}
        self.mesh = fe.RectangleMesh(fe.MPI.comm_self, fe.Point(0,0), fe.Point(malha.L,malha.H), malha.nx,malha.ny)
        # Mesmo refinamento geométrico usado nos volumes finitos.
        self.mesh.coordinates()[:,0] = malha.L*(self.mesh.coordinates()[:,0]/malha.L)**malha.grading
        el = fe.MixedElement([fe.VectorElement('P',self.mesh.ufl_cell(),2),fe.FiniteElement('P',self.mesh.ufl_cell(),1)])
        self.W = fe.FunctionSpace(self.mesh,el)
        self.S = fe.FunctionSpace(self.mesh,'DG',0)

    def passo(self, old, components, reference, dt, force=None, boundary_velocity=None):
        fe,g = self.fe,self.grid
        old = fe.Function(self.W) if old is None else old
        rho_delta = np.zeros(g.shape)
        for k,b in self.coeff.items():
            rho_delta += b*(components[k].reshape(g.shape)-reference[k])
        if np.max(np.abs(rho_delta))/1000 > .01:
            raise ValueError('Hipótese Boussinesq inválida: |delta rho|/rho > 1%')
        density = fe.Function(self.S)
        coords = self.S.tabulate_dof_coordinates().reshape(-1,2)
        ix = np.clip(np.searchsorted(g.xf,coords[:,0])-1,0,g.nx-1)
        iy = np.clip(np.searchsorted(g.yf,coords[:,1])-1,0,g.ny-1)
        density.vector().set_local(rho_delta[iy,ix]); density.vector().apply('insert')
        delta_t=fe.Constant(dt)
        u,p = fe.TrialFunctions(self.W); v,q = fe.TestFunctions(self.W)
        un,_ = fe.split(old)
        a = (1000/delta_t*fe.inner(u,v)+1000*fe.inner(fe.dot(fe.grad(u),un),v)
             +2e-3*fe.inner(fe.sym(fe.grad(u)),fe.sym(fe.grad(v)))) *fe.dx
        a += (-p*fe.div(v)+q*fe.div(u))*fe.dx
        rhs = (1000/delta_t*fe.inner(un,v)+fe.inner(fe.as_vector((0.,-9.81*density)),v))*fe.dx
        if force is not None:
            rhs += fe.inner(force,v)*fe.dx
        bc = fe.DirichletBC(self.W.sub(0),fe.Constant((0.,0.)),lambda x,on:on and fe.near(x[0],0.))
        bcs=[bc]
        if boundary_velocity is not None:
            bcs=[fe.DirichletBC(self.W.sub(0),boundary_velocity,'on_boundary'),
                 fe.DirichletBC(self.W.sub(1),fe.Constant(0.),
                                lambda x,on: fe.near(x[0],0.) and fe.near(x[1],0.),'pointwise')]
        out = fe.Function(self.W)
        fe.solve(a == rhs,out,bcs,solver_parameters={'linear_solver':'mumps'})
        velocity = out.split(deepcopy=True)[0]
        vx = np.array([[velocity(float(x),float(y))[0] for x in g.xf] for y in g.y])
        vy = np.array([[velocity(float(x),float(y))[1] for x in g.x] for y in g.yf])
        umax = max(np.abs(vx).max(),np.abs(vy).max())
        return out,vx,vy,dict(Re=1000*umax*min(g.L,g.H)/1e-3,
                              density_ratio=np.max(np.abs(rho_delta))/1000)
