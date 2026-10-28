import fenics as fe
import numpy as np
import ufl
import csv
import time
import math
import os
import matplotlib.pyplot as plt
import time
fe.set_log_level(fe.LogLevel.WARNING)
from datetime import datetime
from mpi4py import MPI 


tempo_inicio_real = time.time()
comm = fe.MPI.comm_world
rank = comm.rank  

if rank == 0:
    print("=== SIMULADOR 2D FASE 8.2: FÍSICA PURA (SEM FREIOS) E SOLVER BLINDADO ===")
    os.makedirs("resultados_2d/perfis", exist_ok=True)
    os.makedirs("resultados_2d/mapas_2d", exist_ok=True)

# 1. MALHA E GEOMETRIA
L, H = 0.01, 0.02   
nx, ny = 40, 80 
malha = fe.RectangleMesh(fe.Point(0, 0), fe.Point(L, H), nx, ny)

class Contorno_Aberto(fe.SubDomain):
    def inside(self, x, on_boundary): return on_boundary and not fe.near(x[0], 0.0)

def Parede_Metal(x, on_boundary):
    return on_boundary and fe.near(x[0], 0.0)

fronteiras = fe.MeshFunction("size_t", malha, malha.topology().dim() - 1)
fronteiras.set_all(0)
Parede_Esquerda = fe.AutoSubDomain(Parede_Metal)
Parede_Esquerda.mark(fronteiras, 1)

dx = fe.Measure('dx', domain=malha, metadata={'quadrature_degree': 4})
ds = fe.Measure('ds', domain=malha, subdomain_data=fronteiras, metadata={'quadrature_degree': 4})

# 2. ESPAÇOS DE FUNÇÃO (FLUIDO + QUÍMICA)
V_v = fe.VectorElement('P', malha.ufl_cell(), 2); V_p = fe.FiniteElement('P', malha.ufl_cell(), 1)
W_fluido = fe.FunctionSpace(malha, fe.MixedElement([V_v, V_p]))
w = fe.Function(W_fluido); w_n = fe.Function(W_fluido)
v_u, v_p = fe.TestFunctions(W_fluido); u_f, p_f = fe.split(w); u_f_n, p_f_n = fe.split(w_n)

from cinetica import MotorCinetico
from termodinamica import MotorTermodinamico
motor = MotorCinetico(); motor_termo = MotorTermodinamico()

P1 = fe.FiniteElement('P', malha.ufl_cell(), 1); R0 = fe.FiniteElement('R', malha.ufl_cell(), 0)
ME = fe.FunctionSpace(malha, fe.MixedElement([P1, P1, P1, P1, P1, P1, P1, P1, R0]))
V_sub = ME.sub(0).collapse() 

u = fe.Function(ME); u_n = fe.Function(ME)
v_Fe2, v_O2, v_AT, v_DIC, v_Na, v_NO3, v_Fe3, v_Phi, v_E = fe.TestFunctions(ME)
C_Fe2, C_O2, C_AT, C_DIC, C_Na, C_NO3, C_Fe3, Phi, E_metal = fe.split(u)
C_Fe2_n, C_O2_n, C_AT_n, C_DIC_n, C_Na_n, C_NO3_n, C_Fe3_n, Phi_n, E_metal_n = fe.split(u_n)

# FUNÇÕES DE APOIO SNIA E NERNST DINÂMICO
H_eq_ini, OH_eq_ini = motor_termo.calcular_pH_alcalinidade(np.array([1e-4]), np.array([1.4e-2]), 100.0, np.array([1e-12]))
H_fn  = fe.Function(V_sub); H_fn.vector()[:]  = H_eq_ini[0]
OH_fn = fe.Function(V_sub); OH_fn.vector()[:] = OH_eq_ini[0]

E_eq_Fe_fn  = fe.Function(V_sub); E_eq_Fe_fn.vector()[:] = motor.E_eq_Fe_padrao
E_eq_O2_fn  = fe.Function(V_sub); E_eq_O2_fn.vector()[:] = 0.815 
E_eq_H2O_fn = fe.Function(V_sub); E_eq_H2O_fn.vector()[:] = -0.413

# Espessura da camada porosa (metros)
h_f_fn = fe.Function(V_sub); h_f_fn.vector()[:] = 0.0

dof_coords = V_sub.tabulate_dof_coordinates()
dofs_fronteira = np.where(np.isclose(dof_coords[:, 0], 0.0))[0]
dist = (dof_coords[:, 0] - 0.0)**2 + (dof_coords[:, 1] - 0.0)**2
dof_origem = np.argmin(dist) if (len(dist) > 0 and np.min(dist) < 1e-12) else -1

# 3. CONDIÇÕES DE CONTORNO
bc_fluid = fe.DirichletBC(W_fluido.sub(0), fe.Constant((0.0, 0.0)), Parede_Metal)
def ponto_ancora(x, on_boundary): return fe.near(x[0], L) and fe.near(x[1], H)
bc_pressao = fe.DirichletBC(W_fluido.sub(1), fe.Constant(0.0), ponto_ancora, 'pointwise')

bcs_quimica = [
    fe.DirichletBC(ME.sub(7), fe.Constant(0.0),    Contorno_Aberto()),
    fe.DirichletBC(ME.sub(1), fe.Constant(motor.C_O2_sat), Contorno_Aberto()),
    fe.DirichletBC(ME.sub(2), fe.Constant(1e-4),   Contorno_Aberto()),
    fe.DirichletBC(ME.sub(3), fe.Constant(1.4e-2), Contorno_Aberto()),
    fe.DirichletBC(ME.sub(0), fe.Constant(1e-16),  Contorno_Aberto()),
    fe.DirichletBC(ME.sub(6), fe.Constant(1e-16),  Contorno_Aberto()),
    fe.DirichletBC(ME.sub(4), fe.Constant(100.0),  Contorno_Aberto()), 
    fe.DirichletBC(ME.sub(5), fe.Constant(100.0),  Contorno_Aberto())  
]

# 4. FORMULAÇÃO DO FLUIDO (BOUSSINESQ)
dt_val = 0.001
dt = fe.Constant(dt_val) 
M_Fe = 0.0558
F_boyancy = fe.as_vector([0.0, -9.81 * (M_Fe*C_Fe2_n + M_Fe*C_Fe3_n)])
F_fluid_eq = (1000.0/dt)*fe.inner(u_f - u_f_n, v_u)*dx - fe.inner(p_f, fe.div(v_u))*dx + 1e-3*fe.inner(fe.grad(u_f), fe.grad(v_u))*dx - fe.inner(F_boyancy, v_u)*dx + fe.inner(fe.div(u_f), v_p)*dx

problema_fluido = fe.NonlinearVariationalProblem(F_fluid_eq, w, bcs=[bc_fluid, bc_pressao], J=fe.derivative(F_fluid_eq, w))
solver_fluido = fe.NonlinearVariationalSolver(problema_fluido)
# Blindagem do solver de fluidos
solver_fluido.parameters['newton_solver']['linear_solver'] = 'mumps'
solver_fluido.parameters['newton_solver']['error_on_nonconvergence'] = False
solver_fluido.parameters['newton_solver']['maximum_iterations'] = 30

# 5. FORMULAÇÃO QUÍMICA E TRANSPORTE
F_const = motor.F; RT = motor.R * motor.T
D_Fe2, z_Fe2 = 0.71e-9, 2.0; D_O2, z_O2 = 2.10e-9, 0.0
D_AT, z_AT   = 1.18e-9, -1.0; D_DIC, z_DIC = 1.92e-9, 0.0 
D_Na, z_Na = motor.D_Na, motor.z_Na; D_NO3, z_NO3 = motor.D_NO3, motor.z_NO3
D_Fe3, z_Fe3 = 0.60e-9, 3.0; D_H, z_H = 9.31e-9, 1.0; D_OH, z_OH = 5.27e-9, -1.0

def c_safe(C): return fe.conditional(fe.lt(C, 1e-12), 1e-12, C)
a_OH = c_safe(OH_fn)

Ka1_ufl = 10**-6.35
Ka2_ufl = 10**-10.33
H_safe = c_safe(H_fn)
H_molar = H_safe / 1000.0  
denom_c = H_molar**2 + Ka1_ufl * H_molar + Ka1_ufl * Ka2_ufl
C_HCO3_ufl = c_safe(C_DIC) * (Ka1_ufl * H_molar) / denom_c

soma_concs = C_HCO3_ufl + a_OH + H_safe
soma_concs_safe = fe.conditional(fe.lt(soma_concs, 1e-16), 1e-16, soma_concs)
D_AT_efetivo = (D_AT * C_HCO3_ufl + D_OH * a_OH + D_H * H_safe) / soma_concs_safe

kappa_local = (F_const**2 / RT) * ( (z_Fe2**2 * D_Fe2 * c_safe(C_Fe2)) + (z_Fe3**2 * D_Fe3 * c_safe(C_Fe3)) + (z_Na**2 * D_Na * c_safe(C_Na)) + (z_NO3**2 * D_NO3 * c_safe(C_NO3)) + (z_H**2 * D_H * c_safe(H_fn)) + (z_OH**2 * D_OH * a_OH) )

J_Fe2 = -D_Fe2 * fe.grad(C_Fe2) - z_Fe2 * (D_Fe2/RT) * F_const * c_safe(C_Fe2) * fe.grad(Phi) + u_f * C_Fe2
J_O2  = -D_O2  * fe.grad(C_O2)  + u_f * C_O2
J_AT  = -D_AT_efetivo * fe.grad(C_AT) - z_AT * (D_AT_efetivo/RT) * F_const * c_safe(C_AT) * fe.grad(Phi) + u_f * C_AT 
J_NO3 = -D_NO3 * fe.grad(C_NO3) - z_NO3 * (D_NO3/RT) * F_const * c_safe(C_NO3) * fe.grad(Phi) + u_f * C_NO3
J_DIC = -D_DIC * fe.grad(C_DIC) + u_f * C_DIC
J_Na  = -D_Na  * fe.grad(C_Na)  - z_Na  * (D_Na/RT)  * F_const * c_safe(C_Na)  * fe.grad(Phi) + u_f * C_Na
J_Fe3 = -D_Fe3 * fe.grad(C_Fe3) - z_Fe3 * (D_Fe3/RT) * F_const * c_safe(C_Fe3) * fe.grad(Phi) + u_f * C_Fe3

# ==============================================================
# FÍSICA PURA: SOBRETENSÕES REAIS (COM BLINDAGEM NUMÉRICA)
# ==============================================================
h_f_safe = c_safe(h_f_fn)
porosidade = 0.4
kappa_filme = 0.1 
R_f = h_f_safe / kappa_filme 

# BLINDAGEM CONTRA OVERFLOW E TRANSIENTES ABSURDOS:
# Limita a força motriz a +/- 0.6 V. 
# Isso permite correntes extremas (milhares de A/m²), 
# mas impede os trilhões de A/m² que destroem a matriz de massa.
def clamp_eta(eta_raw):
    # Alargado para 2.0V para não achatar a curva de Tafel real
    x = eta_raw / 2.0
    exp_2x = fe.exp(2.0 * x)
    tanh_x = (exp_2x - 1.0) / (exp_2x + 1.0)
    return 2.0 * tanh_x

eta_Fe = clamp_eta(E_metal - Phi - E_eq_Fe_fn)
eta_O2 = clamp_eta(E_metal - Phi - E_eq_O2_fn)
eta_H2O = clamp_eta(E_metal - Phi - E_eq_H2O_fn)
# Potencial de Nernst dinâmico para o Nitrato (depende do pH local)
E_eq_NO3_dinamico = motor.E_eq_NO3_padrao + (RT / (2.0 * F_const)) * fe.ln(1e-14 / (a_OH**2 + 1e-16))
eta_NO3 = clamp_eta(E_metal - Phi - E_eq_NO3_dinamico)

a_OH_lim = fe.conditional(fe.gt(a_OH, 1e-2), 1e-2, a_OH)
fator_bockris = (a_OH_lim / motor.a_OH_ref)**motor.ordem_OH

# Correntes Brutas Ativas (Exp pura da cinética de Butler-Volmer)
# A redução do ferro agora é limitada pela disponibilidade de Fe2+ local
C_Fe2_ref = 1.0 # mol/m3 (referência)
i_Fe_nu = motor.i_0_Fe * fator_bockris * (
    fe.exp((motor.beta_a_Fe * F_const * eta_Fe)/RT) - 
    (c_safe(C_Fe2)/C_Fe2_ref) * fe.exp(-(motor.beta_c_Fe * F_const * eta_Fe)/RT)
)
i_O2_nu = -motor.i_0_O2_ref * (c_safe(C_O2) / motor.C_O2_sat) * fe.exp(-(motor.beta_c_O2 * F_const * eta_O2)/RT)
i_H2O_nu = -motor.i_0_H2O * fe.exp(-(motor.beta_c_H2O * F_const * eta_H2O)/RT)
i_H_nu = -motor.i_0_H * fe.exp(-(motor.beta_c_H2O * F_const * eta_H2O)/RT) 

# Reações Adicionais Mecanicistas
i_NO3_nu = -motor.i_0_NO3 * (c_safe(C_NO3) / 100.0) * fe.exp(-(motor.beta_c_NO3 * F_const * eta_NO3)/RT)

# RESTAURADO O CONSUMO DE FERRO PELA CROSTA (FÍSICA REAL)
i_crust_nu = -motor.k_crust * (h_f_safe / 1e-6) * c_safe(C_Fe2) * fe.exp(-(motor.beta_c_crust * F_const * eta_Fe)/RT)

# Correntes Efetivas Modificadas (IR Drop Físico + Difusão)
b_a_Fe = (motor.beta_a_Fe * F_const) / RT
i_Fe_eff = i_Fe_nu / (1.0 + b_a_Fe * R_f * fe.sqrt(i_Fe_nu**2 + 1e-16))

D_O2_eff = D_O2 * (porosidade**1.5) 
i_lim_O2 = (4.0 * F_const * D_O2_eff * c_safe(C_O2)) / (h_f_safe + 1e-12)
b_c_O2 = (motor.beta_c_O2 * F_const) / RT
abs_i_O2_nu = fe.sqrt(i_O2_nu**2 + 1e-16)
i_O2_eff = i_O2_nu / (1.0 + abs_i_O2_nu / i_lim_O2 + b_c_O2 * R_f * abs_i_O2_nu)

b_c_H2O = (motor.beta_c_H2O * F_const) / RT
i_H2O_eff = i_H2O_nu / (1.0 + b_c_H2O * R_f * fe.sqrt(i_H2O_nu**2 + 1e-16))
i_H_eff = i_H_nu / (1.0 + b_c_H2O * R_f * fe.sqrt(i_H_nu**2 + 1e-16))

b_c_NO3 = (motor.beta_c_NO3 * F_const) / RT
i_NO3_eff = i_NO3_nu / (1.0 + b_c_NO3 * R_f * fe.sqrt(i_NO3_nu**2 + 1e-16))

b_c_crust = (motor.beta_c_crust * F_const) / RT
i_crust_eff = i_crust_nu / (1.0 + b_c_crust * R_f * fe.sqrt(i_crust_nu**2 + 1e-16))

# Balanço de Massas (Crosta volta a consumir Fe2+ do meio, cumprindo a lei de conservação)
Fluxo_Fe2_UFL = (i_Fe_eff + i_crust_eff) / (2.0*F_const)
Fluxo_NO3_UFL = i_NO3_eff / (2.0*F_const)
Fluxo_AT_UFL  = -(i_O2_eff + i_H2O_eff + i_H_eff + i_NO3_eff) / F_const 

i_total_faradaic = i_Fe_eff + i_crust_eff + i_O2_eff + i_H2O_eff + i_H_eff + i_NO3_eff

# STUMM-LEE
OH_safe = c_safe(OH_fn)
k_SL_molar = 1.33e12
f_amortecedor = fe.conditional(fe.gt(OH_safe, 0.1), fe.exp(-5.0*(OH_safe-0.1)), 1.0)
r_SL = k_SL_molar * (c_safe(C_Fe2)/1000) * (c_safe(C_O2)/1000) * (OH_safe/1000)**2 * 1000 * f_amortecedor
taxa_maxima_permitida = c_safe(C_AT) / (8.0 * dt)
r_SL_efetivo = fe.conditional(fe.lt(taxa_maxima_permitida, r_SL), taxa_maxima_permitida, r_SL)

# Saturação e Precipitação Homogênea no Bulk (Física do Fe(OH)2)
Ksp_UFL = 8.0e-16
# Produto iônico Q = [Fe2+] * [OH-]^2 (convertido para mol/L)
Q_UFL = (c_safe(C_Fe2)/1000.0) * (a_OH/1000.0)**2
SI_UFL = fe.ln(Q_UFL / Ksp_UFL + 1e-16)

# Taxa de precipitação (mol/(m³ s)) acionada apenas quando supersaturado (SI > 0)
k_precip_bulk = 1.0  # Constante cinética de precipitação no meio aquoso (1/s)
r_precip_bulk = k_precip_bulk * fe.conditional(fe.gt(SI_UFL, 0.0), SI_UFL * c_safe(C_Fe2), 0.0)

F_total = ( ((C_Fe2 - C_Fe2_n)/dt)*v_Fe2*dx - fe.inner(J_Fe2, fe.grad(v_Fe2))*dx - Fluxo_Fe2_UFL*v_Fe2*ds(1) 
            + 4.0*r_SL_efetivo*v_Fe2*dx + r_precip_bulk*v_Fe2*dx # <-- Consumo de Fe2+ por precipitação volumétrica
          
          + ((C_O2 - C_O2_n)/dt)*v_O2*dx - fe.inner(J_O2, fe.grad(v_O2))*dx - (i_O2_eff/(4.0*F_const))*v_O2*ds(1) + r_SL_efetivo*v_O2*dx
          
          + ((C_AT - C_AT_n)/dt)*v_AT*dx - fe.inner(J_AT, fe.grad(v_AT))*dx - Fluxo_AT_UFL*v_AT*ds(1) 
            + 8.0*r_SL_efetivo*v_AT*dx + 2.0*r_precip_bulk*v_AT*dx # <-- Consumo de Alcalinidade (2 OH- por Fe2+)
          
          + fe.inner(kappa_local * fe.grad(Phi), fe.grad(v_Phi))*dx - i_total_faradaic*v_Phi*ds(1)
          + (fe.Constant(motor.C_dl) * (E_metal - E_metal_n)/dt)*v_E*ds(1) + i_total_faradaic*v_E*ds(1) )

F_DIC = ((C_DIC - C_DIC_n)/dt)*v_DIC*dx - fe.inner(J_DIC, fe.grad(v_DIC))*dx
F_Na  = ((C_Na - C_Na_n)/dt)*v_Na*dx - fe.inner(J_Na, fe.grad(v_Na))*dx
F_NO3 = ((C_NO3 - C_NO3_n)/dt)*v_NO3*dx - fe.inner(J_NO3, fe.grad(v_NO3))*dx - Fluxo_NO3_UFL*v_NO3*ds(1)
F_Fe3 = ((C_Fe3 - C_Fe3_n)/dt)*v_Fe3*dx - fe.inner(J_Fe3, fe.grad(v_Fe3))*dx - 4.0*r_SL_efetivo*v_Fe3*dx

F_total += (F_DIC + F_Na + F_NO3 + F_Fe3)

problema = fe.NonlinearVariationalProblem(F_total, u, bcs=bcs_quimica, J=fe.derivative(F_total, u))
solver = fe.NonlinearVariationalSolver(problema)
# BLINDAGEM DO SOLVER QUÍMICO (Rollback habilitado)
solver.parameters['newton_solver']['linear_solver'] = 'mumps'

# Aumentar o limite para permitir que a descida suave de Newton termine
solver.parameters['newton_solver']['maximum_iterations'] = 50 
solver.parameters['newton_solver']['relaxation_parameter'] = 1.0 # Voltou para 1.0 


# Ajustar tolerâncias para padrões realistas de eletroquímica
solver.parameters['newton_solver']['absolute_tolerance'] = 1e-4
solver.parameters['newton_solver']['relative_tolerance'] = 5e-4

u_n.interpolate(fe.Expression(('1e-12', 'C_sat', '1e-4', '1.4e-2', '100.0', '100.0', '1e-12', '0.0', '-0.210'), C_sat=motor.C_O2_sat, degree=1))
u.assign(u_n)

def extrair_global(vetor, op):
    if len(vetor) > 0: local_val = np.array([vetor.min(), vetor.max()], dtype=np.float64)
    else: local_val = np.array([1e20, -1e20], dtype=np.float64)
    global_val = np.zeros(2, dtype=np.float64)
    if op == "MIN":
        comm.Allreduce(local_val[0:1], global_val[0:1], op=MPI.MIN)
        return global_val[0]
    elif op == "MAX":
        comm.Allreduce(local_val[1:2], global_val[1:2], op=MPI.MAX)
        return global_val[1]

timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M")
csv_filename = f"resultados_2d/Historico_Fase8.2_Puro_{timestamp}.csv"
dados_csv = [["Tempo_h", "OCP_mV", "pH_Min", "pH_Max", "Fe2_Max", "E_eq_Fe_Surf", "h_f_um"]]

tempos_alvo = [0.0, 60.0, 300.0, 600.0, 1800.0, 3600.0]
snapshots_salvos = set()

tempo_atual, tempo_final, passo = 0.0, 3600.0, 0

while tempo_atual <= tempo_final:
    t_start = time.time()  # <--- MARCADOR DE INÍCIO DO PASSO
    
    for t_alvo in tempos_alvo:
        if tempo_atual >= t_alvo and t_alvo not in snapshots_salvos:
            snapshots_salvos.add(t_alvo)
            if rank == 0:
                print(f"📸 Extraindo perfis espaciais para t = {int(t_alvo)}s...")
            
            mask_y = np.abs(dof_coords[:, 1] - H/2.0) < 1e-4
            x_local = dof_coords[mask_y, 0]
            o2_local = u.sub(1, deepcopy=True).vector().get_local()[mask_y]
            h_local = H_fn.vector().get_local()[mask_y]
            
            x_global = comm.gather(x_local, root=0)
            o2_global = comm.gather(o2_local, root=0)
            h_global = comm.gather(h_local, root=0)
            
            if rank == 0:
                x_all = np.concatenate(x_global)
                o2_all = np.concatenate(o2_global)
                h_all = np.concatenate(h_global)
                
                idx = np.argsort(x_all)
                x_sorted = x_all[idx]
                o2_sorted = o2_all[idx]
                ph_sorted = -np.log10(np.maximum(h_all[idx], 1e-16) / 1000.0)
                
                nome_perfil = f"resultados_2d/perfis/Perfil_Linha_t_{int(t_alvo)}s.csv"
                with open(nome_perfil, mode='w', newline='') as f:
                    escritor = csv.writer(f)
                    escritor.writerow(["x_m", "C_O2_mol_m3", "pH"])
                    for x_v, o2_v, ph_v in zip(x_sorted, o2_sorted, ph_sorted):
                        escritor.writerow([x_v, o2_v, ph_v])

    if tempo_atual >= tempo_final:
        break

    dt.assign(dt_val)
    if tempo_atual <= 10.0 or passo % 10 == 0: 
        solver_fluido.solve()
            
    # FORÇAR O SOLVER DIRETO MUMPS A CADA PASSO PARA EVITAR QUE O PETSC REVERTA PARA KRYLOV
    solver.parameters['newton_solver']['linear_solver'] = 'mumps'
    
    # TRATAMENTO DE EXCEÇÃO CASO A MATRIZ DIVERJA NUMERICAMENTE
    try:
        iteracoes, convergiu = solver.solve()
    except RuntimeError:
        convergiu = False
        iteracoes = 999
    
    # -------------------------------------------------------------
    # TRATAMENTO DE ERROS DO NEWTON (ROLLBACK FÍSICO)
    # -------------------------------------------------------------
    if convergiu:
        vetor_Fe2 = u.sub(0, deepcopy=True).vector().get_local()
        vetor_H   = H_fn.vector().get_local()
        vetor_AT  = u.sub(2, deepcopy=True).vector().get_local()
        vetor_DIC = u.sub(3, deepcopy=True).vector().get_local()
        RT_F = (motor.R * motor.T) / motor.F

        vetor_E_Fe  = E_eq_Fe_fn.vector().get_local()
        vetor_E_O2  = E_eq_O2_fn.vector().get_local()
        vetor_E_H2O = E_eq_H2O_fn.vector().get_local()
        
        for dof in range(len(vetor_Fe2)):
            h_local = max(vetor_H[dof], 1e-16) / 1000.0
            vetor_E_O2[dof]  = motor.E_eq_O2_puro + RT_F * math.log(h_local)
            vetor_E_H2O[dof] = motor.E_eq_H2O_puro + RT_F * math.log(h_local)
            
            # Piso contínuo em 1e-9 (Evita o degrau de tensão e a explosão inicial)
            Fe2_safe = max(vetor_Fe2[dof], 1e-9)
            I_local = (2.0 * Fe2_safe + 100.0) / 1000.0 
            g2 = 10**(-0.509 * 4 * (math.sqrt(I_local)/(1+math.sqrt(I_local)) - 0.3*I_local))
            vetor_E_Fe[dof] = motor.E_eq_Fe_padrao + (RT_F/2.0) * math.log(g2 * Fe2_safe / 1000.0)

        E_eq_Fe_fn.vector().set_local(vetor_E_Fe); E_eq_Fe_fn.vector().apply("insert")
        E_eq_O2_fn.vector().set_local(vetor_E_O2); E_eq_O2_fn.vector().apply("insert")
        E_eq_H2O_fn.vector().set_local(vetor_E_H2O); E_eq_H2O_fn.vector().apply("insert")

        novo_H, novo_OH = motor_termo.calcular_pH_alcalinidade(vetor_AT, vetor_DIC, 100.0, vetor_Fe2)
        
        SI_vetor = motor_termo.calcular_SI_FeOH2(vetor_Fe2, novo_H, novo_OH)
        vetor_hf = h_f_fn.vector().get_local()
        k_precip_hf = 1.0e-9 
        k_H_hf = 1.0e-9
        
        for d in dofs_fronteira:
            # A referência passa para 1 mol/m3 e a relação linear evita a explosão quadrática
            f_diss = 0.005 * 1e-9 * max(-SI_vetor[d], 0) + k_H_hf * (novo_H[d] / 1.0)
            dhf = dt_val * (k_precip_hf * max(SI_vetor[d], 0) - f_diss)
            vetor_hf[d] = max(0.0, vetor_hf[d] + dhf)

        H_fn.vector().set_local(novo_H); H_fn.vector().apply("insert")
        OH_fn.vector().set_local(novo_OH); OH_fn.vector().apply("insert")
        h_f_fn.vector().set_local(vetor_hf); h_f_fn.vector().apply("insert")
        
        tempo_atual += dt_val; passo += 1
        u_n.assign(u); w_n.assign(w)
        
        if iteracoes < 100 and dt_val < 0.5: dt_val = min(dt_val * 1.2, 0.5)
        
        ph_min = -math.log10(extrair_global(novo_H, "MAX")/1000.0)
        ph_max = -math.log10(extrair_global(novo_H, "MIN")/1000.0)
        fe2_max = extrair_global(vetor_Fe2, "MAX")
        
        e_metal_integral = fe.assemble(E_metal * ds(1))
        e_metal_val = e_metal_integral / H
        ocp_local = e_metal_val * 1000.0  
        
        if dof_origem != -1:
            vals_local = np.array([vetor_E_Fe[dof_origem], vetor_hf[dof_origem]], dtype=np.float64)
        else:
            vals_local = np.array([-1e20, -1e20], dtype=np.float64)
            
        vals_global_th = np.zeros(2, dtype=np.float64)
        comm.Allreduce(vals_local, vals_global_th, op=MPI.MAX)
        e_eq_superficie = vals_global_th[0]
        hf_superficie_um = vals_global_th[1] * 1e6 
        
        t_end = time.time()  # <--- MARCADOR DE FIM DO PASSO
        
        if rank == 0:
            # Substituímos o print de %10 pelo print em tempo real de cada passo
            print(f"[Progresso] Simulação: {tempo_atual:.4f}s / {tempo_final}s | Passo (dt): {dt_val:.4f}s | Tempo de CPU: {(t_end - t_start):.2f}s | OCP: {ocp_local:.1f} mV | pH max: {ph_max:.2f}")
            dados_csv.append([tempo_atual/3600.0, ocp_local, ph_min, ph_max, fe2_max, e_eq_superficie, hf_superficie_um])
            
    else:
        # ROLLBACK: Se a matriz não convergir ou der erro no PETSc, reduz o dt drasticamente
        dt_val /= 4.0
        u.assign(u_n) 
        w.assign(w_n)
        
        t_end = time.time()  # <--- MARCADOR DE FIM DO PASSO (FALHA)
        
        if rank == 0:
            print(f"⚠️ Instabilidade detectada (CPU: {(t_end - t_start):.2f}s). Restaurando passo anterior e reduzindo dt para {dt_val:.6f}s...")
        if dt_val < 1e-8: 
            if rank == 0: print("❌ Falha crítica: Limite mínimo de dt atingido.")
            break

if rank == 0 and tempo_atual >= tempo_final:
    with open(csv_filename, mode='w', newline='') as f: 
        csv.writer(f).writerows(dados_csv)
        
    tempo_fim_real = time.time()
    tempo_total_seg = tempo_fim_real - tempo_inicio_real
    horas = int(tempo_total_seg // 3600)
    minutos = int((tempo_total_seg % 3600) // 60)
    segundos = tempo_total_seg % 60
    
    print(f"✅ FASE 8.2 FINALIZADA COM SUCESSO. CSV: {csv_filename}")
    print(f"⏱️ Tempo total de execution: {horas}h {minutos}m {segundos:.1f}s")
# =====================================================================
# FASE 2: VARREDURA POTENCIODINÂMICA (3600s a 4200s)
# =====================================================================
if rank == 0:
    print("\n" + "="*60)
    print("=== FASE 2: VARREDURA POTENCIODINÂMICA (POLARIZAÇÃO) ===")
    print("="*60)
    os.makedirs("resultados_2d/graficos", exist_ok=True)

# 1. Configuração da Varredura
v_scan = 0.001  # 1 mV/s
E_start = e_metal_val - 0.300
E_end = e_metal_val + 0.300
tempo_varredura = (E_end - E_start) / v_scan  # 600 segundos
tempo_final_fase2 = tempo_final + tempo_varredura

# 2. Nova Condição de Contorno para E_metal (Potenciostato)
# O potencial no contorno será forçado pela equação: E(t) = E_start + v_scan * (t - t_start)
E_app = fe.Expression('E0 + v*(t - t0)', E0=E_start, v=v_scan, t=tempo_atual, t0=tempo_atual, degree=1)

# RECRIAMOS A FORMULAÇÃO FRACA PARA A FASE 2
# Removemos a dinâmica de capacitor da Fase 1 e forçamos a restrição algébrica: (E_metal - E_app)*v_E*ds(1) = 0
F_total_fase2 = F_total - ((fe.Constant(motor.C_dl) * (E_metal - E_metal_n)/dt)*v_E*ds(1) + i_total_faradaic*v_E*ds(1))
F_total_fase2 += (E_metal - E_app) * v_E * ds(1)

# Recriamos o problema sem o DirichletBC fantasma no potenciostato
problem_fase2 = fe.NonlinearVariationalProblem(F_total_fase2, u, bcs=bcs_quimica, J=fe.derivative(F_total_fase2, u))
solver_fase2 = fe.NonlinearVariationalSolver(problem_fase2)

solver_fase2 = fe.NewtonSolver()
solver_fase2.linear_solver = 'mumps'

# --- NOVOS PARÂMETROS DE ESTABILIZAÇÃO ---
solver_fase2.parameters['newton_solver']['relaxation_parameter'] = 0.5  # Amortece o passo em 50%
solver_fase2.parameters['newton_solver']['absolute_tolerance'] = 1e-6   # Exige maior precisão
solver_fase2.parameters['newton_solver']['relative_tolerance'] = 1e-6
solver_fase2.parameters['newton_solver']['maximum_iterations'] = 50     # Mais tentativas para convergir
# ==========================================

# O laço de tempo da varredura de potencial começa logo depois
while tempo_atual <= tempo_final_fase2:
    
    # Dentro do loop, você deve ter a chamada que usa o solver configurado acima:
    solver_fase2.solve(problem_fase2, u.vector())

while tempo_atual <= tempo_final_fase2:
    t_start_loop = time.time()
    dt.assign(dt_val)
    E_app.t = tempo_atual + dt_val  # Atualiza o potencial aplicado
    
    # ADICIONE ESTA LINHA: Descongela a hidrodinâmica na polarização
    solver_fluido.solve() 
    
    try:
        iteracoes, convergiu = solver_fase2.solve() # Use o solver_fase2 recém criado
    except RuntimeError:
        convergiu = False
        
    if convergiu:
        # Atualizações termodinâmicas (Idêntico à Fase 1)
        vetor_Fe2 = u.sub(0, deepcopy=True).vector().get_local()
        vetor_H   = H_fn.vector().get_local()
        vetor_AT  = u.sub(2, deepcopy=True).vector().get_local()
        vetor_DIC = u.sub(3, deepcopy=True).vector().get_local()
        RT_F = (motor.R * motor.T) / motor.F

        vetor_E_Fe  = E_eq_Fe_fn.vector().get_local()
        vetor_E_O2  = E_eq_O2_fn.vector().get_local()
        vetor_E_H2O = E_eq_H2O_fn.vector().get_local()
        
        for dof in range(len(vetor_Fe2)):
            h_local = max(vetor_H[dof], 1e-16) / 1000.0
            vetor_E_O2[dof]  = motor.E_eq_O2_puro + RT_F * math.log(h_local)
            vetor_E_H2O[dof] = motor.E_eq_H2O_puro + RT_F * math.log(h_local)
            Fe2_safe = max(vetor_Fe2[dof], 1e-9)
            I_local = (2.0 * Fe2_safe + 100.0) / 1000.0 
            g2 = 10**(-0.509 * 4 * (math.sqrt(I_local)/(1+math.sqrt(I_local)) - 0.3*I_local))
            vetor_E_Fe[dof] = motor.E_eq_Fe_padrao + (RT_F/2.0) * math.log(g2 * Fe2_safe / 1000.0)

        E_eq_Fe_fn.vector().set_local(vetor_E_Fe); E_eq_Fe_fn.vector().apply("insert")
        E_eq_O2_fn.vector().set_local(vetor_E_O2); E_eq_O2_fn.vector().apply("insert")
        E_eq_H2O_fn.vector().set_local(vetor_E_H2O); E_eq_H2O_fn.vector().apply("insert")

        novo_H, novo_OH = motor_termo.calcular_pH_alcalinidade(vetor_AT, vetor_DIC, 100.0, vetor_Fe2)
        SI_vetor = motor_termo.calcular_SI_FeOH2(vetor_Fe2, novo_H, novo_OH)
        vetor_hf = h_f_fn.vector().get_local()
        
        for d in dofs_fronteira:
            f_diss = 0.005 * 1e-9 * max(-SI_vetor[d], 0) + 1.0e-9 * (novo_H[d] / 1.0)
            dhf = dt_val * (1.0e-9 * max(SI_vetor[d], 0) - f_diss)
            vetor_hf[d] = max(0.0, vetor_hf[d] + dhf)

        H_fn.vector().set_local(novo_H); H_fn.vector().apply("insert")
        OH_fn.vector().set_local(novo_OH); OH_fn.vector().apply("insert")
        h_f_fn.vector().set_local(vetor_hf); h_f_fn.vector().apply("insert")
        
        # 3. Medição da Corrente (Faradaica + Capacitiva)
        i_farad_integral = fe.assemble(i_total_faradaic * ds(1))
        
        # NOVO CÁLCULO CAPACITIVO ESTÁVEL
        i_cap_constante = motor.C_dl * v_scan
        i_cap_integral = fe.assemble(fe.Constant(i_cap_constante) * ds(1))
        
        i_media = (i_farad_integral + i_cap_integral) / H  # A/m²
        
        tempo_atual += dt_val
        u_n.assign(u)
        
        # Armazenamento Global
        e_imposto = E_app(0.0)
        t_plot.append(tempo_atual)
        E_plot.append(e_imposto * 1000.0) # mV
        I_plot.append(i_media)
        
        # Corrente absoluta para escala logarítmica (evitar erro de domínio)
        logI_plot.append(np.log10(abs(i_media) + 1e-12))
        
        if rank == 0 and int(tempo_atual) % 10 == 0:
            print(f"[Polarização] t = {tempo_atual:.1f}s | E = {e_imposto*1000.0:.1f} mV | i = {i_media:.2e} A/m²")
            
    else:
        dt_val /= 2.0
        u.assign(u_n)
        if dt_val < 1e-3: break

# =====================================================================
# GERAÇÃO DE GRÁFICOS (Rank 0)
# =====================================================================
if rank == 0:
    print("Gerando gráficos da Polarização...")
    import numpy as np
    
    # Recalcula o Log diretamente dos dados puros estabilizados pelo Solver
    logI_plot = np.log10(np.abs(I_plot) + 1e-12)
    
    # 1. Curva de Polarização (Potenciodinâmica)
    plt.figure(figsize=(8, 6))
    plt.plot(logI_plot, E_plot, 'b-', linewidth=2) 
    plt.title("Curva de Polarização Potenciodinâmica", fontsize=14)
    plt.xlabel("Log |i| (A/m²)", fontsize=12)
    plt.ylabel("Potencial vs SHE (mV)", fontsize=12)
    plt.grid(True, which="both", ls="--")
    plt.savefig("resultados_2d/graficos/Curva_Polarizacao.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. Evolução Temporal na Varredura (E vs t)
    fig, ax1 = plt.subplots(figsize=(10, 5))
    ax2 = ax1.twinx()
    
    ax1.plot(t_plot, E_plot, 'g-', linewidth=2, label="Potencial Imposto")
    ax2.plot(t_plot, I_plot, 'r-', linewidth=1.5, label="Corrente Medida")
    
    ax1.set_xlabel("Tempo (s)", fontsize=12)
    ax1.set_ylabel("Potencial vs SHE (mV)", color='g', fontsize=12)
    ax2.set_ylabel("Corrente Total (A/m²)", color='r', fontsize=12)
    plt.title("Estímulo e Resposta durante a Varredura", fontsize=14)
    plt.savefig("resultados_2d/graficos/Evolucao_Temporal_Polarizacao.png", dpi=300, bbox_inches='tight')
    plt.close()