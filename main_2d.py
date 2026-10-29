"""Simulação OCP/condicionamento/polarização com balanços auditáveis.

Execute --help. --smoke percorre as três fases em uma malha pequena;
seus resultados não servem como curva experimental nem estudo de malha.
"""
import argparse
import csv
from dataclasses import asdict
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import sys
import platform
import scipy
import numpy as np
import dolfin  # initialize/finalize MPI through DOLFIN before mpi4py
from mpi4py import MPI
from cinetica import MotorCinetico
from modelo import Modelo, Estado
from transporte import Malha


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--nx',type=int,default=16)
    p.add_argument('--ny',type=int,default=16)
    p.add_argument('--grading',type=float,default=3.)
    p.add_argument('--initial-fe',type=float,default=1e-12,help='Fe(II) total no banho, mol/m³')
    p.add_argument('--initial-no2',type=float,default=1e-12,help='nitrito total inicial, mol/m³; estudar sensibilidade à concentração traço')
    p.add_argument('--ocp',type=float,default=3600.)
    p.add_argument('--conditioning',type=float,default=1.)
    p.add_argument('--span',type=float,default=.3,help='meia largura da varredura, V')
    p.add_argument('--scan-rate',type=float,default=.001,help='V/s')
    p.add_argument('--dt',type=float,default=.001)
    p.add_argument('--dt-max',type=float,default=.5)
    p.add_argument('--dt-min',type=float,default=1e-8)
    p.add_argument('--rtol',type=float,default=None,help='erro temporal relativo; padrão 0.002, smoke 0.05')
    p.add_argument('--profile-times',type=float,nargs='+',help='instantes dos perfis em segundos; padrão: início, 25%, 50%, fim do OCP e fim da varredura')
    p.add_argument('--output',type=Path)
    p.add_argument('--resume',type=Path,help='checkpoint.npz; reutiliza a configuração salva')
    p.add_argument('--density-fe',type=float,default=None,help='d rho/d C_Fe, kg/mol, requer valor medido/justificado')
    p.add_argument('--density-na',type=float,default=None,help='d rho/d C_Na, kg/mol')
    p.add_argument('--smoke',action='store_true')
    return p


def output_snapshot(folder, model, state, phase_index, phase_start, anchor, dt):
    temporary = folder/'checkpoint.tmp.npz'
    fluid_vector=np.array([]) if state.fluid is None else state.fluid.vector().get_local()
    np.savez_compressed(temporary,x=state.x,time=state.time,exchange=state.exchange,
                        initial_delta=model.initial_delta,phase_index=phase_index,
                        phase_start=phase_start,anchor=anchor,dt=dt,fluid=fluid_vector)
    os.replace(temporary,folder/'checkpoint.npz')
    np.savez_compressed(folder/'fields.npz',**field_arrays(model,state,phase_index))


def field_arrays(model,state,phase):
    esp,c,phi,s2,s3,E,h=model.details(state.x)
    names=list(esp.c)
    return dict(x=model.g.x,y=model.g.y,xf=model.g.xf,yf=model.g.yf,
                time_s=state.time,phase=phase,E_SHE_V=E,
                pH=esp.pH.reshape(model.g.shape),phi=phi,
                concentrations=c.reshape((7,*model.g.shape)),
                species_names=np.array(names),
                species_concentrations=np.array([esp.c[name].reshape(model.g.shape) for name in names]),
                FeOH2=s2.reshape(model.g.shape),FeOH3=s3.reshape(model.g.shape))


def save_profiles(folder,model,state,phase,times):
    directory=folder/'snapshots'
    directory.mkdir(exist_ok=True)
    tol=1e-12*max(1.,state.time)
    for index,time in enumerate(times):
        if abs(state.time-time) <= tol:
            np.savez_compressed(directory/f'estado_{index:04d}.npz',**field_arrays(model,state,phase))


def plots(folder,grid):
    from graficos import gerar_graficos
    return gerar_graficos(folder,grid)


def run(args):
    if args.smoke:
        args.nx,args.ny,args.grading = 2,2,1.
        args.ocp,args.conditioning,args.span,args.scan_rate = .002,.002,1e-4,.1
        args.dt,args.dt_max = .001,.001
        if args.rtol is None: args.rtol=.05
        args.initial_fe,args.initial_no2 = 1e-6,1e-6
    if args.resume:
        folder = args.resume.resolve().parent
        saved = json.loads((folder/'config.json').read_text())
        for key,value in saved['arguments'].items():
            if key not in ('resume','output'):
                setattr(args,key,value)
    else:
        folder = args.output or Path('resultados_2d')/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
        if folder.exists() and any(folder.iterdir()):
            raise ValueError('Diretório de saída contém dados; use outro ou --resume checkpoint.npz')
        folder.mkdir(parents=True,exist_ok=True)
    if args.rtol is None: args.rtol=.002
    if args.initial_fe <= 0 or args.initial_no2 <= 0:
        raise ValueError('Concentrações traço Fe/NO2 devem ser positivas para Nernst')
    if not all(np.isfinite(v) for v in (args.ocp,args.conditioning,args.span,args.scan_rate,args.dt,args.dt_max,args.dt_min,args.rtol)):
        raise ValueError('Parâmetros devem ser finitos')
    if min(args.ocp,args.conditioning,args.span,args.scan_rate,args.dt,args.dt_max,args.dt_min,args.rtol) <= 0:
        raise ValueError('Tempos, amplitudes e tolerâncias devem ser positivos')
    if args.dt_min > args.dt or args.dt > args.dt_max:
        raise ValueError('Exige dt_min <= dt <= dt_max')
    total_time=args.ocp+args.conditioning+2*args.span/args.scan_rate
    if args.profile_times is None:
        args.profile_times=[0.,.25*args.ocp,.5*args.ocp,args.ocp,total_time]
    if any(not np.isfinite(t) or t<0 or t>total_time for t in args.profile_times):
        raise ValueError('Instantes dos perfis devem estar entre zero e o fim da simulação')
    args.profile_times=sorted(set([0.,*args.profile_times,total_time]))
    grid = Malha(args.nx,args.ny,grading=args.grading)
    flow = None
    if args.density_fe is not None or args.density_na is not None:
        from escoamento import Escoamento
        coeff = {i:v for i,v in ((0,args.density_fe),(3,args.density_na)) if v is not None}
        flow = Escoamento(grid,coeff)
    motor=MotorCinetico(**saved["kinetics"]) if args.resume else MotorCinetico()
    model = Modelo(grid,flow=flow,motor=motor,bath_fe=args.initial_fe,bath_no2=args.initial_no2)
    state = model.initial()
    inv0 = model.inventory(state)
    pi,start,anchor,dt = 0,0.,model.details(state.x)[5],args.dt
    if args.resume:
        with np.load(args.resume) as f:
            state=Estado(f['x'],float(f['time']),exchange=f['exchange'])
            model.initial_delta=f['initial_delta']
            pi,start,anchor,dt = int(f['phase_index']),float(f['phase_start']),float(f['anchor']),float(f['dt'])
        if flow is not None:
            state.fluid=flow.fe.Function(flow.W)
            with np.load(args.resume) as f:
                state.fluid.vector().set_local(f['fluid'])
                state.fluid.vector().apply('insert')
        inv0=np.array(saved['initial_inventory'])
    else:
        config=dict(arguments={k:(str(v) if isinstance(v,Path) else v) for k,v in vars(args).items()},
                    kinetics=asdict(model.m),initial_inventory=inv0.tolist(),
                    versions=dict(python=platform.python_version(),numpy=np.__version__,
                                  scipy=scipy.__version__,dolfin=dolfin.__version__,
                                  mpi=MPI.Get_library_version().strip()),
                    formulation='electroneutral components / SG / implicit reaction / compact-layer nitrate inventory',
                    potential_reference='SHE, reservoir Phi=0',units='SI; activities relative to 1 mol/L',
                    validation='model parameters inherited/indicative; not experimentally calibrated')
        (folder/'config.json').write_text(json.dumps(config,indent=2))
    # Last-column cell distances and diffusivity set an explicit resolution warning.
    diffusion_length=np.sqrt(2.1e-9*args.dt)
    resolved = grid.dx[0] <= diffusion_length
    if not resolved:
        print('AVISO: primeira célula maior que sqrt(D_O2*dt); refine --nx/--grading antes de interpretar o transiente.',flush=True)
    durations=[args.ocp,args.conditioning,2*args.span/args.scan_rate]
    names=['OCP','condicionamento','polarizacao']
    header=['time_s','phase','dt_s','E_SHE_V','j_far_A_m2','j_cap_A_m2','j_total_A_m2',
            'pH_min','pH_max','Fe_balance_mol_m','N_balance_mol_m','C_balance_mol_m','Na_balance_mol_m',
            'scaled_residual','temporal_error','Re','density_ratio']
    success=False
    save_profiles(folder,model,state,pi,args.profile_times)
    with open(folder/'history.csv','a' if args.resume else 'w',newline='') as history:
        writer=csv.writer(history)
        if not args.resume: writer.writerow(header)
        try:
            while pi < 3:
                end = start+durations[pi]
                if state.time >= end-1e-13*max(1.,end):
                    if pi == 0: anchor=model.details(state.x)[5]
                    start=end; pi+=1
                    output_snapshot(folder,model,state,pi,start,anchor,dt)
                    continue
                next_profile=next((t for t in args.profile_times if t>state.time+1e-13*max(1.,state.time)),end)
                used=min(dt,end-state.time,next_profile-state.time)
                def target(t):
                    if pi == 0: return None
                    if pi == 1: return anchor-args.span*(t-start)/durations[1]
                    return anchor-args.span+args.scan_rate*(t-start)
                try:
                    coarse=model.solve(state,used,target(state.time+used))
                    mid=model.solve(state,used/2,target(state.time+used/2))
                    fine=model.solve(mid,used/2,target(state.time+used))
                    err=model.error_norm(coarse,fine,rtol=args.rtol)
                    if not np.isfinite(err) or err > 1:
                        raise RuntimeError(f'Erro temporal estimado {err:.3g} excede 1')
                except (RuntimeError,ValueError,FloatingPointError) as exc:
                    dt=used/2
                    print(f'Passo rejeitado em t={state.time:.6g}: {exc}; dt={dt:.3g}',flush=True)
                    if dt < args.dt_min:
                        raise RuntimeError('dt mínimo atingido; histórico e checkpoint preservados') from exc
                    continue
                # Log average fluxes over BOTH accepted half steps, not recalculated state currents.
                di,df=mid.info,fine.info
                state=fine
                save_profiles(folder,model,state,pi,args.profile_times)
                esp,_,_,_,_,E,_=model.details(state.x)
                inv=model.inventory(state)
                expected=np.array([state.exchange[0]+state.exchange[7],state.exchange[8],state.exchange[2],state.exchange[3]])
                mass=inv-inv0-expected
                jf=.5*np.mean(di['faradaic']+df['faradaic'])
                jc=.5*np.mean(di['capacitive']+df['capacitive'])
                writer.writerow([state.time,pi,used,E,jf,jc,jf+jc,float(esp.pH.min()),float(esp.pH.max()),
                                 *mass,model.diagnostics['max_residual'],err,
                                 model.diagnostics['Re'],model.diagnostics['density_ratio']])
                history.flush()
                print(f'{names[pi]} t={state.time:.6g}s E={E:.6f} V_SHE j={jf+jc:.3e} A/m² erro={err:.2g}',flush=True)
                dt=min(args.dt_max,used*min(1.5,max(.5,.85/np.sqrt(max(err,1e-8)))))
                output_snapshot(folder,model,state,pi,start,anchor,dt)
            success=True
        finally:
            output_snapshot(folder,model,state,pi,start,anchor,dt)
            (folder/'status.json').write_text(json.dumps(dict(completed=success,time_s=state.time,phase=pi,
                near_wall_initial_resolution=bool(resolved),experimental_validation=False),indent=2))
            plots(folder,grid)
    return folder


def main(argv=None):
    args=parser().parse_args(argv)
    comm=MPI.COMM_WORLD
    code=0
    if comm.rank == 0:
        try:
            folder=run(args)
            print(f'Concluído. Saídas: {folder}',flush=True)
        except Exception as exc:
            import traceback
            traceback.print_exc()
            print(f'ERRO: {exc}',file=sys.stderr,flush=True)
            code=1
    # Single deterministic coordinator; workers do not write or execute solves.
    code=comm.bcast(code,root=0)
    return code


if __name__ == '__main__':
    sys.exit(main())
