"""Gráficos e CSV a partir de resultados salvos, sem repetir a simulação.

Uso: python3 graficos.py DIRETORIO_DOS_RESULTADOS
Concentrações em mol/m³; potenciais vs SHE. Log de corrente é log10(|j|).
"""
import argparse
import json
from html import escape
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

COMPONENTES = ('Fe(II) total', 'O₂', 'DIC', 'Na⁺', 'NO₃⁻', 'Nitrito total', 'Alcalinidade')
CHAVES = ('FeII_total', 'O2', 'DIC', 'Na', 'NO3', 'NO2_total', 'AT')


def log_corrente(j):
    """Zeros não têm logaritmo: NaN no CSV e lacuna no gráfico, sem piso fictício."""
    j = np.asarray(j, float)
    result = np.full(j.shape, np.nan)
    nonzero = np.isfinite(j) & (j != 0)
    result[nonzero] = np.log10(np.abs(j[nonzero]))
    return result


def perfil_meia_altura(array, y, altura):
    """Interpolação entre centros dos volumes em y=H/2, sem escolher uma linha acima."""
    a = np.asarray(array)
    return np.array([np.interp(altura/2, y, a[:,i]) for i in range(a.shape[1])])


def _salvar(fig, folder, nome):
    fig.tight_layout()
    fig.savefig(folder/nome, dpi=180, bbox_inches='tight')
    plt.close(fig)


def _csv(path, columns, header):
    np.savetxt(path, np.column_stack(columns), delimiter=',', header=','.join(header), comments='')


def _campos(path):
    with np.load(path, allow_pickle=False) as f:
        return {k:f[k] for k in f.files}


def _perfil(campos, altura):
    return lambda a: perfil_meia_altura(a, campos['y'], altura)


def gerar_graficos(folder, grid=None):
    folder = Path(folder)
    config_path = folder/'config.json'
    config = json.loads(config_path.read_text()) if config_path.exists() else {}
    fields = _campos(folder/'fields.npz')
    if grid is None and not all(k in fields for k in ('xf','yf')):
        from transporte import Malha
        args = config['arguments']
        grid = Malha(args['nx'], args['ny'], grading=args.get('grading',3.))
    xf = fields['xf'] if 'xf' in fields else grid.xf
    yf = fields['yf'] if 'yf' in fields else grid.yf
    altura = float(yf[-1])
    status_path = folder/'status.json'
    completed = json.loads(status_path.read_text()).get('completed',False) if status_path.exists() else False
    nota = ' (teste curto)' if config.get('arguments',{}).get('smoke',False) else ''
    if not completed: nota+=' (execução parcial)'
    data = np.atleast_1d(np.genfromtxt(folder/'history.csv', delimiter=',', names=True))
    if len(data):
        fig, ax = plt.subplots(figsize=(8,5))
        for phase, label in ((0,'OCP'),(1,'Condicionamento'),(2,'Polarização')):
            d = data[data['phase']==phase]
            if len(d): ax.plot(d['time_s'],1000*d['E_SHE_V'],'.-',label=label)
        ax.set(xlabel='Tempo (s)',ylabel='Potencial vs SHE (mV)',title='Potencial em todas as fases'+nota)
        ax.grid(True,alpha=.3); ax.legend()
        _salvar(fig,folder,'potencial.png')
        ocp = data[data['phase']==0]
        if len(ocp):
            fig,ax = plt.subplots(figsize=(8,5))
            t_ocp,E_ocp=ocp['time_s'],ocp['E_SHE_V']
            initial=folder/'snapshots'/'estado_0000.npz'
            if initial.exists():
                start=_campos(initial)
                if float(start['time_s'])==0. and t_ocp[0]>0:
                    t_ocp=np.r_[0.,t_ocp]; E_ocp=np.r_[float(start['E_SHE_V']),E_ocp]
            ax.plot(t_ocp,1000*E_ocp,'.-')
            ax.set(xlabel='Tempo (s)',ylabel='OCP vs SHE (mV)',title='Potencial de circuito aberto'+nota)
            ax.grid(True,alpha=.3); _salvar(fig,folder,'ocp_vs_tempo.png')
            _csv(folder/'ocp.csv',[t_ocp,E_ocp],['time_s','E_SHE_V'])
        fig,ax=plt.subplots(figsize=(8,5))
        for key,label in (('j_far_A_m2','Faradaica'),('j_cap_A_m2','Capacitiva'),('j_total_A_m2','Total')):
            ax.plot(data['time_s'],data[key],'.-',label=label)
        ax.set(xlabel='Tempo (s)',ylabel='Densidade de corrente (A/m²)',title='Correntes assinadas'+nota)
        ax.legend(); ax.grid(True,alpha=.3); _salvar(fig,folder,'correntes_vs_tempo.png')
        fig,ax=plt.subplots(figsize=(8,5))
        ax.plot(data['time_s'],data['pH_min'],label='Mínimo')
        ax.plot(data['time_s'],data['pH_max'],label='Máximo')
        ax.set(xlabel='Tempo (s)',ylabel='pH de atividade',title='Faixa de pH no domínio'+nota)
        ax.legend(); ax.grid(True,alpha=.3); _salvar(fig,folder,'pH_vs_tempo.png')
        scan=data[data['phase']==2]
        if len(scan):
            E=1000*scan['E_SHE_V']; j=scan['j_total_A_m2']; logs=log_corrente(j)
            farlogs=log_corrente(scan['j_far_A_m2'])
            _csv(folder/'polarizacao.csv',[scan['time_s'],scan['E_SHE_V'],j,scan['j_far_A_m2'],scan['j_cap_A_m2'],logs,farlogs],
                 ['time_s','E_SHE_V','j_total_A_m2','j_far_A_m2','j_cap_A_m2','log10_abs_j_total','log10_abs_j_far'])
            fig,ax=plt.subplots(figsize=(8,5))
            ax.plot(j,E,'.-'); ax.set_xscale('symlog',linthresh=1e-8)
            ax.set(xlabel='Densidade de corrente assinada (A/m²)',ylabel='Potencial vs SHE (mV)',title='Polarização'+nota)
            ax.grid(True,alpha=.3); _salvar(fig,folder,'polarizacao.png')
            fig,ax=plt.subplots(figsize=(8,5))
            for positive,label in ((True,'Total anódica'),(False,'Total catódica')):
                yy=np.where((j>0) if positive else (j<0),logs,np.nan)
                if np.any(np.isfinite(yy)): ax.plot(E,yy,'.-',label=label)
            if np.any(np.isfinite(farlogs)): ax.plot(E,farlogs,'--',label='Faradaica |jF|')
            ax.set(xlabel='Potencial vs SHE (mV)',ylabel='log₁₀(|j| / (A/m²))',title='Log da corrente vs potencial'+nota)
            ax.legend(); ax.grid(True,alpha=.3); _salvar(fig,folder,'log_i_vs_E.png')
            fig,ax=plt.subplots(figsize=(8,5))
            ax.plot(logs,E,'.-',label='Total'); ax.plot(farlogs,E,'--',label='Faradaica')
            ax.set(xlabel='log₁₀(|j| / (A/m²))',ylabel='Potencial vs SHE (mV)',title='Polarização: apresentação tipo Tafel'+nota)
            ax.legend(); ax.grid(True,alpha=.3); _salvar(fig,folder,'E_vs_log_i.png')
    # Perfis ao longo do tempo, nos instantes efetivamente armazenados.
    snapshots=sorted((folder/'snapshots').glob('*.npz'))
    samples=[_campos(p) for p in snapshots] if snapshots else [fields]
    samples.sort(key=lambda f:float(f.get('time_s',0.)))
    fig,axes=plt.subplots(2,3,figsize=(13,8))
    for index,sample in enumerate(samples):
        profile=_perfil(sample,altura)
        profile_dir=folder/'perfis_temporais'
        profile_dir.mkdir(exist_ok=True)
        t=float(sample.get('time_s',np.nan))
        _csv(profile_dir/f'perfil_{index:04d}.csv',
             [np.full(len(sample['x']),t),sample['x'],profile(sample['pH'])]+[profile(c) for c in sample['concentrations']],
             ['time_s','x_m','pH']+[f'C_{key}_mol_m3' for key in CHAVES])
        label=f"t={float(sample.get('time_s',0.)):.6g} s" if 'time_s' in sample else 'Estado final'
        for k,ax in enumerate(axes.flat):
            ax.plot(sample['x']*1000,profile(sample['concentrations'][k]),'.-',label=label)
            ax.set(title=COMPONENTES[k],xlabel='Distância ao eletrodo (mm)',ylabel='Concentração (mol/m³)')
            ax.set_ylim(bottom=0)
            ax.grid(True,alpha=.3)
    axes.flat[0].legend(fontsize=8)
    fig.suptitle('Perfis na meia altura (y=H/2)'+nota)
    _salvar(fig,folder,'perfis_concentracao.png')
    fig,axes=plt.subplots(1,2,figsize=(11,4))
    for sample in samples:
        p=_perfil(sample,altura); label=f"t={float(sample.get('time_s',0.)):.6g} s" if 'time_s' in sample else 'Estado final'
        axes[0].plot(sample['x']*1000,p(sample['pH']),'.-',label=label)
        axes[1].plot(sample['x']*1000,1000*p(sample['phi']),'.-',label=label)
    axes[0].set(ylabel='pH de atividade'); axes[1].set(ylabel='Potencial da solução (mV)')
    for ax in axes: ax.set_xlabel('Distância ao eletrodo (mm)'); ax.grid(True,alpha=.3)
    axes[0].legend(fontsize=8); _salvar(fig,folder,'perfis_pH_potencial.png')
    p=_perfil(fields,altura)
    # Preserve the three original columns and append the remaining components.
    cols=[fields['x'],p(fields['pH']),p(fields['concentrations'][1])]
    names=['x_m','pH','C_O2_mol_m3']
    for k in (0,2,3,4,5,6): cols.append(p(fields['concentrations'][k])); names.append(f'C_{CHAVES[k]}_mol_m3')
    cols.extend((p(fields['phi']),p(fields['FeOH2']),p(fields['FeOH3'])))
    names.extend(('Phi_V','FeOH2_solid_mol_m3','FeOH3_solid_mol_m3'))
    _csv(folder/'perfil.csv',cols,names)
    if 'species_concentrations' in fields:
        names=['x_m']+[f'C_{s}_mol_m3' for s in fields['species_names']]
        _csv(folder/'perfil_especies.csv',[fields['x']]+[p(c) for c in fields['species_concentrations']],names)
        fig,axes=plt.subplots(2,2,figsize=(11,7))
        species={str(s):c for s,c in zip(fields['species_names'],fields['species_concentrations'])}
        for ax,nameset,title in zip(axes.flat,(('Fe','FeOH'),('H','OH'),('CO2','HCO3','CO3'),('NO2','HNO2')),
                                    ('Ferro livre e complexado','H⁺ e OH⁻','Especiação do carbono','Especiação do nitrito')):
            for name in nameset: ax.plot(fields['x']*1000,p(species[name]),'.-',label=name)
            ax.set(title=title,xlabel='Distância ao eletrodo (mm)',ylabel='Concentração (mol/m³)')
            ax.legend(); ax.grid(True,alpha=.3)
        _salvar(fig,folder,'perfis_especies.png')
    for key,title,unit,filename in (('pH','pH no domínio','pH de atividade','mapa_pH.png'),
                                  ('phi','Potencial da solução','V','mapa_potencial.png')):
        fig,ax=plt.subplots(figsize=(6,5))
        im=ax.pcolormesh(xf*1000,yf*1000,fields[key],shading='flat')
        fig.colorbar(im,ax=ax,label=unit); ax.set(title=title+nota,xlabel='x (mm)',ylabel='y (mm)')
        _salvar(fig,folder,filename)
    fig,axes=plt.subplots(2,3,figsize=(13,8))
    for k,ax in enumerate(axes.flat):
        im=ax.pcolormesh(xf*1000,yf*1000,fields['concentrations'][k],shading='flat')
        fig.colorbar(im,ax=ax,label='mol/m³'); ax.set(title=COMPONENTES[k],xlabel='x (mm)',ylabel='y (mm)')
    fig.suptitle('Concentrações no último estado salvo'+nota); _salvar(fig,folder,'mapas_concentracao.png')
    fig,axes=plt.subplots(1,2,figsize=(11,4))
    for ax,key in zip(axes,('FeOH2','FeOH3')):
        im=ax.pcolormesh(xf*1000,yf*1000,fields[key],shading='flat')
        fig.colorbar(im,ax=ax,label='mol/m³'); ax.set(title=key+' sólido',xlabel='x (mm)',ylabel='y (mm)')
    _salvar(fig,folder,'mapas_solidos.png')
    manifest=dict(images=sorted(p.name for p in folder.glob('*.png')),completed=bool(completed),
                  profile_times_s=[float(f['time_s']) for f in samples if 'time_s' in f],
                  log_current='log10(abs(j)/(A/m2)); zero -> NaN',
                  profile_y_m=altura/2)
    (folder/'graficos.json').write_text(json.dumps(manifest,indent=2))
    demonstration=config.get('arguments',{}).get('smoke',False)
    label='Demonstração curta' if demonstration else 'Resultados da simulação'
    if not completed: label+=' — execução parcial'
    cards=''.join(f'<figure><a href="{escape(name)}"><img loading="lazy" src="{escape(name)}" alt="{escape(name)}"></a><figcaption>{escape(name)}</figcaption></figure>' for name in manifest['images'])
    links=''.join(f'<li><a href="{escape(p.name)}">{escape(p.name)}</a></li>' for p in sorted(folder.glob('*.csv')))
    html=f'''<!doctype html><html lang="pt-BR"><meta charset="utf-8">
<title>Gráficos do simulador</title><style>
body{{font-family:system-ui,sans-serif;max-width:1100px;margin:2rem auto;padding:0 1rem;background:#f5f7fa;color:#172234}}
figure{{margin:2rem 0;background:white;padding:1rem;border-radius:8px}}img{{max-width:100%;height:auto}}
a{{color:#1656a6}}figcaption{{font-size:.9rem}}h1{{font-size:1.8rem}}</style>
<h1>{escape(label)}</h1><p>Potenciais vs SHE, concentrações em mol/m³ e densidades de corrente em A/m².
Os gráficos não constituem validação experimental. Log(i) significa log10(|j|); valores zero não têm logaritmo.</p>
<p>Clique em uma figura para abrir a imagem original. Os perfis temporais estão também em CSV na pasta perfis_temporais.</p>
<h2>Dados para download</h2><ul>{links}</ul>{cards}</html>'''
    (folder/'index.html').write_text(html,encoding='utf-8')
    return manifest


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('directory',type=Path)
    args=p.parse_args()
    result=gerar_graficos(args.directory)
    print(f"{len(result['images'])} gráficos salvos em {args.directory}")
