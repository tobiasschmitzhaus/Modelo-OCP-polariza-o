import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
from numpy.testing import assert_allclose
from graficos import gerar_graficos,log_corrente,perfil_meia_altura
from transporte import Malha
from termodinamica import ESPECIES


class Graphics(unittest.TestCase):
    def test_signed_log_zero_and_nonfinite(self):
        actual=log_corrente(np.array([-1e-3,0.,1e2,np.nan,np.inf]))
        assert_allclose(actual,[-3,np.nan,2,np.nan,np.nan],equal_nan=True)

    def test_mid_height_interpolation(self):
        x=np.array([.001,.003,.008]); y=np.array([.004,.015,.018])
        f=2*x[None,:]+3*y[:,None]
        assert_allclose(perfil_meia_altura(f,y,.02),2*x+3*.01)

    def test_exports_match_saved_fields_and_history(self):
        with tempfile.TemporaryDirectory() as directory:
            folder=Path(directory); g=Malha(3,2,grading=1.)
            c=np.array([k+1+g.x[None,:]+2*g.y[:,None] for k in range(7)])
            ph=6+np.zeros(g.shape); phi=np.zeros(g.shape)
            fields=dict(x=g.x,y=g.y,xf=g.xf,yf=g.yf,time_s=3.,E_SHE_V=-.4,
                        pH=ph,phi=phi,concentrations=c,FeOH2=phi,FeOH3=phi,
                        species_names=np.array(list(ESPECIES)),
                        species_concentrations=np.array([np.full(g.shape,1e-6) for _ in ESPECIES]))
            np.savez_compressed(folder/'fields.npz',**fields)
            (folder/'snapshots').mkdir()
            np.savez_compressed(folder/'snapshots'/'estado_0000.npz',**dict(fields,time_s=0.,E_SHE_V=-.5))
            np.savez_compressed(folder/'snapshots'/'estado_0001.npz',**fields)
            (folder/'config.json').write_text(json.dumps({'arguments':{'nx':3,'ny':2,'grading':1.}}))
            (folder/'status.json').write_text(json.dumps({'completed':True}))
            (folder/'history.csv').write_text('time_s,phase,E_SHE_V,j_total_A_m2,j_far_A_m2,j_cap_A_m2,pH_min,pH_max\n'
                '.5,0,-.5,0,0,0,6,7\n1,2,-.6,-.01,-.01,0,6,7\n2,2,-.5,0,0,0,6,7\n3,2,-.4,.1,.1,0,6,7\n')
            result=gerar_graficos(folder)
            self.assertTrue(result['completed'])
            self.assertEqual(result['profile_times_s'],[0.,3.])
            for name in ('ocp_vs_tempo.png','log_i_vs_E.png','E_vs_log_i.png','perfis_concentracao.png',
                         'perfis_especies.png','mapas_concentracao.png','mapas_solidos.png'):
                self.assertIn(name,result['images'])
                self.assertGreater((folder/name).stat().st_size,2000)
                self.assertEqual((folder/name).read_bytes()[:8],b'\x89PNG\r\n\x1a\n')
            curve=np.genfromtxt(folder/'polarizacao.csv',delimiter=',',names=True)
            assert_allclose(curve['E_SHE_V'],[-.6,-.5,-.4])
            assert_allclose(curve['log10_abs_j_total'],[-2,np.nan,-1],equal_nan=True)
            profile=np.genfromtxt(folder/'perfil.csv',delimiter=',',names=True)
            assert_allclose(profile['C_FeII_total_mol_m3'],1+g.x+g.H)
            assert_allclose(profile['C_O2_mol_m3'],2+g.x+g.H)
            self.assertTrue((folder/'perfil_especies.csv').exists())
            ocp=np.genfromtxt(folder/'ocp.csv',delimiter=',',names=True)
            assert_allclose(ocp['time_s'],[0.,.5])


if __name__=='__main__': unittest.main()
