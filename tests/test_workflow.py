import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[1]


class Workflow(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp=tempfile.TemporaryDirectory()
        cls.base=Path(cls.tmp.name)
        cls.serial=cls.base/'serial'
        r=cls.invoke('--smoke','--output',str(cls.serial))
        if r.returncode: raise RuntimeError(r.stdout)

    @classmethod
    def tearDownClass(cls): cls.tmp.cleanup()

    @staticmethod
    def invoke(*args,mpi=False):
        cmd=[sys.executable,str(ROOT/'main_2d.py'),*args]
        if mpi: cmd=['mpirun','--mca','btl','self,tcp','--mca','pml','ob1','-n','2',*cmd]
        env={k:v for k,v in os.environ.items() if not k.startswith(('OMPI_','PMIX_','PMI_','OPAL_'))}
        return subprocess.run(cmd,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=180)

    def test_complete_and_conservative(self):
        folder=self.serial
        self.assertTrue(json.loads((folder/'status.json').read_text())['completed'])
        data=np.genfromtxt(folder/'history.csv',delimiter=',',names=True)
        self.assertEqual(set(data['phase']),{0,1,2})
        self.assertAlmostEqual(data['time_s'][-1],.006,places=12)
        self.assertLess(np.max(data['temporal_error']),1.)
        for key in ('Fe_balance_mol_m','N_balance_mol_m','C_balance_mol_m','Na_balance_mol_m'):
            self.assertLess(np.max(np.abs(data[key])),1e-9)
        for filename in ('checkpoint.npz','polarizacao.png','mapa_pH.png','perfil.csv'):
            self.assertGreater((folder/filename).stat().st_size,0)

    def test_mpi_matches_serial(self):
        out=self.base/'mpi'
        result=self.invoke('--smoke','--output',str(out),mpi=True)
        self.assertEqual(result.returncode,0,result.stdout)
        a=np.load(self.serial/'checkpoint.npz'); b=np.load(out/'checkpoint.npz')
        np.testing.assert_allclose(a['x'],b['x'],rtol=1e-10,atol=1e-11)
        self.assertEqual((self.serial/'history.csv').read_bytes(),(out/'history.csv').read_bytes())

    def test_failure_checkpoint_and_resume(self):
        out=self.base/'recover'
        result=self.invoke('--smoke','--dt-min','.0009','--rtol','.000001','--output',str(out))
        self.assertNotEqual(result.returncode,0,result.stdout)
        self.assertFalse(json.loads((out/'status.json').read_text())['completed'])
        self.assertTrue((out/'checkpoint.npz').exists())
        config=json.loads((out/'config.json').read_text())
        # A recovery changes only a numerical bound, never the physical state.
        config['arguments']['dt_min']=1e-8
        config['arguments']['rtol']=.05
        (out/'config.json').write_text(json.dumps(config))
        result=self.invoke('--resume',str(out/'checkpoint.npz'))
        self.assertEqual(result.returncode,0,result.stdout)
        self.assertTrue(json.loads((out/'status.json').read_text())['completed'])
        data=np.genfromtxt(out/'history.csv',delimiter=',',names=True)
        self.assertAlmostEqual(data['time_s'][-1],.006,places=12)
        self.assertTrue(np.all(np.diff(data['time_s'])>0))


if __name__=='__main__': unittest.main()
