"""Original synthetic regressions: no captured rooms or third-party datasets."""
import copy, csv, io, json, os, random, subprocess, sys, tempfile, unittest
from pathlib import Path
from spatial_witness.core import *
ROOT=Path(__file__).resolve().parents[1]
def cfg():return {'schema':1,'frame':'test','reference_unit':'m','reconstruction_unit':'m','threshold_um':1000000,'minimum_precision_bp':10000,'minimum_recall_bp':10000,'regions':[]}
def points(rows):return ('x,y,z\n'+''.join(','.join(map(str,r))+'\n' for r in rows)).encode()
def run(r,p,c=None):return evaluate(canonical(c or cfg()),points(r),points(p))
class Geometry(unittest.TestCase):
    def test_identical(self):
        r=run([(0,0,0),(1,2,3)],[(0,0,0),(1,2,3)]);self.assertEqual(r['status'],'pass');self.assertEqual(r['fscore'],1)
    def test_strict_boundary(self):
        r=run([(0,0,0)],[(1,0,0)]);self.assertEqual(r['recall']['matched'],0);self.assertEqual(r['status'],'fail_or_unknown')
    def test_inside_micro(self):
        self.assertEqual(run([(0,0,0)],[('0.999999',0,0)])['status'],'pass')
    def test_outside_micro(self):self.assertEqual(run([(0,0,0)],[('1.000001',0,0)])['precision']['matched'],0)
    def test_units(self):
        c=cfg();c['reconstruction_unit']='mm';self.assertEqual(run([(1,2,3)],[(1000,2000,3000)],c)['status'],'pass')
    def test_tie_first(self):self.assertEqual(run([(-1,0,0),(1,0,0)],[(0,0,0)])['reconstruction_to_reference'][0]['nearest_index'],0)
    def test_precision_recall_distinct(self):
        r=run([(0,0,0),(10,0,0)],[(0,0,0)]);self.assertEqual(r['precision']['ratio'],1);self.assertEqual(r['recall']['ratio'],.5);self.assertAlmostEqual(r['fscore'],2/3)
    def test_outlier(self):
        r=run([(0,0,0)],[(0,0,0),(10,0,0)]);self.assertEqual(r['precision']['ratio'],.5);self.assertEqual(r['recall']['ratio'],1)
    def test_duplicates_preserve_weight(self):
        r=run([(0,0,0),(0,0,0),(10,0,0)],[(0,0,0)]);self.assertEqual(r['reference']['duplicate_rows'],1);self.assertEqual(r['recall']['matched'],2)
    def test_empty_reconstruction(self):
        r=run([(0,0,0)],[]);self.assertIsNone(r['precision']['ratio']);self.assertIsNone(r['fscore']);self.assertIsNone(r['reference_to_reconstruction'][0]['nearest_index']);self.assertEqual(r['status'],'fail_or_unknown')
    def test_no_reference(self):
        with self.assertRaises(Invalid):run([],[(0,0,0)])
    def test_zero_gate_not_empty_pass(self):
        c=cfg();c['minimum_precision_bp']=c['minimum_recall_bp']=0;self.assertEqual(run([(0,0,0)],[],c)['status'],'fail_or_unknown')
    def test_region_empty_unknown(self):
        c=cfg();c['regions']=[{'id':'empty','min_um':[10,10,10],'max_um':[20,20,20],'minimum_recall_bp':0}]
        r=run([(0,0,0)],[(0,0,0)],c);self.assertEqual(r['regions'][0]['status'],'unknown_no_reference_samples');self.assertEqual(r['status'],'fail_or_unknown')
    def test_region_inclusive_and_all_candidates(self):
        c=cfg();c['regions']=[{'id':'edge','min_um':[0,0,0],'max_um':[0,0,0],'minimum_recall_bp':10000}]
        r=run([(0,0,0)],[(0.5,0,0)],c);self.assertEqual(r['regions'][0]['reference_indices'],[0]);self.assertEqual(r['regions'][0]['status'],'pass')
    def test_exact_basis_point_gate(self):
        c=cfg();c['minimum_recall_bp']=6666;r=run([(0,0,0),(10,0,0),(20,0,0)],[(0,0,0),(10,0,0)],c);self.assertEqual(r['status'],'pass')
        c['minimum_recall_bp']=6667;self.assertEqual(run([(0,0,0),(10,0,0),(20,0,0)],[(0,0,0),(10,0,0)],c)['status'],'fail_or_unknown')
    def test_extreme_coordinate(self):
        r=run([(-1000000,-1000000,-1000000)],[(1000000,1000000,1000000)]);self.assertEqual(r['reconstruction_to_reference'][0]['distance_squared_um2'],12*10**24)
    def test_seeded_oracle(self):
        rng=random.Random(712)
        for _ in range(100):
            a=[tuple(rng.randint(-5,5) for _ in range(3)) for j in range(rng.randint(1,12))];b=[tuple(rng.randint(-5,5) for _ in range(3)) for j in range(rng.randint(1,12))]
            r=run(a,b)
            for key,left,right in [('reconstruction_to_reference',b,a),('reference_to_reconstruction',a,b)]:
                for p,got in zip(left,r[key]):
                    ds=[((p[0]-q[0])**2+(p[1]-q[1])**2+(p[2]-q[2])**2)*10**12 for q in right]
                    self.assertEqual(got,{'nearest_index':ds.index(min(ds)),'distance_squared_um2':min(ds)})
    def test_max_points(self):
        r=run([(i,0,0) for i in range(1000)],[(0,0,0)]);self.assertEqual(r['recall']['total'],1000);self.assertLess(len(outputs(r)['report.json']),OUTPUT_LIMIT)
class Input(unittest.TestCase):
    def test_coordinate_rejections(self):
        for x in ['nan','inf','1e1','0.0000001','1000001','1'*500,' secret ','=SUM(1)','']:
            with self.subTest(x=x),self.assertRaises(Invalid):cloud(points([(x,0,0)]),'m')
    def test_fixed_decimal(self):self.assertEqual(cloud(points([('-0.000001','+1.0','0')]),'m'),[(-1,1000000,0)])
    def test_csv_shape(self):
        for raw in [b'a,b,c\n',b'x,y,z\n0,0\n',b'x,y,z\n\n',b'x,y,z\n0,0,0,1\n',b'x,y,z\n"oops\n',b'\xff']:
            with self.subTest(raw=raw),self.assertRaises(Invalid):cloud(raw,'m')
    def test_cloud_limit(self):
        with self.assertRaises(Invalid):cloud(points([(0,0,0)]*1001),'m')
    def test_contract_types(self):
        for k,v in [('schema',True),('threshold_um',True),('threshold_um',1.0),('threshold_um',0),('minimum_recall_bp',10001),('frame','<secret>'),('reference_unit','cm'),('regions',{})]:
            c=cfg();c[k]=v
            with self.subTest(k=k,v=v),self.assertRaises(Invalid):contract(canonical(c))
    def test_contract_region(self):
        c=cfg();r={'id':'r','min_um':[0,0,0],'max_um':[1,1,1],'minimum_recall_bp':10000}
        for regs in [[r,r],[{**r,'max_um':[-1,0,0]}],[{**r,'min_um':[True,0,0]}],[r]*17]:
            c['regions']=regs
            with self.assertRaises(Invalid):contract(canonical(c))
    def test_bad_json(self):
        for raw in [b'{"secret":1,"secret":2}',b'{"x":NaN}',b'{"x":1e999}',b'{"x":"\\ud800"}',b'['*20000+b']'*20000,b'x'*(INPUT_LIMIT+1)]:
            with self.subTest(raw=raw[:50]),self.assertRaises(Invalid):contract(raw)
class Bundle(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.p=Path(self.tmp.name);self.c=self.p/'c.json';self.r=self.p/'r.csv';self.n=self.p/'n.csv';self.out=self.p/'bundle'
        self.c.write_bytes(canonical(cfg()));self.r.write_bytes(points([(0,0,0)]));self.n.write_bytes(points([(0,0,0)]))
    def tearDown(self):self.tmp.cleanup()
    def build(self):return build(self.c,self.r,self.n,self.out)
    def test_roundtrip_bytes(self):
        self.build();self.assertEqual(verify(self.out)['status'],'pass')
        for name in ['report.json','review.md','distances.csv']:self.assertNotIn(b'\r',read(self.out/name,OUTPUT_LIMIT))
    def test_tamper_each_output(self):
        self.build()
        for name in ['report.json','review.md','distances.csv','reference.csv','contract.json']:
            p=self.out/name;raw=p.read_bytes();p.write_bytes(raw+b' ')
            with self.assertRaises(Invalid):verify(self.out)
            p.write_bytes(raw)
    def test_no_overwrite(self):
        self.build();raw=(self.out/'report.json').read_bytes()
        with self.assertRaises(FileExistsError):self.build()
        self.assertEqual(raw,(self.out/'report.json').read_bytes())
    def test_invalid_no_directory(self):
        self.r.write_bytes(b'x,y,z\n')
        with self.assertRaises(Invalid):self.build()
        self.assertFalse(self.out.exists())
    def test_symlink(self):
        if os.name=='nt':self.skipTest('Symlink privilege is not assumed on Windows.')
        p=self.p/'link';p.symlink_to(self.r)
        with self.assertRaises(Invalid):read(p)
    def test_fifo(self):
        if not hasattr(os,'mkfifo'):self.skipTest('POSIX FIFO unavailable.')
        p=self.p/'fifo';os.mkfifo(p)
        cmd=subprocess.run([sys.executable,'-m','spatial_witness','inspect',str(p),'--unit','m'],cwd=ROOT,capture_output=True,timeout=3)
        self.assertEqual(cmd.returncode,2)
    def test_cli(self):
        def cli(*args):return subprocess.run([sys.executable,'-m','spatial_witness',*map(str,args)],cwd=ROOT,capture_output=True,text=True,timeout=10)
        self.assertEqual(cli('compare',self.c,self.r,self.n,'--out',self.out,'--require-pass').returncode,0)
        self.assertEqual(cli('verify',self.out,'--require-pass').returncode,0)
        self.assertFalse(json.loads(cli('inspect',self.r,'--unit','m').stdout)['quality_evaluated'])
        self.n.write_bytes(points([(10,0,0)]));self.assertEqual(cli('compare',self.c,self.r,self.n,'--out',self.p/'bad','--require-pass').returncode,3)
        self.assertEqual(cli('verify',self.p/'bad').returncode,0);self.assertEqual(cli('verify',self.p/'bad','--require-pass').returncode,3)
        self.n.write_bytes(b'x,y,z\nSECRET_SENTINEL,0,0\n');q=cli('compare',self.c,self.r,self.n,'--out',self.p/'invalid');self.assertEqual(q.returncode,2);self.assertNotIn('SECRET_SENTINEL',q.stdout+q.stderr);self.assertNotIn('Traceback',q.stderr)
    def test_example_missing_region(self):
        e=ROOT/'examples';r=evaluate(read(e/'contract.json'),read(e/'reference.csv'),read(e/'missing-doorway.csv'));self.assertEqual(r['precision']['ratio'],1);self.assertEqual(r['recall']['ratio'],.5);self.assertEqual(r['regions'][0]['status'],'fail')
if __name__=='__main__':unittest.main()
