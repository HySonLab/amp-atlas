import hashlib,json,math,random,sqlite3,threading,unittest,urllib.request,urllib.error
import numpy as np
from core import describe,clean,distance,similarity,charge_at
from structure import alignment,fit,parse_pdb,shape
from server import ROOT,DB,Handler,ThreadingHTTPServer,get_record,sequence_search,structure_corpus
class Calculations(unittest.TestCase):
 def test_reference_mass(self):
  self.assertAlmostEqual(describe('G')['mass'],75.07,places=2)
  self.assertAlmostEqual(describe('G',True)['mass'],74.08,places=2)
 def test_ph_and_pi(self):
  for s in ['KKDE','GIGKFLHSAKKFGKAFVGEIMNS','DDD']:
   p=describe(s);self.assertLess(abs(charge_at(s,p['isoelectric_point'])),.01)
   self.assertGreater(charge_at(s,2),charge_at(s,12))
  self.assertGreater(describe('KKK')['charge'],2)
 def test_validation(self):
  for s in ['', 'AX', '>FASTA\nAAA','A'*501]:
   with self.assertRaises(ValueError):clean(s)
  self.assertEqual(clean(' gi g\n'),'GIG')
 def test_descriptors(self):
  p=describe('WYK');self.assertEqual(p['extinction_reduced'],6990);self.assertEqual(p['composition']['W'],1)
  self.assertEqual(len(p['charge_curve']),29)
 def test_bit_distance(self):
  rng=random.Random(42)
  for _ in range(200):
   a=''.join(rng.choices('ACDE',k=rng.randrange(1,80)));b=''.join(rng.choices('ACDE',k=rng.randrange(1,80)))
   self.assertEqual(round(100*(1-distance(a,b)/max(len(a),len(b))),1),similarity(a,b))
  self.assertEqual(distance('','AAA'),3)
 def test_alignment(self):
  a=alignment('ACDE','ACE');self.assertEqual(a['query'].replace('-',''),'ACDE');self.assertEqual(a['target'].replace('-',''),'ACE');self.assertEqual(len(a['pairs']),3)
 def test_fit_rotation_translation(self):
  rng=np.random.default_rng(9);a=rng.normal(size=(20,3));rot=np.array([[0,-1,0],[1,0,0],[0,0,1]])
  b=a@rot+np.array([40,-50,70]);rmsd,aligned=fit(a,b,[(i,i) for i in range(20)])
  self.assertLess(rmsd,1e-5);np.testing.assert_allclose(a,aligned,atol=1e-10)
 def test_fit_no_reflection(self):
  a=np.random.default_rng(10).normal(size=(15,3));b=a.copy();b[:,0]*=-1
  rmsd,_=fit(a,b,[(i,i) for i in range(15)]);self.assertGreater(rmsd,.1)
 def test_shape_invariant(self):
  a=np.array([[0,0,0],[1,0,0],[2,0,0]])
  self.assertEqual(shape(a.tolist()),shape((a+20).tolist()))
class Database(unittest.TestCase):
 def test_integrity_budget(self):
  with sqlite3.connect(DB) as db:
   self.assertEqual(db.execute('PRAGMA integrity_check').fetchone()[0],'ok');self.assertEqual(db.execute('SELECT count(*) FROM peptides').fetchone()[0],12784)
  self.assertLess(sum(p.stat().st_size for p in (ROOT/'data').rglob('*') if p.is_file()),1_000_000_000)
 def test_provenance_hash(self):
  with sqlite3.connect(DB) as db:
   for source,sha in db.execute('SELECT source,sha256 FROM provenance WHERE sha256 IS NOT NULL'):
    f=ROOT/'data/sources'/source.removeprefix('DRAMP ');self.assertEqual(hashlib.sha256(f.read_bytes()).hexdigest(),sha)
 def test_all_structure_files(self):
  rows=structure_corpus();self.assertGreater(len(rows),400)
  for r in rows:
   chain=parse_pdb((ROOT/'data'/r['file']).read_text())[0]
   self.assertEqual(chain['sequence'],r['sequence']);np.testing.assert_allclose(chain['coords'],r['coords'],atol=.001)
   with sqlite3.connect(DB) as db:seq=db.execute('SELECT sequence FROM peptides WHERE id=?',(r['peptide_id'],)).fetchone()[0]
   self.assertIn(r['sequence'],seq);self.assertGreaterEqual(len(r['sequence'])/len(seq),.8)
 def test_known_query(self):
  s='GIGKFLHSAKKFGKAFVGEIMNS';hits=sequence_search(s,describe(s));self.assertEqual(hits[0]['score'],100);self.assertEqual(hits[0]['sequence'],s)
  self.assertEqual(len({h['sequence'] for h in hits}),len(hits));self.assertTrue(hits[0]['pubmed_id'])
 def test_property_query(self):
  s='GIGKFLHSAKKFGKAFVGEIMNS';hits=sequence_search(s,describe(s),'properties');self.assertEqual(hits[0]['score'],100)
class HTTP(unittest.TestCase):
 @classmethod
 def setUpClass(cls):
  cls.server=ThreadingHTTPServer(('127.0.0.1',0),Handler);cls.thread=threading.Thread(target=cls.server.serve_forever,daemon=True);cls.thread.start();cls.base='http://127.0.0.1:'+str(cls.server.server_port)
 @classmethod
 def tearDownClass(cls):cls.server.shutdown();cls.server.server_close();cls.thread.join()
 def get(self,path):return urllib.request.urlopen(self.base+path)
 def post(self,path,data):return urllib.request.urlopen(urllib.request.Request(self.base+path,data=json.dumps(data).encode(),headers={'Content-Type':'application/json'}))
 def test_pagination_filters(self):
  a=json.load(self.get('/api/peptides?limit=2&has_structure=1'));b=json.load(self.get('/api/peptides?limit=2&has_structure=1&offset=2'))
  self.assertEqual(len(a['items']),2);self.assertNotEqual(a['items'][0]['id'],b['items'][0]['id']);self.assertTrue(all(x['has_structure'] for x in a['items']))
 def test_bad_input(self):
  for data in [{'sequence':'AAAA','ph':float('nan')},{'sequence':'AX'},{'sequence':'AA','amidated':'false'}]:
   with self.assertRaises(urllib.error.HTTPError) as e:self.post('/api/analyze',data)
   self.assertEqual(e.exception.code,400)
 def test_exports(self):
  with sqlite3.connect(DB) as db:ident=db.execute('SELECT id FROM peptides WHERE smiles IS NOT NULL LIMIT 1').fetchone()[0]
  for fmt in ['json','csv','fasta','fasts']:
   r=self.get('/api/export?format='+fmt+'&id='+ident);self.assertEqual(r.status,200);self.assertIn(b'DRAMP',r.read())
 def test_structure_and_upload(self):
  q=next(r for r in structure_corpus() if r['pdb_id']=='2MAG');d=json.load(self.post('/api/structure-search',{'structure_id':q['id']}));self.assertTrue(d['matches']);self.assertGreaterEqual(d['matches'][0]['coverage'],80)
  raw=(ROOT/'data'/q['file']).read_text();d=json.load(self.post('/api/structure-search',{'pdb_text':raw,'chain':'A'}));self.assertEqual(d['matches'][0]['rmsd'],0)
 def test_static_and_missing(self):
  self.assertEqual(self.get('/').status,200);self.assertEqual(self.get('/vendor/3Dmol-min.js').status,200)
  with self.assertRaises(urllib.error.HTTPError) as e:self.get('/api/peptides/UNKNOWN')
  self.assertEqual(e.exception.code,404)
if __name__=='__main__':unittest.main()
