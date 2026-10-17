#!/usr/bin/env python3
"""Local-only research API; SQLite queries are parameterized."""
import argparse,csv,io,json,math,re,sqlite3,threading,functools
from pathlib import Path
from http.server import ThreadingHTTPServer,BaseHTTPRequestHandler
from urllib.parse import urlparse,parse_qs,unquote
from core import clean,describe,distance
from structure import parse_pdb,shape,alignment,fit
ROOT=Path(__file__).resolve().parent;DB=ROOT/'data/peptides.sqlite'
def connection():
    db=sqlite3.connect(DB);db.row_factory=sqlite3.Row;return db

def unpack(row,detail=False):
    r=dict(row)
    for key in ['properties','metadata','coords','residue_ids']:
        if key in r and isinstance(r[key],str):r[key]=json.loads(r[key])
    return r

def get_record(ident):
    with connection() as db:
        row=db.execute('SELECT * FROM peptides WHERE id=?',(ident,)).fetchone()
        if row is None:raise ValueError('Peptide not found')
        r=unpack(row,True);r['structures']=[unpack(x) for x in db.execute('SELECT * FROM structures WHERE peptide_id=? OR sequence=? ORDER BY coverage DESC,id LIMIT 30',(ident,r['sequence']))]
        return r

@functools.lru_cache(maxsize=1)
def corpus():
    with connection() as db:return [unpack(r) for r in db.execute('SELECT id,name,sequence,properties,source_url,metadata FROM peptides WHERE canonical=1')]
@functools.lru_cache(maxsize=1)
def structure_corpus():
    with connection() as db:return [unpack(r) for r in db.execute('SELECT s.*,p.name FROM structures s JOIN peptides p ON p.id=s.peptide_id ORDER BY s.id')]

def sequence_search(sequence,properties,mode='sequence',limit=15):
    if mode not in ['sequence','properties']:raise ValueError('Choose sequence or properties search')
    scored=[]
    for r in corpus():
        meta=r['metadata']
        has_pubmed=bool(re.search(r'\b\d{6,9}\b',str(meta.get('Pubmed_ID') or '')))
        title=str(meta.get('Title') or '').strip();reference=str(meta.get('Reference') or '').strip()
        if not has_pubmed and not (len(title)>12 and len(reference)>6 and title.lower() not in ['not included yet','not found']):continue
        if mode=='sequence':score=100*(1-distance(sequence,r['sequence'])/max(len(sequence),len(r['sequence'])))
        else:
            # Fixed, documented scales: length 20, mass 2000 Da, charge 3, GRAVY 1, moment 1.
            scales={'length':20,'mass':2000,'charge':3,'gravy':1,'hydrophobic_moment':1}
            d=sum(((properties[k]-r['properties'][k])/scale)**2 for k,scale in scales.items())/len(scales)
            score=100/(1+math.sqrt(d))
        scored.append((score,r))
    scored.sort(key=lambda x:(-x[0],x[1]['id']))
    # Show one representative per exact sequence; duplicate studies remain in the explorer.
    hits=[];seen=set()
    for score,r in scored:
        if r['sequence'] in seen:continue
        seen.add(r['sequence']);a=alignment(sequence,r['sequence']);paired=len(a['pairs']);identical=sum(sequence[i]==r['sequence'][j] for i,j in a['pairs'])
        meta=r['metadata'];hits.append(dict(id=r['id'],name=r['name'],sequence=r['sequence'],score=round(score,2),identity=round(100*identical/max(len(a['query']),1),2),coverage=round(100*paired/len(sequence),2),alignment={k:v for k,v in a.items() if k!='pairs'},properties=r['properties'],source_url=r['source_url'],pubmed_id=meta.get('Pubmed_ID'),reference=meta.get('Reference'),title=meta.get('Title'),activity=meta.get('Activity')))
        if len(hits)>=limit:break
    return hits

COMPARE_LOCK=threading.BoundedSemaphore(2)
def structure_search(data):
    if data.get('structure_id'):
        query=next((r for r in structure_corpus() if r['id']==data['structure_id']),None)
        if query is None:raise ValueError('Structure not found')
        query={**query,'label':query['id']}
    else:
        text=data.get('pdb_text','')
        if not isinstance(text,str) or len(text)>2_000_000:raise ValueError('PDB input must be at most 2 MB')
        chains=parse_pdb(text)
        if data.get('chain'):chains=[r for r in chains if r['chain']==data['chain']]
        if not chains:raise ValueError('No canonical protein C-alpha chain found in first PDB model')
        query=max(chains,key=lambda r:len(r['sequence']));query['label']='Uploaded PDB chain '+query['chain']
    if not 8<=len(query['sequence'])<=150:raise ValueError('Structure search supports 8–150 observed canonical residues')
    hits=[];seen=set()
    for target in structure_corpus():
        if target['id']==query.get('id'):continue
        a=alignment(query['sequence'],target['sequence']);pairs=a['pairs'];coverage=len(pairs)/max(len(query['sequence']),len(target['sequence']))
        if coverage<.8 or len(pairs)<8:continue
        rmsd,transformed=fit(query['coords'],target['coords'],pairs)
        key=(target['pdb_id'],target['sequence'])
        if key in seen:continue
        seen.add(key)
        hits.append(dict(id=target['id'],peptide_id=target['peptide_id'],name=target['name'],pdb_id=target['pdb_id'],chain=target['chain'],sequence=target['sequence'],rmsd=rmsd,coverage=round(coverage*100,1),aligned_residues=len(pairs),aligned_coords=transformed,properties=target['properties'],source_url=target['source_url'],citation=target['citation'],alignment={k:v for k,v in a.items() if k!='pairs'}))
    hits.sort(key=lambda r:(r['rmsd'],-r['coverage'],r['id']))
    return dict(query={'id':query.get('id'),'label':query['label'],'chain':query['chain'],'sequence':query['sequence'],'coords':query['coords'],'properties':shape(query['coords'])},matches=hits[:15],method='Sequence-guided global edit alignment + Kabsch C-alpha fit; minimum 80% coverage of longer chain and 8 aligned residues. Lower RMSD is better. Not fold search, TM-align, or evidence of shared activity.')

class Handler(BaseHTTPRequestHandler):
    def send(self,data,status=200,kind='application/json',filename=None):
        body=json.dumps(data,ensure_ascii=False,allow_nan=False).encode() if kind=='application/json' else data
        self.send_response(status);self.send_header('Content-Type',kind);self.send_header('Content-Length',str(len(body)));self.send_header('Cache-Control','no-cache');self.send_header('X-Content-Type-Options','nosniff')
        if filename:self.send_header('Content-Disposition','attachment; filename="'+filename+'"')
        self.end_headers();self.wfile.write(body)
    def do_GET(self):
        try:self.get()
        except (ValueError,TypeError) as e:self.send({'error':str(e)},400)
    def get(self):
        u=urlparse(self.path);q=parse_qs(u.query)
        if u.path=='/api/stats':
            with connection() as db:
                counts={k:db.execute(sql).fetchone()[0] for k,sql in {'peptides':'SELECT count(*) FROM peptides','sequences':'SELECT count(DISTINCT sequence) FROM peptides','canonical':'SELECT count(*) FROM peptides WHERE canonical=1','structures':'SELECT count(*) FROM structures','smiles':"SELECT count(*) FROM peptides WHERE smiles IS NOT NULL",'publications':"SELECT count(DISTINCT json_extract(metadata,'$.Pubmed_ID')) FROM peptides WHERE json_extract(metadata,'$.Pubmed_ID') IS NOT NULL"}.items()}
                pubmed=set()
                for r in db.execute("SELECT json_extract(metadata,'$.Pubmed_ID') FROM peptides"):
                    pubmed.update(re.findall(r'\b\d{6,9}\b',str(r[0] or '')))
                counts['publications']=len(pubmed)
                counts['sources']=[dict(r) for r in db.execute('SELECT * FROM provenance')]
            counts['data_bytes']=sum(p.stat().st_size for p in (ROOT/'data').rglob('*') if p.is_file());counts['limit_bytes']=1_000_000_000;return self.send(counts)
        if u.path.startswith('/api/peptides/'):
            try:return self.send(get_record(unquote(u.path.rsplit('/',1)[1])))
            except ValueError:return self.send({'error':'Peptide not found'},404)
        if u.path=='/api/peptides':
            limit=max(1,min(int(q.get('limit',['40'])[0]),100));offset=max(0,int(q.get('offset',['0'])[0]));term=q.get('q',[''])[0][:200]
            where=['(p.id LIKE ? OR p.name LIKE ? OR p.sequence LIKE ? OR p.organism LIKE ?)'];args=['%'+term+'%']*4
            if q.get('has_structure',['0'])[0]=='1':where.append('EXISTS(SELECT 1 FROM structures s WHERE s.peptide_id=p.id OR s.sequence=p.sequence)')
            if q.get('short',['0'])[0]=='1':where.append('p.length<=20')
            clause=' AND '.join(where)
            with connection() as db:
                total=db.execute('SELECT count(*) FROM peptides p WHERE '+clause,args).fetchone()[0]
                rows=[unpack(r) for r in db.execute('SELECT p.id,p.name,p.sequence,p.organism,p.length,p.properties,p.canonical,EXISTS(SELECT 1 FROM structures s WHERE s.peptide_id=p.id OR s.sequence=p.sequence) AS has_structure FROM peptides p WHERE '+clause+' ORDER BY p.name,p.id LIMIT ? OFFSET ?',args+[limit,offset])]
            return self.send({'items':rows,'total':total,'offset':offset,'limit':limit})
        if u.path=='/api/structures':return self.send([{k:r[k] for k in ['id','peptide_id','name','pdb_id','chain','sequence','properties']} for r in structure_corpus()])
        if u.path=='/api/export':return self.export(q)
        if u.path.startswith('/structures/'):
            file=unquote(u.path[len('/structures/'):]);path=ROOT/'data/structures'/file
            if path.suffix!='.pdb' or not path.resolve().is_relative_to((ROOT/'data/structures').resolve()) or not path.is_file():return self.send({'error':'Not found'},404)
            return self.send(path.read_bytes(),kind='text/plain',filename=path.name)
        path=ROOT/'static'/('index.html' if u.path=='/' else unquote(u.path.lstrip('/')))
        if not path.resolve().is_relative_to((ROOT/'static').resolve()) or not path.is_file():return self.send({'error':'Not found'},404)
        kind={'.html':'text/html; charset=utf-8','.css':'text/css','.js':'text/javascript','.svg':'image/svg+xml'}.get(path.suffix,'application/octet-stream');return self.send(path.read_bytes(),kind=kind)
    def export(self,q):
        fmt=q.get('format',['json'])[0];ident=q.get('id',[None])[0]
        with connection() as db:
            rows=[unpack(r) for r in db.execute('SELECT * FROM peptides'+(' WHERE id=?' if ident else ''),[ident] if ident else [])]
        if fmt=='json':return self.send(json.dumps(rows,ensure_ascii=False,indent=2).encode(),kind='application/json; charset=utf-8',filename='peptides.json')
        if fmt in ['fasta','fasts']:
            if fmt=='fasts':rows=[r for r in rows if r['smiles']]
            if not rows:return self.send({'error':'No qualifying records'},422)
            text=''.join('>'+r['id']+' '+(r['smiles'] if fmt=='fasts' else r['name'].replace('\n',' '))+'\n'+(str(r['metadata'].get('FASTS_sequence') or r['sequence']) if fmt=='fasts' else r['sequence'])+'\n' for r in rows)
            return self.send(text.encode(),kind='text/plain',filename='peptides.'+fmt)
        if fmt=='csv':
            out=io.StringIO();fields=['id','name','sequence','organism','amidated','pdb_id','source_url','smiles'];w=csv.DictWriter(out,fieldnames=fields);w.writeheader()
            for r in rows:
                # Prevent spreadsheet formula execution in exported untrusted source text.
                w.writerow({k:("'"+str(r[k]) if str(r[k] or '').startswith(('=','+','-','@')) else r[k]) for k in fields})
            return self.send(out.getvalue().encode(),kind='text/csv',filename='peptides.csv')
        return self.send({'error':'Unknown format'},400)
    def do_POST(self):
        if self.path not in ['/api/analyze','/api/structure-search']:return self.send({'error':'Not found'},404)
        try:
            length=int(self.headers.get('Content-Length',0))
            if length<1 or length>2_100_000:raise ValueError('Request must be 1 byte–2.1 MB')
            data=json.loads(self.rfile.read(length))
            if not isinstance(data,dict):raise ValueError('JSON object required')
            if self.path=='/api/structure-search':
                if not COMPARE_LOCK.acquire(blocking=False):return self.send({'error':'Structure search is busy; retry shortly'},429)
                try:return self.send(structure_search(data))
                finally:COMPARE_LOCK.release()
            s=clean(data.get('sequence',''));ph=float(data.get('ph',7))
            if not math.isfinite(ph) or not 0<=ph<=14:raise ValueError('pH must be between 0 and 14')
            amide=data.get('amidated',False)
            if not isinstance(amide,bool):raise ValueError('amidated must be a Boolean')
            properties=describe(s,amide,ph);mode=data.get('mode','sequence')
            # Neighbors use fixed pH 7 for comparability; query calculations use requested pH.
            hits=sequence_search(s,describe(s,amide,7),mode)
            return self.send(dict(sequence=s,properties=properties,matches=hits,mode=mode,model_status='Deterministic calculations; no trained activity model'))
        except (ValueError,TypeError,AttributeError,OverflowError) as e:return self.send({'error':str(e)},400)
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--port',type=int,default=8000);parser.add_argument('--host',default='127.0.0.1');a=parser.parse_args()
    print(f'AMP Atlas v0.2 at http://{a.host}:{a.port}',flush=True);ThreadingHTTPServer((a.host,a.port),Handler).serve_forever()
