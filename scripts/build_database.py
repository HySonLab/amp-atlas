#!/usr/bin/env python3
"""Rebuild local data from official DRAMP downloads and bounded PDB downloads.
Requires openpyxl and numpy. Resume-safe downloads; atomic database replacement.
"""
import argparse,concurrent.futures,hashlib,json,re,sqlite3,sys,time,urllib.request,threading
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from core import describe,AA
from structure import parse_pdb,shape
LIMIT=1_000_000_000
STORAGE_LOCK=threading.Lock()

def bounded_download(url,path,max_bytes=12_000_000):
    if path.exists():return path.read_bytes()
    with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'AMPAtlas-research/0.2'}),timeout=35) as response:
        if int(response.headers.get('Content-Length',0))>max_bytes:raise ValueError('Remote file exceeds download cap')
        data=response.read(max_bytes+1)
    if len(data)>max_bytes:raise ValueError('Remote file exceeds download cap')
    with STORAGE_LOCK:
        if sum(p.stat().st_size for p in (ROOT/'data').rglob('*') if p.is_file())+len(data)>LIMIT-30_000_000:raise ValueError('1 GB storage budget reached')
        path.write_bytes(data)
    return data

def build(max_structures=300,workers=6):
    import openpyxl
    data=ROOT/'data';sources=data/'sources';sources.mkdir(exist_ok=True);(data/'structures').mkdir(exist_ok=True)
    urls={}
    for key,file in [('general_amps','dramp-general.xlsx'),('general_smiles','dramp-general-smiles.xlsx')]:
        url='https://dramp.cpu-bioinfor.org/downloads/download.php?filename=download_data/DRAMP3.0_new/'+key+'.xlsx'
        bounded_download(url,sources/file);urls[file]=url
    values=iter(openpyxl.load_workbook(sources/'dramp-general.xlsx',read_only=True).active.values);headers=next(values);rows=[dict(zip(headers,r)) for r in values]
    smile_values=iter(openpyxl.load_workbook(sources/'dramp-general-smiles.xlsx',read_only=True).active.values);sh=next(smile_values);smile_rows={r[0]:dict(zip(sh,r)) for r in smile_values};smiles={k:r.get('SMILES') for k,r in smile_rows.items()}
    temporary=data/'build.sqlite';temporary.unlink(missing_ok=True);db=sqlite3.connect(temporary)
    db.executescript('''CREATE TABLE peptides(id TEXT PRIMARY KEY,name TEXT NOT NULL,sequence TEXT NOT NULL,organism TEXT,amidated INTEGER,pdb_id TEXT,source_url TEXT,smiles TEXT,coords TEXT NOT NULL,notes TEXT,metadata TEXT NOT NULL,properties TEXT,canonical INTEGER NOT NULL,length INTEGER NOT NULL);
    CREATE INDEX idx_name ON peptides(name); CREATE INDEX idx_length ON peptides(length);CREATE INDEX idx_seq ON peptides(sequence);
    CREATE TABLE structures(id TEXT PRIMARY KEY,peptide_id TEXT NOT NULL,pdb_id TEXT NOT NULL,chain TEXT NOT NULL,sequence TEXT NOT NULL,coords TEXT NOT NULL,residue_ids TEXT NOT NULL,file TEXT NOT NULL,properties TEXT NOT NULL,coverage REAL NOT NULL,source_url TEXT NOT NULL,title TEXT,citation TEXT);
    CREATE INDEX idx_structure_peptide ON structures(peptide_id); CREATE INDEX idx_structure_sequence ON structures(sequence);
    CREATE TABLE provenance(source TEXT PRIMARY KEY,url TEXT,sha256 TEXT,downloaded_at TEXT,license TEXT);
    ''')
    referenced={};canonical_rows={};stats={'source_records':len(rows),'imported':0,'canonical':0,'smiles':0,'structures':0,'download_failures':[],'unmatched_pdb':[]}
    for r in rows:
        seq=str(r.get('Sequence') or '').strip();ident=str(r['DRAMP_ID']);r['FASTS_sequence']=smile_rows.get(ident,{}).get('Sequence');valid=bool(re.fullmatch('[ACDEFGHIKLMNPQRSTVWY]{1,500}',seq));ct=str(r.get('C-terminal_Modification') or '').strip().lower();amide=1 if ct=='amidation' else 0 if ct=='free' else None
        props=describe(seq,bool(amide)) if valid else None
        # Descriptor chemistry is explicitly a sequence-only approximation.
        if props:
            props.pop('charge_curve',None);props.pop('profile',None);props.pop('composition',None)
        pids=sorted(set(x.upper() for x in re.findall(r'\b[0-9][A-Za-z0-9]{3}\b',str(r.get('PDB_ID') or ''))))
        sm=r.get('SMILES') or smiles.get(ident);sm=sm.strip() if isinstance(sm,str) and sm.strip() and sm.strip().lower() not in ['none','not found'] else None
        db.execute('INSERT INTO peptides VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?)',(ident,str(r.get('Name') or ident),seq,str(r.get('Source') or 'Not annotated'),amide,','.join(pids),f'https://dramp.cpu-bioinfor.org/browse/All_Information.php?id={ident}',sm,'[]','Imported DRAMP annotation. Sequence descriptors approximate canonical chemistry; source modifications are not chemically modeled.',json.dumps(r,ensure_ascii=False,default=str),json.dumps(props) if props else None,int(valid),len(seq)))
        stats['imported']+=1;stats['canonical']+=valid;stats['smiles']+=bool(sm)
        if valid:
            canonical_rows[ident]=r
            for pid in pids:referenced.setdefault(pid,[]).append((ident,seq))
    # Prioritize the peptides in the original demonstration, then diverse referenced structures.
    pids=sorted(referenced,key=lambda p:(p not in ['2MAG','2MLT','1G89'],p))[:max_structures]
    cache=sources/'pdb-cache';cache.mkdir(exist_ok=True)
    def fetch(pid):
        try:
            f=cache/(pid+'.pdb');raw=bounded_download('https://files.rcsb.org/download/'+pid+'.pdb',f)
            text=raw.decode('utf-8',errors='replace');chains=parse_pdb(text);lines=text.splitlines();title=' '.join(x[10:].strip() for x in lines if x.startswith('TITLE '));citation='\n'.join(x for x in lines if x.startswith('JRNL'))
            first=[]
            for x in lines:
                if x.startswith('ENDMDL'):break
                first.append(x)
            matched=[]
            for c in chains:
                if not 5<=len(c['sequence'])<=150:continue
                matches=[(i,s) for i,s in referenced[pid] if c['sequence']==s or (c['sequence'] in s and len(c['sequence'])/len(s)>=.8)]
                if not matches:continue
                ident,full=sorted(matches,key=lambda v:(v[1]!=c['sequence'],v[0]))[0]
                atoms=[x for x in first if x.startswith(('ATOM  ','HETATM')) and len(x)>21 and x[21]==c['chain']]
                serials={x[6:11].strip() for x in atoms}
                bonds=[x for x in first if x.startswith('CONECT') and all(x[i:i+5].strip() in serials for i in range(6,len(x.rstrip()),5) if x[i:i+5].strip())]
                annotation=[x for x in first if x.startswith(('HEADER','TITLE ','HELIX ','SHEET ','SSBOND'))]
                payload='\n'.join(annotation+atoms+bonds+['END'])+'\n'
                file=f"structures/{pid}_{ord(c['chain']):02x}.pdb";(data/file).write_text(payload)
                c.update(id=f"{pid}_{ord(c['chain']):02x}",peptide_id=ident,pdb_id=pid,file=file,coverage=len(c['sequence'])/len(full),title=title,citation=citation)
                matched.append(c)
            # Keep no redundant raw ensembles; all downloaded atoms used by the viewer are retained above.
            f.unlink(missing_ok=True)
            return pid,matched,None
        except Exception as e:return pid,[],str(e)
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        for n,(pid,items,error) in enumerate(pool.map(fetch,pids),1):
            if error:stats['download_failures'].append({'pdb_id':pid,'error':error})
            elif not items:stats['unmatched_pdb'].append(pid)
            for c in items:
                db.execute('INSERT OR IGNORE INTO structures VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?)',(c['id'],c['peptide_id'],pid,c['chain'],c['sequence'],json.dumps(c['coords']),json.dumps(c['residue_ids']),c['file'],json.dumps(shape(c['coords'])),c['coverage'],f'https://www.rcsb.org/structure/{pid}',c['title'],c['citation']))
                stats['structures']+=1
            if n%20==0:print(f'{n}/{len(pids)} PDB files processed; {stats["structures"]} chain structures',flush=True)
    for file,url in urls.items():db.execute('INSERT INTO provenance VALUES(?,?,?,?,?)',('DRAMP '+file,url,hashlib.sha256((sources/file).read_bytes()).hexdigest(),time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'CC BY 4.0; cite DRAMP 5.0'))
    db.execute('INSERT INTO provenance VALUES(?,?,?,?,?)',('RCSB PDB','https://files.rcsb.org/download/',None,time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'PDB archive CC0; cite individual depositions'))
    db.commit();db.execute('VACUUM');db.close();temporary.replace(data/'peptides.sqlite')
    stats['unique_sequences']=len(set(str(r.get('Sequence') or '').strip() for r in rows));stats['pdb_files_attempted']=len(pids);stats['data_bytes']=sum(p.stat().st_size for p in data.rglob('*') if p.is_file());stats['limit_bytes']=LIMIT
    (data/'build-report.json').write_text(json.dumps(stats,indent=2));print(json.dumps({k:v for k,v in stats.items() if not isinstance(v,list)},indent=2),flush=True)
    if stats['data_bytes']>=LIMIT:raise RuntimeError('Database budget exceeded')
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--max-structures',type=int,default=300);parser.add_argument('--workers',type=int,default=6);args=parser.parse_args();build(args.max_structures,args.workers)
