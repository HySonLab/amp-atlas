#!/usr/bin/env python3
"""Export a deterministic static GitHub Pages site from the local SQLite snapshot.
No third-party Python packages required. Does not crawl or modify the database.
"""
import csv,hashlib,io,json,shutil,sqlite3,re
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def build():
    docs=ROOT/'docs';docs.mkdir(exist_ok=True);(docs/'data/records').mkdir(parents=True,exist_ok=True);(docs/'downloads').mkdir(exist_ok=True)
    for source in (ROOT/'static').iterdir():
        if source.is_dir():shutil.copytree(source,docs/source.name,dirs_exist_ok=True)
        else:shutil.copy2(source,docs/source.name)
    def write(path,value):path.write_text(json.dumps(value,ensure_ascii=False,separators=(',',':'))+'\n')
    with sqlite3.connect(ROOT/'data/peptides.sqlite') as db:
        db.row_factory=sqlite3.Row
        def unpack(row):
            r=dict(row)
            for key in ['metadata','properties','coords','residue_ids']:
                if key in r and isinstance(r[key],str):r[key]=json.loads(r[key])
            return r
        records=[unpack(r) for r in db.execute('SELECT * FROM peptides ORDER BY name,id')]
        structures=[unpack(r) for r in db.execute('SELECT s.*,p.name FROM structures s JOIN peptides p ON p.id=s.peptide_id ORDER BY s.id')]
        structure_ids={r['peptide_id'] for r in structures};structure_seqs={r['sequence'] for r in structures}
        source=[dict(r) for r in db.execute('SELECT * FROM provenance')]
    shards={};index=[]
    for r in records:
        shard=(f"{int(r['id'][5:])//100:04d}" if re.fullmatch(r"DRAMP[0-9]+",r['id']) else "unassigned")+".json";shards.setdefault(shard,{})[r['id']]=r
        index.append({k:r[k] for k in ['id','name','sequence','organism','length','properties','canonical','source_url']}|{'has_structure':r['id'] in structure_ids or r['sequence'] in structure_seqs,'shard':shard,'study':{k:r['metadata'].get(k) for k in ['Pubmed_ID','Reference','Title','Activity']}})
    for name,values in shards.items():write(docs/'data/records'/name,values)
    write(docs/'data/search-index.json',index);write(docs/'data/structures.json',structures)
    shutil.copytree(ROOT/'data/structures',docs/'structures',dirs_exist_ok=True)
    for ext in ['csv','fasta']:shutil.copy2(ROOT/'data'/('peptides.'+ext),docs/'downloads'/('peptides.'+ext))
    write(docs/'downloads/peptides.json',records)
    (docs/'downloads/peptides.fasts').write_text(''.join('>'+r['id']+' '+r['smiles']+'\n'+str(r['metadata'].get('FASTS_sequence') or r['sequence'])+'\n' for r in records if r['smiles']))
    shutil.copy2(ROOT/'data/peptides.sqlite',docs/'downloads/peptides.sqlite')
    shutil.copy2(ROOT/'LICENSE',docs/'LICENSE');shutil.copy2(ROOT/'data/build-report.json',docs/'data/build-report.json')
    manifest=dict(version='0.3-browser',peptides=len(records),sequences=len({r['sequence'] for r in records}),canonical=sum(r['canonical'] for r in records),structures=len(structures),smiles=sum(bool(r['smiles']) for r in records),publications=len({p for r in records for p in re.findall(r'\b\d{6,9}\b',str(r['metadata'].get('Pubmed_ID') or ''))}),sources=source,limit_bytes=1_000_000_000,data_bytes=0,architecture='Static GitHub Pages; calculations in a browser Web Worker; uploaded sequences and PDB files stay on the device.')
    manifest['data_bytes']=sum(p.stat().st_size for folder in ['data','structures','downloads'] for p in (docs/folder).rglob('*') if p.is_file());write(docs/'data/manifest.json',manifest)
    app=(docs/'app.js').read_text();start=app.index('async function api(');end=app.index('\nfunction errorBox',start)
    app=app[:start]+"async function api(url,data){return window.AMPClient.request(url,data)}"+app[end:]
    app=app.replace("'/api/","'api/").replace('href="/api/','href="api/').replace("fetch('/'+s.file)","fetch(s.file)").replace("download:'/'+s.file","download:s.file")
    (docs/'app.js').write_text(app)
    html=(docs/'index.html').read_text().replace('href="/"','href="./"').replace('href="/style.css"','href="style.css"').replace('href="/api/','href="api/').replace('src="/','src="')
    html=html.replace('<script src="app.js">','<script src="browser-client.js"></script><script src="app.js">').replace('AMP Atlas / v0.2','AMP Atlas / v0.3').replace('Local research edition','HySonLab · Public research edition').replace('SQLite · Offline-ready','Browser-only · No account required').replace('● LOCAL INSTANCE','● BROWSER COMPUTE').replace('Local research edition','Browser research edition').replace('Of 1 GB data budget','Static dataset · Under 1 GB')
    html=html.replace('SQLite search','Browser search').replace('Local build manifest','Dataset & browser manifest').replace('local all-atom PDB','static all-atom PDB').replace('No cloud database,','No cloud database,')
    html=html.replace('<div class="eyebrow">SEQUENCE → STRUCTURE → EVIDENCE</div>','<div class="eyebrow">HYSONLAB · SEQUENCE → STRUCTURE → EVIDENCE</div>')
    html=html.replace('<footer>AMP ATLAS','<div class="privacy-note">Sequence and uploaded PDB calculations run on your device. GitHub serves public data and may log ordinary page/file requests. External study links open their respective providers.</div><footer>AMP ATLAS')
    html=html.replace('<h2>DRAMP collection</h2>','<h2>DRAMP collection</h2><p><a href="downloads/peptides.sqlite" download>↓ Complete SQLite database</a> · <a href="https://github.com/HySonLab/amp-atlas" target="_blank" rel="noreferrer">Source code ↗</a></p>')
    (docs/'index.html').write_text(html);(docs/'.nojekyll').touch()
    with (docs/'style.css').open('a') as f:f.write('\n.privacy-note{margin-top:36px;font-size:10px;line-height:1.8;color:var(--muted);border-top:1px solid var(--line);padding-top:18px}\n')
    total=sum(p.stat().st_size for p in docs.rglob('*') if p.is_file());assert total<1_000_000_000
    assert max(p.stat().st_size for p in docs.rglob('*') if p.is_file())<100_000_000
    print(json.dumps({'records':len(records),'structures':len(structures),'detail_shards':len(shards),'site_bytes':total,'index_bytes':(docs/'data/search-index.json').stat().st_size},indent=2))
if __name__=='__main__':build()
