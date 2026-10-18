"""Validate static data integrity, asset paths, quota, and no backend dependency."""
from pathlib import Path
import json,sqlite3,re
ROOT=Path(__file__).resolve().parents[1];docs=ROOT/'docs'
index=json.loads((docs/'data/search-index.json').read_text());manifest=json.loads((docs/'data/manifest.json').read_text());structures=json.loads((docs/'data/structures.json').read_text())
assert len(index)==manifest['peptides']==12784
assert len(structures)==manifest['structures']==433
assert len({r['id'] for r in index})==len(index)
shards={}
for r in index:
    shard=shards.setdefault(r['shard'],None)
    if shard is None:shard=shards[r['shard']]=json.loads((docs/'data/records'/r['shard']).read_text())
    assert shard[r['id']]['sequence']==r['sequence']
for r in structures:assert (docs/r['file']).is_file()
assert (docs/'downloads/peptides.sqlite').read_bytes()==(ROOT/'data/peptides.sqlite').read_bytes()
for ext in ['json','csv','fasta','fasts']:assert (docs/'downloads'/('peptides.'+ext)).is_file()
html=(docs/'index.html').read_text();app=(docs/'app.js').read_text()
assert not re.search(r'(?:src|href)="/(?!/)',html),'Root-relative asset path breaks GitHub project Pages'
assert "fetch('/api/" not in app and 'AMPClient.request' in app
assert (docs/'.nojekyll').exists()
for file in ['wheel-core.js','wheel-designer.js']:
    assert (docs/file).read_bytes()==(ROOT/'static'/file).read_bytes()
    assert f'src="{file}?v=0.5"' in html
assert "importScripts('browser-core.js?v=0.5','wheel-core.js?v=0.5')" in (docs/'browser-worker.js').read_text()
total=sum(p.stat().st_size for p in docs.rglob('*') if p.is_file());assert total<1_000_000_000
assert max(p.stat().st_size for p in ROOT.rglob('*') if p.is_file() and '.git' not in p.parts)<100_000_000
print(f'Static integrity passed: {len(index):,} records, {len(structures)} chains, {len(shards)} shards; {total/1e6:.1f} MB < 1 GB.')
