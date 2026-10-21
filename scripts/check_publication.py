"""Validate the static scholarly landing page, its metadata and full-text links."""
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import urlparse
import json, xml.etree.ElementTree as ET
ROOT=Path(__file__).resolve().parents[1]
class Page(HTMLParser):
    def __init__(self):
        super().__init__();self.meta={};self.links=[];self.abstract=[];self.in_abstract=False;self.json_ld=[];self.in_json=False
    def handle_starttag(self,tag,attrs):
        a=dict(attrs)
        if tag=='meta' and a.get('name'):self.meta.setdefault(a['name'],[]).append(a.get('content',''))
        if tag=='a':self.links.append(a.get('href',''))
        if tag=='section' and a.get('class')=='abstract':
            assert 'hidden' not in a and 'display:none' not in a.get('style','').replace(' ','');self.in_abstract=True
        if tag=='script' and a.get('type')=='application/ld+json':self.in_json=True
    def handle_endtag(self,tag):
        if tag=='section':self.in_abstract=False
        if tag=='script':self.in_json=False
    def handle_data(self,data):
        if self.in_abstract:self.abstract.append(data)
        if self.in_json:self.json_ld.append(data)
def check():
    docs=ROOT/'docs';paper=docs/'papers/amp-atlas.html';p=Page();p.feed(paper.read_text())
    expected={'citation_title':'AMP Atlas: A Browser-Only Workbench for Antimicrobial Peptide Exploration and Interactive Helical-Wheel Design','citation_author':'Truong-Son Hy','citation_publication_date':'2026/10/18','citation_pdf_url':'https://hysonlab.github.io/amp-atlas/papers/amp-atlas-technical-report.pdf'}
    for key,value in expected.items():assert p.meta.get(key)==[value],key
    pdf_url=urlparse(expected['citation_pdf_url']);abstract_url=urlparse(p.meta['citation_abstract_html_url'][0]);assert pdf_url.netloc==abstract_url.netloc and Path(pdf_url.path).parent==Path(abstract_url.path).parent
    assert len(' '.join(p.abstract).split())>150,'Full abstract must be visible in source HTML'
    assert 'amp-atlas-technical-report.pdf' in p.links and 'amp-atlas.bib' in p.links
    assert 'noindex' not in p.meta.get('robots',[''])[0]
    for name in ['amp-atlas.html','amp-atlas.bib','amp-atlas-technical-report.pdf']:
        source=ROOT/'static/papers'/name;published=docs/'papers'/name;assert published.read_bytes()==source.read_bytes(),name
    pdf=(docs/'papers/amp-atlas-technical-report.pdf').read_bytes();assert pdf.startswith(b'%PDF-') and 0<len(pdf)<5_000_000
    assert pdf==(ROOT/'Technical_Report.pdf').read_bytes(),'Published PDF must preserve the author-supplied report'
    schema=json.loads(''.join(p.json_ld));assert ' '.join(schema['abstract'].split())==' '.join(' '.join(p.abstract).removeprefix('Abstract').split())
    assert schema['headline']==expected['citation_title'] and schema['author'][0]['name']==expected['citation_author'] and schema['encoding']['contentUrl']==expected['citation_pdf_url']
    home=Page();home.feed((docs/'index.html').read_text());assert 'papers/amp-atlas.html' in home.links and 'papers/amp-atlas-technical-report.pdf' in home.links
    assert not any(k.startswith('citation_') for k in home.meta),'Scholar tags belong to the paper page, not the app'
    sitemap=ET.parse(docs/'sitemap.xml');urls={n.text for n in sitemap.findall('.//{http://www.sitemaps.org/schemas/sitemap/0.9}loc')};assert p.meta['citation_abstract_html_url'][0] in urls
    print('Publication checks passed: static full abstract, citation metadata, linked PDF/BibTeX, homepage discovery, JSON-LD and sitemap.')
if __name__=='__main__':check()
