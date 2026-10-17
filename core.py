"""Deterministic sequence descriptors and exact global edit distance.
Canonical residue chemistry only; these functions do not predict AMP activity.
"""
import math,re
AA='ACDEFGHIKLMNPQRSTVWY'
KD=dict(zip(AA,[1.8,2.5,-3.5,-3.5,2.8,-.4,-3.2,4.5,-3.9,3.8,1.9,-3.5,-1.6,-3.5,-4.5,-.8,-.7,4.2,-.9,-1.3]))
MASS=dict(zip(AA,[71.0788,103.1388,115.0886,129.1155,147.1766,57.0519,137.1411,113.1594,128.1741,113.1594,131.1926,114.1038,97.1167,128.1307,156.1875,87.0782,101.1051,99.1326,186.2132,163.1760]))
def clean(sequence):
    if not isinstance(sequence,str):raise ValueError('Sequence must be text')
    s=re.sub(r'\s','',sequence).upper()
    if not s or len(s)>500 or any(a not in AA for a in s):
        raise ValueError('Enter 1–500 canonical amino acids (single-letter codes). FASTA headers and nonstandard residues are not supported in the console.')
    return s

def charge_at(s,ph=7,amidated=False):
    value=1/(1+10**(ph-8.0))-(0 if amidated else 1/(1+10**(3.1-ph)))
    for a,p in [('K',10.5),('R',12.5),('H',6.0)]:value+=s.count(a)/(1+10**(ph-p))
    for a,p in [('D',3.9),('E',4.1),('C',8.3),('Y',10.1)]:value-=s.count(a)/(1+10**(p-ph))
    return value

def describe(sequence,amidated=False,ph=7):
    s=clean(sequence);n=len(s)
    if not math.isfinite(ph) or not 0<=ph<=14:raise ValueError('pH must be between 0 and 14')
    lo,hi=0.,14.
    for _ in range(45):
        mid=(lo+hi)/2
        if charge_at(s,mid,amidated)>0:lo=mid
        else:hi=mid
    x=sum(KD[a]*math.cos(math.radians(i*100)) for i,a in enumerate(s))
    y=sum(KD[a]*math.sin(math.radians(i*100)) for i,a in enumerate(s))
    return dict(length=n,mass=round(sum(MASS[a] for a in s)+18.0153-(.9848 if amidated else 0),2),charge=round(charge_at(s,ph,amidated),3),ph=ph,isoelectric_point=round((lo+hi)/2,3),gravy=round(sum(KD[a] for a in s)/n,3),hydrophobic=round(100*sum(a in 'AVILMFWY' for a in s)/n,1),hydrophobic_moment=round(math.hypot(x,y)/n,3),aromaticity=round(100*sum(s.count(a) for a in 'FWY')/n,2),extinction_reduced=s.count('W')*5500+s.count('Y')*1490,composition={a:s.count(a) for a in AA},profile=[KD[a] for a in s],charge_curve=[dict(ph=v/2,charge=round(charge_at(s,v/2,amidated),3)) for v in range(29)])

def distance(a,b):
    """Myers bit-vector Levenshtein distance, unit insertion/deletion/substitution."""
    if not a:return len(b)
    masks={}
    for i,c in enumerate(a):masks[c]=masks.get(c,0)|(1<<i)
    positive=~0;negative=0;score=len(a);last=1<<(len(a)-1)
    for c in b:
        eq=masks.get(c,0);x=eq|negative;d=(((eq&positive)+positive)^positive)|eq
        ph=negative|~(d|positive);mh=positive&d
        if ph&last:score+=1
        elif mh&last:score-=1
        ph=(ph<<1)|1;mh<<=1;positive=mh|~(x|ph);negative=ph&x
    return score

def similarity(a,b):
    """Slow independent DP reference, retained to test the bit-vector method."""
    if not a and not b:return 100.
    prev=list(range(len(b)+1))
    for i,x in enumerate(a,1):
        row=[i]
        for j,y in enumerate(b,1):row.append(min(row[-1]+1,prev[j]+1,prev[j-1]+(x!=y)))
        prev=row
    return round(100*(1-prev[-1]/max(len(a),len(b))),1)
