"""PDB parsing and sequence-guided C-alpha structural comparison."""
import math
MAP=dict(zip('ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR'.split(),'ACDEFGHIKLMNPQRSTVWY'))
def parse_pdb(text):
    chains={};seen=set()
    for line in text.splitlines():
        if line.startswith('ENDMDL'):break
        if not line.startswith('ATOM  ') or len(line)<54 or line[16] not in ' A':continue
        chain=line[21];key=(chain,line[22:27],line[12:16])
        if key in seen:continue
        seen.add(key)
        if line[12:16].strip()=='CA' and line[17:20] in MAP:
            try:xyz=[float(line[30:38]),float(line[38:46]),float(line[46:54])]
            except ValueError:continue
            if not all(math.isfinite(x) for x in xyz):continue
            c=chains.setdefault(chain,{'sequence':'','coords':[],'residue_ids':[],'chain':chain})
            c['sequence']+=MAP[line[17:20]];c['coords'].append(xyz);c['residue_ids'].append(line[22:27].strip())
    return [c for c in chains.values() if 3<=len(c['sequence'])<=500]
def shape(coords):
    if len(coords)<3:return None
    import numpy as np
    x=np.asarray(coords,dtype=float);centered=x-x.mean(axis=0)
    rg=float(np.sqrt((centered**2).sum(axis=1).mean()))
    return {'radius_gyration':round(rg,3),'end_to_end':round(float(np.linalg.norm(x[0]-x[-1])),3),'observed_residues':len(x)}
def alignment(a,b):
    # Global minimum-edit alignment; diagonal preferred on ties.
    table=[list(range(len(b)+1))]
    for i,x in enumerate(a,1):
        row=[i]
        for j,y in enumerate(b,1):row.append(min(table[-1][j-1]+(x!=y),table[-1][j]+1,row[-1]+1))
        table.append(row)
    i,j=len(a),len(b);pairs=[];qa=[];ta=[]
    while i or j:
        if i and j and table[i][j]==table[i-1][j-1]+(a[i-1]!=b[j-1]):
            pairs.append((i-1,j-1));qa.append(a[i-1]);ta.append(b[j-1]);i-=1;j-=1
        elif i and table[i][j]==table[i-1][j]+1:qa.append(a[i-1]);ta.append('-');i-=1
        else:qa.append('-');ta.append(b[j-1]);j-=1
    return {'pairs':list(reversed(pairs)),'query':''.join(reversed(qa)),'target':''.join(reversed(ta))}
def fit(query,target,pairs):
    import numpy as np
    a=np.asarray([query[i] for i,j in pairs]);b=np.asarray([target[j] for i,j in pairs]);ac=a.mean(0);bc=b.mean(0)
    u,_,vt=np.linalg.svd((b-bc).T@(a-ac));correction=np.eye(3);correction[-1,-1]=np.linalg.det(u@vt);rot=u@correction@vt
    all_b=(np.asarray(target)-bc)@rot+ac
    rmsd=float(np.sqrt(((a-all_b[[j for i,j in pairs]])**2).sum()/len(pairs)))
    return round(rmsd,4),all_b.tolist()
