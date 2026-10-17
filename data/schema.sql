-- Schema only. Full data are in peptides.sqlite.
CREATE TABLE peptides(id TEXT PRIMARY KEY,name TEXT NOT NULL,sequence TEXT NOT NULL,organism TEXT,amidated INTEGER,pdb_id TEXT,source_url TEXT,smiles TEXT,coords TEXT NOT NULL,notes TEXT,metadata TEXT NOT NULL,properties TEXT,canonical INTEGER NOT NULL,length INTEGER NOT NULL);
CREATE TABLE provenance(source TEXT PRIMARY KEY,url TEXT,sha256 TEXT,downloaded_at TEXT,license TEXT);
CREATE TABLE structures(id TEXT PRIMARY KEY,peptide_id TEXT NOT NULL,pdb_id TEXT NOT NULL,chain TEXT NOT NULL,sequence TEXT NOT NULL,coords TEXT NOT NULL,residue_ids TEXT NOT NULL,file TEXT NOT NULL,properties TEXT NOT NULL,coverage REAL NOT NULL,source_url TEXT NOT NULL,title TEXT,citation TEXT);
CREATE INDEX idx_length ON peptides(length);
CREATE INDEX idx_name ON peptides(name);
CREATE INDEX idx_seq ON peptides(sequence);
CREATE INDEX idx_structure_peptide ON structures(peptide_id);
CREATE INDEX idx_structure_sequence ON structures(sequence);
