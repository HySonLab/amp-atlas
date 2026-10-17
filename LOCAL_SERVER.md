# AMP Atlas v0.2

Local antimicrobial-peptide research workspace inspired by Ma et al., **DRAMP 5.0**, Nucleic Acids Research (2026), https://doi.org/10.1093/nar/gkag989. Independent project; not an official DRAMP service. Built 2026-10-17.

## Start

Python **3.10+** and a WebGL-enabled browser are recommended. NumPy is required for 3D comparisons; openpyxl is needed only to rebuild the database. The existing database and viewer work offline after dependencies are installed. No API key or JavaScript build is required.

```sh
cd amp-atlas
python3 -m pip install -r requirements.txt
python3 server.py
```

Open http://localhost:8000. Stop with Ctrl+C. To use another port: `python3 server.py --port 8080`. The server binds to loopback by default.

## What changed

- Official DRAMP General collection: **12,784 records**, **11,718 distinct sequence strings** after trimming outer whitespace, **10,579 canonical sequence records** eligible for descriptor calculations, **5,183 imported source SMILES**.
- **433 experimental PDB chain structures**, linked by exact observed canonical sequence or a contiguous fragment covering at least 80% of a DRAMP sequence. Includes alternative structures/chains for the same sequence. Not 433 distinct peptides. 315 DRAMP records have a sequence-matched local structure in this build.
- Molecular viewer powered by **3Dmol.js**: local all-atom coordinates; cartoon, sticks, spheres, lines, solvent-accessible surface; coloring; residue labels; click-to-identify atoms; rotate/pan/zoom; spinning; PNG and PDB export. This is a browser molecular viewer, not the PyMOL engine.
- Calculation workbench: average molecular mass, charge at chosen pH, estimated pI, GRAVY, hydrophobic fraction, ideal-helix hydrophobic moment, aromaticity, reduced-cysteine extinction coefficient, amino-acid composition, charge curve, hydropathy plot and ideal helical wheel.
- Published neighbor search by sequence similarity or physicochemical descriptor proximity; global alignments and links to DRAMP/PubMed study records.
- 3D neighbor search for a local structure or uploaded PDB (stays local), with Cα RMSD, coverage, radius of gyration, end-to-end distance, and interactive backbone overlays.
- Paginated SQLite search, cached descriptors, source annotations, raw source files, reproducible import script, data budget checks and automated tests.

## Data files and provenance

`data/peptides.sqlite` is the complete runnable database. **Keep the whole `data/` directory**: PDB files live outside SQLite under `data/structures/` and are required for the molecular viewer. `data/schema.sql` contains schema only. Original DRAMP XLSX downloads are retained under `data/sources/`. `data/build-report.json` lists import counts, unavailable PDB IDs and unmatched structures. `data/peptides.csv` and `.fasta` are refreshed exports; complete JSON/FASTS can be downloaded through the UI/API without maintaining additional duplicate copies.

Sources:

- Official General AMPs: https://dramp.cpu-bioinfor.org/downloads/download.php?filename=download_data/DRAMP3.0_new/general_amps.xlsx
- Official SMILES supplement: https://dramp.cpu-bioinfor.org/downloads/download.php?filename=download_data/DRAMP3.0_new/general_smiles.xlsx
- PDB coordinate files: `https://files.rcsb.org/download/{PDB_ID}.pdb`
- DRAMP download index: https://dramp.cpu-bioinfor.org/downloads/

DRAMP annotations include publications, target/strain and activity text, hemolysis, cytotoxicity, structural descriptions and chemical modifications. These are imported source annotations, **not independently re-curated or normalized assay tables**. Some records lack study details, contain noncanonical notation, or have incomplete annotations. Database inclusion does not establish antibacterial efficacy or clinical suitability. Missing information remains missing.

We attempted all **591 PDB IDs linked to canonical records**. Two downloads were unavailable (5HAU, 7P5R); 211 downloaded entries did not yield a qualifying sequence-matched canonical chain. Those were not assigned speculative structures. First models are retained; chain files contain original atoms plus available connectivity/secondary-structure annotations. Unmatched files and redundant full NMR ensembles are discarded to keep storage small. PDB chains shorter than five or longer than 150 canonical observed residues are excluded from this build. PDB format, chain limits, missing residues and chemical modifications limit coverage.

A shared canonical sequence is not proof of shared covalent chemistry: cyclic, branched or modified DRAMP records may display a sequence-matched reference structure of a chemically different molecule. Read the source modifications and the specific PDB study before interpreting the picture. The interface explicitly labels this limitation. The canonical Cα parser excludes nonstandard residue names; the all-atom file may still contain modified residues and ligands.

The complete data directory occupies roughly **71 MB**, far below the **1,000,000,000-byte (1 GB)** ceiling. The exact current count appears in the app. No cloud database, remote inference or online asset CDN is required during ordinary use.

## Calculations and assumptions

All sequence descriptors approximate a canonical, linear peptide with free N terminus and a free or amidated C terminus. Database descriptors use source C-terminal amidation only when explicitly annotated as “Amidation”; unknown annotation assumes a free C terminus **for calculation only**. Console chemistry is controlled by the user.

- Average mass = residue masses + 18.0153 Da water; amidation subtracts 0.9848 Da. No disulfide oxidation, cyclization, branching, noncanonical stereochemistry, N-terminal modification, salt or adduct correction.
- Charge uses Henderson–Hasselbalch: pKa N=8.0, C=3.1, Lys=10.5, Arg=12.5, His=6.0, Asp=3.9, Glu=4.1, Cys=8.3, Tyr=10.1. Amidation removes the terminal carboxyl contribution. pI is a numerical root search over pH 0–14, clipped to those bounds when no root exists there; microscopic pKa and environments are not modeled.
- GRAVY = mean Kyte–Doolittle hydropathy. Hydrophobic fraction uses AVILMFWY; aromatic fraction uses FWY.
- Helical hydrophobic moment = `sqrt((Σhᵢ cos(i×100°))² + (Σhᵢ sin(i×100°))²) / N` on the **unnormalized Kyte–Doolittle scale**. The wheel displays the first 18 residues at 100° increments. Both assume an ideal α-helix; neither predicts structure or proves amphipathicity in an experimental environment.
- ε280 = 5500 × Trp + 1490 × Tyr in M⁻¹cm⁻¹, assuming reduced cysteines. Not a concentration estimate or a measured extinction coefficient.
- Shape descriptors use observed Cα atoms only: radius of gyration is the RMS distance from their centroid; end-to-end distance joins first and last observed residues. Not full-atom or mass-weighted Rg.

Values can differ materially from experimentally measured properties of chemically modified peptides. No trained AMP-activity, toxicity, MIC, secondary-structure or 3D prediction model is included.

## Similarity methods

### Sequence

Exact unit-cost global Levenshtein edit distance, computed with a bit-vector algorithm. Rank score = `100 × (1 − distance / max(query length, target length))`. This is **not BLAST** or BLOSUM-scored sequence alignment.

Only records with a source PubMed ID or meaningful study title/reference are eligible for neighbor results. Top hits are globally aligned using the same unit edit costs (diagonal preference on ties). Alignment identity = exact matching residue columns / all alignment columns, including gaps. Query coverage = paired non-gap columns / query length. Results show one representative per exact sequence; all source records remain accessible in the explorer. Same sequence can have different chemistry and assay results.

### Physicochemical

Fixed-scale RMS descriptor distance using length/20, mass/2000 Da, charge-at-pH-7/3, GRAVY/1, helical moment/1. Rank score = `100 / (1 + RMS difference)`. Scores are heuristic descriptor similarity, not probability, 3D similarity, or calibrated activity evidence. Eligible source-study records and query descriptors are compared at pH 7 even when the calculation report uses another pH.

### Structure

Supports 8–150 observed canonical residues. A global minimum-edit sequence alignment supplies paired Cα atoms; candidates require at least 8 pairs and ≥80% paired coverage of the **longer** observed chain. Kabsch least-squares fit uses a proper rotation (no reflection). RMSD is calculated over the paired atoms, while the overlay shows all observed Cα atoms. Hits are ranked by increasing RMSD and deduplicated by PDB ID plus observed sequence.

This is **sequence-guided structural comparison**, not a sequence-independent fold-search engine, TM-align, flexible alignment, or proof of function. There is no statistical significance estimate. Similar short helices can score well despite unrelated sequences; equally valid alignments can give different fits. NMR first models, solvent/membrane conditions and missing residues affect results. The query is omitted from hits for local structure search. For uploads, the user can specify a single-character chain; otherwise the longest eligible canonical first-model chain is chosen and reported.

## API

- `GET /api/stats` — counts, source URLs/hashes and storage use.
- `GET /api/peptides?q=magainin&has_structure=1&short=0&limit=25&offset=0` — paginated search; max page size 100.
- `GET /api/peptides/DRAMP02271` — full source annotation and linked structures.
- `GET /api/structures` — available structure index.
- `POST /api/analyze` — JSON `{"sequence":"GIGKFLHSAKKFGKAFVGEIMNS","ph":7,"amidated":false,"mode":"sequence"}`. Choose `properties` for descriptor search.
- `POST /api/structure-search` — JSON `{"structure_id":"2MAG_41"}`, or `{"pdb_text":"...","chain":"A"}`. PDB upload ≤2 MB. No upload is stored.
- `GET /api/export?format=json|csv|fasta|fasts` — optional `&id=DRAMP00068`.
- `GET /structures/2MAG_41.pdb` — local first-model all-atom chain.

FASTS uses `>DRAMP_ID SMILES` then the sequence from the official SMILES supplement when available. Original annotations remain in metadata. Imported source SMILES has not been independently validated with a chemistry toolkit. Missing SMILES is not generated. FASTA retains original DRAMP sequence notation, including noncanonical symbols; use canonical filters before machine-learning workflows. CSV neutralizes formula-leading fields for spreadsheet safety.

## Rebuild / expand locally

Back up `data/` before rebuilding. Close the running server during replacement, then restart to clear search caches.

```sh
python3 -m pip install -r requirements.txt
python3 scripts/build_database.py --max-structures 595 --workers 6
```

Existing DRAMP XLSX downloads are reused. To refresh them from the current upstream release, move the two source XLSX files out of `data/sources/` first. Do not silently overwrite a reviewed version. The script imports the General collection, merges SMILES by DRAMP ID and downloads **only** PDB IDs referenced by canonical source records. Download size caps, bounded concurrency and the 1 GB local-data budget are checked; a temporary SQLite database replaces the main database after successful import. The original source filename/hashes and run report make the shipped snapshot inspectable. Upstream changes can change record/structure counts on a future rebuild.

The shipped structure directory includes 433 selected first-model chains; the original source XLSX files are sufficient to reproduce import, but rebuilding structures requires access to RCSB. Actual upstream PDB bytes can evolve. Do not interpret missing files as negative biological evidence. A future true refresh should clear obsolete `data/structures/` in a staging copy before replacement; this script preserves old files rather than deleting user data.

## Tests

```sh
python3 -m unittest -v test_app.py
```

19 tests cover descriptor references, input validation, bit-vector versus independent DP distance, Kabsch rotation/translation invariance and reflection handling, alignment, source hashes, all structure-file mappings, SQLite integrity, budget, neighbor search, pagination, exports, invalid API requests and uploaded PDB comparison. Browser QA checked molecular representations, labels, sequence workbench and structural overlay. PNG generation is available in the viewer; automatic download-event capture was unavailable in the in-app preview, so verify the browser's download behavior when exporting a figure.

## PyMOL and larger deployment

Download any PDB from the viewer and open it directly in PyMOL for ray-traced publication figures. `examples/pymol_magainin.pml` uses an included local structure and supplies a ready-to-run styling example; start PyMOL from the project directory. PyMOL is optional and not bundled.

For much larger collections, move metadata to PostgreSQL, use a proper full-text and sequence index (MMseqs2/BLAST as appropriate), store coordinate objects separately, and normalize assay tables with strain, units, medium and conditions. For general 3D fold search, integrate TM-align/Foldseek with suitable peptide-specific validation. To add AI, use a separate versioned inference service with licensed training data, sequence-cluster splits, external validation and uncertainty reporting. This Python standard-library server is intended for local use; production hosting requires a production server, authentication for write operations, request controls and backups.

## Licenses and attribution

Application source: MIT (LICENSE). DRAMP data: **CC BY 4.0**, retain attribution to DRAMP and the reference paper above; imported annotations, combined SQLite and exports are adaptations/repackaging, not a new official release. PDB archival data: **CC0**, cite each original deposition/research article in scientific work: https://www.rcsb.org/pages/policies. 3Dmol.js is bundled locally with its own license and dependency notices in `static/vendor/`; https://3dmol.org/. The application's MIT license does not replace third-party data or library licenses.
