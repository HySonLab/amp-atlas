# AMP Atlas · HySonLab

A browser-only antimicrobial peptide research workspace: published peptide records, interactive molecular visualization, transparent calculations, and sequence/descriptor/3D neighbor searches.

**Website:** https://hysonlab.github.io/amp-atlas/  
**Source code:** https://github.com/HySonLab/amp-atlas

**Technical report:** https://arxiv.org/pdf/2610.29096

No rented server, API key, sign-in, remote inference or Python runtime is needed for visitors. The static website and its downloadable database are included in this repository. GitHub Pages serves public files; a Web Worker performs calculations on the visitor's device. Uploaded sequences and PDB files are never submitted to a server. GitHub may record ordinary page/data-file requests. External research links open the cited providers.

## Technical Report

The technical report describing the architecture, algorithms, implementation,
and verification of AMP Atlas is publicly available:

**Truong-Son Hy (2026).**  
*AMP Atlas: A Browser-Only Workbench for Antimicrobial Peptide Exploration and Interactive Helical-Wheel Design.*  
Technical Report, University of Alabama at Birmingham, October 2026.

📄 **[Read the Technical Report (PDF)](./Technical_Report.pdf)**

The report presents:
- Integration of antimicrobial peptide records, chemical annotations, and PDB structures.
- Interactive sequence analysis and molecular visualization.
- An inverse helical-wheel designer with constraint-based peptide exploration.
- Sequence-guided structural comparison and physicochemical calculations.
- Reproducible software verification and performance measurements.

### Website hosting and scholarly discovery

- [Read the abstract and citation](https://hysonlab.github.io/amp-atlas/papers/amp-atlas.html)
- [View the report on the website (PDF)](https://hysonlab.github.io/amp-atlas/papers/amp-atlas-technical-report.pdf)

The author-supplied PDF is preserved in `static/papers/` and copied unchanged into `docs/papers/` by the Pages builder. The standalone HTML landing page exposes the full abstract without JavaScript and includes Highwire citation metadata, a same-directory absolute PDF URL, JSON-LD and a downloadable BibTeX entry. Plain homepage links and a sitemap support discovery. Publication checks are part of the existing static-integrity workflow.

These changes support crawler discovery; they do not submit a paper to Google Scholar or guarantee inclusion. See Google's inclusion guidelines: https://scholar.google.com/intl/en/scholar/inclusion.html. An institutional publications page or repository can also link to the report. If an arXiv version is published later, add its actual identifier and link; no arXiv identifier or DOI is claimed here.

## Features

- 12,784 DRAMP source records; 11,718 distinct sequence strings; 10,579 canonical records eligible for sequence calculations; 5,183 source SMILES annotations.
- 433 sequence-matched PDB chain structures, including multiple structures/chains for the same sequence; 315 DRAMP records have a local sequence-matched structure.
- 3Dmol.js all-atom molecular viewer: cartoon, sticks, spheres, lines, solvent-accessible surfaces, atom/residue identification, labels, rotation and PNG/PDB export.
- Interactive inverse helical-wheel designer: choose 1–72 residues, fill numbered positions, read the linear N→C sequence, search partial positional constraints, and export FASTA/SVG.
- Mass, charge versus pH, estimated pI, GRAVY, hydrophobic fraction, ideal-helix hydrophobic moment, aromaticity, reduced-cysteine extinction coefficient, amino-acid composition and plots.
- Published-neighbor searches by global edit similarity or scaled physicochemical descriptors, with alignments and study links.
- 3D comparisons using sequence-guided Cα alignments, proper-rotation least-squares fits, RMSD, coverage, shape descriptors and interactive overlays. Upload a PDB file or select a local structure.
- Bulk and per-record JSON, CSV, FASTA and FASTS downloads; a complete SQLite database.

## Run the browser-only site locally

Serve the existing `docs/` directory with any static HTTP server:

```sh
python3 -m http.server 8000 --directory docs
```

Open http://localhost:8000. Python is only a convenient file server here; the browser does all calculations. Opening `index.html` through a `file://` URL is unsupported because workers and data fetches require HTTP(S).

A current Chrome, Edge, Firefox or Safari with Web Workers and WebGL is recommended. Initial loading downloads a roughly 9 MB uncompressed searchable index; gzip/Brotli-capable hosting reduces transfer. Full annotations load on demand in 145 record shards. PDB files load only when viewed; full database downloads occur only on request. The rendered site, data and downloads total about **139 MB**, below the 1 GB budget. This is an online static site, not a guaranteed fully offline PWA.

## Design a peptide on the helical wheel

Open **Helical-wheel designer** in the sidebar:

1. Set **Number of amino acids** (1–72) and press **Apply length**.
2. Click a numbered wheel position, then select a canonical residue from the palette. You can also focus a wheel position and type letters. Entry advances in sequence order; arrow keys move between positions and Delete clears a position.
3. The **Linear sequence · N → C** field updates immediately in **position-number order**, not clockwise angular order. Empty positions appear as `?`; they are wildcards, not amino acids.
4. With every position filled, copy the sequence, download FASTA, or choose **Analyze & find related peptides** to use the existing calculation/similarity workbench. The workbench's pH, terminal chemistry and similarity mode are preserved.
5. With a partially or fully filled wheel, **Find exact positional matches** searches source-study-linked canonical DRAMP records of exactly the selected length. Every specified residue must match its numbered position; blank positions are unconstrained. The first 100 matching source records are shown, and links open the existing explorer. At least one residue is required. No match is a valid result; it does not imply the designed sequence lacks activity or exists nowhere else.

**Import a linear sequence**, **Use workbench sequence**, or the workbench's **Edit this sequence on the wheel** button provides a round trip for sequences up to 72 residues. Length changes preserve the common prefix. **Undo** restores previous fills, clears, imports and length changes (up to 40 edits). Download **Wheel SVG** for a vector illustration, including partially filled designs.

Positions follow the ideal 100° α-helical projection. After position 18, extra display rings prevent coincident slots; these rings are a layout aid, not covalent cycles, predicted 3D coordinates or a change in sequence order. The designer can construct any canonical helical-peptide design, while database searches are limited to the included DRAMP collection. There is no helicity or activity prediction. The pure model lives in `static/wheel-core.js` and the isolated UI in `static/wheel-designer.js`; existing scientific calculations and PDB data remain unchanged.

## Free GitHub Pages deployment

For the HySonLab organization:

1. Create a **public** repository named `amp-atlas` and push the contents of this project to its `main` branch. Keep the whole project, including `docs/` and the data.
2. In **Settings → Pages**, set Source to **Deploy from a branch**, branch **main**, folder **/docs**, and Save.
3. Wait for the Pages deployment to complete, then visit https://hysonlab.github.io/amp-atlas/.

`.nojekyll` serves the static data without a Jekyll build. All application, worker, coordinate and download paths are relative to the repository base URL. Do not enable Git LFS for data needed directly by Pages; the committed files are real files, each below 100 MB. GitHub Free supports public repository Pages subject to GitHub's usage limits and policies; this project has no server rental or external paid services.

The verification workflow runs Python scientific/database tests plus Node browser-worker tests. The Pages branch deployment is independent of those tests; review failed checks before accepting later scientific changes.

## Architecture and source layout

| Path | Purpose |
| --- | --- |
| `docs/` | Complete deployable static site, including browser data and downloads |
| `static/app.js`, `static/index.html`, `static/style.css` | Shared application UI sources |
| `static/browser-core.js` | Browser descriptor calculations, alignment, PDB parser and quaternion superposition |
| `static/browser-worker.js` | Browser dataset access, searches, local upload handling and record exports |
| `static/browser-client.js` | Worker transport and browser download handling |
| `scripts/build_pages.py` | Rebuild the static site from the supplied SQLite snapshot; standard-library Python only |
| `scripts/check_static.py` | Check the generated site, data consistency, asset paths and size |
| `data/peptides.sqlite` | Complete source database |
| `data/structures/` | All-atom first-model PDB chain files |
| `data/sources/` | Original official DRAMP XLSX source files |
| `scripts/build_database.py` | Rebuild/crawl official DRAMP and referenced RCSB PDB data |
| `core.py`, `structure.py`, `server.py` | Original Python reference implementation and optional local API |
| `LOCAL_SERVER.md` | Detailed v0.2 data provenance, assumptions and original local API documentation |
| `tests/` | Browser implementation tests and fixed NumPy reference fixtures |
| `examples/pymol_magainin.pml` | Optional styling example for an existing PyMOL installation |

Rebuild after changing the shared UI, browser engine, or database:

```sh
python3 scripts/build_pages.py
python3 scripts/check_static.py
pnpm install --frozen-lockfile --ignore-scripts
pnpm test
```

To update the database itself, read `LOCAL_SERVER.md` first, install `requirements.txt`, rebuild the database, refresh CSV/FASTA exports as needed, and then rebuild Pages. The browser export builder copies the authoritative SQLite snapshot and reads its records; it never crawls at deployment time. This repository intentionally includes both source data and generated deployable files for reproducibility and zero-build hosting. All hosted data are public.

## Scientific methods and limits

Sequence descriptors use canonical, linear residue chemistry, free N terminus, and optionally amidated C terminus. Unknown source termini use free termini **only as a calculation assumption**. Cyclic/branched chemistry, disulfide oxidation, nonstandard residues, stereochemistry, N-terminal modifications, salt/adduct state and environment-specific pKa are not chemically modeled. Imported source SMILES has not been independently validated. Numerical descriptors are calculations, not experimental activity measurements.

Mass is average residue mass plus water, optionally adjusted for C-amidation. Charge/pI use the Henderson–Hasselbalch parameters in `core.py`; pI is searched within pH 0–14. GRAVY and moment use the unnormalized Kyte–Doolittle scale. The moment and first-18-residue wheel assume an ideal 100° α-helix; they do not predict secondary structure. ε280 assumes reduced cysteines.

Sequence score is `100 × (1 − Levenshtein distance / longer length)`, with unit edit costs. Global alignment identity includes gap columns. Descriptor score is `100 / (1 + RMS scaled difference)` for length/20, mass/2000, pH-7 charge/3, GRAVY/1 and moment/1. Neighbor search requires a source PubMed ID or meaningful source study title/reference and returns one representative per exact sequence. Same sequence can describe different chemistry.

The browser fits sequence-aligned observed Cα atoms using **Horn's quaternion least-squares method**, numerically equivalent to the Python Kabsch proper-rotation fit. At least 8 pairs and ≥80% coverage of the longer observed chain are required; query lengths are 8–150 canonical residues. Tests compare the browser's fitted coordinates and RMSD against NumPy references. Shape values use Cα atoms, not mass-weighted all-atom coordinates. Structure search is sequence-guided comparison, not sequence-independent fold search, flexible alignment, TM-align or functional validation. Similar short helices may score well despite unrelated sequences. First-model conformations and experimental environments matter.

PDB chains are linked by exact observed canonical sequence or a contiguous source fragment covering ≥80%. A sequence match does not establish identical covalent chemistry. The viewer labels sequence-matched reference structures; read both the DRAMP chemical modifications and the PDB study. Source classifications, activity/hemolysis/cytotoxicity text and references are retained as provided; they are not normalized assay measurements or independently re-curated. One source row lacks a DRAMP identifier and is retained under its original imported `None` identifier in an `unassigned` shard; do not infer a real DRAMP accession for it.

No pretrained AMP, toxicity, MIC, secondary-structure or 3D prediction model is bundled. This tool supports research exploration, not clinical decisions.

## Tests

```sh
python3 -m pip install -r requirements.txt
python3 -m unittest -v test_app.py
pnpm install --frozen-lockfile --ignore-scripts
pnpm test
python3 scripts/check_static.py
```

The Python suite has 23 tests. The browser core/worker suite has 30 checks covering descriptor agreement, alignment, quaternion/Kabsch agreement, reflection handling, dataset loading, subdirectory URLs, search/filter/pagination, sharded records, both neighbor modes, structural comparison, uploaded PDB identity fits, exports and invalid inputs, plus wheel round trips, sequence numbering, all 72 supported geometries, resizing, wildcard matching, and invalid wheel inputs. `check_static.py` checks all record shards, structure paths, database equality and the size budget.

## Attribution and licenses

- Application source: **MIT**, see `LICENSE`.
- DRAMP General and source SMILES data: **CC BY 4.0**. Cite Ma et al., *DRAMP 5.0: an update integrating sequence and structure in a novel FASTS format for AI-driven design of antimicrobial peptides*, Nucleic Acids Research (2026), https://doi.org/10.1093/nar/gkag989. Official download page: https://dramp.cpu-bioinfor.org/downloads/. The SQLite, JSON, CSV, static shards and exports are adaptations/repackaging, not an official DRAMP release. Source URLs/hashes and build outcomes are retained in `data/` and `docs/data/`.
- RCSB PDB archival coordinates: **CC0**; retain deposition identifiers and cite individual research papers for scientific use. https://www.rcsb.org/pages/policies.
- **3Dmol.js** is bundled locally with its license and dependency notices in `static/vendor/` and `docs/vendor/`. https://3dmol.org/.

This is an independent DRAMP-inspired research workspace hosted for HySonLab. Its software license does not replace upstream data/library licenses.

## Live helical-wheel design

Fill residues gradually to search after a 250 ms pause. Exact positional matches use the chosen length and treat empty positions as wildcards. Filter by pH-7 charge, GRAVY, hydrophobic percentage, or linked structures. The selected-position residue counts relax only that position while preserving other constraints; they count source records and are not activity predictions. Live search can be paused; manual search remains available. Complete designs show calculated mass, charge, GRAVY, hydrophobic fraction, pI and ideal-helical hydrophobic moment; incomplete designs show assigned-residue counts.

The DOM integration suite exercises editing, debounce, filter validation, pause/manual search, clear/undo, residue suggestions, complete-design properties and the existing workbench handoff. Development tests require Node.js 22+ and pnpm 11.25.0; deployed Pages has no runtime package dependencies.

## Citation

```bibtex
@misc{hy2026ampatlasbrowseronlyworkbench,
      title={AMP Atlas: A Browser-Only Workbench for Antimicrobial Peptide Exploration and Interactive Helical-Wheel Design}, 
      author={Truong-Son Hy},
      year={2026},
      eprint={2610.29096},
      archivePrefix={arXiv},
      primaryClass={cs.CE},
      url={https://arxiv.org/abs/2610.29096}, 
}
```
