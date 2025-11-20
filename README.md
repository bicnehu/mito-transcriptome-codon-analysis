# Mitonuclear codon-usage dataset (120 parasitic helminths)

**Short description.** Curated coding sequences and derived codon-bias metrics from **120 parasitic helminth species** with matched **mitochondrial genomes** and **nuclear transcriptomes**, plus fully reproducible scripts (R 4.3.x / Python 3.10).

**Repository:** https://github.com/bicnehu/mito-transcriptome-codon-analysis  
**Data DOIs:** Dryad (data) — _ADD DOI_ • Zenodo (code release) — _ADD DOI_  
**Contact:** devbioinfo@gmail.com | devendra@nehu.ac.in

---

## Contents

```
.
├── data_raw/
│   ├── genbank_mito_cds/                 # curated mitochondrial CDS per species
│   └── wormbase_nuclear_cds/             # curated nuclear CDS per species
├── tables/
│   ├── mito_rscu_120x60.csv              # 120×60 RSCU matrix (mitochondria)
│   ├── nuclear_rscu_120x60.csv           # 120×60 RSCU matrix (nuclear)
│   ├── enc_gc3_genelevel.tsv             # ENC and GC3 per gene
│   ├── gc12_gc3_species.tsv              # GC12 and GC3 per species
│   ├── skews_species.tsv                 # AT- and GC-skew per species
│   └── cai_13genes_by_species.csv        # CAI for 13 mt genes × 120 spp
├── code/
│   ├── R/                                # PCA, PERMANOVA, plots
│   └── py/                               # parsing, skews, plotting
├── env/
│   ├── renv.lock                         # R environment (pinned)
│   └── requirements.txt                  # Python environment (pinned)
└── docs/                                 # figure scripts and README assets
```

---

## Quick start

### 1) Recreate environments
```bash
# Python 3.10
python -m venv .venv && source .venv/bin/activate
pip install -r env/requirements.txt

# R 4.3.x
R -q <<'RSCRIPT'
install.packages("renv"); renv::restore()
RSCRIPT
```

### 2) Rebuild key tables
```bash
# Count codons and compute RSCU (mito and nuclear)
python code/py/make_codoncounters_and_rscu.py   --mito data_raw/genbank_mito_cds/   --nuclear data_raw/wormbase_nuclear_cds/   --outdir tables/
```

### 3) Reproduce analyses & figures
```bash
# Skews, ENC–GC3, neutrality plots
python code/py/run_skews_and_enc_gc3.py --in tables/ --out docs/

# PCA, PERMANOVA, biplots
Rscript code/R/run_pca_and_permanova.R   --mito tables/mito_rscu_120x60.csv   --nuclear tables/nuclear_rscu_120x60.csv   --out docs/
```

---

## Methods (how the data were built)

### Data acquisition & curation
- **Sources.** Complete mitochondrial CDS per species from **GenBank**; matched nuclear CDS from **WormBase** (latest stable release per taxon).  
- **QC.** Verified genetic code, removed frameshifts/early stops and sequences with >1% ambiguous bases; minus-strand genes were reverse-complemented.

### Metrics and derived tables
- **Codon counts & RSCU.** Codon counts were generated with **EMBOSS** `cusp`; counts → **RSCU** per species/compartment (60 sense codons; ATG & stops excluded).  
- **Composition & skews.** Per species: GC, **GC3**, and strand-asymmetry **AT-/GC-skews**.  
- **ENC & neutrality inputs.** **ENC** and **GC3** per gene; species-level **GC12 vs GC3** for neutrality plots.  
- **CAI profiling (mitogenome).** **Codon Adaptation Index (CAI)** computed for 13 mitochondrial PCGs using **CodonW**; CAI is the geometric mean of codon adaptiveness against a reference set.

### Multivariate analyses & tests
- **PCA of codon usage.** RSCU matrices were centered/scaled by codon and analyzed with **FactoMineR**; visualization via **factoextra**.  
- **Order-level clustering (PERMANOVA).** PERMANOVA on Euclidean distances of z-scored RSCU (10,000 perms) with **vegan** `adonis2`; dispersion checked with `betadisper`.

### Implementation
- **R (4.3.x):** `vegan 2.6-4`, `FactoMineR 2.8`, `factoextra 1.0.7`.  
- **Python (3.10):** `pandas 2.1`, `numpy 1.26`, `matplotlib 3.8`.

---

## Data dictionary (key files)

| File | Description |
|---|---|
| `tables/mito_rscu_120x60.csv` | Relative synonymous codon usage per mitochondrial dataset (rows: species; cols: 60 codons). |
| `tables/nuclear_rscu_120x60.csv` | RSCU for matched nuclear CDS. |
| `tables/enc_gc3_genelevel.tsv` | Gene-wise ENC and GC3 values for ENC–GC3 analysis. |
| `tables/gc12_gc3_species.tsv` | Species-level GC12 & GC3 for neutrality plots. |
| `tables/skews_species.tsv` | AT- and GC-skew by species and compartment. |
| `tables/cai_13genes_by_species.csv` | CAI for 13 mitogenome genes across species. |

---

## Reuse notes

- Suitable for comparative tests of mutation vs selection, mitonuclear coupling/decoupling, benchmarking codon-bias metrics, phylogenetic or ML analyses, and codon-aware design for expression studies.
- Please cite the **Dryad dataset DOI** (data) and **Zenodo DOI** (code release) in derivative works.

---

## Reproducibility

- Random seeds: `set.seed(12345)` (R); `PYTHONHASHSEED=0` (Python).  
- Environments are pinned (`renv.lock`, `requirements.txt`).  
- We recommend archiving a tagged GitHub **Release** and enabling **Zenodo** to mint a DOI for that release (GitHub↔Zenodo integration).

---

## How to regenerate metrics (CLI examples)

**Codon counts & RSCU via EMBOSS `cusp`:**
```bash
# Example: run on a FASTA of CDS for one species
cusp -sequence data_raw/genbank_mito_cds/SpeciesA_cds.fna -outfile tmp/speciesA_cusp.tab
# Convert to RSCU (script aggregates counts across CDS)
python code/py/cusp_to_rscu.py tmp/speciesA_cusp.tab > tmp/speciesA_rscu.csv
```

**ENC/GC3 and CAI via CodonW:**
```bash
# ENC and GC3 for all CDS in a FASTA
codonw speciesA_cds.fna -enc -gc3s -nomenu -silent

# CAI for mitogenome genes (using species-specific reference)
codonw mito13genes_speciesA.fna -cai -nomenu
```

**PCA & PERMANOVA in R:**
```r
# PCA
library(FactoMineR); library(factoextra)
rscu <- read.csv("tables/mito_rscu_120x60.csv", row.names = 1)
res.pca <- PCA(scale(rscu, scale = TRUE), graph = FALSE)
# PERMANOVA
library(vegan)
distM <- dist(scale(rscu, scale = TRUE), method = "euclidean")
adonis2(distM ~ Order, data = metadata, permutations = 10000)
```

---

## Licensing & ethics

- **Sequences** are from public repositories (GenBank, WormBase).  
- **Data package** (Dryad): recommend **CC0** for maximal reuse; **code** (Zenodo): use a permissive software license (e.g., MIT/GPL) appropriate to your needs.  
- **Ethics:** no human subjects or live-animal work; secondary analysis of public sequences.

---

## Citation

Please cite:  
1) **Dataset (Dryad)** — _ADD DOI once minted_.  
2) **Software release (Zenodo)** — _ADD DOI once minted_.  
3) This repository: `bicnehu/mito-transcriptome-codon-analysis` (add a `CITATION.cff`).

---

## Acknowledgements

We thank the Department of Biotechnology (DBT)–sponsored Bioinformatics Centre at North-Eastern Hill University, Shillong, and the National Network Project @ Biodiversity Informatics Centre, NEHU, for computational resources.# Mitonuclear codon-usage dataset (120 parasitic helminths)

**Short description.** Curated coding sequences and derived codon-bias metrics from **120 parasitic helminth species** with matched **mitochondrial genomes** and **nuclear transcriptomes**, plus fully reproducible scripts (R 4.3.x / Python 3.10).

**Repository:** https://github.com/bicnehu/mito-transcriptome-codon-analysis  
**Data DOIs:** Dryad (data) — _ADD DOI_ • Zenodo (code release) — _ADD DOI_  
**Contact:** devbioinfo@gmail.com | devendra@nehu.ac.in

---

## Contents

```
.
├── data_raw/
│   ├── genbank_mito_cds/                 # curated mitochondrial CDS per species
│   └── wormbase_nuclear_cds/             # curated nuclear CDS per species
├── tables/
│   ├── mito_rscu_120x60.csv              # 120×60 RSCU matrix (mitochondria)
│   ├── nuclear_rscu_120x60.csv           # 120×60 RSCU matrix (nuclear)
│   ├── enc_gc3_genelevel.tsv             # ENC and GC3 per gene
│   ├── gc12_gc3_species.tsv              # GC12 and GC3 per species
│   ├── skews_species.tsv                 # AT- and GC-skew per species
│   └── cai_13genes_by_species.csv        # CAI for 13 mt genes × 120 spp
├── code/
│   ├── R/                                # PCA, PERMANOVA, plots
│   └── py/                               # parsing, skews, plotting
├── env/
│   ├── renv.lock                         # R environment (pinned)
│   └── requirements.txt                  # Python environment (pinned)
└── docs/                                 # figure scripts and README assets
```

---

## Quick start

### 1) Recreate environments
```bash
# Python 3.10
python -m venv .venv && source .venv/bin/activate
pip install -r env/requirements.txt

# R 4.3.x
R -q <<'RSCRIPT'
install.packages("renv"); renv::restore()
RSCRIPT
```

### 2) Rebuild key tables
```bash
# Count codons and compute RSCU (mito and nuclear)
python code/py/make_codoncounters_and_rscu.py   --mito data_raw/genbank_mito_cds/   --nuclear data_raw/wormbase_nuclear_cds/   --outdir tables/
```

### 3) Reproduce analyses & figures
```bash
# Skews, ENC–GC3, neutrality plots
python code/py/run_skews_and_enc_gc3.py --in tables/ --out docs/

# PCA, PERMANOVA, biplots
Rscript code/R/run_pca_and_permanova.R   --mito tables/mito_rscu_120x60.csv   --nuclear tables/nuclear_rscu_120x60.csv   --out docs/
```

---

## Methods (how the data were built)

### Data acquisition & curation
- **Sources.** Complete mitochondrial CDS per species from **GenBank**; matched nuclear CDS from **WormBase** (latest stable release per taxon).  
- **QC.** Verified genetic code, removed frameshifts/early stops and sequences with >1% ambiguous bases; minus-strand genes were reverse-complemented.

### Metrics and derived tables
- **Codon counts & RSCU.** Codon counts were generated with **EMBOSS** `cusp`; counts → **RSCU** per species/compartment (60 sense codons; ATG & stops excluded).  
- **Composition & skews.** Per species: GC, **GC3**, and strand-asymmetry **AT-/GC-skews**.  
- **ENC & neutrality inputs.** **ENC** and **GC3** per gene; species-level **GC12 vs GC3** for neutrality plots.  
- **CAI profiling (mitogenome).** **Codon Adaptation Index (CAI)** computed for 13 mitochondrial PCGs using **CodonW**; CAI is the geometric mean of codon adaptiveness against a reference set.

### Multivariate analyses & tests
- **PCA of codon usage.** RSCU matrices were centered/scaled by codon and analyzed with **FactoMineR**; visualization via **factoextra**.  
- **Order-level clustering (PERMANOVA).** PERMANOVA on Euclidean distances of z-scored RSCU (10,000 perms) with **vegan** `adonis2`; dispersion checked with `betadisper`.

### Implementation
- **R (4.3.x):** `vegan 2.6-4`, `FactoMineR 2.8`, `factoextra 1.0.7`.  
- **Python (3.10):** `pandas 2.1`, `numpy 1.26`, `matplotlib 3.8`.

---

## Data dictionary (key files)

| File | Description |
|---|---|
| `tables/mito_rscu_120x60.csv` | Relative synonymous codon usage per mitochondrial dataset (rows: species; cols: 60 codons). |
| `tables/nuclear_rscu_120x60.csv` | RSCU for matched nuclear CDS. |
| `tables/enc_gc3_genelevel.tsv` | Gene-wise ENC and GC3 values for ENC–GC3 analysis. |
| `tables/gc12_gc3_species.tsv` | Species-level GC12 & GC3 for neutrality plots. |
| `tables/skews_species.tsv` | AT- and GC-skew by species and compartment. |
| `tables/cai_13genes_by_species.csv` | CAI for 13 mitogenome genes across species. |

---

## Reuse notes

- Suitable for comparative tests of mutation vs selection, mitonuclear coupling/decoupling, benchmarking codon-bias metrics, phylogenetic or ML analyses, and codon-aware design for expression studies.
- Please cite the **Dryad dataset DOI** (data) and **Zenodo DOI** (code release) in derivative works.

---

## Reproducibility

- Random seeds: `set.seed(12345)` (R); `PYTHONHASHSEED=0` (Python).  
- Environments are pinned (`renv.lock`, `requirements.txt`).  
- We recommend archiving a tagged GitHub **Release** and enabling **Zenodo** to mint a DOI for that release (GitHub↔Zenodo integration).

---

## How to regenerate metrics (CLI examples)

**Codon counts & RSCU via EMBOSS `cusp`:**
```bash
# Example: run on a FASTA of CDS for one species
cusp -sequence data_raw/genbank_mito_cds/SpeciesA_cds.fna -outfile tmp/speciesA_cusp.tab
# Convert to RSCU (script aggregates counts across CDS)
python code/py/cusp_to_rscu.py tmp/speciesA_cusp.tab > tmp/speciesA_rscu.csv
```

**ENC/GC3 and CAI via CodonW:**
```bash
# ENC and GC3 for all CDS in a FASTA
codonw speciesA_cds.fna -enc -gc3s -nomenu -silent

# CAI for mitogenome genes (using species-specific reference)
codonw mito13genes_speciesA.fna -cai -nomenu
```

**PCA & PERMANOVA in R:**
```r
# PCA
library(FactoMineR); library(factoextra)
rscu <- read.csv("tables/mito_rscu_120x60.csv", row.names = 1)
res.pca <- PCA(scale(rscu, scale = TRUE), graph = FALSE)
# PERMANOVA
library(vegan)
distM <- dist(scale(rscu, scale = TRUE), method = "euclidean")
adonis2(distM ~ Order, data = metadata, permutations = 10000)
```

---

## Licensing & ethics

- **Sequences** are from public repositories (GenBank, WormBase).  
- **Data package** (Dryad): recommend **CC0** for maximal reuse; **code** (Zenodo): use a permissive software license (e.g., MIT/GPL) appropriate to your needs.  
- **Ethics:** no human subjects or live-animal work; secondary analysis of public sequences.

---

## Citation

Please cite:  
1) **Dataset (Dryad)** — _ADD DOI once minted_.  
2) **Software release (Zenodo)** — _ADD DOI once minted_.  
3) This repository: `bicnehu/mito-transcriptome-codon-analysis` (add a `CITATION.cff`).

---

## Acknowledgements

We thank the Department of Biotechnology (DBT)–sponsored Bioinformatics Centre at North-Eastern Hill University, Shillong, and the National Network Project @ Biodiversity Informatics Centre, NEHU, for computational resources.
