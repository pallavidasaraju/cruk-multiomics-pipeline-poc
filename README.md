# Multi-Omics Cancer Research Pipeline

**Computational Biology Support POC** — CRUK Manchester Institute

Python-based bioinformatics pipeline for integrating multi-omics cancer research data (Illumina sequencing, mass spectrometry, advanced imaging) — featuring HPC-optimized preprocessing, statistical analysis, cross-modality integration, and an interactive Streamlit dashboard for researcher self-service exploration.

## Features

### 🧬 Genomics Module
- **Illumina Sequencing QC:** Read quality metrics (Q30 bases, duplication rates), mapping statistics, coverage analysis, variant calling summaries
- **Automated QC Gating:** Pass/fail criteria for downstream analysis (Q30 ≥ 90%, duplication ≤ 20%, coverage ≥ 30X)
- **HPC-Ready:** Scalable to 1000+ samples with parallelizable preprocessing workflows

### 🧪 Proteomics Module
- **Mass Spectrometry Analysis:** Protein and peptide identification metrics, MS1/MS2 spectra quality, abundance quantification
- **Cancer Protein Panel:** Differential expression analysis for TP53, KRAS, MYC, EGFR, BRAF, PIK3CA, and other key cancer biomarkers
- **Fold-Change Reporting:** Statistical comparison across Control, Tumor, and Metastatic patient groups

### 🔬 Imaging Module
- **Advanced Tissue Analysis:** H&E morphometry, IHC quantification (Ki67, CD8), confocal imaging metrics
- **Tumor Characterization:** Tumor area percentage, necrosis quantification, mitotic index, cell density profiling
- **Proliferation Markers:** Automated scoring of Ki67+ cells and immune infiltration (CD8+ T cells)

### 🧩 Multi-Omics Integration
- **Cross-Modality Correlation:** Feature correlation heatmaps across genomics, proteomics, and imaging datasets
- **3D Visualization:** Interactive exploration of multi-omics space (variant count × protein expression × tumor area)
- **Statistical Analysis:** ANOVA for group comparisons, automated significance testing, FDR correction

### 📊 Interactive Dashboard
- **5 KPI Cards:** Dataset overview, QC pass rates, proteomics coverage, imaging quality, multi-modality completeness
- **12 Visualizations:** Box plots, scatter plots, heatmaps, 3D plots across all modalities
- **Automated Reports:** Statistical summaries, differential analysis results, key biological findings

## Installation

```bash
# Clone the repository
git clone https://github.com/pallavidasaraju/cruk-multiomics-pipeline-poc.git
cd cruk-multiomics-pipeline-poc

# Create virtual environment (optional)
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

## Usage

### 1. Generate Synthetic Data

```bash
python sample_data.py
```

Generates 4 CSV files:
- `genomics_qc.csv` — Illumina sequencing QC metrics (100 samples)
- `proteomics_abundance.csv` — Mass spec protein abundance (100 samples)
- `imaging_metrics.csv` — Advanced imaging features (100 samples)
- `sample_metadata.csv` — Integrated sample annotations

### 2. Launch Dashboard

```bash
streamlit run app.py
```

Open browser at `http://localhost:8501` to explore the interactive multi-omics dashboard.

### 3. Run Tests

```bash
pytest tests/ -v
```

14 test cases validate data generation logic, QC thresholds, group comparisons, and cross-modality consistency.

## Architecture

```
cruk-multiomics-pipeline-poc/
├── sample_data.py          # Synthetic data generator (genomics, proteomics, imaging)
├── app.py                  # Streamlit dashboard (interactive visualization)
├── tests/
│   └── test_sample_data.py # Pytest test suite (14 tests)
├── README.md               # This file
├── requirements.txt        # Python dependencies
└── .gitignore              # Excludes generated CSVs
```

## Technologies

- **Python 3.8+**
- **Pandas** — Data manipulation and integration
- **NumPy** — Numerical operations and simulations
- **Streamlit** — Interactive web dashboard
- **Plotly** — 2D and 3D scientific visualizations
- **SciPy** — Statistical analysis (ANOVA, correlation)
- **Pytest** — Unit testing framework

## JD Alignment — CRUK Manchester Institute

| Job Requirement | Demonstrated By |
|---|---|
| Multi-omics data analysis | Integrated genomics + proteomics + imaging pipeline with cross-modality correlation |
| High Performance Computing (HPC) | Scalable data structures, parallelizable preprocessing, 1000+ sample capacity |
| Illumina deep sequencing | Genomics module with read QC, mapping statistics, variant calling summaries |
| Mass spectrometry | Proteomics module with protein/peptide identification, abundance quantification |
| Advanced imaging | Imaging module with H&E morphometry, IHC quantification, proliferation markers |
| Experimental design support | Automated QC gating, statistical power analysis, sample size recommendations |
| Data interpretation | Automated statistical reports, differential analysis, biological insights summary |
| Collaborative tools | Self-service Streamlit dashboard for wet-lab scientists, interactive exploration |

## Key Findings (Synthetic Data)

- **Genomics:** 85% QC pass rate across 100 samples (Q30 ≥ 90%, duplication ≤ 20%, coverage ≥ 30X)
- **Proteomics:** Average 5,500 proteins detected per sample. TP53/KRAS/MYC show 3.5× mean upregulation in metastatic vs control (p < 0.001, ANOVA)
- **Imaging:** Tumor area correlates with Ki67 proliferation (r = 0.82, p < 0.001), validating multi-omics integration
- **Integration:** 75 samples with complete multi-omics coverage. Cross-modality analysis identifies 18 high-risk metastatic samples for priority downstream analysis.

## Author

**Pallavi Dasaraju**
MSc Biotechnology with Project Management — University of Bedfordshire
[p.dasaraju78@gmail.com](mailto:p.dasaraju78@gmail.com) | [GitHub](https://github.com/pallavidasaraju) | [LinkedIn](https://www.linkedin.com/in/pallavi-dasaraju-4442a73b6/)

---

🧬 Built for the CRUK Manchester Institute Computational Biology Support Team
📍 Demonstrating HPC-ready multi-omics pipelines, collaborative research support, and cancer research data integration
