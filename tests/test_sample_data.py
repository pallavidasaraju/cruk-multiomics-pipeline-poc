"""
Tests for multi-omics sample data generation
"""
import pytest
import pandas as pd
import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))
from sample_data import (
    generate_genomics_data,
    generate_proteomics_data,
    generate_imaging_data,
    generate_integrated_metadata
)

def test_genomics_structure():
    df = generate_genomics_data(n_samples=50)
    assert len(df) == 50
    assert "sample_id" in df.columns
    assert "patient_group" in df.columns
    assert set(df["patient_group"].unique()) <= {"Control", "Tumor", "Metastatic"}

def test_genomics_qc_logic():
    df = generate_genomics_data(n_samples=100)
    # Passing samples should have Q30 >= 90, duplication <= 20, coverage >= 30
    passing = df[df["pass_qc"]]
    assert (passing["q30_bases_pct"] >= 90).all()
    assert (passing["duplication_rate_pct"] <= 20).all()
    assert (passing["mean_coverage"] >= 30).all()

def test_genomics_read_counts():
    df = generate_genomics_data(n_samples=100)
    assert (df["total_reads_M"] >= 20).all()
    assert (df["total_reads_M"] <= 80).all()
    assert (df["mapped_reads_pct"] >= 85).all()
    assert (df["mapped_reads_pct"] <= 100).all()

def test_proteomics_structure():
    df = generate_proteomics_data(n_samples=50)
    assert len(df) == 50
    assert "sample_id" in df.columns
    assert "proteins_detected" in df.columns
    protein_fc_cols = [c for c in df.columns if "_fc" in c]
    assert len(protein_fc_cols) == 10  # 10 cancer proteins

def test_proteomics_fold_changes():
    df = generate_proteomics_data(n_samples=100)
    # Metastatic samples should have higher TP53/KRAS/MYC fold-changes on average
    metastatic = df[df["patient_group"] == "Metastatic"]
    control = df[df["patient_group"] == "Control"]
    assert metastatic["TP53_fc"].mean() > control["TP53_fc"].mean()
    assert metastatic["KRAS_fc"].mean() > control["KRAS_fc"].mean()

def test_proteomics_detection_range():
    df = generate_proteomics_data(n_samples=100)
    assert (df["proteins_detected"] >= 3000).all()
    assert (df["proteins_detected"] <= 8000).all()
    assert (df["peptides_identified"] >= 10000).all()

def test_imaging_structure():
    df = generate_imaging_data(n_samples=50)
    assert len(df) == 50
    assert "sample_id" in df.columns
    assert "tumor_area_pct" in df.columns
    assert "mitotic_index" in df.columns

def test_imaging_tumor_logic():
    df = generate_imaging_data(n_samples=200)
    control = df[df["patient_group"] == "Control"]
    tumor = df[df["patient_group"] == "Tumor"]
    metastatic = df[df["patient_group"] == "Metastatic"]

    # Controls should have 0 tumor area
    assert (control["tumor_area_pct"] == 0).all()
    # Metastatic should have higher tumor area than non-metastatic on average
    assert metastatic["tumor_area_pct"].mean() > tumor["tumor_area_pct"].mean()

def test_imaging_quality_scores():
    df = generate_imaging_data(n_samples=100)
    assert (df["image_quality_score"] >= 7).all()
    assert (df["image_quality_score"] <= 10).all()

def test_metadata_structure():
    df = generate_integrated_metadata(n_samples=50)
    assert len(df) == 50
    assert "sample_id" in df.columns
    assert "patient_id" in df.columns
    assert "tissue_type" in df.columns
    assert set(df["genomics_available"].unique()) == {True}

def test_metadata_data_availability():
    df = generate_integrated_metadata(n_samples=200)
    # All samples should have genomics
    assert df["genomics_available"].all()
    # Most should have proteomics and imaging
    assert df["proteomics_available"].mean() >= 0.8
    assert df["imaging_available"].mean() >= 0.7

def test_sample_id_consistency():
    genomics = generate_genomics_data(n_samples=100)
    proteomics = generate_proteomics_data(n_samples=100)
    imaging = generate_imaging_data(n_samples=100)
    metadata = generate_integrated_metadata(n_samples=100)

    # All should have same sample_ids
    assert set(genomics["sample_id"]) == set(proteomics["sample_id"])
    assert set(genomics["sample_id"]) == set(imaging["sample_id"])
    assert set(genomics["sample_id"]) == set(metadata["sample_id"])

def test_patient_group_distribution():
    genomics = generate_genomics_data(n_samples=300)
    group_counts = genomics["patient_group"].value_counts(normalize=True)
    # Should roughly match probabilities [0.3, 0.5, 0.2]
    assert group_counts["Control"] > 0.2
    assert group_counts["Tumor"] > 0.4
    assert group_counts["Metastatic"] > 0.1
