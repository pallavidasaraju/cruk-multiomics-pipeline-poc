"""
Multi-Omics Cancer Research Data Generator
Generates synthetic genomics, proteomics, and imaging datasets for POC demonstration.
"""
import pandas as pd
import numpy as np
from pathlib import Path

np.random.seed(42)
DATA_DIR = Path(__file__).parent

def generate_genomics_data(n_samples=100):
    """Generate synthetic Illumina sequencing QC metrics"""
    samples = []
    for i in range(n_samples):
        sample_id = f"CRUK_{i+1:04d}"
        patient_group = np.random.choice(["Control", "Tumor", "Metastatic"], p=[0.3, 0.5, 0.2])
        total_reads = np.random.randint(20_000_000, 80_000_000)
        mapped_reads = int(total_reads * np.random.uniform(0.85, 0.98))
        q30_bases = np.random.uniform(85, 98)
        gc_content = np.random.uniform(38, 52)
        duplication_rate = np.random.uniform(5, 25)
        coverage_mean = np.random.uniform(30, 150)
        variant_count = np.random.randint(50, 500)

        samples.append({
            "sample_id": sample_id,
            "patient_group": patient_group,
            "total_reads_M": total_reads / 1_000_000,
            "mapped_reads_pct": (mapped_reads / total_reads) * 100,
            "q30_bases_pct": q30_bases,
            "gc_content_pct": gc_content,
            "duplication_rate_pct": duplication_rate,
            "mean_coverage": coverage_mean,
            "variant_count": variant_count,
            "pass_qc": (q30_bases >= 90) and (duplication_rate <= 20) and (coverage_mean >= 30)
        })

    df = pd.DataFrame(samples)
    df.to_csv(DATA_DIR / "genomics_qc.csv", index=False)
    print(f"Generated genomics data: {len(df)} samples")
    return df

def generate_proteomics_data(n_samples=100):
    """Generate synthetic mass spectrometry protein abundance data"""
    samples = []
    proteins = ["TP53", "KRAS", "EGFR", "MYC", "PTEN", "AKT1", "PIK3CA", "BRCA1", "BRCA2", "BRAF"]

    for i in range(n_samples):
        sample_id = f"CRUK_{i+1:04d}"
        patient_group = np.random.choice(["Control", "Tumor", "Metastatic"], p=[0.3, 0.5, 0.2])
        proteins_detected = np.random.randint(3000, 8000)
        peptides_identified = np.random.randint(10000, 25000)
        ms1_intensity = np.random.uniform(1e7, 1e9)
        ms2_spectra = np.random.randint(15000, 50000)

        # Generate fold-changes for key cancer proteins
        protein_fc = {}
        for protein in proteins:
            if patient_group == "Control":
                fc = np.random.uniform(0.9, 1.1)
            elif patient_group == "Tumor":
                fc = np.random.uniform(1.5, 4.0) if protein in ["TP53", "KRAS", "MYC"] else np.random.uniform(0.5, 1.5)
            else:  # Metastatic
                fc = np.random.uniform(2.0, 6.0) if protein in ["TP53", "KRAS", "MYC"] else np.random.uniform(0.3, 2.0)
            protein_fc[f"{protein}_fc"] = fc

        samples.append({
            "sample_id": sample_id,
            "patient_group": patient_group,
            "proteins_detected": proteins_detected,
            "peptides_identified": peptides_identified,
            "ms1_intensity": ms1_intensity,
            "ms2_spectra": ms2_spectra,
            **protein_fc
        })

    df = pd.DataFrame(samples)
    df.to_csv(DATA_DIR / "proteomics_abundance.csv", index=False)
    print(f"Generated proteomics data: {len(df)} samples")
    return df

def generate_imaging_data(n_samples=100):
    """Generate synthetic advanced imaging metrics (H&E, IHC, confocal)"""
    samples = []
    for i in range(n_samples):
        sample_id = f"CRUK_{i+1:04d}"
        patient_group = np.random.choice(["Control", "Tumor", "Metastatic"], p=[0.3, 0.5, 0.2])

        if patient_group == "Control":
            cell_density = np.random.randint(500, 1500)
            tumor_pct = 0
            necrosis_pct = 0
            mitotic_index = np.random.uniform(0, 2)
        elif patient_group == "Tumor":
            cell_density = np.random.randint(1500, 4000)
            tumor_pct = np.random.uniform(30, 80)
            necrosis_pct = np.random.uniform(0, 20)
            mitotic_index = np.random.uniform(5, 20)
        else:  # Metastatic
            cell_density = np.random.randint(2000, 5000)
            tumor_pct = np.random.uniform(50, 95)
            necrosis_pct = np.random.uniform(10, 40)
            mitotic_index = np.random.uniform(15, 50)

        samples.append({
            "sample_id": sample_id,
            "patient_group": patient_group,
            "cell_density_per_mm2": cell_density,
            "tumor_area_pct": tumor_pct,
            "necrosis_area_pct": necrosis_pct,
            "mitotic_index": mitotic_index,
            "ki67_positive_pct": np.random.uniform(10, 80) if patient_group != "Control" else np.random.uniform(1, 10),
            "cd8_infiltration_score": np.random.randint(0, 4),
            "image_quality_score": np.random.uniform(7, 10)
        })

    df = pd.DataFrame(samples)
    df.to_csv(DATA_DIR / "imaging_metrics.csv", index=False)
    print(f"Generated imaging data: {len(df)} samples")
    return df

def generate_integrated_metadata(n_samples=100):
    """Generate sample metadata linking all modalities"""
    samples = []
    for i in range(n_samples):
        sample_id = f"CRUK_{i+1:04d}"
        patient_group = np.random.choice(["Control", "Tumor", "Metastatic"], p=[0.3, 0.5, 0.2])

        samples.append({
            "sample_id": sample_id,
            "patient_id": f"PT_{i//2 + 1:03d}",  # 2 samples per patient on average
            "patient_group": patient_group,
            "tissue_type": np.random.choice(["Breast", "Lung", "Colon", "Prostate", "Ovarian"]),
            "age": np.random.randint(35, 85),
            "sex": np.random.choice(["M", "F"]),
            "stage": np.random.choice(["I", "II", "III", "IV"]) if patient_group != "Control" else "N/A",
            "genomics_available": True,
            "proteomics_available": np.random.choice([True, False], p=[0.9, 0.1]),
            "imaging_available": np.random.choice([True, False], p=[0.85, 0.15]),
            "collection_date": pd.Timestamp("2025-01-01") + pd.Timedelta(days=np.random.randint(0, 365))
        })

    df = pd.DataFrame(samples)
    df.to_csv(DATA_DIR / "sample_metadata.csv", index=False)
    print(f"Generated metadata: {len(df)} samples")
    return df

if __name__ == "__main__":
    print("Generating multi-omics cancer research datasets...")
    generate_genomics_data()
    generate_proteomics_data()
    generate_imaging_data()
    generate_integrated_metadata()
    print("\nAll datasets generated successfully!")
    print(f"  -> genomics_qc.csv")
    print(f"  -> proteomics_abundance.csv")
    print(f"  -> imaging_metrics.csv")
    print(f"  -> sample_metadata.csv")
