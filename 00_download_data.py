"""
Step 0 - Download the public TCGA lung adenocarcinoma (LUAD) data.

Source: cBioPortal DataHub, study "luad_tcga" (TCGA Firehose Legacy).
All files are public, open-access, de-identified summary data.

Run:  python python/00_download_data.py
"""
from pathlib import Path
import urllib.request

BASE = ("https://media.githubusercontent.com/media/cBioPortal/datahub/"
        "master/public/luad_tcga/")

FILES = [
    "data_clinical_patient.txt",   # age, sex, smoking history ...
    "data_clinical_sample.txt",    # sample IDs
    "data_mutations.txt",          # somatic mutations (MAF format)
    "data_cna.txt",                # copy-number calls per gene (GISTIC)
    "data_methylation_hm450.txt",  # DNA methylation, Illumina 450K, per gene
]

out_dir = Path(__file__).resolve().parents[1] / "data" / "raw"
out_dir.mkdir(parents=True, exist_ok=True)

for name in FILES:
    target = out_dir / name
    if target.exists():
        print(f"already present: {name}")
        continue
    print(f"downloading {name} ...")
    urllib.request.urlretrieve(BASE + name, target)
    print(f"  saved {target.stat().st_size / 1e6:.1f} MB")

print("done")
