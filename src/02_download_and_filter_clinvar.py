import argparse
import csv
import gzip
import re
from pathlib import Path
from urllib.request import Request, urlopen

from config import (
    ASSEMBLY,
    CLINVAR_VARIANT_SUMMARY_URL,
    HBB_END,
    HBB_START,
    INTERIM_DIR,
    RAW_DIR,
    ensure_directories,
)


OUTPUT_COLUMNS = [
    "variation_id",
    "name",
    "clinical_significance",
    "review_status",
    "phenotype",
    "chromosome",
    "position",
    "reference",
    "alternate",
    "variant_type",
    "rs_id",
]


def download_with_progress(url: str, destination: Path) -> None:
    request = Request(url, headers={"User-Agent": "UIU-HBB-FYDP/1.0"})
    with urlopen(request, timeout=180) as response, destination.open("wb") as out:
        total = int(response.headers.get("Content-Length", "0"))
        downloaded = 0
        while True:
            chunk = response.read(1024 * 1024)
            if not chunk:
                break
            out.write(chunk)
            downloaded += len(chunk)
            if total:
                print(f"\rDownloaded {downloaded / total:.1%}", end="")
            else:
                print(f"\rDownloaded {downloaded / 1024 / 1024:.1f} MiB", end="")
    print()


def valid_dna_allele(value: str) -> bool:
    return bool(value) and bool(re.fullmatch(r"[ACGT]+", value.upper()))


def exact_pathogenic(significance: str) -> bool:
    normalized = significance.strip().lower()
    return normalized in {
        "pathogenic",
        "likely pathogenic",
        "pathogenic/likely pathogenic",
    }


def exact_benign(significance: str) -> bool:
    normalized = significance.strip().lower()
    return normalized in {
        "benign",
        "likely benign",
        "benign/likely benign",
    }


def row_to_output(row: dict[str, str]) -> dict[str, str] | None:
    position_text = row.get("PositionVCF") or row.get("Start") or ""
    try:
        position = int(position_text)
    except ValueError:
        return None

    reference = (row.get("ReferenceAlleleVCF") or row.get("ReferenceAllele") or "").upper()
    alternate = (row.get("AlternateAlleleVCF") or row.get("AlternateAllele") or "").upper()
    if not valid_dna_allele(reference) or not valid_dna_allele(alternate):
        return None
    if not (HBB_START <= position <= HBB_END):
        return None

    return {
        "variation_id": row.get("VariationID", ""),
        "name": row.get("Name", ""),
        "clinical_significance": row.get("ClinicalSignificance", ""),
        "review_status": row.get("ReviewStatus", ""),
        "phenotype": row.get("PhenotypeList", ""),
        "chromosome": row.get("Chromosome", ""),
        "position": str(position),
        "reference": reference,
        "alternate": alternate,
        "variant_type": row.get("Type", ""),
        "rs_id": row.get("RS# (dbSNP)", ""),
    }


def write_rows(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--redownload",
        action="store_true",
        help="Download variant_summary.txt.gz again even if it already exists.",
    )
    args = parser.parse_args()
    ensure_directories()

    source = RAW_DIR / "variant_summary.txt.gz"
    if args.redownload or not source.exists():
        print("Downloading the official ClinVar variant summary. This is a large file...")
        download_with_progress(CLINVAR_VARIANT_SUMMARY_URL, source)
    else:
        print(f"Using existing file: {source}")

    pathogenic: dict[tuple[str, str, str], dict[str, str]] = {}
    benign: dict[tuple[str, str, str], dict[str, str]] = {}
    hbb_rows = 0

    with gzip.open(source, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row.get("Assembly") != ASSEMBLY:
                continue
            if row.get("GeneSymbol", "").strip() != "HBB":
                continue
            if row.get("Chromosome", "").replace("chr", "") != "11":
                continue

            output = row_to_output(row)
            if output is None:
                continue
            hbb_rows += 1
            key = (output["position"], output["reference"], output["alternate"])
            significance = output["clinical_significance"]
            phenotype = output["phenotype"].lower()

            if exact_pathogenic(significance) and "thalassem" in phenotype:
                pathogenic[key] = output
            elif exact_benign(significance):
                benign[key] = output

    pathogenic_path = INTERIM_DIR / "hbb_beta_thal_pathogenic_clinvar.csv"
    benign_path = INTERIM_DIR / "hbb_benign_clinvar.csv"
    write_rows(pathogenic_path, list(pathogenic.values()))
    write_rows(benign_path, list(benign.values()))

    print(f"Usable GRCh38 HBB rows examined: {hbb_rows}")
    print(f"Unique beta-thalassemia P/LP variants: {len(pathogenic)}")
    print(f"Unique benign/likely-benign HBB variants: {len(benign)}")
    print(f"Saved: {pathogenic_path}")
    print(f"Saved: {benign_path}")
    print("Remember: these counts are variants, not patient samples.")


if __name__ == "__main__":
    main()
