from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = PROJECT_ROOT / "data"
RAW_DIR = DATA_DIR / "raw"
INTERIM_DIR = DATA_DIR / "interim"
PROCESSED_DIR = DATA_DIR / "processed"
RESULTS_DIR = PROJECT_ROOT / "results"
MODELS_DIR = PROJECT_ROOT / "models"

ASSEMBLY = "GRCh38"
CHROMOSOME = "11"
CHROMOSOME_ACCESSION = "NC_000011.10"
HBB_START = 5_225_464  # 1-based, inclusive
HBB_END = 5_227_071    # 1-based, inclusive

REFERENCE_URL = (
    "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi"
    "?id=NC_000011.10&db=nuccore&report=fasta&retmode=text"
    "&from=5225464&to=5227071&strand=1"
)

CLINVAR_VARIANT_SUMMARY_URL = (
    "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/"
    "variant_summary.txt.gz"
)

THOUSAND_GENOMES_CHR11_URL = (
    "https://hgdownload.soe.ucsc.edu/gbdb/hg38/1000Genomes/"
    "ALL.chr11.shapeit2_integrated_snvindels_v2a_27022019."
    "GRCh38.phased.vcf.gz"
)


def ensure_directories() -> None:
    for path in (RAW_DIR, INTERIM_DIR, PROCESSED_DIR, RESULTS_DIR, MODELS_DIR):
        path.mkdir(parents=True, exist_ok=True)


def read_single_fasta(path: Path) -> tuple[str, str]:
    header = ""
    sequence_parts: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    raise ValueError(f"Expected one FASTA record in {path}")
                header = line[1:]
            else:
                sequence_parts.append(line.upper())
    sequence = "".join(sequence_parts)
    if not header or not sequence:
        raise ValueError(f"Invalid FASTA file: {path}")
    if set(sequence) - set("ACGTN"):
        raise ValueError(f"Unexpected nucleotide symbols in {path}")
    return header, sequence
