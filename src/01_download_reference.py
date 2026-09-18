from pathlib import Path
from urllib.request import Request, urlopen

from config import (
    HBB_END,
    HBB_START,
    RAW_DIR,
    REFERENCE_URL,
    ensure_directories,
    read_single_fasta,
)


def download(url: str, destination: Path) -> None:
    request = Request(url, headers={"User-Agent": "UIU-HBB-FYDP/1.0"})
    with urlopen(request, timeout=120) as response:
        data = response.read()
    if data.lstrip().lower().startswith(b"<!doctype html"):
        raise RuntimeError("NCBI returned a web page instead of FASTA data")
    destination.write_bytes(data)


def main() -> None:
    ensure_directories()
    destination = RAW_DIR / "hbb_reference_grch38_plus.fasta"
    if not destination.exists():
        print("Downloading the GRCh38 HBB genomic reference...")
        download(REFERENCE_URL, destination)
    else:
        print(f"Using existing file: {destination}")

    header, sequence = read_single_fasta(destination)
    expected_length = HBB_END - HBB_START + 1
    if len(sequence) != expected_length:
        raise ValueError(
            f"Expected {expected_length} bases, but downloaded {len(sequence)}. "
            "Do not continue until the reference is corrected."
        )
    if "N" in sequence:
        raise ValueError("Reference contains ambiguous N bases")

    print(f"Reference: {header}")
    print(f"Length: {len(sequence)} bp")
    print(f"Saved to: {destination}")


if __name__ == "__main__":
    main()
