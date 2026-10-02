import argparse
import itertools
from pathlib import Path

import joblib
import pandas as pd

from config import MODELS_DIR


EXPECTED_LENGTH = 1608
LABEL_ORDER = ["Normal", "Carrier", "Affected"]


def read_fasta_records(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header = None
    chunks: list[str] = []

    with path.open(encoding="utf-8-sig") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(chunks).upper()))
                header = line[1:].strip()
                chunks = []
                if not header:
                    raise ValueError(f"Empty FASTA header at line {line_number}")
            else:
                if header is None:
                    raise ValueError("FASTA sequence data must follow a >header line")
                chunks.append("".join(line.split()))

    if header is not None:
        records.append((header, "".join(chunks).upper()))
    if len(records) != 2:
        raise ValueError(f"Expected exactly 2 FASTA records (two alleles); found {len(records)}")

    for name, sequence in records:
        if not sequence:
            raise ValueError(f"FASTA record {name!r} has no sequence")
        invalid = set(sequence) - set("ACGTN")
        if invalid:
            raise ValueError(f"FASTA record {name!r} has invalid symbols: {sorted(invalid)}")
        if len(sequence) != EXPECTED_LENGTH:
            raise ValueError(
                f"FASTA record {name!r} is {len(sequence)} bp; expected the "
                f"{EXPECTED_LENGTH} bp GRCh38 HBB interval (NC_000011.10:5225464-5227071)"
            )
    return records


def normalized_kmer_frequencies(sequence: str, k: int) -> dict[str, float]:
    windows = [
        sequence[index : index + k]
        for index in range(len(sequence) - k + 1)
        if set(sequence[index : index + k]) <= set("ACGT")
    ]
    denominator = len(windows) or 1
    counts = {"".join(parts): 0 for parts in itertools.product("ACGT", repeat=k)}
    for window in windows:
        counts[window] += 1
    return {kmer: count / denominator for kmer, count in counts.items()}


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Predict a project class from a two-allele HBB FASTA file."
    )
    parser.add_argument("fasta", type=Path, help="FASTA file containing exactly two HBB allele records")
    parser.add_argument(
        "--model",
        type=Path,
        default=MODELS_DIR / "best_kmer_k4.joblib",
        help="Model bundle created by src/08_compare_models.py",
    )
    args = parser.parse_args()

    if not args.model.exists():
        raise FileNotFoundError(
            f"Model bundle not found: {args.model}. First run "
            "python src/08_compare_models.py --features kmer --k 4"
        )

    bundle = joblib.load(args.model)
    if not isinstance(bundle, dict) or not {"estimator", "feature_columns", "representation", "k"} <= bundle.keys():
        raise ValueError("Expected a model bundle created by src/08_compare_models.py")
    if bundle["representation"] != "kmer":
        raise ValueError("This FASTA predictor supports a kmer model bundle; train with --features kmer")

    k = int(bundle["k"])
    alleles = read_fasta_records(args.fasta)
    freq1 = normalized_kmer_frequencies(alleles[0][1], k)
    freq2 = normalized_kmer_frequencies(alleles[1][1], k)
    all_features = {f"KMER_{name}": (freq1[name] + freq2[name]) / 2 for name in freq1}
    feature_columns = bundle["feature_columns"]
    missing = [column for column in feature_columns if column not in all_features]
    if missing:
        raise ValueError(f"Model expects unsupported feature columns, for example: {missing[:3]}")

    prediction = bundle["estimator"].predict(pd.DataFrame([[all_features[name] for name in feature_columns]], columns=feature_columns))[0]
    print(f"Allele 1: {alleles[0][0]} ({len(alleles[0][1])} bp)")
    print(f"Allele 2: {alleles[1][0]} ({len(alleles[1][1])} bp)")
    print(f"Predicted project class: {prediction}")
    print("This is an educational prediction for simulated genotype proxy classes, not a clinical result.")


if __name__ == "__main__":
    main()
