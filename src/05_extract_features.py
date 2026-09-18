import argparse
import csv
import itertools
import re
from collections import Counter

import pandas as pd

from config import PROCESSED_DIR, ensure_directories


META_COLUMNS = ["sample_id", "base_sample_id", "split", "class_label"]


def all_kmers(k: int) -> list[str]:
    return ["".join(parts) for parts in itertools.product("ACGT", repeat=k)]


def normalized_kmer_counts(sequence: str, k: int) -> Counter:
    sequence = sequence.upper()
    valid_windows = [
        sequence[index : index + k]
        for index in range(max(0, len(sequence) - k + 1))
        if set(sequence[index : index + k]) <= set("ACGT")
    ]
    counts = Counter(valid_windows)
    denominator = len(valid_windows) or 1
    for key in list(counts):
        counts[key] /= denominator
    return counts


def safe_column_name(variant_key: str) -> str:
    return "MF_" + re.sub(r"[^A-Za-z0-9]+", "_", variant_key).strip("_")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--k", type=int, default=4)
    args = parser.parse_args()
    if args.k < 2 or args.k > 6:
        raise ValueError("Use k between 2 and 6 for this starter project")
    ensure_directories()

    input_path = PROCESSED_DIR / "hbb_three_class_genotypes.csv"
    if not input_path.exists():
        raise FileNotFoundError(f"Missing input: {input_path}")

    dataset = pd.read_csv(input_path).fillna("")
    vocabulary = all_kmers(args.k)
    kmer_rows = []
    allele_pair_rows = []
    mutation_keys: set[str] = set()

    for row in dataset.itertuples(index=False):
        count1 = normalized_kmer_counts(row.allele1_sequence, args.k)
        count2 = normalized_kmer_counts(row.allele2_sequence, args.k)
        output = {
            "sample_id": row.sample_id,
            "base_sample_id": row.base_sample_id,
            "split": row.split,
            "class_label": row.class_label,
        }
        pair_output = {key: output[key] for key in META_COLUMNS}
        for kmer in vocabulary:
            frequency1 = count1.get(kmer, 0.0)
            frequency2 = count2.get(kmer, 0.0)
            output[f"KMER_{kmer}"] = (frequency1 + frequency2) / 2
            # Min/max retain allele differences without depending on haplotype order.
            pair_output[f"PAIR_MIN_{kmer}"] = min(frequency1, frequency2)
            pair_output[f"PAIR_MAX_{kmer}"] = max(frequency1, frequency2)
        kmer_rows.append(output)
        allele_pair_rows.append(pair_output)
        if row.pathogenic_variant_1:
            mutation_keys.add(row.pathogenic_variant_1)
        if row.pathogenic_variant_2:
            mutation_keys.add(row.pathogenic_variant_2)

    kmer_frame = pd.DataFrame(kmer_rows)
    allele_pair_frame = pd.DataFrame(allele_pair_rows)
    flag_rows = []
    sorted_keys = sorted(mutation_keys)
    for row in dataset.itertuples(index=False):
        flags = {
            "sample_id": row.sample_id,
            "base_sample_id": row.base_sample_id,
            "split": row.split,
            "class_label": row.class_label,
        }
        observed = Counter(
            key
            for key in (row.pathogenic_variant_1, row.pathogenic_variant_2)
            if key
        )
        for key in sorted_keys:
            flags[safe_column_name(key)] = observed.get(key, 0)
        flag_rows.append(flags)
    flag_frame = pd.DataFrame(flag_rows)

    kmer_path = PROCESSED_DIR / f"features_kmer_k{args.k}.csv"
    allele_pair_path = PROCESSED_DIR / f"features_allele_pair_k{args.k}.csv"
    flags_path = PROCESSED_DIR / "mutation_flags.csv"
    combined_path = PROCESSED_DIR / f"features_combined_k{args.k}.csv"
    kmer_frame.to_csv(kmer_path, index=False)
    allele_pair_frame.to_csv(allele_pair_path, index=False)
    flag_frame.to_csv(flags_path, index=False)
    combined = kmer_frame.merge(
        flag_frame.drop(columns=["base_sample_id", "split", "class_label"]),
        on="sample_id",
        validate="one_to_one",
    )
    combined.to_csv(combined_path, index=False)

    print(f"Rows: {len(dataset)}")
    print(f"k-mer features: {len(vocabulary)}")
    print(f"Mutation-position dosage flags: {len(sorted_keys)}")
    print(f"Primary feature file: {kmer_path}")
    print(f"Allele-pair feature file: {allele_pair_path}")
    print(f"Flags for ablation/explanation: {flags_path}")
    print(f"Combined secondary experiment: {combined_path}")
    print("Do not present combined-feature accuracy without a leakage warning and ablation study.")


if __name__ == "__main__":
    main()
