"""Print human-readable sanity checks without changing any dataset file."""

import pandas as pd

from config import PROCESSED_DIR


def print_table(title: str, table: pd.DataFrame) -> None:
    print(f"\n{title}")
    print(table.to_string())


def main() -> None:
    genotype_path = PROCESSED_DIR / "hbb_three_class_genotypes.csv"
    kmer_path = PROCESSED_DIR / "features_kmer_k4.csv"
    flags_path = PROCESSED_DIR / "mutation_flags.csv"
    combined_path = PROCESSED_DIR / "features_combined_k4.csv"

    for path in (genotype_path, kmer_path, flags_path, combined_path):
        if not path.exists():
            raise FileNotFoundError(f"Missing expected file: {path}")

    data = pd.read_csv(genotype_path).fillna("")
    print("HBB THREE-CLASS DATASET CHECK")
    print(f"Rows: {len(data):,}")
    print(f"Columns: {len(data.columns)}")
    print(f"Unique generated sample IDs: {data['sample_id'].nunique():,}")
    print(f"Unique base population samples: {data['base_sample_id'].nunique():,}")
    print(f"Duplicate sample IDs: {data['sample_id'].duplicated().sum():,}")

    print_table(
        "Class counts:",
        data["class_label"].value_counts().reindex(
            ["Normal", "Carrier", "Affected"], fill_value=0
        ).to_frame("rows"),
    )
    print_table(
        "Split by class:",
        pd.crosstab(data["split"], data["class_label"]).reindex(
            index=["train", "validation", "test"],
            columns=["Normal", "Carrier", "Affected"],
            fill_value=0,
        ),
    )
    print_table(
        "Genotype by label:",
        pd.crosstab(data["class_label"], data["genotype"]),
    )

    length1 = data["allele1_sequence"].str.len()
    length2 = data["allele2_sequence"].str.len()
    print("\nAllele sequence length range:")
    print(f"Allele 1: {length1.min()} to {length1.max()} bp")
    print(f"Allele 2: {length2.min()} to {length2.max()} bp")
    print("Indels can make a reconstructed allele shorter or longer than 1,608 bp.")

    base_split_count = data.groupby("base_sample_id")["split"].nunique()
    print(f"\nBase samples appearing in more than one split: {(base_split_count > 1).sum()}")

    variant_splits: dict[str, set[str]] = {}
    for row in data.itertuples(index=False):
        for variant in (row.pathogenic_variant_1, row.pathogenic_variant_2):
            if variant:
                variant_splits.setdefault(variant, set()).add(row.split)
    cross_split_variants = sum(len(splits) > 1 for splits in variant_splits.values())
    print(f"Pathogenic variants appearing in more than one split: {cross_split_variants}")

    invalid_dna_rows = sum(
        bool(set(sequence.upper()) - set("ACGT"))
        for sequence in data["allele1_sequence"].tolist()
        + data["allele2_sequence"].tolist()
    )
    print(f"Allele sequences containing non-ACGT symbols: {invalid_dna_rows}")

    for path in (kmer_path, flags_path, combined_path):
        header = pd.read_csv(path, nrows=0).columns
        print(f"{path.name}: {len(data):,} rows, {len(header)} columns")

    kmer_columns = [
        column for column in pd.read_csv(kmer_path, nrows=0).columns
        if column.startswith("KMER_")
    ]
    mutation_columns = [
        column for column in pd.read_csv(flags_path, nrows=0).columns
        if column.startswith("MF_")
    ]
    print(f"4-mer feature columns: {len(kmer_columns)} (expected 4^4 = 256)")
    print(f"Mutation-dosage columns used in final rows: {len(mutation_columns)}")

    print("\nInterpretation warning:")
    print("Carrier/Affected rows are simulated genotype proxies, not clinical patients.")
    print("Mutation flags are for explanation/ablation; use k-mer-only as the primary model.")


if __name__ == "__main__":
    main()
