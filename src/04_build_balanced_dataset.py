import argparse
import csv
import random
import re
from dataclasses import dataclass

import pysam

from config import (
    HBB_END,
    HBB_START,
    INTERIM_DIR,
    PROCESSED_DIR,
    RAW_DIR,
    ensure_directories,
    read_single_fasta,
)


@dataclass(frozen=True)
class Variant:
    position: int
    reference: str
    alternate: str
    name: str

    @property
    def key(self) -> str:
        return f"{self.position}:{self.reference}>{self.alternate}"

    @property
    def interval(self) -> tuple[int, int]:
        return self.position, self.position + len(self.reference) - 1


def is_simple_allele(value: str | None) -> bool:
    return bool(value) and bool(re.fullmatch(r"[ACGT]+", value.upper()))


def read_pathogenic_variants(path) -> list[Variant]:
    variants: dict[str, Variant] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            variant = Variant(
                position=int(row["position"]),
                reference=row["reference"].upper(),
                alternate=row["alternate"].upper(),
                name=row["name"],
            )
            variants[variant.key] = variant
    return list(variants.values())


def overlap(a: Variant, b: Variant) -> bool:
    a_start, a_end = a.interval
    b_start, b_end = b.interval
    return not (a_end < b_start or b_end < a_start)


def can_add(existing: list[Variant], candidate: Variant) -> bool:
    return all(not overlap(item, candidate) for item in existing)


def apply_variants(reference: str, variants: list[Variant]) -> str | None:
    sequence = reference
    for variant in sorted(variants, key=lambda item: item.position, reverse=True):
        offset = variant.position - HBB_START
        observed = sequence[offset : offset + len(variant.reference)]
        if observed != variant.reference:
            return None
        sequence = (
            sequence[:offset]
            + variant.alternate
            + sequence[offset + len(variant.reference) :]
        )
    return sequence


def load_population_haplotypes(vcf_path, pathogenic_keys: set[str]):
    vcf = pysam.VariantFile(str(vcf_path))
    samples = list(vcf.header.samples)
    haplotypes: dict[str, list[list[Variant]]] = {
        sample: [[], []] for sample in samples
    }
    has_known_pathogenic = {sample: False for sample in samples}

    for record in vcf.fetch():
        if not is_simple_allele(record.ref) or not record.alts:
            continue
        for sample in samples:
            genotype = record.samples[sample].get("GT")
            if not genotype or len(genotype) != 2 or None in genotype:
                continue
            for haplotype_index, allele_index in enumerate(genotype):
                if allele_index == 0:
                    continue
                if allele_index - 1 >= len(record.alts):
                    continue
                alternate = record.alts[allele_index - 1]
                if not is_simple_allele(alternate):
                    continue
                variant = Variant(
                    position=record.pos,
                    reference=record.ref.upper(),
                    alternate=alternate.upper(),
                    name=record.id or "1000G_variant",
                )
                haplotypes[sample][haplotype_index].append(variant)
                if variant.key in pathogenic_keys:
                    has_known_pathogenic[sample] = True
    vcf.close()

    return {
        sample: ops
        for sample, ops in haplotypes.items()
        if not has_known_pathogenic[sample]
    }


def choose_compatible_variant(
    rng: random.Random,
    pool: list[Variant],
    existing: list[Variant],
) -> Variant:
    candidates = pool[:]
    rng.shuffle(candidates)
    for candidate in candidates:
        if can_add(existing, candidate):
            return candidate
    raise RuntimeError("No compatible pathogenic variant could be added")


def partition(items: list, rng: random.Random):
    shuffled = items[:]
    rng.shuffle(shuffled)
    n = len(shuffled)
    train_end = max(1, int(n * 0.70))
    validation_end = max(train_end + 1, int(n * 0.85))
    validation_end = min(validation_end, n - 1)
    return {
        "train": shuffled[:train_end],
        "validation": shuffled[train_end:validation_end],
        "test": shuffled[validation_end:],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--n-per-class",
        type=int,
        default=300,
        help="Number of rows to build for each class (default: 300).",
    )
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    ensure_directories()

    reference_path = RAW_DIR / "hbb_reference_grch38_plus.fasta"
    pathogenic_path = INTERIM_DIR / "hbb_beta_thal_pathogenic_clinvar.csv"
    population_path = RAW_DIR / "hbb_1000g_grch38.vcf.gz"
    for required in (reference_path, pathogenic_path, population_path):
        if not required.exists():
            raise FileNotFoundError(f"Missing required input: {required}")

    _, reference = read_single_fasta(reference_path)
    pathogenic = read_pathogenic_variants(pathogenic_path)
    usable_pathogenic = [
        variant
        for variant in pathogenic
        if apply_variants(reference, [variant]) is not None
    ]
    if len(usable_pathogenic) < 12:
        raise RuntimeError(
            "Fewer than 12 compatible P/LP variants remain. Check that ClinVar and "
            "the reference both use GRCh38 plus-strand coordinates."
        )

    pathogenic_keys = {variant.key for variant in usable_pathogenic}
    population = load_population_haplotypes(population_path, pathogenic_keys)
    valid_population = {
        sample: haplotypes
        for sample, haplotypes in population.items()
        if apply_variants(reference, haplotypes[0]) is not None
        and apply_variants(reference, haplotypes[1]) is not None
    }
    if len(valid_population) < args.n_per_class:
        raise RuntimeError(
            f"Only {len(valid_population)} usable population samples were found, but "
            f"{args.n_per_class} were requested. Lower --n-per-class."
        )

    rng = random.Random(args.seed)
    selected_samples = list(valid_population)
    rng.shuffle(selected_samples)
    selected_samples = selected_samples[: args.n_per_class]

    sample_splits = partition(selected_samples, rng)
    variant_splits = partition(usable_pathogenic, rng)
    if any(not values for values in variant_splits.values()):
        raise RuntimeError("A train/validation/test pathogenic-variant split is empty")

    output_path = PROCESSED_DIR / "hbb_three_class_genotypes.csv"
    fieldnames = [
        "sample_id",
        "base_sample_id",
        "split",
        "assembly",
        "region",
        "allele1_sequence",
        "allele2_sequence",
        "pathogenic_variant_1",
        "pathogenic_variant_2",
        "genotype",
        "class_label",
        "synthetic_genotype",
        "label_basis",
    ]

    class_counts = {"Normal": 0, "Carrier": 0, "Affected": 0}
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()

        for split_name, sample_ids in sample_splits.items():
            pool = variant_splits[split_name]
            for sample_id in sample_ids:
                hap1_ops, hap2_ops = valid_population[sample_id]
                normal_hap1 = apply_variants(reference, hap1_ops)
                normal_hap2 = apply_variants(reference, hap2_ops)
                if normal_hap1 is None or normal_hap2 is None:
                    continue

                shared = {
                    "base_sample_id": sample_id,
                    "split": split_name,
                    "assembly": "GRCh38",
                    "region": f"11:{HBB_START}-{HBB_END}",
                    "label_basis": "ClinVar P/LP allele dosage; not a clinical diagnosis",
                }

                writer.writerow(
                    {
                        **shared,
                        "sample_id": f"{sample_id}_N",
                        "allele1_sequence": normal_hap1,
                        "allele2_sequence": normal_hap2,
                        "pathogenic_variant_1": "",
                        "pathogenic_variant_2": "",
                        "genotype": "0/0",
                        "class_label": "Normal",
                        "synthetic_genotype": 0,
                    }
                )
                class_counts["Normal"] += 1

                carrier_on_first = bool(rng.getrandbits(1))
                carrier_existing = hap1_ops if carrier_on_first else hap2_ops
                carrier_variant = choose_compatible_variant(rng, pool, carrier_existing)
                carrier_hap1_ops = hap1_ops + ([carrier_variant] if carrier_on_first else [])
                carrier_hap2_ops = hap2_ops + ([] if carrier_on_first else [carrier_variant])
                carrier_hap1 = apply_variants(reference, carrier_hap1_ops)
                carrier_hap2 = apply_variants(reference, carrier_hap2_ops)
                if carrier_hap1 is None or carrier_hap2 is None:
                    raise RuntimeError("Failed to construct a carrier sequence")

                writer.writerow(
                    {
                        **shared,
                        "sample_id": f"{sample_id}_C",
                        "allele1_sequence": carrier_hap1,
                        "allele2_sequence": carrier_hap2,
                        "pathogenic_variant_1": carrier_variant.key,
                        "pathogenic_variant_2": "",
                        "genotype": "0/1",
                        "class_label": "Carrier",
                        "synthetic_genotype": 1,
                    }
                )
                class_counts["Carrier"] += 1

                affected_variant_1 = choose_compatible_variant(rng, pool, hap1_ops)
                affected_variant_2 = (
                    affected_variant_1
                    if rng.random() < 0.5 and can_add(hap2_ops, affected_variant_1)
                    else choose_compatible_variant(rng, pool, hap2_ops)
                )
                affected_hap1 = apply_variants(reference, hap1_ops + [affected_variant_1])
                affected_hap2 = apply_variants(reference, hap2_ops + [affected_variant_2])
                if affected_hap1 is None or affected_hap2 is None:
                    raise RuntimeError("Failed to construct an affected-genotype sequence")

                writer.writerow(
                    {
                        **shared,
                        "sample_id": f"{sample_id}_A",
                        "allele1_sequence": affected_hap1,
                        "allele2_sequence": affected_hap2,
                        "pathogenic_variant_1": affected_variant_1.key,
                        "pathogenic_variant_2": affected_variant_2.key,
                        "genotype": (
                            "1/1"
                            if affected_variant_1.key == affected_variant_2.key
                            else "1/2"
                        ),
                        "class_label": "Affected",
                        "synthetic_genotype": 1,
                    }
                )
                class_counts["Affected"] += 1

    print(f"Compatible ClinVar P/LP variants: {len(usable_pathogenic)}")
    print(f"Usable 1000 Genomes background samples: {len(valid_population)}")
    print(f"Final class counts: {class_counts}")
    print(f"Saved: {output_path}")
    print("Carrier and Affected rows are simulated genotypes and must be reported as such.")


if __name__ == "__main__":
    main()
