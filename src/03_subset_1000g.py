import pysam

from config import (
    HBB_END,
    HBB_START,
    RAW_DIR,
    THOUSAND_GENOMES_CHR11_URL,
    ensure_directories,
)


def main() -> None:
    ensure_directories()
    output_path = RAW_DIR / "hbb_1000g_grch38.vcf.gz"

    print("Opening the remote phased chromosome 11 VCF...")
    print("Only the HBB region will be copied; the full chromosome file is not downloaded.")
    try:
        source = pysam.VariantFile(
            THOUSAND_GENOMES_CHR11_URL,
            index_filename=THOUSAND_GENOMES_CHR11_URL + ".tbi",
        )
    except Exception as exc:
        raise RuntimeError(
            "Could not open the remote 1000 Genomes VCF. Check the internet connection "
            "and try again. The README also gives a browser-based Data Slicer option."
        ) from exc

    if "11" in source.header.contigs:
        contig = "11"
    elif "chr11" in source.header.contigs:
        contig = "chr11"
    else:
        raise RuntimeError("Chromosome 11 was not found in the VCF header")

    header = source.header.copy()
    sample_count = len(header.samples)
    variant_count = 0
    with pysam.VariantFile(str(output_path), "wz", header=header) as output:
        for record in source.fetch(contig, HBB_START - 1, HBB_END):
            output.write(record)
            variant_count += 1
    source.close()
    pysam.tabix_index(str(output_path), preset="vcf", force=True)

    print(f"Samples in VCF header: {sample_count}")
    print(f"Variant sites in the HBB interval: {variant_count}")
    print(f"Saved: {output_path}")


if __name__ == "__main__":
    main()
