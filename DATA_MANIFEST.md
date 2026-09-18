# Data Manifest

Collection and validation date: 2026-09-03

Documentation, inspection code, and included k-mer baseline outputs verified: 2026-09-07

Allele-pair features and validation-based model comparison added and run on the same included 900-row dataset: 2026-09-18 (Python 3.14, scikit-learn 1.9.1). No raw data was downloaded again for this run. Exact new scores are in `docs/EXAMPLE_RESULTS_BN.md` and `results/comparison_*`.

Genome build and interval:

- Assembly: GRCh38.p14
- Chromosome accession: NC_000011.10
- HBB interval: chromosome 11, 5,225,464 to 5,227,071 inclusive
- Reference length: 1,608 bp
- Stored orientation: chromosome plus strand, to match VCF REF/ALT alleles

Source counts obtained by the included scripts:

- 1000 Genomes phased GRCh38 VCF header: 2,548 samples
- 1000 Genomes variant sites in the HBB interval: 71
- ClinVar usable GRCh38 HBB rows examined: 1,745
- Unique ClinVar beta-thalassemia Pathogenic/Likely pathogenic variants retained: 279
- ClinVar P/LP variants compatible with the exact downloaded reference: 271
- Unique ClinVar benign/likely-benign HBB variants retained: 836
- 1000 Genomes samples usable as no-known-P/LP background genotypes: 2,494

Constructed dataset:

- Normal rows: 300
- Carrier rows: 300
- Affected-genotype rows: 300
- Total rows: 900
- Split: 630 train, 135 validation, 135 test
- Unique pathogenic variant sets do not overlap between train, validation and test
- Each base population sample is confined to one split

Important interpretation:

- ClinVar counts are variant counts, not patient counts.
- Normal rows are 1000 Genomes population genotypes with no matched ClinVar P/LP HBB variant in this pipeline.
- Carrier and Affected rows are computationally constructed monoallelic and biallelic genotypes.
- These labels are genotype proxies and are not confirmed clinical diagnoses.
- The dataset is suitable for method development and coursework only, not clinical use.
