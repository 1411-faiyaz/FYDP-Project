# Column dictionary

## `hbb_three_class_genotypes.csv`

| Column | অর্থ |
|---|---|
| `sample_id` | unique generated row ID; suffix `_N`, `_C`, `_A` |
| `base_sample_id` | 1000 Genomes background sample ID |
| `split` | train, validation অথবা test |
| `assembly` | coordinate system; এখানে GRCh38 |
| `region` | HBB genomic interval |
| `allele1_sequence` | reconstructed first HBB allele DNA |
| `allele2_sequence` | reconstructed second HBB allele DNA |
| `pathogenic_variant_1` | first inserted P/LP variant key; format `position:REF>ALT` |
| `pathogenic_variant_2` | second inserted P/LP variant key, যদি থাকে |
| `genotype` | `0/0`, `0/1`, `1/1` বা `1/2` proxy notation |
| `class_label` | supervised target: Normal, Carrier, Affected |
| `synthetic_genotype` | 1 হলে pathogenic genotype computationally constructed |
| `label_basis` | label-এর provenance/limitation statement |

একই `base_sample_id`-এর তিনটি row থাকতে পারে, তবে সবগুলো একই split-এ থাকে।

## `features_kmer_k4.csv`

- 900 row, 260 column
- প্রথম 4টি metadata: `sample_id`, `base_sample_id`, `split`, `class_label`
- বাকি 256টি numeric feature: `KMER_AAAA` থেকে সব possible 4-mer
- value = দুই allele-এর normalized frequency-এর average

## `mutation_flags.csv`

- 900 row, 252 column
- প্রথম 4টি metadata
- বাকি 248টি `MF_...` dosage column
- value 0, 1 বা 2

271 compatible ClinVar P/LP variant থাকলেও final selected 900 rows-এ 248 unique pathogenic variant ব্যবহৃত হয়েছে, তাই 248 mutation feature column।

## `features_combined_k4.csv`

- 900 row, 508 column
- 4 metadata + 256 k-mer + 248 mutation flag
- secondary ablation only; leakage warning ছাড়া main accuracy হিসেবে দেখাবে না

## `hbb_beta_thal_pathogenic_clinvar.csv`

| Column | অর্থ |
|---|---|
| `variation_id` | ClinVar variation identifier |
| `name` | reported HGVS/variant name |
| `clinical_significance` | Pathogenic/Likely pathogenic category |
| `review_status` | ClinVar review level text |
| `phenotype` | associated condition names |
| `chromosome`, `position` | GRCh38 genomic location |
| `reference`, `alternate` | REF এবং ALT allele |
| `variant_type` | SNV/indel ইত্যাদি |
| `rs_id` | dbSNP ID, পাওয়া গেলে |
