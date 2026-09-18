# Code walkthrough

সব code `src/` folder-এ আছে। এই file পড়ার সময় VS Code-এ পাশে সংশ্লিষ্ট `.py` file খুলে রাখো।

## `config.py`

এক জায়গায় project path, HBB coordinate এবং source URL রাখা হয়েছে।

- `PROJECT_ROOT`, `RAW_DIR`, `PROCESSED_DIR`: file কোথায় যাবে
- `HBB_START`, `HBB_END`: HBB genomic interval
- `REFERENCE_URL`, `CLINVAR_VARIANT_SUMMARY_URL`, `THOUSAND_GENOMES_CHR11_URL`: source
- `ensure_directories()`: প্রয়োজনীয় folder বানায়
- `read_single_fasta()`: FASTA header ও sequence পড়ে এবং nucleotide validate করে

## `01_download_reference.py`

NCBI URL open করে FASTA text download করে। তারপর `read_single_fasta()` দিয়ে validate করে length print করে।

বলবে: “I used a fixed GRCh38 genomic coordinate so all later VCF and ClinVar positions refer to the same assembly.”

## `02_download_and_filter_clinvar.py`

1. বড় `variant_summary.txt.gz` download বা existing copy নেয়।
2. `gzip.open(..., 'rt')` দিয়ে compressed table stream করে।
3. GRCh38, HBB, chromosome 11 এবং coordinate filter করে।
4. `exact_pathogenic()` বা `exact_benign()` significance যাচাই করে।
5. `(position, reference, alternate)` key দিয়ে duplicate variant বাদ দেয়।
6. দুইটি filtered CSV save করে।

গুরুত্বপূর্ণ: `Pathogenic` substring দিয়ে loose matching করা হয়নি; exact approved category ব্যবহার করা হয়েছে।

## `03_subset_1000g.py`

`pysam.VariantFile()` দিয়ে remote indexed VCF open করে। `fetch()` শুধু HBB interval আনে। Output bgzip VCF লেখা এবং `pysam.tabix_index()` দিয়ে `.tbi` বানানো হয়।

বলবে: “I did not need the full chromosome file because the index supports region-based retrieval.”

## `04_build_balanced_dataset.py`

এই script সবচেয়ে গুরুত্বপূর্ণ।

- `Variant` dataclass: position, REF, ALT এবং name একসঙ্গে রাখে।
- `apply_variants()`: genomic coordinate-কে reference-string offset-এ convert করে; REF না মিললে variant reject করে; indel-এর জন্য right-to-left apply করে।
- `load_population_haplotypes()`: VCF-এর phased `GT` দেখে প্রতি sample-এর allele 1 এবং allele 2 variant list বানায়।
- `partition()`: sample এবং pathogenic variant আলাদাভাবে train/validation/test-এ ভাগ করে।
- main loop: প্রতি base sample-এর Normal, Carrier ও Affected-proxy row লেখে।

Random selection reproducible রাখতে default `--seed 42`। একই input ও seed হলে একই constructed dataset হবে।

## `05_extract_features.py`

এখানেই k-mer conversion হয়েছে।

```python
def all_kmers(k):
    return ["".join(parts) for parts in itertools.product("ACGT", repeat=k)]
```

এই অংশ A/C/G/T-এর সব possible combination বানায়। `k=4` হলে 256টি।

```python
sequence[index:index + k]
```

এই slicing দিয়ে sequence-এর উপর এক base করে sliding window চলে। Count-কে total valid window দিয়ে ভাগ করে normalized frequency হয়। তারপর:

```python
(count1.get(kmer, 0.0) + count2.get(kmer, 0.0)) / 2
```

দুই allele-এর frequency average করে sample feature বানানো হয়। `safe_column_name()` mutation key-কে CSV-safe `MF_...` column-এ বদলায়।

একই loop-এ `PAIR_MIN_...` ও `PAIR_MAX_...` column-ও বানানো হয়। এগুলো দুই allele-এর frequency-এর ছোট/বড় মান রাখে, তাই allele order না বদলিয়েও পার্থক্য কিছুটা ধরে। Output `features_allele_pair_k4.csv` নতুন sequence-only experiment-এর জন্য।

## `06_train_baseline.py`

1. `class_label`, IDs ও `split` বাদ দিয়ে numeric feature column নেয়।
2. Pre-assigned train/test row আলাদা করে।
3. Random Forest এবং StandardScaler + RBF-SVM fit করে।
4. Test prediction থেকে classification report ও confusion matrix করে।
5. metrics JSON, model JOBLIB এবং matrix PNG save করে।

Validation split এখন future hyperparameter tuning-এর জন্য reserved; baseline code test-এর আগে manual tuning করে না। এটা বললে code-এর বর্তমান সীমা সৎভাবে বোঝানো হয়।

## `08_compare_models.py`

`--features kmer` অথবা `--features allele_pair` দিয়ে sequence-only input নেয়। Dummy, Logistic Regression, KNN, Random Forest, Extra Trees ও RBF-SVM শুধু training split-এ train হয়ে validation macro-F1-তে তুলনা হয়। জয়ী model train+validation-এ আবার fit হয়, তারপর একবার test-এ মাপা হয়। Validation leaderboard CSV, final test JSON/PNG এবং feature-column order-সহ JOBLIB save হয়। পুরো command ও নতুন model যোগ করার উদাহরণ `RUN_AND_EXTEND_MODELS_BN.md`-এ আছে।

## `07_inspect_dataset.py`

Dataset না বদলে sanity check করে: row/column count, class/split distribution, sequence length, genotype-label mapping, synthetic flag, duplicate ID এবং feature dimension print করে। Demo শুরুতে এটি চালাও।

## Code নিজে বুঝেছ কি না যাচাই

নিজে এই তিনটি ছোট experiment করো:

1. `--k 3` দিয়ে feature বানিয়ে দেখো—64 k-mer হওয়ার কথা।
2. `--n-per-class 50` দিলে মোট 150 row কেন হয় ব্যাখ্যা করো।
3. একটি Normal, Carrier, Affected row-এর `genotype`, variant column ও label পাশাপাশি দেখো।

Actual report result-এর file overwrite করতে না চাইলে experiment-এর আগে project folder-এর একটি copy রাখবে।
