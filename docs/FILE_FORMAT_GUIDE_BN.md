# ব্যবহৃত file format guide

## FASTA — `.fasta`

DNA/protein sequence রাখার সাধারণ text format। প্রথম line `>` দিয়ে header; পরের lineগুলো sequence।

```text
>NC_000011.10:5225464-5227071 ...
TTGCAATGAAAATAAATG...
```

এই project-এ: `data/raw/hbb_reference_grch38_plus.fasta`

## VCF — `.vcf.gz`

Variant Call Format। কোন chromosome position-এ reference allele-এর বদলে কোন alternate allele আছে এবং sample-এর genotype কী—তা রাখে।

মূল field:

| Field | অর্থ |
|---|---|
| `CHROM` | chromosome |
| `POS` | 1-based genomic position |
| `REF` | reference allele |
| `ALT` | alternate allele |
| `GT` | genotype; যেমন `0|1` |

`0` = REF allele, `1` = first ALT allele। `|` phased genotype—কোন allele/haplotype-এ variant আছে তা জানা যায়।

এই project-এ VCF gzip-compressed: `data/raw/hbb_1000g_grch38.vcf.gz`

## TBI — `.tbi`

Compressed VCF-এর tabix index। পুরো VCF না পড়ে `11:5225464-5227071`-এর মতো নির্দিষ্ট region দ্রুত read করতে সাহায্য করে। `.vcf.gz` এবং `.tbi` একই folder-এ রাখবে।

## CSV — `.csv`

Comma-separated table। প্রথম row column name, পরের প্রতিটি row একটি record/sample। VS Code-এ Rainbow CSV extension বা Python pandas দিয়ে দেখা যায়। Excel-এ 1,608 bp sequence-এর বড় column দেখতে অসুবিধা হতে পারে; file নষ্ট না করতে save না করে শুধু view করো।

এই project-এ CSV-এর তিন স্তর:

| স্তর | উদাহরণ | অর্থ |
|---|---|---|
| Interim | `hbb_beta_thal_pathogenic_clinvar.csv` | filtered variant catalogue |
| Processed genotype | `hbb_three_class_genotypes.csv` | দুই allele + label |
| ML feature | `features_kmer_k4.csv` | numeric model input |

## GZ — `.gz`

Gzip compression। ClinVar-এর original bulk table এবং VCF compressed রাখতে ব্যবহার হয়। Windows-এ manually unzip করা বাধ্যতামূলক নয়; Python `gzip`/`pysam` সরাসরি পড়তে পারে।

## Python — `.py`

Executable source code। VS Code-এ খুলে পড়া ও terminal-এ `python filename.py` দিয়ে চালানো যায়। `src/` folder-এ পুরো preprocessing এবং ML code আছে।

## Batch — `.bat`

Windows Command Prompt-এ কয়েকটি command ধারাবাহিকভাবে চালায়। এটি নতুন algorithm নয়; আলাদা আলাদা Python command-এর shortcut।

## JSON — `.json`

Model evaluation metrics structured text হিসেবে save হয়। Training চালানোর পর `results/*_metrics.json` পাওয়া যাবে।

## JOBLIB — `.joblib`

Trained scikit-learn model binary file। Training চালানোর পর `models/` folder-এ পাওয়া যাবে। এটি dataset নয়।

## PNG — `.png`

Confusion matrix image। Training script test prediction থেকে এটি বানায়।
