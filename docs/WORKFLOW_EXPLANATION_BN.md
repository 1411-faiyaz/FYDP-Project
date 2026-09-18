# পুরো workflow সহজ ভাষায়

## Research question

HBB genomic sequence-এর pattern দেখে একটি ML model কি `Normal`, `Carrier` এবং `Affected-genotype proxy` আলাদা করতে পারে?

## Data flow

| ধাপ | Input | কাজ | Output |
|---|---|---|---|
| 1 | NCBI GRCh38 | HBB genomic interval download | FASTA reference |
| 2 | NCBI ClinVar | HBB + GRCh38 + thalassemia + P/LP filter | Pathogenic variant CSV |
| 3 | 1000 Genomes | শুধু HBB interval subset | Population VCF |
| 4 | FASTA + VCF + ClinVar | দুই allele reconstruct ও genotype proxy বানানো | 900-row labelled CSV |
| 5 | দুই allele sequence | normalized 4-mer এবং mutation dosage | Feature CSV files |
| 6 | Feature CSV + label | Random Forest ও RBF-SVM training/testing | Metrics, model, confusion matrix |

## ধাপ 1 — reference কেন দরকার

`NC_000011.10:5225464-5227071` হলো GRCh38 assembly-এর chromosome 11-এ HBB genomic region। দৈর্ঘ্য 1,608 bp। Variant apply করার সময় position, REF এবং ALT ঠিক আছে কি না—এটি reference দিয়ে যাচাই হয়।

`NM_000518.5` HBB-এর curated mRNA/transcript accession। HGVS নাম যেমন `c.92+5G>C` বোঝাতে transcript দরকার, কিন্তু এই model-এর sequence input genomic DNA। তাই mRNA dataset দিয়ে model বানানো হয়নি।

## ধাপ 2 — disease তথ্য কোথা থেকে

ClinVar-এর official `variant_summary.txt.gz` থেকে শুধু:

- Assembly = GRCh38
- Gene = HBB
- Chromosome = 11
- HBB region-এর মধ্যে position
- phenotype text-এ thalassemia
- significance = Pathogenic / Likely pathogenic / Pathogenic-Likely pathogenic
- VCF-compatible DNA allele

রাখা হয়েছে। 279 unique P/LP variant পাওয়া গিয়েছিল; exact reference-এর সঙ্গে 271 compatible ছিল। এগুলো রোগীর সংখ্যা নয়।

## ধাপ 3 — normal background কোথা থেকে

1000 Genomes phased VCF-এ 2,548 sample ছিল এবং HBB interval-এ 71 variant site ছিল। প্রত্যেক sample-এর genotype থেকে দুইটি haplotype/allele reconstruct করা হয়। যাদের genotype-এ pipeline-এর matched known ClinVar P/LP HBB variant ছিল, তাদের normal-background pool থেকে বাদ দেওয়া হয়। 2,494 usable background পাওয়া যায়।

`Normal` মানে এখানে **no-known-P/LP-HBB proxy**; কোনো ব্যক্তির সম্পূর্ণ medical health নিশ্চিত করা হয়নি।

## ধাপ 4 — তিনটি class কীভাবে তৈরি

একই selected base population sample থেকে তিনটি row:

| suffix | genotype | pathogenic allele dosage | label |
|---|---:|---:|---|
| `_N` | `0/0` | 0 | Normal proxy |
| `_C` | `0/1` | 1 | Carrier proxy |
| `_A` | `1/1` বা `1/2` | 2 | Affected-genotype proxy |

Carrier row-তে এক allele-এ এবং Affected row-তে দুই allele-এ split-specific ClinVar P/LP variant computationally বসানো হয়। এই কারণে `synthetic_genotype=1` রাখা হয়েছে।

300 base sample × 3 class = 900 row। এটি 900 independent patient নয়; base sample 300টি।

## ধাপ 5 — 4-mer কীভাবে হয়

`k=4` হলে সম্ভব pattern সংখ্যা:

\[
4^4 = 256
\]

কারণ DNA alphabet হলো A, C, G, T। উদাহরণ: `ACGTAC` sequence থেকে 4-mer window হলো `ACGT`, `CGTA`, `GTAC`। প্রতিটি 4-mer-এর count মোট valid window দিয়ে ভাগ করে normalized frequency করা হয়। দুই allele-এর একই k-mer frequency average করা হয়:

\[
f_{sample}(x)=\frac{f_{allele1}(x)+f_{allele2}(x)}{2}
\]

তাই `features_kmer_k4.csv`-এ 4 metadata column + 256 k-mer column = 260 column।

## Mutation flag কী

`MF_...` column নির্দিষ্ট pathogenic variant sample-এর কত allele-এ আছে তা দেখায়:

- 0 = নেই
- 1 = এক allele-এ
- 2 = দুই allele-এ

এটি interpretation/ablation-এর জন্য। যেহেতু P/LP dosage থেকেই label তৈরি, mutation flags primary input দিলে target leakage হতে পারে। তাই main result `features_kmer_k4.csv` থেকে report করবে।

## ধাপ 6 — supervised ML কীভাবে বুঝে

Training row-তে feature-এর সঙ্গে `class_label` দেওয়া আছে। Model feature pattern ও label-এর relationship শেখে। Test row-এর label prediction-এর সময় input হিসেবে দেওয়া হয় না; শেষে true label-এর সঙ্গে prediction compare করে metric বের করা হয়।

Current split:

- Train: 630 rows
- Validation: 135 rows
- Test: 135 rows

একটি base sample শুধু একটি split-এ আছে এবং pathogenic variant set-ও split-এর মধ্যে overlap করে না। এতে data leakage কমে।

## Result বোঝার নিয়ম

শুধু accuracy বলবে না। বলবে:

- macro-F1: তিন class-কে সমান গুরুত্ব দিয়ে overall balance
- precision: একটি label prediction-এর কতগুলো ঠিক
- recall: সত্যিকারের class-এর কতগুলো model ধরেছে
- confusion matrix: কোন class-কে কোন class হিসেবে ভুল করেছে

High score দেখলেই clinical diagnosis দাবি করা যাবে না, কারণ dataset simulated genotype proxy।
