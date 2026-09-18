# HBB FYDP Data Starter

এই project-এর লক্ষ্য হলো একই genome build ব্যবহার করে HBB genomic reference, ClinVar mutation catalogue এবং 1000 Genomes population genotype থেকে একটি reproducible তিন-class research dataset তৈরি করা।

**প্রথমবার হলে আগে `START_HERE_BN.md` পড়ো।** সেখানে quick demo এবং কোন explanation file কোন order-এ পড়বে দেওয়া আছে। `src/` folder-এ ব্যবহৃত সম্পূর্ণ code আছে; `docs/` folder-এ code walkthrough, file formats, viva Q&A এবং report methodology আছে।

**নতুন model train করতে:** `docs/RUN_AND_EXTEND_MODELS_BN.md` পড়ো। `src/08_compare_models.py` validation split-এ কয়েকটি model তুলনা করে, জয়ী model-কে train+validation-এ refit করে এবং শেষে held-out test মাপে। `src/05_extract_features.py` এখন allele-pair feature-ও লেখে।

## প্রথমে যে সিদ্ধান্তটি বুঝতে হবে

Online-এ বড়, open, ready-made `Normal / Carrier / Affected` HBB FASTA dataset সাধারণত পাওয়া যায় না। ClinVar patient sequence dataset নয়; এটি variant এবং clinical interpretation-এর catalogue। 1000 Genomes population genotype দেয়, কিন্তু confirmed beta-thalassemia phenotype দেয় না। তাই এই starter project carrier এবং affected genotype **simulate** করে। এটি clinical patient cohort নয় এবং report-এ অবশ্যই `variant-informed simulated genotype dataset` বলতে হবে।

Real clinical classification দাবি করতে হলে hospital/controlled repository থেকে genotype-এর সঙ্গে confirmed phenotype নিতে হবে, ethical approval ও data-use permission মানতে হবে।

এই package-এ 2026-09-03 তারিখে তৈরি 900-row example dataset ও feature files দেওয়া আছে। Exact source count এবং provenance `DATA_MANIFEST.md`-এ আছে। Research-এর final run-এর আগে scripts আবার চালিয়ে current source version সংগ্রহ করবে।

## ব্যবহৃত official/primary sources

1. NCBI HBB Gene: https://www.ncbi.nlm.nih.gov/datasets/gene/3043/
2. GRCh38 HBB reference region: `NC_000011.10:5225464-5227071`, plus genomic strand used internally so that VCF REF/ALT matches.
3. NCBI ClinVar data access: https://www.ncbi.nlm.nih.gov/clinvar/docs/maintenance_use/
4. ClinVar `variant_summary.txt.gz`: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
5. 1000 Genomes/IGSR data: https://www.internationalgenome.org/data/
6. Remote GRCh38 chromosome 11 VCF mirror: https://hgdownload.soe.ucsc.edu/gbdb/hg38/1000Genomes/
7. IthaGenes HBB page for cross-checking, not automatic labels: https://www.ithanet.eu/db/ithagenes?geneID=10

`NM_000518.5` is the HBB transcript reference used in HGVS names such as `c.92+5G>C`. It is useful for reporting mutation names, but this project uses the genomic HBB region as model input.

## Software

Use VS Code, Python 3.11 or 3.12, and the VS Code Python extension. No GPU is required for the included Random Forest and SVM baselines.

## Folder খোলা ও environment তৈরি

1. Zip extract করে `HBB_FYDP_Starter` folder VS Code-এ open করো।
2. VS Code থেকে `Terminal > New Terminal` খোলো।
3. Windows Command Prompt terminal হলে চালাও:

```bat
py -m venv .venv
.venv\Scripts\activate.bat
python -m pip install --upgrade pip
pip install -r requirements.txt
```

PowerShell activation block হলে VS Code terminal profile থেকে Command Prompt select করলেই execution-policy পরিবর্তন করতে হবে না।

## Step 1: Reference download

```bat
python src\01_download_reference.py
```

Expected output:

```text
data/raw/hbb_reference_grch38_plus.fasta
Length: 1608 bp
```

এটি একটি baseline reference, হাজার sample-এর dataset নয়।

## Step 2: ClinVar download ও filter

```bat
python src\02_download_and_filter_clinvar.py
```

এটি official ClinVar bulk file download করবে এবং শুধু নিচের row রাখবে:

- Assembly: GRCh38
- Gene: HBB
- HBB coordinate: 5225464-5227071
- Phenotype-এ `thalassemia`
- Clinical significance exact `Pathogenic`, `Likely pathogenic`, অথবা `Pathogenic/Likely pathogenic`
- VCF-compatible REF/ALT allele

Output:

```text
data/interim/hbb_beta_thal_pathogenic_clinvar.csv
data/interim/hbb_benign_clinvar.csv
```

Script শেষে exact row count print করবে। এটি **unique variant count**, patient count নয়। ClinVar নিয়মিত update হয়, তাই report-এ collection date ও printed count লিখবে।

## Step 3: 1000 Genomes থেকে শুধু HBB region নেওয়া

```bat
python src\03_subset_1000g.py
```

Script 627 MB chromosome file পুরোটা download না করে indexed VCF থেকে শুধু HBB interval fetch করার চেষ্টা করবে। Output:

```text
data/raw/hbb_1000g_grch38.vcf.gz
data/raw/hbb_1000g_grch38.vcf.gz.tbi
```

এটি VCF header-এর exact sample count এবং HBB interval-এর variant-site count print করবে।

Remote access কাজ না করলে Ensembl Data Slicer ব্যবহার করো:

1. https://www.ensembl.org/Homo_sapiens/Tools/DataSlicer খোলো।
2. VCF URL হিসেবে README-এর remote chromosome 11 VCF link দাও।
3. Region হিসেবে প্রথমে `11:5225464-5227071` দাও। কাজ না করলে `chr11:5225464-5227071` চেষ্টা করো।
4. Output file `hbb_1000g_grch38.vcf.gz` নামে `data/raw/` folder-এ রাখো।

## Step 4: Balanced three-class genotype dataset তৈরি

প্রথম run ছোট রাখো:

```bat
python src\04_build_balanced_dataset.py --n-per-class 100
```

সব ঠিক থাকলে final experimental run:

```bat
python src\04_build_balanced_dataset.py --n-per-class 300
```

Script যা করে:

1. 1000 Genomes sample-এর দুই phased HBB haplotype reconstruct করে।
2. যেসব population sample-এ known ClinVar P/LP HBB variant আছে, সেগুলো normal background থেকে বাদ দেয়।
3. একই population background দিয়ে একটি no-P/LP row, একটি monoallelic P/LP row এবং একটি biallelic P/LP row তৈরি করে।
4. Label দেয় `Normal`, `Carrier`, `Affected`।
5. `synthetic_genotype` এবং `label_basis` column রেখে দেয় যাতে provenance লুকানো না হয়।
6. Pathogenic variants আলাদা train, validation ও test pool-এ ভাগ করে; test mutation training-এ না দেখার ব্যবস্থা করে।

Output:

```text
data/processed/hbb_three_class_genotypes.csv
```

একটি row-তে দুইটি allele sequence থাকবে। কারণ single haploid FASTA দিয়ে Carrier বনাম Affected নির্ধারণ করা যায় না।

## Step 5: k-mer ও mutation flag

Primary experiment-এর জন্য 4-mer:

```bat
python src\05_extract_features.py --k 4
```

Output:

```text
data/processed/features_kmer_k4.csv
data/processed/mutation_flags.csv
data/processed/features_combined_k4.csv
```

- `features_kmer_k4.csv`: দুই allele-এর normalized 4-mer frequency-এর average। এটি primary ML input।
- `mutation_flags.csv`: কোন known mutation কয়টি allele-এ আছে, dosage 0/1/2। এটি explanation এবং ablation study-এর জন্য।
- `features_combined_k4.csv`: k-mer + mutation flags। এটি secondary comparison।

Known ClinVar pathogenic flag দিয়ে label তৈরি করা হয়েছে। তাই একই flag সরাসরি input দিলে label leakage হতে পারে। এজন্য primary result k-mer-only file থেকে report করবে এবং mutation-only/combined result আলাদা ablation হিসেবে দেখাবে।

## Step 6: Baseline ML training

```bat
python src\06_train_baseline.py --features kmer --k 4
```

Secondary ablation:

```bat
python src\06_train_baseline.py --features combined --k 4
```

Random Forest এবং RBF-SVM train হবে। Output `models/` এবং `results/` folder-এ থাকবে। Accuracy-এর সঙ্গে macro-F1, class-wise precision/recall এবং confusion matrix report করবে।

ZIP-এ included example run-এর trained model ও result-ও দেওয়া আছে। Exact score এবং interpretation `docs/EXAMPLE_RESULTS_BN.md`-এ আছে। নিজের environment-এ run করলে package/library version-এর কারণে সামান্য numerical difference হতে পারে।

## Included data inspect করা

```bat
python src\07_inspect_dataset.py
```

এটি কোনো file পরিবর্তন না করে class/split count, genotype mapping, sequence length এবং feature dimensions print করবে। পুরো ছোট demo এক command-এ চালাতে:

```bat
run_quick_demo.bat
```

## Dataset size কীভাবে বুঝবে

প্রতিটি script শেষে count print করে। CSV manually দেখতে VS Code-এ Rainbow CSV extension ব্যবহার করতে পারো। Python দিয়ে count:

```bat
python -c "import pandas as pd; d=pd.read_csv('data/processed/hbb_three_class_genotypes.csv'); print(d['class_label'].value_counts()); print('Total:',len(d))"
```

`--n-per-class 300` হলে expected final size:

```text
Normal      300
Carrier     300
Affected    300
Total       900
```

এগুলো 900 independent clinical patients নয়। Normal backgrounds 1000 Genomes-derived; Carrier/Affected rows variant-informed simulated genotypes। Report-এ unique ClinVar variant count, unique population background count এবং generated genotype count আলাদা করে লিখবে।

## Report-এ ব্যবহারযোগ্য dataset statement

> The study used the GRCh38 HBB genomic interval (NC_000011.10:5225464-5227071) as the reference. Population haplotypes were obtained from phased 1000 Genomes genotypes, while beta-thalassemia-associated pathogenic and likely pathogenic variants were curated from ClinVar. Monoallelic and biallelic variant-informed genotypes were computationally constructed for method development. Therefore, the resulting dataset is a simulated genotype dataset and is not a clinically validated patient cohort.

## Scientific limitations

- `Normal` এখানে no-known-P/LP-HBB genotype proxy; সম্পূর্ণ medically healthy প্রমাণ নয়।
- `Affected` এখানে biallelic P/LP genotype proxy; disease severity শুধু HBB sequence দিয়ে সবসময় নির্ভুলভাবে বলা যায় না।
- Large deletions/CNV, complex structural variants এবং modifier genes এই starter pipeline-এর scope-এর বাইরে।
- এটি clinical diagnosis বা patient treatment-এর জন্য ব্যবহার করা যাবে না।
- Supervisor real patient dataset চাইলে এই simulated workflow final evidence হিসেবে যথেষ্ট নয়; ethics-approved genotype + phenotype cohort লাগবে।
