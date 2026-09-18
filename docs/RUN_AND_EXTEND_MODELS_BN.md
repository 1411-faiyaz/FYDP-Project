# Code দিয়ে পুরো কাজ চালানো এবং নতুন model train করা

এই project-এ preprocessing ও baseline training code আগে থেকেই `src/`-এ ছিল। এই guide-এ কোন command কোন code চালায়, কী output দেয়, এবং নতুন model কীভাবে তুলনা করবে তা এক জায়গায় দেওয়া হলো। Included ফলাফল 2026-09-03/07-এর example run; নতুন run করলে source version ও library version ভেদে ফল বদলাতে পারে।

## 1. প্রথমে environment তৈরি

VS Code-এ project folder খুলে Windows Command Prompt-এ:

```bat
py -3.11 -m venv .venv
.venv\Scripts\python.exe -m pip install -r requirements.txt
```

Python 3.12 থাকলে `py -3.12` ব্যবহার করা যায়। `pysam`-সহ পুরো pipeline চালাতে 3.11/3.12 ব্যবহার করো। পরে সব command-এ `.venv\Scripts\python.exe` ব্যবহার করলে activation প্রয়োজন হয় না। যদি শুধু included CSV থেকে feature/model run করো, download বা `pysam` দরকার হয় না; সেক্ষেত্রে বর্তমান Python দিয়ে virtual environment বানিয়ে `requirements-ml.txt` install করলেই হবে:

```bat
py -m venv .venv
.venv\Scripts\python.exe -m pip install -r requirements-ml.txt
```

## 2. আগে থেকে থাকা data দিয়ে offline পুনরায় চালাও

এই command-গুলো internet ছাড়া included 900-row CSV থেকে feature ও model output তৈরি করে:

```bat
.venv\Scripts\python.exe src\07_inspect_dataset.py
.venv\Scripts\python.exe src\05_extract_features.py --k 4
.venv\Scripts\python.exe src\06_train_baseline.py --features kmer --k 4
.venv\Scripts\python.exe src\08_compare_models.py --features kmer --k 4 --validation-only
.venv\Scripts\python.exe src\08_compare_models.py --features allele_pair --k 4 --validation-only
```

আগের run-এর একই নামের feature/model/result file overwrite হবে। তাই original example artifact রাখতে চাইলে folder copy করে experiment করো। `07_inspect_dataset.py` শুধু পড়ে; অন্য চারটি command output লেখে।

## 3. শুরু থেকে preprocessing: কোন code কী করেছে

| ধাপ | Code ও command | Input | কাজ | Output |
|---|---|---|---|---|
| 1 | `src/01_download_reference.py` | NCBI GRCh38 | HBB genomic interval FASTA download, 1,608 bp ও ACGT validation | `data/raw/hbb_reference_grch38_plus.fasta` |
| 2 | `src/02_download_and_filter_clinvar.py` | ClinVar `variant_summary.txt.gz` | GRCh38, HBB, chr11, region, significance ও phenotype filter; `(position, REF, ALT)` দিয়ে deduplicate | `data/interim/hbb_beta_thal_pathogenic_clinvar.csv`, `hbb_benign_clinvar.csv` |
| 3 | `src/03_subset_1000g.py` | indexed 1000 Genomes chr11 VCF | শুধু HBB interval fetch, local bgzip VCF ও index লেখা | `data/raw/hbb_1000g_grch38.vcf.gz` এবং `.tbi` |
| 4 | `src/04_build_balanced_dataset.py` | FASTA + ClinVar CSV + VCF | phased genotype দিয়ে দুই haplotype reconstruct; REF check; compatible P/LP variant বসিয়ে তিন proxy class; base sample ও mutation আলাদা split | `data/processed/hbb_three_class_genotypes.csv` |
| 5 | `src/05_extract_features.py` | দুই allele sequence | sliding window দিয়ে normalized k-mer; দুই allele-এর mean ও min/max vector; mutation dosage flag | `features_kmer_k4.csv`, `features_allele_pair_k4.csv`, `mutation_flags.csv`, `features_combined_k4.csv` |
| 6 | `src/06_train_baseline.py` | k-mer feature CSV | fixed Random Forest ও scaled RBF-SVM train; test metrics | `models/`-এ JOBLIB, `results/`-এ JSON ও PNG |
| 7 | `src/08_compare_models.py` | k-mer অথবা allele-pair CSV | validation-এ candidate তুলনা; জয়ী model train+validation দিয়ে refit; একবার test | validation CSV, final test JSON/PNG ও JOBLIB |

Full pipeline চালাতে নিচের command-গুলো এই ক্রমে দাও। প্রথম তিন ধাপে internet এবং ClinVar-এর বড় download লাগতে পারে:

```bat
.venv\Scripts\python.exe src\01_download_reference.py
.venv\Scripts\python.exe src\02_download_and_filter_clinvar.py
.venv\Scripts\python.exe src\03_subset_1000g.py
.venv\Scripts\python.exe src\04_build_balanced_dataset.py --n-per-class 300 --seed 42
.venv\Scripts\python.exe src\05_extract_features.py --k 4
.venv\Scripts\python.exe src\07_inspect_dataset.py
.venv\Scripts\python.exe src\06_train_baseline.py --features kmer --k 4
```

`run_full_pipeline.bat` প্রথম ছয় ধাপ ও baseline এক command-এ চালায়। Included source/output-এর exact counts এবং date `DATA_MANIFEST.md`-এ আছে। কোনো নতুন run-এ count বদলালে নতুন date, source এবং count report-এ লিখবে।

### Label তৈরির logic

`04_build_balanced_dataset.py`-এর `apply_variants()` genomic position থেকে sequence offset নেয়, reference allele মিলিয়ে দেখে, এবং indel ঠিক রাখতে ডান থেকে বামে variant বসায়। প্রতি selected population background থেকে একই split-এ তিন row হয়:

| Row | P/LP variant dosage | `class_label` | অর্থ |
|---|---:|---|---|
| `_N` | 0 | `Normal` | no-known-P/LP-HBB proxy |
| `_C` | 1 | `Carrier` | computational monoallelic proxy |
| `_A` | 2 | `Affected` | computational biallelic proxy |

তাই 900 row মানে 900 independent patient নয়; included dataset-এ 300 base population background থেকে তিনটি করে row। `synthetic_genotype` ও `label_basis` column provenance রাখে। Normal মানে medically confirmed healthy নয়, আর Carrier/Affected confirmed diagnosis নয়।

### Feature বের করার logic

`05_extract_features.py`-তে 4-mer হলে `4^4 = 256` possible DNA pattern। প্রতিটি allele-এ 4 base-এর sliding window count করে valid window সংখ্যা দিয়ে ভাগ করা হয়। `KMER_ACGT`-এর মত column দুই allele-এর frequency-এর গড়। নতুন `PAIR_MIN_ACGT` ও `PAIR_MAX_ACGT` একই pattern-এর দুই frequency-এর ছোট ও বড় মান রাখে। ফলে allele order বদলালেও value বদলায় না, কিন্তু দুই allele-এর পার্থক্য কিছুটা থাকে। Pair file-এ 512 numeric column আছে।

`mutation_flags.csv` ও `features_combined_k4.csv` label বানানোর P/LP dosage থেকে এসেছে। এগুলো দিয়ে বেশি score পেলেও primary result হিসেবে দেখানো ঠিক নয়: direct label leakage-এর ঝুঁকি আছে। Primary comparison-এ `kmer` বা `allele_pair` ব্যবহার করো।

## 4. Baseline আর নতুন model comparison-এর তফাৎ

`06_train_baseline.py` আগের example ফল পুনরায় বানায়: fixed Random Forest এবং RBF-SVM training split-এ fit করে test split-এ মাপে। Validation split সেখানে ব্যবহার হয় না। Published example result ব্যাখ্যা `docs/EXAMPLE_RESULTS_BN.md`-এ আছে।

নতুন `08_compare_models.py` ছয় candidate নেয়: Dummy majority-class control, Logistic Regression, k-nearest neighbors (KNN), Random Forest, Extra Trees, RBF-SVM। Scaling প্রয়োজন হয় এমন model-এ `StandardScaler` pipeline-এর ভেতরে আছে; training data থেকেই scaler fit হয়। প্রতিটি model **শুধু train split-এ** fit করে validation split-এর macro-F1 দেখা হয়। সর্বোচ্চ validation macro-F1 (tie হলে accuracy, তারপর নাম) পাওয়া model-টি train+validation data-তে আবার fit হয়। Test split তখন **শুধু সেই জয়ী model-এর জন্য একবার** ব্যবহার হয়।

সঠিক comparison order: প্রথমে দুই representation-ই শুধু validation-এ দেখো, তার ফল দেখে representation ও candidate list ঠিক করো। `--validation-only` দিলে script test score হিসাব করে না এবং final model-ও লেখে না। তারপর বেছে নেওয়া **একটি** representation-এ flag ছাড়া final run করো:

```bat
.venv\Scripts\python.exe src\08_compare_models.py --features kmer --k 4 --validation-only
.venv\Scripts\python.exe src\08_compare_models.py --features allele_pair --k 4 --validation-only
.venv\Scripts\python.exe src\08_compare_models.py --features allele_pair --k 4
```

নির্দিষ্ট candidate নিয়ে experiment করতে:

```bat
.venv\Scripts\python.exe src\08_compare_models.py --features kmer --k 4 --validation-only --models logistic_regression extra_trees svm_rbf
```

Output:

- `results/comparison_kmer_k4_validation.csv`: সব candidate-এর validation macro-F1 ও accuracy; model নির্বাচন এখানেই।
- `results/comparison_kmer_k4_test.json`: জয়ী model-এর test classification report, confusion matrix, library version, split size ও feature-file SHA-256।
- `results/comparison_kmer_k4_confusion_matrix.png`: true label বনাম prediction।
- `models/best_kmer_k4.joblib`: fitted estimator ও feature-column order।

`allele_pair` experiment-এ একই নামের মধ্যে `allele_pair` থাকবে। Accuracy ছাড়াও macro-F1 এবং বিশেষ করে Carrier recall দেখবে। দুই representation-এর test result বারবার দেখে নতুন model/parameter বেছে নিলে test set আর নিরপেক্ষ থাকে না; final report-এর আগে development choice validation-এ ঠিক করো।

## 5. নিজের নতুন model যোগ করা

`src/08_compare_models.py`-তে তিনটি ছোট পরিবর্তন করলেই হয়। উদাহরণ হিসেবে `GradientBoostingClassifier`:

```python
from sklearn.ensemble import GradientBoostingClassifier

MODEL_NAMES = [
    "dummy", "logistic_regression", "knn", "random_forest",
    "extra_trees", "svm_rbf", "gradient_boosting",
]

def build_models(seed: int) -> dict:
    return {
        # বিদ্যমান model-গুলো এখানে থাকবে
        "gradient_boosting": GradientBoostingClassifier(random_state=seed),
    }
```

উপরের snippet-এ পুরো `build_models()` replace করবে না; existing dictionary-তে শুধু নতুন entry যোগ করবে। তারপর চালাও:

```bat
.venv\Scripts\python.exe src\08_compare_models.py --features kmer --k 4 --validation-only --models dummy gradient_boosting
```

নতুন model-এর parameter পরীক্ষা করতে চাইলে একাধিক **candidate name** দিয়ে আলাদা configuration dictionary-তে যোগ করো, validation score দেখে আগে নির্বাচন করো। Test দেখে parameter বদলাবে না। উদাহরণ: `GradientBoostingClassifier(n_estimators=100, learning_rate=0.05, random_state=seed)`। Final chosen model-এর test score একবার report করবে। Deep learning model করতে চাইলে input shape ও training code আলাদা লাগবে; আগে এই reproducible conventional baselines শেষ করো।

## 6. পরের research কাজের বাস্তব ক্রম

1. Included data inspect করে class/split/variant leakage check করো।
2. k=3, 4, 5 feature আলাদা তৈরি করো; প্রতি k-এর candidate selection validation-এ করো। `k=5` হলে 1,024 mean বা 2,048 pair feature হয় এবং memory/time বাড়ে।
3. Carrier recall ও confusion matrix দেখে mean বনাম allele-pair representation তুলনা করো।
4. যেটি আগে ঠিক করেছ সেটির held-out test result report করো। Raw validation score, candidate list, random seed, feature count ও package version রাখো।
5. এরপর real phenotype-validated, ethics-approved independent cohort পেলে external validation করো। এই simulated proxy dataset-এর score clinical performance নয়। Large deletion/CNV ও modifier-gene effect এই pipeline-এ ধরা হয় না।

কোনো saved model দিয়ে prediction করার সময় একই feature schema ও column order লাগবে:

```python
import joblib
import pandas as pd

bundle = joblib.load("models/best_kmer_k4.joblib")  # নিজের তৈরি trusted file
rows = pd.read_csv("data/processed/features_kmer_k4.csv")
predictions = bundle["estimator"].predict(rows[bundle["feature_columns"]])
print(predictions[:5])
```

এই example একই feature CSV-এর row predict করে code দেখায়; এটি independent clinical validation নয়। `joblib` file শুধু trusted source থেকে load করবে।
