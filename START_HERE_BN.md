# এখান থেকে শুরু করো

এই ZIP-এ শুধু তৈরি dataset নেই—raw source, processed data, প্রতিটি processing script, file-format explanation, code walkthrough এবং viva উত্তরও আছে। **সব command ও নতুন model train করার বিস্তারিত guide:** `docs/RUN_AND_EXTEND_MODELS_BN.md`।

## প্রথমে পাঁচটি জিনিস বুঝবে

1. **Reference DNA** হলো GRCh38-এর chromosome 11-এ HBB gene-এর 1,608 bp অঞ্চল। এটি একটি baseline sequence, কোনো patient নয়।
2. **1000 Genomes VCF** থেকে মানুষের population background genotype নেওয়া হয়েছে।
3. **ClinVar CSV** থেকে beta-thalassemia-associated Pathogenic/Likely pathogenic variant নেওয়া হয়েছে। ClinVar row হলো variant, patient নয়।
4. Population background-এর এক allele-এ pathogenic variant বসিয়ে `Carrier` এবং দুই allele-এ বসিয়ে `Affected` genotype proxy তৈরি করা হয়েছে। তাই এগুলো **simulated**, confirmed clinical patient নয়।
5. দুই allele sequence থেকে normalized 4-mer frequency বানিয়ে ML model train করা হয়েছে। `class_label` হলো supervised-learning target।

## VS Code-এ সবচেয়ে সহজ demo

ZIP extract করে folder-টি VS Code-এ open করো। তারপর `Terminal > New Terminal` থেকে Windows Command Prompt-এ চালাও:

```bat
py -m venv .venv
.venv\Scripts\activate.bat
python -m pip install --upgrade pip
pip install -r requirements.txt
run_quick_demo.bat
```

এই demo existing 900-row data inspect করবে, 4-mer feature আবার বানাবে, k-mer-only baseline model train করবে এবং অতিরিক্ত model-গুলো validation-এ তুলনা করবে। বড় ClinVar file আবার download করবে না।

কোনো `.bat` file চালাতে সমস্যা হলে একই command আলাদাভাবে চালাও:

```bat
python src\07_inspect_dataset.py
python src\05_extract_features.py --k 4
python src\06_train_baseline.py --features kmer --k 4
python src\08_compare_models.py --features kmer --k 4 --validation-only
```

## কোন file কোন কাজে

| আগে পড়বে | উদ্দেশ্য |
|---|---|
| `docs/WORKFLOW_EXPLANATION_BN.md` | শুরু থেকে শেষ পর্যন্ত পুরো গবেষণা flow |
| `docs/FILE_FORMAT_GUIDE_BN.md` | FASTA, VCF, CSV, GZ, TBI, PY, JSON, PNG কী |
| `docs/CODE_WALKTHROUGH_BN.md` | প্রতিটি Python script কী করে |
| `docs/RUN_AND_EXTEND_MODELS_BN.md` | preprocessing command, validation দিয়ে নতুন model নির্বাচন, output ও পরের experiment |
| `docs/COLUMN_DICTIONARY_BN.md` | dataset-এর প্রতিটি গুরুত্বপূর্ণ column-এর অর্থ |
| `docs/EXAMPLE_RESULTS_BN.md` | included baseline result কীভাবে ব্যাখ্যা করবে |
| `docs/ONE_MINUTE_PRESENTATION_BN.md` | ma'am-কে এক মিনিটে বলার script |
| `docs/VIVA_QA_BN.md` | সম্ভাব্য viva প্রশ্ন ও সৎ উত্তর |
| `docs/METHODOLOGY_EN.md` | report-এর জন্য English methodology draft |

## একটি বাক্যে project

> We constructed a reproducible, variant-informed simulated HBB genotype dataset from GRCh38, 1000 Genomes and ClinVar, extracted normalized 4-mer features from two alleles, and trained supervised ML baselines to distinguish Normal, Carrier and Affected-genotype proxy classes.

`simulated`, `proxy` এবং `not clinically validated`—এই কথাগুলো কখনো লুকাবে না।
