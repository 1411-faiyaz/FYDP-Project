# Included example result কীভাবে ব্যাখ্যা করবে

Included 900-row dataset-এর predefined test set (135 rows) ব্যবহার করে 4-mer-only baseline আবার চালিয়ে পাওয়া result:

| Model | Test accuracy | Test macro-F1 |
|---|---:|---:|
| Random Forest | 0.6000 | 0.5707 |
| RBF-SVM | 0.7926 | 0.7826 |

RBF-SVM-এর class recall:

| Class | Recall |
|---|---:|
| Normal | 0.9556 |
| Carrier | 0.4889 |
| Affected | 0.9333 |

এখানে Carrier recall সবচেয়ে কম। কারণ দুই allele-এর k-mer vector average করলে একটি mutated allele-এর signal normal allele-এর signal-এর সঙ্গে মিশে যায়; Carrier ও অন্য class-এর boundary কঠিন হতে পারে। এটি result লুকানোর বিষয় নয়—future improvement-এর research direction। উদাহরণ:

- allele 1 ও allele 2-এর feature average না করে concatenate করা
- multiple k values compare করা
- sequence model বা biologically focused variant feature পরীক্ষা করা
- real phenotype-validated external cohort দিয়ে validation

`results/` folder-এ exact JSON metrics এবং confusion-matrix PNG আছে। `models/` folder-এ trained JOBLIB model আছে। এই score simulated genotype proxy dataset-এর জন্য; clinical accuracy হিসেবে দাবি করা যাবে না।

## নতুন validation-based comparison (2026-09-18)

`src/08_compare_models.py` দিয়ে included 900-row data-তে ছয় candidate তুলনা করা হয়েছে। Validation macro-F1 দিয়ে প্রতিটি representation-এর জয়ী model নির্বাচন করে train+validation-এ refit করা হয়। Python 3.14, scikit-learn 1.9.1 ব্যবহার হয়েছে।

| Feature | Validation winner | Validation macro-F1 | Final test macro-F1 | Test accuracy | Test Carrier recall |
|---|---|---:|---:|---:|---:|
| Mean 4-mer | Random Forest | 0.7453 | 0.6066 | 0.6296 | 0.5333 |
| Allele-pair 4-mer | Extra Trees | 0.8253 | 0.6899 | 0.7037 | 0.6667 |

দুটি representation-ই exploration হিসেবে test-এ দেখা হয়েছে; তাই এই table দেখে আরেকটি winner নির্বাচন করে সেটিকে untouched final test result বলা যাবে না। আগের fixed RBF-SVM baseline-এর test macro-F1 0.7826, যা এই দুটি নতুন selected model-এর test score-এর চেয়ে বেশি। Validation থেকে test-এ score কমেছে; variant/background distribution বদলালে performance স্থিতিশীল নাও থাকতে পারে। `results/comparison_*_validation.csv` ও `results/comparison_*_test.json`-এ exact সংখ্যা আছে। Clinical accuracy দাবি করা যাবে না।
