# Ma'am-কে এক মিনিটে বলার script

“আমাদের কাজের target হলো HBB DNA sequence থেকে Normal, Carrier এবং Affected-genotype proxy classify করা। আমরা প্রথমে GRCh38 chromosome 11 থেকে 1,608 base-pair HBB genomic reference নিয়েছি। Population background genotype নিয়েছি phased 1000 Genomes VCF থেকে এবং beta-thalassemia-associated Pathogenic/Likely pathogenic variant নিয়েছি ClinVar থেকে। সব source একই GRCh38 coordinate-এ filter করেছি।

Confirmed open patient-level three-class cohort না পাওয়ায় আমরা data limitation লুকাইনি। No-known-pathogenic 1000 Genomes genotype-কে Normal proxy ধরেছি; এক allele-এ ClinVar pathogenic variant computationally বসিয়ে Carrier এবং দুই allele-এ বসিয়ে Affected-genotype proxy বানিয়েছি। মোট 300 করে 900 labelled row হয়েছে।

এরপর প্রত্যেক sample-এর দুই allele থেকে normalized 4-mer frequency বের করেছি। চারটি nucleotide-এর 4-mer combination 4 to the power 4, অর্থাৎ 256 feature। `class_label` target দিয়ে Random Forest ও RBF-SVM supervised baseline train করেছি। Split করার সময় একই base sample ও pathogenic variant train-test-এ overlap না করতে ব্যবস্থা করেছি। Result accuracy-এর পাশাপাশি macro-F1, class-wise precision/recall এবং confusion matrix দিয়ে evaluate করি। এই dataset method development-এর জন্য; clinical diagnosis claim করি না।”

## যদি ma'am জিজ্ঞেস করেন ‘তুমি code কোথায় লিখেছ?’

বলবে: “`src/` folder-এ numbered Python scripts আছে। `05_extract_features.py` sliding window দিয়ে k-mer বানায় এবং `06_train_baseline.py` model training/evaluation করে। `START_HERE_BN.md`-এ exact run command আছে।”

## যদি live demo চায়

```bat
python src\07_inspect_dataset.py
python src\05_extract_features.py --k 4
python src\06_train_baseline.py --features kmer --k 4
```
