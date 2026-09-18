# সম্ভাব্য viva প্রশ্ন ও উত্তর

## Biology ও data

**Q: HBB chromosome 11-এ কেন?**  
A: মানুষের HBB gene autosomal chromosome 11-এর short arm, 11p15.4 region-এ অবস্থিত। Beta-thalassemia HBB mutation-এর কারণে; এটি sex-linked disease নয়।

**Q: Chromosome 11 মানে 11 নম্বর pair?**  
A: আমরা chromosome type 11 বলি। Diploid মানুষের এই type-এর দুই homologous copy থাকে—একটি মা, একটি বাবা থেকে। তাই sequence model-এ দুই allele রাখা হয়েছে।

**Q: `NM_000518.5` কী? DNA project-এ mRNA কেন?**  
A: এটি HBB RefSeq transcript accession। Coding/HGVS mutation name standard করতে ব্যবহৃত হয়। Model input genomic DNA interval `NC_000011.10` থেকে; transcript sequence input নয়।

**Q: Disease dataset আর mutation dataset কি আলাদা?**  
A: ClinVar এখানে disease-associated mutation catalogue, patient disease cohort নয়। Mutation information দিয়ে simulated carrier/affected genotype তৈরি হয়েছে। তাই তিনটি independent ready-made dataset merge করিনি; three sources-এর ভিন্ন role ছিল।

**Q: Normal কীভাবে চিনেছ?**  
A: 1000 Genomes background genotype-এর মধ্যে pipeline-matched ClinVar P/LP HBB variant না থাকলে no-known-P/LP proxy ধরা হয়েছে। “Confirmed medically normal” বলা হয়নি।

**Q: 900 মানে 900 patient?**  
A: না। 300 base population sample থেকে প্রতি sample-এর Normal, Carrier ও Affected-proxy row বানিয়ে 900 constructed genotype row হয়েছে।

**Q: Carrier আর Affected difference?**  
A: Carrier proxy-তে pathogenic allele dosage 1; Affected-genotype proxy-তে dosage 2। Clinical severity নির্ধারণের জন্য sequence ছাড়াও phenotype/modifier data লাগতে পারে।

## Preprocessing ও feature

**Q: k-mer কেন?**  
A: Variable-length/text DNA-কে fixed-length numeric vector বানায় এবং local sequence pattern ধরে। এতে standard ML model ব্যবহার করা যায়।

**Q: k=4 হলে 256 কেন?**  
A: প্রতিটি position-এ চারটি base-এর একটি হতে পারে, তাই `4 × 4 × 4 × 4 = 256`।

**Q: Raw count না frequency কেন?**  
A: Count-কে valid window দিয়ে ভাগ করলে feature sequence length-এর প্রভাব কমিয়ে comparable হয়।

**Q: দুই allele কীভাবে combine করেছ?**  
A: allele 1 ও allele 2-এর normalized k-mer frequency আলাদাভাবে বের করে element-wise average করেছি।

**Q: Mutation flag-এর use কী?**  
A: specific variant dosage ব্যাখ্যা এবং ablation study। Primary model-এ এটি ব্যবহার করিনি, কারণ label creation-এর একই তথ্য হওয়ায় leakage risk আছে।

**Q: `class_label` কেন দরকার?**  
A: Supervised learning-এ training-এর সময় model-কে correct target দেখাতে হয়। Test prediction-এর সময় label feature হিসেবে দেওয়া হয় না; evaluation-এর জন্য রাখা হয়।

## ML design

**Q: Random Forest ও SVM কেন?**  
A: Structured numeric feature-এর জন্য দুইটি standard baseline: Random Forest nonlinear tree ensemble, RBF-SVM scaled feature space-এ nonlinear boundary শেখে। Deep learning-এর আগে baseline দরকার।

**Q: Split কীভাবে?**  
A: 70% train, 15% validation, 15% test। মোট 630/135/135 row। একই base sample ও pathogenic variant set split-এর বাইরে যায় না।

**Q: Random split করলেই হতো না?**  
A: একই base sample-এর related rows বা একই inserted variant train ও test-এ গেলে model memorization করে inflated score দিতে পারে। Grouped construction এই leakage কমায়।

**Q: Validation set code-এ ব্যবহার হয়নি কেন?**  
A: এটি future hyperparameter tuning-এর জন্য reserved। Included script একটি fixed-parameter baseline; tuning করলে শুধু validation ব্যবহার করে শেষে একবার test evaluate করতে হবে।

**Q: Accuracy যথেষ্ট?**  
A: না। Balanced data হলেও macro-F1, per-class precision/recall এবং confusion matrix report করি।

**Q: Clinical use করা যাবে?**  
A: না। এটি simulated genotype method-development dataset, real phenotype-validated patient cohort নয়। External clinical validation এবং ethics-approved data ছাড়া diagnostic claim করা যাবে না।

## Code ownership নিয়ে সৎ উত্তর

**Q: সব code তুমি একা লিখেছ?**  
A: “আমি AI-assisted development ব্যবহার করেছি, কিন্তু source selection, assumptions, code execution এবং outputs যাচাই করেছি। Reproducibility-এর জন্য complete scripts, manifest এবং limitations জমা দিয়েছি।” নিজের university-এর academic-integrity policy অনুযায়ী AI disclosure করবে।
