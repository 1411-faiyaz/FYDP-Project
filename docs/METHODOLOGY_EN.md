# Methodology draft for the report

## Data sources and reference coordinate system

The analysis used the GRCh38 HBB genomic interval on chromosome 11 (`NC_000011.10:5225464-5227071`), represented on the genomic plus strand to remain consistent with VCF REF and ALT alleles. The HBB RefSeq transcript `NM_000518.5` was used only to interpret transcript-based HGVS names; transcript sequence was not used as the model input. Population genotypes were obtained from the phased GRCh38 1000 Genomes VCF, and clinical variant interpretations were obtained from the NCBI ClinVar variant summary.

## ClinVar filtering

ClinVar records were restricted to GRCh38, gene symbol HBB, chromosome 11, and the selected HBB interval. Beta-thalassemia-associated records with an exact clinical-significance category of Pathogenic, Likely pathogenic, or Pathogenic/Likely pathogenic were retained. Only records with explicit VCF-compatible DNA REF and ALT alleles were used. Records were deduplicated using the tuple of genomic position, REF, and ALT. The collection used for the included example contained 279 unique filtered pathogenic or likely pathogenic variants, of which 271 were compatible with the downloaded genomic reference.

## Construction of the three-class genotype dataset

Phased 1000 Genomes genotypes were used to reconstruct two HBB haplotypes for each background sample. Samples carrying a pipeline-matched ClinVar pathogenic or likely pathogenic HBB variant were excluded from the background pool. For each selected background sample, a no-known-pathogenic genotype was retained as a Normal proxy. A monoallelic pathogenic genotype was computationally constructed as a Carrier proxy, and a biallelic pathogenic genotype was computationally constructed as an Affected-genotype proxy. REF matching and interval-overlap checks were performed before applying each variant. The resulting balanced example dataset contained 300 rows per class and 900 rows in total.

The generated labels are genotype proxies rather than confirmed clinical diagnoses. In particular, the Normal class indicates the absence of a matched known ClinVar pathogenic or likely pathogenic HBB variant within this pipeline, not proof of overall medical health. Carrier and Affected-genotype rows are variant-informed simulations and do not constitute a clinical patient cohort.

## Data partitioning and leakage control

Background samples and pathogenic variants were partitioned into training, validation, and test pools using a fixed random seed of 42. All three rows derived from a base population sample were confined to one split. Pathogenic variant sets were also kept disjoint across splits. The final partition contained 630 training, 135 validation, and 135 test rows. This design reduced direct leakage through shared backgrounds and repeated pathogenic variants.

## Feature extraction

For the primary experiment, each allele sequence was converted into normalized 4-mer frequencies. All 256 possible DNA 4-mers were enumerated, and the count of each valid 4-mer was divided by the number of valid sliding windows in that allele. The two allele-level frequency vectors were averaged to obtain one fixed-length feature vector per sample. Mutation-dosage flags were also generated, with values of 0, 1, or 2 indicating the number of alleles carrying a specific constructed pathogenic variant. Because these flags are closely related to the procedure used to create the class labels, they were reserved for interpretation and secondary ablation analysis rather than the primary model.

## Baseline models and evaluation

Random Forest and radial-basis-function Support Vector Machine classifiers were used as conventional supervised-learning baselines. The SVM pipeline included standard feature scaling. Models were trained using the predefined training partition and evaluated on the held-out test partition. Performance was summarized using accuracy, macro-averaged F1 score, class-wise precision and recall, and confusion matrices. The validation partition was reserved for future hyperparameter selection; the included baseline used fixed hyperparameters.

## Reproducibility and scope

All download, filtering, genotype-construction, feature-extraction, inspection, and training code is provided as numbered Python scripts. Source counts, coordinates, and the collection date are recorded in `DATA_MANIFEST.md`. This workflow is intended for coursework and methodological development. It is not suitable for clinical diagnosis without an independently collected, ethics-approved, phenotype-validated cohort and external validation.
