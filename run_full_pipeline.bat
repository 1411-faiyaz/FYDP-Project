@echo off
setlocal

echo WARNING: The ClinVar step downloads a large current bulk file.
echo The included example data can be demonstrated with run_quick_demo.bat instead.
pause

python src\01_download_reference.py
if errorlevel 1 goto :error
python src\02_download_and_filter_clinvar.py
if errorlevel 1 goto :error
python src\03_subset_1000g.py
if errorlevel 1 goto :error
python src\04_build_balanced_dataset.py --n-per-class 300 --seed 42
if errorlevel 1 goto :error
python src\05_extract_features.py --k 4
if errorlevel 1 goto :error
python src\07_inspect_dataset.py
if errorlevel 1 goto :error
python src\06_train_baseline.py --features kmer --k 4
if errorlevel 1 goto :error

echo Full pipeline complete.
exit /b 0

:error
echo Pipeline stopped because a step failed. Read the error above.
exit /b 1
