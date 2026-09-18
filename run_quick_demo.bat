@echo off
setlocal

echo [1/4] Inspecting the included 900-row dataset...
python src\07_inspect_dataset.py
if errorlevel 1 goto :error

echo [2/4] Rebuilding normalized 4-mer features...
python src\05_extract_features.py --k 4
if errorlevel 1 goto :error

echo [3/4] Training k-mer-only baseline models...
python src\06_train_baseline.py --features kmer --k 4
if errorlevel 1 goto :error

echo [4/4] Comparing additional models on validation data...
python src\08_compare_models.py --features kmer --k 4 --validation-only
if errorlevel 1 goto :error

echo Demo complete. Check the results and models folders.
exit /b 0

:error
echo A step failed. Read the error above and confirm that the virtual environment is active.
exit /b 1
