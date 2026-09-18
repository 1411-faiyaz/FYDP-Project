"""Select a model on validation data, then evaluate the winner once on test data."""

import argparse
import hashlib
import json
import sys

import joblib
import numpy as np
import pandas as pd
import sklearn
from sklearn.base import clone
from sklearn.dummy import DummyClassifier
from sklearn.ensemble import ExtraTreesClassifier, RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    ConfusionMatrixDisplay,
    accuracy_score,
    classification_report,
    confusion_matrix,
    f1_score,
)
from sklearn.neighbors import KNeighborsClassifier
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC

from config import MODELS_DIR, PROCESSED_DIR, RESULTS_DIR, ensure_directories


META_COLUMNS = {"sample_id", "base_sample_id", "split", "class_label"}
LABEL_ORDER = ["Normal", "Carrier", "Affected"]
MODEL_NAMES = [
    "dummy",
    "logistic_regression",
    "knn",
    "random_forest",
    "extra_trees",
    "svm_rbf",
]


def build_models(seed: int) -> dict:
    """All preprocessing that learns from data stays inside each model pipeline."""
    return {
        "dummy": DummyClassifier(strategy="most_frequent"),
        "logistic_regression": make_pipeline(
            StandardScaler(),
            LogisticRegression(max_iter=3000, class_weight="balanced", random_state=seed),
        ),
        "knn": make_pipeline(StandardScaler(), KNeighborsClassifier(n_neighbors=5)),
        "random_forest": RandomForestClassifier(
            n_estimators=300, class_weight="balanced", random_state=seed, n_jobs=-1
        ),
        "extra_trees": ExtraTreesClassifier(
            n_estimators=300, class_weight="balanced", random_state=seed, n_jobs=-1
        ),
        "svm_rbf": make_pipeline(
            StandardScaler(), SVC(C=2.0, kernel="rbf", class_weight="balanced")
        ),
    }


def load_features(path, representation: str):
    if not path.exists():
        raise FileNotFoundError(f"Missing {path}. Run src/05_extract_features.py first.")
    data = pd.read_csv(path)
    if not META_COLUMNS.issubset(data.columns):
        raise ValueError(f"Required metadata columns are missing from {path}")
    if data.empty or data["sample_id"].isna().any() or data["sample_id"].duplicated().any():
        raise ValueError("Feature rows need unique, nonempty sample IDs")
    if set(data["split"]) != {"train", "validation", "test"}:
        raise ValueError("Expected nonempty train, validation, and test splits")
    if set(data["class_label"]) != set(LABEL_ORDER):
        raise ValueError(f"Expected exactly these labels: {LABEL_ORDER}")
    if (data.groupby("base_sample_id")["split"].nunique() > 1).any():
        raise ValueError("A base sample occurs in more than one split")

    feature_columns = [name for name in data.columns if name not in META_COLUMNS]
    prefix = "KMER_" if representation == "kmer" else ("PAIR_MIN_", "PAIR_MAX_")
    if not feature_columns or not all(name.startswith(prefix) for name in feature_columns):
        raise ValueError("Unexpected feature columns; use a sequence-only feature file")
    try:
        values = data[feature_columns].to_numpy(dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError("All model features must be numeric") from exc
    if not np.isfinite(values).all():
        raise ValueError("Feature file has missing or nonfinite values")
    return data, feature_columns


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--k", type=int, default=4, help="k value used in feature extraction")
    parser.add_argument("--features", choices=["kmer", "allele_pair"], default="kmer")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--validation-only", action="store_true",
        help="Save validation scores without fitting a final model or evaluating test.",
    )
    parser.add_argument(
        "--models", nargs="+", choices=MODEL_NAMES, default=MODEL_NAMES,
        help="Candidate models to compare on validation data",
    )
    args = parser.parse_args()
    ensure_directories()

    feature_path = PROCESSED_DIR / f"features_{args.features}_k{args.k}.csv"
    data, feature_columns = load_features(feature_path, args.features)
    train = data.loc[data["split"] == "train"]
    validation = data.loc[data["split"] == "validation"]
    models = build_models(args.seed)
    leaderboard = []

    for name in dict.fromkeys(args.models):
        model = clone(models[name])
        model.fit(train[feature_columns], train["class_label"])
        predicted = model.predict(validation[feature_columns])
        leaderboard.append({
            "model": name,
            "validation_macro_f1": f1_score(
                validation["class_label"], predicted, labels=LABEL_ORDER,
                average="macro", zero_division=0,
            ),
            "validation_accuracy": accuracy_score(validation["class_label"], predicted),
        })

    leaderboard.sort(
        key=lambda row: (-row["validation_macro_f1"], -row["validation_accuracy"], row["model"])
    )
    stem = f"comparison_{args.features}_k{args.k}"
    validation_path = RESULTS_DIR / f"{stem}_validation.csv"
    pd.DataFrame(leaderboard).to_csv(validation_path, index=False)
    winner_name = leaderboard[0]["model"]
    if args.validation_only:
        print(pd.DataFrame(leaderboard).to_string(index=False))
        print(f"Validation winner: {winner_name}")
        print(f"Validation scores: {validation_path}")
        print("Test split was not evaluated.")
        return

    # The winner is fixed using validation data before any test prediction is made.
    final_model = clone(models[winner_name])
    train_and_validation = data.loc[data["split"].isin(["train", "validation"])]
    final_model.fit(train_and_validation[feature_columns], train_and_validation["class_label"])
    test = data.loc[data["split"] == "test"]
    test_prediction = final_model.predict(test[feature_columns])
    test_matrix = confusion_matrix(test["class_label"], test_prediction, labels=LABEL_ORDER)
    test_report = classification_report(
        test["class_label"], test_prediction, labels=LABEL_ORDER,
        output_dict=True, zero_division=0,
    )

    feature_sha256 = hashlib.sha256(feature_path.read_bytes()).hexdigest()
    result = {
        "dataset": str(feature_path.relative_to(PROCESSED_DIR.parent.parent)),
        "feature_sha256": feature_sha256,
        "representation": args.features,
        "k": args.k,
        "seed": args.seed,
        "python_version": sys.version.split()[0],
        "sklearn_version": sklearn.__version__,
        "label_order": LABEL_ORDER,
        "feature_count": len(feature_columns),
        "split_rows": {part: int((data["split"] == part).sum()) for part in ("train", "validation", "test")},
        "validation_leaderboard": leaderboard,
        "selected_model": winner_name,
        "final_fit_rows": len(train_and_validation),
        "test_report": test_report,
        "test_confusion_matrix": test_matrix.tolist(),
        "interpretation": "Variant-informed simulated genotype proxies; not clinical accuracy.",
    }
    metrics_path = RESULTS_DIR / f"{stem}_test.json"
    metrics_path.write_text(json.dumps(result, indent=2), encoding="utf-8")

    model_path = MODELS_DIR / f"best_{args.features}_k{args.k}.joblib"
    joblib.dump({
        "estimator": final_model,
        "feature_columns": feature_columns,
        "label_order": LABEL_ORDER,
        "feature_sha256": feature_sha256,
        "representation": args.features,
        "k": args.k,
    }, model_path)

    import matplotlib.pyplot as plt

    display = ConfusionMatrixDisplay(test_matrix, display_labels=LABEL_ORDER)
    display.plot(cmap="Blues", values_format="d", colorbar=False)
    plt.tight_layout()
    matrix_path = RESULTS_DIR / f"{stem}_confusion_matrix.png"
    plt.savefig(matrix_path, dpi=180)
    plt.close()

    print(pd.DataFrame(leaderboard).to_string(index=False))
    print(f"Selected on validation: {winner_name}")
    print(f"Final test macro-F1: {test_report['macro avg']['f1-score']:.4f}")
    print(f"Validation scores: {validation_path}")
    print(f"Test metrics: {metrics_path}")
    print(f"Trained model: {model_path}")
    print(f"Confusion matrix: {matrix_path}")


if __name__ == "__main__":
    main()
