import argparse
import json

import joblib
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC

from config import MODELS_DIR, PROCESSED_DIR, RESULTS_DIR, ensure_directories


META_COLUMNS = {"sample_id", "base_sample_id", "split", "class_label"}
LABEL_ORDER = ["Normal", "Carrier", "Affected"]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--k", type=int, default=4)
    parser.add_argument(
        "--features",
        choices=["kmer", "combined"],
        default="kmer",
        help="Use kmer for the primary experiment; combined is only an ablation comparison.",
    )
    args = parser.parse_args()
    ensure_directories()

    filename = (
        f"features_kmer_k{args.k}.csv"
        if args.features == "kmer"
        else f"features_combined_k{args.k}.csv"
    )
    path = PROCESSED_DIR / filename
    data = pd.read_csv(path)
    feature_columns = [column for column in data.columns if column not in META_COLUMNS]

    train = data[data["split"] == "train"]
    validation = data[data["split"] == "validation"]
    test = data[data["split"] == "test"]
    if train.empty or validation.empty or test.empty:
        raise RuntimeError("Train, validation and test splits must all be non-empty")

    x_train = train[feature_columns]
    y_train = train["class_label"]
    x_test = test[feature_columns]
    y_test = test["class_label"]

    models = {
        "random_forest": RandomForestClassifier(
            n_estimators=300,
            class_weight="balanced",
            random_state=42,
            n_jobs=-1,
        ),
        "svm_rbf": make_pipeline(
            StandardScaler(),
            SVC(C=2.0, kernel="rbf", class_weight="balanced"),
        ),
    }

    for model_name, model in models.items():
        model.fit(x_train, y_train)
        prediction = model.predict(x_test)
        report = classification_report(
            y_test,
            prediction,
            labels=LABEL_ORDER,
            output_dict=True,
            zero_division=0,
        )
        matrix = confusion_matrix(y_test, prediction, labels=LABEL_ORDER)

        result_stem = f"{model_name}_{args.features}_k{args.k}"
        (RESULTS_DIR / f"{result_stem}_metrics.json").write_text(
            json.dumps(report, indent=2),
            encoding="utf-8",
        )
        joblib.dump(model, MODELS_DIR / f"{result_stem}.joblib")

        plt.figure(figsize=(6, 5))
        sns.heatmap(
            matrix,
            annot=True,
            fmt="d",
            cmap="Blues",
            xticklabels=LABEL_ORDER,
            yticklabels=LABEL_ORDER,
        )
        plt.xlabel("Predicted")
        plt.ylabel("Actual")
        plt.tight_layout()
        plt.savefig(RESULTS_DIR / f"{result_stem}_confusion_matrix.png", dpi=180)
        plt.close()

        print(
            f"{model_name}: test macro-F1="
            f"{report['macro avg']['f1-score']:.4f}"
        )


if __name__ == "__main__":
    main()
