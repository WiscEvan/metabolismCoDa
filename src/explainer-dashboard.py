#!/usr/bin/env python


import argparse
import pandas as pd
from explainerdashboard import ExplainerDashboard, ClassifierExplainer
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier


def read_kofam_matrix(matrix: str) -> pd.DataFrame:
    mat_df = pd.read_table(matrix, index_col="ko_num")
    sample_cols = [col for col in mat_df.columns if ".tsv" in col]
    mat_df = mat_df[sample_cols].T.fillna(0)
    mat_df.index.name = "filename"
    return mat_df


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--matrix", help="Path to kofamscan_results_matrix.tsv", required=True
    )
    parser.add_argument(
        "--ko-data",
        help="Path to metabolism_feature_analysis.tsv (downloaded from feature-analysis-app.py)",
        required=True,
    )
    parser.add_argument(
        "--factor-name",
        help="Factor to use for modeling feature analysis",
    )
    parser.add_argument(
        "--n-estimators",
        "-T",
        help="Number of trees to use for training RandomForestClassifier",
        default=50,
        type=int,
    )
    parser.add_argument(
        "--n-jobs",
        help="Parallelizes jobs using joblib. For now only used for calculating permutation importances.",
        default=None,
        type=int,
    )
    parser.add_argument(
        "--host",
        help="Host address to use for dashboard",
        type=str,
        default="0.0.0.0",
    )
    parser.add_argument(
        "--port",
        help="Port number to use for dashboard",
        type=str,
        default="8855",
    )
    args = parser.parse_args()
    ko_matrix = read_kofam_matrix(args.matrix)
    ko_data = pd.read_table(args.ko_data)
    if not args.factor_name:
        factor_names = ", ".join(ko_data.set_index("filename").columns.tolist())
        print("A --factor-name must be provided...")
        print(f"Factor names available: {factor_names}")
        exit()
    factor_values = ko_data.set_index("filename")[args.factor_name]
    X = ko_matrix.loc[factor_values.index].copy()
    le = LabelEncoder()
    # le.fit(factor_values)
    # y = le.transform(factor_values)
    y = le.fit_transform(factor_values)
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=42)

    model = RandomForestClassifier(
        criterion="gini",
        random_state=42,
        n_estimators=args.n_estimators,
    )
    model.fit(X_train, y_train)
    explainer = ClassifierExplainer(
        model,
        X_test,
        y_test,
        n_jobs=args.n_jobs,
    )
    dashboard = ExplainerDashboard(explainer)
    dashboard.run(host=args.host, port=args.port)


if __name__ == "__main__":
    main()
