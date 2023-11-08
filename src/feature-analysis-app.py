#!/usr/bin/env python

import argparse
from dash.exceptions import PreventUpdate
from dash_extensions.enrich import (
    DashProxy,
    ServersideOutputTransform,
    Serverside,
    State,
    html,
    Input,
    Output,
    dcc,
)
import dash_ag_grid as dag
import dash_mantine_components as dmc
from dash_iconify import DashIconify

from functools import partial
from enum import auto
from strenum import StrEnum
import plotly.express as px

import pandas as pd
import numpy as np

from typing import Dict, List, Optional, Tuple

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split

from plotly import graph_objects as go


class Ids(StrEnum):
    ko_data_store = auto()
    ko_scatterplot = auto()
    factor_name_select = auto()
    feature_importance_btn = auto()
    feature_importance_barchart = auto()
    factor_name_text_input = auto()
    factor_value_text_input = auto()
    factor_unclassified_text_input = auto()
    feature_importance_data_store = auto()
    feature_importance_clf_store = auto()
    forest_criterion = auto()
    feature_data_download = auto()
    feature_data_download_button = auto()
    ko_data_download = auto()
    ko_data_download_button = auto()
    clf_score_badge = auto()
    feature_table = auto()


def read_kofam_matrix(matrix: str) -> pd.DataFrame:
    mat_df = pd.read_table(matrix, index_col="ko_num")
    sample_cols = [col for col in mat_df.columns if ".tsv" in col]
    mat_df = mat_df[sample_cols].T.fillna(0)
    mat_df.index.name = "filename"
    return mat_df


def read_kofam_table(table: str) -> pd.DataFrame:
    ko_table = pd.read_table(table, usecols=["ko_num", "ko_def"])
    ko_table = ko_table.groupby("ko_num")[["ko_def"]].agg(",".join).unstack()
    ko_table = (
        ko_table.to_frame(name="ko_def").reset_index(level=0).drop(columns=["level_0"])
    )
    ko_table.ko_def = ko_table.ko_def.map(lambda x: ",".join(set(x.split(","))))
    return ko_table


def read_kofam_embedding(embedding: str) -> pd.DataFrame:
    embed_df = pd.read_table(embedding)
    embed_df["Sample Name"] = embed_df.filename.map(
        lambda fname: fname.replace(".kofamscan.tsv", "")
    )
    return embed_df


def render_ko_data_store(app: DashProxy, ko_embed_df: pd.DataFrame) -> dcc.Store:
    @app.callback(
        Output(Ids.ko_data_store, "data"),
        Input(Ids.ko_scatterplot, "selectedData"),
        Input(Ids.factor_name_text_input, "value"),
        Input(Ids.factor_value_text_input, "value"),
        Input(Ids.factor_unclassified_text_input, "value"),
        State(Ids.ko_data_store, "data"),
    )
    def update_ko_data(
        selected_data: Dict[str, List[Dict[str, str]]],
        factor_name: str,
        factor_value: str,
        unclassified_value: str,
        ko_data: pd.DataFrame,
    ) -> pd.DataFrame:
        if not isinstance(ko_data, pd.DataFrame) and ko_data == None:
            ko_data = ko_embed_df.copy()
            return Serverside(ko_data)
        if not factor_name or not factor_value:
            raise PreventUpdate
        if not selected_data or not selected_data["points"]:
            raise PreventUpdate
        # Get points based on filename (index as these are unique)
        selected_points_filenames = [
            point["customdata"][2] for point in selected_data["points"]
        ]
        selected_points = ko_data.filename.isin(selected_points_filenames)
        if factor_name not in ko_data.columns:
            ko_data[factor_name] = unclassified_value
        ko_data.loc[selected_points, factor_name] = factor_value
        return Serverside(ko_data)

    return dcc.Store(Ids.ko_data_store, storage_type="session")


def render_ko_scatterplot(app: DashProxy) -> html.Div:
    @app.callback(
        Output(Ids.ko_scatterplot, "figure"),
        Input(Ids.factor_name_select, "value"),
        Input(Ids.ko_data_store, "data"),
    )
    def plot_by_factor(factor_name: Optional[str], ko_data: pd.DataFrame) -> go.Figure:
        layout = go.Layout(
            xaxis=go.layout.XAxis(title=dict(text=r"$\text{X}_{1}$")),
            yaxis=go.layout.YAxis(title=dict(text=r"$\text{X}_{2}$")),
            height=550,
            uirevision=factor_name,
            title="Genome KEGG ortholog embedding",
            template="plotly_white",
        )
        fig = px.scatter(
            ko_data,
            x="x_1",
            y="x_2",
            color=factor_name,
            hover_data={
                "x_1": False,
                "x_2": False,
                "Organism Name": True,
                "Sample Name": True,
                "filename": True,
            },
            labels=dict(x_1=r"$\text{X}_{1}$", x_2=r"$\text{X}_{2}$"),
            template="plotly_white",
        )
        return go.Figure(data=fig.data, layout=layout)

    return html.Div(dcc.Graph(id=Ids.ko_scatterplot, mathjax=True))


def render_factor_name_text_input() -> html.Div:
    return dmc.TextInput(id=Ids.factor_name_text_input, label="Factor name")


def render_factor_value_text_input() -> html.Div:
    return dmc.TextInput(id=Ids.factor_value_text_input, label="Factor value")


def render_factor_unclassified_text_input() -> html.Div:
    return dmc.TextInput(
        id=Ids.factor_unclassified_text_input,
        label="Factor's unclassified value",
        value="unclassified",
    )


def render_factor_dropdown(app: DashProxy) -> html.Div:
    @app.callback(
        Output(Ids.factor_name_select, "data"),
        Input(Ids.ko_data_store, "data"),
    )
    def update_factor_names(ko_data: pd.DataFrame) -> List[str]:
        return ko_data.columns.tolist()

    return html.Div(
        dmc.Select(
            id=Ids.factor_name_select,
            label="Select a factor",
            clearable=True,
            persistence="session",
        )
    )


def render_forest_criterion_dropdown() -> html.Div:
    # https://scikit-learn.org/stable/modules/generated/sklearn.ensemble.RandomForestClassifier.html#sklearn.ensemble.RandomForestClassifier
    return html.Div(
        dmc.Select(
            id=Ids.forest_criterion,
            label="Select the classifier criterion",
            description="The function to measure the quality of a split",
            clearable=False,
            data=[
                dict(label="Gini", value="gini"),
                dict(label="Entropy", value="entropy"),
                dict(label="Log loss", value="log_loss"),
            ],
            value="gini",
            persistence="session",
        )
    )


def render_feature_importance_button(app: DashProxy) -> html.Div:
    @app.callback(
        Output(Ids.feature_importance_btn, "disabled"),
        Input(Ids.factor_name_select, "value"),
    )
    def toggle_disabled(factor_name: str) -> bool:
        return not factor_name

    return html.Div(
        dmc.Button(
            "Get feature importances",
            id=Ids.feature_importance_btn,
            fullWidth=True,
            variant="light",
            color="dark",
        )
    )


def render_feature_importance_download_data_button(app: DashProxy) -> html.Div:
    @app.callback(
        Output(Ids.feature_data_download, "data"),
        Output(Ids.feature_data_download_button, "n_clicks"),
        Input(Ids.feature_importance_data_store, "data"),
        Input(Ids.forest_criterion, "value"),
        Input(Ids.factor_name_select, "value"),
        Input(Ids.feature_data_download_button, "n_clicks"),
        prevent_initial_call=True,
    )
    def retrieve_download_data(
        feature_importance: pd.DataFrame,
        criterion: str,
        factor_name: str,
        n_clicks: int,
    ) -> pd.DataFrame:
        if not n_clicks:
            raise PreventUpdate
        to_tsv_method = partial(feature_importance.to_csv, sep="\t")
        to_tsv_method = partial(to_tsv_method, header=True)
        to_tsv_method = partial(to_tsv_method, index=False)
        to_tsv_method.__name__ = "to_csv"  # Necessary hack to get around explicit 'to_csv' check in dash_extensions.enrich
        factor_name = factor_name.replace(" ", "_")
        return (
            dcc.send_data_frame(
                to_tsv_method, f"{factor_name}_{criterion}_feature_importances.tsv"
            ),
            0,
        )

    return html.Div(
        [
            dmc.Tooltip(
                label="Download feature importance data",
                children=[
                    dmc.ActionIcon(
                        DashIconify(icon="basil:download-outline"),
                        variant="light",
                        color="orange",
                        size="lg",
                        id=Ids.feature_data_download_button,
                    ),
                    dcc.Download(id=Ids.feature_data_download),
                ],
            )
        ]
    )


def render_ko_data_download_data_button(app: DashProxy) -> html.Div:
    @app.callback(
        Output(Ids.ko_data_download, "data"),
        Output(Ids.ko_data_download_button, "n_clicks"),
        Input(Ids.ko_data_store, "data"),
        Input(Ids.ko_data_download_button, "n_clicks"),
        prevent_initial_call=True,
    )
    def retrieve_download_data(
        ko_data: pd.DataFrame,
        n_clicks: int,
    ) -> pd.DataFrame:
        if not n_clicks:
            raise PreventUpdate
        to_tsv_method = partial(ko_data.to_csv, sep="\t")
        to_tsv_method = partial(to_tsv_method, header=True)
        to_tsv_method = partial(to_tsv_method, index=False)
        to_tsv_method.__name__ = "to_csv"  # Necessary hack to get around explicit 'to_csv' check in dash_extensions.enrich
        return dcc.send_data_frame(to_tsv_method, "metabolism_feature_analysis.tsv"), 0

    return html.Div(
        [
            dmc.Tooltip(
                label="Download feature analysis data",
                children=[
                    dmc.ActionIcon(
                        DashIconify(icon="basil:download-outline"),
                        variant="light",
                        color="yellow",
                        size="lg",
                        id=Ids.ko_data_download_button,
                    ),
                    dcc.Download(id=Ids.ko_data_download),
                ],
            )
        ]
    )


def render_feature_importance_clf_store(
    app: DashProxy, ko_matrix: pd.DataFrame
) -> dcc.Store:
    @app.callback(
        Output(Ids.feature_importance_clf_store, "data"),
        Input(Ids.ko_data_store, "data"),
        Input(Ids.factor_name_select, "value"),
        Input(Ids.forest_criterion, "value"),
        Input(Ids.feature_importance_btn, "n_clicks"),
        prevent_initial_call=True,
    )
    def update_classifier_and_train_test_split_data(
        ko_data: pd.DataFrame,
        factor_name: str,
        forest_criterion: str,
        feature_importance_btn: int,
    ) -> Tuple[
        RandomForestClassifier,
        pd.DataFrame,
        pd.DataFrame,
        pd.DataFrame,
        pd.DataFrame,
        List[str],
    ]:
        if not feature_importance_btn or not factor_name:
            raise PreventUpdate
        # Prepare features, target labels and train classifier
        factor_values = ko_data.set_index("filename")[factor_name]
        X = ko_matrix.loc[factor_values.index].copy().to_numpy()
        X_names = ko_matrix.columns.tolist()
        y = factor_values.to_numpy()
        X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=42)
        clf = RandomForestClassifier(
            criterion=forest_criterion,
            random_state=42,
            n_estimators=1000,
        )
        clf.fit(X_train, y_train)
        return Serverside((clf, X_train, X_test, y_train, y_test, X_names))

    return dcc.Store(id=Ids.feature_importance_clf_store, storage_type="session")


def render_feature_importance_store(
    app: DashProxy, ko_table: pd.DataFrame
) -> dcc.Store:
    @app.callback(
        Output(Ids.feature_importance_data_store, "data"),
        Input(Ids.feature_importance_clf_store, "data"),
    )
    def update_feature_importance_data(
        clf_store_data: Tuple[
            RandomForestClassifier,
            pd.DataFrame,
            pd.DataFrame,
            pd.DataFrame,
            pd.DataFrame,
            List[str],
        ],
    ):
        if not clf_store_data:
            raise PreventUpdate
        clf, *__, X_names = clf_store_data
        importances_df = pd.Series(
            clf.feature_importances_,
            index=X_names,
            name="feature_importance",
        )
        # TODO: Incorporate other relevant feature analysis info
        # forest.classes_
        # forest.get_params()
        # zero_importance = importances_df.loc[importances_df.eq(0)]
        # nonzero_importance = importances_df.loc[importances_df.ne(0)].sort_values(
        #     ascending=False
        # )
        ko_features_importance_stddev = np.std(
            [tree.feature_importances_ for tree in clf.estimators_], axis=0
        )
        importances_df = pd.merge(
            importances_df,
            ko_table,
            how="left",
            left_index=True,
            right_index=True,
        ).reset_index(names=["ko_num"])
        ko_features_importance_stddev = pd.Series(
            ko_features_importance_stddev,
            index=X_names,
            name="feature_importance_std",
        )
        importances_df = pd.merge(
            importances_df,
            ko_features_importance_stddev,
            how="left",
            left_on="ko_num",
            right_index=True,
        ).sort_values(by="feature_importance", ascending=False)
        return Serverside(importances_df)

    return dcc.Store(Ids.feature_importance_data_store)


def render_clf_score_badge(app: DashProxy) -> html.Div:
    @app.callback(
        Output(Ids.clf_score_badge, "children"),
        Input(Ids.feature_importance_clf_store, "data"),
    )
    def update_classifier_score(
        clf_store_data: Tuple[
            RandomForestClassifier,
            pd.DataFrame,
            pd.DataFrame,
            pd.DataFrame,
            pd.DataFrame,
            List[str],
        ]
    ) -> str:
        if not clf_store_data:
            raise PreventUpdate
        clf, __, X_test, __, y_test, __ = clf_store_data
        score = clf.score(X_test, y_test)
        return f"Classifier baseline accuracy on test data: {score:.3f}"

    return html.Div(
        dmc.Badge(
            "Classifier score", variant="outline", color="gray", id=Ids.clf_score_badge
        ),
    )


def render_feature_importances_table(app: DashProxy) -> html.Div:
    @app.callback(
        Output(Ids.feature_table, "rowData"),
        Input(Ids.feature_importance_data_store, "data"),
        Input(Ids.factor_name_select, "value"),
        prevent_initial_call=True,
    )
    def get_feature_importance_data(
        importances_df: pd.DataFrame,
        factor_name: str,
    ) -> List[str]:
        if not factor_name:
            raise PreventUpdate
        if not isinstance(importances_df, pd.DataFrame):
            raise PreventUpdate
        if isinstance(importances_df, pd.DataFrame) and importances_df.empty:
            raise PreventUpdate
        q = importances_df.query("feature_importance > 0").feature_importance.quantile(
            0.5
        )
        return (
            importances_df.query("feature_importance >= @q")
            .sort_values(by="feature_importance", ascending=False)
            .round(decimals=10)
            .to_dict("records")
        )

    columnDefs = [
        dict(
            field="ko_num",
            headerName="KEGG ortholog ID",
            width=175,
            minWidth=80,
            maxWidth=175,
        ),
        dict(
            field="feature_importance",
            headerName="feature importance",
            width=185,
            maxWidth=185,
        ),
        dict(
            field="feature_importance_std",
            headerName="feature importance std.",
            width=205,
            maxWidth=205,
        ),
        dict(field="ko_def", headerName="Definition"),
    ]
    return html.Div(
        dag.AgGrid(
            id=Ids.feature_table,
            columnDefs=columnDefs,
            defaultColDef=dict(resizable=True, sortable=True, filter=True),
            columnSize="sizeToFit",
            style=dict(height=1200),
        )
    )


def render_feature_importance_barchart(app: DashProxy) -> html.Div:
    @app.callback(
        Output(Ids.feature_importance_barchart, "figure"),
        Input(Ids.feature_importance_data_store, "data"),
        Input(Ids.factor_name_select, "value"),
        prevent_initial_call=True,
    )
    def get_feature_importance_data(
        importances_df: pd.DataFrame,
        factor_name: str,
    ) -> go.Figure:
        if not factor_name:
            raise PreventUpdate
        quantile = importances_df.query(
            "feature_importance > 0"
        ).feature_importance.quantile(0.5)
        data_frame = importances_df.query(
            "feature_importance >= @quantile"
        ).sort_values(by="feature_importance", ascending=True)
        fig = px.bar(
            data_frame,
            x="feature_importance",
            y="ko_num",
            hover_data=["ko_num", "ko_def", "feature_importance"],
            title=f"{factor_name} feature importance",
            template="plotly_white",
            orientation="h",
            labels=dict(
                ko_num="KEGG Ortholog ID",
                ko_def="KO definition",
                feature_importance="Feature importance",
            ),
            height=1200,
            text="ko_def",
            text_auto=True,
        )
        fig.update_traces(textposition="outside")
        return fig

    return html.Div(
        dcc.Loading(
            dcc.Graph(
                id=Ids.feature_importance_barchart,
                figure=px.bar(
                    title=f"Factor feature importance",
                    template="plotly_white",
                    orientation="h",
                    height=1200,
                ),
            )
        )
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--matrix", help="path to kofamscan_results_matrix.tsv", required=True
    )
    parser.add_argument(
        "--table", help="path to kofamscan_results_table.tsv", required=True
    )
    parser.add_argument(
        "--embedding", help="path to kofamscan_results_embedding.tsv", required=True
    )
    parser.add_argument(
        "--debug",
        help="Set app.debug to True",
        required=False,
        action="store_true",
        default=False,
    )
    args = parser.parse_args()

    ko_matrix = read_kofam_matrix(args.matrix)  # features
    ko_table = read_kofam_table(args.table)  # ko_nums to ko_defs
    ko_embedding = read_kofam_embedding(args.embedding)

    app = DashProxy(
        title="Metabolism CoDa feature analysis",
        transforms=[ServersideOutputTransform()],
    )
    app.layout = dmc.Container(
        [
            render_ko_data_store(app, ko_embedding),
            render_feature_importance_clf_store(app, ko_matrix),
            render_feature_importance_store(app, ko_table),
            html.H1("Metabolism CoDa feature analysis"),
            render_ko_scatterplot(app),
            dmc.Group(
                [
                    render_factor_name_text_input(),
                    render_factor_value_text_input(),
                    render_factor_unclassified_text_input(),
                    render_factor_dropdown(app),
                ],
                position="center",
                grow=True,
                spacing="sm",
                align="end",
            ),
            dmc.Group(
                [
                    render_forest_criterion_dropdown(),
                    render_feature_importance_button(app),
                    dmc.Group(
                        [
                            render_clf_score_badge(app),
                            render_feature_importance_download_data_button(app),
                            render_ko_data_download_data_button(app),
                        ],
                        spacing="xs",
                        position="right",
                    ),
                ],
                position="center",
                grow=True,
                spacing="sm",
                align="end",
            ),
            dmc.Space(h=20),
            dmc.Group(
                [
                    render_feature_importance_barchart(app),
                    render_feature_importances_table(app),
                ],
                grow=True,
                align="apart",
                spacing="md",
                position="center",
            ),
        ],
        fluid=True,
    )

    app.run(debug=args.debug)


if __name__ == "__main__":
    main()
