"""Make the layout of HomologyViz GUI.

License
-------
This file is part of HomologyViz
BSD 3-Clause License
Copyright (c) 2024, Ivan Munoz Gutierrez
"""

import dash_ag_grid as dag
from dash import Dash, html, dcc
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc
from dash.development.base_component import Component
from dash_iconify import DashIconify
import plotly.express as px


TAB_LABEL_STYLE = {
    "fontSize": "14px",
    "padding": "0.6rem 1rem",
}


def make_dmc_select(**kwargs) -> Component:
    return dmc.Select(
        w=200,
        size="md",
        style={"padding": "0"},
        styles={
            "input": {"fontSize": "14px"},
            "label": {"fontSize": "14px"},
            "option": {"fontSize": "14px"},
        },
        **kwargs,
    )


def list_sequential_color_scales() -> list[str]:
    """Provide a list of Ploly's sequential color scales"""
    sequential_color_scales = [
        name for name in dir(px.colors.sequential) if not name.startswith("_")
    ]
    return sequential_color_scales


def make_tab_main() -> dbc.Tab:
    """Make tab Main."""
    tab_main = dbc.Tab(
        label="Main",
        tab_id="tab-main",
        label_style=TAB_LABEL_STYLE,
        children=[
            dbc.Row(  # ==== UPLOAD FILES SECTION ====================================== #
                [
                    dcc.Upload(
                        id="upload",
                        children=dmc.Button(
                            "Drag & Drop or Browse Files",
                            color="#3a7ebf",
                            leftSection=DashIconify(
                                icon="bytesize:upload",
                                width=25,
                            ),
                            variant="outline",
                            size="md",
                            style={
                                "fontSize": "14px",
                                "borderStyle": "dashed",
                                "borderWidth": "2px",
                                "width": "100%",
                                "height": "60px",
                            },
                        ),
                        multiple=True,
                        accept=".gb, .gbk",
                        className="d-flex justify-content-center",
                    ),
                    html.Div(  # Div to center AgGrid
                        [
                            dag.AgGrid(  # ==== TABLE TO DISPLAY FILE NAMES AND PATHS == #
                                id="files-table",
                                columnDefs=[
                                    {
                                        "headerName": "File Name",
                                        "field": "filename",
                                        "rowDrag": True,
                                        "sortable": True,
                                        "editable": False,
                                        "checkboxSelection": True,
                                        "headerCheckboxSelection": True,
                                        "cellStyle": {"fontSize": "12px"},
                                    },
                                ],
                                defaultColDef={"resizable": True},
                                dashGridOptions={
                                    "rowDragManaged": True,
                                    "localeText": {"noRowsToShow": "No Uploaded Files"},
                                    "rowSelection": "multiple",
                                },
                                rowData=[],  # Empty at start
                                columnSize="sizeToFit",
                                style={
                                    "height": "250px",
                                    "width": "100%",
                                    "fontSize": "14px",
                                },
                                className="ag-theme-alpine-dark",
                            ),
                        ],
                        style={"margin": "10px"},
                        className="d-flex justify-content-center",
                    ),
                ],
                className="d-flex justify-content-center mt-3",
                style={
                    "margin": "2px",
                },
            ),
            dbc.Row(  # ==== PLOT SECTION ============================================== #
                [
                    dbc.Row(
                        [
                            dmc.Button(
                                "Trash Selected Files",
                                id="delete-selected-files-button",
                                leftSection=DashIconify(
                                    icon="material-symbols-light:delete-outline-rounded",
                                    width=25,
                                ),
                                color="#3a7ebf",
                                size="md",
                                style={"fontSize": "14px", "width": "200px"},
                            ),
                        ],
                        className="d-flex justify-content-evenly mb-2",
                    ),
                    dbc.Row(
                        [
                            dmc.Button(  # RESET
                                "Reset",
                                id="reset-button",
                                leftSection=DashIconify(
                                    icon="material-symbols-light:reset-settings-rounded",
                                    width=25,
                                ),
                                color="#3a7ebf",
                                size="md",
                                style={"fontSize": "14px", "width": "200px"},
                            ),
                        ],
                        className="d-flex justify-content-evenly mb-2",
                    ),
                    dbc.Row(
                        [
                            dmc.Button(  # ERASE
                                "Erase Plot",
                                id="erase-button",
                                leftSection=DashIconify(
                                    icon="clarity:eraser-line",
                                    width=20,
                                ),
                                color="#3a7ebf",
                                size="md",
                                style={"fontSize": "14px", "width": "200px"},
                            ),
                        ],
                        className="d-flex justify-content-evenly mb-2",
                    ),
                    dbc.Row(
                        [
                            dmc.Button(  # DRAW PLOT
                                "Plot",
                                id="draw-button",
                                leftSection=DashIconify(
                                    icon="stash:pencil-writing-light",
                                    width=25,
                                ),
                                color="#b303b3",
                                size="md",
                                style={"fontSize": "14px", "width": "200px"},
                            ),
                        ],
                        className="d-flex justify-content-evenly mb-2",
                    ),
                ],
                className="d-flex justify-content-center mt-2",
                style={"margin": "2px"},
            ),
        ],
    )
    return tab_main


def make_tab_view() -> dbc.Tab:
    """Make tab view."""
    tab_view = dbc.Tab(
        label="View",
        tab_id="tab-view",
        label_style=TAB_LABEL_STYLE,
        children=[
            dbc.Row(
                [
                    dbc.Row(
                        make_dmc_select(
                            label="Align Plot",
                            id="align-plot",
                            value="left",
                            data=[
                                {"value": "left", "label": "Left"},
                                {"value": "center", "label": "Center"},
                                {"value": "right", "label": "Right"},
                            ],
                        ),
                        className="d-flex justify-content-evenly my-1",
                    ),
                    dbc.Row(
                        make_dmc_select(
                            label="Get Genes Info From",
                            id="use-genes-info-from",
                            value="gene",
                            data=[
                                {"value": "gene", "label": "CDS Gene"},
                                {"value": "product", "label": "CDS Product"},
                            ],
                        ),
                        className="d-flex justify-content-evenly mb-1",
                    ),
                    dbc.Row(
                        make_dmc_select(
                            id="annotate-genes",
                            label="Annotate Genes",
                            value="no",
                            data=[
                                {"value": "no", "label": "No"},
                                {"value": "top", "label": "Top genes"},
                                {"value": "bottom", "label": "Bottom genes"},
                                {
                                    "value": "top-bottom",
                                    "label": "Top and bottom genes",
                                },
                            ],
                        ),
                        className="d-flex justify-content-evenly mb-1",
                    ),
                    dbc.Row(
                        make_dmc_select(
                            id="annotate-sequences",
                            label="Annotate Sequences",
                            value="no",
                            data=[
                                {"value": "no", "label": "No"},
                                {"value": "accession", "label": "Accession"},
                                {"value": "name", "label": "Sequence name"},
                                {"value": "fname", "label": "File name"},
                            ],
                        ),
                        className="d-flex justify-content-evenly mb-1",
                    ),
                    dbc.Row(
                        make_dmc_select(
                            id="scale-bar",
                            label="Scale Bar",
                            value="yes",
                            data=[
                                {"value": "no", "label": "No"},
                                {"value": "yes", "label": "Yes"},
                            ],
                        ),
                        className="d-flex justify-content-evenly mb-1",
                    ),
                    dbc.Row(
                        dmc.NumberInput(  # Minimun homology lenght to plot
                            id="minimum-homology-length",
                            label="Min Homology Length",
                            value=100,
                            min=1,
                            step=50,
                            w=200,
                            suffix=" bp",
                            size="md",
                            style={"padding": "0"},
                            styles={
                                "input": {"fontSize": "14px"},
                                "label": {"fontSize": "14px"},
                            },
                        ),
                        className="d-flex justify-content-evenly mb-1",
                    ),
                    dbc.Row(
                        dmc.Button(
                            "Update View",
                            id="update-annotations",
                            leftSection=DashIconify(
                                icon="radix-icons:update",
                                width=25,
                            ),
                            color="#b303b3",
                            size="md",
                            style={
                                "fontSize": "14px",
                                "width": "200px",
                                "padding": "0",
                            },
                        ),
                        className="d-flex justify-content-evenly mt-2",
                    ),
                ],
                className="d-flex justify-content-center mt-2",
                style={"margin": "5px"},
            ),
        ],
        style={"margin": "5px"},
    )
    return tab_view


def make_tab_edit() -> dbc.Tab:
    """Make tab edit."""
    tab_edit = dbc.Tab(
        label="Edit",
        tab_id="tab-edit",
        label_style=TAB_LABEL_STYLE,
        children=[
            dbc.Row(  # ==== SELECT AND CHANGE COLOR SECTION =========================== #
                [
                    dmc.Divider(className="mt-2 mb-1"),
                    dbc.Row(
                        [
                            dmc.ColorInput(
                                id="color-input",
                                label="Edit Color of Selected Items",
                                value="rgb(0, 255, 255)",
                                w=200,
                                format="rgb",
                                swatches=[
                                    "rgb(255,0,255)",
                                    "rgb(0,255,255)",
                                    "rgb(255,26,0)",
                                    "rgb(255,116,0)",
                                    "rgb(255,255,0)",
                                    "rgb(0,255,0)",
                                    "rgb(151,59,255)",
                                    "rgb(0,0,0)",
                                ],
                                size="md",
                                style={"padding": "0"},
                                styles={
                                    "input": {"fontSize": "14px"},
                                    "label": {"fontSize": "14px"},
                                },
                            ),
                        ],
                        className="d-flex justify-content-evenly my-2",
                    ),
                    dbc.Row(
                        [
                            dmc.Button(
                                "Select Items",
                                id="select-change-color-button",
                                leftSection=DashIconify(
                                    icon="material-symbols-light:arrow-selector-tool-outline",
                                    width=30,
                                ),
                                color="#3a7ebf",
                                size="md",
                                variant="outline",
                                disabled=True,
                                style={"fontSize": "12px", "width": "200px"},
                            ),
                            dcc.Store(id="select-button-state-store", data=False),
                        ],
                        className="d-flex justify-content-evenly mb-2",
                    ),
                    dbc.Row(
                        [
                            dmc.Button(
                                "Change Color",
                                id="change-gene-color-button",
                                leftSection=DashIconify(
                                    icon="oui:color",
                                    width=20,
                                ),
                                color="#b303b3",
                                size="md",
                                style={"fontSize": "12px", "width": "200px"},
                            ),
                        ],
                        className="d-flex justify-content-evenly mb-2",
                    ),
                ],
                className="d-flex justify-content-center my-1",
                style={"margin": "2px"},
            ),
            dbc.Row(  # ==== SELECT COLORMAP FOR HOMOLOGY REGIONS ====================== #
                [
                    dmc.Divider(className="my-1"),
                    dbc.Row(
                        make_dmc_select(
                            id="color-scale",
                            label="Change Colormap of Homologies",
                            value="Greys",
                            data=list_sequential_color_scales(),
                        ),
                        className="d-flex justify-content-evenly mt-2 mb-2",
                    ),
                    dbc.Row(
                        html.Div(
                            dcc.Graph(
                                id="color-scale-display",
                                config={"displayModeBar": False, "staticPlot": True},
                                style={"width": "100%"},
                                className="border",
                            ),
                            style={"width": "90%"},
                        ),
                        className="d-flex justify-content-center mt-2 mb-1",
                    ),
                    dbc.Row(
                        html.Div(
                            dmc.RangeSlider(
                                id="range-slider",
                                value=[0, 75],
                                marks=[
                                    {"value": 25, "label": "25%"},
                                    {"value": 50, "label": "50%"},
                                    {"value": 75, "label": "75%"},
                                ],
                                size="md",
                                style={"width": "90%", "fontSize": "14px"},
                            ),
                            className="d-flex justify-content-center mt-1 mb-3",
                        ),
                    ),
                    dbc.Row(
                        [
                            html.Span("Truncate or Set Colormap to Extreme Homologies"),
                        ],
                        className="d-flex justify-content-center text-left mt-3 mb-0",
                        style={"fontSize": "14px", "width": "90%"},
                    ),
                    dbc.Row(
                        dmc.ButtonGroup(
                            [
                                dmc.Button(
                                    "Truncate",
                                    id="truncate-colorscale-button",
                                    variant="filled",
                                    size="md",
                                    style={
                                        "pointer-events": "none",
                                    },
                                    styles={
                                        "root": {
                                            "fontSize": "14px",
                                        }
                                    },
                                ),
                                dmc.Button(
                                    "Extreme Homol.",
                                    id="extreme-homologies-button",
                                    variant="subtle",
                                    size="md",
                                    styles={
                                        "root": {
                                            "fontSize": "14px",
                                        }
                                    },
                                ),
                                dcc.Store(
                                    id="is_set_to_extreme_homologies",
                                    data=False,
                                ),
                            ],
                            style={
                                "padding": "0px",
                                "borderWidth": "solid",
                                "borderStyle": "solid",
                                "borderColor": "#424242",
                                "borderRadius": "5px",
                                "backgroundColor": "#2e2e2e",
                            },
                        ),
                        style={"width": "85%"},
                        className="d-flex justify-content-evenly mb-2",
                    ),
                    dbc.Row(
                        [
                            dmc.Button(
                                "Update Homologies",
                                id="change-homology-color-button",
                                leftSection=DashIconify(
                                    icon="radix-icons:update",
                                    width=25,
                                ),
                                color="#b303b3",
                                size="md",
                                style={"fontSize": "14px", "width": "200px"},
                            ),
                        ],
                        className="d-flex justify-content-evenly my-2",
                    ),
                ],
                className="d-flex justify-content-center mt-2",
                style={"margin": "2px"},
            ),
        ],
    )
    return tab_edit


def make_tab_save() -> dbc.Tab:
    """Make tab save."""
    tab_save = dbc.Tab(
        label="Save",
        tab_id="tab-save",
        label_style=TAB_LABEL_STYLE,
        children=[
            dbc.Row(
                [
                    dbc.Row(
                        make_dmc_select(
                            label="Format",
                            id="figure-format",
                            value="png",
                            data=[
                                {"value": "png", "label": "png"},
                                {"value": "jpg", "label": "jpg"},
                                {"value": "pdf", "label": "pdf"},
                                {"value": "svg", "label": "svg"},
                                {"value": "html", "label": "html"},
                            ],
                        ),
                        className="d-flex justify-content-evenly my-1",
                    ),
                    dbc.Row(
                        dmc.NumberInput(
                            label="Width",
                            id="figure-width",
                            value=1200,
                            step=10,
                            w=200,
                            size="md",
                            suffix=" px",
                            style={"padding": "0"},
                            styles={
                                "input": {"fontSize": "14px"},
                                "label": {"fontSize": "14px"},
                            },
                        ),
                        className="d-flex justify-content-evenly mb-1",
                    ),
                    dbc.Row(
                        dmc.NumberInput(
                            label="Height",
                            id="figure-height",
                            value=1000,
                            step=10,
                            w=200,
                            size="md",
                            suffix=" px",
                            style={"padding": "0"},
                            styles={
                                "input": {"fontSize": "14px"},
                                "label": {"fontSize": "14px"},
                            },
                        ),
                        className="d-flex justify-content-evenly mb-1",
                    ),
                    dbc.Row(
                        dmc.NumberInput(
                            label="Scale",
                            id="figure-scale",
                            value=1,
                            step=0.2,
                            min=1,
                            max=10,
                            w=200,
                            size="md",
                            style={"padding": "0"},
                            styles={
                                "input": {"fontSize": "14px"},
                                "label": {"fontSize": "14px"},
                            },
                        ),
                        className="d-flex justify-content-evenly mb-3",
                    ),
                    dbc.Row(
                        [
                            dmc.Button(
                                "Download",
                                id="download-plot-button",
                                leftSection=DashIconify(
                                    icon="bytesize:download",
                                    width=25,
                                ),
                                variant="outline",
                                color="#3a7ebf",
                                size="md",
                                style={
                                    "fontSize": "14px",
                                    "borderWidth": "2px",
                                    "width": "200px",
                                },
                            ),
                            dcc.Download(id="download-plot-component"),
                        ],
                        className="d-flex justify-content-evenly mb-1",
                    ),
                ],
                className="d-flex justify-content-center mt-2",
                style={"margin": "5px"},
            ),
        ],
    )
    return tab_save


def create_layout(app: Dash) -> Dash:
    """Create app layout."""
    # Wrap layout with dmc.MantineProvider
    app.layout = dmc.MantineProvider(
        dmc.Grid(
            children=[
                dcc.Location(id="url", refresh=True),  # Allows refreshing app
                dmc.GridCol(
                    html.Div(  # ==== PLOT CONTROL ===================================== #
                        children=[
                            html.Img(
                                src="/assets/logo.png",
                                className="mx-auto my-2 d-block text-white fw-bold text-center",
                                alt="HomologyViz",
                                style={
                                    "height": "40px",
                                    "fontSize": "24px",
                                },
                            ),
                            html.Div(  # Tabs menu
                                dbc.Tabs(
                                    [
                                        make_tab_main(),
                                        make_tab_view(),
                                        make_tab_edit(),
                                        make_tab_save(),
                                    ],
                                    id="tabs",
                                ),
                                className="mt-1",
                                style={
                                    "height": "90%",
                                    "width": "100%",
                                    "overflow": "auto",
                                },
                            ),
                        ],
                        style={
                            "backgroundColor": "#242424",
                            "height": "96vh",
                            "overflow": "auto",
                        },
                    ),
                    span=3,
                ),
                dmc.GridCol(
                    html.Div(  # ==== GRAPH ============================================ #
                        children=[
                            dmc.Skeleton(
                                id="plot-skeleton",
                                visible=False,
                                children=dcc.Graph(
                                    id="plot",
                                    style={"height": "100%"},
                                ),
                                height="100%",
                            ),
                        ],
                        style={
                            "border": "1px solid black",
                            "height": "96vh",
                        },
                    ),
                    span=9,
                ),
            ],
            align="center",
            justify="flex-start",
            gutter="xs",
            style={"padding": "8px"},
        ),
        forceColorScheme="dark",
    )

    return app
