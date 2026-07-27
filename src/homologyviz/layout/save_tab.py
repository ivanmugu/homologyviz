from dash import dcc
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc
from dash_iconify import DashIconify

from homologyviz.layout.constants import TAB_LABEL_STYLE, make_dmc_select


def make_save_tab() -> dbc.Tab:
    """
    Create the 'Save' tab layout for exporting the plotted figure.

    This tab allows users to customize export settings and download the current plot
    in various formats. It provides controls to define output dimensions and scale.

    UI Elements:

    - Format selector (PNG, JPG, PDF, SVG, or HTML).
    - Numeric inputs for specifying figure width, height, and scale.
    - Download button that triggers file generation and download.
    - Dash `dcc.Download` component to handle file delivery.

    Returns
    -------
    dbc.Tab
        A Dash Bootstrap Component Tab containing the UI layout for the "Save" tab.
    """
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
