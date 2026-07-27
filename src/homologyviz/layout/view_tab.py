import dash_bootstrap_components as dbc
import dash_mantine_components as dmc
from dash_iconify import DashIconify

from homologyviz.layout.constants import TAB_LABEL_STYLE, make_dmc_select


def make_view_tab() -> dbc.Tab:
    """
    Create the 'View' tab layout for the HomologyViz interface.

    This tab allows users to customize how the DNA sequences and homology regions
    are displayed in the plot. Users can adjust layout alignment, annotations,
    and minimum homology length threshold.

    Features included:

    - Dropdowns for:
      - Aligning sequences (left, center, right)
      - Choosing gene info source (gene or product)
      - Annotating genes (none, top, bottom, or both)
      - Annotating DNA sequences (accession, name, or file name)
      - Toggling the scale bar
    - Number input to set the minimum homology length to display
    - Button to apply view updates to the plot

    Returns
    -------
    dbc.Tab
        A Dash Bootstrap Component Tab containing the UI layout for the "View" tab.
    """
    tab_view = dbc.Tab(
        label="View",
        tab_id="tab-view",
        label_style=TAB_LABEL_STYLE,
        children=[
            dbc.Row(
                [
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
                        className="d-flex justify-content-evenly mt-3 mb-1",
                    ),
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
                            label="Homology Connector Style",
                            id="homology-style",
                            value="straight",
                            data=[
                                {"value": "straight", "label": "Straight"},
                                {"value": "curve", "label": "Curve"},
                            ],
                        ),
                        className="d-flex justify-content-evenly my-1",
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
                ],
                className="d-flex justify-content-center mt-2",
                style={"margin": "5px"},
            ),
        ],
        style={"margin": "5px"},
    )
    return tab_view
