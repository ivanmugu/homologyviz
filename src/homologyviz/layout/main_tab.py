import dash_ag_grid as dag
import dash_bootstrap_components as dbc
from dash import html, dcc
import dash_mantine_components as dmc
from dash_iconify import DashIconify

from homologyviz.layout.constants import TAB_LABEL_STYLE


def make_main_tab() -> dbc.Tab:
    """
    Create the 'Main' tab layout for the HomologyViz interface.

    This tab provides users with the interface to upload GenBank files,
    manage them in a table, and control the main plotting functions. It includes:

    - A drag-and-drop upload area for `.gb` or `.gbk` files.
    - An AG Grid table to display and manage uploaded file names.
    - Buttons for:
        - Deleting selected files
        - Resetting the app
        - Erasing the plot
        - Drawing the plot

    Returns
    -------
    dbc.Tab
        A Dash Bootstrap Component Tab containing the UI layout for the "Main" tab.
    """
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
                                id="trash-selected-files-button",
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
                                id="plot-button",
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
