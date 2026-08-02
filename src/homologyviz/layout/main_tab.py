import dash_ag_grid as dag
import dash_bootstrap_components as dbc
from dash import html, dcc
import dash_mantine_components as dmc
from dash_iconify import DashIconify

from homologyviz.layout.common import TAB_LABEL_STYLE


# ==== UPLOAD FILES SECTION ============================================================ #
def make_upload_button() -> dcc.Upload:
    """
    Create a drag-and-drop upload button for `.gb` or `.gbk` files.

    Returns
    -------
    dcc.Upload
        A Dash Core Component Upload element configured for GenBank files.
    """
    return dcc.Upload(
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
        accept=".gb, .gbk, .bgff",
        className="d-flex justify-content-center",
    )


def make_upload_table() -> dag.AgGrid:
    """
    Create an AG Grid table to display uploaded file names.

    Returns
    -------
    dag.AgGrid
        An AG Grid component configured to display uploaded file names.
    """
    return dag.AgGrid(
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
    )


def upload_files_section() -> dbc.Row:
    """
    Create the upload files section of the 'Main' tab layout.

    This section provides users with a drag-and-drop upload area for `.gb` or `.gbk`
    files, and an AG Grid table to display and manage uploaded file names.

    Returns
    -------
    dbc.Row
        A Dash Bootstrap Component Row containing the upload files section layout.
    """
    return dbc.Row(
        children=[
            make_upload_button(),
            html.Div(  # Div to center AgGrid
                make_upload_table(),
                style={"margin": "10px"},
                className="d-flex justify-content-center",
            ),
        ],
        className="d-flex justify-content-center mt-3",
        style={
            "margin": "2px",
        },
    )


# ==== PLOT SECTION ==================================================================== #
def make_trash_button() -> dmc.Button:
    """
    Create a button to trigger the deletion of selected files.

    Returns
    -------
    dmc.Button
        A Dash Mantine Component Button configured to delete selected files.
    """
    return dmc.Button(
        "Trash Selected Files",
        id="trash-selected-files-button",
        leftSection=DashIconify(
            icon="material-symbols-light:delete-outline-rounded",
            width=25,
        ),
        color="#3a7ebf",
        size="md",
        style={"fontSize": "14px", "width": "200px"},
    )


def make_reset_button() -> dmc.Button:
    """
    Create a button to reset the application.

    Returns
    -------
    dmc.Button
        A Dash Mantine Component Button configured to reset the app.
    """
    return dmc.Button(
        "Reset",
        id="reset-button",
        leftSection=DashIconify(
            icon="material-symbols-light:reset-settings-rounded",
            width=25,
        ),
        color="#3a7ebf",
        size="md",
        style={"fontSize": "14px", "width": "200px"},
    )


def make_erase_button() -> dmc.Button:
    """
    Create a button to erase the plot.

    Returns
    -------
    dmc.Button
        A Dash Mantine Component Button configured to erase the plot.
    """
    return dmc.Button(
        "Erase Plot",
        id="erase-button",
        leftSection=DashIconify(
            icon="clarity:eraser-line",
            width=20,
        ),
        color="#3a7ebf",
        size="md",
        style={"fontSize": "14px", "width": "200px"},
    )


def make_plot_button() -> dmc.Button:
    """
    Create a button to trigger the plotting function.

    Returns
    -------
    dmc.Button
        A Dash Mantine Component Button configured to trigger the plot.
    """
    return dmc.Button(
        "Plot",
        id="plot-button",
        leftSection=DashIconify(
            icon="stash:pencil-writing-light",
            width=25,
        ),
        color="#b303b3",
        size="md",
        style={"fontSize": "14px", "width": "200px"},
    )


def plot_section_buttons() -> dbc.Row:
    """
    Create the plot section of the 'Main' tab layout.

    This section provides users with buttons to manage the plotting functions, including
    deleting selected files, resetting the app, erasing the plot, and drawing the plot.

    Returns
    -------
    dbc.Row
        A Dash Bootstrap Component Row containing the plot section layout.
    """
    classname_row = "d-flex justify-content-center mt-2"
    return dbc.Row(
        children=[
            dbc.Row(
                children=make_trash_button(),
                className=classname_row,
            ),
            dbc.Row(
                children=make_reset_button(),
                className=classname_row,
            ),
            dbc.Row(
                children=make_erase_button(),
                className=classname_row,
            ),
            dbc.Row(
                children=make_plot_button(),
                className=classname_row,
            ),
        ],
        className=classname_row,
        style={"margin": "2px"},
    )


# === MAIN TAB ========================================================================= #
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
            upload_files_section(),
            plot_section_buttons(),
        ],
    )
    return tab_main
