"""
Define the layout for the HomologyViz graphical user interface (GUI).

This module builds the entire front-end layout of the HomologyViz Dash application using
Dash Mantine Components (DMC), Dash Bootstrap Components (DBC), and Plotly Graphs. The GUI
includes interactive controls to upload files, edit plots, adjust views, and export
figures. It is structured into multiple tabs—Main, View, Edit, and Save—and integrates
seamlessly with Dash callbacks.

License
-------
This file is part of HomologyViz
BSD 3-Clause License
Copyright (c) 2024, Iván Muñoz Gutiérrez
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
    """
    Create a styled Dash Mantine Components (DMC) Select element.

    This utility function returns a DMC Select component with predefined styling,
    including fixed width, padding, and consistent font size across input, label,
    and options. Additional keyword arguments are passed directly to the `dmc.Select`.

    Parameters
    ----------
    **kwargs : dict
        Additional properties to customize the Select component (e.g., `data`, `value`, `label`).

    Returns
    -------
    Component
        A configured `dmc.Select` component ready to be used in the Dash layout.
    """
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
    """
    List all Plotly sequential color scales.

    This function returns the names of all sequential color scale options available
    in `plotly.express.colors.sequential`. These color scales are typically used for
    gradient-style visualizations such as heatmaps or homology identity shading.

    Returns
    -------
    list of str
        A list of sequential color scale names (e.g., "Viridis", "Blues", "Greys").
    """
    sequential_color_scales = [
        name for name in dir(px.colors.sequential) if not name.startswith("_")
    ]
    return sequential_color_scales


def make_tab_main() -> dbc.Tab:
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


def make_tab_view() -> dbc.Tab:
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
                                {"value": "all-above", "label": "All genes above"},
                                {"value": "all-below", "label": "All genes below"},
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
                ],
                className="d-flex justify-content-center mt-2",
                style={"margin": "5px"},
            ),
        ],
        style={"margin": "5px"},
    )
    return tab_view


def make_accordion_item_edit_color() -> dmc.AccordionItem:
    """
    Create a Dash Mantine Components AccordionItem for editing the color of selected
    items.

    This UI component includes:

    - A `ColorInput` widget for selecting a color (HEX format) from predefined swatches
      or custom values.
    - A "Select Items" button to enable item selection mode within the plot.
    - A "Change Color" button to apply the selected color to the currently selected items.
    - A hidden `dcc.Store` to keep track of the selection mode state (enabled/disabled).

    Returns
    -------
    dmc.AccordionItem
        A fully constructed AccordionItem containing the color editing UI for selected
        plot items.

    Notes
    -----
    - The component assumes that callbacks elsewhere in the app handle selection logic and
      color application.
    - Styling is handled using Bootstrap classes (`d-flex`, `justify-content-evenly`,
      `my-2`, etc.) and inline styles.
    - Color swatches include commonly used HEX values to improve usability.

    Component IDs
    -------------
    - "color-input": The HEX color selector input.
    - "select-items-button": Triggers selection mode for interactive elements.
    - "select-items-button-store": A hidden Store tracking whether selection mode is
      active.
    - "change-gene-color-button": Applies the selected color to all currently selected
      items.
    """
    return dmc.AccordionItem(
        [
            dmc.AccordionControl("Edit Color of Selected Items"),
            dmc.AccordionPanel(
                dbc.Row(
                    [
                        dbc.Row(
                            [
                                dmc.ColorInput(
                                    id="color-input",
                                    value="#00FFFF",
                                    w=200,
                                    format="hex",
                                    swatches=[
                                        "#FF00FF",
                                        "#00FFFF",
                                        "#FF1A00",
                                        "#FF7400",
                                        "#FFFF00",
                                        "#00FF00",
                                        "#973BFF",
                                        "#000000",
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
                                    id="select-items-button",
                                    leftSection=DashIconify(
                                        icon="material-symbols-light:arrow-selector-tool-outline",
                                        width=30,
                                    ),
                                    color="#3a7ebf",
                                    size="md",
                                    variant="outline",
                                    style={
                                        "fontSize": "12px",
                                        "width": "200px",
                                    },
                                ),
                                dcc.Store(
                                    id="select-items-button-store",
                                    data=False,
                                ),
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
                                    style={
                                        "fontSize": "12px",
                                        "width": "200px",
                                    },
                                ),
                            ],
                            className="d-flex justify-content-evenly mb-2",
                        ),
                    ],
                    className="d-flex justify-content-center my-1",
                ),
            ),
        ],
        value="edit-color",
    )


def make_accordion_item_homology() -> dmc.AccordionItem:
    """
    Create a Dash Mantine Components AccordionItem for customizing homology region colors.

    This UI component allows users to:

    - Select a sequential color scale for homology identity shading.
    - Preview the selected colormap in a static Plotly graph.
    - Adjust the effective identity range using a range slider.
    - Choose between truncating the colormap or setting it to the full (extreme) homology
      range.
    - Apply changes to the visualization with a button click.

    Returns
    -------
    dmc.AccordionItem
        A fully constructed AccordionItem containing UI controls for modifying the
        homology color mapping in the plot.

    Notes
    -----
    - The dropdown menu (`make_dmc_select`) uses available sequential color scales.
    - A small preview of the current color scale is shown via a static `dcc.Graph`.
    - The range slider allows users to limit the range of identity values visualized
      (e.g., 0-75%).
    - Two buttons ("Truncate" and "Extreme") toggle how the color scale range is handled.
    - The "Update Homologies" button triggers a callback to re-render regions with the
      selected color scale and identity thresholds.

    Component IDs
    -------------
    - "color-scale": Dropdown for selecting a colormap.
    - "color-scale-display": Plotly graph displaying a preview of the colormap.
    - "range-slider": Slider to adjust the visible range of homology identity.
    - "truncate-colorscale-button": Button indicating colormap is truncated.
    - "extreme-homologies-button": Button for stretching the colormap to extremes.
    - "is-set-to-extreme-homologies": Hidden Store tracking colormap state.
    - "change-homology-color-button": Button to apply updated color mapping.
    """
    return dmc.AccordionItem(
        [
            dmc.AccordionControl("Change Homology Colormap"),
            dmc.AccordionPanel(
                dbc.Row(
                    [
                        dbc.Row(
                            make_dmc_select(
                                id="color-scale",
                                value="Greys",
                                data=list_sequential_color_scales(),
                            ),
                            className="d-flex justify-content-evenly mt-2 mb-2",
                        ),
                        dbc.Row(
                            html.Div(
                                dcc.Graph(
                                    id="color-scale-display",
                                    config={
                                        "displayModeBar": False,
                                        "staticPlot": True,
                                    },
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
                                    style={
                                        "width": "90%",
                                        "fontSize": "14px",
                                    },
                                ),
                                className="d-flex justify-content-center mt-1 mb-3",
                            ),
                        ),
                        dbc.Row(
                            [
                                html.Span(
                                    "Truncate or Set Colormap to Extreme Homologies"
                                ),
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
                                        "Extreme",
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
                                        id="is-set-to-extreme-homologies",
                                        data=False,
                                    ),
                                ],
                                style={
                                    "padding": "0px",
                                    "borderWidth": "1px",
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
                                    style={
                                        "fontSize": "14px",
                                        "width": "200px",
                                    },
                                ),
                            ],
                            className="d-flex justify-content-evenly my-2",
                        ),
                    ],
                    className="d-flex justify-content-center mt-2",
                ),
            ),
        ],
        value="edit-homology-regions",
    )


def make_tab_edit() -> dbc.Tab:
    """
    Create the 'Edit' tab layout for the HomologyViz interface.

    This tab allows users to customize visual aspects of the plot, including:

    - Selecting specific gene or homology traces and applying custom colors.
    - Picking from a list of predefined colors using a color input.
    - Changing the colormap used for homology identity shading.
    - Adjusting the colormap range (e.g., truncating or setting extreme bounds).
    - Previewing the selected colormap via a horizontal colorbar.
    - Updating the plot to reflect all visual changes.

    UI Elements:

    - Color input with swatches and RGB support.
    - Buttons for selecting items and applying color changes.
    - Dropdown to choose a Plotly sequential colorscale.
    - Static plot to preview the colorscale.
    - Range slider to control truncation percentage.
    - Button group to toggle between truncating or fixing homology value bounds.
    - Button to apply the updated homology colormap.

    Returns
    -------
    dbc.Tab
        A Dash Bootstrap Component Tab containing the UI layout for the "Edit" tab.
    """
    tab_edit = dbc.Tab(
        label="Edit",
        tab_id="tab-edit",
        label_style=TAB_LABEL_STYLE,
        children=[
            dmc.Accordion(
                children=[
                    make_accordion_item_edit_color(),
                    make_accordion_item_homology(),
                ],
                variant="separated",
                className="mt-3",
            ),
        ],
    )
    return tab_edit


def make_tab_save() -> dbc.Tab:
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


def create_layout(app: Dash) -> Dash:
    """
    Construct the full layout for the HomologyViz Dash app.

    This function defines the GUI structure, including the control panel and plot display.
    It uses Dash Mantine Components for styling and layout organization. The layout is
    composed of two primary columns:

    - Left Column: Control panel with the HomologyViz logo and tabbed interface for
      uploading files, customizing views, editing plots, and saving outputs.
    - Right Column: Main plotting area displaying the generated figure using
      `dcc.Graph`, wrapped in a `dmc.Skeleton` for loading effects.

    Parameters
    ----------
    app : dash.Dash
        The Dash application instance to which the layout will be assigned.

    Returns
    -------
    dash.Dash
        The Dash app with its layout fully configured and assigned.
    """
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
                                className="mx-auto my-4 d-block text-white fw-bold text-center",
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
                                    "height": "85%",
                                    "width": "100%",
                                    "overflow": "auto",
                                },
                            ),
                        ],
                        style={
                            "backgroundColor": "#242424",
                            "height": "95vh",
                            "overflow": "auto",
                        },
                    ),
                    span="auto",
                    style={"maxWidth": "340px", "minWidth": "200px"},
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
