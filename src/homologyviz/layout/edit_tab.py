from dash import html, dcc
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc
from dash_iconify import DashIconify
import dash_ag_grid as dag
import plotly.express as px

from homologyviz.layout.common import TAB_LABEL_STYLE, make_dmc_select


# ==== EDIT TITLE ==== #
def make_update_button_insert_title() -> dmc.Button:
    return dmc.Button(
        "Update Title",
        id="update-title-button",
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
    )


def make_accordion_item_insert_title() -> dmc.AccordionItem:
    """
    Create a Dash Mantine Components AccordionItem for adding a title to the plot.

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
            dmc.AccordionControl("Title"),
            dmc.AccordionPanel(
                dbc.Row(
                    [
                        dbc.Row(
                            [
                                dmc.TextInput(
                                    id="title-input",
                                    placeholder="Plot Title",
                                    style={"width": "100%"},
                                )
                            ],
                            className="d-flex justify-content-evenly my-2",
                        ),
                        dbc.Row(
                            [
                                make_update_button_insert_title(),
                            ],
                            className="d-flex justify-content-evenly my-2",
                        ),
                    ],
                    className="d-flex justify-content-center my-1",
                ),
            ),
        ],
        value="insert-title",
    )


# ==== EDIT SEQUENCES AND GENES ==== #
def make_accordion_item_edit_button(name: str, id: str) -> dmc.Button:
    """
    Create a Dash Mantine Components AccordionItem for editing the annotations of
    sequences.

    Parameters
    ----------
    name : str
        The name to display on the button.
    id : str
        The id to assign to the button.

    Returns
    -------
    dmc.Button
        A Dash Mantine Components Button that opens the offcanvas for editing annotations.
    """
    return dmc.Button(
        name,
        id=id,
        variant="light",
        color="gray",
        size="md",
        fullWidth=True,
        justify="space-between",
        rightSection=DashIconify(icon="lucide:chevron-right", width=16),
        style={
            "fontSize": "14px",
        },
    )


def make_accordion_item_edit_button_group() -> dmc.ButtonGroup:
    """
    Create a Dash Mantine Components ButtonGroup for editing the annotations of
    sequences and genes.

    Returns
    -------
    dmc.ButtonGroup
        A Dash Mantine Components ButtonGroup that contains buttons to open the offcanvas
        for editing sequence and gene annotations.
    """
    return dmc.ButtonGroup(
        children=[
            make_accordion_item_edit_button(
                name="Sequences", id="open-offcanvas-edit-sequence-annotations"
            ),
            make_accordion_item_edit_button(
                name="Genes", id="open-offcanvas-edit-gene-annotations"
            ),
        ],
        style={
            "width": "200px",
            "padding": "0px",
            "borderWidth": "1px",
            "borderStyle": "solid",
            "borderColor": "#424242",
            "borderRadius": "5px",
        },
        orientation="vertical",
    )


# ==== Edit Gene Annotations ==== #
def make_offcanvas_update_gene_annotations_positions() -> html.Div:
    """Make offcanvas segment control"""
    return html.Div(
        children=[
            dbc.Label(
                "Annotate Genes Positions:",
                style={
                    "fontSize": "16px",
                },
            ),
            dmc.SegmentedControl(
                id="gene-annotation-positions-choice",
                data=[
                    {"value": "no", "label": "No"},
                    {"value": "top", "label": "Top"},
                    {"value": "bottom", "label": "Bottom"},
                    {
                        "value": "top-bottom",
                        "label": "Top and bottom",
                    },
                    {"value": "all-above", "label": "All-above"},
                    {"value": "all-below", "label": "All-below"},
                ],
                value="no",
                color="#3a7ebf",
                radius="xl",
                size="md",
                withItemsBorders=True,
                styles={"fontSize": "16px"},
                fullWidth=False,
                transitionDuration=500,
            ),
        ]
    )


def make_offcanvas_update_gene_segment_control() -> html.Div:
    """Make offcanvas segment control"""
    return html.Div(
        children=[
            dbc.Label(
                "Annotate Genes From:",
                style={
                    "fontSize": "16px",
                },
            ),
            dmc.SegmentedControl(
                id="gene-annotation-from-choice",
                data=[
                    {"label": "Gene", "value": "gene"},
                    {"label": "Product", "value": "product"},
                    {"label": "Custom Name", "value": "custom_name"},
                ],
                value="gene",
                color="#3a7ebf",
                radius="xl",
                size="md",
                withItemsBorders=True,
                styles={"fontSize": "16px"},
                fullWidth=False,
                transitionDuration=500,
            ),
        ]
    )


def make_offcanvas_update_annotations_text() -> dmc.Text:
    return dmc.Text(
        "You can type a custom name in the ✎ Custom name column for annotation"
    )


def make_offcanvas_update_button(name: str, id: str) -> dmc.Button:
    return dmc.Button(
        name,
        id=id,
        leftSection=DashIconify(
            icon="radix-icons:update",
            width=25,
        ),
        color="#b303b3",
        size="md",
        style={
            "fontSize": "14px",
            "width": "150px",
        },
    )


def make_offcanvas_edit_gene_paper() -> dmc.Paper:
    """Make offcanvas paper for gene table"""
    return dmc.Paper(
        children=[
            dbc.Row(
                children=[
                    make_offcanvas_update_gene_annotations_positions(),
                ]
            ),
            dbc.Row(
                children=[
                    make_offcanvas_update_gene_segment_control(),
                ]
            ),
            dbc.Row(
                children=[
                    dbc.Col(
                        make_offcanvas_update_annotations_text(),
                    ),
                    dbc.Col(
                        make_offcanvas_update_button(
                            name="Update",
                            id="offcanvas-update-gene-annotations-button",
                        ),
                    ),
                ],
                justify="center",
                align="center",
                className="mt-1",
            ),
        ],
        shadow="md",
        radius="md",
        p="md",
        withBorder=True,
    )


def make_offcanvas_update_sequence_segment_control() -> html.Div:
    """Make offcanvas segment control"""
    return html.Div(
        children=[
            dbc.Label(
                "Annotate Sequences From:",
                style={
                    "fontSize": "16px",
                },
            ),
            dmc.SegmentedControl(
                id="annotation-column-choice",
                data=[
                    {"label": "No", "value": "no"},
                    {"label": "Accession", "value": "accession"},
                    {"label": "Record Name", "value": "record_name"},
                    {"label": "File Name", "value": "file_name"},
                    {"label": "Custom Name", "value": "custom_name"},
                ],
                value="no",
                color="#3a7ebf",
                radius="xl",
                size="md",
                withItemsBorders=True,
                styles={"fontSize": "16px"},
                fullWidth=False,
                transitionDuration=500,
            ),
        ]
    )


def make_offcanvas_edit_sequence_paper() -> dmc.Paper:
    """Make offcanvas paper for sequence table"""
    return dmc.Paper(
        children=[
            dbc.Row(
                children=[
                    make_offcanvas_update_sequence_segment_control(),
                ]
            ),
            dbc.Row(
                children=[
                    dbc.Col(
                        make_offcanvas_update_annotations_text(),
                    ),
                    dbc.Col(
                        make_offcanvas_update_button(
                            name="Update",
                            id="offcanvas-update-sequence-annotations-button",
                        ),
                    ),
                ],
                justify="center",
                align="center",
                className="mt-1",
            ),
        ],
        shadow="md",
        radius="md",
        p="md",
        withBorder=True,
    )


def make_offcanvas_gene_table() -> dag.AgGrid:
    """
    Make the offcanvas gene table

    Returns
    -------
    dag.AgGrid
    """
    return dag.AgGrid(
        id="gene-table",
        columnDefs=[
            {
                "headerName": "CDS",
                "field": "cds_number",
                "width": 80,
                "maxWidth": 100,
                "minWidth": 40,
                "suppressSizeToFit": True,
                "rowDrag": False,
                "sortable": False,
                "editable": False,
                "checkboxSelection": False,
                "headerCheckboxSelection": False,
                "cellStyle": {"fontSize": "12px"},
            },
            {
                "headerName": "Accession",
                "field": "accession",
                "width": 140,
                "maxWidth": 160,
                "minWidth": 100,
                "suppressSizeToFit": True,
                "rowDrag": False,
                "sortable": False,
                "editable": False,
                "checkboxSelection": False,
                "headerCheckboxSelection": False,
                "cellStyle": {"fontSize": "12px"},
            },
            {
                "headerName": "Gene",
                "field": "gene",
                "width": 140,
                "maxWidth": 160,
                "minWidth": 100,
                "suppressSizeToFit": True,
                "rowDrag": False,
                "sortable": False,
                "editable": True,
                "checkboxSelection": False,
                "headerCheckboxSelection": False,
                "cellStyle": {"fontSize": "12px"},
            },
            {
                "headerName": "Product",
                "field": "product",
                "width": 200,
                "minWidth": 100,
                "suppressSizeToFit": True,
                "rowDrag": False,
                "sortable": False,
                "editable": True,
                "checkboxSelection": False,
                "headerCheckboxSelection": False,
                "cellStyle": {"fontSize": "12px"},
            },
            {
                "headerName": "Custom name",
                "field": "custom_name",
                "headerTooltip": "Type a custom name for annotation",
                "sortable": False,
                "editable": True,
                "checkboxSelection": False,
                "headerCheckboxSelection": False,
                "cellStyle": {
                    "fontSize": "12px",
                    "fontStyle": "italic",
                },
            },
        ],
        className="ag-theme-alpine-dark",
    )


def make_offcanvas_edit_gene_annotations() -> dbc.Offcanvas:
    return dbc.Offcanvas(
        style={"width": "80%", "fontSize": "12px"},
        id="offcanvas-edit-gene-annotations",
        title="Edit Gene Annotations",
        is_open=False,
        backdrop="static",
        children=[
            dbc.Row(
                children=[
                    dbc.Col(
                        make_offcanvas_edit_gene_paper(),
                        width=10,
                    ),
                ],
                justify="center",
                align="center",
                className="mb-3",
            ),
            dbc.Row(
                children=[
                    dbc.Col(
                        make_offcanvas_gene_table(),
                        width=10,
                    ),
                ],
                className="mb-3",
                justify="center",
                align="center",
            ),
        ],
    )


# ==== Edit Sequence Annotations ==== #
def make_offcanvas_sequence_table() -> dag.AgGrid:
    """
    Make the offcanvas sequence table

    Returns
    -------
    dag.AgGrid
    """
    return dag.AgGrid(
        id="sequence-table",
        columnDefs=[
            {
                "headerName": "File",
                "field": "file_number",
                "width": 60,
                "maxWidth": 80,
                "minWidth": 40,
                "suppressSizeToFit": True,
                "rowDrag": False,
                "sortable": False,
                "editable": False,
                "checkboxSelection": False,
                "headerCheckboxSelection": False,
                "cellStyle": {"fontSize": "12px"},
            },
            {
                "headerName": "Accession",
                "field": "accession",
                "width": 120,
                "maxWidth": 140,
                "minWidth": 100,
                "suppressSizeToFit": True,
                "rowDrag": False,
                "sortable": False,
                "editable": False,
                "checkboxSelection": False,
                "headerCheckboxSelection": False,
                "cellStyle": {"fontSize": "12px"},
            },
            {
                "headerName": "Record name",
                "field": "record_name",
                "rowDrag": False,
                "sortable": False,
                "editable": False,
                "checkboxSelection": False,
                "headerCheckboxSelection": False,
                "cellStyle": {"fontSize": "12px"},
            },
            {
                "headerName": "File name",
                "field": "file_name",
                "rowDrag": False,
                "sortable": False,
                "editable": False,
                "checkboxSelection": False,
                "headerCheckboxSelection": False,
                "cellStyle": {"fontSize": "12px"},
            },
            {
                "headerName": "Custom name",
                "field": "custom_name",
                "headerTooltip": "Type a custom name for annotation",
                "rowDrag": False,
                "sortable": False,
                "editable": True,
                "checkboxSelection": False,
                "headerCheckboxSelection": False,
                "cellStyle": {
                    "fontSize": "12px",
                    "fontStyle": "italic",
                },
            },
        ],
        dashGridOptions={
            "tooltipShowDelay": 100,
        },
        defaultColDef={"resizable": True},
        columnSize="sizeToFit",
        style={
            "height": "400px",
            "width": "100%",
            "fontSize": "12px",
        },
        className="ag-theme-alpine-dark",
    )


def make_offcanvas_edit_sequence_annotations() -> dbc.Offcanvas:
    return dbc.Offcanvas(
        style={"width": "80%", "fontSize": "12px"},
        id="offcanvas-edit-sequence-annotations",
        title="Edit Sequence Annotations",
        is_open=False,
        backdrop="static",
        children=[
            dbc.Row(
                children=[
                    dbc.Col(
                        make_offcanvas_edit_sequence_paper(),
                        width=10,
                    ),
                ],
                className="mb-3",
                justify="center",
                align="center",
            ),
            dbc.Row(
                children=[
                    dbc.Col(
                        make_offcanvas_sequence_table(),
                        width=10,
                    ),
                ],
                justify="center",
            ),
        ],
    )


def make_accordion_item_edit_sequence_and_gene_annotations() -> dmc.AccordionItem:
    """
    Create the Dash Mantine Components AccordionItem for editing the annotations of
    sequences and genes.
    """
    return dmc.AccordionItem(
        children=[
            dmc.AccordionControl("Annotations"),
            dmc.AccordionPanel(
                children=[
                    make_accordion_item_edit_button_group(),
                    make_offcanvas_edit_sequence_annotations(),
                    make_offcanvas_edit_gene_annotations(),
                ],
                className="d-flex justify-content-evenly mb-2",
            ),
        ],
        value="edit-annotations",
    )


# ==== EDIT COLOR ==== #
def make_accordion_item_edit_color_input() -> dmc.ColorInput:
    """
    Create a Dash Mantine Components ColorInput for editing the color of selected items.

    Returns
    -------
    dmc.ColorInput
        A Dash Mantine Components ColorInput that allows users to select a color (HEX
        format) from predefined swatches or custom values.
    """
    return dmc.ColorInput(
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
    )


def make_accordion_item_edit_color_button(
    name: str, id: str, icon: str, color: str
) -> dmc.Button:
    return dmc.Button(
        name,
        id=id,
        leftSection=DashIconify(
            icon=icon,
            width=25,
        ),
        color=color,
        size="md",
        variant="outline",
        style={
            "fontSize": "14px",
            "width": "200px",
        },
    )


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
        children=[
            dmc.AccordionControl("Color of Selected Items"),
            dmc.AccordionPanel(
                dbc.Row(
                    children=[
                        dbc.Row(
                            children=make_accordion_item_edit_color_input(),
                            className="d-flex justify-content-evenly my-2",
                        ),
                        dbc.Row(
                            children=[
                                make_accordion_item_edit_color_button(
                                    name="Select Items",
                                    id="select-items-button",
                                    icon="mdi:cursor-default-outline",
                                    color="#3a7ebf",
                                ),
                                dcc.Store(
                                    id="select-items-button-store",
                                    data=False,
                                ),
                            ],
                            className="d-flex justify-content-evenly mb-2",
                        ),
                        dbc.Row(
                            children=[
                                make_accordion_item_edit_color_button(
                                    name="Change Color",
                                    id="change-gene-color-button",
                                    icon="mdi:color",
                                    color="#b303b3",
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


# ==== EDIT HOMOLOGY ==== #
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
            dmc.AccordionControl("Homology Colormap"),
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


def make_edit_tab() -> dbc.Tab:
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
                    make_accordion_item_insert_title(),
                    make_accordion_item_edit_sequence_and_gene_annotations(),
                    make_accordion_item_edit_color(),
                    make_accordion_item_homology(),
                ],
                variant="default",
                chevronPosition="left",
                className="mt-3",
            ),
        ],
    )
    return tab_edit
