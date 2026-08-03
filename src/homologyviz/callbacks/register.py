"""
Register Dash callback functions for the HomologyViz graphical interface.

This module wires together the core interactive components of the app, including:
- File upload and deletion for GenBank files
- Plot generation using BLASTn alignments
- UI controls for adjusting annotations, homology colors, and visibility
- Custom color selection and trace selection logic
- Application reset and download features
- Heartbeat monitoring to shut down the app when the browser tab is closed

All callbacks are registered through the `register_callbacks(app)` function.

Notes
-----
- This file is part of HomologyViz
- BSD 3-Clause License
- Copyright (c) 2024, Iván Muñoz Gutiérrez
"""

import base64
import tempfile
import atexit
from pathlib import Path
from io import BytesIO
import os
import signal
import time
import threading
from flask import request, jsonify, Response
import json

import dash
from dash import Input, Output, State
from plotly.graph_objects import Figure

from homologyviz.parameters import PlotParameters
from homologyviz import plotter as plt
from homologyviz.callbacks.electrocardiograph import HeartBeatsParameters
from homologyviz.callbacks.files_manipulation import save_uploaded_file
from homologyviz.callbacks.clicks import (
    handle_plot_button_click,
    handle_update_view_click,
    handle_update_homologies_click,
    handle_update_title_click,
    handle_change_color_click,
    handle_select_traces_click,
)
from homologyviz.callbacks.updates import (
    update_genes_annotations,
    update_dna_sequence_annotations,
)

# import logging
# logging.basicConfig(level=logging.INFO)
# logger = logging.getLogger(__name__)


def register_callbacks(app: dash.Dash) -> dash.Dash:
    """
    Register all Dash callbacks for the app, including plotting logic, UI interactins,
    and server shutdown monitoring.

    This function sets up the full interactivity of the Dash app, including:
        - Handling file uploads and deletions.
        - Executing BLASTn alignments and plotting homology regions.
        - Updating annotations, colors, layout, and display options.
        - Managing UI elements like buttons, skeleton loaders, and input states.
        - Generating downloadable figures in various formats.
        - Monitoring heartbeat pings from the frontend to detect tab closure and
          gracefully shut down the app server when inactive.

    Parameters
    ----------
    app : dash.Dash
        The Dash app instance to which all callback functions and server routes will
        be attached.

    Returns
    -------
    dash.Dash
        The same Dash app instance, now with all callbacks registered.
    """
    # Create the tmp directory and ensure it's deleted when the app stops
    tmp_directory = tempfile.TemporaryDirectory()
    atexit.register(tmp_directory.cleanup)
    tmp_directory_path = Path(tmp_directory.name)
    # Monitor the Dash app tab status
    heartbeat_parameters = HeartBeatsParameters()

    # Class to store alignments data
    dash_parameters = PlotParameters()

    # = Files-table for selected GenBank files ========================================= #
    @app.callback(
        Output("files-table", "rowData"),
        [
            Input("upload", "filename"),
            Input("upload", "contents"),
            Input("trash-selected-files-button", "n_clicks"),
        ],
        [
            State("files-table", "rowData"),
            State("files-table", "selectedRows"),
        ],
    )
    def update_file_table(
        filenames: list | None,
        contents: list | None,
        n_clicks: int | None,
        current_row_data: list[dict] | None,
        selected_rows: list[dict] | None,
    ) -> list[dict[str, str]]:
        """
        Update the GenBank files table based on uploaded files or deletion actions.

        This callback populates the table with uploaded files by decoding their content
        and saving them temporarily. It also supports removing selected rows when the
        "Trash Selected Files" button is clicked.

        Parameters
        ----------
        filenames : list[str] or None
            List of filenames uploaded via the upload component.
        contents : list[str] or None
            Corresponding list of base64-encoded file contents.
        n_clicks : int or None
            Number of items the delete button has been clicked.
        current_row_data : list[dict] or None
            Current content of the table (`files-table`) as a list of dictionaries.
        selected_rows : list[dict] or None
            Subset of `current_row_data` that the user has selected for deletion.

        Returns
        -------
        list[dict]
            Uploaded list of table rows reflecting uploaded files or deletions.
        """
        ctx = dash.callback_context
        ctx_id = ctx.triggered[0]["prop_id"].split(".")[0]
        print(f"clicked from update file table: {ctx_id}")
        # Update table with uploaded files.
        if (ctx_id == "upload") and filenames and contents:
            new_rows = []
            # Simulate saving each file and creating a temporary file path
            for name, content in zip(filenames, contents):
                file_path = save_uploaded_file(name, content, tmp_directory_path)
                new_rows.append({"filename": name, "filepath": file_path})

            # Append new filenames and file paths to the table data
            return current_row_data + new_rows if current_row_data else new_rows

        # Delete selected rows
        if ctx_id == "trash-selected-files-button":
            updated = [row for row in current_row_data if row not in selected_rows]
            return updated

        return current_row_data if current_row_data else []

    # = MAIN CALLBACK FUNCTION ========================================================= #
    @app.callback(
        [
            Output("plot", "figure"),
            Output("plot", "clickData"),
            Output("plot-skeleton", "visible"),
        ],
        [
            Input("plot-button", "n_clicks"),
            Input("erase-button", "n_clicks"),
            Input("plot", "clickData"),
            Input("change-homology-color-button", "n_clicks"),
            Input("change-gene-color-button", "n_clicks"),
            Input("update-annotations", "n_clicks"),
            Input("update-title-button", "n_clicks"),
            Input("offcanvas-update-sequence-annotations-button", "n_clicks"),
            Input("offcanvas-update-gene-annotations-button", "n_clicks"),
        ],
        [
            State("files-table", "virtualRowData"),
            State("tabs", "active_tab"),
            State("plot", "figure"),
            State("color-input", "value"),
            State("select-items-button-store", "data"),
            State("color-scale", "value"),
            State("range-slider", "value"),
            State("align-plot", "value"),
            State("homology-style", "value"),
            State("minimum-homology-length", "value"),
            State("is-set-to-extreme-homologies", "data"),
            State("scale-bar", "value"),
            State("gene-annotation-positions-choice", "value"),
            State("gene-annotation-from-choice", "value"),
            State("gene-table", "virtualRowData"),
            State("title-input", "value"),
            State("sequence-table", "virtualRowData"),
            State("annotation-column-choice", "value"),
        ],
        prevent_initial_call=True,
    )
    def main_plot(
        plot_button_clicks: int | None,
        clear_button_clicks: int | None,
        click_data: dict | None,
        change_homology_color_button_clicks: int | None,
        change_gene_color_button_clicks: int | None,
        update_annotations_clicks: int | None,
        update_title_button_clicks: int | None,
        update_sequence_annotations: int | None,
        update_gene_annotations: int | None,
        virtual: list[dict[str, str]],
        active_tab: str,
        figure_state: dict,
        color_input_state: str,
        select_items_state: bool,
        color_scale_state: str,
        range_slider_state: list[int, int],
        align_plot_state: str,
        homology_style_state: str,
        minimum_homology_length_state: int,
        is_set_to_extreme_homologies: bool,
        scale_bar_state: str,
        gene_annotation_positions_state: str,
        gene_annotation_from_state: str,
        gene_table_state: list[dict] | None,
        title_input_state: str,
        sequence_table_state: list[dict] | None,
        annotation_column_choice_state: str,
    ) -> tuple[Figure | dict, None, bool]:
        """
        Master callback function for generating and modifying the alignment plot.

        This function coordinates user interactions accross multiple tabs:
        - In the **Main** tab, it triggers the BLASTn alignments and generates the plot.
        - In the **Edit** tab, it enables trace selection and color modifications.
        - In the **View** tab, it updates gene/sequence annotations, scale bar visibility,
          alignment position, and homology filtering.

        The function also resets the plot when the user clicks the "Erase" button.

        Parameters
        ----------
        plot_button_clicks : int | None
            Number of times the "Plot" button has been clicked.
        clear_button_clicks : int | None
            Number of times the "Erase" button has been clicked.
        click_data : dict | None
            Data from clicking a trace on the plot (used for selecting/deselecting
            traces).
        change_homology_color_button_clicks : int | None
            Click count for the button that changes homology colors and legend.
        change_gene_color_button_clicks : int | None
            Click count for the button that updates selected gene trace colors.
        update_annotations_clicks : int | None
            Click count for the button that updates annotations and plot layout.
        update_title_button_clicks : int | None
            Click count for the button that updates the figure's title.
        update_sequence_annotations : int | None
            Click count for the button that updates DNA sequence annotations.
        update_gene_annotations : int | None
            Click count for the button that updates gene annotations.
        virtual : list[dict[str, str]] | None
            Virtual row data from the GenBank file upload table.
        active_tab : str
            The currently active tab in the UI (e.g., "tab-main", "tab-edit", "tab-view").
        figure_state : dict
            The current figure state stored in Dash, used to rebuild or modify the plot.
        color_input_state : str
            Color selected for updating gene trace colors (hex string).
        select_items_state : bool
            Whether the "Select Items" button is active for toggling trace selection.
        color_scale_state : str
            The selected colorscale name used for identity-based homology coloring.
        range_slider_state : list[int, int]
            The selected range of identity percentages used for color scaling.
        align_plot_state : str
            Alignment layout setting (e.g., "left", "center", "right").
        homology_style_state : str
            Whether the style of homology connections are straight or curve (Bezier).
        minimum_homology_length_state : int
            Minimum homology length (in bp) to display in the plot.
        is_set_to_extreme_homologies : bool
            Whether to stretch the color scale to the dataset min/max identity values.
        annotate_genes_state : str
            Whether and where to annotate gene features (e.g., "top", "bottom", "no").
        scale_bar_state : str
            Whether to show the scale bar ("yes" or "no").
        use_genes_info_from_state : str
            Source of gene labels used for annotation (e.g., "CDS product", "CDS gene").
        title_input_state : str
            String holding the figure's title.
        sequence_table_state : list[dict[str, str, str, str, str]] | None
            Table from 'Edit Sequence Annotations' with the following columns: File,
            Accession, Record Name, File name, and Custom name.
        annotation_column_choice_state : str
            Whether and how to annotate DNA sequences (e.g., "accession", "name", "no",
            "custom").

        Returns
        -------
        fig : plotly.graph_objects.Figure
            The updated Plotly figure, either newly created or modified.
        None
            Placeholder to reset `clickData` in Dash (prevents stuck selections).
        bool
            A flag (`False`) to hide the dmc.Skeleton loading component after plot
            rendering.

        Notes
        -----
        - Uses `dash.callback_context` to determine which button triggered the callback.
        - This function is central to all updates affecting the alignment plot.
        """
        # Use context to find the button that triggered the callback.
        ctx = dash.callback_context
        button_id = ctx.triggered[0]["prop_id"].split(".")[0]
        print(f"button_id: {button_id}")

        # = TAB MAIN =================================================================== #
        # Perform Alignments & Plot
        if (button_id == "plot-button") and virtual:
            print("Plot button clicked, starting alignment and plotting...")
            return handle_plot_button_click(
                dash_parameters=dash_parameters,
                virtual=virtual,
                tmp_directory_path=tmp_directory_path,
                align_plot_state=align_plot_state,
                color_scale_state=color_scale_state,
                range_slider_state=range_slider_state,
                is_set_to_extreme_homologies=is_set_to_extreme_homologies,
                annotation_column_choice_state=annotation_column_choice_state,
                annotate_genes_from_state=gene_annotation_from_state,
                annotate_genes_positions_state=gene_annotation_positions_state,
                homology_style_state=homology_style_state,
                minimum_homology_length_state=minimum_homology_length_state,
                scale_bar_state=scale_bar_state,
                title_input_state=title_input_state,
            )

        # Erase Plot & Reset All Parameters
        if button_id == "erase-button":
            dash_parameters.reset()
            # Return an empty figure, None for clickdata, and False for skeleton
            return {}, None, False

        # = TAB VIEW =================================================================== #
        # Update Annotations and View
        if button_id == "update-annotations":
            return handle_update_view_click(
                figure_state=figure_state,
                dash_parameters=dash_parameters,
                align_plot_state=align_plot_state,
                homology_style_state=homology_style_state,
                scale_bar_state=scale_bar_state,
                minimum_homology_length_state=minimum_homology_length_state,
            )

        # = TAB EDIT =================================================================== #
        # Change Homology Color Regions and Colorscale Bar Legend
        if button_id == "change-homology-color-button":
            return handle_update_homologies_click(
                figure_state,
                dash_parameters,
                color_scale_state,
                range_slider_state,
                is_set_to_extreme_homologies,
            )

        # Insert title to plot
        if button_id == "update-title-button":
            return handle_update_title_click(
                figure_state=figure_state,
                dash_parameters=dash_parameters,
                title_input_state=title_input_state,
            )

        # Update sequence annotations
        if button_id == "offcanvas-update-sequence-annotations-button":
            fig = Figure(data=figure_state["data"], layout=figure_state["layout"])
            fig = update_dna_sequence_annotations(
                fig=fig,
                dash_parameters=dash_parameters,
                annotation_column_choice_state=annotation_column_choice_state,
                table=sequence_table_state,
            )
            return fig, None, False

        # Update gene annotations
        if button_id == "offcanvas-update-gene-annotations-button":
            print("Updating gene annotations based on user preferences...")
            fig = Figure(data=figure_state["data"], layout=figure_state["layout"])
            fig = update_genes_annotations(
                fig=fig,
                dash_parameters=dash_parameters,
                annotate_genes_from_state=gene_annotation_from_state,
                annotate_genes_positions_state=gene_annotation_positions_state,
                table=gene_table_state,
            )
            return fig, None, False

        # Change Color of Selected Traces
        if button_id == "change-gene-color-button":
            return handle_change_color_click(
                figure_state,
                dash_parameters,
                color_input_state,
            )

        # Select Traces for Changing Color
        if (
            (active_tab == "tab-edit")
            and (click_data is not None)
            and select_items_state
        ):
            return handle_select_traces_click(
                figure_state,
                dash_parameters,
                click_data,
            )

        return figure_state, None, False

    # = Activate all update buttons only when there is a plot ========================== #
    @app.callback(
        [
            Output("erase-button", "disabled"),
            Output("update-annotations", "disabled"),
            Output("change-gene-color-button", "disabled"),
            Output("change-homology-color-button", "disabled"),
            Output("select-items-button", "disabled"),
            Output("update-title-button", "disabled"),
            Output("offcanvas-update-sequence-annotations-button", "disabled"),
            Output("offcanvas-update-gene-annotations-button", "disabled"),
        ],
        Input("plot", "figure"),
    )
    def toggle_update_buttons(figure: dict) -> list[bool]:
        """
        Enable or disable editing buttons based on whether a plot is currently displayed.

        This callback disables the erase, update view, select items, change color, update
        homologies, and update title buttons when no figure has been generated (i.e., the
        figure is empty). It re-enables them when valid plot data is available.

        Parameters
        ----------
        figure : dict
            The current Plotly figure dictionary from the graph component.

        Returns
        -------
        list[bool]
            A list of six boolean values corresponding to the disabled states of:
            [erase-button, update-annotations, change-gene-color-button,
             change-homology-color-button, select-items-button, update-title-button,
             offcanvas-update-sequence-annotations-button].
        """
        if figure and figure.get("data", []):
            return (False,) * 8
        return (True,) * 8

    # = Activate the plot button only if there are files in the input table ============ #
    @app.callback(
        [
            Output("plot-button", "disabled"),
            Output("trash-selected-files-button", "disabled"),
        ],
        Input("files-table", "rowData"),
    )
    def toggle_plot_button(row_data: list[dict[str, str]] | None) -> bool:
        """
        Enable or disable the plot and trash buttons based on the number of uploaded
        GenBank files.

        Both buttons are enabled only when the upload table contains at least two files.
        Otherwise, they remain disabled.

        Parameters
        ----------
        row_data : list[dict[str, str]] or None
            The current contents of the GenBank file upload table.

        Returns
        -------
        tuple[bool, bool]
            The disabled state of the plot and trash buttons, respectively.
        """
        enabled = row_data is not None and len(row_data) >= 2
        return (not enabled, not enabled)  # Disable both buttons if not enough files

    # = Activate the select items button when the user clicks on it ==================== #
    @app.callback(
        [
            Output("select-items-button", "variant"),
            Output("select-items-button-store", "data"),
        ],
        Input("select-items-button", "n_clicks"),
        State("select-items-button-store", "data"),
    )
    def toggle_select_items_button(
        n_clicks: int | None,
        is_active: bool,
    ) -> tuple[str, bool]:
        """
        Toggle the state and appearance of the 'Select Items' button when clicked.

        This callback switches the internal selection mode on or off and updates
        the button's visual style (`variant`) accordingly. When active, the button
        appears filled; when inactive, it appears outlined.

        Parameters
        ----------
        n_clicks : int or None
            Number of times the 'Select Items' button has been clicked.
        is_active : bool
            The current selection mode state stored in Dash.

        Returns
        -------
        variant : str
            The button style variant ("filled" if active, "outline" if inactive).
        is_active : bool
            The updated state of the selection mode.
        """
        if n_clicks:
            # Toggle the active state on click
            is_active = not is_active

        # Set button style based on the active state
        if is_active:
            button_style = "filled"
        else:
            button_style = "outline"
        return button_style, is_active

    # = Toggle the colorscale buttons depending on users preference ==================== #
    @app.callback(
        [
            Output("extreme-homologies-button", "variant"),
            Output("extreme-homologies-button", "style"),
            Output("truncate-colorscale-button", "variant"),
            Output("truncate-colorscale-button", "style"),
            Output("is-set-to-extreme-homologies", "data"),
        ],
        [
            Input("extreme-homologies-button", "n_clicks"),
            Input("truncate-colorscale-button", "n_clicks"),
        ],
    )
    def toggle_colorscale_buttons(
        extreme_clicks: int | None,
        truncate_clicks: int | None,
    ) -> tuple[str, dict, str, dict, bool]:
        """
        Toggle the state and appearance of the 'Extreme Homologies' and
        'Truncate Colorscale' buttons.

        This callback ensures that only one of the two color scale options is active
        at a time. The active button is visually styled as "filled" and interactive;
        the inactive button is styled as "subtle" and disabled (via CSS pointer-events).
        The corresponding state value (`is_set_to_extreme_homologies`) is also updated.

        Parameters
        ----------
        extreme_clicks : int or None
            Number of times the 'Extreme Homologies' button has been clicked.
        truncate_clicks : int or None
            Number of times the 'Truncate Colorscale' button has been clicked.

        Returns
        -------
        extreme_variant : str
            Variant for the 'Extreme Homologies' button ("filled" or "subtle").
        extreme_style : dict
            CSS style dictionary for the 'Extreme Homologies' button.
        truncate_variant : str
            Variant for the 'Truncate Colorscale' button ("filled" or "subtle").
        truncate_style : dict
            CSS style dictionary for the 'Truncate Colorscale' button.
        is_set_to_extreme : bool
            Whether the extreme homology range setting is active.
        """
        ctx = dash.callback_context

        option1 = (
            "subtle",
            {"width": "280px", "padding": "5px"},
            "filled",
            {"width": "280px", "padding": "5px", "pointer-events": "none"},
            False,
        )
        option2 = (
            "filled",
            {"width": "280px", "padding": "5px", "pointer-events": "none"},
            "subtle",
            {"width": "280px", "padding": "5px"},
            True,
        )

        if not ctx.triggered:
            return option1

        triggered_id = ctx.triggered[0]["prop_id"].split(".")[0]

        if triggered_id == "extreme-homologies-button":
            return option2
        elif triggered_id == "truncate-colorscale-button":
            return option1

        return option1

    # = Update the color scale gradient after a user selection ========================= #
    @app.callback(
        Output("color-scale-display", "figure"),
        Input("color-scale", "value"),
    )
    def update_color_scale(value: str) -> Figure:
        """
        Update the horizontal color gradient display based on the selected colorscale.

        This callback is triggered when the user selects a new colorscale from the
        dropdown menu in the `Edit` tab. It passes the selected value to the
        `create_color_line` function to generate a smooth gradient for visual feedback.

        Parameters
        ----------
        value : str
            The name of the selected Plotly sequential colorscale (e.g., "Greys",
            "Blues").

        Returns
        -------
        figure : plotly.graph_objects.Figure
            A Plotly figure displaying a horizontal gradient representing the selected colorscale.
        """
        return plt.create_color_line(value.capitalize())

    # = Update edit sequence annotations table ========================================= #
    @app.callback(
        Output("sequence-table", "rowData"),
        Input("open-offcanvas-edit-sequence-annotations", "n_clicks"),
        State("sequence-table", "rowData"),
        prevent_initial_call=True,
    )
    def update_edit_sequence_annotations_table(figure: dict, table) -> list[dict]:
        if dash_parameters.gb_df is None:
            return []
        if figure:
            rows = []
            for _, gb_file in dash_parameters.gb_df.iterrows():
                rows.append(
                    {
                        "file_number": (gb_file["file_number"] + 1),
                        "accession": gb_file["accession"],
                        "record_name": gb_file["record_name"],
                        "file_name": gb_file["file_name"],
                        "custom_name": gb_file["custom_name"],
                    }
                )
            return rows
        else:
            return []

    # = Update edit gene annotations table ============================================= #
    @app.callback(
        Output("gene-table", "rowData"),
        Input("open-offcanvas-edit-gene-annotations", "n_clicks"),
        State("gene-table", "rowData"),
        prevent_initial_call=True,
    )
    def update_edit_gene_annotations_table(figure: dict, table) -> list[dict]:
        if dash_parameters.cds_df is None:
            return []
        if figure:
            rows = []
            for _, cds in dash_parameters.cds_df.iterrows():
                rows.append(
                    {
                        "cds_number": (cds["cds_number"] + 1),
                        "gene": cds["gene"],
                        "product": cds["product"],
                        "accession": cds["accession"],
                        "custom_name": cds["custom_name"],
                        "file_number": cds["file_number"],
                        "start": cds["start"],
                        "end": cds["end"],
                        "strand": cds["strand"],
                    }
                )
            return rows
        else:
            return []

    # = Toggle the offcanvas for the edit sequence annotations options ================= #
    @app.callback(
        Output("offcanvas-edit-sequence-annotations", "is_open"),
        Input("open-offcanvas-edit-sequence-annotations", "n_clicks"),
        [State("offcanvas-edit-sequence-annotations", "is_open")],
    )
    def toggle_offcanvas_edit_sequence_annotations(n1, is_open):
        if n1:
            return not is_open
        return is_open

    # = Toggle the offcanvas for the edit gene annotations options ===================== #
    @app.callback(
        Output("offcanvas-edit-gene-annotations", "is_open"),
        Input("open-offcanvas-edit-gene-annotations", "n_clicks"),
        [State("offcanvas-edit-gene-annotations", "is_open")],
    )
    def toggle_offcanvas_edit_gene_annotations(n1, is_open):
        if n1:
            return not is_open
        return is_open

    # = RESET THE APP ================================================================== #
    @app.callback(
        Output("url", "href"),
        Input("reset-button", "n_clicks"),
        prevent_initial_call=True,
    )
    def reset_app(n_clicks: int | None) -> str:
        """
        Reload the app when the "Reset" button is clicked.

        This callback returns the current URL path ("/"), which triggers a full page
        reload in Dash. It serves as a way to reset the interface and clear any stored
        state.

        Parameters
        ----------
        n_clicks : int or None
            Number of times the "Reset" button has been clicked.

        Returns
        -------
        str
            The URL path ("/") to trigger a browser reload of the app.
        """
        print("clicked Reset and I am reseting...")
        if n_clicks:
            # Return the current URL to trigger a reload
            return "/"

    # = Download the plot according to user preference format ========================== #
    @app.callback(
        Output("download-plot-component", "data"),
        Input("download-plot-button", "n_clicks"),
        [
            State("plot", "figure"),
            State("figure-format", "value"),
            State("figure-scale", "value"),
            State("figure-width", "value"),
            State("figure-height", "value"),
        ],
        prevent_initial_call=True,
    )
    def download_plot(
        n_clicks: int | None,
        figure: dict,
        figure_format: str,
        scale: int,
        width: int,
        height: int,
    ) -> dict:
        """
        Generate downloadable plot data in the selected format when the user clicks the "Download" button.

        This callback converts the current Plotly figure into either an HTML string or
        a static image (PNG, JPEG, SVG, etc.), encodes it in base64, and returns the data
        in a format compatible with the `dmc.Download` component.

        Parameters
        ----------
        n_clicks : int or None
            Number of times the "Download" button has been clicked.
        figure : dict
            The current Plotly figure as a dictionary (from `dcc.Graph`).
        figure_format : str
            The desired output format ("html", "png", "jpeg", "svg", etc.).
        scale : int
            Scaling factor for image resolution (used for static exports).
        width : int
            Width of the exported figure in pixels.
        height : int
            Height of the exported figure in pixels.

        Returns
        -------
        dict
            A dictionary containing the base64-encoded content, filename, MIME type,
            and `base64=True` flag for download via `dmc.Download`.
        """
        # Convert figure dictionary into a Figure object
        fig = Figure(data=figure["data"], layout=figure["layout"])

        if figure_format == "html":
            html_content = fig.to_html(full_html=True, include_plotlyjs="cdn")
            figure_name = "plot.html"

            # Encode the HTML content to base64 for download
            encoded = base64.b64encode(html_content.encode()).decode()

            # Return data for dmc.Download to prompt a download
            return dict(
                base64=True, content=encoded, filename=figure_name, type="text/html"
            )

        # If user didn't select html convert Figure object into an image in the
        # chosen format and DPI
        else:
            # Create an in-memory bytes buffer
            buffer = BytesIO()

            fig.write_image(
                buffer,
                format=figure_format,
                width=width,
                height=height,
                scale=scale,
                engine="kaleido",
            )

            # Encode the buffer as a base64 string
            encoded = base64.b64encode(buffer.getvalue()).decode()
            figure_name = f"plot.{figure_format}"

            # Return data for dmc.Download to prompt a download
            return dict(
                base64=True, content=encoded, filename=figure_name, type=figure_format
            )

    # ↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓ CHECKING IF TAB WAS CLOSED TO KILL SERVER ↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓ #
    @app.server.route("/heartbeat", methods=["POST"])
    def heartbeat() -> tuple[Response, int]:
        """
        Receive heartbeat pings from the frontend to monitor whether the app tab is open.

        This route is periodically called by the frontend to indicate that the app is
        still active. It parses the POST payload (JSON or raw) and updates the internal
        heartbeat counter and timestamp. If no data is received, it returns a warning.

        Returns
        -------
        tuple
            A Flask response with a JSON payload indicating success or failure,
            and an HTTP status code (200 or 500).
        """
        try:
            data = None

            # Attempt to parse the JSON payload
            if request.is_json:
                data = request.get_json()
            elif request.data:
                data = json.loads(request.data.decode("utf-8"))

            # Handle cases where no data is received
            if not data:
                print("Warning: No data received in the heartbeat request.", flush=True)
                return jsonify(success=False, message="No data received"), 200

            counter = data.get("counter", 0)
            heartbeat_parameters.last_heartbeat["timestamp"] = time.time()
            heartbeat_parameters.last_heartbeat["counter"] = counter

            return jsonify(success=True), 200
        except Exception as e:
            print(f"Error in /heartbeat route: {e}", flush=True)
            return jsonify(success=False, error=str(e)), 500

    def monitor_heartbeats() -> None:
        """
        Continuously monitor heartbeat timestamps to detect tab closure and shut down the server.

        This function runs in a background thread and checks whether the most recent
        heartbeat has timed out (based on `heartbeat_parameters.timeout_seconds`).
        If no new heartbeat is detected for a set period, and the heartbeat counter
        remains unchanged, the server is shut down gracefully.
        """
        counter = 0
        while True:
            now = time.time()
            elapsed_time = now - heartbeat_parameters.last_heartbeat["timestamp"]
            counter += 1
            # If timeout occurs, shut down the server
            if elapsed_time > heartbeat_parameters.timeout_seconds:
                print("Timeout: No heartbeats. Checking if counter has stopped...")
                # Check if the counter has stopped increasing
                initial_counter = heartbeat_parameters.last_heartbeat["counter"]
                time.sleep(5)  # Wait to see if the counter increases
                if heartbeat_parameters.last_heartbeat["counter"] == initial_counter:
                    shutdown_server()
            time.sleep(1)  # Regular monitoring interval

    # STARTING HEARTBEATS!
    if not heartbeat_parameters.heartbeat_monitor_started:
        heartbeat_parameters.heartbeat_monitor_started = True
        print("Initiating heartbeat_monitor_started")
        # Start the monitoring thread
        threading.Thread(target=monitor_heartbeats, daemon=True).start()

    @app.server.route("/shutdown", methods=["POST"])
    def shutdown_server() -> tuple[str, int]:
        """
        Shut down the Dash server when triggered.

        This endpoint is called by `monitor_heartbeats` when the app tab is closed
        and no heartbeats are received for a prolonged period. It sends a SIGINT
        signal to terminate the current process.

        Returns
        -------
        tuple
            A string message and HTTP status code 200 indicating shutdown.
        """
        os.kill(os.getpid(), signal.SIGINT)  # Send a signal to terminate the process
        print("Server shutting down...")
        return "Server shutting down...", 200

    # ↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑ CHECKING IF TAB WAS CLOSED TO KILL SERVER ↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑ #

    return app
