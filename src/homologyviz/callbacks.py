"""Register callbacks for HomologyViz GUI.

License
-------
This file is part of HomologyViz
BSD 3-Clause License
Copyright (c) 2024, Ivan Munoz Gutierrez
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
from flask import request, jsonify
import json

import dash
from dash import Input, Output, State
from plotly.graph_objects import Figure

from homologyviz import plotter as plt
from homologyviz.gb_files_manipulation import get_longest_sequence_dataframe


class HeartBeatsParameters:
    """Parameters to monitor heart beats of the Dash app.

    The monitoring of the heart beats allows to stop the server when the app tab is closed
    in the browser.

    Attributes
    ----------
    last_heartbeat : dict
        A dictionary storing the timestamp of the last heartbeat and a counter.
    timeout_seconds : int
        The number of seconds before a timeout occurs if no heartbeat is received.
    heartbeat_monitor_started : bool
        Whether the heartbeat monitor has been started
    """

    def __init__(
        self,
        last_heartbeat: dict | None = None,
        timeout_seconds: int = 5,
        heartbeat_monitor_started: bool = False,
    ) -> None:
        """Initialize HeartBeatsParameters

        Parameters
        ----------
        last_heartbeat : dict, optional
            Initial dictionary storing the timestamp and counter. Defaults to current time.
        timeout_seconds : int, optional
            Timeout duration in seconds. Default is 5 seconds.
        heartbeat_monitor_started : bool, optional
            Whether the monitor is started. Default is False.
        """
        self.last_heartbeat = (
            last_heartbeat
            if last_heartbeat is not None
            else {"timestamp": time.time(), "counter": 0}
        )
        self.timeout_seconds = timeout_seconds
        self.heartbeat_monitor_started = heartbeat_monitor_started


def save_uploaded_file(file_name, content, temp_folder_path: Path) -> str:
    """Decode the content and write it to a temporary file."""
    # Decode content
    data = content.split(";base64,")[1]
    decoded_data = base64.b64decode(data)

    # Save uploaded file
    output_path = temp_folder_path / file_name
    with open(output_path, "wb") as f:
        f.write(decoded_data)
    # Dash doesn't like Path; hence, we need to cast Path to str.
    return str(output_path)


def register_callbacks(app: dash.Dash) -> dash.Dash:
    """Callbacks for Dash app."""
    # Create the tmp directory and ensure it's deleted when the app stops
    tmp_directory = tempfile.TemporaryDirectory()
    atexit.register(tmp_directory.cleanup)
    tmp_path = Path(tmp_directory.name)

    # Monitor the Dash app tab status
    heartbeat_parameters = HeartBeatsParameters()

    # Class to store alignments data
    dash_parameters = plt.PlotParameters()

    # ==== files-table for selected GenBank files ====================================== #
    @app.callback(
        Output("files-table", "rowData"),
        [
            Input("upload", "filename"),
            Input("upload", "contents"),
            Input("delete-selected-files-button", "n_clicks"),
        ],
        [State("files-table", "rowData"), State("files-table", "selectedRows")],
    )
    def update_file_table(
        filenames, contents, n_clicks, current_row_data, selected_rows
    ):
        ctx = dash.callback_context
        ctx_id = ctx.triggered[0]["prop_id"].split(".")[0]
        # Update table with uploaded files.
        if (ctx_id == "upload") and filenames and contents:
            new_rows = []
            # Simulate saving each file and creating a temporary file path
            for name, content in zip(filenames, contents):
                file_path = save_uploaded_file(name, content, tmp_path)
                new_rows.append({"filename": name, "filepath": file_path})

            # Append new filenames and file paths to the table data
            return current_row_data + new_rows if current_row_data else new_rows

        # Delete selected rows
        if ctx_id == "delete-selected-files-button":
            updated = [row for row in current_row_data if row not in selected_rows]
            return updated

        return current_row_data if current_row_data else []

    # ==== plot the alignments -> MAIN FUNCTION ======================================== #
    @app.callback(
        [
            Output("plot", "figure"),
            Output("plot", "clickData"),
            Output("plot-skeleton", "visible"),
        ],
        [
            Input("draw-button", "n_clicks"),
            Input("erase-button", "n_clicks"),
            Input("plot", "clickData"),
            Input("change-homology-color-button", "n_clicks"),
            Input("change-gene-color-button", "n_clicks"),
            Input("update-annotations", "n_clicks"),
            Input("update-align-sequences-button", "n_clicks"),
        ],
        [
            State("files-table", "virtualRowData"),
            State("tabs", "active_tab"),
            State("plot", "figure"),
            State("color-input", "value"),
            State("select-button-state-store", "data"),
            State("color-scale", "value"),
            State("range-slider", "value"),
            State("align-plot", "value"),
            # State("homology-lines", "value"),
            State("minimum-homology-length", "value"),
            State("is_set_to_extreme_homologies", "data"),
            State("annotate-sequences", "value"),
            State("annotate-genes", "value"),
            State("scale-bar", "value"),
            State("use-genes-info-from", "value"),
        ],
        prevent_initial_call=True,
    )
    def main_plot(
        plot_button_clicks,  # input plot button click
        clear_button_clicks,  # input clear button click
        click_data,  # input click data form plot
        change_homology_color_button_clicks,  # input selected color scale value
        change_gene_color_button_clicks,  # input selected color value
        update_annotations_clicks,  # input to update annotate sequences
        update_align_sequences_clicks,  # input to update annotate sequences
        virtual,  # state of table with path to GenBank files
        active_tab,  # state activet tab
        figure_state,  # state output Figure object
        color_input_state,  # state color input
        select_button_state,  # state select button state store
        color_scale_state,  # state color scale
        range_slider_state,  # state range slider for color scale
        align_plot_state,  # state align plot
        # homology_lines_state,  # state homology lines
        minimum_homology_length_state,  # state miminum homology length
        is_set_to_extreme_homologies,  # state button colorscale range
        annotate_sequences_state,  # state annotate sequences
        annotate_genes_state,  # state annotate sequences
        scale_bar_state,  # state scale bar
        use_genes_info_from_state,  # state genes info
    ) -> Figure:
        """Main function controling the plot.

        Notes
        -----
        The Output to Plot—clickData to None in all the returns allows selecting and
        deselecting traces in the plot.
        """
        # Use context to find the button that triggered the callback.
        ctx = dash.callback_context
        button_id = ctx.triggered[0]["prop_id"].split(".")[0]

        # ============================================================================== #
        #                             MAIN TAB -> Plot                                   #
        # ============================================================================== #
        # ==== Plot Sequences ========================================================== #
        if (button_id == "draw-button") and virtual:
            # Store all the metadata for plotting into the dash_parameters object. The
            # dash_parameters object is a PlotParameters class declared at the top of
            # the `create_dash_app` function

            # Drawing plot from draw-button
            dash_parameters.draw_from_button = button_id

            # Convert paths to input_files and output_folder into Path objects
            input_files = [Path(row["filepath"]) for row in virtual]
            output_folder = Path(tmp_directory.name)
            # Store input_files and output_folder Paths in dash_parameters
            dash_parameters.input_files = input_files
            dash_parameters.output_folder = output_folder
            # Save number of files
            dash_parameters.number_gb_records = len(input_files)
            # Make alignments
            gb_df, cds_df, alignments_df, regions_df = plt.make_alignments(
                input_files, output_folder
            )
            # Add alignments to dash_parameters
            dash_parameters.gb_df = gb_df
            dash_parameters.cds_df = cds_df
            dash_parameters.alignments_df = alignments_df
            dash_parameters.alignments_regions_df = regions_df
            # Find longest sequence
            dash_parameters.longest_sequence = get_longest_sequence_dataframe(gb_df)

            # Store the rest of information provided by the user for plotting
            dash_parameters.alignments_position = align_plot_state
            dash_parameters.identity_color = color_scale_state
            dash_parameters.colorscale_vmin = range_slider_state[0] / 100
            dash_parameters.colorscale_vmax = range_slider_state[1] / 100
            dash_parameters.set_colorscale_to_extreme_homologies = (
                is_set_to_extreme_homologies
            )
            dash_parameters.annotate_sequences = annotate_sequences_state
            dash_parameters.annotate_genes = annotate_genes_state
            dash_parameters.annotate_genes_with = use_genes_info_from_state
            # For now set it to "straight". Change it after fixing the bezier storage data
            dash_parameters.straight_homology_regions = "straight"
            dash_parameters.minimum_homology_length = minimum_homology_length_state
            dash_parameters.add_scale_bar = scale_bar_state
            dash_parameters.selected_traces = []

            # Set separation of sequences in the Y axis
            dash_parameters.y_separation = 10

            # Make figure
            fig = plt.make_figure(dash_parameters)

            fig.update_layout(clickmode="event+select")
            print("figure is displayed")
            return fig, None, False

        # ==== Erase plot ============================================================== #
        if button_id == "erase-button":
            # Reset all attributes of dash_parameters.
            dash_parameters.reset()
            # Return an empty figure, None for clickdata, and False for skeleton
            return {}, None, False

        # ============================================================================== #
        #                           EDIT TAB -> colors                                   #
        # ============================================================================== #
        # = Change homology color and colorscale bar legend = #
        if button_id == "change-homology-color-button":
            # If color_scale is different, update dash_paremeters.identity_color
            if color_scale_state != dash_parameters.identity_color:
                dash_parameters.identity_color = color_scale_state
            # Change homology color traces
            fig = plt.change_homoloy_color_traces(
                figure=figure_state,
                colorscale_name=color_scale_state,
                vmin_truncate=range_slider_state[0] / 100,
                vmax_truncate=range_slider_state[1] / 100,
                set_colorscale_to_extreme_homologies=is_set_to_extreme_homologies,
                lowest_homology=dash_parameters.lowest_identity,
                highest_homology=dash_parameters.highest_identity,
            )
            # Remove old colorscale bar legend
            fig = plt.remove_traces_by_name(fig, "colorbar legend")
            # Convert the fig dictionary return by remove_traces_by_name into a Figure
            # object
            fig = Figure(data=fig["data"], layout=fig["layout"])
            # Make new colorscale bar legend
            fig = plt.plot_colorbar_legend(
                fig=fig,
                colorscale=plt.get_truncated_colorscale(
                    colorscale_name=color_scale_state,
                    vmin=range_slider_state[0] / 100,
                    vmax=range_slider_state[1] / 100,
                ),
                min_value=dash_parameters.lowest_identity,
                max_value=dash_parameters.highest_identity,
                set_colorscale_to_extreme_homologies=is_set_to_extreme_homologies,
            )
            return fig, None, False

        # = Change color of selected traces = #
        if button_id == "change-gene-color-button":
            curve_numbers = set(dash_parameters.selected_traces)
            # Iterate over selected curve numbers and change color
            for curve_number in curve_numbers:
                figure_state["data"][curve_number]["fillcolor"] = color_input_state
                figure_state["data"][curve_number]["line"]["color"] = color_input_state
                figure_state["data"][curve_number]["line"]["width"] = 1
            # Empty "selected_traces" list.
            dash_parameters.selected_traces.clear()
            return figure_state, None, False

        # = Select traces for changing color = #
        if (
            (active_tab == "tab-edit")
            and (click_data is not None)
            and select_button_state
        ):
            # Get curve_number (selected trace)
            curve_number = click_data["points"][0]["curveNumber"]
            # If curve_number already in "selected_traces", remove it from the list and
            # restore trace to its previous state; this creates the effect of deselecting.
            if curve_number in dash_parameters.selected_traces:
                dash_parameters.selected_traces.remove(curve_number)
                fillcolor = figure_state["data"][curve_number]["fillcolor"]
                figure_state["data"][curve_number]["line"]["color"] = fillcolor
                figure_state["data"][curve_number]["line"]["width"] = 1
                return figure_state, None, False
            # Save the curve number in "selected_traces" for future modification
            dash_parameters.selected_traces.append(curve_number)
            # Make selection effect by changing line color of selected trace
            fig = plt.make_selection_effect(figure_state, curve_number)
            return fig, None, False

        # ==== Align sequences in the plot ============================================= #
        if button_id == "update-align-sequences-button":
            if align_plot_state != dash_parameters.alignments_position:
                # Change the value of dash_parameters -> alignments_position
                dash_parameters.alignments_position = align_plot_state
                # Make figure and get lowest and highest identities
                fig = plt.make_figure(dash_parameters)
                return fig, None, False

        # ============================================================================== #
        #                                  VIEW TAB                                      #
        # ============================================================================== #
        if figure_state and button_id == "update-annotations":
            # Convert the figure_state dictionary into a Figure object
            fig = Figure(data=figure_state["data"], layout=figure_state["layout"])

            # ==== Genes Annotations =================================================== #
            # check if user changed use_genes_info_from_state and if user wants to
            # to annotate genes
            if (
                use_genes_info_from_state != dash_parameters.annotate_genes_with
                and annotate_genes_state != "no"
            ):
                # Update dash_parameters.
                dash_parameters.annotate_genes_with = use_genes_info_from_state
                # Remove any gene annotations
                fig = plt.remove_annotations_by_name(fig, "Gene annotation:")
                # Annotate with the new parameter
                fig = plt.annotate_genes(fig, dash_parameters)
            # check if value of annotate_genes_state is different in dash_parameters
            if annotate_genes_state != dash_parameters.annotate_genes:
                # change value of dash_parameters -> annotate_genes
                dash_parameters.annotate_genes = annotate_genes_state
                # Remove any gene annotations
                fig = plt.remove_annotations_by_name(fig, "Gene annotation:")
                # If asked add new annotations
                if annotate_genes_state != "no":
                    fig = plt.annotate_genes(fig, dash_parameters)
            # ==== DNA Sequences Annotations =========================================== #
            # check if value of annotate_sequences_state is different in dash_parameters
            if annotate_sequences_state != dash_parameters.annotate_sequences:
                # Change value of dash_parameters -> annotate_sequences
                dash_parameters.annotate_sequences = annotate_sequences_state
                # Remove any dna sequence annotations
                fig = plt.remove_annotations_by_name(fig, "Sequence annotation:")
                # If annotate_sequences_state is not "no" add annotations.
                if annotate_sequences_state != "no":
                    fig = plt.annotate_dna_sequences(
                        fig=fig,
                        gb_records=dash_parameters.gb_df,
                        longest_sequence=dash_parameters.longest_sequence,
                        number_gb_records=dash_parameters.number_gb_records,
                        annotate_with=dash_parameters.annotate_sequences,
                        y_separation=dash_parameters.y_separation,
                    )
            # ==== Toggle scale bar ==================================================== #
            # check if value of scale_bar_state is different in dash_parameters
            if scale_bar_state != dash_parameters.add_scale_bar:
                # change value of dash_parameters -> add_cale_bar
                dash_parameters.add_scale_bar = scale_bar_state
                # toggle scale bar
                fig = plt.toggle_scale_bar(
                    fig, True if scale_bar_state == "yes" else False
                )
            # ==== Minimum Homology Length ============================================= #
            # check if minimum homology length is different from dash_parameters
            if minimum_homology_length_state != dash_parameters.minimum_homology_length:
                # change value of dash_parameters -> minimum_homology_length
                dash_parameters.minimum_homology_length = minimum_homology_length_state
                # Update homology regions.
                fig = plt.hide_homology(fig, int(minimum_homology_length_state))

            return fig, None, False

        return figure_state, None, False

    # ==== activate update buttons only when there is a figure ========================= #
    @app.callback(
        [
            Output("erase-button", "disabled"),
            Output("update-annotations", "disabled"),
            Output("change-gene-color-button", "disabled"),
            Output("change-homology-color-button", "disabled"),
            Output("select-change-color-button", "disabled"),
            Output("update-align-sequences-button", "disabled"),
        ],
        Input("plot", "figure"),
    )
    def toggle_update_buttons(figure) -> bool:
        if figure and figure.get("data", []):
            return [False] * 6
        return [True] * 6

    # ==== activate Draw button when files in upload table ============================= #
    @app.callback(
        Output("draw-button", "disabled"),
        Input("files-table", "rowData"),
    )
    def toggle_draw_button(row_data) -> bool:
        return False if row_data else True

    # ==== activate Select button ====================================================== #
    @app.callback(
        [
            Output("select-change-color-button", "variant"),
            Output("select-button-state-store", "data"),
        ],
        Input("select-change-color-button", "n_clicks"),
        State("select-button-state-store", "data"),
    )
    def toggle_select_button(n_clicks, is_active):
        if n_clicks:
            # Toggle the active state on click
            is_active = not is_active

        # Set button style based on the active state
        if is_active:
            button_style = "filled"
        else:
            button_style = "outline"
        return button_style, is_active

    # ==== toggle between set colorscale buttons ======================================= #
    @app.callback(
        [
            Output("extreme-homologies-button", "variant"),
            Output("extreme-homologies-button", "style"),
            Output("truncate-colorscale-button", "variant"),
            Output("truncate-colorscale-button", "style"),
            Output("is_set_to_extreme_homologies", "data"),
        ],
        [
            Input("extreme-homologies-button", "n_clicks"),
            Input("truncate-colorscale-button", "n_clicks"),
        ],
    )
    def toggle_colorscale_buttons(extreme_clicks, truncate_clicks):
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

    # ==== update the color scale display ============================================== #
    @app.callback(
        Output("color-scale-display", "figure"),
        Input("color-scale", "value"),
    )
    def update_color_scale(value):
        return plt.create_color_line(value.capitalize())

    # ==== reset app (page) ============================================================ #
    @app.callback(
        Output("url", "href"),
        Input("reset-button", "n_clicks"),
        prevent_initial_call=True,
    )
    def reset_page(n_clicks):
        if n_clicks:
            # Return the current URL to trigger a reload
            return "/"

    # ==== save file =================================================================== #
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
        n_clicks,
        figure,
        figure_format,
        scale,
        width,
        height,
    ):
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

    # ↓↓↓↓ CHECKING IF TAB WAS CLOSED TO KILL SERVER ↓↓↓↓ #
    @app.server.route("/heartbeat", methods=["POST"])
    def heartbeat():
        """Receive heartbeat pings from the client."""
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

    def monitor_heartbeats():
        """Monitor the heartbeats and detect tab closure."""
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

    if not heartbeat_parameters.heartbeat_monitor_started:
        heartbeat_parameters.heartbeat_monitor_started = True
        print("Initiating heartbeat_monitor_started")
        # Start the monitoring thread
        threading.Thread(target=monitor_heartbeats, daemon=True).start()

    # Endpoint to shut down the server
    @app.server.route("/shutdown", methods=["POST"])
    def shutdown_server():
        os.kill(os.getpid(), signal.SIGINT)  # Send a signal to terminate the process
        print("Server shutting down...")
        return "Server shutting down...", 200

    # ↑↑↑↑ CHECKING IF TAB WAS CLOSED TO KILL SERVER ↑↑↑↑ #

    return app
