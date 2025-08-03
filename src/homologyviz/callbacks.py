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
import binascii
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
from pandas import DataFrame

import dash
from dash import Input, Output, State
from plotly.graph_objects import Figure

from homologyviz.parameters import PlotParameters
from homologyviz import plotter as plt
from homologyviz.gb_files_manipulation import get_longest_sequence_dataframe

# import logging
# logging.basicConfig(level=logging.INFO)
# logger = logging.getLogger(__name__)


# TODO: If we change a color in the edit tab, the changes are return to the
# original colors when performing any change in the view tab.


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


def save_uploaded_file(
    file_name: str, content: str, temp_folder_path: Path
) -> str | None:
    """Decode the content and write it to a temporary file.

    Returns the file path as a string if successful, otherwise returns None.
    """
    try:
        # Ensure content is properly formatted
        if ";base64," not in content:
            raise ValueError("Content is not base64-encoded or improperly formatted.")

        # Decode content
        data = content.split(";base64,")[1]
        decoded_data = base64.b64decode(data)

        # Save uploaded file
        output_path = temp_folder_path / file_name
        with open(output_path, "wb") as f:
            f.write(decoded_data)

        # Dash doesn't like Path; hence, we need to cast Path to str.
        return str(output_path)

    except (IndexError, ValueError, binascii.Error) as e:
        print(f"Failed to decode and save uplaoded file: {e}")
        return None


def handle_plot_button_click(
    dash_parameters: PlotParameters,
    virtual: list[dict[str, str]],
    tmp_directory_path: Path,
    align_plot_state: str,
    color_scale_state: str,
    range_slider_state: list[int, int],
    is_set_to_extreme_homologies: bool,
    annotate_sequences_state: str,
    annotate_genes_state: str,
    use_genes_info_from_state: str,
    homology_style_state: str,
    minimum_homology_length_state: int,
    scale_bar_state: str,
    title_input_state: str,
) -> tuple[Figure, None, bool]:
    """
    Perform BLASTn alignments and generate a homology plot for Dash.

    This function prepares alignments data from input files, sets plotting parameters, and
    generates a Plotly figure representing sequence alignments and homologies.

    Parameters
    ----------
    dash_parameters : PlotParameters
        Object holding all plotting configuration and data.
    virtual : list[dict[str, str]]
        Metadata for uploaded files, including file names and file paths.
    tmp_directory_path : Path
        Path to the temporary folder for storing alignments results.
    align_plot_state : str
        Layout preference for positioning the alignments in the plot (e.g. "left",
        "center", "right").
    color_scale_state : str
        Name of the color scale used to represent homology identity levels.
    range_slider_state : list[int, int]
        Percent identity range (e.g. [50, 100]) selected by the used to define color
        scalling.
    is_set_to_extreme_homologies : bool
        Whether to stretch the color scale to the min/max homology identity values in the
        data.
    annotate_sequences_state : str
        Whether and how to annotate sequence names.
    annotate_genes_state : str
        Whether gene features shold be annotated.
    use_genes_info_from_state : str
        Indicate source for genes annotations (e.g. "CDS product", "CDS gene").
    homology_style_state : str
        Whether the connections between homology regions are straight or curved (Bezier).
    minimum_homology_length_state : int
        Minimum length of homology region to be displayed.
    scale_bar_state : str
        Whether to include a scale bar in the plot.
    title_input_state : str
        String holding the figure's title.

    Returns
    -------
    fig : plotly.graph_objects.Figure
    None
        Placeholder to reset 'clickData' in Dash callbacks.
    bool
        A flag (`False`) to indicate that the dmc.Skeleton loading component should be
        hidden.
    """
    print("clicking plot-button")
    print(f"tmp directory path: {tmp_directory_path}")
    dash_parameters.draw_from_button = "plot-button"

    input_files = [Path(row["filepath"]) for row in virtual]
    dash_parameters.input_files = input_files
    dash_parameters.output_folder = tmp_directory_path
    dash_parameters.number_gb_records = len(input_files)

    gb_df, cds_df, alignments_df, regions_df = plt.make_alignments(
        input_files, tmp_directory_path
    )
    dash_parameters.gb_df = gb_df
    dash_parameters.cds_df = cds_df
    dash_parameters.alignments_df = alignments_df
    dash_parameters.alignments_regions_df = regions_df
    dash_parameters.longest_sequence = get_longest_sequence_dataframe(gb_df)

    dash_parameters.alignments_position = align_plot_state
    dash_parameters.identity_color = color_scale_state
    dash_parameters.colorscale_vmin = range_slider_state[0] / 100
    dash_parameters.colorscale_vmax = range_slider_state[1] / 100
    dash_parameters.set_colorscale_to_extreme_homologies = is_set_to_extreme_homologies
    dash_parameters.annotate_sequences = annotate_sequences_state
    dash_parameters.annotate_genes = annotate_genes_state
    dash_parameters.annotate_genes_with = use_genes_info_from_state
    dash_parameters.style_homology_regions = homology_style_state
    dash_parameters.minimum_homology_length = minimum_homology_length_state
    dash_parameters.add_scale_bar = scale_bar_state
    dash_parameters.selected_traces = []
    dash_parameters.y_separation = 10
    dash_parameters.plot_title = title_input_state

    fig = plt.make_figure(dash_parameters)
    fig.update_layout(clickmode="event+select")
    print("figure is displayed")
    return fig, None, False


def check_plot_parameters_for_update_homologies(
    dash_parameters: PlotParameters,
    color_scale_state: str,
    range_slider_state: list[int, int],
    is_set_to_extreme_homologies: bool,
) -> bool:
    """Check if plotting parameters for homology regions provided by the user are the same
    as the current stored values in the PlotParameters Object

    If values are the same return False. Otherwise, update PlotParameters values and
    return True

    Parameters
    ----------
    dash_parameters : PlotParameters
        Object holding all plotting configuration and data
    color_scale_state : str
        Name of the color scale used to represent homology identity levels
    range_slider_state : list[int, int]
        Percent identity range (e.g. [50, 100]) selected by the used to define color
        scalling.
    is_set_to_extreme_homologies : bool
        Whether to stretch the color scale to the min/max homology identity values in the
        data

    Returns
    -------
    bool
        A flag to indicate if values are the same (`True`) or not (`False`)
    """
    vmin = range_slider_state[0] / 100
    vmax = range_slider_state[1] / 100
    # if all values are the same as in dash_parameter, then return False
    if (
        color_scale_state == dash_parameters.identity_color
        and vmin == dash_parameters.colorscale_vmin
        and vmax == dash_parameters.colorscale_vmax
        and is_set_to_extreme_homologies
        == dash_parameters.set_colorscale_to_extreme_homologies
    ):
        return False
    else:
        dash_parameters.identity_color = color_scale_state
        dash_parameters.colorscale_vmin = vmin
        dash_parameters.colorscale_vmax = vmax
        dash_parameters.set_colorscale_to_extreme_homologies = (
            is_set_to_extreme_homologies
        )
        return True


def handle_update_homologies_click(
    figure_state: dict,
    dash_parameters: PlotParameters,
    color_scale_state: str,
    range_slider_state: list[int, int],
    is_set_to_extreme_homologies: bool,
) -> tuple[Figure, None, bool]:
    """Update the homology trace colors and regenerate the colorscale bar legend.

    This function updates the figure based on a new colorscale or identity range, and
    regenerates the corresponding colorbar legend for homology visualization.

    Parameters
    ----------
    figure_state : dict
        Dictionary representing the current Plotly figure, retrieved from dcc.Graph via
        Dash State
    dash_parameters : PlotParameters
        Object holding all plotting configuration and data
    color_scale_state : str
        Name of the color scale used to represent homology identity levels
    range_slider_state : list[int, int]
        Percent identity range (e.g. [50, 100]) selected by the used to define color
        scalling.
    is_set_to_extreme_homologies : bool
        Whether to stretch the color scale to the min/max homology identity values in the
        data

    Returns
    -------
    fig : plotly.graph_objects.Figure
    None
        Placeholder to reset 'clickData' in Dash callbacks
    bool
        A flag (`False`) to indicate that the dmc.Skeleton loading component should be
        hidden
    """
    if not check_plot_parameters_for_update_homologies(
        dash_parameters,
        color_scale_state,
        range_slider_state,
        is_set_to_extreme_homologies,
    ):
        return figure_state, None, False

    fig = plt.change_homology_color(
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

    # Convert the fig dictionary return by remove_traces_by_name into a Figure object
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


def change_color_cell_cds_dataframe(
    cds_dataframe: DataFrame,
    file_number: int,
    cds_number: int,
    new_color: str,
) -> None:
    """
    Update the color value for a specific coding sequence in the DataFrame.

    This function locates the row in the `cds_dataframe` corresponding to the given
    `file_number` and `cds_number`, and updates the value in the "color" column to the
    specified `new_color`. The function modifies the DataFrame in place.

    Parameters
    ----------
    cds_dataframe : pandas.DataFrame
        The DataFrame containig coding sequence data, including columns "file_number",
        "cds_number", and "color".
    file_number : int
        The file identifier used to locate the target row.
    cds_number : int
        The CDS idenfifier used to locate the target row.
    new_color : str
        The new color value to assign, typically in hexadecimal format.

    Returns
    -------
    None
        The input DataFrame is modified in place.
    """
    # change the value of color in DataFrame
    cds_dataframe.loc[
        (cds_dataframe["file_number"] == file_number)
        & (cds_dataframe["cds_number"] == cds_number),
        "color",
    ] = new_color


def handle_change_color_click(
    figure_state: dict, dash_parameters: PlotParameters, color_input_state: str
) -> tuple[Figure, None, bool]:
    """Change color of selected traces.

    Applies the chosen color to all traces currently marked as selected in
    `dash_parameters.selected_traces`, then clears the selection list.

    Parameters
    ----------
    figure_state : dict
        Dictionary representing the current Plotly figure, retrieved from dcc.Graph via
        Dash State.
    dash_parameters : PlotParameters
        Object holding all plotting configuration and data, including selected traces.
    color_input_state : str
        Hex color code (e.g., "#FF0000) selected by the user to apply to the selected
        traces

    Returns
    -------
    fig : plotly.graph_objects.Figure
        The updated Plotly figure with modified trace colors.
    None
        Placeholder to reset 'clickData' in Dash callbacks
    bool
        A flag (`False`) to indicate that the dmc.Skeleton loading component should be
        hidden
    """
    # Iterate over selected curve numbers and change color
    for curve_number in set(dash_parameters.selected_traces):
        customdata = figure_state["data"][curve_number]["customdata"]
        # Change the value of color in the DataFrame
        change_color_cell_cds_dataframe(
            dash_parameters.cds_df, customdata[0], customdata[1], color_input_state
        )
        figure_state["data"][curve_number]["fillcolor"] = color_input_state
        figure_state["data"][curve_number]["line"]["color"] = color_input_state
        figure_state["data"][curve_number]["line"]["width"] = 1
    # Empty "selected_traces" list.
    dash_parameters.selected_traces.clear()
    return figure_state, None, False


def handle_update_title_click(
    figure_state: dict,
    dash_parameters: PlotParameters,
    title_input_state: str,
) -> tuple[Figure, None, bool]:
    """
    Handle update title button click and update the Plotly figure accordingly.

    If the provided title is unchanged or only whitespace, the figure title is removed.
    Otherwise, the new title is set and centered.

    Parameters
    ----------
    figure_state : dict
        Current Plotly figure state, retrieved from dcc.Graph via Dash callback.
    dash_parameters : PlotParameters
        Object holding plotting configuration and metadata. This function updates its
        `plot_title` attribute.
    title_input_state : str
        New title input from the user.

    Returns
    -------
    fig : plotly.graph_objects.Figure
        The updated Plotly figure.
    None
        Placeholder to reset 'clickData' in Dash callbacks.
    bool
        `False` to hide the dmc.Skeleton loading indicator.
    """
    # Convert figure_state dictionary into a Figure object
    fig = Figure(data=figure_state["data"], layout=figure_state["layout"])
    # If title is the same, do nothing and return
    if title_input_state == dash_parameters.plot_title:
        return fig, None, False

    # update dash_parametres.plot_title
    dash_parameters.plot_title = title_input_state
    fig = plt.add_or_remove_title(fig, title_input_state)

    return fig, None, False


def handle_select_traces_click(
    figure_state: dict,
    dash_parameters: PlotParameters,
    click_data: dict,
) -> tuple[Figure, None, bool]:
    """
    Handle click events on traces to toggle selection and update the figure.

    This function stores the selected trace index from `click_data`, applies a visual
    selection effect (e.g., line color/width change), and allows toggling the selection
    on repeated clicks.

    Parameters
    ----------
    figure_state : dict
        Dictionary representing the current Plotly figure, retrieved from dcc.Graph via
        Dash State.
    dash_parameters : PlotParameters
        Object holding all plotting configuration and data.
    click_data : dict
        Dictionary representing data about the clicked point, as returned by Dash's
        `clickData`. Must contain a "points" list with "curveNumber" to identify the
        clicked trace.

    Returns
    -------
    fig : plotly.graph_objects.Figure
    None
        Placeholder to reset 'clickData' in Dash callbacks.
    bool
        A flag (`False`) to indicate that the dmc.Skeleton loading component should be
        hidden.
    """
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


def align_plot(
    figure_state: dict, dash_parameters: PlotParameters, align_plot_state: str
) -> Figure:
    """Align the homology plot to the left, center, or right based on user preference.

    If the selected alignment differs from the current one stored in `dash_parameters`,
    a new figure is generated. Otherwise, the existing figure state is converted back
    to a Plotly Figure object.

    Parameters
    ----------
    figure_state : dict
        Dictionary representing the current Plotly figure, retrieved from dcc.Graph via
        Dash State.
    dash_parameters : PlotParameters
        Object holding all plotting configuration and data.
    align_plot_state : str
        Layout preference for positioning the alignments in the plot (e.g. "left",
        "center", "right").

    Returns
    -------
    fig : plotly.graph_objects.Figure
        The updated or restored Plotly figure.
    """
    # ==== Align sequences in the plot ========================================= #
    # Check if user wants to change the plot position
    if align_plot_state != dash_parameters.alignments_position:
        # Change the value of dash_parameters -> alignments_position
        dash_parameters.alignments_position = align_plot_state
        # Make figure with new plot position
        fig = plt.make_figure(dash_parameters)
    # Otherwise, convert figure_state dictionary into a Figure object
    else:
        fig = Figure(data=figure_state["data"], layout=figure_state["layout"])
    return fig


def update_homology_regions(
    figure_state: dict,
    dash_parameters: PlotParameters,
    align_plot_state: str,
    homology_style_state: str,
) -> Figure:
    """
    Update the homology region style and alignment position in the plot if user
    preferences change.

    This function checks whether the user has requested a change to the style of homology
    region shading (e.g., Bezier vs. straight) or to the alignment position of the
    homology region plot (e.g., left, center, or right). If either preference has changed,
    the figure is redrawn using the updated parameters. Otherwise, the figure is
    reconstructed from the current figure state.

    Parameters
    ----------
    figure_state : dict
        A dictionary containing the current Plotly figure's 'data' and 'layout'.
    dash_parameters : PlotParameters
        An object storing the current state of plotting parameters, including alignment
        position and homology style.
    align_plot_state : str
        The desired position of the homology alignment plot in the figure ('left',
        'center', or 'right')
    homology_style_state : str
        The desired visual style of homology regions ('Bezier' or 'straight').

    Returns
    -------
    plotly.graph_objects.Figure
        The updated Plotly figure with the appropriate homology style and alignment
        position.
    """
    print(f"homology style state: {homology_style_state}")
    # Check if user wants to change the plot location and homology style
    if (
        align_plot_state != dash_parameters.alignments_position
        or homology_style_state != dash_parameters.style_homology_regions
    ):
        # Update dash_parameters
        dash_parameters.alignments_position = align_plot_state
        dash_parameters.style_homology_regions = homology_style_state
        # Make figure with new plot loation and style
        fig = plt.make_figure(dash_parameters)

    # Otherwise, convert figure_state dictionary into a Figure object
    else:
        fig = Figure(data=figure_state["data"], layout=figure_state["layout"])

    return fig


def update_genes_annotations(
    fig: Figure,
    dash_parameters: PlotParameters,
    use_genes_info_from_state: str,
    annotate_genes_state: str,
) -> Figure:
    """
    Update genes annotations in the plot based on user preferences.

    If the user changes either the annotation source (e.g., product or gene) or the
    visibility of annotation (top, bottom, both, or none), the function updates the
    figure accordingly. If no change are neede, the input figure is returned unchanged.

    Parameters
    ----------
    fig : plotly.graph_objects.Figure.
        The current Plotly figure to update.
    dash_parameters : PlotParameters
        Object holding all plotting configuration and data.
    use_genes_info_from_state : str
        Source of gene annotation labels (e.g., "CDS product", "CDS gene").
    annotate_genes_state : str
        Desired gene annotation display setting (e.g., "top", "bottom", "both", "no").

    Returns
    -------
    fig : plotly.graph_objects.Figure.
        The updated or original Plotly figure, depending on whether changes are needed.
    """

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
    return fig


def update_dna_sequences_annotations(
    fig: Figure,
    dash_parameters: PlotParameters,
    annotate_sequences_state: str,
) -> Figure:
    """
    Update DNA sequences annotations in the plot based on user preferences.

    If the user changes either the annotation source (e.g., accession number, sequence
    name, or file name) or the visibility of annotation (no), the function updates the
    figure accordingly. If no changes are neede, the input figure is returned unchanged.

    Parameters
    ----------
    fig : plotly.graph_objects.Figure.
        The current Plotly figure to update.
    dash_parameters : PlotParameters
        Object holding all plotting configuration and data.
    annotate_sequences_state : str
        Desired sequence annotation display setting (e.g., "accession", "name", "fname",
        "no").

    Returns
    -------
    fig : plotly.graph_objects.Figure
        The updated or original Plotly figure, depending on whether changes are needed.
    """
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
    return fig


def update_gb_dataframe_custom_name(table: list[dict], gb_df: DataFrame):
    """
    Update the 'custom_name' column of a DataFrame using values from a table of
    dictionaries.

    Parameters
    ----------
    table : list of dict
        A list where each dictionary contain a 'custom_name' key.
    gb_df : pandas.DataFrame
        The DataFrame whose 'custom_name' column will be updated.

    Returns
    -------
    pandas.DataFrame
        The updated DataFrame with the 'custom_name' column set from `table`.
    """
    custom_names = [row["custom_name"] for row in table]
    gb_df["custom_name"] = custom_names
    return gb_df


def update_dna_sequence_annotations(
    fig: Figure,
    dash_parameters: PlotParameters,
    annotation_column_choice_state: str,
    table: list[dict],
):
    """
    Update DNA sequence annotations in the Plotly figure based on user-selected
    annotations options.

    This function is triggered when the user edits the `Annotations/Sequences` dropdown in
    the `Edit` tab of the app. It removes existing DNA sequence annotations and re-adds
    them based on the current user selection. It also updates the 'custom_name' column in
    the GenBank DataFrame.

    Parameters
    ----------
    fig : plotly.graph_objects.Figure
        The Plotly figure containing the DNA sequence tracks and annotations.
    dash_parameters : PlotParameters
        An object containing state variables for plotting, such as GenBank data and
        annotation preferences.
    annotation_column_choice_state : str
        The column name selected by the user to annotate sequences with. Use "no" to
        disable annotations.
    table : list of dict
        A list of dictionaries representing the rows in the editable sequence table.
        Each dictionary must include a "custom_name" field.

    Returns
    -------
    plotly.graph_objects.Figure
        The updated figure with DNA sequence annotations added or removed based on user
        input.
    """
    # if annotation_column_choice_state != dash_parameters.annotate_sequences:
    # Update dash_parameters
    dash_parameters.annotate_sequences = annotation_column_choice_state
    # Update custom name field of GenBank DataFrame
    update_gb_dataframe_custom_name(table, dash_parameters.gb_df)
    # Remove any DNA sequence annotations
    fig = plt.remove_annotations_by_name(fig, "Sequence annotation:")
    # If user selected a different value than "no" add annotations.
    if annotation_column_choice_state != "no":
        fig = plt.annotate_dna_sequences(
            fig=fig,
            gb_records=dash_parameters.gb_df,
            longest_sequence=dash_parameters.longest_sequence,
            number_gb_records=dash_parameters.number_gb_records,
            annotate_with=dash_parameters.annotate_sequences,
            y_separation=dash_parameters.y_separation,
        )
    return fig


def update_scale_bar(
    fig: Figure,
    dash_parameters: PlotParameters,
    scale_bar_state: str,
) -> Figure:
    """
    Update the visibility of the scale bar in the plot based on user preferences.

    If the user changes the scale bar setting, the function updates the figure
    accordingly. If no changes are needed, the input figure is returned unchanged.

    Parameters
    ----------
    fig : plotly.graph_objects.Figure.
        The current Plotly figure to update.
    dash_parameters : PlotParameters
        Object holding all plotting configuration and data.
    scale_bar_state : str
        Desired scale bar annotation display setting ("yes" to show, "no" to hide).

    Returns
    -------
    fig : plotly.graph_objects.Figure
        The updated or original Plotly figure, depending on whether changes are needed.
    """
    # check if value of scale_bar_state is different in dash_parameters
    if scale_bar_state != dash_parameters.add_scale_bar:
        # change value of dash_parameters -> add_cale_bar
        dash_parameters.add_scale_bar = scale_bar_state
        # toggle scale bar
        fig = plt.toggle_scale_bar(fig, True if scale_bar_state == "yes" else False)
    return fig


def update_minimum_homology_length(
    fig: Figure,
    dash_parameters: PlotParameters,
    minimum_homology_length_state: int,
) -> Figure:
    """
    Update minimum homology length displayed in the plot based on user preferences.

    If the user changes the minimum homology length setting, the function updates the
    figure accordingly by hidding homology regions shorter than the specified length.
    If no changes are needed, the input figure is returned unchanged.

    Parameters
    ----------
    fig : plotly.graph_objects.Figure.
        The current Plotly figure to update.
    dash_parameters : PlotParameters
        Object holding all plotting configuration and data.
    minimum_homology_length_state : int
        The new minimum homology length to display in the plot.

    Returns
    -------
    fig : plotly.graph_objects.Figure
        The updated or original Plotly figure, depending on whether changes are needed.
    """
    # check if minimum homology length is different from dash_parameters
    if minimum_homology_length_state != dash_parameters.minimum_homology_length:
        # change value of dash_parameters -> minimum_homology_length
        dash_parameters.minimum_homology_length = minimum_homology_length_state
        # Update homology regions.
        fig = plt.hide_homology(fig, int(minimum_homology_length_state))
    return fig


def handle_update_view_click(
    figure_state: dict,
    dash_parameters: PlotParameters,
    align_plot_state: str,
    homology_style_state: str,
    use_genes_info_from_state: str,
    annotate_genes_state: str,
    annotate_sequences_state: str,
    scale_bar_state: str,
    minimum_homology_length_state: int,
) -> tuple[Figure, None, bool]:
    """
    Handle the 'update view' button click event.

    This function updates the current figure layout and annotations based on user
    preferences, including alignment positioning, gene/sequence annotations, scale bar
    visibility, and minimum homology length.

    Parameters
    ----------
    figure_state : dict
        Dictionary representing the current Plotly figure, retrieved from dcc.Graph via
        Dash State.
    dash_parameters : PlotParameters
        Object holding all plotting configuration and data.
    align_plot_state : str
        Layout preference for positioning the alignments in the plot (e.g. "left",
        "center", "right").
    homology_style_state : str
        Whether the connections between homology regions are straight or curved (Bezier).
    use_genes_info_from_state : str
        Indicate source for genes annotations (e.g. "CDS product", "CDS gene").
    annotate_genes_state : str
        whether gene features shold be annotated.
    annotate_sequences_state : str
        Whether and how to annotate sequence names.
    scale_bar_state : str
        Whether to display the scale bar ("yes" or "no").
    minimum_homology_length_state : int
        Minimum length (in bp) of homology region to be displayed.

    Returns
    -------
    fig : plotly.graph_objects.Figure
        The updated Plotly figure with applied user preferences.
    None
        Placeholder to reset 'clickData' in Dash callbacks.
    bool
        A flag (`False`) to indicate that the dmc.Skeleton loading component should be
        hidden.
    """
    # Update homology connector style and/or positio of the alignment in the figure
    fig = update_homology_regions(
        figure_state,
        dash_parameters,
        align_plot_state,
        homology_style_state,
    )

    # Update genes annotations
    fig = update_genes_annotations(
        fig,
        dash_parameters,
        use_genes_info_from_state,
        annotate_genes_state,
    )

    # # Update DNA sequences annotations
    # fig = update_dna_sequences_annotations(
    #     fig, dash_parameters, annotate_sequences_state
    # )

    # Toggle scale bar
    fig = update_scale_bar(fig, dash_parameters, scale_bar_state)

    # Update the minimum homology displayed
    fig = update_minimum_homology_length(
        fig, dash_parameters, minimum_homology_length_state
    )

    return fig, None, False


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
            State("annotate-sequences", "value"),
            State("annotate-genes", "value"),
            State("scale-bar", "value"),
            State("use-genes-info-from", "value"),
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
        annotate_sequences_state: str,
        annotate_genes_state: str,
        scale_bar_state: str,
        use_genes_info_from_state: str,
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
        annotate_sequences_state : str
            Whether and how to annotate DNA sequences (e.g., "accession", "name", "no").
        annotate_genes_state : str
            Whether and where to annotate gene features (e.g., "top", "bottom", "no").
        scale_bar_state : str
            Whether to show the scale bar ("yes" or "no").
        use_genes_info_from_state : str
            Source of gene labels used for annotation (e.g., "CDS product", "CDS gene").
        title_input_state : str
            String holding the figure's title.

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
        # TODO: fix annotate_squences_state with annotation_column_choice_state
        # Use context to find the button that triggered the callback.
        ctx = dash.callback_context
        button_id = ctx.triggered[0]["prop_id"].split(".")[0]
        print(f"button_id: {button_id}")

        # = TAB MAIN =================================================================== #
        # Perform Alignments & Plot
        if (button_id == "plot-button") and virtual:
            return handle_plot_button_click(
                dash_parameters,
                virtual,
                tmp_directory_path,
                align_plot_state,
                color_scale_state,
                range_slider_state,
                is_set_to_extreme_homologies,
                annotate_sequences_state,
                annotate_genes_state,
                use_genes_info_from_state,
                homology_style_state,
                minimum_homology_length_state,
                scale_bar_state,
                title_input_state,
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
                figure_state,
                dash_parameters,
                align_plot_state,
                homology_style_state,
                use_genes_info_from_state,
                annotate_genes_state,
                annotate_sequences_state,
                scale_bar_state,
                minimum_homology_length_state,
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
                figure_state,
                dash_parameters,
                title_input_state,
            )

        # Update sequence annotations
        if button_id == "offcanvas-update-sequence-annotations-button":
            fig = Figure(data=figure_state["data"], layout=figure_state["layout"])
            fig = update_dna_sequence_annotations(
                fig,
                dash_parameters,
                annotation_column_choice_state,
                sequence_table_state,
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
            return [False] * 7
        return [True] * 7

    # = Activate the plot button only if there are files in the input table ============ #
    # TODO: make sure the button is actevited only if there are two or more files in the
    # table. Also, deactivate the button if the user erases all files.
    @app.callback(
        Output("plot-button", "disabled"),
        Input("files-table", "rowData"),
    )
    def toggle_plot_button(row_data: list[dict[str, str]] | None) -> bool:
        """
        Enable or disable the 'Plot' button based on the presence of uploaded files.

        This callback activates the 'Plot' button only when the upload table contains
        at least one file. If no files are uploaded, the button is disabled to prevent
        plotting without input data.

        Parameters
        ----------
        row_data : list[dict[str, str]] or None
            The current contents of the GenBank file upload table.

        Returns
        -------
        bool
            True if the button should be disabled, False otherwise.
        """
        return False if row_data else True

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
            # if table is not None:
            #     custom = [x["custom_name"] for x in table]
            #     print(f"value custom name: \n{custom}")
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
