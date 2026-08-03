from pathlib import Path

from plotly.graph_objects import Figure

from homologyviz.parameters import PlotParameters
from homologyviz.gb_files_manipulation import (
    get_longest_sequence_dataframe,
)
from homologyviz import plotter as plt
from homologyviz.callbacks.updates import (
    check_plot_parameters_for_update_homologies,
    update_homology_regions,
    update_scale_bar,
    update_minimum_homology_length,
    change_color_cell_cds_dataframe,
)


def handle_plot_button_click(
    *,
    dash_parameters: PlotParameters,
    virtual: list[dict[str, str]],
    tmp_directory_path: Path,
    align_plot_state: str,
    color_scale_state: str,
    range_slider_state: list[int, int],
    is_set_to_extreme_homologies: bool,
    annotation_column_choice_state: str,
    annotate_genes_from_state: str,
    annotate_genes_positions_state: str,
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
    annotation_column_choice_state : str
        Wheter and how to annotate sequence names.
    annotate_genes_from_state : str
        Whether gene representations should be annotated and from where to take the
        information (e.g. "gene", "product", or "custom").
    annotate_genes_positions_state : str
        Whether gene representations should be annotated at the top, bottom, top-bottom,
        all-above, or all-below.
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
    dash_parameters.annotate_sequences = annotation_column_choice_state
    dash_parameters.annotate_genes_from = annotate_genes_from_state
    dash_parameters.annotate_genes_positions = annotate_genes_positions_state
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


def handle_change_color_click(
    figure_state: dict,
    dash_parameters: PlotParameters,
    color_input_state: str,
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
    *,
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


def handle_update_view_click(
    *,
    figure_state: dict,
    dash_parameters: PlotParameters,
    align_plot_state: str,
    homology_style_state: str,
    scale_bar_state: str,
    minimum_homology_length_state: int,
) -> tuple[Figure, None, bool]:
    """
    Handle the 'update view' button click event.

    This function updates the current figure layout and annotations based on user
    preferences, including alignment positioning, scale bar visibility, and minimum
    homology length.

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
        figure_state=figure_state,
        dash_parameters=dash_parameters,
        align_plot_state=align_plot_state,
        homology_style_state=homology_style_state,
    )

    # Toggle scale bar
    fig = update_scale_bar(fig, dash_parameters, scale_bar_state)

    # Update the minimum homology displayed
    fig = update_minimum_homology_length(
        fig, dash_parameters, minimum_homology_length_state
    )

    return fig, None, False
