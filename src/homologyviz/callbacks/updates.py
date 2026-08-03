from typing import Final, Literal

import pandas as pd
from pandas import DataFrame
from plotly.graph_objects import Figure

from homologyviz.parameters import PlotParameters
from homologyviz import plotter as plt


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


def update_homology_regions(
    *,
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
        # Make figure with new plot location and style
        fig = plt.make_figure(dash_parameters)

    # Otherwise, convert figure_state dictionary into a Figure object
    else:
        fig = Figure(data=figure_state["data"], layout=figure_state["layout"])

    return fig


def update_genes_annotations(
    *,
    fig: Figure,
    dash_parameters: PlotParameters,
    annotate_genes_from_state: str,
    annotate_genes_positions_state: str,
    table: list[dict],
) -> Figure:
    """
    Update genes annotations in the plot based on user preferences.

    If the user changes either the annotation source (e.g., gene, product, or custom) or
    the visibility of annotation (no, top, bottom, top-bottom, all-below or all-above),
    the function updates the figure accordingly. If no change are needed, the input figure
    is returned unchanged.

    Parameters
    ----------
    fig : plotly.graph_objects.Figure.
        The current Plotly figure to update.
    dash_parameters : PlotParameters
        Object holding all plotting configuration and data.
    annotate_genes_from_state : str
        Source of gene annotation labels (e.g., "CDS product", "CDS gene", "custom").
    annotate_genes_positions_state : str
        Desired gene annotation display setting (e.g., "no", "top", "bottom",
        "top-bottom").
    table : list of dict
        A list where each dictionary represents a row in the editable gene annotation
        table, containing a "custom_name" key for custom gene labels.

    Returns
    -------
    fig : plotly.graph_objects.Figure.
        The updated or original Plotly figure, depending on whether changes are needed.
    """
    # If user wants to remove gene annotations, update dash_parameters and remove any
    # existing gene annotations from the figure.
    if annotate_genes_positions_state == "no":
        # Update dash_parameters.
        dash_parameters.annotate_genes_positions = "no"
        # Remove any gene annotations
        fig = plt.remove_annotations_by_name(fig, "Gene annotation:")
        return fig
    # If user wants to add or change gene annotations, update dash_parameters and
    # re-annotate the figure accordingly. Additionally, update CDS DataFrame.
    dash_parameters.annotate_genes_positions = annotate_genes_positions_state
    dash_parameters.annotate_genes_from = annotate_genes_from_state
    # Update CDS DataFrame.
    dash_parameters.cds_df = update_cds_df(
        table=table,
        cds_df=dash_parameters.cds_df,
        header=annotate_genes_from_state,
    )
    # Remove any gene annotations
    fig = plt.remove_annotations_by_name(fig, "Gene annotation:")
    # Annotate with the new parameter
    fig = plt.annotate_genes(fig, dash_parameters)
    return fig


def update_dna_sequence_annotations(
    *,
    fig: Figure,
    dash_parameters: PlotParameters,
    annotation_column_choice_state: str,
    table: list[dict],
) -> Figure:
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
    fig : plotly.graph_objects.Figure
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


def update_gb_dataframe_custom_name(table: list[dict], gb_df: DataFrame) -> DataFrame:
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


AnnotationHeader = Literal["gene", "product", "custom_name"]

ANNOTATION_HEADERS: Final[set[str]] = {
    "gene",
    "product",
    "custom_name",
}

CDS_ID_COLUMNS: Final[tuple[str, ...]] = (
    "file_number",
    "accession",
    "start",
    "end",
    "strand",
)


def update_cds_df(
    table: list[dict],
    cds_df: DataFrame,
    header: AnnotationHeader,
) -> DataFrame:
    """
    Update an annotation column in a CDS DataFrame from an AG Grid table.

    The input table is validated to ensure it contains the required identifier columns and
    the selected annotation column. Rows are matched to the CDS DataFrame using the CDS
    identifier columns rather than their row positions, allowing the table to be sorted
    without affecting the correctness of the update.

    Parameters
    ----------
    table : list of dict
        Table data from AG Grid. Each dictionary must contain the CDS identifier columns
        and the annotation column specified by `header`.
    cds_df : pandas.DataFrame
        CDS DataFrame to update.
    header : {"gene", "product", "custom_name"}
        Annotation column to update.

    Returns
    -------
    pandas.DataFrame
        A copy of `cds_df` with the selected annotation column updated.

    Raises
    ------
    ValueError
        If `header` is invalid, if duplicate CDS identifiers are found, or if the table
        and DataFrame do not contain the same CDS identifiers.
    KeyError
        If the table or the CDS DataFrame is missing one or more required columns.

    Notes
    -----
    CDS rows are matched using the columns defined in `CDS_ID_COLUMNS`, allowing the input
    table to be reordered without assigning annotations to the wrong CDS.
    """
    if header not in ANNOTATION_HEADERS:
        raise ValueError(
            f"Invalid header {header!r}. "
            f"Expected one of {sorted(ANNOTATION_HEADERS)}."
        )

    # Make a DataFrame from the table list of dictionaries
    table_df = DataFrame(table)

    required_columns = {*CDS_ID_COLUMNS, header}

    # Check required columns exists (set arithmetic)
    missing_table_columns = required_columns - set(table_df.columns)
    if missing_table_columns:
        raise KeyError(
            "The table is missing required columns: " f"{sorted(missing_table_columns)}"
        )
    missing_df_columns = required_columns - set(cds_df.columns)
    if missing_df_columns:
        raise KeyError(
            "The CDS DataFrame is missing required columns: "
            f"{sorted(missing_df_columns)}"
        )

    # Check for duplicates in table_df and cds_df
    if table_df.duplicated(list(CDS_ID_COLUMNS)).any():
        raise ValueError("The table contains duplicated CDS identifiers.")
    if cds_df.duplicated(list(CDS_ID_COLUMNS)).any():
        raise ValueError("The CDS DataFrame contains duplicated CDS identifiers.")

    # Create tables with the indexes of table_df and cds_df for comparison.
    table_ids = pd.MultiIndex.from_frame(table_df.loc[:, list(CDS_ID_COLUMNS)])
    dataframe_ids = pd.MultiIndex.from_frame(cds_df.loc[:, list(CDS_ID_COLUMNS)])
    # Check if the indexes of table_df and cds_df are equal
    if not table_ids.equals(dataframe_ids):
        if set(table_ids) != set(dataframe_ids):
            raise ValueError("The table rows do not match the CDS DataFrame rows.")

    indexed_table = table_df.set_index(list(CDS_ID_COLUMNS))
    updates = indexed_table[header]

    updated_df = cds_df.copy()
    updated_index = pd.MultiIndex.from_frame(updated_df.loc[:, list(CDS_ID_COLUMNS)])
    # This will align the updates with the order of cds_df. This will ensure that the
    # updates are applied to the correct rows in the original DataFrame in the case the
    # table from the Ag-Grid is not in the same order as the original DataFrame.
    ordered_updates = updates.reindex(updated_index)
    updated_df[header] = ordered_updates.to_numpy()

    return updated_df


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
