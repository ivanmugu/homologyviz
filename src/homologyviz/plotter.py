"""
Functions and classes to generate graphical representations of BLASTn alignments.

Homology Visualization (HomologyViz) uses GenBank files (.gb) to align sequences, plot
genes, and visualize homology regions. Gene features are extracted from the `CDS` section
of each GenBank file.

Gene color customization is supported by adding a `/Color` qualifier to a `CDS` feature.
For example:
    /Color="#00ff00"
will render the gene in green. If no color is provided, HomologyViz defaults to yellow.

License
-------
This file is part of HomologyViz
BSD 3-Clause License
Copyright (c) 2024, Ivan Muñoz Gutiérrez
"""

import math
from pathlib import Path
from pandas import DataFrame

import plotly.graph_objects as go
from plotly.graph_objects import Figure, Heatmap
import plotly.colors as colors
import numpy as np
import matplotlib.colors as mcolors

from homologyviz.parameters import PlotParameters
from homologyviz.arrow import Arrow
from homologyviz.rectangle_bezier import RectangleCurveHeight
from homologyviz import gb_files_manipulation as genbank
from homologyviz import miscellaneous as misc


# TODO: I have to check how to store the data of bezier lines in a dataframe.
# TODO: The remove_annotations_by_name needs checking because Plotly annotations can be
# stored as Layout.annotation.Annotation objects, not just dicts. We may get a
# TypeError on annotation["name"] so we might need to switch to annotation.name


def create_color_line(colors: str) -> Figure:
    """
    Create a continuous horizontal color gradient to display a selected colorscale.

    This function takes input from the dropdown menu in the `Edit` tab of the app, which
    lists all Plotly sequential colorscales (e.g., "Greys", "Greens", "GnBu", "Blues_r").
    It generates a slim heatmap that visually represents the full range of the selected
    colorscale, helping users preview how the color mapping will appear when applied
    to homology identity values.

    Parameters
    ----------
    colors : str
        A Plotly sequential colorscale name (e.g., "Greys", "Greens", "GnBu", and
        "Blues_r").

    Returns
    -------
    figure : plotly.graph_objects.Figure
        A Plotly figure showing a smooth horizontal color gradient without axes or labels.
    """
    # Create a large number of z values for a smooth gradient (e.g., 1000 values)
    z_values = np.linspace(0, 1, 1000).reshape(1, -1)  # 1 row of 1000 values

    # Create heatmap with fine z-values
    figure = Figure(
        Heatmap(
            z=z_values,  # Fine-grained z-values
            colorscale=colors,
            showscale=False,  # Disable the color bar
            xgap=0,
            ygap=0,
        )
    )

    figure.update_layout(
        xaxis=dict(visible=False),  # hide x-axis
        yaxis=dict(visible=False),  # Hide y-axis
        height=40,  # Adjust height to display the line clearly
        margin=dict(l=0, r=0, t=0, b=0),  # Remove margins around the figure
        plot_bgcolor="white",  # Set background color to white
    )

    return figure


def get_color_from_colorscale(value: float, colorscale_name: str = "Greys") -> str:
    """
    Retrieve an RGB color from a Plotly colorscale based on a normalized value.

    This function samples a color from the specified Plotly colorscale using a
    value between 0 and 1. It's useful for mapping numerical data (e.g., identity
    percentages) to a corresponding color.

    Parameters
    ----------
    value : float
        A normalized value between 0 and 1 indicating the position within the colorscale.
    colorscale_name : str, default="Greys"
        The name of the Plotly colorscale to use (e.g., "Greys", "Blues", "Viridis").

    Returns
    -------
    str
        A string representing the RGB color (e.g., "rgb(200, 200, 200)").
    """
    # Sample the color from the colorscale
    rgb_color = colors.sample_colorscale(colorscale_name, [value])[0]
    return rgb_color


def get_truncated_colorscale(
    colorscale_name: str = "Greys",
    vmin: float = 0,
    vmax: float = 0.75,
    n_samples: int = 256,
) -> list[tuple[float, str]]:
    """
    Generate a truncated Plotly colorscale between two normalized values.

    This function samples a subset of a Plotly colorscale between `vmin` and `vmax`,
    returning a list of (position, color) tuples. It is useful for focusing a colorscale
    on a specific range of values (e.g., homology identity percentages).

    For the "Greys" colorscale, the `vmax` is capped at 0.99 to prevent rendering issues
    that may cause the darkest values to appear incorrectly when `vmax=1`.

    Parameters
    ----------
    colorscale_name : str, default="Greys"
        The name of the Plotly colorscale to sample from.
    vmin : float, default=0
        The lower bound of the normalized range (between 0 and 1).
    vmax : float, default=0.75
        The upper bound of the normalized range (between 0 and 1).
    n_samples : int, default=256
        Number of samples to generate across the specified range.

    Returns
    -------
    list of tuple[float, str]
        A list of (normalized_position, color_string) tuples representing the truncated
        colorscale.
    """
    # IMPORTANT: This fix a problem with Greys set to 100% (i.e. vmax 1) that shows
    # shadows with off colors.
    if colorscale_name == "Greys" and vmax == 1:
        vmax = 0.99
    values = np.linspace(vmin, vmax, n_samples)
    truncated_colors = colors.sample_colorscale(colorscale_name, values)
    return truncated_colors


def sample_from_truncated_colorscale(
    truncated_colorscale: list[tuple[float, str]], homology_value: float
) -> str:
    """
    Sample a color from a truncated colorscale based on a normalized homology value.

    This function interpolates between the first and last colors in a truncated
    Plotly colorscale using the given homology value (between 0 and 1), and returns
    the corresponding RGB color.

    Parameters
    ----------
    truncated_colorscale : list of tuple[float, str]
        A colorscale defined by a list of (normalized_position, color_string) tuples.
        Typically produced by `get_truncated_colorscale`.
    homology_value : float
        A value between 0 and 1 indicating the relative position to sample within the
        colorscale.

    Returns
    -------
    str
        The interpolated RGB color string (e.g., "rgb(100, 150, 200)").
    """
    sampled_color = colors.find_intermediate_color(
        truncated_colorscale[0],  # The first color in the truncated colorscale
        truncated_colorscale[-1],  # The last color in the truncated colorscale
        homology_value,  # The input value between 0 and 1
        colortype="rgb",  # Return the color in RGB format
    )
    return sampled_color


def sample_colorscale_setting_lowest_and_highest_homologies(
    truncated_colorscale: list[tuple[float, str]],
    homology_value: float,
    lowest_homology: float,
    highest_homology: float,
) -> str:
    """
    Sample a color from a truncated colorscale based on actual homology value range.

    This function maps a homology identity value to a position within a truncated
    colorscale, taking into account the actual minimum and maximum identity values
    in the dataset. This ensures consistent color mapping even when the colorscale
    has been stretched or truncated.

    Parameters
    ----------
    truncated_colorscale : list of tuple[float, str]
        A colorscale defined by a list of (normalized_position, color_string) tuples.
    homology_value : float
        The homology identity value to map (e.g., 0.83 for 83% identity).
    lowest_homology : float
        The lowest identity value present in the dataset (used to normalize the range).
    highest_homology : float
        The highest identity value present in the dataset.

    Returns
    -------
    str
        The interpolated RGB color string corresponding to the input homology value.
    """
    delta_highest_to_lowest_homology = highest_homology - lowest_homology
    delta_highest_to_value_homology = highest_homology - homology_value
    if delta_highest_to_lowest_homology == 0:
        value = 1.0
    else:
        value = (
            1.0
            - (delta_highest_to_value_homology * 1.0) / delta_highest_to_lowest_homology
        )
    sampled_color = colors.find_intermediate_color(
        truncated_colorscale[0],  # The first color in the truncated colorscale
        truncated_colorscale[-1],  # The last color in the truncated colorscale
        value,  # The input value between 0 and 1
        colortype="rgb",  # Return the color in RGB format
    )
    return sampled_color


def plot_colorbar_legend(
    fig: Figure,
    colorscale: list[tuple[float, str]],
    min_value: float,
    max_value: float,
    set_colorscale_to_extreme_homologies: bool = False,
) -> Figure:
    """
    Add a horizontal colorbar legend to a Plotly figure to indicate homology identity
    range.

    This function creates a dummy scatter trace with a customized colorbar to serve
    as a legend for the homology identity colors. It adjusts the color scale display
    and tick labels based on whether the full identity range or a truncated scale
    is being used.

    Parameters
    ----------
    fig : plotly.graph_objects.Figure
        The Plotly figure to which the colorbar legend should be added.
    colorscale : list of tuple[float, str]
        A Plotly-compatible colorscale representing identity values.
    min_value : float
        The minimum identity value represented in the current color scale (normalized 0-1).
    max_value : float
        The maximum identity value represented in the current color scale (normalized 0-1).
    set_colorscale_to_extreme_homologies : bool, default=False
        If True, the color scale is stretched to the actual min and max identity values.
        If False, it uses a truncated range (e.g., based on a user-defined threshold).

    Returns
    -------
    fig : plotly.graph_objects.Figure
        The updated Plotly figure with a horizontal identity colorbar legend.
    """
    colorbar_len: float = 0.3
    title_position: str = "bottom"
    # Check if plot was set to set colorscale to extreme homologies
    if min_value != max_value and set_colorscale_to_extreme_homologies:
        updated_colorscale = colorscale
        tickvals = [min_value, max_value]
        ticktext = [f"{min_value*100:.2f}%", f"{max_value*100:.2f}%"]
    if min_value != max_value and not set_colorscale_to_extreme_homologies:
        updated_colorscale = get_truncated_colorscale(colorscale, min_value, max_value)
        tickvals = [min_value, max_value]
        ticktext = [f"{min_value*100:.2f}%", f"{max_value*100:.2f}%"]
    # If min and max values are the same add only one tick value.
    if min_value == max_value:
        updated_colorscale = get_truncated_colorscale(colorscale, min_value, max_value)
        tickvals = [max_value]
        ticktext = [f"{max_value * 100:.2f}%"]
        colorbar_len = 0.15
        title_position = "right"

    fig.add_trace(
        go.Scatter(
            y=[None],  # Dummy values
            x=[None],  # Dummy values
            customdata=[
                {
                    "min_identity": min_value,
                    "max_identity": max_value,
                }
            ],
            name="colorbar legend",
            mode="markers",
            marker=dict(
                colorscale=updated_colorscale,
                cmin=min_value,
                cmax=max_value,
                color=[min_value, max_value],
                colorbar=dict(
                    title=dict(
                        text="Identity", font=dict(size=18), side=title_position
                    ),
                    orientation="h",
                    x=0.75,
                    y=-0.02,
                    tickfont=dict(size=18),
                    tickvals=tickvals,
                    ticktext=ticktext,
                    len=colorbar_len,
                ),
            ),
            hoverinfo="none",
        )
    )
    return fig


def plot_polygon(
    fig: Figure,
    x_values: list,
    y_values: list,
    color: str = "blue",
    name: str = "",
    customdata: list = [],
    visible: bool = True,
) -> None:
    """
    Add a filled polygon to the Plotly figure to represent a gene or homology region.

    This function draws a closed shape by connecting the points defined by `x_values`
    and `y_values`, and fills it with the specified color. It is used to visually
    represent genes, coding regions, or homology blocks on the plot.

    Parameters
    ----------
    fig : plotly.graph_objects.Figure
        The Plotly figure to which the polygon should be added.
    x_values : list
        The x-coordinates of the polygon vertices.
    y_values : list
        The y-coordinates of the polygon vertices.
    color : str, default="blue"
        The color used for both the polygon outline and fill (RGB or hex format).
    name : str, optional
        A label used for hovering (displayed as `%{text}`).
    customdata : list, optional
        Extra data attached to each vertex, accessible in callbacks or hover events.
    visible : bool, default=True
        Whether the polygon is initially visible in the plot.

    Returns
    -------
    None
        The function modifies the input figure in-place and returns nothing.
    """
    fig.add_trace(
        go.Scatter(
            x=x_values,
            y=y_values,
            text=name,
            fill="toself",
            mode="lines",
            line=dict(color=color, width=1),
            fillcolor=color,
            name="",
            customdata=customdata,
            hoverlabel=dict(font_size=14),
            hovertemplate="%{text}<extra></extra>",
            visible=visible,
        )
    )


def plot_line(
    fig: Figure,
    x_values: list,
    y_values: list,
    name: str | None = None,
    customdata: list = [],
    color: str = "black",
) -> None:
    """
    Add a thick horizontal line to the Plotly figure to represent a DNA sequence.

    This function draws a straight line between the specified x and y coordinates.
    It is typically used to visualize the backbone of each DNA sequence in the alignment.

    Parameters
    ----------
    fig : plotly.graph_objects.Figure
        The Plotly figure to which the line should be added.
    x_values : list
        The x-coordinates defining the start and end of the line.
    y_values : list
        The y-coordinates defining the vertical position of the line (usually constant).
    name : str or None, optional
        Name label used for hover display and trace identification.
    customdata : list, optional
        Extra data attached to the line, accessible via hover or callbacks.
    color : str, default="black"
        The color of the line (e.g., "black", "#FF0000", or "rgb(0,0,0)").

    Returns
    -------
    None
        The function modifies the input figure in-place and returns nothing.
    """
    fig.add_trace(
        go.Scatter(
            x=x_values,
            y=y_values,
            mode="lines",
            name=name,
            customdata=customdata,
            line=dict(color=color, width=4),
            hoverlabel=dict(font_size=14),
            hovertemplate=f"{name}<extra></extra>",
        )
    )


def plot_dna_sequences(
    fig: Figure,
    gb_records: DataFrame,
    y_separation: int = 10,
) -> Figure:
    """
    Plot horizontal lines representing DNA sequences using metadata from a DataFrame.

    This function iterates over a DataFrame of GenBank metadata and draws one horizontal
    line per sequence, spaced vertically by `y_separation`. Each line represents a DNA
    sequence from the dataset and is labeled with the sequence name.

    Parameters
    ----------
    fig : plotly.graph_objects.Figure
        The Plotly figure to which the sequence lines will be added.
    gb_records : pandas.DataFrame
        A DataFrame containing metadata for each sequence, including:
        - 'sequence_start'
        - 'sequence_end'
        - 'record_name'
    y_separation : int, default=10
        The vertical distance between stacked sequence lines.

    Returns
    -------
    fig : plotly.graph_objects.Figure
        The updated Plotly figure with all DNA sequences plotted as horizontal lines.
    """
    y_distance = len(gb_records) * y_separation
    for _, row in gb_records.iterrows():
        x1 = row["sequence_start"]
        x2 = row["sequence_end"]
        record_name = row["record_name"]
        x_values = np.array([x1, x2])
        y_values = np.array([y_distance, y_distance])
        trace_name = f"Sequence: {record_name}"
        plot_line(
            fig,
            x_values,
            y_values,
            name=trace_name,
        )
        y_distance -= y_separation
    return fig


def plot_genes(
    fig: Figure,
    number_gb_records: int,
    longest_sequence: int,
    cds_records: DataFrame,
    name_from: str = "product",
    y_separation: int = 10,
) -> Figure:
    """
    Plot arrows representing genes using metadata from a CDS DataFrame.

    This function iterates over coding sequence (CDS) entries grouped by GenBank file.
    Each gene is rendered as a directional polygon (arrow) positioned on the y-axis
    according to the file order. Arrowhead size is scaled relative to the longest sequence
    to maintain proportionality.

    Parameters
    ----------
    fig : plotly.graph_objects.Figure
        The Plotly figure to which gene arrows will be added.
    number_gb_records : int
        The total number of GenBank files, used to calculate vertical placement.
    longest_sequence : int
        The length of the longest DNA sequence, used to scale arrowhead size.
    cds_records : pandas.DataFrame
        A DataFrame containing gene metadata. Expected columns include:
        - 'file_number'
        - 'start_plot'
        - 'end_plot'
        - 'color'
        - 'product'
        - 'gene'
    name_from : str, default="product"
        Determines the label shown for each gene: either "product" or "gene".
    y_separation : int, default=10
        Vertical spacing between rows of gene arrows for different GenBank files.

    Returns
    -------
    fig : plotly.graph_objects.Figure
        The updated Plotly figure with gene arrows plotted for each sequence.
    """
    # Position of the first DNA sequence in the y axis. Plotting starts at the top.
    y = number_gb_records * y_separation
    # Ratio head_height vs lenght of longest sequence
    ratio = 0.02
    head_height = longest_sequence * ratio
    # Iterate over gb_records dataframe to plot genes.
    for _, cds_group in cds_records.groupby(["file_number"]):
        for _, row in cds_group.iterrows():
            x1 = row["start_plot"]
            x2 = row["end_plot"]
            color = row["color"]
            arrow = Arrow(x1=x1, x2=x2, y=y, head_height=head_height)
            x_values, y_values = arrow.get_coordinates()
            # Get name and check if is None
            name = row["product"] if name_from == "product" else row["gene"]
            name = name if name is not None else "no name"
            plot_polygon(
                fig,
                x_values,
                y_values,
                color=color,
                name=name,
            )
        y -= y_separation
    return fig


def plot_homology_regions_with_dataframe(
    fig: Figure,
    alignments_df: DataFrame,
    regions_df: DataFrame,
    y_separation: int = 10,
    homology_padding: float = 1.1,
    colorscale: str = "Greys",
    straight_heights: bool = True,
    minimum_homology_length: int = 1,
    set_colorscale_to_extreme_homologies: bool = False,
    lowest_homology: None | float = None,
    highest_homology: None | float = None,
) -> Figure:
    """
    Plot homology regions as filled polygons using metadata from alignment dataframes.

    This function visualizes BLASTn alignment regions as colored polygons between
    query and subject sequences. The color represents the identity percentage, based
    on a colorscale that may be truncated or stretched to fit user-defined identity
    bounds.

    Parameters
    ----------
    fig : plotly.graph_objects.Figure
        The Plotly figure to which homology polygons will be added.
    alignments_df : pandas.DataFrame
        A DataFrame with metadata about BLAST alignments. Expected columsn include:
        - 'alignment_number', 'query_name', 'hit_name', 'query_len', and 'hit_len'.
    regions_df : pandas.DataFrame
        A DataFrame containing the coordinates of individual homology blocks to plot.
        Expected columns include:
        - 'alignment_number', 'query_from_plot', 'query_to_plot',
        - 'hit_from_plot', 'hit_to_plot', and 'homology'.
    y_separation : int, default=10
        Vertical distance between stacked alignment rows.
    homology_padding : float, default=1.1
        Padding added to the y-coordinates to offset polygon positioning for visual
        clarity.
    colorscale : str
        Name of the Plotly colorscale used to color-code identity values.
    straight_heights : bool, default=True
        If True, polygons are drawn with flat tops and bottoms.
        If False, curved edges are used for a smoother appearance.
    minimum_homology_length : int, default=1
        Minimum length (in bp) for a homology region to be displayed.
    set_colorscale_to_extreme_homologies : bool, default=False
        If True, colors are scaled using the dataset's actual min/max identity values.
    lowest_homology : float or None, optional
        The minimum identity value used for scaling (required if
        `set_colorscale_to_extreme_homologies` is True).
    highest_homology : float or None, optional
        The maximum identity value used for scaling (required if
        `set_colorscale_to_extreme_homologies` is True).

    Returns
    -------
    fig : plotly.graph_objects.Figure
        The updated Plotly figure with homology regions plotted as polygons.
    """
    # Get length of alignments and add 1
    alignments_len = len(alignments_df) + 1
    # Get the y distance to start plotting at the top of the graph
    y_distances = (alignments_len) * y_separation
    # Iterate over homologous regions for plotting
    for i, region_group in regions_df.groupby(["alignment_number"]):
        for _, row in region_group.iterrows():
            # Get region coordinates
            x1 = row["query_from_plot"]
            x2 = row["query_to_plot"]
            x3 = row["hit_to_plot"]
            x4 = row["hit_from_plot"]
            y1 = y_distances - homology_padding
            y2 = y_distances - homology_padding
            y3 = y_distances - y_separation + homology_padding
            y4 = y_distances - y_separation + homology_padding
            homology_length = x2 - x1
            # If homology length is less or equalts to the minimun required, ignore it
            if homology_length <= minimum_homology_length:
                visible = False
            else:
                visible = True
            # If user requested straight lines convert coordinates to np.array
            if straight_heights:
                xpoints = np.array([x1, x2, x3, x4, x1])
                ypoints = np.array([y1, y2, y3, y4, y1])
            # Otherwise, convert coordinates to bezier coordinates.
            else:
                xpoints, ypoints = RectangleCurveHeight(
                    x_coordinates=[x1, x2, x3, x4],
                    y_coordinates=[y1, y2, y3, y4],
                    proportions=[0, 0.2, 0.8, 1],
                ).coordinates_rectangle_height_bezier()
            # Get the identity to match with the correct color.
            homology = row["homology"]
            # Sample color depending on how the user set the colorscale
            if set_colorscale_to_extreme_homologies:
                color = sample_colorscale_setting_lowest_and_highest_homologies(
                    truncated_colorscale=colorscale,
                    homology_value=homology,
                    lowest_homology=lowest_homology,
                    highest_homology=highest_homology,
                )
            else:
                color = sample_from_truncated_colorscale(
                    truncated_colorscale=colorscale,
                    homology_value=homology,
                )
            customdata = ["identity", homology, homology_length]
            plot_polygon(
                fig,
                xpoints,
                ypoints,
                color=color,
                name=f"Identity: {homology*100:.2f}%",
                customdata=customdata,
                visible=visible,
            )
        y_distances -= y_separation
    return fig


def annotate_dna_sequences(
    fig: Figure,
    gb_records: DataFrame,
    longest_sequence: int,
    number_gb_records: int,
    annotate_with: str = "accession",
    y_separation: int = 10,
    padding: int = 10,
) -> Figure:
    """
    Add text annotations to DNA sequence lines using metadata from GenBank records.

    This function appends labels to the right of each sequence line, based on user
    preference. Labels can include the sequence's accession number, internal name, or
    original filename. The annotations are positioned with consistent vertical spacing.

    Parameters
    ----------
    fig : plotly.graph_objects.Figure
        The Plotly figure to which annotations will be added.
    gb_records : pandas.DataFrame
        A DataFrame containing metadata for each GenBank sequence. Expected columns
        include:
        - 'accession', 'record_name', and 'file_name'.
    longest_sequence : int
        Length of the longest sequence, used to offset the annotation on the x-axis.
    number_gb_records : int
        Total number of sequences to be annotated (used to calculate vertical placement).
    annotate_with : str, default="accession"
        Field to use for annotation. Must be one of:
        - "accession" (e.g., NCBI ID)
        - "name" (internal record name)
        - "fname" (original file name)
    y_separation : int, default=10
        Vertical spacing between each sequence line.
    padding : int, default=10
        Horizontal space between the end of the longest sequence and its annotation.

    Returns
    -------
    fig : plotly.graph_objects.Figure
        The updated Plotly figure with sequence annotations added.
    """
    y = y_separation * number_gb_records
    options = {"accession": "accession", "name": "record_name", "fname": "file_name"}
    option = options.get(annotate_with, "accession")
    for _, row in gb_records.iterrows():
        name = row[option]
        fig.add_annotation(
            x=longest_sequence + padding,
            xref="x",
            y=y,
            name=f"Sequence annotation: {name}",
            text=name,
            font=dict(size=18),
            showarrow=False,
            xanchor="left",
            yanchor="middle",
        )
        y -= y_separation
    return fig


def annotate_top_genes(
    fig: Figure,
    annotate_genes_with: str,
    number_gb_records: int,
    cds_records: DataFrame,
    y_separation: int = 10,
) -> Figure:
    """
    Add vertical text annotations for genes on the top sequence row.

    This function annotates genes from the first GenBank file (`file_number == 0`)
    using a specified metadata field (e.g., gene name or product). Annotations are
    positioned above the top row of gene arrows and are rotated vertically.

    Parameters
    ----------
    fig : plotly.graph_objects.Figure
        The Plotly figure to which annotations will be added.
    annotate_genes_with : str
        The metadata field to use for labeling each gene (e.g., "gene", "product").
    number_gb_records : int
        Total number of GenBank records, used to calculate the top y-axis position.
    cds_records : pandas.DataFrame
        A DataFrame containing CDS metadata. Must include the following columns:
        - 'file_number', 'start_plot', 'end_plot',
        - and the column specified by `annotate_genes_with` ('gene' and 'product').
    y_separation : int, default=10
        Vertical spacing between sequence rows; controls annotation height.

    Returns
    -------
    fig : plotly.graph_objects.Figure
        The updated Plotly figure with gene annotations added above the top sequence.
    """
    # Filter rows to get the ones that belong to the top sequence
    df_file_number_0 = cds_records.loc[cds_records["file_number"] == 0]
    y = y_separation * number_gb_records
    # Iterate over rows of file_number 0 to annotate genes
    for _, row in df_file_number_0.iterrows():
        x_start = row["start_plot"]
        x_end = row["end_plot"]
        x = (x_start + x_end) / 2
        name = row[annotate_genes_with]
        fig.add_annotation(
            x=x,
            y=y + 1.1,
            text=name,
            name=f"Gene annotation: {name}",
            showarrow=False,
            textangle=270,
            font=dict(size=16),
            xanchor="center",
            yanchor="bottom",
        )
    return fig


def annotate_bottom_genes(
    fig: Figure,
    annotate_genes_with: str,
    number_gb_records: int,
    cds_records: DataFrame,
    y_separation: int = 10,
) -> Figure:
    """
    Add vertical text annotations for genes on the bottom sequence row.

    This function annotates genes from the last GenBank file (`file_number ==
    number_gb_records - 1`) using a specified metadata field (e.g., gene name or product).
    Annotations are positioned below the bottom row of gene arrows and rotated vertically
    for readability.

    Parameters
    ----------
    fig : plotly.graph_objects.Figure
        The Plotly figure to which annotations will be added.
    annotate_genes_with : str
        The metadata field to use for labeling each gene (e.g., "gene" or "product").
    number_gb_records : int
        Total number of GenBank records, used to identify the bottom sequence.
    cds_records : pandas.DataFrame
        A DataFrame containing CDS metadata. Must include the following columns:
        - 'file_number', 'start_plot', 'end_plot'
        - and the column specified by `annotate_genes_with` ('gene' or 'product').
    y_separation : int, default=10
        Vertical spacing between sequence rows; controls annotation height.

    Returns
    -------
    fig : plotly.graph_objects.Figure
        The updated Plotly figure with gene annotations added below the bottom sequence.
    """
    # Filter rows to get the ones that belong to the bottom sequence
    df_file_number_0 = cds_records.loc[
        cds_records["file_number"] == number_gb_records - 1
    ]
    y = y_separation
    # Iterate over rows of file_number 0 to annotate genes
    for _, row in df_file_number_0.iterrows():
        x_start = row["start_plot"]
        x_end = row["end_plot"]
        x = (x_start + x_end) / 2
        name = row[annotate_genes_with]
        fig.add_annotation(
            x=x,
            y=y - 1.1,
            text=name,
            name=f"Gene annotation: {name}",
            showarrow=False,
            textangle=270,
            font=dict(size=16),
            xanchor="center",
            yanchor="top",
        )
    return fig


def annotate_genes(fig: Figure, plot_parameters: PlotParameters) -> Figure:
    """
    Annotate gene features on the top, bottom, or both sequence rows based on user
    preferences.

    This function uses the `PlotParameters` object to determine how and where gene
    annotations should be applied. It delegates to `annotate_top_genes` and/or
    `annotate_bottom_genes` depending on the value of `plot_parameters.annotate_genes`.

    Parameters
    ----------
    fig : plotly.graph_objects.Figure
        The Plotly figure to which gene annotations will be added.
    plot_parameters : PlotParameters
        An object containing configuration and metadata, including:
        - `annotate_genes`: str, one of "top", "bottom", or "top-bottom"
        - `annotate_genes_with`: str, column name to use for labels (e.g., "gene",
          "product")
        - `cds_df`: DataFrame with CDS metadata
        - `number_gb_records`: int, number of GenBank files
        - `y_separation`: int, spacing between sequence rows

    Returns
    -------
    fig : plotly.graph_objects.Figure
        The updated Plotly figure with gene annotations applied according to settings.
    """
    annotate_genes = plot_parameters.annotate_genes
    if annotate_genes == "top" or annotate_genes == "top-bottom":
        fig = annotate_top_genes(
            fig=fig,
            annotate_genes_with=plot_parameters.annotate_genes_with,
            number_gb_records=plot_parameters.number_gb_records,
            cds_records=plot_parameters.cds_df,
            y_separation=plot_parameters.y_separation,
        )
    if annotate_genes == "bottom" or annotate_genes == "top-bottom":
        fig = annotate_bottom_genes(
            fig=fig,
            annotate_genes_with=plot_parameters.annotate_genes_with,
            number_gb_records=plot_parameters.number_gb_records,
            cds_records=plot_parameters.cds_df,
            y_separation=plot_parameters.y_separation,
        )
    return fig


def remove_traces_by_name(figure: dict, name: str) -> dict:
    """
    Remove all traces from a Plotly figure dictionary whose 'name' contains a given
    substring.

    This function is useful for removing dynamically generated traces (e.g., colorbar
    legends, annotations, or temporary highlights) based on partial name matching.

    Parameters
    ----------
    figure : dict
        A dictionary representing a Plotly figure (as returned by dcc.Graph).
    name : str
        A substring to match within the 'name' field of each trace. Any trace whose
        name contains this string will be removed.

    Returns
    -------
    dict
        The updated figure dictionary with the matching traces removed.
    """
    data = []
    for trace in figure["data"]:
        if ("name" not in trace) or (name not in trace["name"]):
            data.append(trace)
    figure["data"] = data
    return figure


def hide_homology(figure: Figure, min_homology: int) -> Figure:
    """
    Toggle visibility of homology traces based on their length.

    This function inspects each trace in the Plotly figure and hides those representing
    homology regions whose length is less than or equal to `min_homology`. Length
    information is expected in the trace's `customdata`.

    The expected structure of `customdata` is:
        [ "identity", identity_value: float, homology_length: int ]

    Parameters
    ----------
    figure : plotly.graph_objects.Figure
        A Plotly figure object containing homology region traces with `customdata`.
    min_homology : int
        Minimum homology length required for a region to remain visible.

    Returns
    -------
    figure : plotly.graph_objects.Figure
        The updated figure with trace visibility toggled based on homology length.
    """
    figure_type = type(figure)
    print(f"the type of figure in hide_homology is: {figure_type}")
    for trace in figure["data"]:
        if trace["customdata"] is not None:
            if "identity" in trace["customdata"]:
                homology_length = int(trace["customdata"][2])
                if homology_length <= min_homology:
                    trace.visible = False
                if homology_length > min_homology:
                    trace.visible = True
    return figure


def remove_annotations_by_name(figure: Figure, name: str) -> Figure:
    """
    Remove annotations from a Plotly figure whose 'name' contains the given substring.

    This function filters `figure.layout['annotations']` and removes any annotations
    where the 'name' field contains the specified string. It is useful for dynamically
    updating or clearing labeled elements such as gene or sequence annotations.

    Parameters
    ----------
    figure : plotly.graph_objects.Figure
        A Plotly figure object that may contain annotations.
    name : str
        A substring to search for in the 'name' field of each annotation.
        Annotations matching this string will be removed.

    Returns
    -------
    figure : plotly.graph_objects.Figure
        The updated Plotly figure with matching annotations removed.
    """
    annotations = []
    for annotation in figure.layout["annotations"]:
        if name in annotation["name"]:
            continue
        annotations.append(annotation)
    annotations = tuple(annotations)
    figure.layout["annotations"] = annotations
    return figure


def make_selection_effect(figure: Figure, curve_number: int) -> Figure:
    """
    Visually highlight a selected trace by changing its line color and thickness.

    This function adjusts the line color and width of a specific trace (based on its
    `curve_number`) to create a selection effect. It uses the trace's current color to
    determine whether to use a dark or light outline for optimal contrast.

    Parameters
    ----------
    figure : plotly.graph_objects.Figure
        The Plotly figure containing the trace to modify.
    curve_number : int
        The index of the trace to update (corresponds to 'curveNumber' from Dash
        clickData).

    Returns
    -------
    figure : plotly.graph_objects.Figure
        The updated figure with the specified trace visually emphasized.
    """
    # Make effect of selection by changing the color of the line
    default_black = "rgb(30, 30, 30)"
    default_light = "rgb(230, 230, 230)"
    color_curve = figure["data"][curve_number]["line"]["color"]
    # if color is hex change to rgb list
    if "#" in color_curve:
        color_curve = mcolors.to_rgb(color_curve)
    # otherwise convert rgb to list
    else:
        color_curve = color_curve.replace("rgb(", "").replace(")", "").split(",")
        # normalize rgb
        color_curve = [float(n) / 255 for n in color_curve]
    # define lightness
    lightness = (
        0.2126 * color_curve[0] + 0.7152 * color_curve[1] + 0.0722 * color_curve[2]
    )
    # if color is too dark use make the color line light
    if lightness < 0.3:
        color_line = default_light
    # else black
    else:
        color_line = default_black
    # Update color line
    figure["data"][curve_number]["line"]["color"] = color_line
    figure["data"][curve_number]["line"]["width"] = 6
    return figure


def change_homology_color(
    figure: dict,
    colorscale_name: str,
    vmin_truncate: float,
    vmax_truncate: float,
    set_colorscale_to_extreme_homologies: bool = False,
    lowest_homology: None | float = None,
    highest_homology: None | float = None,
) -> dict:
    """
    Update the color of homology region traces based on identity values and a selected colorscale.

    This function iterates over traces in a Plotly figure and updates the `fillcolor` and
    `line.color` properties based on the homology identity value stored in `customdata`.
    Color mapping is done using either a truncated or full-range colorscale.

    Expected structure of `customdata` for homology regions:
        [ "identity", identity_value: float, homology_length: int ]

    Parameters
    ----------
    figure : dict
        A dictionary representation of a Plotly figure (typically from Dash state).
    colorscale_name : str
        Name of the Plotly colorscale to use (e.g., "Greys", "Viridis").
    vmin_truncate : float
        Lower bound (0–1) of the normalized identity range for color mapping.
    vmax_truncate : float
        Upper bound (0–1) of the normalized identity range for color mapping.
    set_colorscale_to_extreme_homologies : bool, default=False
        If True, stretch the color scale based on dataset-wide min/max identity values.
    lowest_homology : float or None
        The minimum identity value in the dataset, required if using extreme homology
        scaling.
    highest_homology : float or None
        The maximum identity value in the dataset, required if using extreme homology
        scaling.

    Returns
    -------
    figure : dict
        The updated Plotly figure dictionary with modified homology trace colors.
    """
    # Get new colorscale
    colorscale = get_truncated_colorscale(
        colorscale_name=colorscale_name,
        vmin=vmin_truncate,
        vmax=vmax_truncate,
    )
    for trace in figure["data"]:
        if "customdata" not in trace:
            continue
        if trace["customdata"] is None:
            continue
        if "identity" in trace["customdata"]:
            # Get identity information from customdata
            identity = trace["customdata"][1]
            identity = float(identity)
            # Sample colorscale with identity value.
            if set_colorscale_to_extreme_homologies:
                color = sample_colorscale_setting_lowest_and_highest_homologies(
                    truncated_colorscale=colorscale,
                    homology_value=identity,
                    lowest_homology=lowest_homology,
                    highest_homology=highest_homology,
                )
            else:
                color = sample_from_truncated_colorscale(
                    truncated_colorscale=colorscale, homology_value=identity
                )
            trace["fillcolor"] = color
            trace["line"]["color"] = color

    return figure


def round_up_to_nearest_significant_digit(number: float) -> int:
    """
    Round up a float to the nearest significant digit.

    This function rounds a given number up to the next multiple of the most significant
    power of ten. It is useful for generating clean scale values or axis limits in plots.
    This function is used in HomologyViz to show the scale bar value.

    Examples
    --------
    >>> round_up_to_nearest_significant_digit(142)
    200

    Parameters
    ----------
    number : float
        The number to round up.

    Returns
    -------
    int
        The input number rounded up to the nearest significant digit.
    """
    # Determine the nearest power of ten (e.g., 1000, 100, 10, etc.)
    power_of_ten = 10 ** math.floor(math.log10(number))
    # Round up to the next multiple of that power
    return math.ceil(number / power_of_ten) * power_of_ten


def plot_scale(
    figure: Figure, length_longest_sequence: int, add_scale: bool = True
) -> Figure:
    """
    Add a scale bar to the plot representing distance in base pairs (bp).

    The scale length is calculated as one-fifth of the longest sequence, rounded up to the
    nearest significant digit for visual clarity. A line and label are drawn at the bottom
    of the figure to indicate scale. The scale can be hidden by setting `add_scale=False`.

    Parameters
    ----------
    figure : plotly.graph_objects.Figure
        The Plotly figure to which the scale bar will be added.
    length_longest_sequence : int
        Length of the longest DNA sequence in the plot (used to size the scale).
    add_scale : bool, default=True
        Whether to display the scale. If False, a fully transparent line and label are
        added.

    Returns
    -------
    figure : plotly.graph_objects.Figure
        The updated figure with the scale line and annotation.
    """
    # TODO: Dynamically position scale bar annotation based on y_separation to avoid
    # overlapping with bottom traces.

    scale_length: int = round_up_to_nearest_significant_digit(
        length_longest_sequence / 5
    )
    color: str = "rgba(0, 0, 0, 1)" if add_scale else "rgba(0,0,0,0)"
    # add line representing scale
    figure.add_trace(
        go.Scatter(
            x=[0, scale_length],
            y=[0, 0],
            mode="lines",
            name="Scale trace",
            line=dict(color=color, width=4),
            showlegend=False,
            hoverinfo="skip",
        )
    )
    # add annotation to scale
    figure.add_annotation(
        x=scale_length / 2,
        y=-1,
        text=f"{scale_length:,.0f} bp",
        name="Scale annotation",
        showarrow=False,
        yshift=-10,
        font=dict(size=18, color=color),
        hovertext=None,
        hoverlabel=None,
    )
    return figure


def toggle_scale_bar(figure: Figure, show: bool) -> Figure:
    """
    Toggle the visibility of the scale bar by adjusting its alpha channel.

    This function searches the figure for the scale bar trace and annotation using
    their names ("Scale trace" and "Scale annotation"). It modifies the color's
    alpha value to either show or fully hide them without removing the elements.

    Parameters
    ----------
    figure : plotly.graph_objects.Figure
        The Plotly figure containing the scale bar trace and annotation.
    show : bool
        If True, make the scale bar visible. If False, hide it using transparent colors.

    Returns
    -------
    figure : plotly.graph_objects.Figure
        The updated figure with the scale bar toggled on or off.
    """
    color: str = "rgba(0, 0, 0, 1)" if show else "rgba(0,0,0,0)"
    for trace in figure["data"]:
        if ("name" in trace) and ("Scale trace" in trace["name"]):
            trace["line"]["color"] = color
    for annotation in figure.layout["annotations"]:
        if "Scale annotation" in annotation["name"]:
            annotation["font"]["color"] = color
    return figure


def make_alignments(
    input_files: list[Path], output_folder: Path
) -> tuple[DataFrame, DataFrame, DataFrame, DataFrame]:
    """
    Perform BLASTn alignments between input sequences and return metadata for plotting.

    This function processes GenBank files to extract DNA sequences, runs BLASTn locally
    to identify alignments, and parses the results into structured Pandas DataFrames
    for visualization.

    Parameters
    ----------
    input_files : list of pathlib.Path
        List of file paths pointing to input GenBank files.
    output_folder : pathlib.Path
        Directory where temporary BLAST-related files will be written.

    Returns
    -------
    gb_df : pandas.DataFrame
        DataFrame with GenBank record metadata, including:
        - file path, file name, sequence length, file number, and accession number.
        Used to relate to CDS data (`cds_df`) by file number or accession.

    cds_df : pandas.DataFrame
        DataFrame with coding sequence (CDS) metadata for each GenBank file, including:
        - gene name, product name, strand, start/end coordinates (both raw and for
        plotting), file number, and accession number.

    alignments_df : pandas.DataFrame
        DataFrame summarizing BLASTn results per alignment, including:
        - alignment number, query name, hit name, query length, and hit length.

    regions_df : pandas.DataFrame
        DataFrame describing individual regions of sequence similarity (homology)
        between aligned sequences. Includes:
        - alignment number (to relate to `alignments_df`) and region-level coordinates.
        See `parse_blast_record` in the `gb_files_manipulation` module for details.
    """
    # Create fasta files for BLASTing using the gb files
    faa_files = genbank.make_fasta_files(input_files, output_folder)
    # Run blastn locally to make alignments.
    blast_xml_results = genbank.run_blastn(faa_files, output_folder)
    # Make alignments and regions dataframes from blast results
    alignments_df, regions_df = genbank.get_blast_metadata(blast_xml_results)
    # Make GenBank records and coding sequences dataframes
    gb_df, cds_df = genbank.genbank_files_metadata_to_dataframes(input_files)
    # Delete the documents used for genereting BLASTn results.
    misc.delete_files(faa_files)
    misc.delete_files(blast_xml_results)
    return gb_df, cds_df, alignments_df, regions_df


def make_figure(plot_parameters: PlotParameters) -> Figure:
    """
    Create a multiple sequence alignment plot using BLASTn and GenBank metadata.

    This function builds a complete Plotly figure for visualizing DNA sequences, gene
    annotations, and homology regions between multiple GenBank files. It uses the
    configuration and metadata provided in a fully populated `PlotParameters` object.

    Components included in the plot (depending on user settings):
    - DNA sequences as horizontal lines
    - Genes as directional polygons (arrows)
    - Homology regions as color-coded polygons
    - Annotations for sequences and genes
    - Optional scale bar (in base pairs)
    - Optional color scale legend for identity values

    Parameters
    ----------
    plot_parameters : PlotParameters
        Object containing all configuration settings and metadata needed to construct
        the plot. Must include parsed GenBank records, BLAST alignment results, color
        settings, visibility options, and layout preferences.

    Returns
    -------
    fig : plotly.graph_objects.Figure
        A fully assembled Plotly figure ready for display in a Dash app.
    """
    # TODO: Validate that all required attributes in plot_parameters are set before
    # plotting.

    # Before plotting, check if `alignments_position` option is not selected to the
    # `left` to adjust the coordinates of the sequences and genes accordingly.
    if (
        plot_parameters.draw_from_button == "plot-button"
        and plot_parameters.alignments_position != "left"
    ):
        genbank.adjust_positions_sequences_and_alignments_df_for_plotting(
            gb_records=plot_parameters.gb_df,
            cds=plot_parameters.cds_df,
            alignments=plot_parameters.alignments_df,
            regions=plot_parameters.alignments_regions_df,
            size_longest_sequence=plot_parameters.longest_sequence,
            position=plot_parameters.alignments_position,
        )

    # Create a blank figure
    fig = go.Figure()

    # Customize layout
    fig.update_layout(
        showlegend=False,
        xaxis=dict(
            showline=False,
            showgrid=False,
            showticklabels=False,
            zeroline=False,
        ),
        yaxis=dict(
            showline=False, showticklabels=False, showgrid=False, zeroline=False
        ),
        plot_bgcolor="white",
        hovermode="closest",
    )

    # Get lowest and hightest homologies.
    lowest_identity, highest_identity = (
        genbank.find_lowest_and_highest_homology_dataframe(
            plot_parameters.alignments_regions_df
        )
    )
    # Add lowest and highest identities to plot_parameters
    plot_parameters.lowest_identity = lowest_identity
    plot_parameters.highest_identity = highest_identity
    # Check if user set the colorscale to extreme homologies
    set_colorscale_to_extreme_homologies = (
        plot_parameters.set_colorscale_to_extreme_homologies
    )

    # Select colormap and its range to plot the homology regions
    colorscale = get_truncated_colorscale(
        colorscale_name=plot_parameters.identity_color,
        vmin=plot_parameters.colorscale_vmin,
        vmax=plot_parameters.colorscale_vmax,
    )

    # Plot the DNA sequences
    fig = plot_dna_sequences(
        fig=fig,
        gb_records=plot_parameters.gb_df,
        y_separation=plot_parameters.y_separation,
    )

    # Plot the homology regions
    straight_homology_regions = (
        True if plot_parameters.straight_homology_regions == "straight" else False
    )

    fig = plot_homology_regions_with_dataframe(
        fig=fig,
        alignments_df=plot_parameters.alignments_df,
        regions_df=plot_parameters.alignments_regions_df,
        y_separation=plot_parameters.y_separation,
        colorscale=colorscale,
        straight_heights=straight_homology_regions,
        minimum_homology_length=plot_parameters.minimum_homology_length,
        set_colorscale_to_extreme_homologies=set_colorscale_to_extreme_homologies,
        lowest_homology=lowest_identity,
        highest_homology=highest_identity,
    )

    # Plot genes
    fig = plot_genes(
        fig=fig,
        number_gb_records=plot_parameters.number_gb_records,
        longest_sequence=plot_parameters.longest_sequence,
        cds_records=plot_parameters.cds_df,
        y_separation=plot_parameters.y_separation,
    )
    # Annotate genes
    if plot_parameters.annotate_genes != "no":
        fig = annotate_genes(fig, plot_parameters)

    # Annotate DNA sequences
    if plot_parameters.annotate_sequences != "no":
        fig = annotate_dna_sequences(
            fig=fig,
            gb_records=plot_parameters.gb_df,
            longest_sequence=plot_parameters.longest_sequence,
            number_gb_records=plot_parameters.number_gb_records,
            annotate_with=plot_parameters.annotate_sequences,
            y_separation=plot_parameters.y_separation,
        )

    # Plot DNA scale
    if plot_parameters.add_scale_bar == "yes":
        fig = plot_scale(
            fig, plot_parameters.longest_sequence, plot_parameters.add_scale_bar
        )

    # Plot colorscale legend
    fig = plot_colorbar_legend(
        fig,
        colorscale,
        lowest_identity,
        highest_identity,
        set_colorscale_to_extreme_homologies=set_colorscale_to_extreme_homologies,
    )

    return fig
