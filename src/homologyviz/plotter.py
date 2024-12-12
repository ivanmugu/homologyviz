"""Make a graphical representation of BLASTn alignments.

Homology Visualization (HomologyViz) uses GenBank files (.gb) to align the sequences and
plot the genes. HomologyViz uses the information from the `CDS` features section to plot
the genes. To customize the colors for plotting genes, you can add a `Color` tag in the
`CDS` features with a color in hexadecimal. For example, add the tag `/Color="#00ff00"`
to show a green gene.

License
-------
This file is part of HomologyViz
BSD 3-Clause License
Copyright (c) 2024, Ivan Munoz Gutierrez
"""

import math
from pathlib import Path
from pandas import DataFrame

import plotly.graph_objects as go
from plotly.graph_objects import Figure, Heatmap
import plotly.colors as colors
import numpy as np
import matplotlib.colors as mcolors

from homologyviz.arrow import Arrow
from homologyviz.rectangle_bezier import RectangleCurveHeight
from homologyviz import gb_files_manipulation as genbank
from homologyviz.gb_files_manipulation import PlotParameters
from homologyviz import miscellaneous as misc


# TODO: I have to check how to store the data of bezier lines in a dataframe.
# ====================================================================================== #
#                             FUNCTIONS FOR PLOTTING                                     #
# ====================================================================================== #
def create_color_line(colors) -> Figure:
    """Create a continues color line based on selected color scale.

    Function to show a colormap in the `Edit` tab to set the colorscale range.
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
    """Get the RGB value from a Plotly colorscale given a value between 0 and 1."""
    # Sample the color from the colorscale
    rgb_color = colors.sample_colorscale(colorscale_name, [value])[0]
    return rgb_color


def get_truncated_colorscale(
    colorscale_name: str = "Greys",
    vmin: float = 0,
    vmax: float = 0.75,
    n_samples: int = 256,
) -> list[tuple[float, str]]:
    """Get truncated colorscale between vmin and vmax"""
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
    """Sample a color from a truncated colorscale given a value between 0 and 1."""
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
    """Sample colorscale by setting lowest and highest homologies"""
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
    """Plot colorbar legend."""
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
    customdata: list = [{}],
) -> None:
    """Plot polygon representing genes or homology regions"""
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
        )
    )


def plot_line(
    fig: Figure,
    x_values: list,
    y_values: list,
    name: str | None = None,
    customdata: list = [{}],
    color: str = "black",
) -> None:
    """Plot line representing DNA sequences"""
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


def plot_dna_sequences_with_dataframe(
    fig: Figure,
    gb_records: DataFrame,
    longest_sequence: int,
    y_separation: int = 10,
) -> Figure:
    """Plot DNA sequences"""
    y_distance = len(gb_records) * y_separation
    for _, row in gb_records.iterrows():
        x1 = row["sequence_start"]
        x2 = row["sequence_end"]
        accession = row["accession"]
        record_name = row["record_name"]
        file_name = row["file_name"]
        x_values = np.array([x1, x2])
        y_values = np.array([y_distance, y_distance])
        trace_name = f"Sequence: {record_name}"
        plot_line(
            fig,
            x_values,
            y_values,
            name=trace_name,
            customdata=[
                {
                    "type": "sequence_info",
                    "accession": accession,
                    "name": record_name,
                    "fname": file_name,
                    "longest_sequence": longest_sequence,
                },
            ],
        )
        y_distance -= y_separation
    return fig


def plot_genes_with_dataframe(
    fig: Figure,
    gb_records: DataFrame,
    cds_records: DataFrame,
    name_from: str = "product",
    y_separation: int = 10,
) -> Figure:
    """Plot arrows representing genes using metadata stored in a pandas DataFrame.

    The customdata has information for annotating the plot dinamically in Dash.
    """
    # Number of gb_records (rows)
    number_gb_records = len(gb_records)
    # Position of the first DNA sequence in the y axis. Plotting starts at the top.
    y = number_gb_records * y_separation
    # Find longest sequence
    longest_sequence = genbank.get_longest_sequence_dataframe(gb_records)
    # Ratio head_height vs lenght of longest sequence
    ratio = 0.02
    head_height = longest_sequence * ratio
    # cds_records_groups = cds_records.groupby(["file_number"])
    # Iterate over gb_records dataframe to plot genes.
    for i, cds_group in cds_records.groupby(["file_number"]):
        for _, row in cds_group.iterrows():
            file_numer = row["file_number"]
            cds_number = row["cds_number"]
            plot_id = f"{file_numer}.{cds_number}"
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
                customdata=[
                    {
                        "type": "gene_info",
                        "gb_record_number": i[0],
                        "number_gb_records": number_gb_records,
                        "x_start": row["start"],
                        "x_end": row["end"],
                        "x_start_plot": row["start_plot"],
                        "x_end_plot": row["end_plot"],
                        "y": y,
                        "plot_id": plot_id,
                    }
                ],
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
    minimum_homology_length: int = 0,
    set_colorscale_to_extreme_homologies: bool = False,
    lowest_homology: None | float = None,
    highest_homology: None | float = None,
) -> Figure:
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
                continue
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
            # Plot the homology region
            plot_polygon(
                fig,
                xpoints,
                ypoints,
                color=color,
                name=f"Identity: {homology*100:.2f}%",
                customdata=[
                    {
                        "type": "homology_info",
                        "identity": homology,
                        "homology_length": homology_length,
                    },
                ],
            )
        y_distances -= y_separation
    return fig


def annotate_dna_sequences_using_trace_customdata(
    figure: Figure, annotate_with: str, padding: int = 10
) -> Figure:
    """Annotate dna sequences using trace customdata

    The trace customdata is created when plotting using the plot_dna_sequences function.
    annotate_with options are "accession", "name", "fname".
    """
    for trace in figure["data"]:
        if "customdata" not in trace:
            continue
        if trace["customdata"] is None:
            continue
        if "sequence_info" in trace["customdata"][0].get("type", ""):
            text = trace["customdata"][0][annotate_with]
            figure.add_annotation(
                x=trace["customdata"][0]["longest_sequence"] + padding,
                xref="x",
                y=trace["y"][0],
                name=f"Sequence annotation: {text}",
                text=text,
                font=dict(size=18),
                showarrow=False,
                xanchor="left",
                yanchor="middle",
            )
    return figure


def annotate_genes_top_using_trace_customdata(figure: Figure) -> Figure:
    """Annotate top genes using trace customdata

    The trace customdata is created when plotting using the plot_genes function.
    """
    for trace in figure["data"]:
        if "customdata" not in trace:
            continue
        if trace["customdata"] is None:
            continue
        if "gene_info" in trace["customdata"][0].get("type", ""):
            gb_record = trace["customdata"][0]["gb_record_number"]
            # If gb_record is not the top (i. e. zero), continue.
            if gb_record != 0:
                continue
            x_start = float(trace["customdata"][0]["x_start_plot"])
            x_end = float(trace["customdata"][0]["x_end_plot"])
            y = float(trace["customdata"][0]["y"])
            x_position = (x_start + x_end) / 2
            y_position = y
            name = trace["text"]
            figure.add_annotation(
                x=x_position,
                y=y_position + 1.1,
                text=name,
                name=f"Gene annotation: {name}",
                showarrow=False,
                textangle=270,
                font=dict(size=16),
                xanchor="center",
                yanchor="bottom",
            )
    return figure


def annotate_genes_bottom_using_trace_customdata(figure: Figure) -> Figure:
    """Annotate bottom genes using trace customdata

    The trace customdata is created when plotting using the plot_genes function.
    """
    for trace in figure["data"]:
        if "customdata" not in trace:
            continue
        if trace["customdata"] is None:
            continue
        if "gene_info" in trace["customdata"][0].get("type", ""):
            gb_record = trace["customdata"][0]["gb_record_number"]
            number_gb_records = trace["customdata"][0]["number_gb_records"]
            # if gb_record is not the last one, continue
            if gb_record != (number_gb_records - 1):
                continue
            # print(trace)
            x_start = float(trace["customdata"][0]["x_start_plot"])
            x_end = float(trace["customdata"][0]["x_end_plot"])
            y = float(trace["customdata"][0]["y"])
            x_position = (x_start + x_end) / 2
            y_position = y
            name = trace["text"]
            figure.add_annotation(
                x=x_position,
                y=y_position - 1.1,
                text=name,
                name=f"Gene annotation: {name}",
                showarrow=False,
                textangle=270,
                font=dict(size=16),
                xanchor="center",
                yanchor="top",
            )
    return figure


def remove_traces_by_name(figure: Figure, name: str) -> dict:
    """Remove traces that have the any instance of the str name in trace['name']"""
    data = []
    for trace in figure["data"]:
        if ("name" not in trace) or (name not in trace["name"]):
            data.append(trace)
    figure["data"] = data
    return figure


def remove_annotations_by_name(figure: Figure, name: str) -> dict:
    """
    Remove annotations that have any instance of the str name in
    figure.layout['annotations']
    """
    annotations = []
    for annotation in figure.layout["annotations"]:
        if name in annotation["name"]:
            continue
        annotations.append(annotation)
    annotations = tuple(annotations)
    figure.layout["annotations"] = annotations
    return figure


def make_selection_effect(figure: Figure, curve_number: int):
    """Change border line of trace to make a selection efect."""
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


def change_homoloy_color_traces(
    figure: Figure,
    colorscale_name: str,
    vmin_truncate: float,
    vmax_truncate: float,
    set_colorscale_to_extreme_homologies: bool = False,
    lowest_homology: None | float = None,
    highest_homology: None | float = None,
) -> Figure:
    """Change the Figure's homology color traces

    This function works using the traces' customdata.
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
        if "homology_info" in trace["customdata"][0].get("type", ""):
            print(trace["customdata"][0].get("type", ""))
            # Get identity information from customdata
            identity = trace["customdata"][0]["identity"]
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
    # Determine the nearest power of ten (e.g., 1000, 100, 10, etc.)
    power_of_ten = 10 ** math.floor(math.log10(number))
    # Round up to the next multiple of that power
    return math.ceil(number / power_of_ten) * power_of_ten


def plot_scale(
    figure: Figure, length_longest_sequence: int, add_scale: bool = True
) -> Figure:
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
    """BLAST nucleotide sequences to make alignments."""
    # Create fasta files for BLASTing using the gb files
    faa_files = genbank.make_fasta_files(input_files, output_folder)
    # Run blastn locally to make alignments.
    blast_xml_results = genbank.run_blastn(faa_files, output_folder)
    # Make alignments and regions dataframes from blast results
    alignments_df, regions_df = genbank.blast_alignments_to_dataframe(blast_xml_results)
    # Make GenBank records and coding sequences dataframes
    gb_df, cds_df = genbank.genbank_files_metadata_to_dataframes(input_files)
    # Delete the documents used for genereting BLASTn results.
    misc.delete_files(faa_files)
    misc.delete_files(blast_xml_results)
    return gb_df, cds_df, alignments_df, regions_df


# ====================================================================================== #
#                                  MAIN FUNCTION GUI                                     #
# ====================================================================================== #
def make_figure(plot_parameters: PlotParameters) -> Figure:
    """Make a multiple sequence alignment plot

    Return
    ------
    Figure : plotly.graph_objects.Figure
    """
    # Before plotting, we must check the alignments position option to adjust the
    # coordinates of the sequences and genes.
    if (
        plot_parameters.draw_from_button == "draw-button"
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
    else:
        genbank.adjust_positions_sequences_and_alignments_df_for_plotting(
            gb_records=plot_parameters.gb_df,
            cds=plot_parameters.cds_df,
            alignments=plot_parameters.alignments_df,
            regions=plot_parameters.alignments_regions_df,
            size_longest_sequence=plot_parameters.longest_sequence,
            position=plot_parameters.alignments_position,
        )

    # = PLOT ALIGNMENTS =================================================================
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
    fig = plot_dna_sequences_with_dataframe(
        fig,
        plot_parameters.gb_df,
        plot_parameters.longest_sequence,
    )
    # Plot the homology regions
    straight_homology_regions = (
        True if plot_parameters.straight_homology_regions == "straight" else False
    )

    fig = plot_homology_regions_with_dataframe(
        fig,
        alignments_df=plot_parameters.alignments_df,
        regions_df=plot_parameters.alignments_regions_df,
        colorscale=colorscale,
        straight_heights=straight_homology_regions,
        minimum_homology_length=plot_parameters.minimum_homology_length,
        set_colorscale_to_extreme_homologies=set_colorscale_to_extreme_homologies,
        lowest_homology=lowest_identity,
        highest_homology=highest_identity,
    )

    fig = plot_genes_with_dataframe(
        fig=fig,
        gb_records=plot_parameters.gb_df,
        cds_records=plot_parameters.cds_df,
    )
    # Annotate genes
    if plot_parameters.annotate_genes == "top":
        fig = annotate_genes_top_using_trace_customdata(fig)
    if plot_parameters.annotate_genes == "bottom":
        fig = annotate_genes_bottom_using_trace_customdata(fig)
    if plot_parameters.annotate_genes == "top-bottom":
        fig = annotate_genes_top_using_trace_customdata(fig)
        fig = annotate_genes_bottom_using_trace_customdata(fig)
    # Annotate DNA sequences
    if (sequence_name := plot_parameters.annotate_sequences) != "no":
        # fig = annotate_dna_sequences(fig, gb_records, longest_sequence, sequence_name)
        fig = annotate_dna_sequences_using_trace_customdata(fig, sequence_name)

    # Plot DNA scale
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
