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

import plotly.graph_objects as go
from plotly.graph_objects import Figure, Heatmap
import plotly.colors as colors
import numpy as np
import matplotlib.colors as mcolors

from homologyviz.arrow import Arrow
from homologyviz.rectangle_bezier import RectangleCurveHeight
from homologyviz import gb_files_manipulation as genbank
from homologyviz.gb_files_manipulation import GenBankRecord, BlastnAlignment
from homologyviz import miscellaneous as misc


# TODO: when adjusting alignments to the center or right, don't modify the cds data.
# Instead, create a dictionary that stores the positions of the genes and sequence.


# = FUNCTIONS FOR ALIGNMENTS ============================================================
def get_longest_sequence(gb_records: list[GenBankRecord]) -> int:
    """Find the longest sequence in the gb records."""
    longest = 0
    for record in gb_records:
        if record.length > longest:
            longest = record.length
    return longest


def find_lowest_and_highest_homology(alignments: list[BlastnAlignment]) -> tuple:
    """Find the lowest and highest homology in a list of BlastnAlignmet objects."""
    lowest = 1
    highest = 0
    for alignment in alignments:
        for region in alignment.regions:
            if region.homology < lowest:
                lowest = region.homology
            if region.homology > highest:
                highest = region.homology
    return lowest, highest


def adjust_positions_sequences_right(
    gb_records: list[GenBankRecord], size_longest_sequence: int
) -> None:
    """Adjust the position of the sequences to the right including data in CDSs."""
    for record in gb_records:
        delta = size_longest_sequence - record.length
        record.sequence_start = record.sequence_start + delta
        record.sequence_end = record.sequence_end + delta
        for sequence in record.cds:
            sequence.start = sequence.start + delta
            sequence.end = sequence.end + delta


def adjust_positions_alignments_right(
    alignments: list[BlastnAlignment], size_longest_sequence: int
) -> None:
    """Adjust position of the alignments to the right."""
    for alignment in alignments:
        delta_query = size_longest_sequence - alignment.query_len
        delta_hit = size_longest_sequence - alignment.hit_len
        for region in alignment.regions:
            region.query_from = region.query_from + delta_query
            region.query_to = region.query_to + delta_query
            region.hit_from = region.hit_from + delta_hit
            region.hit_to = region.hit_to + delta_hit


def adjust_positions_sequences_center(
    gb_records: list[GenBankRecord], size_longest_sequence: int
) -> None:
    """Adjust position of sequences to the center including data in CDSs."""
    for record in gb_records:
        shift = (size_longest_sequence - record.length) / 2
        record.sequence_start = record.sequence_start + shift
        record.sequence_end = record.sequence_end + shift
        for sequence in record.cds:
            sequence.start = sequence.start + shift
            sequence.end = sequence.end + shift


def adjust_positions_alignments_center(
    alignments: list[BlastnAlignment], size_longest_sequence: int
) -> None:
    """Adjust position of alignmets to the center."""
    for alignment in alignments:
        shift_q = (size_longest_sequence - alignment.query_len) / 2
        shift_h = (size_longest_sequence - alignment.hit_len) / 2
        for region in alignment.regions:
            region.query_from = region.query_from + shift_q
            region.query_to = region.query_to + shift_q
            region.hit_from = region.hit_from + shift_h
            region.hit_to = region.hit_to + shift_h


# = FUNCTIONS FOR PLOTTING ==============================================================
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
            customdata=[f"min_identity={min_value}", f"max_identity={max_value}"],
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
    color="blue",
    name=None,
    custom_data=None,
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
            customdata=custom_data,
            hoverlabel=dict(font_size=14),
            hovertemplate="%{text}<extra></extra>",
        )
    )


def plot_line(
    fig: Figure,
    x_values: list,
    y_values: list,
    name=None,
    custom_data=[],
    color="black",
) -> None:
    """Plot line representing DNA sequences"""
    fig.add_trace(
        go.Scatter(
            x=x_values,
            y=y_values,
            mode="lines",
            name=name,
            customdata=custom_data,
            line=dict(color=color, width=4),
            hoverlabel=dict(font_size=14),
            hovertemplate=f"{name}<extra></extra>",
        )
    )


def plot_dna_sequences(
    fig: Figure,
    gb_records: list[GenBankRecord],
    longest_sequence: int,
    y_separation: int = 10,
) -> Figure:
    """Plot DNA sequences"""
    y_distance = len(gb_records) * y_separation
    for gb_record in gb_records:
        x1 = gb_record.sequence_start
        x2 = gb_record.sequence_end
        accession = gb_record.accession
        record_name = gb_record.name
        file_name = gb_record.file_name
        x_values = np.array([x1, x2])
        y_values = np.array([y_distance, y_distance])
        trace_name = f"Sequence: {record_name}"
        plot_line(
            fig,
            x_values,
            y_values,
            name=trace_name,
            custom_data=[accession, record_name, file_name, longest_sequence],
        )
        y_distance -= y_separation
    return fig


def plot_genes(
    fig: Figure,
    gb_records: list[GenBankRecord],
    name_from: str = "product",
    y_separation: int = 10,
) -> Figure:
    """Plot arrows representing genes.

    The custom_data has information for annotating the plot dinamically in Dash.
    """
    # Number of gb_records
    number_gb_records = len(gb_records)
    # Position of the first DNA sequence in the y axis. Plotting starts at the top.
    y = number_gb_records * y_separation
    # Find longest sequence
    longest_sequence = get_longest_sequence(gb_records)
    # Ratio head_height vs lenght of longest sequence
    ratio = 0.02
    head_height = longest_sequence * ratio
    # Iterate over GenBankRecords to plot genes.
    for i, record in enumerate(gb_records):
        for cds in record.cds:
            x1 = cds.start
            x2 = cds.end
            color = cds.color
            arrow = Arrow(x1=x1, x2=x2, y=y, head_height=head_height)
            x_values, y_values = arrow.get_coordinates()
            # Get name and check if is None
            name = cds.product if name_from == "product" else cds.gene
            name = name if name is not None else "no name"
            plot_polygon(
                fig,
                x_values,
                y_values,
                color=color,
                name=name,
                custom_data=[
                    "type=gene_info",
                    f"gb_record={i}",
                    f"number_gb_records={number_gb_records}",
                    f"x_start={cds.start}",
                    f"x_end={cds.end}",
                    f"y={y}",
                ],
            )
        y -= y_separation
    return fig


def plot_homology_regions(
    fig: Figure,
    alignments: list[BlastnAlignment],
    y_separation: int = 10,
    homology_padding: float = 1.1,
    colorscale: str = "Greys",
    straight_heights: bool = True,
    minimum_homology_length: int = 0,
    set_colorscale_to_extreme_homologies: bool = False,
    lowest_homology: None | float = None,
    highest_homology: None | float = None,
) -> Figure:
    """Plot homology regions of aligned sequences

    The customdata has information to change the colors dinamically in Dash
    """
    # Get the y distance to start plotting at the top of the graph
    y_distances = (len(alignments) + 1) * y_separation
    # Iterate over alignmets
    for alignment in alignments:
        # On each alignment iterate over the homologous regions for plotting
        for region in alignment.regions:
            # Get region coordinates
            x1 = region.query_from
            x2 = region.query_to
            x3 = region.hit_to
            x4 = region.hit_from
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
            homology = region.homology
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
                custom_data=[
                    "type=homology_info",
                    f"Identity={homology}",
                    f"homology_length={homology_length}",
                ],
            )
        y_distances -= y_separation
    return fig


def annotate_genes_top(
    fig: Figure, gb_records: list[GenBankRecord], y_separation: int = 10
) -> Figure:
    """Annotate genes of top sequence."""
    gb_record = gb_records[0]
    y_distance = len(gb_records) * y_separation
    for gene in gb_record.cds:
        x_position = (gene.start + gene.end) / 2
        y_position = y_distance
        name = gene.product
        fig.add_annotation(
            x=x_position,
            y=y_position + 1.1,
            text=name,
            showarrow=False,
            textangle=270,
            font=dict(size=16),
            xanchor="center",
            yanchor="bottom",
        )
    return fig


def annotate_genes_bottom(
    fig: Figure, gb_records: list[GenBankRecord], y_separation: int = 10
) -> Figure:
    """Annotate genes of top sequence."""
    gb_record = gb_records[len(gb_records) - 1]
    for gene in gb_record.cds:
        x_position = (gene.start + gene.end) / 2
        y_position = y_separation
        name = gene.product
        fig.add_annotation(
            x=x_position,
            y=y_position - 1.1,
            text=name,
            showarrow=False,
            textangle=270,
            font=dict(size=16),
            xanchor="center",
            yanchor="top",
        )
    return fig


def annotate_dna_sequences(
    fig: Figure,
    gb_records: list[GenBankRecord],
    longest_sequence: int,
    sequence_name: str = "accession",
    y_separation: int = 10,
    padding: int = 10,
) -> Figure:
    """Annotate DNA sequences."""
    y_distance = len(gb_records) * y_separation
    for record in gb_records:
        if sequence_name == "accession":
            sequence_name = record.accession
        if sequence_name == "name":
            sequence_name = record.name
        if sequence_name == "fname":
            sequence_name = record.file_name
        fig.add_annotation(
            x=longest_sequence + padding,
            # x=0.95,
            # xref="paper",
            xref="x",
            y=y_distance,
            text=sequence_name,
            font=dict(size=18),
            showarrow=False,
            xanchor="left",
            yanchor="middle",
        )
        y_distance -= y_separation
    return fig


def annotate_dna_sequences_using_trace_customdata(
    figure: Figure, annotate_with: str, padding: int = 10
) -> Figure:
    """Annotate dna sequences using trace customdata

    The trace customdata is created when plotting dna sequences and has the following
    information:
    custom_data: list[str, str, str, int] = [
        accession, record_name, file_name, longest_sequence
    ]
    """
    option = ["accession", "name", "fname"].index(annotate_with)
    for trace in figure["data"]:
        if ("name" in trace) and ("Sequence:" in trace["name"]):
            text = trace["customdata"][option]
            figure.add_annotation(
                x=trace["customdata"][3] + padding,
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
    for trace in figure["data"]:
        if trace["customdata"] is None:
            continue
        if "gene_info" in trace["customdata"][0]:
            gb_record = int(trace["customdata"][1].split("=")[1])
            # If gb_record is not the top (i. e. zero), continue.
            if gb_record != 0:
                continue
            x_start = float(trace["customdata"][3].split("=")[1])
            x_end = float(trace["customdata"][4].split("=")[1])
            y = float(trace["customdata"][5].split("=")[1])
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
    for trace in figure["data"]:
        if trace["customdata"] is None:
            continue
        if "gene_info" in trace["customdata"][0]:
            gb_record = int(trace["customdata"][1].split("=")[1])
            number_gb_records = int(trace["customdata"][2].split("=")[1])
            # if gb_record is not the last one, continue
            if gb_record != (number_gb_records - 1):
                continue
            # print(trace)
            x_start = float(trace["customdata"][3].split("=")[1])
            x_end = float(trace["customdata"][4].split("=")[1])
            y = float(trace["customdata"][5].split("=")[1])
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
        if ("customdata" in trace) and ("homology_info" in trace["customdata"][0]):
            # Get identity information from customdata
            identity = trace["customdata"][1]
            identity = float(identity.split("=")[1])
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


def make_alignments(plot_parameters: dict) -> tuple[BlastnAlignment, GenBankRecord]:
    """BLAST nucleotide sequences to make alignments."""
    # Convert paths of input_files and ouput_folder into Path objects
    input_files = [Path(input_file) for input_file in plot_parameters["input_files"]]
    output_folder = Path(plot_parameters["output_folder"])
    # Create fasta files for BLASTing using the gb files
    faa_files = genbank.make_fasta_files(input_files, output_folder)
    # Run blastn locally to make alignments.
    blast_xml_results = genbank.run_blastn(faa_files, output_folder)
    # Delete the fasta files used for BLASTing.
    misc.delete_files(faa_files)
    # Make a list of dict representing `BlastnAlignment` classes from the xml blatn
    # results.
    alignments = genbank.get_alignment_records(blast_xml_results)
    # Delete the xml documents.
    misc.delete_files(blast_xml_results)
    # Make a list of `GenBankRecord` classes using the gb files.
    gb_records = genbank.get_gb_records(input_files)
    return alignments, gb_records


# ====================================================================================== #
#                                  MAIN FUNCTION GUI                                     #
# ====================================================================================== #
def make_figure(plot_parameters: dict) -> tuple[Figure, float, float]:
    """Make a multiple sequence alignment plot

    Return
    ------
    Figure : plotly.graph_objects.Figure
    lowest_identity : float
    highest_identity : float
    """
    # = Get alignments and gb_records ===================================================
    alignments, gb_records = make_alignments(plot_parameters)

    # Before plotting, we must check the alignments position option to adjust the
    # coordinates of the sequences and genes.
    longest_sequence = get_longest_sequence(gb_records)
    if plot_parameters["alignments_position"] == "right":
        adjust_positions_alignments_right(alignments, longest_sequence)
        adjust_positions_sequences_right(gb_records, longest_sequence)
    if plot_parameters["alignments_position"] == "center":
        adjust_positions_alignments_center(alignments, longest_sequence)
        adjust_positions_sequences_center(gb_records, longest_sequence)

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

    # Find lowest and highest homologies
    lowest_identity, highest_identity = find_lowest_and_highest_homology(alignments)
    set_colorscale_to_extreme_homologies = plot_parameters[
        "set_colorscale_to_extreme_homologies"
    ]

    # Select colormap and its range to plot the homology regions
    colorscale = get_truncated_colorscale(
        colorscale_name=plot_parameters["identity_color"],
        vmin=plot_parameters["colorscale_vmin"],
        vmax=plot_parameters["colorscale_vmax"],
    )

    # Plot the DNA sequences
    fig = plot_dna_sequences(fig, gb_records, longest_sequence)
    # Plot the homology regions
    straight_homology_regions = (
        True if plot_parameters["straight_homology_regions"] == "straight" else False
    )
    fig = plot_homology_regions(
        fig,
        alignments,
        colorscale=colorscale,
        straight_heights=straight_homology_regions,
        minimum_homology_length=plot_parameters["minimum_homology_length"],
        set_colorscale_to_extreme_homologies=set_colorscale_to_extreme_homologies,
        lowest_homology=lowest_identity,
        highest_homology=highest_identity,
    )
    # Plot the genes
    fig = plot_genes(fig, gb_records, name_from=plot_parameters["annotate_genes_with"])
    # Annotate genes
    if plot_parameters["annotate_genes"] == "top":
        fig = annotate_genes_top_using_trace_customdata(fig)
    if plot_parameters["annotate_genes"] == "bottom":
        fig = annotate_genes_bottom_using_trace_customdata(fig)
    if plot_parameters["annotate_genes"] == "top-bottom":
        fig = annotate_genes_top_using_trace_customdata(fig)
        fig = annotate_genes_bottom_using_trace_customdata(fig)
    # Annotate DNA sequences
    if (sequence_name := plot_parameters["annotate_sequences"]) != "no":
        # fig = annotate_dna_sequences(fig, gb_records, longest_sequence, sequence_name)
        fig = annotate_dna_sequences_using_trace_customdata(fig, sequence_name)

    # Plot DNA scale
    fig = plot_scale(fig, longest_sequence, plot_parameters["add_scale_bar"])

    # Plot colorscale legend
    fig = plot_colorbar_legend(
        fig,
        colorscale,
        lowest_identity,
        highest_identity,
        set_colorscale_to_extreme_homologies=set_colorscale_to_extreme_homologies,
    )

    return fig, lowest_identity, highest_identity
