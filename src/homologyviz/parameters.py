"""
Class for storing plot configuration and metadata during Dash callbacks.

This module defines the `PlotParameters` class, which centralizes all relevant data
needed for plotting in the HomologyViz app. It is used to share state between
callbacks and to retain user-defined options and metadata derived from GenBank
and BLASTn results.

License
-------
This file is part of HomologyViz
BSD 3-Clause License
Copyright (c) 2024, Iván Muñoz Gutiérrez
"""

from pathlib import Path
from pandas import DataFrame


class PlotParameters:
    """Store alignments information and user input for plotting in Dash.

    This class was designed to collect information during the Dash callbacks.

    Parameters
    ----------
    input_files : Path
        List of genbank files to BLAST.
    number_gb_records : int
        Number of gb files to BLAST.
    output_folder : Path
        Path to output any result. This can be the path to a temporary directory.
    alignments_position : str
        Position of the alignemts in the plot. Options are left, center, or right.
    identity_color : str
        Selected colormap to show the different shades that represent identities. For
        example `Greys`, `Greens`, and `Blues`.
    colorscale_vmin : float
        Minimum value to use in the colormap to represent identities. Values can go from
        0 to 1.0; for example, a value of 0.5 represents the shade at the center of the
        colormap.
    colorscale_vmax : float
        Maximum value to use in the colormap to represent identities. Values can go from
        0 to 1.0; for example, a value of 1.0 represents the shade with the highest value
        in the colormap.
    set_colorscale_to_extreme_homologies : bool
        If this parameter is set to True, the lowest and highest homologies will be
        represented by the values used in colorscale_vmin and colorscale_vmax,
        respectively. Otherwise, the lowest and highest homologies will be represented by
        0 and 1.0 in the colorsle by, respectively.
    annotate_sequences : str
        Annotate DNA sequences. Options:
        - no: no annotations
        - accession: use the accesion number
        - name: use the sequence name
        - fname: use the file name
    annotate_genes : str
        Annotate genes on the DNA sequences. Options:
        - no: no annotations
        - top: annotate only the genes at the top sequence.
        - bottom: annotate only the genes at the bottom sequence.
        - top-bottom: annotate only the genes at the top and bottom sequences.
    annotate_genes_with : str
        Annotate genes using GenBank file metadata stored in `CDS gene` or `CDS product`.
        Options are `gene` and `product`.
    stright_homology_regions : bool
        It this parameter is the to False, the shadows reprenting homologies with have a
        Bezier shape.
    minimium_homology_lenght : int
        This number represent the lenght of the minimum homology shown in the plot. For
        example, if it is set to 500, all homologies spanning 500 or more nucleotides are
        shown.
    add_scale_bar : str
        Show the scale bar in plot. Option are `yes` or `no`.
    selected_traces : list
        List of `curve_numbers` of selected traces. This list is used when the user selec
        genes in the `edit` tab to change their colors.
    lowest_identity : float
        Lowest identity in the BLASTn analysis.
    highest_identity : float
        Highest identity in the BLASTn analysis.
    longest_sequence : int
        Lenght of the longest sequence during the BLASTn analysis.
    gb_df : pandas.DataFrame
        Pandas DataFrame storing metadata of the GenBank files for plotting.
    cds_df : pandas.DataFrame
        Pandas DataFrame storing metadata of the GenBank CDS files for plotting.
    alignments_df : pandas.DataFrame
        Pandas DataFrame storing metadata of the BLASTn results for plotting.
    alignments_regions_df : pandas.DataFrame
        Pandas DataFrame stroing metadata from the homology regions found after BLASTning
        for plotting.
    draw_from_button : str
        Stores the name id of the button that triggered the callback for plotting. This
        parameter is import to distinguish between the `Draw` button and the rest of
        buttons used to update the plot. The id for the `Draw` button is `draw-button`.
    y_separation float
        Number to plot the sequences in the y-axis. The values to plot in the x-axis are
        stored in the different Pandas DataFrames.
    """

    def __init__(
        self,
        input_files: None | list[Path] = None,
        number_gb_records: None | int = None,
        output_folder: None | Path = None,
        alignments_position: None | str = None,
        identity_color: None | str = None,
        colorscale_vmin: None | float = None,
        colorscale_vmax: None | float = None,
        set_colorscale_to_extreme_homologies: None | bool = None,
        annotate_sequences: None | str = None,
        annotate_genes: None | str = None,
        annotate_genes_with: None | str = None,
        straight_homology_regions: None | bool = None,
        minimum_homology_length: None | int = None,
        add_scale_bar: None | str = None,
        selected_traces: None | list = None,
        lowest_identity: None | float = None,
        highest_identity: None | float = None,
        longest_sequence: None | int = None,
        gb_df: None | DataFrame = None,
        cds_df: None | DataFrame = None,
        alignments_df: None | DataFrame = None,
        alignments_regions_df: None | DataFrame = None,
        draw_from_button: None | str = None,
        y_separation: None | float = None,
    ):
        self.input_files = input_files
        self.number_gb_records = number_gb_records
        self.output_folder = output_folder
        self.alignments_position = alignments_position
        self.identity_color = identity_color
        self.colorscale_vmin = colorscale_vmin
        self.colorscale_vmax = colorscale_vmax
        self.set_colorscale_to_extreme_homologies = set_colorscale_to_extreme_homologies
        self.annotate_sequences = annotate_sequences
        self.annotate_genes = annotate_genes
        self.annotate_genes_with = annotate_genes_with
        self.straight_homology_regions = straight_homology_regions
        self.minimum_homology_length = minimum_homology_length
        self.add_scale_bar = add_scale_bar
        self.selected_traces = selected_traces
        self.lowest_identity = lowest_identity
        self.highest_identity = highest_identity
        self.longest_sequence = longest_sequence
        self.gb_df = gb_df
        self.cds_df = cds_df
        self.alignments_df = alignments_df
        self.alignments_regions_df = alignments_regions_df
        self.draw_from_button = draw_from_button
        self.y_separation = y_separation

    def reset(self):
        """Reset all attributes to their default values."""
        self.__init__()
