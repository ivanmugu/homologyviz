"""
Utilities for processing GenBank files and BLASTn results in HomologyViz.

This module provides functions to:
    - Convert GenBank (.gb) files to FASTA format for BLASTn (`make_fasta_files`).
    - Run local BLASTn alignments and capture XML results (`run_blastn`,
      `blastn_command_line`).
    - Parse BLASTn XML into structured DataFrames (`get_blast_metadata`,
      `parse_blast_record`).
    - Extract sequence- and feature-level metadata from GenBank records
      (`genbank_files_metadata_to_dataframes`, `parse_genbank_cds_to_df`).
    - Determine the longest sequence and homology bounds for plotting
      (`get_longest_sequence_dataframe`, `find_lowest_and_highest_homology_dataframe`).

These utilities underpin the data preparation pipeline for visualizing homology and gene
annotations.

Notes
-----
- This file is part of HomologyViz
- BSD 3-Clause License
- Copyright (c) 2024, Iván Muñoz Gutiérrez
"""

from pathlib import Path
import subprocess

import pandas as pd
from pandas import DataFrame
from _plotly_utils.basevalidators import ColorValidator

from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.SeqRecord import SeqRecord
from Bio.Blast import Record


def make_fasta_files(gb_files: list[Path], output_path: Path) -> list[Path]:
    """
    Convert GenBank files to FASTA format for downstream processing (e.g., BLASTn).

    Each input GenBank file is parsed to extract its sequence and metadata,
    and a corresponding FASTA (.faa) file is written to the specified output folder.

    Parameters
    ----------
    gb_files : list of pathlib.Path
        List of paths to GenBank (.gb) files to be converted.
    output_path : pathlib.Path
        Directory where the resulting FASTA files will be saved.

    Returns
    -------
    list of pathlib.Path
        List of paths to the generated FASTA (.faa) files.
    """
    # Initiate list to store paths to fasta files.
    faa_files = []
    # Iterate over paths of gb files.
    for gb_file in gb_files:
        # Read gb files and make a new record
        record = SeqIO.read(gb_file, "genbank")
        new_record = SeqRecord(record.seq, id=record.id, description=record.description)
        # Get name of gb file without extension
        name = gb_file.name.split(".")[0]
        faa_name = name + ".faa"
        # Make otuput path
        output_file = output_path / faa_name
        # Create fasta file
        SeqIO.write(new_record, output_file, "fasta")
        # Append path of fasta file to faa_files list.
        faa_files.append(output_file)

    return faa_files


def run_blastn(faa_files: list[Path], output_path: Path) -> list[Path]:
    """
    Run local BLASTn alignments between consecutive FASTA files and save results in XML
    format.

    For a given list of FASTA files, this function performs pairwise comparisons in order:
    file[0] vs file[1], file[1] vs file[2], and so on. The results are saved as XML files
    using BLAST output format 5.

    Parameters
    ----------
    faa_files : list of pathlib.Path
        List of paths to nucleotide FASTA (.faa) files to be compared.
    output_path : pathlib.Path
        Directory where the resulting BLASTn XML output files will be saved.

    Returns
    -------
    list of pathlib.Path
        List of paths to the BLASTn result files in XML format.
    """
    # TODO: Consider logging to file instead of printing directly

    # Initiate list to store paths to xml results.
    results = []
    # Iterate over paths of fasta files.
    for i in range(len(faa_files) - 1):
        # Make path to outpu file
        output_file_name = "result" + str(i) + ".xml"
        output_file = output_path / output_file_name
        # Run blastn
        std = blastn_command_line(
            query=faa_files[i], subject=faa_files[i + 1], outfmt=5, out=output_file
        )
        # Append path to xlm results to the result list
        results.append(output_file)
        print(f"BLASTing {faa_files[i]} (query) and {faa_files[i + 1]} (subject)\n")
        print(std)

    return results


def blastn_command_line(query: Path, subject: Path, out: Path, outfmt: int = 5) -> str:
    """
    Run a local BLASTn alignment between two nucleotide sequences using the command line.

    Executes BLASTn with the given query and subject FASTA files, writes results to the
    specified output file, and returns the standard output or error message.

    Parameters
    ----------
    query : pathlib.Path
        Path to the query FASTA file.
    subject : pathlib.Path
        Path to the subject FASTA file.
    out : pathlib.Path
        Path to the file where BLASTn output will be written.
    outfmt : int, default=5
        BLAST output format (5 = XML). HomologyViz requires XML format for parsing.

    Returns
    -------
    str
        The standard output from the BLASTn command if successful, or the error message
        if the command fails.

    Notes
    -----
    - Both `query` and `subject` must be valid nucleotide FASTA files.
    - The default output format (XML) is required for compatibility with HomologyViz.
    """
    # TODO: Add a raise instead of print() for critical failures if you're building a CLI
    # or GUI around this.

    # Define the BLASTn command
    command = [
        "blastn",
        "-query",
        str(query),
        "-subject",
        str(subject),
        "-outfmt",
        str(outfmt),
        "-out",
        str(out),
    ]

    # Run BLASTn
    try:
        result = subprocess.run(command, capture_output=True, text=True, check=True)
        return result.stdout
    except subprocess.CalledProcessError as e:
        print("Error running BLASTn:", e)
        return e.stderr


def genbank_files_metadata_to_dataframes(
    gb_files: list[Path],
) -> tuple[DataFrame, DataFrame]:
    """
    Parse GenBank files and return sequence and CDS metadata as structured DataFrames.

    This function reads a list of GenBank files and extracts relevant metadata for
    downstream plotting or analysis. It separates the data into two related tables:
    one for general sequence information (`gb_df`) and one for coding sequences
    (`cds_df`).

    Parameters
    ----------
    gb_files : list of pathlib.Path
        List of paths to GenBank (.gb) files.

    Returns
    -------
    gb_df : pandas.DataFrame
        DataFrame with GenBank record-level metadata, including:
            - file number, file path, file name, custom_name, record name, accession,
              sequence length,
              and plotting coordinates (`sequence_start`, `sequence_end`).

    cds_df : pandas.DataFrame
        DataFrame with CDS (gene) feature metadata from all GenBank files, including:
            - file number, accession, gene name, product name, strand, color (if
              available), and plotting coordinates (`start_plot`, `end_plot`).

    Notes
    -----
    These two DataFrames are linked via the `file_number` and `accession` fields.
    """
    # headers gb_files_df
    headers_gb_files_df = [
        "file_number",
        "file_path",
        "file_name",
        "custom_name",
        "record_name",
        "accession",
        "length",
        "sequence_start",
        "sequence_end",
        "sequence_start_plot_left",
        "sequence_end_plot_left",
        "sequence_start_plot_center",
        "sequence_end_plot_center",
        "sequence_start_plot_right",
        "sequence_end_plot_right",
    ]
    # Initiate dictionary to store data
    gb_files_data = dict(
        file_number=[],
        file_path=[],
        file_name=[],
        custom_name=[],
        record_name=[],
        accession=[],
        length=[],
        sequence_start=[],
        sequence_end=[],
        sequence_start_plot_left=[],
        sequence_end_plot_left=[],
        sequence_start_plot_center=[],
        sequence_end_plot_center=[],
        sequence_start_plot_right=[],
        sequence_end_plot_right=[],
    )
    # Initiate a list of cds DataFrames
    cds_dataframes = []

    # Iterate over GenBank files
    for i, gb_file in enumerate(gb_files):
        # fill data related to the file
        gb_files_data["file_number"].append(i)
        gb_files_data["file_path"].append(gb_file)
        gb_files_data["file_name"].append(gb_file.stem)

        # Add an empty string to custom name for future manipulation in the GUI
        gb_files_data["custom_name"].append("")

        # Read the file into a temporary variable
        record = SeqIO.read(gb_file, "genbank")

        # fill data related to the GenBank record
        gb_files_data["record_name"].append(record.name)
        gb_files_data["accession"].append(record.id)
        seq_length = len(record)
        gb_files_data["length"].append(float(seq_length))
        gb_files_data["sequence_start"].append(0.0)
        gb_files_data["sequence_end"].append(float(seq_length))
        gb_files_data["sequence_start_plot_left"].append(0.0)
        gb_files_data["sequence_end_plot_left"].append(float(seq_length))

        # Placeholders for future alignment logic; these will be updated later in the
        # pipeline
        gb_files_data["sequence_start_plot_center"].append(None)
        gb_files_data["sequence_end_plot_center"].append(None)
        gb_files_data["sequence_start_plot_right"].append(None)
        gb_files_data["sequence_end_plot_right"].append(None)

        # Get a DataFrame from the cds
        cds_dataframes.append(
            parse_genbank_cds_to_df(
                record=record,
                file_number=i,
                accession=record.id,
            )
        )

    # Create the GenBank files DataFrame
    gb_df = DataFrame(gb_files_data, columns=headers_gb_files_df)
    # Concatenate the cds_dataframes list into a single DataFrame
    cds_df = pd.concat(cds_dataframes, ignore_index=True)

    gb_df, cds_df = include_coordinates_to_center_align_sequences(
        gb_records=gb_df,
        cds=cds_df,
    )
    gb_df, cds_df = include_coordinates_to_right_align_sequences(
        gb_records=gb_df,
        cds=cds_df,
    )

    return gb_df, cds_df


def get_color_qualifier(qualifiers: dict) -> str | None:
    for key, value in qualifiers.items():
        if key.lower() == "color" and value:
            return value[0]

    return None


def is_valid_color(color: str) -> bool:
    validator = ColorValidator("color", "homologyviz")

    try:
        validator.validate_coerce(color)
    except ValueError:
        return False

    return True


def parse_genbank_cds_to_df(
    record: SeqRecord, file_number: int, accession: str
) -> DataFrame:
    """
    Extract CDS feature metadata from a GenBank record and return it as a DataFrame.

    This function parses a `Bio.SeqRecord` GenBank object and compiles information
    from all its `CDS` features.

    Parameters
    ----------
    record : Bio.SeqRecord.SeqRecord
        The parsed GenBank record to extract CDS data from.
    file_number : int
        Index of the GenBank file, used for relational tracking.
    accession : str
        Accession ID of the sequence, used for relational grouping.

    Returns
    -------
    pandas.DataFrame
        A DataFrame with one row per CDS part, containing the following columns:
        - file_number: index of source file
        - cds_number: index of the CDS feature within the record
        - accession: accession ID of the sequence
        - gene: gene name (if available)
        - product: protein product name (if available)
        - start, end: sequence coordinates (1-based)
        - strand: strand orientation (+1 or -1)
        - color: gene color (from `/Color` qualifier, or default "#ffff00")
        - start_plot, end_plot: adjusted coordinates for plotting purposes

    Notes
    -----
    - Compound CDS locations are treated as a single CDS feature and represented
      by one row using the outermost coordinates of the full location.
    - Start/end coordinates are stored as floats for consistency with plotting tools.
    """
    # DataFrame headers
    headers = [
        "file_number",
        "cds_number",
        "accession",
        "gene",
        "product",
        "custom_name",
        "start",
        "end",
        "strand",
        "color",
        "start_plot_left",
        "end_plot_left",
        "start_plot_center",
        "end_plot_center",
        "start_plot_right",
        "end_plot_right",
    ]
    # Initiate dictionary to store data
    data = dict(
        file_number=[],
        cds_number=[],
        accession=[],
        gene=[],
        product=[],
        custom_name=[],
        start=[],
        end=[],
        strand=[],
        color=[],
        start_plot_left=[],
        end_plot_left=[],
        start_plot_center=[],
        end_plot_center=[],
        start_plot_right=[],
        end_plot_right=[],
    )
    # Initialize counter to track cds; enumerate will not give continues numbers.
    counter = 0

    # Iterate over features to extract data. Make sure that if there is no metadata,
    # then add None.
    for feature in record.features:
        if feature.type != "CDS":
            continue

        data["file_number"].append(file_number)
        data["cds_number"].append(counter)
        counter += 1

        data["accession"].append(accession)

        # Add an empty string to custom name for future manipulation in the GUI
        data["custom_name"].append("")

        if gene := feature.qualifiers.get("gene", None):
            data["gene"].append(gene[0])
        else:
            data["gene"].append(None)
        if product := feature.qualifiers.get("product", None):
            data["product"].append(product[0])
        else:
            data["product"].append(None)

        # CDS color obtained from a case-insensitive `color` qualifier.
        # Invalid or missing colors default to yellow (#ffff00).
        color = get_color_qualifier(feature.qualifiers)
        if color is not None and is_valid_color(color):
            data["color"].append(color)
        else:
            data["color"].append("#ffff00")

        # Treat each CDS feature as one DataFrame row, including CDSs represented by
        # compound locations such as join(...).
        location = feature.location
        strand = location.strand

        data["strand"].append(strand)

        if strand == -1:
            start = float(location.end)
            end = float(location.start + 1)
        else:
            start = float(location.start + 1)
            end = float(location.end)

        data["start"].append(start)
        data["start_plot_left"].append(start)
        data["end"].append(end)
        data["end_plot_left"].append(end)

        # Placeholders for future alignment logic; these will be updated later in the
        # pipeline
        data["start_plot_center"].append(None)
        data["start_plot_right"].append(None)
        data["end_plot_center"].append(None)
        data["end_plot_right"].append(None)

    return DataFrame(data, columns=headers)


def get_blast_metadata(
    xml_alignment_result: list[Path],
    size_longest_sequence: int,
) -> tuple[DataFrame, DataFrame]:
    """
    Parse BLASTn XML result files into structured Pandas DataFrames.

    This function extracts both summary and detailed region metadata from a list of
    BLASTn XML result files, returning two linked DataFrames:
    - `alignments_df`: High-level metadata for each BLAST alignment.
    - `regions_df`: Local matching regions for each alignment.

    Parameters
    ----------
    xml_alignment_result : list of pathlib.Path
        List of paths to XML-formatted BLASTn result files (outfmt=5).
    size_longest_sequence : int
        Length of the longest sequence in the alignment set. Used to calculate plotting
        coordinates for center- and right-aligned sequences.

    Returns
    -------
    alignments_df : pandas.DataFrame
        Summary table with one row per alignment. Columns include:
            - alignment_number (int): Unique alignment index
            - query_name (str): Query sequence ID
            - hit_name (str): Subject sequence ID
            - query_len (int): Query sequence length
            - hit_len (int): Subject sequence length

    regions_df : pandas.DataFrame
        Detailed region-level metadata for all matching regions across alignments.
        Includes start/end positions and identity metrics. Each row corresponds to
        one HSP (high-scoring pair). The `alignment_number` field links this table
        to `alignments_df`.

    Notes
    -----
    For additional region-level metadata, see the `parse_blast_record` function.
    """
    headers = ["alignment_number", "query_name", "hit_name", "query_len", "hit_len"]
    data = dict(
        alignment_number=[],
        query_name=[],
        hit_name=[],
        query_len=[],
        hit_len=[],
    )
    regions = []
    # Iterate over xml files containing alignment results
    for i, xml_file in enumerate(xml_alignment_result):
        with open(xml_file, "r") as result_handle:
            blast_record = NCBIXML.read(result_handle)
            # Add alignment number for a relational database
            data["alignment_number"].append(i)
            # Get metadata
            data["query_name"].append(blast_record.query)
            data["hit_name"].append(blast_record.alignments[0].hit_def)
            data["query_len"].append(int(blast_record.query_length))
            data["hit_len"].append(int(blast_record.alignments[0].length))
            regions.append(
                parse_blast_record(
                    blast_record=blast_record,
                    alignment_number=i,
                )
            )
    # Create DataFrame
    alignments_df = DataFrame(data, columns=headers)
    regions_df = pd.concat(regions, ignore_index=True)

    alignments_df, regions_df = include_coordinates_to_center_align_alignments(
        alignments=alignments_df,
        regions=regions_df,
        size_longest_sequence=size_longest_sequence,
    )
    alignments_df, regions_df = include_coordinates_to_right_align_alignments(
        alignments=alignments_df,
        regions=regions_df,
        size_longest_sequence=size_longest_sequence,
    )

    return alignments_df, regions_df


def parse_blast_record(blast_record: Record, alignment_number: int) -> DataFrame:
    """
    Parse HSP metadata from a pairwise BLAST record.

    This function extracts coordinate and similarity information for each high-scoring
    pair (HSP) associated with the first BLAST hit in the record. In HomologyViz, each
    BLAST record represents a pairwise comparison between two adjacent sequences in the user-defined sequence order.

    Raw BLAST coordinates are preserved, and plotting coordinates are initialized for
    left-aligned visualization. Center- and right-aligned plotting coordinates are left
    empty for later calculation.

    Parameters
    ----------
    blast_record : Bio.Blast.Record
        Parsed BLAST record produced by ``Bio.Blast.NCBIXML.read``. The record must
        contain at least one alignment.
    alignment_number : int
        Identifier assigned to the alignment. The same value is stored for every HSP
        extracted from the record. This number is used to link the HSPs to their parent
        alignment in a relational DataFrame.

    Returns
    -------
    pandas.DataFrame
        One row per HSP from the first BLAST alignment, with the following
        columns:

        alignment_number : int
            Identifier linking each HSP to its parent alignment.
        query_from, query_to : float
            Start and end coordinates of the HSP in the query sequence.
        hit_from, hit_to : float
            Start and end coordinates of the HSP in the subject sequence.
        query_from_plot_left, query_to_plot_left : float
            Query plotting coordinates for left-aligned visualization.
        hit_from_plot_left, hit_to_plot_left : float
            Subject plotting coordinates for left-aligned visualization.
        query_from_plot_center, query_to_plot_center : float or None
            Query plotting coordinates for center alignment. Initialized as
            missing and calculated later in the pipeline.
        hit_from_plot_center, hit_to_plot_center : float or None
            Subject plotting coordinates for center alignment. Initialized as
            missing and calculated later in the pipeline.
        query_from_plot_right, query_to_plot_right : float or None
            Query plotting coordinates for right alignment. Initialized as
            missing and calculated later in the pipeline.
        hit_from_plot_right, hit_to_plot_right : float or None
            Subject plotting coordinates for right alignment. Initialized as
            missing and calculated later in the pipeline.
        identity : int
            Number of identical positions in the HSP.
        positive : int
            Number of positive-scoring positions in the HSP.
        align_len : int
            Length of the HSP alignment.
        homology : float
            Fraction of identical positions, calculated as
            ``identity / align_len``.
    """
    headers = [
        "alignment_number",
        "query_from",
        "query_to",
        "query_from_plot_left",
        "query_to_plot_left",
        "query_from_plot_center",
        "query_to_plot_center",
        "query_from_plot_right",
        "query_to_plot_right",
        "hit_from",
        "hit_to",
        "hit_from_plot_left",
        "hit_to_plot_left",
        "hit_from_plot_center",
        "hit_to_plot_center",
        "hit_from_plot_right",
        "hit_to_plot_right",
        "identity",
        "positive",
        "align_len",
        "homology",
    ]
    data = dict(
        alignment_number=[],
        query_from=[],
        query_to=[],
        query_from_plot_left=[],
        query_to_plot_left=[],
        query_from_plot_center=[],
        query_to_plot_center=[],
        query_from_plot_right=[],
        query_to_plot_right=[],
        hit_from=[],
        hit_to=[],
        hit_from_plot_left=[],
        hit_to_plot_left=[],
        hit_from_plot_center=[],
        hit_to_plot_center=[],
        hit_from_plot_right=[],
        hit_to_plot_right=[],
        identity=[],
        positive=[],
        align_len=[],
        homology=[],
    )
    for region in blast_record.alignments[0].hsps:
        data["alignment_number"].append(alignment_number)
        data["query_from"].append(float(region.query_start))
        data["query_to"].append(float(region.query_end))
        data["query_from_plot_left"].append(float(region.query_start))
        data["query_to_plot_left"].append(float(region.query_end))
        data["hit_from"].append(float(region.sbjct_start))
        data["hit_to"].append(float(region.sbjct_end))
        data["hit_from_plot_left"].append(float(region.sbjct_start))
        data["hit_to_plot_left"].append(float(region.sbjct_end))
        data["identity"].append(int(region.identities))
        data["positive"].append(int(region.positives))
        data["align_len"].append(int(region.align_length))
        homology = int(region.identities) / int(region.align_length)
        data["homology"].append(homology)

        # Placeholders for future alignment logic; these will be updated later in the
        # pipeline
        data["query_from_plot_center"].append(None)
        data["query_to_plot_center"].append(None)
        data["query_from_plot_right"].append(None)
        data["query_to_plot_right"].append(None)
        data["hit_from_plot_center"].append(None)
        data["hit_to_plot_center"].append(None)
        data["hit_from_plot_right"].append(None)
        data["hit_to_plot_right"].append(None)

    return DataFrame(data, columns=headers)


def include_coordinates_to_center_align_sequences(
    gb_records: DataFrame,
    cds: DataFrame,
) -> tuple[DataFrame, DataFrame]:
    """
    Center-align sequences and their CDS features for plotting.

    Each sequence is shifted horizontally relative to the longest sequence so that
    all sequences are centered on the same plotting range. The same shift is
    applied to every CDS belonging to that sequence.

    Parameters
    ----------
    gb_records : pandas.DataFrame
        Sequence metadata. Must contain ``file_number``, ``length``,
        ``sequence_start``, and ``sequence_end``.
    cds : pandas.DataFrame
        CDS metadata. Must contain ``file_number``, ``start``, and ``end``.

    Returns
    -------
    tuple[pandas.DataFrame, pandas.DataFrame]
        Updated sequence and CDS DataFrames. Sequence coordinates are stored in
        ``sequence_start_plot_center`` and ``sequence_end_plot_center``. CDS
        coordinates are stored in ``start_plot_center`` and ``end_plot_center``.
    """
    size_longest_sequence = gb_records["length"].max()
    # Iterate over gb_records rows to find the shift value
    for i, row in gb_records.iterrows():
        # Get value to shift sequences to the center
        shift = (size_longest_sequence - row["length"]) / 2
        file_number = row["file_number"]
        cds_mask = cds["file_number"] == file_number

        # Change the values of the sequence_start and sequence_end of gb_records
        gb_records.loc[i, "sequence_start_plot_center"] = row["sequence_start"] + shift
        gb_records.loc[i, "sequence_end_plot_center"] = row["sequence_end"] + shift

        # Change the values of start_plot and end_plot of the cds DataFrame
        cds.loc[cds_mask, "start_plot_center"] = cds.loc[cds_mask, "start"] + shift
        cds.loc[cds_mask, "end_plot_center"] = cds.loc[cds_mask, "end"] + shift

    return gb_records, cds


def include_coordinates_to_right_align_sequences(
    gb_records: DataFrame,
    cds: DataFrame,
) -> tuple[DataFrame, DataFrame]:
    """
    Right-align sequences and their CDS features for plotting.

    Each sequence is shifted horizontally relative to the longest sequence so that
    all sequences share the same right boundary. The same shift is applied to
    every CDS belonging to that sequence.

    Parameters
    ----------
    gb_records : pandas.DataFrame
        Sequence metadata. Must contain ``file_number``, ``length``,
        ``sequence_start``, and ``sequence_end``.
    cds : pandas.DataFrame
    CDS metadata. Must contain ``file_number``, ``start``, and ``end``.

    Returns
    -------
    tuple[pandas.DataFrame, pandas.DataFrame]
        Updated sequence and CDS DataFrames. Sequence coordinates are stored in
        ``sequence_start_plot_right`` and ``sequence_end_plot_right``. CDS
        coordinates are stored in ``start_plot_right`` and ``end_plot_right``.
    """
    size_longest_sequence = gb_records["length"].max()

    # Iterate over gb_records rows to find the shift value
    for i, row in gb_records.iterrows():
        # Get value to shift sequences to the center
        shift = size_longest_sequence - row["length"]
        file_number = row["file_number"]
        cds_mask = cds["file_number"] == file_number

        # Change the values of the sequence_start and sequence_end of gb_records
        gb_records.loc[i, "sequence_start_plot_right"] = row["sequence_start"] + shift
        gb_records.loc[i, "sequence_end_plot_right"] = row["sequence_end"] + shift

        # Change the values of start_plot and end_plot of the cds DataFrame
        cds.loc[cds_mask, "start_plot_right"] = cds.loc[cds_mask, "start"] + shift
        cds.loc[cds_mask, "end_plot_right"] = cds.loc[cds_mask, "end"] + shift

    return gb_records, cds


def include_coordinates_to_center_align_alignments(
    alignments: DataFrame,
    regions: DataFrame,
    size_longest_sequence: int,
) -> tuple[DataFrame, DataFrame]:
    """
    Center-align BLAST alignment regions for plotting.

    For each alignment, this function independently shifts the query and hit
    coordinates relative to the longest sequence so that both sequences are
    centered within the same plotting range. The calculated shifts are applied
    to all HSP regions associated with the alignment.

    Parameters
    ----------
    alignments : pandas.DataFrame
        BLAST alignment summary metadata. Must contain ``alignment_number``,
        ``query_len``, and ``hit_len``.
    regions : pandas.DataFrame
        BLAST HSP metadata. Must contain ``alignment_number``, ``query_from``,
        ``query_to``, ``hit_from``, and ``hit_to``.
    size_longest_sequence : int
        Length of the longest sequence in the sequence set. Used as the common
        plotting width for calculating the query and hit shifts.

    Returns
    -------
    tuple[pandas.DataFrame, pandas.DataFrame]
        Updated alignment and region DataFrames. Center-aligned HSP coordinates
        are stored in ``query_from_plot_center``, ``query_to_plot_center``,
        ``hit_from_plot_center``, and ``hit_to_plot_center``.
    """
    # Iterate over alignments to find the shift value
    for _, alignment in alignments.iterrows():
        # Find the amount to add to shift the alignments to the center.
        shift_q = (size_longest_sequence - alignment["query_len"]) / 2
        shift_h = (size_longest_sequence - alignment["hit_len"]) / 2
        alignment_number = alignment["alignment_number"]
        alignment_mask = regions["alignment_number"] == alignment_number

        # Change the values of the regions used for plotting.
        regions.loc[alignment_mask, "query_from_plot_center"] = (
            regions.loc[alignment_mask, "query_from"] + shift_q
        )
        regions.loc[alignment_mask, "query_to_plot_center"] = (
            regions.loc[alignment_mask, "query_to"] + shift_q
        )
        regions.loc[alignment_mask, "hit_from_plot_center"] = (
            regions.loc[alignment_mask, "hit_from"] + shift_h
        )
        regions.loc[alignment_mask, "hit_to_plot_center"] = (
            regions.loc[alignment_mask, "hit_to"] + shift_h
        )

    return alignments, regions


def include_coordinates_to_right_align_alignments(
    alignments: DataFrame,
    regions: DataFrame,
    size_longest_sequence: int,
) -> tuple[DataFrame, DataFrame]:
    """
    Right-align BLAST alignment regions for plotting.

    For each alignment, this function independently shifts the query and hit
    coordinates relative to the longest sequence so that both sequences share
    the same right plotting boundary. The calculated shifts are applied to all
    HSP regions associated with the alignment.

    Parameters
    ----------
    alignments : pandas.DataFrame
        BLAST alignment summary metadata. Must contain ``alignment_number``,
        ``query_len``, and ``hit_len``.
    regions : pandas.DataFrame
        BLAST HSP metadata. Must contain ``alignment_number``, ``query_from``,
        ``query_to``, ``hit_from``, and ``hit_to``.
    size_longest_sequence : int
        Length of the longest sequence in the sequence set. Used as the common
        plotting width for calculating the query and hit shifts.

    Returns
    -------
    tuple[pandas.DataFrame, pandas.DataFrame]
        Updated alignment and region DataFrames. Right-aligned HSP coordinates
        are stored in ``query_from_plot_right``, ``query_to_plot_right``,
        ``hit_from_plot_right``, and ``hit_to_plot_right``.
    """
    # Iterate over alignments to find the shift value
    for _, alignment in alignments.iterrows():
        # Find the amount to add to shift the alignments to the right.
        delta_query = size_longest_sequence - alignment["query_len"]
        delta_hit = size_longest_sequence - alignment["hit_len"]
        alignment_number = alignment["alignment_number"]
        alignment_mask = regions["alignment_number"] == alignment_number

        # Change the values fo the regions used for plotting.
        regions.loc[alignment_mask, "query_from_plot_right"] = (
            regions.loc[alignment_mask, "query_from"] + delta_query
        )
        regions.loc[alignment_mask, "query_to_plot_right"] = (
            regions.loc[alignment_mask, "query_to"] + delta_query
        )
        regions.loc[alignment_mask, "hit_from_plot_right"] = (
            regions.loc[alignment_mask, "hit_from"] + delta_hit
        )
        regions.loc[alignment_mask, "hit_to_plot_right"] = (
            regions.loc[alignment_mask, "hit_to"] + delta_hit
        )

    return alignments, regions


if __name__ == "__main__":
    gb_records = pd.DataFrame(
        {
            "file_number": [5, 8, 12],
            "length": [1000, 600, 800],
            "sequence_start": [1, 1, 1],
            "sequence_end": [1000, 600, 800],
        }
    )

    cds = pd.DataFrame(
        {
            "file_number": [5, 8, 12],
            "start": [100, 100, 100],
            "end": [200, 200, 200],
        }
    )

    for i, row in gb_records.iterrows():
        cds_mask = cds["file_number"] == row["file_number"]
        print(f"i: {i}, row: \n{row}")
        print(f"cds_mask: \n{cds_mask}")
        print(f"type of cds_mask: {type(cds_mask)}")
    # _, cds = include_coordinates_to_center_align_sequences(
    #     gb_records=gb_records,
    #     cds=cds,
    # )
    # print(cds)
