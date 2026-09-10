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
    - Each part of a multi-segment CDS is treated as a separate row.
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

        # TODO: allow user to provide "color", "Color", or "COLOR"
        if feature.qualifiers.get("Color", None):
            data["color"].append(feature.qualifiers["Color"][0])
        else:
            data["color"].append("#ffff00")  # Make yellow default color

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
    size_longest_sequence: int = None,
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

    output_folder = Path(
        "/Users/msp/Documents/coding/python_projects/HomologyViz/data/SW4848_paper"
    )
    alignments_df.to_csv(output_folder / "alignments_df.csv", index=False)
    regions_df.to_csv(output_folder / "regions_df.csv", index=False)

    return alignments_df, regions_df


def parse_blast_record(blast_record: Record, alignment_number: int) -> DataFrame:
    """
    Parse a BLAST record and extract metadata for all matching regions (HSPs).

    This function processes the first alignment in a BLAST record and extracts key
    information about each high-scoring pair (HSP), including coordinate ranges,
    identity metrics, and computed homology. It returns a DataFrame with one row
    per region, ready for downstream plotting or filtering.

    Parameters
    ----------
    blast_record : Bio.Blast.Record
        A parsed BLAST record object from Bio.Blast.NCBIXML.read().
        Must contain at least one alignment with HSPs.
    alignment_number : int
        Unique index for this alignment, used to link with the summary DataFrame.

    Returns
    -------
    pandas.DataFrame
        DataFrame where each row represents a BLAST high-scoring pair (HSP).
        Includes both raw start/end coordinates and pre-scaled values for plotting.

        Columns:
        - alignment_number : int
        - query_from, query_to : float
        - query_from_plot, query_to_plot : float
        - hit_from, hit_to : float
        - hit_from_plot, hit_to_plot : float
        - identity : int (number of identical matches)
        - positive : int (number of positive-scoring matches)
        - align_len : int (alignment length)
        - homology : float (identity / alignment length)
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
    Adjust plotting coordinates to center-align each sequence and its CDS features.

    This function horizontally centers all sequences relative to the longest sequence.
    It appends ``sequence_start_center`` and ``sequence_end_center`` columns to
    `gb_records`, and ``start_plot_center`` and ``end_plot_center`` columns to `cds`,
    which represent the adjusted coordinates for plotting.

    Parameters
    ----------
    gb_records : pandas.DataFrame
        DataFrame containing metadata for GenBank sequences. Must include:
        - 'length', 'sequence_start', 'sequence_end', and 'file_number'.

    cds : pandas.DataFrame
        DataFrame containing CDS feature metadata. Must include:
        - 'file_number', 'start_plot', and 'end_plot'.
    """
    size_longest_sequence = gb_records["length"].max()
    # Iterate over gb_records rows to find the shift value
    for i, row in gb_records.iterrows():
        # Get value to shift sequences to the center
        shift = (size_longest_sequence - row["length"]) / 2
        # Change the values of the sequence_start and sequence_end of gb_records
        gb_records.loc[i, "sequence_start_plot_center"] = row["sequence_start"] + shift
        gb_records.loc[i, "sequence_end_plot_center"] = row["sequence_end"] + shift
        # Change the values of start_plot and end_plot of the cds DataFrame
        cds.loc[cds["file_number"] == i, "start_plot_center"] = (
            cds.loc[cds["file_number"] == i, "start"] + shift
        )
        cds.loc[cds["file_number"] == i, "end_plot_center"] = (
            cds.loc[cds["file_number"] == i, "end"] + shift
        )

    return gb_records, cds


def include_coordinates_to_right_align_sequences(
    gb_records: DataFrame,
    cds: DataFrame,
) -> tuple[DataFrame, DataFrame]:
    """
    Adjust plotting coordinates to right-align sequences and CDS features.

    This function horizontally right-aligns each sequence relative to the longest
    sequence. It updates the `sequence_start_right` and `sequence_end_right` columns in
    `gb_records`, and adjusts the CDS plotting coordinates (`start_plot_right`,
    `end_plot_right`) in `cds`.

    Parameters
    ----------
    gb_records : pandas.DataFrame
        DataFrame containing GenBank sequence metadata. Must include:
        - 'length', 'sequence_start', 'sequence_end', and 'file_number'.

    cds : pandas.DataFrame
        DataFrame containing CDS feature metadata. Must include:
        - 'file_number', 'start_plot', and 'end_plot'.
    """
    size_longest_sequence = gb_records["length"].max()

    # Iterate over gb_records rows to find the shift value
    for i, row in gb_records.iterrows():
        # Get value to shift sequences to the center
        shift = size_longest_sequence - row["length"]
        # Change the values of the sequence_start and sequence_end of gb_records
        gb_records.loc[i, "sequence_start_plot_right"] = row["sequence_start"] + shift
        gb_records.loc[i, "sequence_end_plot_right"] = row["sequence_end"] + shift
        # Change the values of start_plot and end_plot of the cds DataFrame
        cds.loc[cds["file_number"] == i, "start_plot_right"] = (
            cds.loc[cds["file_number"] == i, "start"] + shift
        )
        cds.loc[cds["file_number"] == i, "end_plot_right"] = (
            cds.loc[cds["file_number"] == i, "end"] + shift
        )

    return gb_records, cds


def include_coordinates_to_center_align_alignments(
    alignments: DataFrame,
    regions: DataFrame,
    size_longest_sequence: int = None,
) -> tuple[DataFrame, DataFrame]:
    """
    Center-align alignment regions for plotting relative to the longest sequence.

    This function adjusts the plotting coordinates of each alignment region to center
    both the query and hit sequences. It shifts the `*_plot` columns
    (`query_from_plot_center`, `query_to_plot_center`, `hit_from_plot_center`
    `hit_to_plot_center`) based on the difference between each alignment's sequence length
    and the longest sequence in the dataset.

    Parameters
    ----------
    alignments : pandas.DataFrame
        DataFrame containing BLAST alignment summary metadata. Must include:
            - 'alignment_number', 'query_len', 'hit_len'.

    regions : pandas.DataFrame
        DataFrame containing BLAST alignment region metadata. Must include:
            - 'alignment_number', 'query_from_plot', 'query_to_plot',
              'hit_from_plot', 'hit_to_plot'.
    """
    # Iterate over alignments to find the shift value
    for i, alignment in alignments.iterrows():
        # Find the amount to add to shift the alignments the the center.
        shift_q = (size_longest_sequence - alignment["query_len"]) / 2
        shift_h = (size_longest_sequence - alignment["hit_len"]) / 2
        # Change the values of the regions used for plotting.
        regions.loc[regions["alignment_number"] == i, "query_from_plot_center"] = (
            regions.loc[regions["alignment_number"] == i, "query_from"] + shift_q
        )
        regions.loc[regions["alignment_number"] == i, "query_to_plot_center"] = (
            regions.loc[regions["alignment_number"] == i, "query_to"] + shift_q
        )
        regions.loc[regions["alignment_number"] == i, "hit_from_plot_center"] = (
            regions.loc[regions["alignment_number"] == i, "hit_from"] + shift_h
        )
        regions.loc[regions["alignment_number"] == i, "hit_to_plot_center"] = (
            regions.loc[regions["alignment_number"] == i, "hit_to"] + shift_h
        )

    return alignments, regions


def include_coordinates_to_right_align_alignments(
    alignments: DataFrame,
    regions: DataFrame,
    size_longest_sequence: int,
) -> tuple[DataFrame, DataFrame]:
    """
    Right-align alignment regions for plotting relative to the longest sequence.

    This function shifts the `*_plot` coordinates (`query_from_plot`, `query_to_plot`,
    `hit_from_plot`, `hit_to_plot`) of each region so that both query and hit
    alignments appear right-aligned in the plot.

    Parameters
    ----------
    alignments : pandas.DataFrame
        DataFrame containing summary metadata for each alignment. Must include:
            - 'alignment_number', 'query_len', and 'hit_len'.

    regions : pandas.DataFrame
        DataFrame with alignment region metadata. Must include:
            - 'alignment_number', 'query_from_plot', 'query_to_plot',
              'hit_from_plot', and 'hit_to_plot'.

    size_longest_sequence : int
        The length of the longest sequence in the dataset. Used to compute the
        right-shift offset for alignment display.
    """
    # Iterate over alignments to find the shift value
    for i, alignment in alignments.iterrows():
        # Find the amount to add to shift the alignments the the right.
        delta_query = size_longest_sequence - alignment["query_len"]
        delta_hit = size_longest_sequence - alignment["hit_len"]
        # Change the values fo the regions used for plotting.
        regions.loc[regions["alignment_number"] == i, "query_from_plot_right"] = (
            regions.loc[regions["alignment_number"] == i, "query_from"] + delta_query
        )
        regions.loc[regions["alignment_number"] == i, "query_to_plot_right"] = (
            regions.loc[regions["alignment_number"] == i, "query_to"] + delta_query
        )
        regions.loc[regions["alignment_number"] == i, "hit_from_plot_right"] = (
            regions.loc[regions["alignment_number"] == i, "hit_from"] + delta_hit
        )
        regions.loc[regions["alignment_number"] == i, "hit_to_plot_right"] = (
            regions.loc[regions["alignment_number"] == i, "hit_to"] + delta_hit
        )

    return alignments, regions
