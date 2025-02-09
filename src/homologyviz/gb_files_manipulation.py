"""Functions to manipulate GenBank files and BLASTn results for HomologyViz.

License
-------
This file is part of HomologyViz
BSD 3-Clause License
Copyright (c) 2024, Ivan Munoz Gutierrez
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
    """Make fasta files from GenBank files.

    Parameters
    ----------
    gb_files : list[Path]
        Paths' list of GenBank files.
    output_path : Path
        Path to folder that will store the fasta files.

    Returns
    -------
    faa_files : list[Path]
        Paths' list of fasta files names.
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
    """Run blastn locally and create xml result file(s).

    Parameters
    ----------
    faa_files : list[Path]
        Paths' list of fasta files.
    output_path : Path
        Path to save files produced by blastn

    Returns
    -------
    results : list[Path]
        Paths' list of xml files with blastn results.
    """
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
        print(f"BLASTing {faa_files[i]} (query) and {faa_files[i+1]} (subject)\n")
        print(std)
    return results


def blastn_command_line(query: Path, subject: Path, out: Path, outfmt: int = 5) -> str:
    """Run BLASTn locally to compare two DNA sequences.

    Notes
    -----
    query and subject must be FASTA files.
    output format must be xml for HomologyViz to work, therefore, outfmt input is 5.
    """
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
    """Compile GenBank files metadata into Pandas DataFrames.

    This function is designed to parse and store metadata and features of GenBank files
    for downstream applications such as visualization or analysis. It extracts attributes
    like the sequence name, coding sequences, and sequence length, providing a structured
    representation of the data using Pandas DataFrames.

    Returns
    -------
    gb_df : pandas.DataFrame
        DataFrame storing GenBank record data such as file path, file name, and sequence
        length. The DataFrame includes a file number and an acession number colum to make
        a relational database with the cds_df
    cds_df : pandas.DataFrame
        DataFrame storing the coding sequences information from each GenBank file. The
        DataFrame includes a file number and an acession number colum to make a relational
        database with the gb_df. This dataframe compiles gene name, product name, start of
        the gene, end of the gene, strand, color, start of the gene for ploting and end of
        the gene for plotting.
    """
    # headers gb_files_df
    headers_gb_files_df = [
        "file_number",
        "file_path",
        "file_name",
        "record_name",
        "accession",
        "length",
        "sequence_start",
        "sequence_end",
    ]
    # Initiate dictionary to stroe data
    gb_files_data = dict(
        file_number=[],
        file_path=[],
        file_name=[],
        record_name=[],
        accession=[],
        length=[],
        sequence_start=[],
        sequence_end=[],
    )
    # Initiate a list of cds DataFrames
    cds_dataframes = []
    # Iterate over GenBank files
    for i, gb_file in enumerate(gb_files):
        # fill data related to the file
        gb_files_data["file_number"].append(i)
        gb_files_data["file_path"].append(gb_file)
        gb_files_data["file_name"].append(gb_file.stem)
        # Read the file into a temporary variable
        record = SeqIO.read(gb_file, "genbank")
        # fill data related to the GenBank record
        gb_files_data["record_name"].append(record.name)
        gb_files_data["accession"].append(record.id)
        seq_length = len(record)
        gb_files_data["length"].append(float(seq_length))
        gb_files_data["sequence_start"].append(0.0)
        gb_files_data["sequence_end"].append(float(seq_length))

        # Get a DataFrame from the cds
        cds_dataframes.append(
            parse_genbank_cds_to_df(record=record, file_number=i, accession=record.id)
        )
    # Create the GenBank files DataFrame
    gb_df = DataFrame(gb_files_data, columns=headers_gb_files_df)
    # Concatenate the cds_dataframes list into a single DataFrame
    cds_df = pd.concat(cds_dataframes, ignore_index=True)

    return gb_df, cds_df


def parse_genbank_cds_to_df(
    record: SeqRecord, file_number: int, accession: str
) -> DataFrame:
    """Parse a Bio.SeqRecord from a GenBank file and make a Pandas DataFrame.

    Parameters
    ----------
    record : Bio.SeqRecord object

    Returns
    -------
    pandas.DataFrame
        DataFrame holding CDSs' information.
    """
    # DataFrame headers
    headers = [
        "file_number",
        "cds_number",
        "accession",
        "gene",
        "product",
        "start",
        "end",
        "strand",
        "color",
        "start_plot",
        "end_plot",
    ]
    # Initiate dictionary to store data
    data = dict(
        file_number=[],
        cds_number=[],
        accession=[],
        gene=[],
        product=[],
        start=[],
        end=[],
        strand=[],
        color=[],
        start_plot=[],
        end_plot=[],
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
        if gene := feature.qualifiers.get("gene", None):
            data["gene"].append(gene[0])
        else:
            data["gene"].append(None)
        if product := feature.qualifiers.get("product", None):
            data["product"].append(product[0])
        else:
            data["product"].append(None)
        if feature.qualifiers.get("Color", None):
            data["color"].append(feature.qualifiers["Color"][0])
        else:
            data["color"].append("#ffff00")  # Make yellow default color
        # Some CDS are composed of more than one parts, like introns, or,
        # in the case of some bacteria, some genes have frameshifts as a
        # regulatory function (some transposase genes have frameshifts as
        # a regulatory function).
        for part in feature.location.parts:
            strand = part._strand
            data["strand"].append(strand)
            if strand == -1:
                data["start"].append(float(part._end))
                data["start_plot"].append(float(part._end))
                data["end"].append(float(part._start + 1))
                data["end_plot"].append(float(part._start + 1))
            else:
                data["start"].append(float(part._start + 1))
                data["start_plot"].append(float(part._start + 1))
                data["end"].append(float(part._end))
                data["end_plot"].append(float(part._end))
    # Create DataFrame
    df = DataFrame(data, columns=headers)
    return df


def blast_alignments_to_dataframe(
    xml_alignment_result: list[Path],
) -> list[DataFrame, DataFrame]:
    """Store a BLASTn alignments results into Pandas DataFrames.

    Return
    ------
    alignments_df : pandas.DataFrame
        DataFrame storing BLAST metadata results such as alignments number, query name,
        hit name, query length, and hit length.
    regions_df : pandas.DataFrame
        DataFrame storing the regions that math between two sequences during BLASTing. The
        metadata includes an alignment number to make a relational database with the
        alignments_df. For more details of metadata stored in regions_df check the
        `parse_blast_record` function in the gb_files_manipulation module.
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
                parse_blast_record(blast_record=blast_record, alignment_number=i)
            )
    # Create DataFrame
    alignments_df = DataFrame(data, columns=headers)
    regions_df = pd.concat(regions, ignore_index=True)
    return alignments_df, regions_df


def parse_blast_record(blast_record: Record, alignment_number: int) -> DataFrame:
    """Parse blast records to store metadata related the regions that match during
    BLASTing two sequences.

    Check the metadata in the `headers` variable.
    """
    headers = [
        "alignment_number",
        "query_from",
        "query_to",
        "query_from_plot",
        "query_to_plot",
        "hit_from",
        "hit_to",
        "hit_from_plot",
        "hit_to_plot",
        "identity",
        "positive",
        "align_len",
        "homology",
    ]
    data = dict(
        alignment_number=[],
        query_from=[],
        query_to=[],
        query_from_plot=[],
        query_to_plot=[],
        hit_from=[],
        hit_to=[],
        hit_from_plot=[],
        hit_to_plot=[],
        identity=[],
        positive=[],
        align_len=[],
        homology=[],
    )
    for region in blast_record.alignments[0].hsps:
        data["alignment_number"].append(alignment_number)
        data["query_from"].append(float(region.query_start))
        data["query_to"].append(float(region.query_end))
        data["query_from_plot"].append(float(region.query_start))
        data["query_to_plot"].append(float(region.query_end))
        data["hit_from"].append(float(region.sbjct_start))
        data["hit_to"].append(float(region.sbjct_end))
        data["hit_from_plot"].append(float(region.sbjct_start))
        data["hit_to_plot"].append(float(region.sbjct_end))
        data["identity"].append(int(region.identities))
        data["positive"].append(int(region.positives))
        data["align_len"].append(int(region.align_length))
        homology = int(region.identities) / int(region.align_length)
        data["homology"].append(homology)
    regions_df = pd.DataFrame(data, columns=headers)
    return regions_df


def get_longest_sequence_dataframe(gb_records: DataFrame) -> int:
    """Finde the longest sequence in the GenBank metadata DataFrame."""
    longest = gb_records["length"].max()
    return longest


def find_lowest_and_highest_homology_dataframe(regions_df: DataFrame) -> tuple:
    lowest = regions_df["homology"].min()
    highest = regions_df["homology"].max()
    return lowest, highest


def adjust_positions_sequences_df_left(gb_records: DataFrame, cds: DataFrame) -> None:
    """Adjust the position of the sequences to the left including the CDSs.

    The start_plot and end_plot columns of the cds DataFrame are used for plotting.
    Therefore, their values are modified.
    """
    # Reset the values of gb_records and cds to the left
    gb_records["sequence_start"] = 0.0
    gb_records["sequence_end"] = gb_records["length"]
    cds["start_plot"] = cds["start"]
    cds["end_plot"] = cds["end"]


def adjust_positions_sequences_df_center(
    gb_records: DataFrame, cds: DataFrame, size_longest_sequence: int
) -> None:
    """Adjust position of sequences to the center including the CDSs.

    The start_plot and end_plot columns of the cds DataFrame are used for plotting.
    Therefore, this function modifies their values.
    """
    # Check if sequences are at the left. If not, reset the values to the left
    if not check_if_sequences_are_at_left(cds):
        adjust_positions_sequences_df_left(gb_records, cds)
    # Iterate over gb_records rows to find the shift value
    for i, row in gb_records.iterrows():
        # Get value to shift sequences to the center
        shift = (size_longest_sequence - row["length"]) / 2
        # Change the values of the sequence_start and sequence_end of gb_records
        gb_records.loc[i, "sequence_start"] = row["sequence_start"] + shift
        gb_records.loc[i, "sequence_end"] = row["sequence_end"] + shift
        # Change the values of start_plot and end_plot of the cds DataFrame
        cds.loc[cds["file_number"] == i, "start_plot"] += shift
        cds.loc[cds["file_number"] == i, "end_plot"] += shift


def adjust_positions_sequences_df_right(
    gb_records: DataFrame, cds: DataFrame, size_longest_sequence: int
) -> None:
    """Adjust position of sequences to the right including the CDSs.

    The start_plot and end_plot columns of the cds DataFrame are used for plotting.
    Therefore, their values are modified.
    """
    # Check if sequences are at the left. If not, reset the values to the left
    if not check_if_sequences_are_at_left(cds):
        adjust_positions_sequences_df_left(gb_records, cds)
    # Iterate over gb_records rows to find the shift value
    for i, row in gb_records.iterrows():
        # Get value to shift sequences to the center
        shift = size_longest_sequence - row["length"]
        # Change the values of the sequence_start and sequence_end of gb_records
        gb_records.loc[i, "sequence_start"] += shift
        gb_records.loc[i, "sequence_end"] += shift
        # Change the values of start_plot and end_plot of the cds DataFrame
        cds.loc[cds["file_number"] == i, "start_plot"] += shift
        cds.loc[cds["file_number"] == i, "end_plot"] += shift


def adjust_positions_alignments_df_left(regions: DataFrame) -> None:
    """Adjust position of alignmets to the left.

    The query_from_plot, query_to_plot, hit_from_plot, and hit_to_plot columns of the
    alignments DataFramme are used for plotting. Therefore, this function modifies their
    values.
    """
    # Reset values
    regions["query_from_plot"] = regions["query_from"]
    regions["query_to_plot"] = regions["query_to"]
    regions["hit_from_plot"] = regions["hit_from"]
    regions["hit_to_plot"] = regions["hit_to"]


def adjust_positions_alignments_df_center(
    alignments: DataFrame, regions: DataFrame, size_longest_sequence: int
) -> None:
    """Adjust position of alignmets to the center.

    The query_from_plot, query_to_plot, hit_from_plot, and hit_to_plot columns of the
    alignments DataFramme are used for plotting. Therefore, this function modifies their
    values.
    """
    # Check if alignments are at the left. If not, reset the values to the left
    if not check_if_alignments_are_at_left(regions):
        adjust_positions_alignments_df_left(regions)
    # Iterate over alignments to find the shift value
    for i, alignment in alignments.iterrows():
        # Find the amount to add to shift the alignments the the center.
        shift_q = (size_longest_sequence - alignment["query_len"]) / 2
        shift_h = (size_longest_sequence - alignment["hit_len"]) / 2
        # Change the values of the regions used for plotting.
        regions.loc[regions["alignment_number"] == i, "query_from_plot"] += shift_q
        regions.loc[regions["alignment_number"] == i, "query_to_plot"] += shift_q
        regions.loc[regions["alignment_number"] == i, "hit_from_plot"] += shift_h
        regions.loc[regions["alignment_number"] == i, "hit_to_plot"] += shift_h


def adjust_positions_alignments_df_right(
    alignments: DataFrame, regions: DataFrame, size_longest_sequence: int
) -> None:
    """Adjust position of the alignments to the right.

    The query_from_plot, query_to_plot, hit_from_plot, and hit_to_plot columns of the
    alignments DataFramme are used for plotting. Therefore, this function modifies their
    values.
    """
    # Check if alignments are at the left. If not, reset the values to the left
    if not check_if_alignments_are_at_left(regions):
        adjust_positions_alignments_df_left(regions)
    # Iterate over alignments to find the shift value
    for i, alignment in alignments.iterrows():
        # Find the amount to add to shift the alignments the the right.
        delta_query = size_longest_sequence - alignment["query_len"]
        delta_hit = size_longest_sequence - alignment["hit_len"]
        # Change the values fo the regions used for plotting.
        regions.loc[regions["alignment_number"] == i, "query_from_plot"] += delta_query
        regions.loc[regions["alignment_number"] == i, "query_to_plot"] += delta_query
        regions.loc[regions["alignment_number"] == i, "hit_from_plot"] += delta_hit
        regions.loc[regions["alignment_number"] == i, "hit_to_plot"] += delta_hit


def adjust_positions_sequences_and_alignments_df_for_plotting(
    gb_records: DataFrame,
    cds: DataFrame,
    alignments: DataFrame,
    regions: DataFrame,
    size_longest_sequence: None | int = None,
    position: str = "left",
) -> None:
    if position == "left":
        adjust_positions_sequences_df_left(gb_records=gb_records, cds=cds)
        adjust_positions_alignments_df_left(regions=regions)
    if position == "center":
        adjust_positions_sequences_df_center(
            gb_records=gb_records,
            cds=cds,
            size_longest_sequence=size_longest_sequence,
        )
        adjust_positions_alignments_df_center(
            alignments=alignments,
            regions=regions,
            size_longest_sequence=size_longest_sequence,
        )
    if position == "right":
        adjust_positions_sequences_df_right(
            gb_records=gb_records,
            cds=cds,
            size_longest_sequence=size_longest_sequence,
        )
        adjust_positions_alignments_df_right(
            alignments=alignments,
            regions=regions,
            size_longest_sequence=size_longest_sequence,
        )


def check_if_alignments_are_at_left(regions: DataFrame) -> bool:
    left = regions["query_from_plot"].equals(regions["query_from"])
    return left


def check_if_sequences_are_at_left(cds: DataFrame):
    left = cds["start_plot"].equals(cds["start"])
    return left


if __name__ == "__main__":
    # # test
    # xml1 = Path(
    #     "/Users/msp/Documents/Coding/python_projects/HomologyViz/data/SW4848_paper/result0.xml"
    # )
    # xml2 = Path(
    #     "/Users/msp/Documents/Coding/python_projects/HomologyViz/data/SW4848_paper/result1.xml"
    # )
    # alignments_df, regions_df = blast_alignments_to_dataframe([xml1, xml2])

    # print(alignments_df)
    # print(regions_df["homology"].min())

    gb1 = Path(
        "/Users/msp/Documents/Coding/python_projects/HomologyViz/data/SW4848_paper/Tn21.gb"
    )
    gb_df, cds_df = genbank_files_metadata_to_dataframes([gb1])
    print(cds_df)
    # cds_df_groups = cds_df.groupby(["file_number"])
    # print(f"the length of the groups is: {len(cds_df_groups)}")
    # for file_number, group in cds_df_groups:
    #     print(file_number)
    #     print(group)
