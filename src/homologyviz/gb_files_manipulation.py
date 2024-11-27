"""Functions and classes to manipulate genbank files for msplotly.

License
-------
This file is part of MSPloter
BSD 3-Clause License
Copyright (c) 2024, Ivan Munoz Gutierrez
"""

from pathlib import Path
import subprocess

from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.SeqRecord import SeqRecord


class GenBankRecord:
    """Store relevant info from a GenBank file.

    Attributes
    ----------
    file_name : str
        file name.
    name : str
        Sequence name shown next to the `LOCUS` tag.
    accession : str
        Sequence accession number with version.
    description : str
        Description shown next to the `DEFINITION` tag.
    length : int
        Sequence length.
    sequence_start : int
        Start coordinate used for plotting. Default value is zero but changes
        if `alignments_position` of `MakeFigure` class is set to center or
        left.
    sequence_end : int
        End coordinate used for plotting. Default value is zero but changes if
        `alignments_position` of `MakeFigure` class is set to center or left.
    cds : list
        List of `CodingSequence` classes with info of `CDS` tags as product,
        start, end, strand, and color.
    num_cds : int
        Number of CDSs
    """

    def __init__(self, file_name):
        """
        file_name : Path oject
            Path to file.
        """
        record = SeqIO.read(file_name, "genbank")
        self.file_name = file_name.stem
        self.name = record.name
        self.accession = record.id
        self.description = record.description
        self.length = len(record)
        self.sequence_start = 0
        self.sequence_end = self.length
        self.cds = self.parse_gb(record)
        self.num_cds = len(self.cds)

    def parse_gb(self, record):
        """Parse gb file and make a list of `CodingSequence` classes.

        Parameters
        ----------
        record : Bio SeqIO.read object.

        Returns
        -------
        coding_sequences : list
            List of `CodingSequence` classes holding CDSs' information.
        """
        coding_sequences = []
        for feature in record.features:
            if feature.type != "CDS":
                continue
            if gene := feature.qualifiers.get("gene", None):
                gene = gene[0]
            if product := feature.qualifiers.get("product", None):
                product = product[0]
            if color := feature.qualifiers.get("Color", None):
                color = feature.qualifiers["Color"][0]
            else:
                color = "#ffff00"  # Make yellow default color
            # Some CDS are composed of more than one parts, like introns, or,
            # in the case of some bacteria, some genes have frameshifts as a
            # regulatory function (some transposase genes have frameshifts as
            # a regulatory function).
            for part in feature.location.parts:
                strand = part._strand
                if strand == -1:
                    start = part._end
                    end = part._start + 1
                else:
                    start = part._start + 1
                    end = part._end
                # Append cds.
                coding_sequences.append(
                    CodingSequence(gene, product, start, end, strand, color)
                )
        return coding_sequences


class CodingSequence:
    """Store Coding Sequence (CDS) information from gb file."""

    def __init__(self, gene, product, start, end, strand, color):
        self.gene = gene
        self.product = product
        self.start = int(start)
        self.end = int(end)
        self.strand = int(strand)
        self.color = color


class BlastnAlignment:
    """Store blastn alignment results.

    Attributes
    ----------
    query_name : str
        Name of query sequence.
    hit_name : str
        Name of subject sequence.
    query_len : int
        Length of query sequence.
    hit_len : int
        Length of subject sequence.
    regions : list
        List of `RegionAlignmentResult` classes with info of aligned region as
        query_from, query_to, hit_from, hit_to, and identity.

    Methods
    -------
    blastn_alignment_to_dict
        Makes a Dict representing `BlastnAlignmet` class.
    """

    def __init__(self, xml_alignment_result):
        with open(xml_alignment_result, "r") as result_handle:
            blast_record = NCBIXML.read(result_handle)
            self.query_name = blast_record.query
            self.hit_name = blast_record.alignments[0].hit_def
            self.query_len = int(blast_record.query_length)
            self.hit_len = int(blast_record.alignments[0].length)
            self.regions = self.parse_blast_regions(blast_record)

    def parse_blast_regions(self, blast_record):
        """Parse blastn aligned regions to store the information.

        Parameters
        ----------
        blast_record : NCBIXML object
            Harbors blastn alignment results in xml format.

        Returns
        -------
        regions : list
            List of `RegionAlignmentResult` classes with info from alignment.
        """
        regions = []
        for region in blast_record.alignments[0].hsps:
            regions.append(
                RegionAlignmentResult(
                    query_from=int(region.query_start),
                    query_to=int(region.query_end),
                    hit_from=int(region.sbjct_start),
                    hit_to=int(region.sbjct_end),
                    identity=int(region.identities),
                    positive=int(region.positives),
                    align_len=int(region.align_length),
                )
            )
        return regions


class RegionAlignmentResult:
    """Save blastn results of region that aligned."""

    def __init__(
        self, query_from, query_to, hit_from, hit_to, identity, positive, align_len
    ):
        self.query_from = query_from
        self.query_to = query_to
        self.hit_from = hit_from
        self.hit_to = hit_to
        self.identity = identity
        self.positive = positive
        self.align_len = align_len
        self.homology = identity / align_len


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


def get_alignment_records(alignment_files: list[Path]) -> list:
    """Parse xml alignment files and make list of `BlastnAlignment` classes."""
    alignments = [BlastnAlignment(alignment) for alignment in alignment_files]
    return alignments


def get_gb_records(gb_files: list[Path]) -> list:
    """Parse gb files and make list of `GenBankRecord` classes."""
    gb_records = [GenBankRecord(gb_file) for gb_file in gb_files]
    return gb_records


if __name__ == "__main__":
    # test
    print("printing from gb_files_manipulation.")
