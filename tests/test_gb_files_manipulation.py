import pytest
from pathlib import Path

import pandas as pd
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import CompoundLocation
from Bio.Blast import NCBIXML

from homologyviz.gb_files_manipulation import (
    get_blast_metadata,
    parse_genbank_cds_to_df,
    parse_blast_record,
    include_coordinates_to_center_align_sequences,
    include_coordinates_to_right_align_sequences,
    include_coordinates_to_center_align_alignments,
    include_coordinates_to_right_align_alignments,
)

TEST_DATA = Path(__file__).parent / "data"


def test_parse_genbank_cds_to_df_forward_strand():
    record = SeqRecord(Seq("A" * 1000))

    feature = SeqFeature(
        location=SimpleLocation(99, 200, strand=1),
        type="CDS",
        qualifiers={
            "gene": ["test_gene"],
            "product": ["test protein"],
        },
    )

    record.features.append(feature)

    result = parse_genbank_cds_to_df(
        record=record,
        file_number=1,
        accession="TEST001",
    )

    assert len(result) == 1

    row = result.iloc[0]

    assert row["file_number"] == 1
    assert row["cds_number"] == 0
    assert row["accession"] == "TEST001"
    assert row["gene"] == "test_gene"
    assert row["product"] == "test protein"
    assert row["custom_name"] == ""
    assert row["start"] == 100.0
    assert row["end"] == 200.0
    assert row["strand"] == 1
    assert row["color"] == "#ffff00"


def test_parse_genbank_cds_to_df_reverse_strand():
    record = SeqRecord(Seq("A" * 1000))

    feature = SeqFeature(
        location=SimpleLocation(99, 200, strand=-1),
        type="CDS",
        qualifiers={
            "gene": ["reverse_gene"],
            "product": ["reverse protein"],
        },
    )

    record.features.append(feature)

    result = parse_genbank_cds_to_df(
        record=record,
        file_number=1,
        accession="TEST001",
    )

    assert len(result) == 1

    row = result.iloc[0]

    assert row["start"] == 200.0
    assert row["end"] == 100.0
    assert row["strand"] == -1


def test_parse_genbank_cds_to_df_compound_location():
    record = SeqRecord(Seq("A" * 40000))

    location = CompoundLocation(
        [
            SimpleLocation(33669, 33931, strand=1),
            SimpleLocation(33930, 34832, strand=1),
        ],
        operator="join",
    )

    feature = SeqFeature(
        location=location,
        type="CDS",
        qualifiers={
            "gene": ["insD"],
            "product": ["IS3 family transposase"],
        },
    )

    record.features.append(feature)

    result = parse_genbank_cds_to_df(
        record=record,
        file_number=1,
        accession="TEST001",
    )

    assert len(result) == 1

    row = result.iloc[0]

    assert row["start"] == 33670.0
    assert row["end"] == 34832.0
    assert row["strand"] == 1
    assert row["gene"] == "insD"


def test_parse_genbank_cds_to_df_missing_qualifiers():
    record = SeqRecord(Seq("A" * 1000))

    feature = SeqFeature(
        location=SimpleLocation(99, 200, strand=1),
        type="CDS",
        qualifiers={},
    )

    record.features.append(feature)

    result = parse_genbank_cds_to_df(
        record=record,
        file_number=1,
        accession="TEST001",
    )

    row = result.iloc[0]

    assert row["gene"] is None
    assert row["product"] is None
    assert row["custom_name"] == ""
    assert row["color"] == "#ffff00"


def test_parse_genbank_cds_to_df_color_qualifier():
    record = SeqRecord(Seq("A" * 1000))

    feature = SeqFeature(
        location=SimpleLocation(99, 200, strand=1),
        type="CDS",
        qualifiers={
            "gene": ["test_gene"],
            "Color": ["#ff0000"],
        },
    )

    record.features.append(feature)

    result = parse_genbank_cds_to_df(
        record=record,
        file_number=1,
        accession="TEST001",
    )

    row = result.iloc[0]

    assert row["color"] == "#ff0000"


@pytest.mark.parametrize("color_key", ["Color", "color", "COLOR"])
def test_parse_genbank_cds_to_df_color_qualifier_case_insensitive(color_key):
    record = SeqRecord(Seq("A" * 1000))

    feature = SeqFeature(
        location=SimpleLocation(99, 200, strand=1),
        type="CDS",
        qualifiers={
            color_key: ["#ff0000"],
        },
    )

    record.features.append(feature)

    result = parse_genbank_cds_to_df(
        record=record,
        file_number=1,
        accession="TEST001",
    )

    assert result.iloc[0]["color"] == "#ff0000"


@pytest.mark.parametrize(
    "invalid_color",
    [
        "not-a-color",
        "#gg0000",
        "#12345",
        "",
    ],
)
def test_parse_genbank_cds_to_df_invalid_color_uses_default(invalid_color):
    record = SeqRecord(Seq("A" * 1000))

    feature = SeqFeature(
        location=SimpleLocation(99, 200, strand=1),
        type="CDS",
        qualifiers={
            "Color": [invalid_color],
        },
    )

    record.features.append(feature)

    result = parse_genbank_cds_to_df(
        record=record,
        file_number=1,
        accession="TEST001",
    )

    assert result.iloc[0]["color"] == "#ffff00"


def test_parse_genbank_cds_to_df_named_color():
    record = SeqRecord(Seq("A" * 1000))

    feature = SeqFeature(
        location=SimpleLocation(99, 200, strand=1),
        type="CDS",
        qualifiers={
            "Color": ["red"],
        },
    )

    record.features.append(feature)

    result = parse_genbank_cds_to_df(
        record=record,
        file_number=1,
        accession="TEST001",
    )

    assert result.iloc[0]["color"] == "red"


def test_get_blast_metadata():
    xml_file = TEST_DATA / "blast0_test.xml"

    alignments_df, regions_df = get_blast_metadata(
        [xml_file],
        size_longest_sequence=22600,
    )

    assert len(alignments_df) == 1
    assert len(regions_df) == 5

    alignment = alignments_df.iloc[0]

    assert alignment["alignment_number"] == 0
    assert alignment["query_name"] == "CP107285.1 Escherichia coli SW4848"
    assert alignment["query_len"] == 22600
    assert alignment["hit_len"] == 11139

    region = regions_df.iloc[0]

    assert region["alignment_number"] == 0
    assert region["query_from"] == 13358
    assert region["query_to"] == 22585
    assert region["hit_from"] == 1900
    assert region["hit_to"] == 11139
    assert region["identity"] == 9210
    assert region["align_len"] == 9243
    assert region["homology"] == pytest.approx(9210 / 9243)


def test_parse_blast_record():
    xml_file = TEST_DATA / "blast0_test.xml"

    with open(xml_file) as handle:
        blast_record = NCBIXML.read(handle)

    result = parse_blast_record(
        blast_record=blast_record,
        alignment_number=3,
    )

    assert len(result) == 5
    assert (result["alignment_number"] == 3).all()

    region = result.iloc[0]

    assert region["query_from_plot_left"] == 13358
    assert region["query_to_plot_left"] == 22585
    assert region["hit_from_plot_left"] == 1900
    assert region["hit_to_plot_left"] == 11139

    assert pd.isna(region["query_from_plot_center"])
    assert pd.isna(region["query_to_plot_center"])
    assert pd.isna(region["query_from_plot_right"])
    assert pd.isna(region["query_to_plot_right"])

    assert pd.isna(region["hit_from_plot_center"])
    assert pd.isna(region["hit_to_plot_center"])
    assert pd.isna(region["hit_from_plot_right"])
    assert pd.isna(region["hit_to_plot_right"])


def test_include_coordinates_to_center_align_sequences():
    gb_records = pd.DataFrame(
        {
            "file_number": [0, 1],
            "length": [1000, 600],
            "sequence_start": [1, 1],
            "sequence_end": [1000, 600],
        }
    )

    cds = pd.DataFrame(
        {
            "file_number": [0, 1],
            "start": [100, 100],
            "end": [200, 200],
        }
    )

    gb_records, cds = include_coordinates_to_center_align_sequences(
        gb_records=gb_records,
        cds=cds,
    )

    assert gb_records.loc[0, "sequence_start_plot_center"] == 1
    assert gb_records.loc[0, "sequence_end_plot_center"] == 1000

    assert gb_records.loc[1, "sequence_start_plot_center"] == 201
    assert gb_records.loc[1, "sequence_end_plot_center"] == 800

    assert cds.loc[0, "start_plot_center"] == 100
    assert cds.loc[0, "end_plot_center"] == 200

    assert cds.loc[1, "start_plot_center"] == 300
    assert cds.loc[1, "end_plot_center"] == 400


def test_center_align_sequences_uses_file_number_for_cds():
    gb_records = pd.DataFrame(
        {
            "file_number": [5, 8],
            "length": [1000, 600],
            "sequence_start": [1, 1],
            "sequence_end": [1000, 600],
        }
    )

    cds = pd.DataFrame(
        {
            "file_number": [5, 8],
            "start": [100, 100],
            "end": [200, 200],
        }
    )

    _, cds = include_coordinates_to_center_align_sequences(
        gb_records=gb_records,
        cds=cds,
    )

    assert cds.loc[0, "start_plot_center"] == 100
    assert cds.loc[0, "end_plot_center"] == 200

    assert cds.loc[1, "start_plot_center"] == 300
    assert cds.loc[1, "end_plot_center"] == 400


def test_include_coordinates_to_right_align_sequences():
    gb_records = pd.DataFrame(
        {
            "file_number": [5, 8],
            "length": [1000, 600],
            "sequence_start": [1, 1],
            "sequence_end": [1000, 600],
        }
    )

    cds = pd.DataFrame(
        {
            "file_number": [5, 8],
            "start": [100, 100],
            "end": [200, 200],
        }
    )

    gb_records, cds = include_coordinates_to_right_align_sequences(
        gb_records=gb_records,
        cds=cds,
    )

    assert gb_records.loc[0, "sequence_start_plot_right"] == 1
    assert gb_records.loc[0, "sequence_end_plot_right"] == 1000

    assert gb_records.loc[1, "sequence_start_plot_right"] == 401
    assert gb_records.loc[1, "sequence_end_plot_right"] == 1000

    assert cds.loc[0, "start_plot_right"] == 100
    assert cds.loc[0, "end_plot_right"] == 200

    assert cds.loc[1, "start_plot_right"] == 500
    assert cds.loc[1, "end_plot_right"] == 600


def test_include_coordinates_to_center_align_alignments():
    alignments = pd.DataFrame(
        {
            "alignment_number": [5],
            "query_len": [1000],
            "hit_len": [600],
        }
    )

    regions = pd.DataFrame(
        {
            "alignment_number": [5],
            "query_from": [100],
            "query_to": [200],
            "hit_from": [100],
            "hit_to": [200],
        }
    )

    _, regions = include_coordinates_to_center_align_alignments(
        alignments=alignments,
        regions=regions,
        size_longest_sequence=1000,
    )

    region = regions.iloc[0]

    assert region["query_from_plot_center"] == 100
    assert region["query_to_plot_center"] == 200

    assert region["hit_from_plot_center"] == 300
    assert region["hit_to_plot_center"] == 400


def test_include_coordinates_to_right_align_alignments():
    alignments = pd.DataFrame(
        {
            "alignment_number": [5],
            "query_len": [1000],
            "hit_len": [600],
        }
    )

    regions = pd.DataFrame(
        {
            "alignment_number": [5],
            "query_from": [100],
            "query_to": [200],
            "hit_from": [100],
            "hit_to": [200],
        }
    )

    _, regions = include_coordinates_to_right_align_alignments(
        alignments=alignments,
        regions=regions,
        size_longest_sequence=1000,
    )

    region = regions.iloc[0]

    assert region["query_from_plot_right"] == 100
    assert region["query_to_plot_right"] == 200

    assert region["hit_from_plot_right"] == 500
    assert region["hit_to_plot_right"] == 600
