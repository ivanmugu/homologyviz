import dash
from dash import Input, Output

from homologyviz.parameters import PlotParameters


def register_annotation_callbacks(
    app: dash.Dash,
    dash_parameters: PlotParameters,
) -> None:
    """Register callbacks related to sequence and gene annotation tables."""

    # = Update edit sequence annotations table ========================================= #
    @app.callback(
        Output("sequence-table", "rowData"),
        Input("open-offcanvas-edit-sequence-annotations", "n_clicks"),
        prevent_initial_call=True,
    )
    def populate_sequence_annotations_table(
        n_clicks: int | None,
    ) -> list[dict]:
        """
        Populate the sequence annotation table when the annotation editor is opened.
        The table is built from the GenBank sequence information stored in
        ``dash_parameters.gb_df``. It displays the available sequence annotation
        fields, including accession, record name, file name, and custom name.
        This callback only initializes the table contents. User edits made to the
        ``custom_name`` column are read later by the plot callback when the user
        applies the sequence annotation changes.

        Parameters
        ----------
        n_clicks : int or None
            Number of times the button for opening the sequence annotation editor
            has been clicked.

        Returns
        -------
        list[dict]
            Rows used to populate the sequence annotation table. Returns an empty
            list if no GenBank sequence data are available.
        """
        if dash_parameters.gb_df is None:
            return []

        if n_clicks:
            rows = []
            for _, gb_file in dash_parameters.gb_df.iterrows():
                rows.append(
                    {
                        "file_number": (gb_file["file_number"] + 1),
                        "accession": gb_file["accession"],
                        "record_name": gb_file["record_name"],
                        "file_name": gb_file["file_name"],
                        "custom_name": gb_file["custom_name"],
                    }
                )
            return rows

        return []

    # = Update edit gene annotations table ============================================= #
    @app.callback(
        Output("gene-table", "rowData"),
        Input("open-offcanvas-edit-gene-annotations", "n_clicks"),
        prevent_initial_call=True,
    )
    def populate_gene_annotations_table(
        n_clicks: int | None,
    ) -> list[dict]:
        """
        Populate the gene annotation table when the annotation editor is opened.

        The table is built from the CDS information stored in ``dash_parameters.cds_df``.
        It displays gene annotation fields together with the identifiers and genomic
        coordinates needed to associate each table row with its corresponding CDS.

        This callback only initializes the table contents. User edits made to editable
        annotation columns are read later by the plot callback when the user applies the
        gene annotation changes.

        Parameters
        ----------
        n_clicks : int or None
            Number of times the button for opening the gene annotation editor has been
            clicked.

        Returns
        -------
        list[dict]
            Rows used to populate the gene annotation table. Returns an empty list if no
            CDS data are available.
        """
        if dash_parameters.cds_df is None:
            return []

        if n_clicks:
            rows = []
            for _, cds in dash_parameters.cds_df.iterrows():
                rows.append(
                    {
                        "cds_number": (cds["cds_number"] + 1),
                        "gene": cds["gene"],
                        "product": cds["product"],
                        "accession": cds["accession"],
                        "custom_name": cds["custom_name"],
                        "file_number": cds["file_number"],
                        "start": cds["start"],
                        "end": cds["end"],
                        "strand": cds["strand"],
                    }
                )
            return rows

        return []
