from pathlib import Path
import base64
import binascii

import dash
from dash import Input, Output, State


def save_uploaded_file(
    file_name: str,
    content: str,
    tmp_directory: Path,
) -> str | None:
    """Decode the content and write it to a temporary file.

    Returns the file path as a string if successful, otherwise returns None.
    """
    try:
        # Ensure content is properly formatted
        if ";base64," not in content:
            raise ValueError("Content is not base64-encoded or improperly formatted.")

        # Decode content
        data = content.split(";base64,")[1]
        decoded_data = base64.b64decode(data)

        # Save uploaded file
        output_path = tmp_directory / file_name
        with open(output_path, "wb") as f:
            f.write(decoded_data)

        # Dash doesn't like Path; hence, we need to cast Path to str.
        return str(output_path)

    except (IndexError, ValueError, binascii.Error) as e:
        print(f"Failed to decode and save uplaoded file: {e}")
        return None


def register_file_callbacks(
    app: dash.Dash,
    tmp_directory: Path,
) -> None:
    @app.callback(
        Output("files-table", "rowData"),
        [
            Input("upload", "filename"),
            Input("upload", "contents"),
            Input("trash-selected-files-button", "n_clicks"),
        ],
        [
            State("files-table", "rowData"),
            State("files-table", "selectedRows"),
        ],
    )
    def update_file_table(
        filenames: list | None,
        contents: list | None,
        n_clicks: int | None,
        current_row_data: list[dict] | None,
        selected_rows: list[dict] | None,
    ) -> list[dict[str, str]]:
        """
        Update the GenBank files table based on uploaded files or deletion actions.

        This callback populates the table with uploaded files by decoding their content
        and saving them temporarily. It also supports removing selected rows when the
        "Trash Selected Files" button is clicked.

        Parameters
        ----------
        filenames : list[str] or None
            List of filenames uploaded via the upload component.
        contents : list[str] or None
            Corresponding list of base64-encoded file contents.
        n_clicks : int or None
            Number of items the delete button has been clicked.
        current_row_data : list[dict] or None
            Current content of the table (`files-table`) as a list of dictionaries.
        selected_rows : list[dict] or None
            Subset of `current_row_data` that the user has selected for deletion.

        Returns
        -------
        list[dict]
            Uploaded list of table rows reflecting uploaded files or deletions.
        """
        ctx_id = dash.ctx.triggered_id

        current_rows = current_row_data or []
        selected = selected_rows or []

        # Update table with uploaded files.
        if ctx_id == "upload" and filenames and contents:
            new_rows = []

            # Simulate saving each file and creating a temporary file path
            for name, content in zip(filenames, contents):
                file_path = save_uploaded_file(
                    file_name=name,
                    content=content,
                    tmp_directory=tmp_directory,
                )

                new_rows.append(
                    {
                        "filename": name,
                        "filepath": file_path,
                    }
                )

            return (current_row_data or []) + new_rows

        # Delete selected rows
        if ctx_id == "trash-selected-files-button":
            return [row for row in current_rows if row not in selected]

        return current_row_data or []

    @app.callback(
        [
            Output("plot-button", "disabled"),
            Output("trash-selected-files-button", "disabled"),
        ],
        Input("files-table", "rowData"),
    )
    def toggle_file_buttons(
        row_data: list[dict[str, str]] | None,
    ) -> tuple[bool, bool]:
        """
        Enable or disable file-related buttons based on the number of uploaded files.

        The Plot button is enabled when at least two files are available. The Trash
        button is enabled when at least one file is available.

        Parameters
        ----------
        row_data : list[dict[str, str]] or None
            Current rows in the uploaded-files table.

        Returns
        -------
        tuple[bool, bool]
            Disabled states for the Plot and Trash buttons, respectively.
        """
        file_count = len(row_data) if row_data else 0

        plot_disabled = file_count < 2
        trash_disabled = file_count < 1

        return plot_disabled, trash_disabled
