import base64
from io import BytesIO

import dash
from dash import Input, Output, State
from plotly.graph_objects import Figure


def create_download_data(
    content: bytes,
    filename: str,
    mime_type: str,
) -> dict:
    """
    Encode file content as base64 and build Dash download data.

    Parameters
    ----------
    content : bytes
        File content to encode.
    filename : str
        Name of the downloaded file.
    mime_type : str
        MIME type of the downloaded file.

    Returns
    -------
    dict
        Download data compatible with the Dash Download component.
    """
    encoded = base64.b64encode(content).decode()

    return {
        "base64": True,
        "content": encoded,
        "filename": filename,
        "type": mime_type,
    }


def register_download_callbacks(app: dash.Dash) -> None:
    """Register callbacks related to figure downloads."""

    @app.callback(
        Output("download-plot-component", "data"),
        Input("download-plot-button", "n_clicks"),
        [
            State("plot", "figure"),
            State("figure-format", "value"),
            State("figure-scale", "value"),
            State("figure-width", "value"),
            State("figure-height", "value"),
        ],
        prevent_initial_call=True,
    )
    def download_plot(
        n_clicks: int | None,
        figure: dict,
        figure_format: str,
        scale: int,
        width: int,
        height: int,
    ) -> dict:
        """
        Generate downloadable verison of the current plot.

        The figure is exported as HTML or as a static image using the format and
        dimensions already selected by the user.

        Parameters
        ----------
        n_clicks : int or None
            Number of times the "Download" button has been clicked.
        figure : dict
            The current Plotly figure as a dictionary (from `dcc.Graph`).
        figure_format : str
            The desired output format ("html", "png", "jpeg", "svg", etc.).
        scale : int
            Scaling factor for image resolution (used for static exports).
        width : int
            Width of the exported figure in pixels.
        height : int
            Height of the exported figure in pixels.

        Returns
        -------
        dict
            Download data containing the encoded figure, filename, MIME type, and base64
            flag.
        """
        # Convert figure dictionary into a Figure object
        fig = Figure(
            data=figure["data"],
            layout=figure["layout"],
        )

        if figure_format == "html":
            html_content = fig.to_html(
                full_html=True,
                include_plotlyjs="cdn",
            )

            return create_download_data(
                content=html_content.encode(),
                filename="plot.html",
                mime_type="text/html",
            )

        buffer = BytesIO()

        fig.write_image(
            buffer,
            format=figure_format,
            width=width,
            height=height,
            scale=scale,
            engine="kaleido",
        )
        return create_download_data(
            content=buffer.getvalue(),
            filename=f"plot.{figure_format}",
            mime_type=figure_format,
        )
