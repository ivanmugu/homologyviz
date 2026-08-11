import dash
from dash import Input, Output, State
from plotly.graph_objects import Figure

from homologyviz import plotter as plt


def register_ui_callbacks(app: dash.Dash) -> None:
    @app.callback(
        [
            Output("select-items-button", "variant"),
            Output("select-items-button-store", "data"),
        ],
        Input("select-items-button", "n_clicks"),
        State("select-items-button-store", "data"),
    )
    def toggle_select_items_button(
        n_clicks: int | None,
        is_active: bool,
    ) -> tuple[str, bool]:
        """
        Toggle trace-selection mode and update the button appearance.

        Parameters
        ----------
        n_clicks : int or None
            Number of times the Select Items button has been clicked.
        is_active : bool
            Current selection-mode state stored in Dash.

        Returns
        -------
        tuple[str, bool]
            The button variant and the updated selection-mode state.
        """
        if n_clicks:
            is_active = not is_active

        variant = "filled" if is_active else "outline"

        return variant, is_active

    @app.callback(
        [
            Output("extreme-homologies-button", "variant"),
            Output("extreme-homologies-button", "style"),
            Output("truncate-colorscale-button", "variant"),
            Output("truncate-colorscale-button", "style"),
            Output("is-set-to-extreme-homologies", "data"),
        ],
        [
            Input("extreme-homologies-button", "n_clicks"),
            Input("truncate-colorscale-button", "n_clicks"),
        ],
    )
    def toggle_colorscale_buttons(
        extreme_clicks: int | None,
        truncate_clicks: int | None,
    ) -> tuple[str, dict, str, dict, bool]:
        """
        Toggle between extreme-homology and truncated colorscale modes.

        The selected mode is shown as the filled button and is made non-interactive until
        the other mode is selected.

        Parameters
        ----------
        extreme_clicks : int or None
            Number of times the Extreme Homologies button has been clicked.

        truncate_clicks : int or None
            Number of times the Truncate Colorscale button has been clicked.

        Returns
        -------
        tuple[str, dict, str, dict, bool]
            Variants and styles for both buttons, followed by whether the extreme-homology
            mode is active.
        """
        truncate_mode = (
            "subtle",
            {"width": "280px", "padding": "5px"},
            "filled",
            {"width": "280px", "padding": "5px", "pointer-events": "none"},
            False,
        )
        extreme_mode = (
            "filled",
            {"width": "280px", "padding": "5px", "pointer-events": "none"},
            "subtle",
            {"width": "280px", "padding": "5px"},
            True,
        )

        triggered_id = dash.ctx.triggered_id

        if triggered_id == "extreme-homologies-button":
            return extreme_mode

        return truncate_mode

    @app.callback(
        Output("color-scale-display", "figure"),
        Input("color-scale", "value"),
    )
    def update_color_scale(value: str) -> Figure:
        """
        Update the horizontal color gradient display based on the selected colorscale.

        This callback is triggered when the user selects a new colorscale from the
        dropdown menu in the `Edit` tab. It passes the selected value to the
        `create_color_line` function to generate a smooth gradient for visual feedback.

        Parameters
        ----------
        value : str
            The name of the selected Plotly sequential colorscale (e.g., "Greys", "Blues")

        Returns
        -------
        figure : plotly.graph_objects.Figure
            A Plotly figure displaying a horizontal gradient representing the selected
            colorscale.
        """
        return plt.create_color_line(value.capitalize())

    @app.callback(
        Output("offcanvas-edit-sequence-annotations", "is_open"),
        Input("open-offcanvas-edit-sequence-annotations", "n_clicks"),
        State("offcanvas-edit-sequence-annotations", "is_open"),
    )
    def toggle_sequence_annotations_offcanvas(
        n_clicks: int | None,
        is_open: bool,
    ):
        """
        Toggle the sequence annotation editor when its open button is clicked.

        Parameters
        ----------
        n_clicks : int or None
            Number of times the button for opening the sequence annotation editor
            has been clicked.

        is_open : bool
            Current open state of the sequence annotation offcanvas.

        Returns
        -------
        bool
            The updated open state of the offcanvas.
        """
        if n_clicks:
            return not is_open
        return is_open

    @app.callback(
        Output("offcanvas-edit-gene-annotations", "is_open"),
        Input("open-offcanvas-edit-gene-annotations", "n_clicks"),
        State("offcanvas-edit-gene-annotations", "is_open"),
    )
    def toggle_gene_annotations_offcanvas(
        n_clicks: int | None,
        is_open: bool,
    ):
        """
        Toggle the gene annotation editor when its open button is clicked.

        Parameters
        ----------
        n_clicks : int or None
            Number of times the button for opening the gene annotation editor
            has been clicked.

        is_open : bool
            Current open state of the gene annotation offcanvas.

        Returns
        -------
        bool
            The updated open state of the offcanvas.
        """
        if n_clicks:
            return not is_open
        return is_open

    @app.callback(
        Output("url", "href"),
        Input("reset-button", "n_clicks"),
        prevent_initial_call=True,
    )
    def reset_app(n_clicks: int) -> str:
        """
        Reload the app when the "Reset" button is clicked.

        This callback returns the current URL path ("/"), which triggers a full page
        reload in Dash. It serves as a way to reset the interface and clear any stored
        state.

        Parameters
        ----------
        n_clicks : int
            Number of times the "Reset" button has been clicked.

        Returns
        -------
        str
            The URL path ("/") to trigger a browser reload of the app.
        """
        print("Reseting application...")
        return "/"
