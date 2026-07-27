from dash import html, dcc
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc

from homologyviz.layout.main_tab import make_main_tab
from homologyviz.layout.view_tab import make_view_tab
from homologyviz.layout.edit_tab import make_edit_tab
from homologyviz.layout.save_tab import make_save_tab


def make_layout_control_panel() -> html.Div:
    """
    Create the left-side control panel layout for the HomologyViz Dash app.

    This panel includes:
    - The HomologyViz logo at the top.
    - A vertically scrollable `dbc.Tabs` menu for navigation between key UI sections:
        - Main (file input, sequence selection)
        - View (customization of layout and display)
        - Edit (annotation and color editing)
        - Save (export options)

    Returns
    -------
    html.Div
        A styled Dash HTML Div component containing the control panel layout.
    """
    return html.Div(
        children=[
            html.Img(
                src="/assets/logo.png",
                className="mx-auto my-4 d-block text-white fw-bold text-center",
                alt="HomologyViz",
                style={
                    "height": "40px",
                    "fontSize": "24px",
                },
            ),
            html.Div(  # Tabs menu
                dbc.Tabs(
                    [
                        make_main_tab(),
                        make_view_tab(),
                        make_edit_tab(),
                        make_save_tab(),
                    ],
                    id="tabs",
                ),
                className="mt-1",
                style={
                    "height": "85%",
                    "width": "100%",
                    "overflow": "auto",
                },
            ),
        ],
        style={
            "backgroundColor": "#242424",
            "height": "95vh",
            "overflow": "auto",
        },
    )


def make_layout_plot_panel() -> html.Div:
    """
    Create the right-side plot panel layout for the HomologyViz Dash app.

    This panel includes:
    - A `dcc.Graph` component where the main BLASTn alignment figure is rendered.
    - A `dmc.Skeleton` component used as a loading placeholder while the figure is
      updating.

    The layout is styled to occupy nearly full vertical height and has a visible border.

    Returns
    -------
    html.Div
        A styled Dash HTML Div containing the main plot area and loading skeleton.
    """
    return html.Div(
        children=[
            dmc.Skeleton(
                id="plot-skeleton",
                visible=False,
                children=dcc.Graph(
                    id="plot",
                    style={"height": "100%"},
                ),
                height="100%",
            ),
        ],
        style={
            "border": "1px solid black",
            "height": "96vh",
        },
    )
