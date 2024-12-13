"""Make HomologyViz GUI.

License
-------
This file is part of HomologyViz
BSD 3-Clause License
Copyright (c) 2024, Ivan Munoz Gutierrez
"""

import webbrowser

import dash
from dash import _dash_renderer
import dash_bootstrap_components as dbc

from homologyviz.callbacks import register_callbacks
from homologyviz.layout import create_layout

from homologyviz.cli import parse_command_line_input


def main() -> None:
    """Make the app."""
    # Parse command line input
    parse_command_line_input()

    # This variable must be set according to the Dash Mantine Components
    _dash_renderer._set_react_version("18.2.0")

    # Initialize the Dash app with a Bootstrap theme
    app = dash.Dash(__name__, external_stylesheets=[dbc.themes.DARKLY])

    # Create the app layout
    app = create_layout(app)

    # Register callbacks
    app = register_callbacks(app=app)

    # Open the app in the default web browser
    webbrowser.open("http://127.0.0.1:8050")

    # Run the app
    app.run(
        # debug=True,
        # use_reloader=False,
    )
