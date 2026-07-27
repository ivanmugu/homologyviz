"""
Define the layout for the HomologyViz graphical user interface (GUI).

This module builds the entire front-end layout of the HomologyViz Dash application using
Dash Mantine Components (DMC), Dash Bootstrap Components (DBC), and Plotly Graphs. The GUI
includes interactive controls to upload files, edit plots, adjust views, and export
figures. It is structured into multiple tabs—Main, View, Edit, and Save—and integrates
seamlessly with Dash callbacks.

Notes
-----
- This file is part of HomologyViz
- BSD 3-Clause License
- Copyright (c) 2024, Iván Muñoz Gutiérrez
"""

from dash import Dash, dcc
import dash_mantine_components as dmc

from homologyviz.layout.panels import make_layout_control_panel, make_layout_plot_panel


def create_layout(app: Dash) -> Dash:
    """
    Construct the full layout for the HomologyViz Dash app.

    This function defines the GUI structure, including the control panel and plot display.
    It uses Dash Mantine Components for styling and layout organization. The layout is
    composed of two primary columns:

    - Left Column: Control panel with the HomologyViz logo and tabbed interface for
      uploading files, customizing views, editing plots, and saving outputs.
    - Right Column: Main plotting area displaying the generated figure using
      `dcc.Graph`, wrapped in a `dmc.Skeleton` for loading effects.

    Parameters
    ----------
    app : dash.Dash
        The Dash application instance to which the layout will be assigned.

    Returns
    -------
    dash.Dash
        The Dash app with its layout fully configured and assigned.
    """
    # Wrap layout with dmc.MantineProvider
    app.layout = dmc.MantineProvider(
        dmc.Grid(
            children=[
                dcc.Location(id="url", refresh=True),
                dmc.GridCol(
                    make_layout_control_panel(),
                    span="auto",
                    style={"maxWidth": "340px", "minWidth": "200px"},
                ),
                dmc.GridCol(
                    make_layout_plot_panel(),
                    span=9,
                ),
            ],
            align="center",
            justify="flex-start",
            gutter="xs",
            style={"padding": "8px"},
        ),
        forceColorScheme="dark",
    )

    return app
