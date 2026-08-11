"""
Register Dash callback functions for the HomologyViz graphical interface.

This module wires together the core interactive components of the app, including:
- File upload and deletion for GenBank files
- Plot generation using BLASTn alignments
- UI controls for adjusting annotations, homology colors, and visibility
- Custom color selection and trace selection logic
- Application reset and download features
- Heartbeat monitoring to shut down the app when the browser tab is closed

All callbacks are registered through the `register_callbacks(app)` function.

Notes
-----
- This file is part of HomologyViz
- BSD 3-Clause License
- Copyright (c) 2024, Iván Muñoz Gutiérrez
"""

import tempfile
import atexit
from pathlib import Path

import dash

from homologyviz.parameters import PlotParameters
from homologyviz.callbacks.electrocardiograph import HeartBeatsParameters
from homologyviz.callbacks.ui import register_ui_callbacks
from homologyviz.callbacks.files import register_file_callbacks
from homologyviz.callbacks.download import register_download_callbacks
from homologyviz.callbacks.annotations import register_annotation_callbacks
from homologyviz.callbacks.plot import register_plot_callbacks
from homologyviz.callbacks.lifecycle import register_server_lifecycle

# import logging
# logging.basicConfig(level=logging.INFO)
# logger = logging.getLogger(__name__)


def register_callbacks(app: dash.Dash) -> dash.Dash:
    """
    Register all Dash callbacks for the app, including plotting logic, UI interactins,
    and server shutdown monitoring.

    This function sets up the full interactivity of the Dash app, including:
        - Handling file uploads and deletions.
        - Executing BLASTn alignments and plotting homology regions.
        - Updating annotations, colors, layout, and display options.
        - Managing UI elements like buttons, skeleton loaders, and input states.
        - Generating downloadable figures in various formats.
        - Monitoring heartbeat pings from the frontend to detect tab closure and
          gracefully shut down the app server when inactive.

    Parameters
    ----------
    app : dash.Dash
        The Dash app instance to which all callback functions and server routes will
        be attached.

    Returns
    -------
    dash.Dash
        The same Dash app instance, now with all callbacks registered.
    """
    # Create the tmp directory and ensure it's deleted when the app stops
    tmp_directory = tempfile.TemporaryDirectory()
    atexit.register(tmp_directory.cleanup)
    tmp_directory_path = Path(tmp_directory.name)

    # Monitor the Dash app tab status
    heartbeat_parameters = HeartBeatsParameters()

    # Class to store alignments data
    dash_parameters = PlotParameters()

    # Registers
    register_file_callbacks(
        app=app,
        tmp_directory=tmp_directory_path,
    )

    register_ui_callbacks(
        app=app,
    )

    register_download_callbacks(
        app=app,
    )

    register_annotation_callbacks(
        app=app,
        dash_parameters=dash_parameters,
    )

    register_plot_callbacks(
        app=app,
        dash_parameters=dash_parameters,
        tmp_directory=tmp_directory_path,
    )

    register_server_lifecycle(
        server=app.server,
        heartbeat_parameters=heartbeat_parameters,
    )

    return app
