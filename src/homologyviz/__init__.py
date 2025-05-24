"""
HomologyViz: Visualize BLASTn alignments between annotated DNA sequences.

This package powers the HomologyViz Dash application, designed for bioinformatics
researchers to create graphical representations of pairwise BLASTn alignments. It supports
GenBank files, gene annotations, customizable color scales, and interactive plot editing.

Modules
-------
- app : Initialize and run the HomologyViz graphical user interface
- arrow : Class for plotting arrows to represent genes
- callbacks : Register Dash callback functions for the HomologyViz graphical interface
- cli : Command-line interface (CLI) utilities for HomologyViz
- gb_files_manipulation : Utilities for processing GenBank files and BLASTn results
- layout : Define the layout for the HomologyViz graphical user interface (GUI)
- logger : Provide a reusable function for setting up and retrieving a logger instance
- miscellaneous : General-purpose utility functions used throughout the app
- parameters : Class for storing plot configuration and metadata during Dash callbacks
- plotter : Functions to render DNA sequences, genes, and homology regions
- rectangle_bezier : Generate coordinates for drawing rectangles with curve sides

Example
-------
To launch the app from the command line:

    $ homologyviz

License
-------
This file is part of HomologyViz
BSD 3-Clause License
Copyright (c) 2024, Iván Muñoz Gutiérrez
"""
