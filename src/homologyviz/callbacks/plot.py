from pathlib import Path

import dash
from dash import Input, Output, State
from plotly.graph_objects import Figure

from homologyviz.parameters import PlotParameters
from homologyviz.callbacks.clicks import (
    handle_plot_button_click,
    handle_update_view_click,
    handle_update_homologies_click,
    handle_update_title_click,
    handle_change_color_click,
    handle_select_traces_click,
)
from homologyviz.callbacks.updates import (
    update_genes_annotations,
    update_dna_sequence_annotations,
)


def register_plot_callbacks(
    app: dash.Dash,
    dash_parameters: PlotParameters,
    tmp_directory: Path,
) -> None:
    """Register callbacks related to the alignment plot."""

    @app.callback(
        [
            Output("erase-button", "disabled"),
            Output("update-annotations", "disabled"),
            Output("change-gene-color-button", "disabled"),
            Output("change-homology-color-button", "disabled"),
            Output("select-items-button", "disabled"),
            Output("update-title-button", "disabled"),
            Output("offcanvas-update-sequence-annotations-button", "disabled"),
            Output("offcanvas-update-gene-annotations-button", "disabled"),
        ],
        Input("plot", "figure"),
    )
    def toggle_update_buttons(figure: dict) -> tuple[bool, ...]:
        """
        Enable editing controls when a plot is available.

        Parameters
        ----------
        figure : dict
            Current Plotly figure stored in the graph component.

        Returns
        -------
        tuple[bool, ...]
            Disabled states for the eight controls associated with plot editing.
            All controls are enabled when the figure contains plot data and
            disabled otherwise.
        """
        has_plot = bool(figure and figure.get("data"))
        disabled = not has_plot
        return (disabled,) * 8

    @app.callback(
        [
            Output("plot", "figure"),
            Output("plot", "clickData"),
            Output("plot-skeleton", "visible"),
        ],
        [
            Input("plot-button", "n_clicks"),
            Input("erase-button", "n_clicks"),
            Input("plot", "clickData"),
            Input("change-homology-color-button", "n_clicks"),
            Input("change-gene-color-button", "n_clicks"),
            Input("update-annotations", "n_clicks"),
            Input("update-title-button", "n_clicks"),
            Input("offcanvas-update-sequence-annotations-button", "n_clicks"),
            Input("offcanvas-update-gene-annotations-button", "n_clicks"),
        ],
        [
            State("files-table", "virtualRowData"),
            State("tabs", "active_tab"),
            State("plot", "figure"),
            State("color-input", "value"),
            State("select-items-button-store", "data"),
            State("color-scale", "value"),
            State("range-slider", "value"),
            State("align-plot", "value"),
            State("homology-style", "value"),
            State("minimum-homology-length", "value"),
            State("is-set-to-extreme-homologies", "data"),
            State("scale-bar", "value"),
            State("gene-annotation-positions-choice", "value"),
            State("gene-annotation-from-choice", "value"),
            State("gene-table", "virtualRowData"),
            State("title-input", "value"),
            State("sequence-table", "virtualRowData"),
            State("annotation-column-choice", "value"),
        ],
        prevent_initial_call=True,
    )
    def main_plot(
        plot_button_clicks: int | None,
        clear_button_clicks: int | None,
        click_data: dict | None,
        change_homology_color_button_clicks: int | None,
        change_gene_color_button_clicks: int | None,
        update_annotations_clicks: int | None,
        update_title_button_clicks: int | None,
        update_sequence_annotations: int | None,
        update_gene_annotations: int | None,
        virtual: list[dict[str, str]],
        active_tab: str,
        figure_state: dict,
        color_input_state: str,
        select_items_state: bool,
        color_scale_state: str,
        range_slider_state: list[int],
        align_plot_state: str,
        homology_style_state: str,
        minimum_homology_length_state: int,
        is_set_to_extreme_homologies: bool,
        scale_bar_state: str,
        gene_annotation_positions_state: str,
        gene_annotation_from_state: str,
        gene_table_state: list[dict] | None,
        title_input_state: str,
        sequence_table_state: list[dict] | None,
        annotation_column_choice_state: str,
    ) -> tuple[Figure | dict, None, bool]:
        """
        Master callback function for generating and modifying the alignment plot.

        This function coordinates user interactions accross multiple tabs:
        - In the **Main** tab, it triggers the BLASTn alignments and generates the plot.
        - In the **Edit** tab, it enables trace selection and color modifications.
        - In the **View** tab, it updates gene/sequence annotations, scale bar visibility,
          alignment position, and homology filtering.

        The function also resets the plot when the user clicks the "Erase" button.

        Parameters
        ----------
        plot_button_clicks : int | None
            Number of times the "Plot" button has been clicked.
        clear_button_clicks : int | None
            Number of times the "Erase" button has been clicked.
        click_data : dict | None
            Data from clicking a trace on the plot (used for selecting/deselecting
            traces).
        change_homology_color_button_clicks : int | None
            Click count for the button that changes homology colors and legend.
        change_gene_color_button_clicks : int | None
            Click count for the button that updates selected gene trace colors.
        update_annotations_clicks : int | None
            Click count for the button that updates annotations and plot layout.
        update_title_button_clicks : int | None
            Click count for the button that updates the figure's title.
        update_sequence_annotations : int | None
            Click count for the button that updates DNA sequence annotations.
        update_gene_annotations : int | None
            Click count for the button that updates gene annotations.
        virtual : list[dict[str, str]] | None
            Virtual row data from the GenBank file upload table.
        active_tab : str
            The currently active tab in the UI (e.g., "tab-main", "tab-edit", "tab-view").
        figure_state : dict
            The current figure state stored in Dash, used to rebuild or modify the plot.
        color_input_state : str
            Color selected for updating gene trace colors (hex string).
        select_items_state : bool
            Whether the "Select Items" button is active for toggling trace selection.
        color_scale_state : str
            The selected colorscale name used for identity-based homology coloring.
        range_slider_state : list[int, int]
            The selected range of identity percentages used for color scaling.
        align_plot_state : str
            Alignment layout setting (e.g., "left", "center", "right").
        homology_style_state : str
            Whether the style of homology connections are straight or curve (Bezier).
        minimum_homology_length_state : int
            Minimum homology length (in bp) to display in the plot.
        is_set_to_extreme_homologies : bool
            Whether to stretch the color scale to the dataset min/max identity values.
        annotate_genes_state : str
            Whether and where to annotate gene features (e.g., "top", "bottom", "no").
        scale_bar_state : str
            Whether to show the scale bar ("yes" or "no").
        use_genes_info_from_state : str
            Source of gene labels used for annotation (e.g., "CDS product", "CDS gene").
        title_input_state : str
            String holding the figure's title.
        sequence_table_state : list[dict[str, str, str, str, str]] | None
            Table from 'Edit Sequence Annotations' with the following columns: File,
            Accession, Record Name, File name, and Custom name.
        annotation_column_choice_state : str
            Whether and how to annotate DNA sequences (e.g., "accession", "name", "no",
            "custom").

        Returns
        -------
        fig : plotly.graph_objects.Figure
            The updated Plotly figure, either newly created or modified.
        None
            Placeholder to reset `clickData` in Dash (prevents stuck selections).
        bool
            A flag (`False`) to hide the dmc.Skeleton loading component after plot
            rendering.

        Notes
        -----
        - Uses `dash.callback_context` to determine which button triggered the callback.
        - This function is central to all updates affecting the alignment plot.
        """
        # Use context to find the button that triggered the callback.
        button_id = dash.callback_context.triggered_id
        print(f"button_id: {button_id}")

        # = TAB MAIN =================================================================== #
        # Perform Alignments & Plot
        if (button_id == "plot-button") and virtual:
            print("Plot button clicked, starting alignment and plotting...")
            return handle_plot_button_click(
                dash_parameters=dash_parameters,
                virtual=virtual,
                tmp_directory_path=tmp_directory,
                align_plot_state=align_plot_state,
                color_scale_state=color_scale_state,
                range_slider_state=range_slider_state,
                is_set_to_extreme_homologies=is_set_to_extreme_homologies,
                annotation_column_choice_state=annotation_column_choice_state,
                annotate_genes_from_state=gene_annotation_from_state,
                annotate_genes_positions_state=gene_annotation_positions_state,
                homology_style_state=homology_style_state,
                minimum_homology_length_state=minimum_homology_length_state,
                scale_bar_state=scale_bar_state,
                title_input_state=title_input_state,
            )

        # Erase Plot & Reset All Parameters
        if button_id == "erase-button":
            dash_parameters.reset()
            # Return an empty figure, None for clickdata, and False for skeleton
            return {}, None, False

        # = TAB VIEW =================================================================== #
        # Update Annotations and View
        if button_id == "update-annotations":
            return handle_update_view_click(
                figure_state=figure_state,
                dash_parameters=dash_parameters,
                align_plot_state=align_plot_state,
                homology_style_state=homology_style_state,
                scale_bar_state=scale_bar_state,
                minimum_homology_length_state=minimum_homology_length_state,
            )

        # = TAB EDIT =================================================================== #
        # Change Homology Color Regions and Colorscale Bar Legend
        if button_id == "change-homology-color-button":
            return handle_update_homologies_click(
                figure_state,
                dash_parameters,
                color_scale_state,
                range_slider_state,
                is_set_to_extreme_homologies,
            )

        # Insert title to plot
        if button_id == "update-title-button":
            return handle_update_title_click(
                figure_state=figure_state,
                dash_parameters=dash_parameters,
                title_input_state=title_input_state,
            )

        # Update sequence annotations
        if button_id == "offcanvas-update-sequence-annotations-button":
            fig = Figure(data=figure_state["data"], layout=figure_state["layout"])
            fig = update_dna_sequence_annotations(
                fig=fig,
                dash_parameters=dash_parameters,
                annotation_column_choice_state=annotation_column_choice_state,
                table=sequence_table_state,
            )
            return fig, None, False

        # Update gene annotations
        if button_id == "offcanvas-update-gene-annotations-button":
            print("Updating gene annotations based on user preferences...")
            fig = Figure(data=figure_state["data"], layout=figure_state["layout"])
            fig = update_genes_annotations(
                fig=fig,
                dash_parameters=dash_parameters,
                annotate_genes_from_state=gene_annotation_from_state,
                annotate_genes_positions_state=gene_annotation_positions_state,
                table=gene_table_state,
            )
            return fig, None, False

        # Change Color of Selected Traces
        if button_id == "change-gene-color-button":
            return handle_change_color_click(
                figure_state,
                dash_parameters,
                color_input_state,
            )

        # Select Traces for Changing Color
        if (
            (active_tab == "tab-edit")
            and (click_data is not None)
            and select_items_state
        ):
            return handle_select_traces_click(
                figure_state,
                dash_parameters,
                click_data,
            )

        return figure_state, None, False
