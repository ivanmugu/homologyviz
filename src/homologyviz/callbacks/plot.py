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


def figure_from_state(figure_state: dict) -> Figure:
    """
    Rebuild a Plotly Figure from its Dash figure state.

    Parameters
    ----------
    figure_state : dict
        Serialized Plotly figure containing ``data`` and ``layout``.

    Returns
    -------
    plotly.graph_objects.Figure
        Reconstructed Plotly figure.
    """
    return Figure(
        data=figure_state["data"],
        layout=figure_state["layout"],
    )


def handle_sequence_annotations_update(
    figure_state: dict,
    dash_parameters: PlotParameters,
    annotation_column_choice_state: str,
    sequence_table_state: list[dict] | None,
) -> tuple[Figure, None, bool]:
    """Apply sequence annotation edits to the current figure."""
    fig = figure_from_state(figure_state)
    fig = update_dna_sequence_annotations(
        fig=fig,
        dash_parameters=dash_parameters,
        annotation_column_choice_state=annotation_column_choice_state,
        table=sequence_table_state,
    )
    return fig, None, False


def handle_gene_annotations_update(
    figure_state: dict,
    dash_parameters: PlotParameters,
    gene_annotation_from_state: str,
    gene_annotation_positions_state: str,
    gene_table_state: list[dict] | None,
) -> tuple[Figure, None, bool]:
    """Apply gene annotation edits to the current figure."""
    fig = figure_from_state(figure_state)
    fig = update_genes_annotations(
        fig=fig,
        dash_parameters=dash_parameters,
        annotate_genes_from_state=gene_annotation_from_state,
        annotate_genes_positions_state=gene_annotation_positions_state,
        table=gene_table_state,
    )
    return fig, None, False


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
        output=[
            Output("plot", "figure"),
            Output("plot", "clickData"),
            Output("plot-skeleton", "visible"),
        ],
        inputs=dict(
            _actions=dict(
                plot=Input("plot-button", "n_clicks"),
                erase=Input("erase-button", "n_clicks"),
                change_homology_color=Input("change-homology-color-button", "n_clicks"),
                change_gene_color=Input("change-gene-color-button", "n_clicks"),
                update_view=Input("update-annotations", "n_clicks"),
                update_title=Input("update-title-button", "n_clicks"),
                update_sequence_annotations=Input(
                    "offcanvas-update-sequence-annotations-button", "n_clicks"
                ),
                update_gene_annotations=Input(
                    "offcanvas-update-gene-annotations-button", "n_clicks"
                ),
            ),
            click_data=Input("plot", "clickData"),
        ),
        state=dict(
            virtual=State("files-table", "virtualRowData"),
            active_tab=State("tabs", "active_tab"),
            figure_state=State("plot", "figure"),
            title_input_state=State("title-input", "value"),
            edit=dict(
                color=State("color-input", "value"),
                select_item=State("select-items-button-store", "data"),
            ),
            view=dict(
                align_plot=State("align-plot", "value"),
                homology_style=State("homology-style", "value"),
                minimum_homology_length=State("minimum-homology-length", "value"),
                scale_bar=State("scale-bar", "value"),
            ),
            annotations=dict(
                gene_positions=State("gene-annotation-positions-choice", "value"),
                gene_from=State("gene-annotation-from-choice", "value"),
                gene_table=State("gene-table", "virtualRowData"),
                sequence_table=State("sequence-table", "virtualRowData"),
                sequence_column=State("annotation-column-choice", "value"),
            ),
            homology_colors=dict(
                color_scale=State("color-scale", "value"),
                identity_range=State("range-slider", "value"),
                use_extreme=State("is-set-to-extreme-homologies", "data"),
            ),
        ),
        prevent_initial_call=True,
    )
    def main_plot(
        _actions: dict,
        click_data: dict | None,
        virtual: list[dict[str, str]] | None,
        active_tab: str,
        figure_state: dict,
        title_input_state: str,
        edit: dict,
        view: dict,
        annotations: dict,
        homology_colors: dict,
    ) -> tuple[Figure | dict, None, bool]:
        """
        Dispatch user actions that create or modify the alignment plot.

        This callback is the main controller for interactions with the HomologyViz
        alignment figure. It determines which input triggered the callback using
        ``dash.ctx.triggered_id`` and delegates the requested operation to the
        appropriate handler.

        Supported actions include generating and erasing the plot, updating the
        plot view, changing homology and gene colors, editing the title, applying
        sequence and gene annotations, and selecting traces from the figure.

        Parameters
        ----------
        _actions : dict
            Grouped button inputs that trigger plot-related actions. Their click
            counts are not used directly; the triggering component is identified
            with ``dash.ctx.triggered_id``.
        click_data : dict or None
            Data generated when the user clicks a trace in the Plotly figure.
            Used for trace selection when selection mode is active.
        virtual : list[dict[str, str]]
            Current virtual rows of the uploaded-files table, including the files
            used to generate the alignment plot.
        active_tab : str
            Identifier of the currently active application tab.
        figure_state : dict
            Current serialized state of the Plotly alignment figure.
        title_input_state : str
            Current value of the figure title input.
        edit : dict
            Settings related to interactive editing. Contains the selected color
            and the current trace-selection state.
        homology_colors : dict
            Homology coloring settings, including the colorscale, identity range,
            and whether the colorscale uses the extreme homology values.
        view : dict
            Plot display settings, including sequence alignment, homology style,
            minimum homology length, and scale-bar visibility.
        annotations : dict
            Sequence and gene annotation settings and their editable table data.

        Returns
        -------
        tuple[Figure | dict, None, bool]
            The updated Plotly figure (or serialized figure state), ``None`` to
            clear ``clickData``, and ``False`` to hide the plot skeleton loader.

        Notes
        -----
        The callback primarily acts as a dispatcher. Plot generation and editing
        operations are delegated to dedicated handler functions rather than being
        implemented directly in this callback.
        """
        # Use context to find the button that triggered the callback.
        button_id = dash.callback_context.triggered_id

        if (button_id == "plot-button") and virtual:
            return handle_plot_button_click(
                dash_parameters=dash_parameters,
                virtual=virtual,
                tmp_directory_path=tmp_directory,
                align_plot_state=view["align_plot"],
                color_scale_state=homology_colors["color_scale"],
                range_slider_state=homology_colors["identity_range"],
                is_set_to_extreme_homologies=homology_colors["use_extreme"],
                annotation_column_choice_state=annotations["sequence_column"],
                annotate_genes_from_state=annotations["gene_from"],
                annotate_genes_positions_state=annotations["gene_positions"],
                homology_style_state=view["homology_style"],
                minimum_homology_length_state=view["minimum_homology_length"],
                scale_bar_state=view["scale_bar"],
                title_input_state=title_input_state,
            )

        if button_id == "erase-button":
            dash_parameters.reset()
            return {}, None, False

        if button_id == "update-annotations":
            return handle_update_view_click(
                figure_state=figure_state,
                dash_parameters=dash_parameters,
                align_plot_state=view["align_plot"],
                homology_style_state=view["homology_style"],
                scale_bar_state=view["scale_bar"],
                minimum_homology_length_state=view["minimum_homology_length"],
            )

        if button_id == "change-homology-color-button":
            return handle_update_homologies_click(
                figure_state=figure_state,
                dash_parameters=dash_parameters,
                color_scale_state=homology_colors["color_scale"],
                range_slider_state=homology_colors["identity_range"],
                is_set_to_extreme_homologies=homology_colors["use_extreme"],
            )

        if button_id == "update-title-button":
            return handle_update_title_click(
                figure_state=figure_state,
                dash_parameters=dash_parameters,
                title_input_state=title_input_state,
            )

        if button_id == "offcanvas-update-sequence-annotations-button":
            return handle_sequence_annotations_update(
                figure_state=figure_state,
                dash_parameters=dash_parameters,
                annotation_column_choice_state=annotations["sequence_column"],
                sequence_table_state=annotations["sequence_table"],
            )

        if button_id == "offcanvas-update-gene-annotations-button":
            return handle_gene_annotations_update(
                figure_state=figure_state,
                dash_parameters=dash_parameters,
                gene_annotation_from_state=annotations["gene_from"],
                gene_annotation_positions_state=annotations["gene_positions"],
                gene_table_state=annotations["gene_table"],
            )

        if button_id == "change-gene-color-button":
            return handle_change_color_click(
                figure_state=figure_state,
                dash_parameters=dash_parameters,
                color_input_state=edit["color"],
            )

        if (
            (active_tab == "tab-edit")
            and (click_data is not None)
            and edit["select_item"]
        ):
            return handle_select_traces_click(
                figure_state=figure_state,
                dash_parameters=dash_parameters,
                click_data=click_data,
            )

        return figure_state, None, False
