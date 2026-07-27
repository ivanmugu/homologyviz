import dash_mantine_components as dmc

TAB_LABEL_STYLE = {
    "fontSize": "14px",
    "padding": "0.6rem 1rem",
}


def make_dmc_select(**kwargs) -> dmc.Select:
    """
    Create a styled Dash Mantine Components (DMC) Select element.

    This utility function returns a DMC Select component with predefined styling,
    including fixed width, padding, and consistent font size across input, label,
    and options. Additional keyword arguments are passed directly to the `dmc.Select`.

    Parameters
    ----------
    **kwargs : dict
        Additional properties to customize the Select component (e.g., `data`, `value`,
        `label`).

    Returns
    -------
    A configured `dmc.Select` component ready to be used in the Dash layout.
    """
    return dmc.Select(
        w=200,
        size="md",
        style={"padding": "0"},
        styles={
            "input": {"fontSize": "14px"},
            "label": {"fontSize": "14px"},
            "option": {"fontSize": "14px"},
        },
        **kwargs,
    )
