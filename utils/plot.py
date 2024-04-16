import matplotlib.pyplot as plt
import seaborn as sns


def clear_axes(
    ax=None, top=True, right=True, left=False, bottom=False, minorticks_off=True
):
    """A more forcing version of sns.despine.
    ax : matplotlib axes, optional
        Specific axes object to despine. Ignored if fig is provided.
    top, right, left, bottom : boolean, optional
        If True, remove that spine.
    minorticks_off: boolean, optional
        If True, remove all minor ticks
    """
    if ax is None:
        axes = plt.gcf().axes
    else:
        axes = [ax]

    for ax_i in axes:
        sns.despine(ax=ax_i, top=top, right=right, left=left, bottom=bottom)
        if minorticks_off:
            ax_i.minorticks_off()
        ax_i.tick_params(axis="x", which="both", top=not top)
        ax_i.tick_params(axis="y", which="both", right=not right)
        ax_i.tick_params(axis="y", which="both", left=not left)
        ax_i.tick_params(axis="x", which="both", bottom=not bottom)


def set_size(width, layout="h", fraction=1):
    """Set figure dimensions to avoid scaling in LaTeX.

    Heavily inspired by: https://jwalton.info/Embed-Publication-Matplotlib-Latex/

    Parameters
    ----------
    width: float
        Document textwidth or columnwidth in pts
        Report: 390 pt
    layout: string
        h: horizontal layout
        v: vertical layout
        s: square layout
    fraction: float, optional
        Fraction of the width which you wish the figure to occupy

    Returns
    -------
    fig_dim: tuple
        Dimensions of figure in inches (width, height)
    """
    # Width of figure (in pts)
    fig_width_pt = width * fraction

    # Convert from pt to inches
    inches_per_pt = 1 / 72.27

    # Golden ratio for aesthetic figures
    # https://disq.us/p/2940ij3
    golden_ratio = (5**0.5 - 1) / 2

    # Figure width in inches
    fig_width_in = fig_width_pt * inches_per_pt
    # Figure height in inches
    if layout == "h":
        fig_height_in = fig_width_in * golden_ratio
    elif layout == "v":
        fig_height_in = fig_width_in / golden_ratio
    elif layout == "s":
        fig_height_in = fig_width_in

    fig_dim = (fig_width_in, fig_height_in)

    return fig_dim
