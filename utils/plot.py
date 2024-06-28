"""Plotting utility functions."""

import itertools

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from matplotlib.collections import LineCollection


def clear_axes(
    ax=None, top=True, right=True, left=False, bottom=False, minorticks_off=True
):
    """A more forcing version of sns.despine.

    Parameters
    ----------
    ax : matplotlib axes, optional
        Specific axes object to despine. Ignored if fig is provided.
    top, right, left, bottom : bool, optional
        If True, remove that spine.
    minorticks_off : bool, optional
        If True, remove all minor ticks, by default True.
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
    width : float
        Document textwidth or columnwidth in pts
        Report: 390 pt
    layout : string
        h: horizontal layout
        v: vertical layout
        s: square layout
    fraction : float, optional
        Fraction of the width which you wish the figure to occupy

    Returns
    -------
    fig_dim : tuple
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


# neil IFNDEF 04.09
def plot_tree(
    tree,
    align_names=False,
    name_offset=None,
    font_size=9,
    label_with_node_style=True,
    ax=None,
):
    """
    Plots a ete3.Tree object using matploltib.

    Adapted from: https://gist.github.com/jolespin/5d90deff552138d73de7ed4bdd9ac57a

    Parameters
    ----------
    tree : ete Tree object
    align_names: bool
        If True names will be aligned vertically, by default False
    name_offset : float, optional
        Offset relative to tips to write leaf_names. In BL scale, by default None
    font_size : int, optional
        Text font size, by default 9
    label_with_node_style : bool, optional
        If True, color the node label with the node's style color, by default True
    ax : matplotlib.Axes object, optional
        Object on which the tree will be plotted, by default None

    Returns
    -------
    ax : matplotlib.Axes object
        The matplotlib axes containing the plot.
    """
    shape_dict = {"circle": "o", "square": "s", "sphere": "o"}
    linestyle_dict = dict(enumerate(("-", "--", ":")))

    if ax is None:
        ax = plt.gca()

    aligned_lines = []

    max_x = max(n.get_distance(tree) for n in tree.iter_leaves())

    if name_offset is None:
        name_offset = max_x / 50.0

    node_pos = {n2: i for i, n2 in enumerate(tree.get_leaves()[::-1])}
    node_list = itertools.chain(tree.iter_descendants(strategy="postorder"), [tree])

    # draw tree
    for node in node_list:
        # Parent style
        pstyle = node.img_style

        x = sum(n2.dist for n2 in node.iter_ancestors()) + node.dist

        if node.is_leaf():
            y = node_pos[node]
            if align_names:
                x = max_x
                aligned_lines.append(((x, y), (max_x + name_offset, y)))

        else:
            y = np.mean([node_pos[n2] for n2 in node.children])
            node_pos[node] = y

            # draw vertical line
            ax.plot(
                [x, x],
                [node_pos[node.children[0]], node_pos[node.children[-1]]],
                c=pstyle["vt_line_color"],
                linestyle=linestyle_dict[pstyle["vt_line_type"]],
                linewidth=0.5 * (pstyle["vt_line_width"] + 1),
            )

            # draw horizontal lines
            for child in node.children:
                # Child style
                cstyle = child.img_style
                ax.plot(
                    [x, x + child.dist],
                    [node_pos[child], node_pos[child]],
                    c=cstyle["hz_line_color"],
                    linestyle=linestyle_dict[cstyle["hz_line_type"]],
                    linewidth=0.5 * (cstyle["hz_line_width"] + 1),
                )

        # Node label
        ax.text(
            x + name_offset,
            y,
            node.name,
            va="center",
            size=font_size,
            c=pstyle["fgcolor"] if label_with_node_style else "k",
        )

        # Node point
        ax.scatter(
            x,
            y,
            s=pstyle["size"] ** 2 / 2,
            marker=shape_dict[pstyle["shape"]],
            c=pstyle["fgcolor"],
            zorder=10,
        )

    ali_line_col = LineCollection(aligned_lines, colors="k")

    ax.add_collection(ali_line_col)

    ax.set_axis_off()
    return ax


# neil ENDIF 04.09
