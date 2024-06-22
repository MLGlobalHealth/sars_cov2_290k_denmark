"""Utility functions for plotting"""

# Author: Neil Scheidwasser (neil.clow@sund.ku.dk)

import calendar
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from utils.config import VARIANT_PALETTE
from utils.plot import clear_axes, set_size


def _get_name_and_label(x, nnz):
    if x == "n_changes":
        xname = "bl"
        xlabel = "# Subsitutions"
    elif x == "rate":
        xname = "rate"
        xlabel = "Substitution rate \n [N/(site * year)]"
    else:
        raise ValueError

    if nnz:
        xname = "nnz_" + xname

    return xname, xlabel


def pointplot(
    df,
    x,
    y,
    nnz=False,
    xlabel=None,
    rotate=False,
    ax=None,
    # show=False,
    **kwargs,
):
    ax = sns.pointplot(
        x=x,
        y=y,
        data=df.query(f"{y} > 0") if nnz else df,
        order=sorted(df.loc[:, x].unique()),
        linestyles="none",
        errwidth=2,
        scale=0.8,
        ax=ax,
        **kwargs,
    )

    if y == "rate":
        ax.ticklabel_format(style="sci", axis="y", scilimits=(0, 0), useOffset=False)

    if rotate:
        ax.tick_params(axis="x", rotation=45)

    if xlabel is None:
        xlabel = x.capitalize()

    _, ylabel = _get_name_and_label(y, nnz)

    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)

    clear_axes()

    return ax


def month_plot(df, y, nnz=False):
    """Point plot between-month differences in y

    Parameters
    ----------
    df : pandas.DataFrame
        Data
    y : str
        regressand
    nnz : bool, optional
        If True, only keep non-zero values, by default False
    """
    _, ax = plt.subplots(figsize=set_size(240, "h"))

    sns.pointplot(
        x="month",
        y=y,
        data=df.query(f"{y} > 0") if nnz else df,
        order=calendar.month_abbr[1:],
        n_boot=100,
        color="k",
        linestyles="none",
        ax=ax,
    )

    yname, ylabel = _get_name_and_label(y, nnz)

    if y == "rate":
        ax.ticklabel_format(style="sci", axis="y", scilimits=(0, 0), useOffset=False)

    ax.set_xlabel("")
    ax.set_ylabel(rf"{ylabel} \n (mean +/- 95\% CI)")
    clear_axes()

    ax.tick_params(axis="x", rotation=45)
    ymin, ymax = 0, 8
    alpha = 0.25

    spans = [
        (0, 1, "Wild type"),
        (1, 5.5, "Alpha"),
        (5.5, 10.5, "Delta"),
        (10.5, 11, "Omicron"),
    ]

    for span in spans:
        xmin, xmax, key = span

        ax.axvspan(
            xmin=xmin,
            xmax=xmax,
            ymin=ymin,
            ymax=ymax,
            alpha=alpha,
            color=VARIANT_PALETTE[key],
            label=key,
        )

    ax.legend(loc="center", bbox_to_anchor=(0.5, 1.21), ncol=2)

    plt.savefig(
        f"img/{yname}_vs_month.pdf", format="pdf", bbox_inches="tight", pad_inches=0.1
    )
    plt.show()


def reg_coef_plot(params, group, ref, x="var_group", y="coef", ax=None, **kwargs):
    """Plot regression plots for groups within a variable

    Parameters
    ----------
    params : pandas.DataFrame
        Table with coefs, p-value, CIs etc. for each parametrer
    group : str
        Grouping variable
    x : str, optional
        Group column, by default "var_group"
    y : str, optional
        Coef column, by default "coef"
    palette : str, list, or dict, optional
        Colors to use for the different levels, by default None
    ax : matplotlib.axes.Axes
        Axes object to draw the plot onto
    model_name : str, optional
        Regression model name, by default ""
    """
    # Prepare data
    sub_params = params.query("var_name == @group")

    # Main plot
    ax = sns.pointplot(
        x=x,
        y=y,
        data=sub_params,
        linestyles="none",
        color="k",
        scale=0.8,
        ax=ax,
        **kwargs,
    )

    # Error bars
    x_coords = np.arange(sub_params.shape[0])
    y_coords = sub_params.coef
    y_err = sub_params.error
    ax.errorbar(
        x_coords,
        y_coords,
        yerr=y_err,
        color="k",
        fmt="none",
        elinewidth=1.5,
        alpha=0.8,
    )

    # Labels
    x_name = group.split("_")[0]
    ax.set_xlabel(
        f"{x_name.capitalize().replace('Vacc', 'Vaccination status')}\n(Ref. group: {ref})"
    )

    ax.ticklabel_format(style="sci", axis="y", scilimits=(0, 0), useOffset=False)
    clear_axes()

    return ax
