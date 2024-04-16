"""Script to recreate Figure 4b"""

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import pandas as pd

# pylint: disable=unused-import
import scienceplots

# pylint: enable=unused-import
import seaborn as sns

plt.style.use("nature")
plt.rcParams["font.sans-serif"] = "Arial"
FONT_SIZE = 11
plt.rcParams["font.size"] = FONT_SIZE
plt.rcParams["legend.title_fontsize"] = FONT_SIZE
plt.rcParams["legend.fontsize"] = FONT_SIZE
plt.rcParams["axes.labelsize"] = FONT_SIZE
plt.rcParams["xtick.labelsize"] = FONT_SIZE
plt.rcParams["ytick.labelsize"] = FONT_SIZE


def parse_args():
    """Parse input-related options."""
    parser = ArgumentParser(
        description="Arguments for Figure 4b",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-f",
        "--date_filepath",
        type=str,
        help="Date file path",
    )
    return parser.parse_args()


def get_dates(path_to_date_file):
    df = pd.read_csv(path_to_date_file, header=[0, 1, 2]).droplevel(
        level=[1, 2], axis=1
    )

    dates = df["Sampling Dates"].rename(index="dates")

    dates = pd.DataFrame(
        dates.apply(lambda x: x[1:-1].split(",")).to_list(),
        columns=["first_date", "last_date"],
    )

    for col in dates.columns:
        dates[col] = pd.to_datetime(
            dates[col].apply(lambda x: "2021-" + x.strip()), format="%Y-%m-%d"
        )

    dates["Clade"] = df["Clade"]

    dates["variant"] = df["Clade"].apply(lambda x: x.split(" ")[0])

    return dates


def main():
    args = parse_args()
    dates = get_dates(args.date_filepath)

    palettes = {
        "Alpha": "Greens_r",
        "Delta": "autumn",
        "Epsilon": "Blues_r",
        "Omicron": "PiYG",
    }

    colors = {"Beta": "gray", "Eta": "y", "Gamma": "r", "Mu": "purple", "Zeta": "k"}

    counts = {}

    _, ax = plt.subplots(1, 1, figsize=(4.43, 8))

    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b\n%Y"))

    # scale = 1

    for i in range(dates.shape[0]):
        row = dates.iloc[i]
        variant = row.variant

        if variant in counts:
            counts[variant] += 1
        else:
            counts[variant] = 0

        if variant in palettes:
            color = sns.color_palette(palettes[variant])[counts[variant]]
        else:
            color = colors[variant]

        dx = (row.last_date - row.first_date).days

        ax.arrow(
            x=row.first_date,
            dx=dx,
            y=dates.shape[0] - i,
            dy=0,
            head_width=0.25,
            head_length=10,
            color=color,
        )

        ax.arrow(
            x=row.last_date,
            dx=-dx,
            y=dates.shape[0] - i,
            dy=0,
            head_width=0.25,
            head_length=10,
            color=color,
        )

        ax.text(
            x=row.first_date + (row.last_date - row.first_date) / 2,
            y=dates.shape[0] - i + 0.25,
            s=row.Clade,
            fontdict={"color": color, "size": 11.5},
            ha="center",
        )

    plt.ylim((dates.shape[0] - i - 0.8, dates.shape[0] + 0.2))
    sns.despine(left=True)
    ax.set_yticks([], [])
    plt.savefig("img/timelines.pdf", format="pdf", pad_inches=0, bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    main()
