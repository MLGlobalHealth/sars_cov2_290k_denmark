"""Constants, encodings, and colour palettes"""

from typing import Final

N_DAYS: Final = 365
DK_POP: Final = 5840045  # Population of Denmark
REF_LEN: Final = 29891  # Reference SARS-CoV 2 sequence length

# Abbreviations for regions
REGION_MAPPING: Final = {
    "Hovedstaden": "H",
    "Midtjylland": "M",
    "Nordjylland": "N",
    "Sjælland": "SJ",
    "Syddanmark": "SY",
}

# Encoding for vaccination status
VACC_MAPPING: Final = {"0": "Unvacc.", "1": "Partial", "2": "Full"}

# Viridis
REGION_PALETTE: Final = {
    "H": "#440154FF",
    "M": "#3B528BFF",
    "N": "#21908CFF",
    "SJ": "#5DC863FF",
    "SY": "#FDE725FF",
}

# gnuplot
VACC_PALETTE: Final = {"Unvacc.": "#6001c6", "Partial": "#a71470", "Full": "#d75d00"}

# colorblind
VARIANT_PALETTE: Final = {
    "Wild type": "#CC78BC",
    "Alpha": "#0173B2",
    "Delta": "#DE8F05",
    "Omicron": "#029E73",
    "Other": "k",
}
