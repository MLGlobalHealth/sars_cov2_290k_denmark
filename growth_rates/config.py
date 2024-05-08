from typing import Final

REGION_DICT: Final = {
    "Nordjylland": 1081,
    "Midtjylland": 1082,
    "Syddanmark": 1083,
    "Hovedstaden": 1084,
    "Sjælland": 1085,
}

STEM_LINEAGES: Final = {
    "Alpha - B.1.1.7": ["B.1.1.7", "Q"],
    "Delta - B.1.617.2": [
        "B.1.617.2",
        "AY",
    ],
    "Omicron - BA.1": [
        "B.1.1.529.1",
        "BA.1",
    ],
    "Omicron - BA.2": [
        "B.1.1.529.2",
        "BA.2",
    ],
    "Epsilon": ["B.1.429", "B.1.427"],
    "Eta": ["B.1.525"],
    "Gamma": [
        "P",
    ],
}
