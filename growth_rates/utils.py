# pylint: disable=invalid-name
import numpy as np


def get_lineage_counts(data, dict_lineage, dict_time):
    day_col, lineage_col = data.columns

    day_counts = data.groupby(day_col).count()
    column_name = {lineage_col: "All"}
    day_counts.rename(columns=column_name, inplace=True)

    for key, value in dict_lineage.items():
        # extract rows belonging to each lineage from dataframe
        # requires a grouping of sublineages into major lineages using dict_lineage
        mask = data[lineage_col].isin(value)

        # Group by each day
        df_grouped = data[mask].groupby(day_col)
        # Total counts on each day for full dataframe - this is incidence
        Z = df_grouped.count()
        # Start and end date for each lineage
        t0, t1 = dict_time[key]
        # Incidence from start to end of each lineage
        Z = Z.loc[t0:t1, :]
        columns_name = {lineage_col: key}
        Z.rename(columns=columns_name, inplace=True)
        # Create query for dataframe - not the NA values
        query = f"`{key}`.notna()"

        # Find first and last non-NA values
        indices = Z.query(query, engine="python").index
        xmin, xmax = indices[[0, -1]]

        # Create new time frame, filling in missing values with zeros
        Z = Z.loc[xmin : xmax + 1, :]
        Z = Z.reindex(np.arange(xmin, xmax + 1))
        Z.fillna(0, inplace=True)

        # Final df to return
        day_counts = day_counts.join(Z, how="outer")
    day_counts.fillna(0, inplace=True)
    day_columns = day_counts.columns
    day_counts["Other"] = day_counts["All"] - day_counts[day_columns[1:]].sum(axis=1)
    return day_counts.astype(int)


# def get_lineage_probabilities(data, dict_lineage, dict_time):
#     df_probabilities = pd.DataFrame()
#     df = get_lineage_counts(data, dict_lineage, dict_time)
#     for lineage in df.columns[1:]:
#         p_lineage = df.eval(f"`{lineage}` / All")
#         p_lineage.name = f"P_{lineage}"
#         df_probabilities = df_probabilities.join(p_lineage, how="outer")
#     return df_probabilities


def get_lineage_occurences(data, dict_lineage):
    day_col, lineage_col = data.columns
    occurences = {}

    for lineage in dict_lineage.keys():

        query = f"{lineage_col} in @lineage_list"

        # pylint: disable=unused-variable
        x_min, x_max = data.query(query)[day_col].agg([min, max]).T
        # pylint: enable=unused-variable

        df_lineage = data.query(f"@x_min <= {day_col} <= @x_max")

        eval_ = f"{lineage_col} = {query}"
        df_lineage.eval(eval_, inplace=True)
        df_lineage.reset_index(drop=True, inplace=True)
        occurences[lineage] = df_lineage

    return occurences
