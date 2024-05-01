"""Utility functions for (regression) modelling"""

# Author: Neil Scheidwasser (neil.clow@sund.ku.dk)

import pandas as pd
import statsmodels.api as sm

BASE_RHS = "age_groups + C(major_major_variant) + regions2 + sex + C(vacc_status2)"


def _load_regressor(model_name):
    """Load a statsmodels model for regression

    Parameters
    ----------
    model_name : str
        Model name

    Returns
    -------
    function
        A regressor from statsmodels.api
    """
    if model_name == "ols":
        return sm.OLS
    elif model_name == "zip":
        return sm.ZeroInflatedPoisson
    elif model_name == "zinb":
        return sm.ZeroInflatedNegativeBinomialP
    elif model_name == "zigp":
        return sm.ZeroInflatedGeneralizedPoisson
    else:
        raise ValueError(f"No model found under name {model_name}.")


def regress(
    model_name, df, y, extra_formula=None, extra_factors=None, regularized=False
):
    """Regression of branch length information vs. other common variables

    Parameters
    ----------
    model_name : str
        Regression model name
    df : pandas.DataFrame
        Data
    y : str
        Regressand
    extra_formula : str, optional
        A string to add extra variables to the base formula (in patsy), by default None
    extra_factors : list, optional
        A list of extra factors to add to the base formula, by default None
    regularized : bool, optional
        if True, perform a regularized fit (using regularized MLE), by default False

    Returns
    -------
    res : statsmodels.iolab.summary.Summary
        Summary of the model fit
    params : pandas.DataFrame
        Table with coefs, p-value, CIs etc. for each parametrer
    """
    rhs = BASE_RHS

    # Adjust formula depending on extra variables
    if extra_formula is not None:
        rhs += f"+ {extra_formula}"
    if extra_factors is not None:
        rhs += f"+ {'+'.join(extra_factors)}"

    # Adjust formula depending on model
    if model_name in ["zip", "zinb", "zigp"]:
        # ZIP/ZINB need an int regressand (count data)
        df[f"{y}_int"] = df.loc[:, y].astype(int)
        formula = f"{y}_int ~ {rhs}"
    else:
        formula = f"{y} ~ {rhs}"

    model = _load_regressor(model_name).from_formula(formula, data=df)

    # Fit model
    if regularized:
        res = model.fit_regularized().summary()
    else:
        res = model.fit().summary()

    # Create param table
    params = pd.read_html(res.tables[1].as_html(), header=0, index_col=0)[0]

    params.reset_index(names="variable", inplace=True)

    params["error"] = 1.96 * params["std err"]

    params[["var_name", "var_group"]] = params["variable"].str.split("[", expand=True)

    params["var_group"] = params["var_group"].apply(
        lambda x: str(x).replace("T.", "").replace("]", "")
    )

    params["var_name"] = params["var_name"].apply(
        lambda x: x.split(",")[0].replace("C(", "").replace(")", "")
    )

    return res, params
