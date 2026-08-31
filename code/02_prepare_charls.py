from pathlib import Path
import json
import os

import numpy as np
import pandas as pd
from patsy import dmatrices
from scipy.optimize import minimize
from scipy.stats import norm


ROOT = Path(os.environ.get("PEF_PROJECT_ROOT", ".")).resolve()
DATA = Path(os.environ["CHARLS_DATA_FILE"]).expanduser().resolve()
OUT = ROOT / "outputs" / "charls_preparation"
OUT.mkdir(parents=True, exist_ok=True)

COLS = [
    "ID", "communityID", "inw3", "inw4", "r3agey", "ragender",
    "raeduc_c", "h3rural", "h3atotb", "r3mheight", "r3mbmi",
    "r3puff", "r3puffeff", "r3lunge", "r4lunge", "r3asthmae",
    "r3smokev", "r3smoken", "r3vgact_c", "r3vgactx_c",
    "r3mdact_c", "r3mdactx_c", "r3hibpe", "r3diabe", "r3hearte",
    "r3stroke", "r3cancre", "r3cesd10", "r3wtrespb",
]


def binary(s):
    x = pd.to_numeric(s, errors="coerce")
    return x.where(x.isin([0, 1]))


def poisson_cluster(x, y, weights, clusters):
    weights = np.asarray(weights, float)
    weights = weights / np.mean(weights)
    y = np.asarray(y, float)
    x = np.asarray(x, float)

    def objective(beta):
        eta = np.clip(x @ beta, -20, 10)
        mu = np.exp(eta)
        return float(np.sum(weights * (mu - y * eta)))

    def gradient(beta):
        eta = np.clip(x @ beta, -20, 10)
        return x.T @ (weights * (np.exp(eta) - y))

    start = np.zeros(x.shape[1])
    start[0] = np.log(max(np.average(y, weights=weights), 1e-5))
    fit = minimize(
        objective, start, jac=gradient, method="L-BFGS-B",
        options={"maxiter": 3000, "ftol": 1e-12, "gtol": 1e-8},
    )
    beta = fit.x
    mu = np.exp(np.clip(x @ beta, -20, 10))
    bread = np.linalg.pinv(x.T @ ((weights * mu)[:, None] * x))
    scores = x * (weights * (y - mu))[:, None]
    unique, inv = np.unique(np.asarray(clusters).astype(str), return_inverse=True)
    summed = np.zeros((len(unique), x.shape[1]))
    np.add.at(summed, inv, scores)
    meat = summed.T @ summed
    n, p, g = len(y), x.shape[1], len(unique)
    correction = (g / (g - 1)) * ((n - 1) / (n - p)) if g > 1 and n > p else 1.0
    cov = bread @ meat @ bread * correction
    return beta, cov, bool(fit.success), int(fit.nit), len(unique)


def rank_safe_design(xdf, forced_terms):
    forced = [c for c in xdf.columns if c in forced_terms]
    selected = forced.copy()
    current_rank = np.linalg.matrix_rank(xdf[selected].to_numpy(float))
    for c in [x for x in xdf.columns if x not in selected]:
        candidate = xdf[selected + [c]].to_numpy(float)
        new_rank = np.linalg.matrix_rank(candidate)
        if new_rank > current_rank:
            selected.append(c)
            current_rank = new_rank
    return xdf[selected]


raw = pd.read_stata(DATA, columns=COLS, convert_categoricals=False)

risk = (
    raw["inw3"].eq(1) & raw["inw4"].eq(1)
    & raw["r3agey"].between(45, 80)
    & raw["r3puff"].between(30, 900)
    & raw["r3mheight"].between(1.2, 2.2)
    & raw["ragender"].isin([1, 2])
    & binary(raw["r3lunge"]).eq(0)
    & binary(raw["r3asthmae"]).eq(0)
    & binary(raw["r4lunge"]).notna()
    & raw["r3wtrespb"].gt(0)
    & raw["communityID"].notna()
)
df = raw.loc[risk].copy()
df["event"] = binary(df["r4lunge"]).astype(int)
df["age"] = pd.to_numeric(df["r3agey"], errors="coerce")
df["sex"] = df["ragender"]
df["height_m"] = pd.to_numeric(df["r3mheight"], errors="coerce")
df["pef"] = pd.to_numeric(df["r3puff"], errors="coerce")
df["race_context"] = binary(df["h3rural"])
df["education"] = pd.to_numeric(df["raeduc_c"], errors="coerce")
df["bmi"] = pd.to_numeric(df["r3mbmi"], errors="coerce").where(lambda x: x.between(12, 70))
df["hypertension"] = binary(df["r3hibpe"])
df["diabetes"] = binary(df["r3diabe"])
df["heart"] = binary(df["r3hearte"])
df["stroke"] = binary(df["r3stroke"])
df["cancer"] = binary(df["r3cancre"])
df["cesd"] = pd.to_numeric(df["r3cesd10"], errors="coerce").where(lambda x: x.between(0, 30))
df["weight"] = pd.to_numeric(df["r3wtrespb"], errors="coerce")
df["age_c10"] = (df["age"] - 60) / 10
df["age_c10_sq"] = df["age_c10"] ** 2
df["height_c10"] = (df["height_m"] - 1.65) / 0.10
df["bmi_c5"] = (df["bmi"] - 24) / 5
df["cesd_per5"] = df["cesd"] / 5

# Smoking: never/former/current.
ever = binary(df["r3smokev"])
now = binary(df["r3smoken"])
df["smoking"] = np.select(
    [now.eq(1), now.eq(0) & ever.eq(1), now.eq(0) & ever.eq(0)],
    ["current", "former", "never"], default=None,
)

# Physical activity aligned approximately to the HRS frequency categories.
vig_any = binary(df["r3vgact_c"])
mod_any = binary(df["r3mdact_c"])
vig_days = pd.to_numeric(df["r3vgactx_c"], errors="coerce").where(lambda x: x.between(0, 7))
mod_days = pd.to_numeric(df["r3mdactx_c"], errors="coerce").where(lambda x: x.between(0, 7))
vig_days = vig_days.where(vig_any.ne(0), 0)
mod_days = mod_days.where(mod_any.ne(0), 0)
best_days = pd.concat([vig_days, mod_days], axis=1).max(axis=1, skipna=False)
df["activity"] = np.select(
    [best_days.ge(2), best_days.eq(1), best_days.eq(0)],
    ["frequent", "some", "inactive"], default=None,
)

# Baseline wealth quintile (rank based, robust to negative wealth).
wealth = pd.to_numeric(df["h3atotb"], errors="coerce")
df["wealth_q"] = np.nan
ok = wealth.notna()
df.loc[ok, "wealth_q"] = np.ceil(wealth[ok].rank(method="average", pct=True) * 5).clip(1, 5)

# Harmonized PEF definition: sex-specific log(PEF) residual adjusted for age,
# age squared and measured height; bottom residual quintile = low PEF.
df["pef_residual"] = np.nan
df["low_pef_adj"] = np.nan
df["pef_lower_1sd"] = np.nan
pef_definition = {}
for sex, g in df.groupby("sex"):
    x = np.column_stack([np.ones(len(g)), g["age"], g["age"] ** 2, g["height_m"]])
    y = np.log(g["pef"].to_numpy())
    beta, *_ = np.linalg.lstsq(x, y, rcond=None)
    residual = y - x @ beta
    cutoff = float(np.quantile(residual, .20))
    sd = float(np.std(residual, ddof=1))
    df.loc[g.index, "pef_residual"] = residual
    df.loc[g.index, "low_pef_adj"] = (residual <= cutoff).astype(int)
    df.loc[g.index, "pef_lower_1sd"] = -residual / sd
    pef_definition[str(int(sex))] = {"residual_q20": cutoff, "residual_sd": sd, "n": len(g)}

FORMULAS = {
    "M0": "event ~ exposure",
    "M1": "event ~ exposure + age_c10 + age_c10_sq + sex + height_c10",
    "M2": (
        "event ~ exposure + age_c10 + age_c10_sq + sex + height_c10 + C(race_context) "
        "+ C(education) + C(smoking) + bmi_c5"
    ),
    "M3": (
        "event ~ exposure + age_c10 + age_c10_sq + sex + height_c10 + C(race_context) "
        "+ C(education) + C(smoking) + bmi_c5 "
        "+ hypertension + diabetes + heart + stroke + cancer + cesd_per5"
    ),
    "M4_expanded": (
        "event ~ exposure + age_c10 + age_c10_sq + sex + height_c10 + C(race_context) "
        "+ C(education) + C(smoking) + bmi_c5 + hypertension + diabetes + heart + stroke "
        "+ cancer + cesd_per5 + C(wealth_q) + C(activity)"
    ),
}


def fit_set(data, tag, weighted=True):
    output = []
    for exposure in ["low_pef_adj", "pef_lower_1sd"]:
        for model, template in FORMULAS.items():
            formula = template.replace("exposure", exposure)
            ydf, xdf = dmatrices(formula, data=data, return_type="dataframe", NA_action="drop")
            xdf = rank_safe_design(xdf, ["Intercept", exposure])
            use = xdf.index
            weights = data.loc[use, "weight"].to_numpy() if weighted else np.ones(len(use))
            beta, cov, converged, iterations, clusters = poisson_cluster(
                xdf.to_numpy(), ydf.iloc[:, 0].to_numpy(), weights,
                data.loc[use, "communityID"].astype(str).to_numpy(),
            )
            pos = list(xdf.columns).index(exposure)
            b = float(beta[pos])
            se = float(np.sqrt(max(cov[pos, pos], 0)))
            output.append({
                "analysis": tag,
                "weighted": weighted,
                "exposure": exposure,
                "model": model,
                "n": len(use),
                "events": int(ydf.iloc[:, 0].sum()),
                "clusters": clusters,
                "rr": float(np.exp(b)),
                "ci_low": float(np.exp(b - 1.96 * se)),
                "ci_high": float(np.exp(b + 1.96 * se)),
                "p": float(2 * norm.sf(abs(b / se))) if se > 0 else np.nan,
                "converged": converged,
                "iterations": iterations,
                "formula": formula,
            })
    return output


results = []
results += fit_set(df, "primary", weighted=True)
results += fit_set(df, "primary_unweighted", weighted=False)
results += fit_set(df[df["r3puffeff"].eq(1)], "good_effort_only", weighted=True)
results += fit_set(df[df["pef"].between(50, 800)], "pef_50_to_800", weighted=True)
results = pd.DataFrame(results)

variables = [
    "low_pef_adj", "pef_lower_1sd", "race_context", "education", "wealth_q",
    "smoking", "bmi", "activity", "hypertension", "diabetes", "heart",
    "stroke", "cancer", "cesd",
]
missing = pd.DataFrame({
    "variable": variables,
    "nonmissing_n": [int(df[v].notna().sum()) for v in variables],
    "missing_n": [int(df[v].isna().sum()) for v in variables],
    "missing_pct": [float(df[v].isna().mean() * 100) for v in variables],
})

groups = df.groupby("low_pef_adj")["event"].agg(n="size", events="sum", risk="mean").reset_index()

df.to_csv(OUT / "charls_2015_2018_harmonized_dataset.csv", index=False, encoding="utf-8-sig")
results.to_csv(OUT / "modified_poisson_results.csv", index=False, encoding="utf-8-sig")
missing.to_csv(OUT / "missingness.csv", index=False, encoding="utf-8-sig")
groups.to_csv(OUT / "event_counts_by_low_pef.csv", index=False, encoding="utf-8-sig")

primary_m3 = results.query("analysis == 'primary' and model == 'M3'")
summary = {
    "risk_set_n": int(len(df)),
    "events": int(df["event"].sum()),
    "event_risk": float(df["event"].mean()),
    "primary_m3_n": int(primary_m3.iloc[0]["n"]),
    "primary_m3_events": int(primary_m3.iloc[0]["events"]),
    "pef_definition": pef_definition,
    "weight": "r3wtrespb individual weight with household and individual nonresponse adjustment",
    "cluster": "communityID",
    "note": "r3wtrespbiob was not used because it is a biomarker-subsample weight, whereas PEF is a physical measure.",
}
(OUT / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")

print(json.dumps(summary, ensure_ascii=False, indent=2))
print("\nEvent counts:\n", groups.to_string(index=False))
print("\nMissingness:\n", missing.to_string(index=False))
print("\nModels:\n", results.to_string(index=False))
