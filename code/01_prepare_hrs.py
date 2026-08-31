from pathlib import Path
import json
import os
import warnings

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf


ROOT = Path(os.environ.get("PEF_PROJECT_ROOT", ".")).resolve()
DATA = Path(os.environ["HRS_DATA_FILE"]).expanduser().resolve()
OUT = ROOT / "outputs" / "hrs_preparation"
OUT.mkdir(parents=True, exist_ok=True)

BASE_WAVES = list(range(8, 15))
YEAR = {w: 1990 + 2 * w for w in range(1, 17)}

fixed = ["hhidpn", "ragender", "raracem", "raedyrs", "raestrat", "raehsamp"]
wave_stems = [
    "agey_e", "lunge", "lungs", "puff", "pmhght", "bmi", "smokev",
    "smoken", "vgactx", "mdactx", "hibpe", "diabe", "hearte",
    "stroke", "cancre", "cesd", "wtresp",
]
cols = fixed[:]
for w in range(8, 17):
    cols += [f"inw{w}", f"r{w}lunge", f"r{w}lungs"]
for w in BASE_WAVES:
    cols += [f"r{w}{stem}" for stem in wave_stems]
    cols += [f"h{w}atotb"]
cols = list(dict.fromkeys(cols))

raw = pd.read_stata(DATA, columns=cols, convert_categoricals=False)


def binary(s):
    return pd.to_numeric(s, errors="coerce").where(lambda x: x.isin([0, 1]))


wave_frames = []
for w in BASE_WAVES:
    age = pd.to_numeric(raw[f"r{w}agey_e"], errors="coerce")
    pef = pd.to_numeric(raw[f"r{w}puff"], errors="coerce")
    height = pd.to_numeric(raw[f"r{w}pmhght"], errors="coerce")
    lung0 = binary(raw[f"r{w}lunge"])
    lung1 = binary(raw[f"r{w+1}lunge"])
    risk = (
        raw[f"inw{w}"].eq(1) & raw[f"inw{w+1}"].eq(1)
        & age.between(45, 80) & pef.between(30, 900)
        & lung0.eq(0) & lung1.notna()
    )
    ix = np.flatnonzero(risk.to_numpy())
    smoke_ever = binary(raw[f"r{w}smokev"]).iloc[ix].to_numpy()
    smoke_now = binary(raw[f"r{w}smoken"]).iloc[ix].to_numpy()
    smoke_complete = ~pd.isna(smoke_ever) & ~pd.isna(smoke_now)
    smoking = np.full(len(ix), None, dtype=object)
    smoking[smoke_complete & (smoke_now == 1)] = "current"
    smoking[smoke_complete & (smoke_now == 0) & (smoke_ever == 1)] = "former"
    smoking[smoke_complete & (smoke_now == 0) & (smoke_ever == 0)] = "never"

    vig = pd.to_numeric(raw[f"r{w}vgactx"], errors="coerce").iloc[ix].to_numpy()
    mod = pd.to_numeric(raw[f"r{w}mdactx"], errors="coerce").iloc[ix].to_numpy()
    activity = np.full(len(ix), None, dtype=object)
    act_complete = np.isin(vig, [1, 2, 3, 4, 5]) & np.isin(mod, [1, 2, 3, 4, 5])
    best = np.fmin(vig, mod)
    activity[act_complete & (best <= 2)] = "frequent"
    activity[act_complete & (best > 2) & (best <= 4)] = "some"
    activity[act_complete & (best == 5)] = "inactive"

    confirm_observed = np.zeros(len(ix), dtype=bool)
    event_confirmed = np.full(len(ix), np.nan)
    lung_at_next2 = np.full(len(ix), np.nan)
    if w + 2 <= 16:
        later = binary(raw[f"r{w+2}lunge"]).iloc[ix].to_numpy()
        confirm_observed = raw[f"inw{w+2}"].iloc[ix].eq(1).to_numpy() & ~pd.isna(later)
        lung_at_next2[confirm_observed] = later[confirm_observed]
        event_confirmed[confirm_observed] = (
            (lung1.iloc[ix].to_numpy()[confirm_observed] == 1)
            & (later[confirm_observed] == 1)
        ).astype(int)

    wave_frames.append(pd.DataFrame({
        "row": ix,
        "hhidpn": raw["hhidpn"].iloc[ix].to_numpy(),
        "wave": w,
        "year": YEAR[w],
        "event": (lung1.iloc[ix].to_numpy() == 1).astype(int),
        "confirm_observed": confirm_observed,
        "event_confirmed": event_confirmed,
        "lung_at_next2": lung_at_next2,
        "pef": pef.iloc[ix].to_numpy(),
        "age": age.iloc[ix].to_numpy(),
        "sex": raw["ragender"].iloc[ix].to_numpy(),
        "height_m": height.iloc[ix].to_numpy(),
        "race": raw["raracem"].iloc[ix].to_numpy(),
        "survey_stratum": raw["raestrat"].iloc[ix].to_numpy(),
        "survey_half_sample": raw["raehsamp"].iloc[ix].to_numpy(),
        "education_years": raw["raedyrs"].iloc[ix].to_numpy(),
        "wealth": pd.to_numeric(raw[f"h{w}atotb"], errors="coerce").iloc[ix].to_numpy(),
        "smoking": smoking,
        "bmi": pd.to_numeric(raw[f"r{w}bmi"], errors="coerce").iloc[ix].to_numpy(),
        "activity": activity,
        "hypertension": binary(raw[f"r{w}hibpe"]).iloc[ix].to_numpy(),
        "diabetes": binary(raw[f"r{w}diabe"]).iloc[ix].to_numpy(),
        "heart": binary(raw[f"r{w}hearte"]).iloc[ix].to_numpy(),
        "stroke": binary(raw[f"r{w}stroke"]).iloc[ix].to_numpy(),
        "cancer": binary(raw[f"r{w}cancre"]).iloc[ix].to_numpy(),
        "cesd": pd.to_numeric(raw[f"r{w}cesd"], errors="coerce").iloc[ix].to_numpy(),
        "respondent_weight": pd.to_numeric(raw[f"r{w}wtresp"], errors="coerce").iloc[ix].to_numpy(),
    }))

df = pd.concat(wave_frames, ignore_index=True).sort_values(["row", "wave"]).drop_duplicates("row", keep="first")
df = df.drop(columns="row").reset_index(drop=True)
df.loc[~df["height_m"].between(1.2, 2.2), "height_m"] = np.nan
df.loc[~df["bmi"].between(12, 70), "bmi"] = np.nan
df.loc[~df["cesd"].between(0, 8), "cesd"] = np.nan
df.loc[~df["education_years"].between(0, 17), "education_years"] = np.nan
df.loc[~df["race"].isin([1, 2, 3]), "race"] = np.nan

# Wave-specific wealth quintiles limit distortion from nominal dollar changes.
df["wealth_q"] = np.nan
for w, g in df.groupby("wave"):
    ok = g["wealth"].notna()
    ranks = g.loc[ok, "wealth"].rank(method="average", pct=True)
    df.loc[ranks.index, "wealth_q"] = np.ceil(ranks * 5).clip(1, 5)

# Residualized PEF: sex-specific log(PEF) adjusted for age, age^2, height, wave.
df["pef_residual"] = np.nan
df["low_pef_adj"] = np.nan
df["pef_lower_1sd"] = np.nan
pef_definition = {}
for sex, g0 in df.groupby("sex"):
    g = g0.dropna(subset=["height_m"]).copy()
    wave_dummies = pd.get_dummies(g["wave"].astype(int), drop_first=True, dtype=float)
    x = np.column_stack([
        np.ones(len(g)), g["age"], g["age"] ** 2, g["height_m"], wave_dummies.to_numpy()
    ])
    y = np.log(g["pef"].to_numpy())
    beta, *_ = np.linalg.lstsq(x, y, rcond=None)
    residual = y - x @ beta
    cutoff = float(np.quantile(residual, .20))
    sd = float(np.std(residual, ddof=1))
    df.loc[g.index, "pef_residual"] = residual
    df.loc[g.index, "low_pef_adj"] = (residual <= cutoff).astype(int)
    df.loc[g.index, "pef_lower_1sd"] = -residual / sd
    pef_definition[str(int(sex))] = {"residual_q20": cutoff, "residual_sd": sd, "n": len(g)}

variables = [
    "event", "pef", "height_m", "race", "education_years", "wealth_q",
    "smoking", "bmi", "activity", "hypertension", "diabetes", "heart",
    "stroke", "cancer", "cesd", "low_pef_adj", "pef_lower_1sd",
]
missing = pd.DataFrame({
    "variable": variables,
    "nonmissing_n": [int(df[v].notna().sum()) for v in variables],
    "missing_n": [int(df[v].isna().sum()) for v in variables],
    "missing_pct": [float(df[v].isna().mean() * 100) for v in variables],
})

formulas = {
    "M0": "event ~ exposure",
    "M1": "event ~ exposure + age + I(age ** 2) + sex + height_m + C(wave)",
    "M2_core": (
        "event ~ exposure + age + I(age ** 2) + sex + height_m + C(wave) "
        "+ C(race) + education_years + C(smoking) + bmi"
    ),
    "M3_core": (
        "event ~ exposure + age + I(age ** 2) + sex + height_m + C(wave) "
        "+ C(race) + education_years + C(smoking) + bmi "
        "+ hypertension + diabetes + heart + stroke + cancer + cesd"
    ),
    "M4_expanded": (
        "event ~ exposure + age + I(age ** 2) + sex + height_m + C(wave) "
        "+ C(race) + education_years + C(smoking) + bmi "
        "+ hypertension + diabetes + heart + stroke + cancer + cesd "
        "+ C(wealth_q) + C(activity)"
    ),
}


def fit_models(data, outcome="event", tag="primary"):
    results = []
    for exposure in ["low_pef_adj", "pef_lower_1sd"]:
        for model, template in formulas.items():
            formula = template.replace("event", outcome, 1).replace("exposure", exposure)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                fit = smf.glm(
                    formula=formula,
                    data=data,
                    family=sm.families.Poisson(link=sm.families.links.Log()),
                ).fit(cov_type="HC0")
            b = float(fit.params[exposure])
            se = float(fit.bse[exposure])
            results.append({
                "analysis": tag,
                "exposure": exposure,
                "model": model,
                "n": int(fit.nobs),
                "events": int(np.asarray(fit.model.endog).sum()),
                "rr": float(np.exp(b)),
                "ci_low": float(np.exp(b - 1.96 * se)),
                "ci_high": float(np.exp(b + 1.96 * se)),
                "p": float(fit.pvalues[exposure]),
                "formula": formula,
            })
    return results


results = []
results += fit_models(df, tag="primary_all_waves")
results += fit_models(df[df["wave"] != 14], tag="exclude_2018_to_2020")
confirmed = df[df["confirm_observed"]].copy()
results += fit_models(confirmed, outcome="event_confirmed", tag="confirmed_at_next_interview")
delayed = df[df["confirm_observed"] & df["event"].eq(0)].copy()
delayed["delayed_event"] = delayed["lung_at_next2"]
results += fit_models(delayed, outcome="delayed_event", tag="four_year_delayed_incidence")
results = pd.DataFrame(results)

df.to_csv(OUT / "hrs_first_eligible_analysis_dataset.csv", index=False, encoding="utf-8-sig")
missing.to_csv(OUT / "missingness.csv", index=False, encoding="utf-8-sig")
results.to_csv(OUT / "modified_poisson_results.csv", index=False, encoding="utf-8-sig")

summary = {
    "risk_set_n": int(len(df)),
    "risk_set_events": int(df["event"].sum()),
    "pef_height_complete_n": int(df["low_pef_adj"].notna().sum()),
    "pef_height_complete_events": int(df.loc[df["low_pef_adj"].notna(), "event"].sum()),
    "full_model_complete_case_n": int(results.query("analysis == 'primary_all_waves' and model == 'M3_core'").iloc[0]["n"]),
    "full_model_complete_case_events": int(results.query("analysis == 'primary_all_waves' and model == 'M3_core'").iloc[0]["events"]),
    "delayed_risk_set_n": int(len(delayed)),
    "delayed_events": int(delayed["delayed_event"].sum()),
    "pef_definition": pef_definition,
    "weight_note": "RAND respondent weight is retained but not used; enhanced face-to-face physical-measure weight is absent locally.",
}
(OUT / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")

print(json.dumps(summary, ensure_ascii=False, indent=2))
print("\nMissingness:\n", missing.to_string(index=False))
print("\nModel results:\n", results.to_string(index=False))
