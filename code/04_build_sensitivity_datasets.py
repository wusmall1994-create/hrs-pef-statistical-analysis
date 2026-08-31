from pathlib import Path
import json
import os
import numpy as np
import pandas as pd

ROOT = Path(os.environ.get("PEF_PROJECT_ROOT", ".")).resolve()
HRS_DATA = Path(os.environ["HRS_DATA_FILE"]).expanduser().resolve()
CHARLS_DATA = Path(os.environ["CHARLS_DATA_FILE"]).expanduser().resolve()
OUT = ROOT / "outputs" / "sensitivity"
OUT.mkdir(parents=True, exist_ok=True)


def binary(s):
    x = pd.to_numeric(s, errors="coerce")
    return x.where(x.isin([0, 1]))


def smoking_from(ever, now):
    ever = binary(ever)
    now = binary(now)
    return np.select(
        [now.eq(1), now.eq(0) & ever.eq(1), now.eq(0) & ever.eq(0)],
        ["current", "former", "never"], default=None,
    )


def hrs_datasets():
    waves = list(range(8, 15))
    fixed = ["hhidpn", "ragender", "raracem", "raedyrs", "raestrat", "raehsamp", "radyear"]
    cols = fixed[:]
    stems = [
        "agey_e", "lunge", "lunge2", "puff", "pmhght", "bmi", "smokev", "smoken",
        "vgactx", "mdactx", "hibpe", "diabe", "hearte", "stroke", "cancre", "cesd", "wtresp", "iwstat",
    ]
    for w in range(8, 17):
        cols += [f"inw{w}"]
        cols += [f"r{w}{s}" for s in ["lunge", "lunge2", "iwstat"]]
    for w in waves:
        cols += [f"r{w}{s}" for s in stems if f"r{w}{s}" not in cols]
    cols = list(dict.fromkeys(cols))
    raw = pd.read_stata(HRS_DATA, columns=cols, convert_categoricals=False)

    frames = []
    lagged_frames = []
    for w in waves:
        year = 1990 + 2 * w
        age = pd.to_numeric(raw[f"r{w}agey_e"], errors="coerce")
        pef = pd.to_numeric(raw[f"r{w}puff"], errors="coerce")
        height = pd.to_numeric(raw[f"r{w}pmhght"], errors="coerce")
        lung0 = binary(raw[f"r{w}lunge"])
        wt = pd.to_numeric(raw[f"r{w}wtresp"], errors="coerce")
        eligible = (
            raw[f"inw{w}"].eq(1) & age.between(50, 80) & pef.between(30, 900)
            & height.between(1.2, 2.2) & raw["ragender"].isin([1, 2]) & lung0.eq(0)
            & wt.gt(0) & raw["raestrat"].notna() & raw["raehsamp"].isin([1, 2])
        )
        ix = np.flatnonzero(eligible.to_numpy())
        next_lung = binary(raw[f"r{w+1}lunge"]).iloc[ix].to_numpy()
        next_inw = raw[f"inw{w+1}"].iloc[ix].eq(1).to_numpy()
        observed = next_inw & ~pd.isna(next_lung)
        iwstat = pd.to_numeric(raw[f"r{w+1}iwstat"], errors="coerce").iloc[ix].to_numpy()
        death = (~observed) & np.isin(iwstat, [5, 6])
        radyear = pd.to_numeric(raw["radyear"], errors="coerce").iloc[ix].to_numpy()
        death |= (~observed) & np.isfinite(radyear) & (radyear >= year) & (radyear <= year + 2)

        smoke = smoking_from(raw[f"r{w}smokev"].iloc[ix], raw[f"r{w}smoken"].iloc[ix])
        frame = pd.DataFrame({
            "hhidpn": raw["hhidpn"].iloc[ix].to_numpy(), "wave": w, "year": year,
            "followup_observed": observed.astype(int), "death_before_followup": death.astype(int),
            "followup_status": np.where(observed, "observed", np.where(death, "death", "other_nonresponse")),
            "event": np.where(observed, next_lung, np.nan),
            "pef": pef.iloc[ix].to_numpy(), "age": age.iloc[ix].to_numpy(),
            "sex": raw["ragender"].iloc[ix].to_numpy(), "height_m": height.iloc[ix].to_numpy(),
            "race": raw["raracem"].iloc[ix].to_numpy(), "survey_stratum": raw["raestrat"].iloc[ix].to_numpy(),
            "survey_half_sample": raw["raehsamp"].iloc[ix].to_numpy(),
            "education_years": raw["raedyrs"].iloc[ix].to_numpy(), "smoking": smoke,
            "bmi": pd.to_numeric(raw[f"r{w}bmi"], errors="coerce").iloc[ix].to_numpy(),
            "hypertension": binary(raw[f"r{w}hibpe"]).iloc[ix].to_numpy(),
            "diabetes": binary(raw[f"r{w}diabe"]).iloc[ix].to_numpy(),
            "heart": binary(raw[f"r{w}hearte"]).iloc[ix].to_numpy(),
            "stroke": binary(raw[f"r{w}stroke"]).iloc[ix].to_numpy(),
            "cancer": binary(raw[f"r{w}cancre"]).iloc[ix].to_numpy(),
            "cesd": pd.to_numeric(raw[f"r{w}cesd"], errors="coerce").iloc[ix].to_numpy(),
            "respondent_weight": wt.iloc[ix].to_numpy(),
        })
        frames.append(frame)

        if w + 2 <= 16:
            mid_lung = binary(raw[f"r{w+1}lunge"]).iloc[ix].to_numpy()
            late_lung = binary(raw[f"r{w+2}lunge"]).iloc[ix].to_numpy()
            mid_obs = raw[f"inw{w+1}"].iloc[ix].eq(1).to_numpy() & ~pd.isna(mid_lung)
            late_obs = raw[f"inw{w+2}"].iloc[ix].eq(1).to_numpy() & ~pd.isna(late_lung)
            lag_ok = mid_obs & late_obs & (mid_lung == 0)
            lag = frame.loc[lag_ok].copy()
            lag["event"] = late_lung[lag_ok]
            lag["followup_years"] = 4
            lagged_frames.append(lag)

    all_first = pd.concat(frames, ignore_index=True).sort_values(["hhidpn", "wave"]).drop_duplicates("hhidpn", keep="first")
    all_first.loc[~all_first["bmi"].between(12, 70), "bmi"] = np.nan
    all_first.loc[~all_first["cesd"].between(0, 8), "cesd"] = np.nan
    all_first.loc[~all_first["education_years"].between(0, 17), "education_years"] = np.nan
    all_first.loc[~all_first["race"].isin([1, 2, 3]), "race"] = np.nan
    all_first.to_csv(OUT / "hrs_first_baseline_all_followup_status.csv", index=False, encoding="utf-8-sig")

    lagged = pd.concat(lagged_frames, ignore_index=True).sort_values(["hhidpn", "wave"]).drop_duplicates("hhidpn", keep="first")
    lagged.loc[~lagged["bmi"].between(12, 70), "bmi"] = np.nan
    lagged.loc[~lagged["cesd"].between(0, 8), "cesd"] = np.nan
    lagged.loc[~lagged["education_years"].between(0, 17), "education_years"] = np.nan
    lagged.loc[~lagged["race"].isin([1, 2, 3]), "race"] = np.nan
    lagged.to_csv(OUT / "hrs_four_year_lagged.csv", index=False, encoding="utf-8-sig")
    return {
        "eligible_first_baseline": int(len(all_first)),
        "observed": int(all_first.followup_observed.sum()),
        "deaths": int(all_first.death_before_followup.sum()),
        "other_nonresponse": int((all_first.followup_status == "other_nonresponse").sum()),
        "lagged_n": int(len(lagged)), "lagged_events": int(lagged.event.sum()),
    }


def charls_datasets():
    fixed = ["ID", "communityID", "ragender", "raeduc_c", "radyear"]
    cols = fixed + ["inw2", "inw3", "inw4", "r3iwstat", "r4iwstat"]
    stems = [
        "agey", "lunge", "asthmae", "mheight", "mbmi", "puff", "puff1", "puff2", "puff3",
        "puffeff", "puffcomp", "smokev", "smoken", "hibpe", "diabe", "hearte", "stroke", "cancre",
        "cesd10", "wtrespb", "rxlung", "rxlung_c",
    ]
    for w in [2, 3]:
        cols += [f"r{w}{s}" for s in stems]
        cols += [f"h{w}rural"]
    cols += ["r4lunge"]
    cols = list(dict.fromkeys(cols))
    raw = pd.read_stata(CHARLS_DATA, columns=cols, convert_categoricals=False)

    age = pd.to_numeric(raw["r3agey"], errors="coerce")
    pef = pd.to_numeric(raw["r3puff"], errors="coerce")
    height = pd.to_numeric(raw["r3mheight"], errors="coerce")
    weight = pd.to_numeric(raw["r3wtrespb"], errors="coerce")
    eligible = (
        raw["inw3"].eq(1) & age.between(50, 80) & pef.between(30, 900) & height.between(1.2, 2.2)
        & raw["ragender"].isin([1, 2]) & binary(raw["r3lunge"]).eq(0) & binary(raw["r3asthmae"]).eq(0)
        & weight.gt(0) & raw["communityID"].notna()
    )
    ix = np.flatnonzero(eligible.to_numpy())
    next_lung = binary(raw["r4lunge"]).iloc[ix].to_numpy()
    observed = raw["inw4"].iloc[ix].eq(1).to_numpy() & ~pd.isna(next_lung)
    iwstat = pd.to_numeric(raw["r4iwstat"], errors="coerce").iloc[ix].to_numpy()
    death = (~observed) & np.isin(iwstat, [5, 6])
    radyear = pd.to_numeric(raw["radyear"], errors="coerce").iloc[ix].to_numpy()
    death |= (~observed) & np.isfinite(radyear) & (radyear >= 2015) & (radyear <= 2018)
    attempts = [pd.to_numeric(raw[f"r3puff{i}"], errors="coerce").iloc[ix].to_numpy() for i in [1, 2, 3]]
    valid_attempts = np.column_stack([np.where((a >= 30) & (a <= 900), a, np.nan) for a in attempts])
    sorted_attempts = np.sort(valid_attempts, axis=1)
    top1, top2 = sorted_attempts[:, 2], sorted_attempts[:, 1]
    attempt_n = np.sum(np.isfinite(valid_attempts), axis=1)
    repeat40 = (attempt_n >= 2) & np.isfinite(top1) & np.isfinite(top2) & ((top1 - top2) <= 40)
    repeat10pct = (attempt_n >= 2) & np.isfinite(top1) & np.isfinite(top2) & ((top1 - top2) <= 0.10 * top1)
    smoke = smoking_from(raw["r3smokev"].iloc[ix], raw["r3smoken"].iloc[ix])
    frame = pd.DataFrame({
        "ID": raw["ID"].iloc[ix].to_numpy(), "communityID": raw["communityID"].iloc[ix].to_numpy(),
        "followup_observed": observed.astype(int), "death_before_followup": death.astype(int),
        "followup_status": np.where(observed, "observed", np.where(death, "death", "other_nonresponse")),
        "event": np.where(observed, next_lung, np.nan), "age": age.iloc[ix].to_numpy(),
        "sex": raw["ragender"].iloc[ix].to_numpy(), "height_m": height.iloc[ix].to_numpy(),
        "pef": pef.iloc[ix].to_numpy(), "pef1": attempts[0], "pef2": attempts[1], "pef3": attempts[2],
        "valid_attempts_n": attempt_n, "repeatable_40": repeat40.astype(int), "repeatable_10pct": repeat10pct.astype(int),
        "full_effort": (pd.to_numeric(raw["r3puffeff"], errors="coerce").iloc[ix].to_numpy() == 1).astype(int),
        "lung_medication": binary(raw["r3rxlung_c"]).iloc[ix].to_numpy(),
        "race_context": binary(raw["h3rural"]).iloc[ix].to_numpy(), "education": raw["raeduc_c"].iloc[ix].to_numpy(),
        "smoking": smoke, "bmi": pd.to_numeric(raw["r3mbmi"], errors="coerce").iloc[ix].to_numpy(),
        "hypertension": binary(raw["r3hibpe"]).iloc[ix].to_numpy(), "diabetes": binary(raw["r3diabe"]).iloc[ix].to_numpy(),
        "heart": binary(raw["r3hearte"]).iloc[ix].to_numpy(), "stroke": binary(raw["r3stroke"]).iloc[ix].to_numpy(),
        "cancer": binary(raw["r3cancre"]).iloc[ix].to_numpy(),
        "cesd": pd.to_numeric(raw["r3cesd10"], errors="coerce").iloc[ix].to_numpy(), "weight": weight.iloc[ix].to_numpy(),
    })
    frame.loc[~frame["bmi"].between(12, 70), "bmi"] = np.nan
    frame.loc[~frame["cesd"].between(0, 30), "cesd"] = np.nan
    frame.to_csv(OUT / "charls_baseline_all_followup_status.csv", index=False, encoding="utf-8-sig")

    # Five-year delayed-incidence sensitivity: 2013 baseline, lung disease-free in 2013 and 2015,
    # outcome assessed in 2018.
    age2 = pd.to_numeric(raw["r2agey"], errors="coerce")
    pef2 = pd.to_numeric(raw["r2puff"], errors="coerce")
    height2 = pd.to_numeric(raw["r2mheight"], errors="coerce")
    wt2 = pd.to_numeric(raw["r2wtrespb"], errors="coerce")
    lag_ok = (
        raw["inw2"].eq(1) & raw["inw3"].eq(1) & raw["inw4"].eq(1)
        & age2.between(50, 80) & pef2.between(30, 900) & height2.between(1.2, 2.2)
        & raw["ragender"].isin([1, 2]) & binary(raw["r2lunge"]).eq(0) & binary(raw["r2asthmae"]).eq(0)
        & binary(raw["r3lunge"]).eq(0) & binary(raw["r4lunge"]).notna()
        & wt2.gt(0) & raw["communityID"].notna()
    )
    jx = np.flatnonzero(lag_ok.to_numpy())
    lag = pd.DataFrame({
        "ID": raw["ID"].iloc[jx].to_numpy(), "communityID": raw["communityID"].iloc[jx].to_numpy(),
        "event": binary(raw["r4lunge"]).iloc[jx].to_numpy(), "age": age2.iloc[jx].to_numpy(),
        "sex": raw["ragender"].iloc[jx].to_numpy(), "height_m": height2.iloc[jx].to_numpy(), "pef": pef2.iloc[jx].to_numpy(),
        "race_context": binary(raw["h2rural"]).iloc[jx].to_numpy(), "education": raw["raeduc_c"].iloc[jx].to_numpy(),
        "smoking": smoking_from(raw["r2smokev"].iloc[jx], raw["r2smoken"].iloc[jx]),
        "bmi": pd.to_numeric(raw["r2mbmi"], errors="coerce").iloc[jx].to_numpy(),
        "hypertension": binary(raw["r2hibpe"]).iloc[jx].to_numpy(), "diabetes": binary(raw["r2diabe"]).iloc[jx].to_numpy(),
        "heart": binary(raw["r2hearte"]).iloc[jx].to_numpy(), "stroke": binary(raw["r2stroke"]).iloc[jx].to_numpy(),
        "cancer": binary(raw["r2cancre"]).iloc[jx].to_numpy(),
        "cesd": pd.to_numeric(raw["r2cesd10"], errors="coerce").iloc[jx].to_numpy(), "weight": wt2.iloc[jx].to_numpy(),
        "followup_years": 5,
    })
    lag.loc[~lag["bmi"].between(12, 70), "bmi"] = np.nan
    lag.loc[~lag["cesd"].between(0, 30), "cesd"] = np.nan
    lag.to_csv(OUT / "charls_five_year_lagged.csv", index=False, encoding="utf-8-sig")
    return {
        "eligible_baseline": int(len(frame)), "observed": int(frame.followup_observed.sum()),
        "deaths": int(frame.death_before_followup.sum()),
        "other_nonresponse": int((frame.followup_status == "other_nonresponse").sum()),
        "full_effort_n": int(frame.full_effort.sum()), "repeatable_40_n": int(frame.repeatable_40.sum()),
        "repeatable_10pct_n": int(frame.repeatable_10pct.sum()),
        "baseline_lung_medication_yes": int(frame.lung_medication.eq(1).sum()),
        "lagged_n": int(len(lag)), "lagged_events": int(lag.event.sum()),
    }


summary = {"HRS": hrs_datasets(), "CHARLS": charls_datasets()}
(OUT / "dataset_build_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
print(json.dumps(summary, ensure_ascii=False, indent=2))
