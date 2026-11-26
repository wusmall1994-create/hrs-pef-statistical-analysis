"""Extract HRS physical-measure weights and peak-flow trial quality indicators."""
from pathlib import Path
import sys
import numpy as np
import pandas as pd

root = Path(sys.argv[1]) / "HRS data"
dest = Path(sys.argv[2])
tracker, = (root / "04_tracker").rglob("trk2022tr_r.dta")
prefixes = "klmnopq"
with pd.read_stata(tracker, iterator=True) as reader:
    columns = {name.lower(): name for name in reader.variable_labels()}
names = ["hhid", "pn"] + [p + "pmwgtr" for p in prefixes]
data = pd.read_stata(tracker, columns=[columns[k] for k in names], convert_categoricals=False)
data.columns = names
ids = data.hhid.astype(int) * 1000 + data.pn.astype(int)
weights = [pd.DataFrame(dict(id=ids, wave=i+8, pm=data[p+"pmwgtr"])) for i,p in enumerate(prefixes)]
pd.concat(weights).to_csv(dest/"pm.csv", index=False)
quality = []
for year, prefix in zip((2006,2010,2014,2018), "kmoq"):
    source, = (root/"02_rand_fat"/str(year)).rglob("*.dta")
    trials = [f"{prefix}i{i}" for i in range(807,810)]
    effort = prefix + "i810"
    data = pd.read_stata(source, columns=["hhid","pn",*trials,effort], convert_categoricals=False)
    blows = data[trials]
    ordered = np.sort(blows.where((blows >= 30) & (blows <= 900)).fillna(-9999), axis=1)
    quality.append(pd.DataFrame(dict(id=data.hhid.astype(int)*1000+data.pn.astype(int),
        wave=prefixes.index(prefix)+8, eff=data[effort],
        rep=(ordered[:,1]>0) & (ordered[:,2]-ordered[:,1]<=40))))
pd.concat(quality).to_csv(dest/"qual.csv", index=False)
