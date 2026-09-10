"""Build locked-baseline cohorts. Raw and participant data stay local."""
from pathlib import Path
import argparse, json
import numpy as np
import pandas as pd

ap=argparse.ArgumentParser()
ap.add_argument('--data-root',required=True)
ap.add_argument('--out',required=True)
a=ap.parse_args(); root=Path(a.data_root); out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
flows=[]; inventory=[]
def bin01(x): return pd.to_numeric(x,errors='coerce').where(lambda s:s.isin([0,1]))
def yn(x): return x.map({1:1.,5:0.})
def read(p,cols): return pd.read_stata(p,columns=list(dict.fromkeys(cols)),convert_categoricals=False)
def smoking(ev,no):
 return pd.Series(np.select([no.eq(1),no.eq(0)&ev.eq(1),no.eq(0)&ev.eq(0)],['current','former','never'],default='missing'),index=ev.index).replace('missing',np.nan)
def clean(d):
 d.loc[~d.bmi.between(12,70),'bmi']=np.nan
 d.loc[~d.cesd.between(0,8 if d.cohort.iloc[0]=='HRS' else 30),'cesd']=np.nan
 d.loc[~d.grip.between(0,100),'grip']=np.nan
 if 'education_years' in d:d.loc[~d.education_years.between(0,17),'education_years']=np.nan
 if 'race' in d:d.loc[~d.race.isin([1,2,3]),'race']=np.nan
 d['grip10']=d.grip/10
 d['cig10']=np.where(d.smoking.isin(['never','former']),0,d.cigs/10)
 d.loc[d.smoking.isna(),'cig10']=np.nan
 return d
def sequential(raw,conditions,cohort):
 mask=pd.Series(True,index=raw.index); before=len(raw)
 for label,cond in conditions:
  mask &= cond.fillna(False)
  n=int(mask.sum());flows.append(dict(cohort=cohort,stage=label,remaining=n,excluded=before-n));before=n
 return mask

hpath=next((root/'HRS data/01_rand_longitudinal').rglob('randhrs1992_2022v1.dta'))
stems=['agey_e','lunge','puff','pmhght','bmi','smokev','smoken','hibpe','diabe','hearte','stroke','cancre','cesd','wtresp','grp','doctor','hosp','iwstat']
cols=['hhidpn','ragender','raracem','raedyrs','raestrat','raehsamp','radyear']
for w in range(8,17):cols += [f'inw{w}',f'r{w}lunge',f'r{w}iwstat']
for w in range(8,15):cols += [f'r{w}{s}' for s in stems]
h=read(hpath,cols); frames=[]
for w in range(8,15):
 g=lambda s:h[f'r{w}{s}']
 d=pd.DataFrame({'id':h.hhidpn.astype('int64').astype(str),'wave':w,'year':1990+2*w,'cohort':'HRS','inw':h[f'inw{w}'],'age':g('agey_e'),'pef':g('puff'),'height_m':g('pmhght'),'lung0':bin01(g('lunge')),'sex':h.ragender,'race':h.raracem,'education_years':h.raedyrs,'stratum':h.raestrat,'psu':h.raehsamp,'weight':g('wtresp'),'bmi':g('bmi'),'cesd':g('cesd'),'smoking':smoking(bin01(g('smokev')),bin01(g('smoken'))),'grip':g('grp'),'doctor':bin01(g('doctor')),'hospital':bin01(g('hosp')),'cigs':np.nan,'breathless':np.nan,'cough':np.nan,'raw_available':0})
 for dest,src in [('hypertension','hibpe'),('diabetes','diabe'),('heart','hearte'),('stroke','stroke'),('cancer','cancre')]:d[dest]=bin01(g(src))
 for j in [1,2]:
  lung=bin01(h[f'r{w+j}lunge']);obs=h[f'inw{w+j}'].eq(1)&lung.notna()
  d[f'observed{j}']=obs.astype(int);d[f'event{j}']=lung.where(obs)
 d['death']=((d.observed1==0)&(h[f'r{w+1}iwstat'].isin([5,6])|h.radyear.between(1990+2*w,1992+2*w))).astype(int)
 frames.append(d)
pool=pd.concat(frames,ignore_index=True)
conds=[('Interviewed baseline wave',pool.inw.eq(1)),('Age 50 to 80',pool.age.between(50,80)),('PEF 30 to 900 L/min',pool.pef.between(30,900)),('Height 1.2 to 2.2 m',pool.height_m.between(1.2,2.2)),('Recorded sex',pool.sex.isin([1,2])),('No baseline reported lung disease',pool.lung0.eq(0)),('Positive weight and survey design',pool.weight.gt(0)&pool.stratum.notna()&pool.psu.isin([1,2]))]
mask=sequential(pool,conds,'HRS person-wave');hd=pool[mask].sort_values(['id','wave']).drop_duplicates('id').copy()
flows.append(dict(cohort='HRS',stage='First eligible baseline locked before follow-up',remaining=len(hd),excluded=int(mask.sum())-len(hd)))
# Supplemental raw HRS files are locally available in these waves only.
for year,w,prefix in [(2006,8,'k'),(2010,10,'m'),(2014,12,'o'),(2018,14,'q')]:
 fs=list((root/f'HRS data/02_rand_fat/{year}').rglob('*.dta'))
 if not fs:continue
 with pd.read_stata(fs[0],iterator=True) as rr:labs=rr.variable_labels()
 lookup={k.lower():k for k in labs}; names=['hhid','pn']+[prefix+s for s in ['c118','c119','c144','c149']]
 names=[lookup[n] for n in names if n in lookup];r=read(fs[0],names);r.columns=r.columns.str.lower()
 r['id']=(pd.to_numeric(r.hhid)*1000+pd.to_numeric(r.pn)).astype('int64').astype(str)
 r['cigs']=pd.to_numeric(r[prefix+'c118'],errors='coerce').where(lambda s:s.between(0,100))
 # Missing cigarettes can be supplied by a valid packs/day answer; DK/RF excluded.
 packs=pd.to_numeric(r.get(prefix+'c119',pd.Series(np.nan,index=r.index)),errors='coerce').where(lambda s:s.between(0,5))
 r['cigs']=r.cigs.fillna(packs*20)
 for dest,src in [('breathless','c144'),('cough','c149')]:r[dest]=yn(r.get(prefix+src,pd.Series(np.nan,index=r.index)))
 r=r.set_index('id');ix=hd.wave.eq(w)
 for v in ['cigs','breathless','cough']:hd.loc[ix,v]=hd.loc[ix,'id'].map(r[v])
 hd.loc[ix,'raw_available']=1
 for v in ['cigs','breathless','cough']:inventory.append(dict(cohort='HRS',wave=w,variable=v,eligible=int(ix.sum()),nonmissing=int(hd.loc[ix,v].notna().sum())))
hd=clean(hd);hd['event']=hd.event1;hd['observed']=hd.observed1
assert hd.id.is_unique
hd.to_csv(out/'HRS_baseline.csv',index=False)

cpath=next((root/'CHARLS date/Harmonized CHARLS/H_CHARLS_D_Data').rglob('*.dta'))
cs=['agey','puff','mheight','mbmi','lunge','asthmae','smokev','smoken','smokef','hibpe','diabe','hearte','stroke','cancre','cesd10','wtrespb','lgrip','rgrip','doctor1m','hosp1y','puff1','puff2','puff3','puffeff']
cols=['ID','communityID','ragender','raeduc_c','radyear','r4lunge','r4iwstat','inw4']
for w in [2,3]:cols += [f'inw{w}',f'h{w}rural']+[f'r{w}{s}' for s in cs]
c=read(cpath,cols)
def charls_frame(w):
 g=lambda s:c[f'r{w}{s}']; d=pd.DataFrame({'id':c.ID.astype(str).str.zfill(12),'cohort':'CHARLS','wave':w,'year':2011+2*(w-1),'inw':c[f'inw{w}'],'age':g('agey'),'pef':g('puff'),'height_m':g('mheight'),'lung0':bin01(g('lunge')),'asthma0':bin01(g('asthmae')),'sex':c.ragender,'education':c.raeduc_c,'rural':bin01(c[f'h{w}rural']),'psu':c.communityID,'stratum':1,'weight':g('wtrespb'),'bmi':g('mbmi'),'cesd':g('cesd10'),'smoking':smoking(bin01(g('smokev')),bin01(g('smoken'))),'grip':pd.concat([g('lgrip'),g('rgrip')],axis=1).max(axis=1),'doctor':bin01(g('doctor1m')),'hospital':bin01(g('hosp1y')),'cigs':g('smokef').where(lambda s:s.between(0,100)),'full_effort':g('puffeff').eq(1).astype(int)})
 for dest,src in [('hypertension','hibpe'),('diabetes','diabe'),('heart','hearte'),('stroke','stroke'),('cancer','cancre')]:d[dest]=bin01(g(src))
 vals=np.column_stack([g('puff1'),g('puff2'),g('puff3')]);vals=np.where((vals>=30)&(vals<=900),vals,np.nan)
 # Sort missing attempts below valid values so two valid attempts remain usable.
 sortedvals=np.sort(np.where(np.isfinite(vals),vals,-np.inf),axis=1)
 d['repeat40']=((np.isfinite(vals).sum(axis=1)>=2)&((sortedvals[:,-1]-sortedvals[:,-2])<=40)).astype(int)
 obs=c.inw4.eq(1)&bin01(c.r4lunge).notna();d['observed']=obs.astype(int);d['event']=bin01(c.r4lunge).where(obs)
 d['death']=((~obs)&(c.r4iwstat.isin([5,6])|c.radyear.between(2015,2018))).astype(int)
 d['mid_clear']=c.inw3.eq(1)&bin01(c.r3lunge).eq(0)
 return clean(d)
def cconditions(d):return [('Interviewed baseline wave',d.inw.eq(1)),('Age 50 to 80',d.age.between(50,80)),('PEF 30 to 900 L/min',d.pef.between(30,900)),('Height 1.2 to 2.2 m',d.height_m.between(1.2,2.2)),('Recorded sex',d.sex.isin([1,2])),('No reported lung disease or asthma',d.lung0.eq(0)&d.asthma0.eq(0)),('Positive weight and survey design',d.weight.gt(0)&d.psu.notna())]
cd=charls_frame(3);mask=sequential(cd,cconditions(cd),'CHARLS');cd=cd[mask].copy();assert cd.id.is_unique
cd.to_csv(out/'CHARLS_baseline.csv',index=False)
lag=charls_frame(2);mask=pd.Series(True,index=lag.index)
for _,m in cconditions(lag):mask &= m.fillna(False)
lag=lag[mask].copy();lag.to_csv(out/'CHARLS_2013_baseline.csv',index=False)
pd.DataFrame(flows).to_csv(out/'selection_flow.csv',index=False)
pd.DataFrame(inventory).to_csv(out/'raw_variable_coverage.csv',index=False)
for d in [cd,hd]:print(d.cohort.iloc[0], 'baseline',len(d),'observed',int(d.observed.sum()),'events',int(d.event.sum()),'deaths',int(d.death.sum()),flush=True)
print('BUILD COMPLETE',flush=True)
