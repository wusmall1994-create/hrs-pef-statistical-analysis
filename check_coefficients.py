from pathlib import Path
import numpy as np,pandas as pd
import sys
res=Path(sys.argv[1]);rows=[]
expected=pd.read_csv(res/'associations.csv')
for co in ['CHARLS','HRS']:
 d=pd.read_csv(res/f'{co}_validation_local.csv')
 continuous=['low','age','height_m','bmi','hypertension','diabetes','heart','stroke','cancer','cesd']
 cats=['sex_f','smoking_f']+(['race_f','wave_f'] if co=='HRS' else ['education_f','rural_f'])
 if co=='HRS':continuous.append('education_years')
 X=d[continuous].copy();X['age2']=d.age**2
 for c in cats:X=pd.concat([X,pd.get_dummies(d[c],prefix=c,drop_first=True,dtype=float)],axis=1)
 # Centring/scaling improves numerical conditioning without changing the model space.
 mean=X.mean();sd=X.std();Z=(X-mean)/sd;Z.insert(0,'intercept',1.)
 xx=Z.to_numpy(float);y=d.event.to_numpy(float);w=d.weight.to_numpy(float);w=w/w.mean()
 beta=np.zeros(xx.shape[1]);beta[0]=np.log(np.average(y,weights=w))
 for iteration in range(300):
  eta=xx@beta;mu=np.exp(np.clip(eta,-30,50));target=eta+(y-mu)/mu;ww=np.sqrt(w*mu)
  new=np.linalg.lstsq(xx*ww[:,None],target*ww,rcond=None)[0]
  old_obj=np.sum(w*(mu-y*eta))
  for shrink in range(40):
   proposed=xx@new
   if np.max(proposed)<100 and np.sum(w*(np.exp(proposed)-y*proposed))<=old_obj+1e-10:break
   new=(new+beta)/2
  new_obj=np.sum(w*(np.exp(np.clip(xx@new,-30,50))-y*(xx@new)))
  low_index=list(Z.columns).index('low')
  if np.max(np.abs(new-beta))<1e-10 or (abs(new_obj-old_obj)<1e-9 and abs(new[low_index]-beta[low_index])<1e-9):beta=new;break
  beta=new
 ratio=np.exp(beta[list(Z.columns).index('low')]/sd['low'])
 rr=expected.query('cohort==@co and analysis=="Core model" and exposure=="low"').rr.iloc[0]
 assert abs(ratio-rr)<1e-7,(ratio,rr)
 rows.append({'cohort':co,'python_weighted_poisson_rr':ratio,'r_survey_rr':rr,'absolute_difference':abs(ratio-rr),'iterations':iteration+1})
pd.DataFrame(rows).to_csv(res/'independent_coefficient_check.csv',index=False)
print('PASS independent NumPy IRLS reproduces both primary coefficients')
