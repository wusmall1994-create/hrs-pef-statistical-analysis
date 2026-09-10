options(survey.lonely.psu='adjust',contrasts=c('contr.treatment','contr.poly'))
suppressPackageStartupMessages({library(survey);library(mice)})
args<-commandArgs(TRUE);input<-args[1];out<-args[2];dir.create(out,recursive=TRUE,showWarnings=FALSE)
writeout<-function(x,n)write.csv(x,file.path(out,paste0(n,'.csv')),row.names=FALSE,na='')
readone<-function(n)read.csv(file.path(input,paste0(n,'.csv')),colClasses=c(id='character'),check.names=FALSE)
common<-c('age','I(age^2)','sex_f','height_m','smoking_f','bmi','hypertension','diabetes','heart','stroke','cancer','cesd')
covars<-function(cohort)c(common,if(cohort=='HRS')c('race_f','education_years','wave_f') else c('education_f','rural_f'))
ref_fit<-function(d){
 lapply(split(d,d$sex),function(g){f<-if(g$cohort[1]=='HRS')log(pef)~age+I(age^2)+height_m+factor(wave) else log(pef)~age+I(age^2)+height_m
 fit<-lm(f,data=g);r<-residuals(fit);list(fit=fit,cut=quantile(r,.2),sd=sd(r))})
}
ref_apply<-function(d,refs){
 d$residual<-d$z<-d$low<-NA_real_
 for(s in names(refs)){ix<-which(d$sex==as.numeric(s));rr<-log(d$pef[ix])-predict(refs[[s]]$fit,newdata=d[ix,]);d$residual[ix]<-rr;d$z[ix]<-rr/refs[[s]]$sd;d$low[ix]<-as.integer(rr<=refs[[s]]$cut)}
 d$lower_sd<--d$z;d
}
prepare<-function(d){d$sex_f<-factor(d$sex);d$smoking_f<-factor(d$smoking,levels=c('never','former','current'));d$wave_f<-factor(d$wave)
 if(d$cohort[1]=='HRS')d$race_f<-factor(d$race) else {d$education_f<-factor(d$education);d$rural_f<-factor(d$rural)}
 d$cluster<-interaction(d$stratum,d$psu,drop=TRUE);d
}
des<-function(d){if(d$cohort[1]=='HRS')svydesign(ids=~psu,strata=~stratum,weights=~weight,nest=TRUE,data=d) else svydesign(ids=~psu,weights=~weight,data=d)}
complete_for<-function(d,terms,y='event'){v<-all.vars(as.formula(paste(y,'~',paste(terms,collapse='+'))));d[complete.cases(d[,v,drop=FALSE]),,drop=FALSE]}
fitrr<-function(d,cohort,label,x='low',extra=character(),y='event',core=covars(cohort)){
 cc<-complete_for(d,c(x,core,extra),y);f<-as.formula(paste(y,'~',paste(c(x,core,extra),collapse='+')))
 fit<-svyglm(f,des(cc),family=quasipoisson('log'));b<-coef(fit)[x];se<-sqrt(vcov(fit)[x,x]);mf<-model.frame(fit)
 row<-data.frame(cohort=cohort,analysis=label,exposure=x,n=nrow(mf),events=sum(model.response(mf)),rr=exp(b),lo=exp(b-1.96*se),hi=exp(b+1.96*se),p=2*pnorm(-abs(b/se)))
 list(row=row,fit=fit,data=cc)
}
marginal<-function(obj){fit<-obj$fit;d<-obj$data;w<-d$weight/sum(d$weight);b<-coef(fit);V<-vcov(fit);risks<-grad<-list()
 for(i in 0:1){z<-d;z$low<-i;X<-model.matrix(delete.response(terms(fit)),z);mu<-exp(drop(X%*%b));risks[[i+1]]<-sum(w*mu);grad[[i+1]]<-colSums(X*(w*mu))}
 vv<-function(g)sqrt(drop(t(g)%*%V%*%g));se0<-vv(grad[[1]]);se1<-vv(grad[[2]]);rd<-risks[[2]]-risks[[1]];serd<-vv(grad[[2]]-grad[[1]])
 data.frame(cohort=d$cohort[1],n=nrow(d),events=sum(d$event),risk0=risks[[1]],risk0_lo=risks[[1]]-1.96*se0,risk0_hi=risks[[1]]+1.96*se0,risk1=risks[[2]],risk1_lo=risks[[2]]-1.96*se1,risk1_hi=risks[[2]]+1.96*se1,rd=rd,rd_lo=rd-1.96*serd,rd_hi=rd+1.96*serd)
}
rrs<-risks<-flows<-missing<-states<-subgroups<-splines<-spline_tests<-table1<-mi_counts<-data.frame();data_list<-list()
for(cohort in c('CHARLS','HRS')){
 cat('START',cohort,'\n');flush.console()
 all<-prepare(readone(paste0(cohort,'_baseline')));refs<-ref_fit(all);all<-ref_apply(all,refs);saveRDS(refs,file.path(out,paste0(cohort,'_reference.rds')))
 d<-all[all$observed==1,];cv<-covars(cohort);cc<-complete_for(d,c('low',cv));data_list[[cohort]]<-list(all=all,observed=d,cc=cc)
 for(stage in c('baseline','observed','complete_case')){z<-switch(stage,baseline=all,observed=d,complete_case=cc);flows<-rbind(flows,data.frame(cohort=cohort,stage=stage,n=nrow(z),events=sum(z$event,na.rm=TRUE)))}
 for(v in unique(c(all.vars(as.formula(paste('~',paste(cv,collapse='+')))),'grip10','doctor','hospital','cig10'))){missing<-rbind(missing,data.frame(cohort=cohort,variable=v,n=nrow(d),missing=sum(is.na(d[[v]]))))}
 for(x in c('low','lower_sd')){res<-fitrr(d,cohort,'Core model',x);rrs<-rbind(rrs,res$row);if(x=='low')risks<-rbind(risks,marginal(res))}
 for(nm in c('Grip','Healthcare contact','Grip and healthcare contact','Smoking intensity')){
  ex<-switch(nm,Grip='grip10','Healthcare contact'=c('doctor','hospital'),'Grip and healthcare contact'=c('grip10','doctor','hospital'),'Smoking intensity'='cig10')
  z<-complete_for(d,c('low',cv,ex));if(nm=='Smoking intensity'&&cohort=='HRS')z<-z[z$raw_available==1,]
  for(expanded in c(FALSE,TRUE)){r<-fitrr(z,cohort,paste(nm,if(expanded)'expanded' else 'matched core'),extra=if(expanded)ex else character());rrs<-rbind(rrs,r$row)}
 }
 if(cohort=='HRS'){
  z<-complete_for(d,c('low',cv,'breathless','cough'));writeout(data.frame(n=nrow(z),events=sum(z$event)), 'symptom_eligibility')
  if(nrow(z)>=500&&sum(z$event)>=100){for(expanded in c(FALSE,TRUE))rrs<-rbind(rrs,fitrr(z,cohort,paste('Symptoms',if(expanded)'expanded' else 'matched core'),extra=if(expanded)c('breathless','cough') else character())$row)}
 }
 # Frozen exposure definition; alternatives are fitted on same eligible population.
 for(p in c(.1,.25)){
  alt<-d;for(s in unique(all$sex)){cut<-quantile(all$residual[all$sex==s],p);alt$low[alt$sex==s]<-as.integer(alt$residual[alt$sex==s]<=cut)}
  rrs<-rbind(rrs,fitrr(alt,cohort,paste0('Residual lowest ',100*p,'%'))$row)
 }
 rrs<-rbind(rrs,fitrr(d[d$pef>=60&d$pef<=800,],cohort,'PEF range 60 to 800')$row)
 if(cohort=='CHARLS'){
  for(v in c('full_effort','repeat40'))rrs<-rbind(rrs,fitrr(d[d[[v]]==1,],cohort,v)$row)
  ref<-d;ref$low<-as.integer(ifelse(ref$sex==1,ref$pef<138.64,ref$pef<91.94));rrs<-rbind(rrs,fitrr(ref,cohort,'Published Ji threshold')$row)
  lall<-prepare(readone('CHARLS_2013_baseline'));lall<-ref_apply(lall,ref_fit(lall));lag<-lall[lall$observed==1 & lall$mid_clear %in% c(TRUE,'True','TRUE',1),];rrs<-rbind(rrs,fitrr(lag,cohort,'2013 PEF 2015 clear 2018 outcome')$row)
 }else{
  ref<-d[d$age>=65,];pred<-ifelse(ref$sex==1,-130.15+5.67*ref$height_m*100-5.66*ref$age,213.41+2.80*ref$height_m*100-4.91*ref$age);sdv<-ifelse(ref$sex==1,118.31,73.86);ref$low<-as.integer((ref$pef-pred)/sdv< -1.645);rrs<-rbind(rrs,fitrr(ref,cohort,'Published Donahue age 65 to 80')$row)
  lag<-all[all$observed1==1&all$observed2==1&all$event1==0,];lag$event<-lag$event2;rrs<-rbind(rrs,fitrr(lag,cohort,'Locked baseline middle clear later outcome')$row)
  confirm<-all[all$observed1==1&all$observed2==1,];confirm$event<-as.integer(confirm$event1==1&confirm$event2==1);rrs<-rbind(rrs,fitrr(confirm,cohort,'Two consecutive positive reports')$row)
 }
 # Follow-up states are exhaustive, with nonresponse never coded as non-diagnosis.
 all$state<-ifelse(all$observed==1,ifelse(all$event==1,'Diagnosis','No diagnosis'),ifelse(all$death==1,'Death','Other nonresponse'))
 for(g in 0:1){z<-all[all$low==g,];dd<-des(z);for(st in c('Diagnosis','No diagnosis','Death','Other nonresponse')){est<-svymean(as.formula(paste0('~I(as.numeric(state=="',st,'"))')),dd);states<-rbind(states,data.frame(cohort=cohort,low=g,state=st,n=sum(z$state==st),prob=as.numeric(est),se=as.numeric(SE(est))))}}
 comp<-all[all$observed==1|all$death==1,];comp$event<-as.integer(comp$death==1|(!is.na(comp$event)&comp$event==1));rrs<-rbind(rrs,fitrr(comp,cohort,'Diagnosis or death observed status')$row)
 # IPCW on same locked baseline, covariate-complete population.
 ac<-complete_for(all,c('low',cv),'observed');obfit<-svyglm(as.formula(paste('observed~',paste(c('low',cv),collapse='+'))),des(ac),family=quasibinomial());pr<-as.numeric(predict(obfit,type='response'));iw<-weighted.mean(ac$observed,ac$weight)/pmax(pr,.01);limits<-quantile(iw,c(.01,.99));iw<-pmin(pmax(iw,limits[1]),limits[2]);ac$weight<-ac$weight*iw;rrs<-rbind(rrs,fitrr(ac[ac$observed==1,],cohort,'IPCW')$row)
 writeout(data.frame(cohort=cohort,minimum_probability=min(pr),maximum_weight=max(iw),trim_low=limits[1],trim_high=limits[2]),paste0(cohort,'_ipcw_diagnostics'))
 # Subgroup associations with a joint interaction test.
 cc$age_group<-factor(ifelse(cc$age<65,'50-64','65-80'))
 for(mod in c('sex_f','age_group','smoking_f')){
  intfit<-svyglm(as.formula(paste('event~',paste(unique(c(cv,mod)),collapse='+'),'+low+low:',mod)),des(cc),family=quasipoisson('log'))
  ip<-regTermTest(intfit,as.formula(paste('~low:',mod)),method='Wald')$p
  for(lev in levels(cc[[mod]])){z<-cc[cc[[mod]]==lev,];core<-setdiff(cv,mod);r<-fitrr(z,cohort,paste(mod,lev),core=core)$row;r$interaction_p<-ip;subgroups<-rbind(subgroups,r)}
 }
 # Natural cubic spline, equivalent restricted cubic shape with four knots.
 ks<-as.numeric(svyquantile(~z,des(cc),quantiles=c(.05,.35,.65,.95),ci=FALSE)[[1]])
 bas<-splines::ns(cc$z,knots=ks[2:3],Boundary.knots=ks[c(1,4)]);cc$s1<-bas[,1];cc$s2<-bas[,2];cc$s3<-bas[,3]
 sf<-svyglm(as.formula(paste('event~s1+s2+s3+',paste(cv,collapse='+'))),des(cc),family=quasipoisson('log'));pv<-regTermTest(sf,~s1+s2+s3)$p
 xs<-seq(quantile(cc$z,.01),quantile(cc$z,.99),length.out=150);bb<-predict(bas,xs);b0<-predict(bas,0);delta<-sweep(bb,2,b0[1,],'-');bet<-coef(sf)[c('s1','s2','s3')];VV<-vcov(sf)[c('s1','s2','s3'),c('s1','s2','s3')];eta<-drop(delta%*%bet);se<-sqrt(rowSums((delta%*%VV)*delta))
 splines<-rbind(splines,data.frame(cohort=cohort,z=xs,rr=exp(eta),lo=exp(eta-1.96*se),hi=exp(eta+1.96*se)));spline_tests<-rbind(spline_tests,data.frame(cohort=cohort,k1=ks[1],k2=ks[2],k3=ks[3],k4=ks[4],overall_p=pv))
 # Weighted descriptive values and analysis denominator.
 for(g in c('Overall','Non-low','Low')){z<-if(g=='Overall')d else d[d$low==as.integer(g=='Low'),];dd<-des(z)
  for(v in c('age','pef','height_m','bmi','grip')){m<-as.numeric(svymean(as.formula(paste0('~',v)),dd,na.rm=TRUE));sdv<-sqrt(as.numeric(svyvar(as.formula(paste0('~',v)),dd,na.rm=TRUE)));table1<-rbind(table1,data.frame(cohort=cohort,group=g,variable=v,value=m,sd=sdv,n=nrow(z)))}
  for(v in c('sex==2','smoking=="current"','doctor==1','hospital==1')){m<-as.numeric(svymean(as.formula(paste0('~I(as.numeric(',v,'))')),dd,na.rm=TRUE));table1<-rbind(table1,data.frame(cohort=cohort,group=g,variable=v,value=m,sd=NA,n=nrow(z)))}
 }
 saveRDS(data_list[[cohort]],file.path(out,paste0(cohort,'_analysis_local.rds')))
 cat('MI START',cohort,'\n');flush.console()
 vars<-unique(c('event','low','lower_sd',all.vars(as.formula(paste('~',paste(cv,collapse='+')))),'weight','psu','stratum'))
 dat<-d[,vars];methods<-make.method(dat);protected<-c('event','low','lower_sd','age','sex_f','height_m','wave_f','weight','psu','stratum');methods[intersect(protected,names(methods))]<-''
 pred<-make.predictorMatrix(dat);pred[,intersect(c('weight','psu','stratum'),colnames(pred))]<-0;pred[intersect(protected,rownames(pred)),]<-0
 imp<-mice(dat,m=20,maxit=10,method=methods,predictorMatrix=pred,printFlag=FALSE,seed=20260909+match(cohort,c('CHARLS','HRS')))
 saveRDS(imp,file.path(out,paste0(cohort,'_imputation_local.rds')))
 if(!is.null(imp$loggedEvents))writeout(imp$loggedEvents,paste0(cohort,'_mi_logged_events'))
 for(x in c('low','lower_sd')){qs<-us<-numeric(20)
  for(i in 1:20){z<-complete(imp,i);z$cohort<-cohort;r<-fitrr(z,cohort,'MI20',x);qs[i]<-log(r$row$rr);us[i]<-vcov(r$fit)[x,x];mi_counts<-rbind(mi_counts,data.frame(cohort=cohort,exposure=x,imputation=i,n=r$row$n,events=r$row$events))}
  total<-mean(us)+(1+1/20)*var(qs);se<-sqrt(total);df<-if(var(qs)>0)19*(1+mean(us)/((1+1/20)*var(qs)))^2 else Inf;crit<-qt(.975,df);b<-mean(qs)
  rrs<-rbind(rrs,data.frame(cohort=cohort,analysis='MI20',exposure=x,n=nrow(dat),events=sum(dat$event),rr=exp(b),lo=exp(b-crit*se),hi=exp(b+crit*se),p=2*pt(-abs(b/se),df)))
 }
 cat('DONE',cohort,'\n');flush.console()
 writeout(rrs,'associations');writeout(risks,'absolute_risks');writeout(flows,'analysis_flow');writeout(missing,'missingness');writeout(states,'followup_states');writeout(subgroups,'subgroups');writeout(splines,'spline_curves');writeout(spline_tests,'spline_tests');writeout(table1,'table1');writeout(mi_counts,'mi_actual_counts')
}
# QBA retains all scenarios, including incompatible corrected probabilities.
qba<-data.frame()
for(cohort in c('CHARLS','HRS')){r<-risks[risks$cohort==cohort,];grid<-expand.grid(sensitivity=c(.1,.13,.2,.4),specificity=c(.95,.989,.995,.999))
 for(i in 1:nrow(grid)){se<-grid$sensitivity[i];sp<-grid$specificity[i];p0<-(r$risk0+sp-1)/(se+sp-1);p1<-(r$risk1+sp-1)/(se+sp-1);ok<-p0>=0&&p0<=1&&p1>=0&&p1<=1;qba<-rbind(qba,data.frame(cohort=cohort,sensitivity=se,specificity=sp,admissible=ok,corrected_risk0=if(ok)p0 else NA,corrected_risk1=if(ok)p1 else NA,corrected_rr=if(ok)p1/p0 else NA))}}
writeout(qba,'misclassification_scenarios');writeLines(capture.output(sessionInfo()),file.path(out,'R_session.txt'))
cat('ANALYSIS COMPLETE\n')
