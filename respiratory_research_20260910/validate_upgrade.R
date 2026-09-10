suppressPackageStartupMessages({library(survey);library(mice)})
options(survey.lonely.psu='adjust')
args<-commandArgs(TRUE);out<-args[1];archive<-args[2]
checks<-list()
# Independently refit the old analysis from the archived person-level datasets.
for(co in c('CHARLS','HRS')){
 d<-read.csv(file.path(archive,paste0(tolower(co),'_formal_dataset.csv')),check.names=FALSE)
 names(d)[1]<-if(co=='HRS')'hhidpn' else 'ID';d$sex_f<-factor(d$sex);d$smoking_f<-factor(d$smoking,levels=c('never','former','current'))
 common<-c('age','I(age^2)','sex_f','height_m','smoking_f','bmi','hypertension','diabetes','heart','stroke','cancer','cesd')
 if(co=='HRS'){d$race_f<-factor(d$race);d$wave_f<-factor(d$wave);cv<-c(common,'race_f','education_years','wave_f');des<-svydesign(ids=~survey_half_sample,strata=~survey_stratum,weights=~respondent_weight,nest=TRUE,data=d)}else{d$education_f<-factor(d$education);d$rural_f<-factor(d$race_context);cv<-c(common,'education_f','rural_f');des<-svydesign(ids=~communityID,weights=~weight,data=d)}
 fit<-svyglm(as.formula(paste('event~low_pef_adj+',paste(cv,collapse='+'))),des,family=quasipoisson('log'));b<-coef(fit)['low_pef_adj'];se<-sqrt(vcov(fit)['low_pef_adj','low_pef_adj']);mf<-model.frame(fit)
 checks[[co]]<-data.frame(cohort=co,n=nrow(mf),events=sum(model.response(mf)),rr=exp(b),lo=exp(b-1.96*se),hi=exp(b+1.96*se))
 stopifnot(nrow(mf)==if(co=='CHARLS')8087 else 17030,round(exp(b),2)==if(co=='CHARLS')1.63 else 3.15)
 obj<-readRDS(file.path(out,paste0(co,'_analysis_local.rds')))
 stopifnot(!anyDuplicated(obj$all$id),all(is.finite(obj$all$z)),all(obj$all$height_m>=1.2 & obj$all$height_m<=2.2),all(is.na(obj$all$event[obj$all$observed==0])))
 # Save only local rows for independent coefficient validation; never distribute these files.
 write.csv(obj$cc,file.path(out,paste0(co,'_validation_local.csv')),row.names=FALSE)
}
write.csv(do.call(rbind,checks),file.path(out,'archived_result_reproduction.csv'),row.names=FALSE)
m<-read.csv(file.path(out,'mi_actual_counts.csv'));f<-read.csv(file.path(out,'analysis_flow.csv'))
for(co in unique(m$cohort)){z<-m[m$cohort==co,];ff<-f[f$cohort==co & f$stage=='observed',];stopifnot(length(unique(z$n))==1,all(z$n==ff$n),all(z$events==ff$events))}
s<-read.csv(file.path(out,'followup_states.csv'));tot<-aggregate(s$prob,s[,c('cohort','low')],sum);stopifnot(all(abs(tot$x-1)<1e-10))
r<-read.csv(file.path(out,'associations.csv'));stopifnot(all(r$lo<r$rr&r$rr<r$hi),all(r$n>=r$events),all(is.finite(r$rr)))
reps<-read.csv(file.path(out,'prediction_repeats.csv'));fold<-read.csv(file.path(out,'prediction_fold_audit.csv'))
stopifnot(nrow(reps)==60,all(reps$auc>=0&reps$auc<=1),all(reps$brier>=0&reps$brier<=1))
for(co in unique(fold$cohort))for(su in unique(fold$suite))for(rp in 1:5){z<-subset(fold,cohort==co & suite==su & repeat_id==rp);n<-subset(reps,cohort==co & suite==su & repeat_id==rp)$n[1];stopifnot(sum(z$test_n)==n,all(z$train_n+z$test_n==n))}
writeLines(c('PASS archived estimates reproduced','PASS unique locked baselines and valid exposures','PASS imputation actual sample/event counts','PASS exhaustive follow-up probabilities','PASS association intervals and finite estimates','PASS paired grouped cross-validation fold totals'),file.path(out,'validation_checks.txt'))
cat('VALIDATION COMPLETE\n')
