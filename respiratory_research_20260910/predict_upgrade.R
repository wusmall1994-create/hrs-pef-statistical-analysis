options(survey.lonely.psu='adjust')
args<-commandArgs(TRUE);out<-args[1]
dir.create(out,recursive=TRUE,showWarnings=FALSE)
ref_fit<-function(d){lapply(split(d,d$sex),function(g){f<-if(g$cohort[1]=='HRS')log(pef)~age+I(age^2)+height_m+factor(wave) else log(pef)~age+I(age^2)+height_m;fit<-lm(f,data=g);list(fit=fit,sd=sd(residuals(fit)))})}
ref_apply<-function(d,refs){d$z<-NA_real_;for(s in names(refs)){ix<-which(d$sex==as.numeric(s));d$z[ix]<-(log(d$pef[ix])-predict(refs[[s]]$fit,newdata=d[ix,]))/refs[[s]]$sd};d}
wauc<-function(y,p,w){a<-aggregate(cbind(cases=w*y,controls=w*(1-y)),list(p=p),sum);a<-a[order(a$p),];sum(a$cases*(cumsum(a$controls)-a$controls/2))/(sum(a$cases)*sum(a$controls))}
metrics<-function(y,p,w){p<-pmin(pmax(p,1e-8),1-1e-8);lp<-qlogis(p);ww<-w/mean(w);cal<-glm(y~lp,weights=ww,family=quasibinomial());citl<-glm(y~1+offset(lp),weights=ww,family=quasibinomial());c(auc=wauc(y,p,w),brier=weighted.mean((y-p)^2,w),calibration_intercept=coef(citl)[1],calibration_slope=coef(cal)[2])}
res<-calpoints<-coefs<-data.frame();fold_audit<-data.frame()
set.seed(20260909)
for(cohort in c('CHARLS','HRS')){
 obj<-readRDS(file.path(out,paste0(cohort,'_analysis_local.rds')));d0<-obj$observed
 if(cohort=='CHARLS')d0$education_pred<-cut(d0$education,breaks=c(-Inf,3,4,5,Inf),labels=c('Below primary','Primary','Middle','High or above'))
 for(suite in c('Basic','Expanded','Full clinical')){
  terms<-c('age','I(age^2)','sex_f','height_m','smoking_f','bmi',if(cohort=='HRS')'wave_f',if(suite!='Basic')c('grip10','doctor','hospital'),if(suite=='Full clinical')c('hypertension','diabetes','heart','stroke','cancer','cesd',if(cohort=='HRS')c('race_f','education_years') else c('education_pred','rural_f')))
  vv<-unique(c('event','pef','sex','wave',all.vars(as.formula(paste('~',paste(terms,collapse='+'))))))
  d<-d0[complete.cases(d0[,vv]),];d<-droplevels(d);preds<-array(NA_real_,c(nrow(d),2,5));clusters<-unique(as.character(d$cluster))
  for(rep in 1:5){map<-setNames(sample(rep(1:5,length.out=length(clusters))),clusters);fold<-unname(map[as.character(d$cluster)])
   for(k in 1:5){tr<-d[fold!=k,];te<-d[fold==k,];stopifnot(length(intersect(unique(tr$cluster),unique(te$cluster)))==0)
    ref<-ref_fit(tr);tr<-ref_apply(tr,ref);te<-ref_apply(te,ref)
    for(j in 1:2){f<-as.formula(paste('event~',paste(c(terms,if(j==2)'z'),collapse='+')));fit<-glm(f,data=tr,weights=weight/mean(weight),family=quasibinomial(),control=glm.control(maxit=50));stopifnot(fit$converged,!anyNA(coef(fit)));preds[fold==k,j,rep]<-predict(fit,newdata=te,type='response')}
    fold_audit<-rbind(fold_audit,data.frame(cohort=cohort,suite=suite,repeat_id=rep,fold=k,train_n=nrow(tr),train_events=sum(tr$event),test_n=nrow(te),test_events=sum(te$event),train_psu=length(unique(tr$cluster)),test_psu=length(unique(te$cluster))))
   }
   for(j in 1:2){m<-metrics(d$event,preds[,j,rep],d$weight);res<-rbind(res,data.frame(cohort=cohort,suite=suite,repeat_id=rep,model=if(j==1)'Without PEF' else 'With PEF',n=nrow(d),events=sum(d$event),auc=unname(m[1]),brier=unname(m[2]),calibration_intercept=unname(m[3]),calibration_slope=unname(m[4])))}
   cat('CV',cohort,suite,rep,'done\n');flush.console()
  }
  stopifnot(all(is.finite(preds)))
  # Calibration points use averaged held-out predictions, not in-sample fits.
  for(j in 1:2){p<-rowMeans(preds[,j,]);ord<-order(p);grp<-integer(length(p));grp[ord]<-pmin(10,ceiling(seq_along(ord)/length(ord)*10))
   for(g in 1:10){ix<-grp==g;calpoints<-rbind(calpoints,data.frame(cohort=cohort,suite=suite,model=if(j==1)'Without PEF' else 'With PEF',bin=g,n=sum(ix),events=sum(d$event[ix]),predicted=weighted.mean(p[ix],d$weight[ix]),observed=weighted.mean(d$event[ix],d$weight[ix])))}
  }
  # Full-data fit and reference transformations are local reusable model objects.
  ref<-ref_fit(d);z<-ref_apply(d,ref)
  for(j in 1:2){fit<-glm(as.formula(paste('event~',paste(c(terms,if(j==2)'z'),collapse='+'))),data=z,weights=weight/mean(weight),family=quasibinomial());coefs<-rbind(coefs,data.frame(cohort=cohort,suite=suite,model=if(j==1)'Without PEF' else 'With PEF',term=names(coef(fit)),coefficient=unname(coef(fit))))}
  refrows<-do.call(rbind,lapply(names(ref),function(s)data.frame(cohort=cohort,suite=suite,sex=s,term=names(coef(ref[[s]]$fit)),coefficient=unname(coef(ref[[s]]$fit)),residual_sd=ref[[s]]$sd)))
  write.csv(refrows,file.path(out,paste0('prediction_reference_',cohort,'_',suite,'.csv')),row.names=FALSE)
 }
}
write.csv(res,file.path(out,'prediction_repeats.csv'),row.names=FALSE)
summary<-aggregate(res[,c('auc','brier','calibration_intercept','calibration_slope')],res[,c('cohort','suite','model','n','events')],mean)
write.csv(summary,file.path(out,'prediction_summary.csv'),row.names=FALSE)
write.csv(calpoints,file.path(out,'calibration_points.csv'),row.names=FALSE)
write.csv(coefs,file.path(out,'prediction_coefficients.csv'),row.names=FALSE)
write.csv(fold_audit,file.path(out,'prediction_fold_audit.csv'),row.names=FALSE)
cat('PREDICTION COMPLETE\n')
