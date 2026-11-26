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
