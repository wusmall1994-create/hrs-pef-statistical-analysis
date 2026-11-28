library(survey)
options(survey.lonely.psu="adjust")
args <- commandArgs(TRUE)
input <- args[1]
out <- args[2]
intervals <- data.frame()
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE)[1])), "survey_models.R"))
rows <- list()
set.seed(17834)
for(co in c("CHARLS","HRS")){
 o=readRDS(file.path(out,paste0(co,"_analysis_local.rds"))); a=o$all
 print(table(o$cc$low));refs=ref_fit(a)
 for(s in names(refs)){v=refs[[s]];write.csv(data.frame(term=names(coef(v$fit)),coef=coef(v$fit),cut=as.numeric(v$cut),sd=v$sd),file.path(out,paste0(co,s,"ref.csv")),row.names=F)}
 wr=lapply(split(a,a$sex),function(g){m=lm(formula(refs[[1]]$fit),g,weights=weight);r=resid(m);i=order(r);list(fit=m,cut=r[i][which(cumsum(g$weight[i])>=.2*sum(g$weight))[1]],sd=sqrt(weighted.mean(r^2,g$weight)))})
 z=ref_apply(a,wr);rows[[co]] <- fitrr(z[z$observed==1,],co,"Weighted reference")$row
 if(co=="CHARLS"){b=prepare(readone("CHARLS_2013_baseline"));b=ref_apply(b,ref_fit(b));z=merge(o$observed,b[,c("id","low")],by="id",suffixes=c("",".early"));z$pattern=factor(z$low+2*z$low.early,levels=0:3);fit=svyglm(as.formula(paste("event~pattern+",paste(covars(co),collapse="+"))),des(complete_for(z,c("pattern",covars(co)))),family=quasipoisson("log"));counts <- table(model.frame(fit)$pattern); z_history <- z
 estimates <- exp(cbind(coef(fit),confint(fit)))[2:4,,drop=FALSE]
 history <- data.frame(pattern=c("Neither","2015 only","2013 only","Both"),
                       n=as.integer(counts),rr=c(1,estimates[,1]),
                       lower=c(NA,estimates[,2]),upper=c(NA,estimates[,3]))
 write.csv(history,file.path(out,"exposure_history.csv"),row.names=FALSE)}
 for(suite in c("Basic","Expanded","Full clinical")){
 x=readRDS(file.path(out,paste0(co,suite,".rds")));d=x$d;p=x$p;y=d$event;w=d$weight
 groups=split(unique(data.frame(cl=as.character(d$cluster),st=d$stratum))$cl,unique(data.frame(cl=as.character(d$cluster),st=d$stratum))$st)
 cl=unique(as.character(d$cluster));idx=match(as.character(d$cluster),cl)
 # Assign tie groups once for exact weighted concordance
 forms=lapply(1:10,function(k){v=matrix(p,nrow=nrow(d))[,k];i=order(v);list(i=i,g=cumsum(c(TRUE,diff(v[i])!=0)),sq=(y-v)^2)})
 score=function(wt){au=br=numeric(10);for(k in 1:10){f=forms[[k]];i=f$i;g=f$g;ca=as.vector(rowsum((wt*y)[i],g,reorder=F));cn=as.vector(rowsum((wt*(1-y))[i],g,reorder=F));au[k]=sum(ca*(cumsum(cn)-cn/2))/sum(ca)/sum(cn);br[k]=weighted.mean(f$sq,wt)};c(mean(diff(matrix(au,2))),mean(diff(matrix(br,2))))}
 reps=replicate(499,{mult=setNames(rep(0,length(cl)),cl);for(g in groups){n=length(g);mult[g]=if(n==1)1 else tabulate(sample.int(n,n-1,T),n)*n/(n-1)};score(w*mult[idx])})
 vals=c(score(w),apply(reps,1,quantile,.025),apply(reps,1,quantile,.975));intervals <- rbind(intervals, data.frame(cohort=co,model=suite,delta_auc=vals[1],delta_brier=vals[2],auc_lower=vals[3],brier_lower=vals[4],auc_upper=vals[5],brier_upper=vals[6],row.names=NULL))
 }}

write.csv(intervals,file.path(out,"paired_intervals.csv"),row.names=FALSE)

hrs <- readRDS(file.path(out,"HRS_analysis_local.rds"))$observed
physical <- merge(hrs,read.csv(file.path(input,"pm.csv")),by=c("id","wave"))
physical <- subset(physical,pm>0)
physical$respondent_weight <- physical$weight
physical$weight <- physical$pm
fit <- fitrr(physical,"HRS","Physical measure weights")
rows$physical <- fit$row
write.csv(marginal(fit),file.path(out,"physical_weight_risks.csv"),row.names=FALSE)
quality <- merge(hrs,read.csv(file.path(input,"qual.csv")),by=c("id","wave"))
rows$effort <- fitrr(subset(quality,eff==1),"HRS","Full effort")$row
rows$repeatability <- fitrr(subset(quality,rep=="True"),"HRS","Repeatability")$row
write.csv(do.call(rbind,rows),file.path(out,"measurement_sensitivity.csv"),row.names=FALSE)

# Direct contrasts within a common analysis sample.
history_sample <- complete_for(z_history, c("pattern", covars("CHARLS")))
history_sample$pattern <- relevel(history_sample$pattern, ref = "1")
history_fit <- svyglm(reformulate(c("pattern", covars("CHARLS")), "event"),
                      des(history_sample), family = quasipoisson("log"))
history_contrast <- exp(cbind(coef(history_fit), confint(history_fit)))["pattern3", ]
write.csv(as.data.frame(t(history_contrast)), file.path(out, "persistent_low_contrast.csv"),
          row.names = FALSE)

physical$weight <- physical$respondent_weight
respondent_fit <- fitrr(physical, "HRS", "Respondent weights in physical-measure sample")
write.csv(respondent_fit$row, file.path(out, "matched_respondent_weight.csv"), row.names = FALSE)
