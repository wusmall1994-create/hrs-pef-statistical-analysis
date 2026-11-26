out <- commandArgs(TRUE)[1]
for (co in c("CHARLS","HRS")) {
  obj <- readRDS(file.path(out,paste0(co,"_analysis_local.rds")))
  stopifnot(!anyDuplicated(obj$all$id),all(is.finite(obj$all$z)),
            all(is.na(obj$all$event[obj$all$observed==0])))
  write.csv(obj$cc,file.path(out,paste0(co,"_validation_local.csv")),row.names=FALSE)
}
m<-read.csv(file.path(out,'mi_actual_counts.csv'));f<-read.csv(file.path(out,'analysis_flow.csv'))
for(co in unique(m$cohort)){z<-m[m$cohort==co,];ff<-f[f$cohort==co & f$stage=='observed',];stopifnot(length(unique(z$n))==1,all(z$n==ff$n),all(z$events==ff$events))}
s<-read.csv(file.path(out,'followup_states.csv'));tot<-aggregate(s$prob,s[,c('cohort','low')],sum);stopifnot(all(abs(tot$x-1)<1e-10))
r<-read.csv(file.path(out,'associations.csv'));stopifnot(all(r$lo<r$rr&r$rr<r$hi),all(r$n>=r$events),all(is.finite(r$rr)))
reps<-read.csv(file.path(out,'prediction_repeats.csv'));fold<-read.csv(file.path(out,'prediction_fold_audit.csv'))
stopifnot(nrow(reps)==60,all(reps$auc>=0&reps$auc<=1),all(reps$brier>=0&reps$brier<=1))
for(co in unique(fold$cohort))for(su in unique(fold$suite))for(rp in 1:5){z<-subset(fold,cohort==co & suite==su & repeat_id==rp);n<-subset(reps,cohort==co & suite==su & repeat_id==rp)$n[1];stopifnot(sum(z$test_n)==n,all(z$train_n+z$test_n==n))}
