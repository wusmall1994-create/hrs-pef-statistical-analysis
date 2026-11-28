# Run after prediction.R: Rscript refit.R <RDS directory> <output directory> [workers]
args <- commandArgs(TRUE)
out <- args[2]
dir.create(out,showWarnings=FALSE,recursive=TRUE)
auc <- function(y,p,w) {
 i <- order(p); g <- cumsum(c(TRUE,diff(p[i])!=0))
 a <- rowsum(cbind(w[i]*y[i],w[i]*(1-y[i])),g,reorder=FALSE)
 sum(a[,1]*(cumsum(a[,2])-a[,2]/2))/prod(colSums(a))
}
one <- function(b) {
 set.seed(seed+b)
 mult <- rep(1,length(groups0))
 for (g in groups) if(length(g)>1) mult[g] <- tabulate(sample(g,length(g)-1,replace=TRUE),length(mult))[g]
 active <- which(mult>0)
 fold <- integer(length(mult));fold[active] <- sample(rep(1:5,length.out=length(active)))
 w <- d$weight*mult[cid]*rescale[cid]
 f <- fold[cid];p <- matrix(NA_real_,nrow(d),2)
 for(k in 1:5) {
 tr <- f>0 & f!=k; te <- f==k; z <- numeric(nrow(d))
 for (s in unique(d$sex)) {
 i <- which(tr & d$sex==s); j <- which(f>0 & d$sex==s)
 fit <- lm.wfit(ref[i,,drop=FALSE],log(d$pef[i]),mult[cid[i]])
 if(anyNA(fit$coefficients)) stop("reference rank")
 sd0 <- sqrt(sum(mult[cid[i]]*fit$residuals^2)/(sum(mult[cid[i]])-1))
 z[j] <- (log(d$pef[j])-ref[j,,drop=FALSE]%*%fit$coefficients)/sd0
 }
 for(j in 1:2) {
 xx <- if(j==1) x else cbind(x,z)
 fit <- suppressWarnings(glm.fit(xx[tr,,drop=FALSE],d$event[tr],weights=w[tr]/mean(w[tr]),family=quasibinomial(),control=glm.control(maxit=60)))
 if(!fit$converged || anyNA(fit$coefficients)) stop("logistic fit")
 p[te,j] <- plogis(xx[te,,drop=FALSE]%*%fit$coefficients)
 }
 }
 i <- f>0;y <- d$event[i];w <- w[i]
 c(diff(apply(p[i,,drop=FALSE],2,function(v)auc(y,v,w))),diff(apply(p[i,,drop=FALSE],2,function(v)weighted.mean((y-v)^2,w))))
}
cl <- parallel::makeCluster(if(length(args)>2) as.integer(args[3]) else 2)
for(co in c("CHARLS","HRS")){
 d <- readRDS(file.path(args[1],paste0(co,"Full clinical.rds")))$d
 seed <- if(co=="HRS") 700000 else 600000
 terms <- c("age","I(age^2)","sex_f","height_m","smoking_f","bmi","grip10","doctor","hospital","hypertension","diabetes","heart","stroke","cancer","cesd",if(co=="HRS")c("wave_f","race_f","education_years") else c("education_pred","rural_f"))
 x <- model.matrix(reformulate(terms),d)
 ref <- model.matrix(reformulate(c("age","I(age^2)","height_m",if(co=="HRS")"factor(wave)")),d)
 groups0 <- unique(as.character(d$cluster));cid <- match(d$cluster,groups0)
 groups <- split(seq_along(groups0),d$stratum[match(groups0,d$cluster)])
 rescale <- rep(1,length(groups0));for(g in groups) if(length(g)>1) rescale[g] <- length(g)/(length(g)-1)
 parallel::clusterExport(cl,c("d","x","ref","cid","groups","groups0","rescale","seed","auc","one"))
 v <- parallel::parLapplyLB(cl,1:499,function(b)tryCatch(one(b),error=function(e)stop(conditionMessage(e))))
 v <- as.data.frame(do.call(rbind,v)); names(v) <- c("delta_auc","delta_brier")
 write.csv(data.frame(draw=seq_len(nrow(v)),v),file.path(out,paste0(co,"_refit.csv")),row.names=FALSE)
 print(co);print(colSums(is.na(v)));print(sapply(v,quantile,c(.025,.975),na.rm=TRUE));flush.console()
}
parallel::stopCluster(cl)
