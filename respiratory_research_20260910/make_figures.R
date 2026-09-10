suppressPackageStartupMessages({library(ggplot2);library(grid)})
args<-commandArgs(TRUE);input<-args[1];dest<-args[2];dir.create(dest,recursive=TRUE,showWarnings=FALSE)
rd<-function(n)read.csv(file.path(input,paste0(n,'.csv')))
theme_set(theme_classic(base_size=9,base_family='Arial')+theme(strip.background=element_blank(),strip.text=element_text(face='bold'),legend.position='bottom',axis.title=element_text(size=9),plot.margin=margin(8,9,8,8)))
cols<-c(CHARLS='#236B8E',HRS='#B76527')
saveplot<-function(p,name,w=180,h=125){
 ggsave(file.path(dest,paste0(name,'.pdf')),p,width=w,height=h,units='mm',device=cairo_pdf)
 ggsave(file.path(dest,paste0(name,'.tiff')),p,width=w,height=h,units='mm',dpi=600,compression='lzw')
 ggsave(file.path(dest,paste0(name,'.png')),p,width=w,height=h,units='mm',dpi=160)
}
flow<-rd('analysis_flow');states<-rd('followup_states')
flowdraw<-function(){grid.newpage()
 for(i in 1:2){co<-c('CHARLS','HRS')[i];f<-flow[flow$cohort==co,];cx<-if(i==1).25 else .75;color<-cols[co]
  grid.text(co,x=cx,y=.95,gp=gpar(fontsize=13,fontface='bold',col=color))
  n<-function(s)f$n[f$stage==s];ev<-f$events[f$stage=='observed'];st<-states[states$cohort==co,];dead<-sum(st$n[st$state=='Death']);other<-sum(st$n[st$state=='Other nonresponse'])
  labels<-c(paste0('Eligible baseline\nn = ',format(n('baseline'),big.mark=','),'\n',if(co=='CHARLS')'2015' else 'First eligible visit 2006 to 2018'),paste0('Observed primary outcome\nn = ',format(n('observed'),big.mark=','),'\n',ev,' diagnosis reports'),paste0('Core complete-case analysis\nn = ',format(n('complete_case'),big.mark=','),'\n',f$events[f$stage=='complete_case'],' diagnosis reports'))
  ys<-c(.78,.48,.18)
  for(j in 1:3){grid.roundrect(x=cx,y=ys[j],width=.41,height=.18,r=unit(.025,'snpc'),gp=gpar(fill='white',col=color,lwd=1));grid.text(labels[j],x=cx,y=ys[j],gp=gpar(fontsize=10,lineheight=1.25))}
  for(j in 1:2)grid.lines(x=c(cx,cx),y=c(ys[j]-.095,ys[j+1]+.095),arrow=arrow(length=unit(2,'mm')),gp=gpar(col=color))
  grid.text(paste0('Death: ',dead,'; other nonresponse: ',other),x=cx,y=.63,gp=gpar(fontsize=8),just='centre')
  grid.text(paste0('Missing core covariates: ',n('observed')-n('complete_case')),x=cx,y=.33,gp=gpar(fontsize=8))
 }
}
for(fmt in c('pdf','tiff','png')){
 p<-file.path(dest,paste0('Figure1_selection.',fmt))
 if(fmt=='pdf')cairo_pdf(p,width=180/25.4,height=135/25.4,family='Arial')
 if(fmt=='tiff')tiff(p,width=180,height=135,units='mm',res=600,compression='lzw',type='cairo',family='Arial')
 if(fmt=='png')png(p,width=180,height=135,units='mm',res=160,type='cairo',family='Arial')
 flowdraw();dev.off()
}
a<-rd('associations');z<-a[a$exposure=='low' & a$analysis %in% c('Core model','Grip and healthcare contact matched core','Grip and healthcare contact expanded'),]
lab<-c('Core model'='Core model\nAll core-complete participants','Grip and healthcare contact matched core'='Core model\nSame sample as expanded model','Grip and healthcare contact expanded'='Core + grip + healthcare contact\nSame-sample expanded model')
z$label<-factor(lab[z$analysis],levels=rev(unname(lab)));z$text<-sprintf('%.2f (%.2f–%.2f)',z$rr,z$lo,z$hi)
p<-ggplot(z,aes(rr,label,color=cohort))+geom_vline(xintercept=1,linetype=2,color='grey55')+geom_errorbar(aes(xmin=lo,xmax=hi),orientation='y',width=.12,linewidth=.55)+geom_point(size=2.5)+facet_wrap(~cohort,nrow=1)+scale_color_manual(values=cols,guide='none')+scale_x_log10()+labs(x='Risk ratio (95% CI)',y=NULL)
saveplot(p,'Figure2_expanded_associations',180,105)
s<-rd('spline_curves');p<-ggplot(s,aes(z,rr,color=cohort,fill=cohort))+geom_ribbon(aes(ymin=lo,ymax=hi),alpha=.15,color=NA)+geom_line(linewidth=.75)+geom_hline(yintercept=1,linetype=2,color='grey55')+facet_wrap(~cohort,scales='free_y')+scale_color_manual(values=cols,guide='none')+scale_fill_manual(values=cols,guide='none')+scale_y_log10()+labs(x='Standardised PEF residual',y='Risk ratio (95% CI)')
saveplot(p,'Figure3_dose_response',180,110)
cal<-subset(rd('calibration_points'),suite=='Full clinical');p<-ggplot(cal,aes(predicted,observed,color=model,shape=model))+geom_abline(slope=1,intercept=0,linetype=2,color='grey55')+geom_line(linewidth=.45)+geom_point(size=1.8)+facet_wrap(~cohort,scales='free',nrow=1)+scale_color_manual(values=c('Without PEF'='#777777','With PEF'='#236B8E'))+scale_x_continuous(labels=function(x)paste0(round(100*x),'%'),expand=expansion(mult=c(.02,.1)))+scale_y_continuous(labels=function(x)paste0(round(100*x),'%'),expand=expansion(mult=c(.02,.1)))+labs(x='Mean held-out predicted risk',y='Observed risk',color=NULL,shape=NULL)
saveplot(p,'Figure4_predictive_calibration',180,110)
z<-a[a$exposure=='low' & !(grepl('matched core|expanded',a$analysis)),];z$analysis[z$analysis=='full_effort']<-'Full recorded effort';z$analysis[z$analysis=='repeat40']<-'Two best attempts within 40 L/min';z$analysis[z$analysis=='MI20']<-'Multiple imputation';z$analysis<-factor(z$analysis,levels=rev(unique(z$analysis)))
p<-ggplot(z,aes(rr,analysis,color=cohort))+geom_vline(xintercept=1,linetype=2,color='grey55')+geom_errorbar(aes(xmin=lo,xmax=hi),orientation='y',width=.2,linewidth=.45,position=position_dodge(width=.45))+geom_point(size=2,position=position_dodge(width=.45))+scale_color_manual(values=cols)+scale_x_log10()+labs(x='Risk ratio (95% CI)',y=NULL,color=NULL)+theme(axis.text.y=element_text(size=8))
saveplot(p,'FigureS1_sensitivity',180,165)
cat('FIGURES COMPLETE\n')
