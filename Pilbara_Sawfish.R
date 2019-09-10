# Analysis of sawfish catch and effort in Pilbara trawl
#notes:
#       variable HC represents a shot,there is one row per shot

rm(list=ls(all=TRUE))

# DATA SECTION
library(lubridate)
library(chron)
library(plotrix)
library(ggplot2)
library(caret)    #names(getModelInfo())  all caret models
library(doSNOW)  #do training in paralel
library(parallel)
library(gbm)
library(pROC)    #for AUC
library(DMwR)    #for SMOTE (Synthetic Minority Over-sampling Technique  Chawla et al 2002)
library(pscl)  #zero-inflated models
library('dismo')
library(glmmADMB)
library(nlme)
library(MuMIn)  #R2 glmm Nakagawa, S, Schielzeth, H. (2013). A general and simple method for obtaining R? 
                #from Generalized Linear Mixed-effects Models. Methods in Ecology and Evolution 4: 133-142
library(brglm)  #biased corrected glm
library(dplyr)
#any(grepl("rpart", installed.packages()))  #check if rpart is installed

source("C:/Matias/Analyses/SOURCE_SCRIPTS/Git_Population.dynamics/fn.fig.R")
Do.tiff="NO"    #select figure extension
Do.jpeg="YES"
source("C:/Matias/Analyses/SOURCE_SCRIPTS/Git_other/Plot.Map.R") 

setwd('C:/Matias/Analyses/Sawfish/Pilbara')
Data=read.csv('data.csv',stringsAsFactors=F)


Current.date=as.POSIXlt(Sys.Date())

# PARAMETERS SECTION



# PROCEDURE SECTION

#1. Data manipulation
#turn times to time variable
Data$StartTime=chron(times=paste(Data$StartTime,":00",sep=''))
Data$EndTime=chron(times=paste(Data$EndTime,":00",sep=''))

#turn dates into date variable
Data$UNLOAD=as.POSIXlt(as.character(Data$UNLOAD),format='%d/%m/%Y')
Data$StartDate=as.POSIXlt(as.character(Data$StartDate),format='%d/%m/%Y')
Data$EndDate=as.POSIXlt(as.character(Data$EndDate),format='%d/%m/%Y')

#SuMERY=summary(Data)
#write.csv(SuMERY,"summary.csv")

#extract year, month and day
Data$UNLOAD.yr=with(Data,year(UNLOAD))
Data$UNLOAD.mn=with(Data,month(UNLOAD))
Data$UNLOAD.dy=with(Data,day(UNLOAD))

Data$StartDate.yr=with(Data,year(StartDate))
Data$StartDate.mn=with(Data,month(StartDate))
Data$StartDate.dy=with(Data,day(StartDate))

Data$EndDate.yr=with(Data,year(EndDate))
Data$EndDate.mn=with(Data,month(EndDate))
Data$EndDate.dy=with(Data,day(EndDate))

Data$StartDate.Time=ymd_hms(with(Data,paste(as.character(StartDate),as.character(StartTime)))) 
Data$EndDate.Time=ymd_hms(with(Data,paste(as.character(EndDate),as.character(EndTime)))) 

Data$Start.hour=hour(Data$StartDate.Time)
Data$End.hour=hour(Data$EndDate.Time)

Data$StartDate.yr=with(Data,ifelse(StartDate.yr==2107,2017,StartDate.yr))
Data$EndDate.yr=with(Data,ifelse(EndDate.yr==2107,2017,StartDate.yr))

Data$hrs.trawld=as.numeric(difftime(Data$EndDate.Time,Data$StartDate.Time,units="hours"))


#Sum dolphins and sawsharks
Data$Dolphin=Data$DolphinALIVE+Data$DolphinDEAD
Data$SawfishNarrow=Data$SawfishNarrowALIVE+Data$SawfishNarrowDEAD
Data$SawfishGreen=Data$SawfishGreenALIVE+Data$SawfishGreenDEAD

smart.par=function(n.plots,MAR,OMA,MGP) return(par(mfrow=n2mfrow(n.plots),mar=MAR,oma=OMA,las=1,mgp=MGP))

WD=getwd()

Species=c("SawfishNarrow","SawfishGreen")
Ktch.vars=c("TotalCatch","TotalDiscard","TotalBEEmp","TotalBSSnapper","TotalThreadfin",      
            "TotalRedEmp","TotalRankinCod","TotalFrypan","TotalSaddletail",     
            "TotalSpangledEmp","TotalCrimSnapper","TotalGoldband")
Effort.vars=c("SWEEPS","HEADROPE","BRIDLE","HRGridDistance")
varlist=c(Effort.vars,"AREA","SDEPTH",Ktch.vars)
Predictors=c("Skipper", "VESSEL","StartDate.yr","StartDate.mn","Start.hour",varlist,"SLAT","SLONG")
covars=c("long","lat","SLAT","SLONG",Effort.vars,"depth",Ktch.vars)


#mid points (lat,long depth)
Data$long <- apply(Data[,c("SLONG","ELONG")], 1, mean,na.rm=T)  
Data$lat <- apply(Data[,c("SLAT","ELAT")], 1, mean,na.rm=T)  
Data$depth <- apply(Data[,c("SDEPTH","EDEPTH")], 1, mean,na.rm=T)  



#flag and fix typos
Data=subset(Data,StartDate.yr>2005)

Data$depth=with(Data,ifelse(SDEPTH>125 | SDEPTH<5 ,EDEPTH,SDEPTH))
Data$depth=with(Data,ifelse(is.na(depth),EDEPTH,depth))
Data$LFB=with(Data,ifelse(LFB=="f550","F550",LFB))
Data$Skipper=with(Data,paste(FirstName,Surname))
Data$Skipper=with(Data,ifelse(Skipper==" ","Unknown",Skipper))

wrong.total.catch=c(154222, 125202, 120177, 154921, 154109, 134415)
Data=subset(Data,!HC%in%wrong.total.catch)

Wrong.hours.combined.shots=c(124582,124583,124584,124585,118776,118777,118778,118779,118780,118781,
                             118782,118783,113895,113896,113897,118785,118786,118787,118788,118789,118790,118791, 
                             150522,124586,124587,124588,124589) 
Data=subset(Data,!HC%in%Wrong.hours.combined.shots)
Data=subset(Data,hrs.trawld>=10/60) #keep records trawling for at least 10 mins (no sawfish caught in less than 10 mins)
Data=subset(Data,!HC==124586)  #combined shots into one record

do.exploratory="NO"
if(do.exploratory=="YES")
{
  setwd(paste(WD,"Preliminary",sep='/'))
  fn.date.chck=function(x,y)
  {
    id=which(Data[,match(x,names(Data))]>y)
    if(length(id>0)) return(Data[id,match(c('HC',x),names(Data))])
  }
  
  typo.UNLOAD.yr=fn.date.chck("UNLOAD.yr",year(Current.date))
  typo.StartDate.yr=fn.date.chck("StartDate.yr",year(Current.date))
  typo.EndDate.yr=fn.date.chck("EndDate.yr",year(Current.date))
  
  typo.UNLOAD.mn=fn.date.chck("UNLOAD.mn",12)
  typo.StartDate.mn=fn.date.chck("StartDate.mn",12)
  typo.EndDate.mn=fn.date.chck("EndDate.mn",12)
  
  typo.UNLOAD.dy=fn.date.chck("UNLOAD.dy",31)
  typo.StartDate.dy=fn.date.chck("StartDate.dy",31)
  typo.EndDate.dy=fn.date.chck("EndDate.dy",31)
  
  
  fun1=function(What)
  {
    x=Data[,match(What,names(Data))]
    RANG=range(x,na.rm=T)
    QUANT=quantile(x,probs=seq(0,1,.01),na.rm=T)
    plot(density(x,adjust = 2,na.rm=T),main=What,xlab="",yaxt='n',ylab="")
    legend('top',paste(RANG,collapse="-"),bty='n')
    n=length(QUANT)
    text(QUANT[n-1],0,"    99%",col=2,srt=90,font=2,pos=3,cex=1.25)
    text(QUANT[n],0,"    100%",col=2,srt=90,font=2,pos=3,cex=1.25)
  }
  fn.fig("predictor_dist",1600,2400)
  par(mfcol=c(7,3),mai=c(.2,.2,.2,.1),oma=c(1,1.25,1,.1),las=1,mgp=c(1.5,.5,0))
  for(i in 1:length(varlist))fun1(What=varlist[i])
  mtext("Density",2,outer=T,las=3,line=-1,cex=1.75)
  dev.off()
  
  fun2=function(What)
  {
    x=Data[,match(What,names(Data))]
    QUANT=quantile(x,probs=seq(0,1,.01),na.rm=T)
    plot(x,main="",xlab="",ylab="")
    mtext(What,2,las=3,line=2.35,cex=.8)
    abline(h=QUANT[match("75%",names(QUANT))],col="green")
    abline(h=QUANT[match("99%",names(QUANT))],col="orange")
    abline(h=QUANT[match("100%",names(QUANT))],col="red")
    
  }
  fn.fig("predictor_dist2",1600,2400)
  par(mfcol=c(7,3),mai=c(.2,.4,.1,.1),oma=c(.1,1.25,.1,.1),las=1,mgp=c(1,.5,0))
  for(i in 1:length(varlist))fun2(What=varlist[i])
  dev.off()
  
  
  #HOURS TRAWLED
  fn.fig("Hours_tralwed",1600,2400)
  par(mfcol=c(2,1),mai=c(.2,.2,.2,.1),oma=c(1,1.25,1,.1),las=1,mgp=c(1.5,.5,0))
  fun1(What="hrs.trawld")
  mtext("Density",2,las=3,line=1,cex=1.75)
  
  fun2(What="hrs.trawld")
  dev.off()
  
  
  #spatial distribution of species-gear interactions
  #effort as number of shots
  fn.1=function(species,VAR)
  {
    id=match(species,names(Data))
    a=Data[Data[,id]>0,]
    XLIM=range(a$SLONG)
    YLIM=range(a$SLAT)
    vaR=sort(unique(Data[,match(VAR,names(Data))]))
    #vaR=sort(unique(a[,match(VAR,names(a))]))
    CL=rgb(.1,.1,.2,alpha=0.2)
    smart.par(n.plots=length(vaR),MAR=c(2,2,1,1),OMA=c(1.75,2,.5,.1),MGP=c(1,.5,0))
    for(y in 1:length(vaR))
    {
      #with(Data[Data[,match(VAR,names(Data))]==vaR[y],],smoothScatter(SLONG,SLAT, nrpoints = 0,ylab="",xlab="",ylim=YLIM,xlim=XLIM))  #number of shots
      with(Data[Data[,match(VAR,names(Data))]==vaR[y],],plot(SLONG,SLAT,ylab="",xlab="",pch=21,col=CL,bg=CL,ylim=YLIM,xlim=XLIM))
      aa=a[a[,match(VAR,names(a))]==vaR[y],]
      if(nrow(aa)>0)with(aa,points(SLONG,SLAT,pch=21,col='white',bg='deepskyblue2',cex=1.5))   
      #if(y==1)legend('topleft',c("shots","shots with interaction"),pch=21,col=c(CL,1),pt.bg=c(CL,'deepskyblue2'),bty='n',cex=0.9)
      legend('bottomright',paste(vaR[y]," (",sum(aa[,id])," inter.)",sep=""),bty='n',cex=1.1,text.col="firebrick")
    }
    mtext("Longitude",1,outer=T)
    mtext("Latitude",2,line=0.5,las=3,outer=T)
    mtext(species,3,line=-1,outer=T)
  }
  #year
  for(s in 1:length(Species))
  {
    fn.fig(paste("spatial_inter_year_n_shots",Species[s],sep="_"),2400,2400)
    fn.1(species=Species[s],VAR='StartDate.yr')
    dev.off()
  }
  #month
  for(s in 1:length(Species))
  {
    fn.fig(paste("spatial_inter_month_n_shots",Species[s],sep="_"),2400,2400)
    fn.1(species=Species[s],VAR='StartDate.mn')
    dev.off()
  }
  #hour
  for(s in 1:length(Species))
  {
    fn.fig(paste("spatial_inter_hour_n_shots",Species[s],sep="_"),2400,2400)
    fn.1(species=Species[s],VAR='Start.hour')
    dev.off()
  }
  #Skipper
  for(s in 1:length(Species))
  {
    fn.fig(paste("spatial_inter_Skipper_n_shots",Species[s],sep="_"),2400,2400)
    fn.1(species=Species[s],VAR='Skipper')
    dev.off()
  }
  #VESSEL
  for(s in 1:length(Species))
  {
    fn.fig(paste("spatial_inter_VESSEL_n_shots",Species[s],sep="_"),2400,2400)
    fn.1(species=Species[s],VAR='VESSEL')
    dev.off()
  }
  
  
  #effort as sum of effort
  fn.effort.plot=function(DATA,species,VAR,numInt) 
  {
    id=match(species,names(Data))
    a=Data[Data[,id]>0,]
    vaR=sort(unique(Data[,match(VAR,names(Data))]))
    CL=rgb(.1,.1,.2,alpha=0.2)
    DATA$LAT=as.numeric(substr(DATA$SLAT,1,5))     #6 minute blocks
    DATA$LONG=as.numeric(substr(DATA$SLONG,1,5))  
    A=aggregate(hrs.trawld~DATA[,match(VAR,names(DATA))]+LONG+LAT,DATA,sum)
    Ymax=max(A$hrs.trawld)
    Ymin=min(A$hrs.trawld)
    Breaks=c(0,seq(Ymin,Ymax,length.out=(numInt)))
    b=range(DATA$SLAT)
    ab=range(DATA$SLONG)
    Colfunc <- colorRampPalette(c("yellow","red"))
    Couleurs=c("white",Colfunc(numInt-1))
    numberLab=10
    colLeg=(rep(c("black",rep("transparent",numberLab-1)),(numInt+1)/numberLab))
    smart.par(n.plots=length(vaR),MAR=c(2,2,1,1),OMA=c(1.75,2,.5,.1),MGP=c(1,.5,0))
    for(y in 1:length(vaR))
    {
      AA=subset(A,A[,1]==vaR[y])
      AA$LAT.cen=AA$LAT-.05
      AA$LONG.cen=AA$LONG+.05 
      AA=AA[order(AA$LAT.cen),]
      lat=unique(AA$LAT.cen)
      Reshaped=as.matrix(reshape(subset(AA,select=c(LONG.cen,LAT.cen,hrs.trawld)),idvar="LONG.cen",timevar="LAT.cen",v.names="hrs.trawld", direction="wide"))	
      Reshaped=Reshaped[order(Reshaped[,1]),]
      lon=Reshaped[,1]
      Reshaped=Reshaped[,-1]	
      image(lon,lat,z=Reshaped,xlab="",ylab="",col =Couleurs,breaks=Breaks,xlim=ab,ylim=b)
      if(y==1)color.legend(quantile(ab,probs=.9),quantile(b,probs=.5),quantile(ab,probs=.975),quantile(b,probs=.05),
                           paste(round(Breaks,0),"hrs"),rect.col=Couleurs,gradient="y",col=colLeg,cex=.7)
      box()
      aa=a[a[,match(VAR,names(a))]==vaR[y],]
      with(aa,points(SLONG,SLAT,pch=21,col='white',bg='deepskyblue2',cex=1.5))   
      legend('topleft',paste(vaR[y]," (",sum(aa[,id])," inter.)",sep=""),bty='n',cex=1.1,text.col="firebrick")
    }
    mtext("Longitude",1,outer=T)
    mtext("Latitude",2,line=0.5,las=3,outer=T)
    mtext(species,3,line=-1,outer=T)
  }
  #year
  for(s in 1:length(Species))
  {
    fn.fig(paste("spatial_inter_year_sum_effort",Species[s],sep="_"),2400,2400)
    fn.effort.plot(DATA=Data,species=Species[s],VAR='StartDate.yr',numInt=20)
    dev.off()
  }
  #month
  for(s in 1:length(Species))
  {
    fn.fig(paste("spatial_inter_month_sum_effort",Species[s],sep="_"),2400,2400)
    fn.effort.plot(DATA=Data,species=Species[s],VAR='StartDate.mn',numInt=20)
    dev.off()
  }
  #hour
  for(s in 1:length(Species))
  {
    fn.fig(paste("spatial_inter_hour_sum_effort",Species[s],sep="_"),2400,2400)
    fn.effort.plot(DATA=Data,species=Species[s],VAR='Start.hour',numInt=20)
    dev.off()
  }
  
  
  
  #barplots
  fn.2=function(species,what)
  {
    id=match(species,names(Data))
    if(what=="cpue")
    {
      a=Data
      TabYr=aggregate(a[,id]/a$hrs.trawld~StartDate.yr,a,mean)
      TabMn=aggregate(a[,id]/a$hrs.trawld~StartDate.mn,a,mean)
      TabHr=aggregate(a[,id]/a$hrs.trawld~Start.hour,a,mean)
      Tabdpth=aggregate(a[,id]/a$hrs.trawld~SDEPTH,a,mean)
    }else
    {
      a=Data[Data[,id]>0,]
      TabYr=aggregate(a[,id]~StartDate.yr,a,sum)
      TabMn=aggregate(a[,id]~StartDate.mn,a,sum)
      TabHr=aggregate(a[,id]~Start.hour,a,sum)
      Tabdpth=aggregate(a[,id]~SDEPTH,a,sum)
    }
    Efrt.yr=aggregate(hrs.trawld~StartDate.yr,Data,sum)
    Efrt.yr=subset(Efrt.yr,StartDate.yr%in%sort(TabYr$StartDate.yr))
    Efrt.Mn=aggregate(hrs.trawld~StartDate.mn,Data,sum)
    Efrt.Mn=subset(Efrt.Mn,StartDate.mn%in%sort(TabMn$StartDate.mn))
    Efrt.Hr=aggregate(hrs.trawld~Start.hour,Data,sum)
    Efrt.Hr=subset(Efrt.Hr,Start.hour%in%sort(TabHr$Start.hour))
    Efrt.dpth=aggregate(hrs.trawld~SDEPTH,Data,sum)
    Efrt.dpth=subset(Efrt.dpth,SDEPTH%in%sort(Tabdpth$SDEPTH))
    
    
    smart.par(n.plots=4,MAR=c(2.5,2,1.5,2),OMA=c(1.75,2,1,2),MGP=c(1.5,.5,0))
    
    plot(TabYr[,1],TabYr[,2],type='h',xlab='year',cex.lab=1.5,ylab='',ylim=c(0,max(TabYr[,2])))
    par(new = T)
    plot(Efrt.yr[,1],Efrt.yr[,2],type='l',lwd=2,col=2,axes=F, xlab=NA, ylab=NA,ylim=c(0,max(Efrt.yr[,2])))
    axis(side = 4)
    
    plot(TabMn[,1],TabMn[,2],type='h',xlab='month',cex.lab=1.5,ylab='',ylim=c(0,max(TabMn[,2])))
    par(new = T)
    plot(Efrt.Mn[,1],Efrt.Mn[,2],type='l',lwd=2,col=2,axes=F, xlab=NA, ylab=NA,ylim=c(0,max(Efrt.Mn[,2])))
    axis(side = 4)
    
    plot(TabHr[,1],TabHr[,2],type='h',xlab='hour',cex.lab=1.5,ylab='',ylim=c(0,max(TabHr[,2])))
    par(new = T)
    plot(Efrt.Hr[,1],Efrt.Hr[,2],type='l',lwd=2,col=2,axes=F, xlab=NA, ylab=NA,ylim=c(0,max(Efrt.Hr[,2])))
    axis(side = 4)
    
    plot(Tabdpth[,1],Tabdpth[,2],type='h',xlab='depth',cex.lab=1.5,ylab='',xlim=c(50,120),ylim=c(0,max(Tabdpth[,2])))
    par(new = T)
    plot(Efrt.dpth[,1],Efrt.dpth[,2],type='l',lwd=2,col=2,axes=F, xlab=NA, ylab=NA,xlim=c(50,120),ylim=c(0,max(Efrt.dpth[,2])))
    axis(side = 4)
    
    
    mtext(side = 4, line = 0, 'Hours trawled',outer=T,las=3,col=2,cex=1.5)
    mtext(species,3,line=-1,outer=T,cex=1.5)
    if(what=="cpue")mtext("Number per hour trawled",2,line=0.5,las=3,outer=T,cex=1.5) else mtext("Interactions",2,line=0.5,las=3,outer=T,cex=1.5)
  }
  for(s in 1:length(Species))
  {
    fn.fig(paste("year_month_hour_",Species[s],sep="_"),2400,2400)
    fn.2(species=Species[s],what='catch')
    dev.off()
  }
  for(s in 1:length(Species))
  {
    fn.fig(paste("year_month_hour_cpue_",Species[s],sep="_"),2400,2400)
    fn.2(species=Species[s],what='cpue')
    dev.off()
  }
  
  
  #boxplot
  fn.5=function(species,what)
  {
    id=match(species,names(Data))
    iid=match("hrs.trawld",names(Data))
    a=Data[Data[,id]>0,c(id,match(Predictors,names(Data)),iid)]
    smart.par(n.plots=length(Predictors),MAR=c(1,1,1,1),OMA=c(2,2.5,.5,.1),MGP=c(2,.5,0))
    for(p in 1:length(Predictors))
    {
      d=a[,c(1,match(Predictors[p],names(a)),match("hrs.trawld",names(a)))]
      if(!Predictors[p]%in%covars) d[,2]=factor(d[,2])
      if(Predictors[p]%in%covars) d[,2]=cut(d[,2],10,include.lowest=F)
      if(what=='cpue')boxplot(d[,1]/d$hrs.trawld~d[,2],xlab='',ylab="")else boxplot(d[,1]~d[,2],xlab='',ylab="")
      mtext(Predictors[p],3,-1.5,cex=.8)
    }
    if(what=='cpue')mtext(paste('Positive shots: number of',species,'caught/ hour'),2,outer=T,line=1,cex=1.25,las=3) else
      mtext(paste('Positive shots: number of',species,'caught'),2,outer=T,line=1,cex=1.25,las=3)
  }
  for(s in 1:length(Species))
  {
    fn.fig(paste("boxplot_cpue",Species[s],sep="_"),2400,2400)
    fn.5(species=Species[s],what='cpue')
    dev.off()
  }
  for(s in 1:length(Species))
  {
    fn.fig(paste("boxplot_catch",Species[s],sep="_"),2400,2400)
    fn.5(species=Species[s],what='catch')
    dev.off()
  }
  
  #plot total number by covariate
  fn.5.1=function(species)
  {
    id=match(species,names(Data))
    iid=match("hrs.trawld",names(Data))
    a=Data[,c(id,match(Predictors,names(Data)),iid)]
    smart.par(n.plots=length(Predictors),MAR=c(1,1,1,1),OMA=c(2,2.5,.5,.1),MGP=c(2,.5,0))
    for(p in 1:length(Predictors))
    {
      d=a[,c(1,match(Predictors[p],names(a)),match("hrs.trawld",names(a)))]
      if(!Predictors[p]%in%covars) d[,2]=factor(d[,2])
      if(Predictors[p]%in%covars) d[,2]=cut(d[,2],10,include.lowest=F)
      d$Pos=ifelse(d[,1]>0,1,0)
      x=aggregate(Pos~d[,2],d,sum)
      plot(x[,1],x$Pos,xlab='',ylab="",col=2)
      mtext(Predictors[p],3,-1.5,cex=.8)
    }
    mtext(paste('Number of',species,'caught'),2,outer=T,line=1,cex=1.25,las=3)
  }
  for(s in 1:length(Species))
  {
    fn.fig(paste("catch_by_predictor",Species[s],sep="_"),2400,2400)
    fn.5.1(species=Species[s])
    dev.off()
  }
  
  #catch target species
  fn.3=function(species,VAR,VAR1)
  {
    id=match(species,names(Data))
    a=Data[Data[,id]>0,]
    XLIM=range(a$SLONG)
    YLIM=range(a$SLAT)
    vaR=sort(unique(a[,match(VAR,names(a))]))
    CL=rgb(.1,.4,.1,alpha=0.2)
    smart.par(n.plots=length(vaR),MAR=c(2,2,1,1),OMA=c(1.75,2,.5,.1),MGP=c(1,.5,0))
    for(y in 1:length(vaR))
    {
      d=Data[Data[,match(VAR,names(Data))]==vaR[y],]
      d=subset(d,hrs.trawld>0)
      d$catch.rate=d[,match(VAR1,names(d))]/d$hrs.trawld
      with(d,plot(SLONG,SLAT,ylab="",xlab="",pch=21,cex=3*(d$catch.rate/max(d$catch.rate,na.rm=T)),col=CL,bg=CL,ylim=YLIM,xlim=XLIM))
      aa=a[a[,match(VAR,names(a))]==vaR[y],]
      with(aa,points(SLONG,SLAT,pch=21,col='white',bg='deepskyblue2',cex=1.5))   
      if(y==1)legend('topleft',c(paste(VAR1,"kg/hr"),"interaction"),pch=21,col=c(CL,1),pt.bg=c(CL,'deepskyblue2'),bty='n',cex=0.9)
      legend('bottomright',paste(vaR[y]," (",sum(aa[,id])," inter.)",sep=""),bty='n',cex=1.1,text.col="firebrick")
    }
    mtext("Longitude",1,outer=T)
    mtext("Latitude",2,line=0.5,las=3,outer=T)
    mtext(species,3,line=-1,outer=T)
  }
  setwd(paste(getwd(),'Catch_target',sep='/'))
  for(s in 1:length(Species)) for(v in 1:length(Ktch.vars))
  {
    fn.fig(paste("spatial",Species[s],'StartDate.yr',Ktch.vars[v],sep="_"),2400,2400)
    fn.3(species=Species[s],VAR='StartDate.yr',VAR1=Ktch.vars[v])
    dev.off()
    
    fn.fig(paste("spatial",Species[s],'StartDate.mn',Ktch.vars[v],sep="_"),2400,2400)
    fn.3(species=Species[s],VAR='StartDate.mn',VAR1=Ktch.vars[v])
    dev.off()
    
    fn.fig(paste("spatial",Species[s],'Start.hour',Ktch.vars[v],sep="_"),2400,2400)
    fn.3(species=Species[s],VAR='Start.hour',VAR1=Ktch.vars[v])
    dev.off()
  }
  
  
  fn.4=function(species,VAR,VAR2)
  {
    id=match(species,names(Data))
    id2=match(VAR,names(Data))
    a=Data[Data[,id]>0,]
    if(!VAR2=='cpue')
    {
      a$cpue=(a[,id2]/1000)    #tons 
      hist(a$cpue,main="",xlab=paste(VAR,'tons'),cex.lab=1.5,col=3)
    }
    if(VAR2=='cpue')
    {
      a$cpue=(a[,id2]/1000)/a$hrs.trawld    #tons per hour
      hist(a$cpue,main="",xlab=paste(VAR,'tons/hr'),cex.lab=1.5,col=3)
    }
    
    box()
  }
  for(s in 1:length(Species))
  {
    fn.fig(paste("interactions_by_cpue",Species[s],sep="_"),2400,2400)
    smart.par(n.plots=length(Ktch.vars),MAR=c(3.5,3.5,1,1),OMA=c(1.75,2,.5,.1),MGP=c(2,.5,0))
    for(v in 1:length(Ktch.vars)) fn.4(species=Species[s],VAR=Ktch.vars[v],VAR2='cpue') 
    dev.off()
    
    fn.fig(paste("interactions_by_catch",Species[s],sep="_"),2400,2400)
    smart.par(n.plots=length(Ktch.vars),MAR=c(3.5,3.5,1,1),OMA=c(1.75,2,.5,.1),MGP=c(2,.5,0))
    for(v in 1:length(Ktch.vars)) fn.4(species=Species[s],VAR=Ktch.vars[v],VAR2='catch') 
    dev.off()
    
  }
  
}

#Look at temperature
Temp=read.csv("C:/Matias/Data/Oceanography/SST.csv")
Temp=subset(Temp,Lat>(-21) & Lat<(-17) &Long<121 & Long>114)
plot(Temp$Long,Temp$Lat)
Temp$Temperature=Temp$value
Ag.T=aggregate(Temperature~year+month,Temp,mean)
head(Ag.T)
plot(1,col="transparent",xlim=c(1,12),ylim=c(21,31))
Yr=unique(Ag.T$Year);CL=1:12;for(i in 1:12) with(subset(Ag.T,year==Yr[i]),lines(month,Temperature,col=CL[i]))

Avg.ann.T=aggregate(Temperature~month,Temp,mean)
Data=merge(Data,Avg.ann.T,by.x="StartDate.mn",by.y="month",all.x=T)

#Add annual effort
Ann.eff=aggregate(cbind(SawfishGreen,SawfishNarrow,hrs.trawld)~StartDate.yr,Data,sum)

fn.fig("Number caught by trawled hours",2400,2400)
with(Ann.eff,plot(hrs.trawld,SawfishGreen,ylab="Number caught",pch=19))
with(Ann.eff,points(hrs.trawld,SawfishNarrow,pch=19,col=2))
legend("topleft",c("SawfishGreen",'SawfishNarrow'),pch=19,col=1:2,bty='n')
abline(v=5000)
dev.off()

names(Ann.eff)[match('hrs.trawld',names(Ann.eff))]='Ann.eff'   #in hours trawled
Data=merge(Data,subset(Ann.eff,select=c(StartDate.yr,Ann.eff)),by='StartDate.yr')


#Change vars to two-leveled factor to avoid overfitting due to small sample size
Data$Day_night=with(Data,ifelse(Start.hour%in%6:18,"Day","Night"))
Data$Season=with(Data,ifelse(StartDate.mn%in%3:8,"AuWin","SprSu")) #aggregation related to changes in Temp and algal blooms (Corey)


#Remove unvalidated records for unknown vessel
Data=subset(Data,!VESSEL=="")

#2. Data analyses
setwd(WD)

#Calculate annual catch rates
    ## Arithmetic Mean CPUE
Effort.scaler=1000  #in hours
fn.cpue=function(d,var)
{
  out = d %>%
    mutate(cpue=Effort.scaler*get(var)/hrs.trawld,
           year=as.numeric(StartDate.yr))%>%
    group_by(year) %>%
    summarise(mean = mean(cpue),
              n = length(cpue),
              sd = sd(cpue)) %>%
    mutate(lowCL = mean - 1.96*sd/sqrt(n),
           uppCL = mean + 1.96*sd/sqrt(n)) %>%
    as.data.frame
  return(out)
}
Narrow.cpue=fn.cpue(d=Data%>%select(StartDate.yr,hrs.trawld,SawfishNarrow,SawfishGreen),var="SawfishNarrow")
Green.cpue=fn.cpue(d=Data%>%select(StartDate.yr,hrs.trawld,SawfishNarrow,SawfishGreen),var="SawfishGreen")

fn.fig("Annual cpue",2400,2400)
plot(Narrow.cpue$year,Narrow.cpue$mean,ylab=paste('Numbers caught per',Effort.scaler,'trawled hour',sep=" "),xlab="Year",
     ylim=c(0,max(c(Narrow.cpue$uppCL,Green.cpue$uppCL))),pch=19,cex=1.5)
with(Narrow.cpue,segments(year,lowCL,year,uppCL,lwd=2))
points(Green.cpue$year+.2,Green.cpue$mean,cex=1.5,pch=19,col="grey60")
with(Green.cpue,segments(year+.2,lowCL,year+.2,uppCL,lwd=2,col="grey60"))
dev.off()


#Correlations among predictors
Chck.Corr="NO"
if(Chck.Corr=="YES")
{
  library(GGally)
  d=Data[,match(covars,names(Data))]   
  
  fn.fig("Covariate_correlation_catch",2400,2400)
  ggpairs(d[,-match(Effort.vars,names(d))])
  dev.off()
  
  fn.fig("Covariate_correlation_efort",2400,2400)
  ggpairs(d[,match(Effort.vars,names(d))])
  dev.off()
  
  rm(d)
}



#MODEL              
#always check there's no unique predictor levels
#dropped LAT due to strong correlation with LONG and DEPTH

TERMS=c("Season",'hrs.trawld',"VESSEL","Day_night","long","depth","TotalCatch")
# rules of thumb: to avoid over-fitting for LR models, have at least 10 "events" per 
#       covariate degree of freedom. Hence, only selected these terms (based on "catch_by_species" plot)
#       A continuous covariate would be 1 df, an N-level factor would be N-1 df. 
# Cross-validation can detect overfit models by determining how well your model generalizes 
#       to other data sets by partitioning your data. If the model does a poor job at predicting 
#       the removed observations, this indicates that the model is probably tailored to the specific
#       data points that are included in the sample and not generalizable outside the sample.


COVS=c(covars,'Ann.eff','hrs.trawld')


#Impute missing values of depth
Data$MissingDepth <- ifelse(is.na(Data$depth),"Y", "N")  #to allow tracking what records were imputed
dummy.vars <- dummyVars(~ ., data = Data[,c("HC","SLAT","SLONG","depth")])    
dummy <- predict(dummy.vars, Data[,c("HC","SLAT","SLONG","depth")])
pre.process <- preProcess(dummy, method = "bagImpute") 
imputed.data <- predict(pre.process, dummy)
Data$depth.ori=Data$depth
Data=Data[,-match("depth",names(Data))]
Data=merge(Data,imputed.data,by=c("HC","SLAT","SLONG"),all.x=T)
Data$depth=round(Data$depth)




#normalise covariates so influence on same scale
Model.d=Data[,c(Species,TERMS)]

#remove VESSEL=="Australia Bay I" as only started fishing in 2017
Model.d=subset(Model.d,!VESSEL=="Australia Bay I")

numerics=TERMS[which(TERMS%in%COVS)]

fn.scale=function(D,Numeric)
{
  procValues=preProcess(D[,Numeric],method=c('center','scale')) 
  PRED=predict(procValues,D[,Numeric])
  return(cbind(D[,-match(Numeric,names(D))],PRED)) 
}
Model.d.scaled=fn.scale(D=Model.d,Numeric=numerics)



#1. Define how to use effort in model
do.effort.cut="NO"
if(do.effort.cut=="YES")
{
  fn.effort.levels=function(species)
  {
    #select if doing analyses on normalised data or not
    d=Model.d[,c(species,TERMS)]
    
    d=subset(d,hrs.trawld<6)  #remove noise
    
    #set character variables to factors
    idd=which(!names(d)%in%c(COVS,species))
    for(f in 1:length(idd))d[,idd[f]]=as.factor(d[,idd[f]])
    d$Occurrence=d[,species]
    d=d[,-match(species,names(d))]
    d$Occurrence=ifelse(d$Occurrence>0,1,0)
    fit.gbm <- gbm.step(data=d,   #data frame
                        gbm.x = match(TERMS,names(d)),          #predictors
                        gbm.y = which(!names(d)%in%TERMS),             #response
                        family = "bernoulli",  #error structure for presence/abscence
                        tree.complexity = 5,   
                        learning.rate = 0.01,
                        bag.fraction = 0.5)
    return(fit.gbm)
  }
  
  Effort.out=vector('list',length(Species))
  names(Effort.out)=Species
  system.time({for(s in 1:length(Species))Effort.out[[s]]=fn.effort.levels(species=Species[s])})  
  
  fun.plot.gbm=function(gbm.object,var,Xlim)
  {
    gbm.call <- gbm.object$gbm.call
    k=match(var,gbm.call$predictor.names)
    response.matrix <- gbm::plot.gbm(gbm.object, k, return.grid = TRUE)
    response <- response.matrix[, 2] - mean(response.matrix[,2])
    plot(response.matrix[,1],response,type='l',lwd=2,xlim=Xlim,xlab='',
         ylab="Marginal effect on logit(p)",cex.lab=1.5)
    
  }
  
  fn.fig("Effort_bin",2400,2400)
  smart.par(n.plots=length(Species),MAR=c(2,4,1,1),OMA=rep(1,4),MGP=c(2.5,.7,0))
  for(s in 1:length(Species))
  {
    fun.plot.gbm(gbm.object=Effort.out[[s]],var="hrs.trawld",Xlim=c(.5,6))
    legend("topleft",Species[s],bty='n',cex=1.25)
  }
  mtext("Hours trawled",1,outer=T,line=-0.15,cex=1.5)
  dev.off()
  
  
  #Distribution of effort and interactions
  d=Model.d
  d$Sum=d$SawfishNarrow+d$SawfishGreen  #remove noise
  id=subset(d,Sum>0)
  d=subset(d,hrs.trawld<=max(id$hrs.trawld)) 
  
  d$hrs.trawld.bin=cut(d$hrs.trawld,breaks=seq(0,round(max(d$hrs.trawld)),.5))#30 min bins (the way the fishers operate)
  d$dummy=1
  Ag.eff=aggregate(dummy~hrs.trawld.bin,d,sum)
  
  d$Pos=d$SawfishNarrow+d$SawfishGreen
  Pos=subset(d,Pos>0 )
  
  Ag=aggregate(Pos~hrs.trawld.bin,Pos,sum)
  Add=Ag.eff$hrs.trawld.bin[which(!Ag.eff$hrs.trawld.bin%in%Ag$hrs.trawld.bin)]
  Ag=rbind(Ag,data.frame(hrs.trawld.bin=Add,Pos=0))
  Ag=Ag[order(Ag$hrs.trawld.bin),]
  

    
  fn.fig("Effort.dist_interactions_dist",2400,2400)
  smart.par(n.plots=2,MAR=c(2,4,1,1),OMA=rep(1,4),MGP=c(3,.7,0))
  barplot(Ag.eff$dummy,names.arg=Ag.eff$hrs.trawld.bin,ylab="Number of records",cex.lab=1.5)
  box()
  barplot(Ag$Pos,names.arg=Ag$hrs.trawld.bin,ylab="Number of interactions",cex.lab=1.5)
  box()
  mtext("Hours trawled",1,outer=T,line=-0.15,cex=1.5)
  
  #inset
  par(fig=c(0.575,0.975,.15,0.45), new = T,mgp=c(1,.5,0),las=1)
#  with(Pos,plot(hrs.trawld,Pos,pch=19,xlab="",ylab="",cex.axis=.7,cex=.8,ylim=c(0,max(Pos))))
  # mtext("Hours trawled",1,line=1.5,cex=1)
  # mtext("Number of interactions",2,las=3,line=1.5,cex=1)
  
  
  Test=merge(Ag,Ag.eff,by="hrs.trawld.bin")
  Test$bin.mid=seq(.25,5.75,.5)
  Test$Prop.pos.rec=with(Test,Pos/dummy)
  plot(Prop.pos.rec~hrs.trawld.bin,Test,ylab="",xlab="",cex.axis=.7)
  mtext("Hours trawled",1,line=1.5,cex=.8)
  mtext("Proportion of records with interaction",2,las=3,line=1.5,cex=.8)
  
  dev.off()
  
  Test.same.dist=merge(Ag,Ag.eff,"hrs.trawld.bin")
  Chi.test=chisq.test(Test.same.dist$Pos, Test.same.dist$dummy)
  Pvalue=Chi.test$p.value
  
  # Bin.mid=seq(.25,5.75,.5)
  # y=Bin.mid*Ag$Pos
  # x=Bin.mid*Ag.eff$dummy
  # plot(x,y)
}


#2. Apply binomial and boosted regresssion trees
#TERMS=c(TERMS[-match("hrs.trawld",TERMS)],"Eff.bin")

Model.d$Eff.bin=factor(with(Model.d,ifelse(hrs.trawld<2,"Short",
                ifelse(hrs.trawld>=2 & hrs.trawld<4,"Standard","Long"))))


#Check if factors more than 32 levels (Random Forest won't work)
too.many.levels <- function(x)  is.factor(x) == TRUE & length(levels(x)) > 32
delete <- lapply(Model.d, too.many.levels)


#run models
#notes: 
#1. apply SMOTE for oversampling rare events (i.e. imbalanced data)  see http://amunategui.github.io/smote/
#smote() uses bootstrapping and k-n nearest neighbour to synthetically create more 
# observations of the rare event
# Synthetic minority sampling technique (SMOTE): down samples the majority class and synthesizes 
# new minority instances by interpolating between existing one

#2. AUC is best way to understand outcomes of rare event modelling (more powerful than accuracy 
# because if the model predicts correctly all the 0s, i.e. high accuracy, I don't care, I need
# to predict the 1s)

#3. classification models (response is factor: yes/no)

#parallel processing
CORES=detectCores()-1
cat(CORES, " cores detected.")


#Fit models
Model.out=vector('list',length(Species))
names(Model.out)=Species
Best.MDL=FinalModel=Model.out

fn.models=function(species,METRIC)
{
  #select if doing analyses on normalised data or not
  d=Model.d[,c(species,TERMS)]
  
  #remove noise
  id=d[d[,match(species,names(d))]>0,]
  d=subset(d,hrs.trawld<=max(id$hrs.trawld)) 
  
  #set character variables to factors   
  idd=which(!names(d)%in%c(COVS,species))
  for(f in 1:length(idd)) d[,idd[f]]=as.factor(d[,idd[f]])
  d$Occurrence=d[,species]
  d=d[,-match(species,names(d))]
  d$Occurrence=ifelse(d$Occurrence>0,1,0)
  d$Occurrence=factor(ifelse(d$Occurrence==1,"yes","no"),levels=c("yes","no"))    #have presence as first levels
  
  
  #weight   
  # d$hrs.trawld.bin=cut(d$hrs.trawld,breaks=seq(0,round(max(d$hrs.trawld)),.5))#30 min bins (the way the fishers operate)
  # Tab.weight=as.matrix(table(d$hrs.trawld.bin))
  # Tab.weight=data.frame(hrs.trawld.bin=row.names(Tab.weight),weights=1/Tab.weight)
  # d=merge(d,Tab.weight,by="hrs.trawld.bin")
  
  
  #Partition data
  ind <- createDataPartition(d$Occurrence,times = 1,p = 0.7,list = FALSE)
  train <- d[ind,]
  test <- d[-ind,]
  
  #Smote data
  train.smote <- SMOTE(Occurrence~ ., train, perc.over = 100, perc.under=200)
  #train.smote$Occurrence=as.character(train.smote$Occurrence)
  #train.smote$Occurrence=with(train.smote,ifelse(Occurrence=="1",1,0))
  #prop.table(table(train.smote$Occurrence))
  
  
  #Fit models
  
  #Formula=formula(paste("Occurrence",paste(TERMS,collapse="+"),sep="~"))
  
  #Binomial
  BIN=glm(Occurrence~.,data =train, family="binomial", maxit=500)
  
  #Binomial to smote data
  BIN.smote=glm(Occurrence~.,data =train.smote, family="binomial", maxit=500)
  
  #Binomial brglm
  Brglm <- brglm(Occurrence~., data = train)
  
  
  #GBM
  #register cluster to train in parallel (only works for Xgboost)
  cl <-makeCluster(CORES, type = "SOCK") 
  registerDoSNOW(cl)
  
  #Control
  ctrl <- trainControl(method = "repeatedcv",
                       number = 10,     #split data in ten ways to build 10 models
                       repeats = 3,     #repeat 10 fold cross validation 3 times
                       classProbs = TRUE,     #needed to score models using AUC
                       summaryFunction = twoClassSummary,
                       search = "grid")     #find optimal collection of parameters using default grid
  
  #gbm
  fit.gbm <- train(Occurrence~.,data=train, method='gbm', metric = METRIC,trControl=ctrl)
  #gbm smote
  fit.gbm.smote <- train(Occurrence~.,data=train.smote, method='gbm', metric = METRIC,trControl=ctrl)
  

  #clear parameters for parallel computation
  stopCluster(cl)
  registerDoSEQ()
  remove(cl)
  
  
  return(list(DATA=d,train=train,test=test,BIN=BIN,BIN.smote=BIN.smote,
              Brglm=Brglm,GBM=fit.gbm,GBM.smote=fit.gbm.smote))
}
system.time({for(s in 1:length(Species))Model.out[[s]]=fn.models(species=Species[s],METRIC="Sens")})

#fit selected model
final.mod=function(d,MOD)
{
  if(MOD=="BIN")
  {
    Modl=glm(Occurrence~.,data =d, family="binomial", maxit=500)
    
    T.table=summary(Modl)$coefficients
    Percent.dev.exp=100*(Modl$null.deviance-Modl$deviance)/Modl$null.deviance
    
    Anova.tab=anova(Modl,test="Chisq")
    
    n=2:length(Anova.tab$Deviance)
    Term.dev.exp=100*(Anova.tab$Deviance[n]/Modl$null.deviance)
    names(Term.dev.exp)=rownames(Anova.tab)[n]
    Anov.tab=as.data.frame.matrix(Anova.tab)
    Term.tab=data.frame(Percent.dev.exp=Term.dev.exp)
    Anova.tab=Anova.tab[-1,match(c("Deviance","Pr(>Chi)"),names(Anova.tab))]
    Anova.tab=cbind(Anova.tab,Term.tab)
    Anova.tab=Anova.tab[,-match("Deviance",names(Anova.tab))]
    Anova.tab$"Pr(>Chi)"=ifelse(Anova.tab$"Pr(>Chi)"<0.001,"<0.001",round(Anova.tab$"Pr(>Chi)",3))
    Total=Anova.tab[1,]
    Total$"Pr(>Chi)"=""
    Total$Percent.dev.exp=sum(Anova.tab$Percent.dev.exp)
    rownames(Total)="Total"
    Anova.tab=rbind(Anova.tab,Total)
    Anova.tab$Percent.dev.exp=round(Anova.tab$Percent.dev.exp,2)
    return(list(Modl=Modl,Anova.tab=Anova.tab,COEFS=T.table))
  }
  
}
for(s in 1:length(Species))FinalModel[[s]]=final.mod(d=Model.out[[s]]$DATA,MOD="BIN")





# REPORT SECTION
setwd(paste(WD,"paper",sep="/"))

#Figure 1. map
add.depth="NO"
if(add.depth=="YES")
{
  Bathymetry_120=read.table("C:/Matias/Data/Mapping/get_data112_120.cgi")
  Bathymetry_138=read.table("C:/Matias/Data/Mapping/get_data120.05_138.cgi")
  Bathymetry=rbind(Bathymetry_120,Bathymetry_138)
  Bathymetry=Bathymetry[order(Bathymetry$V1,Bathymetry$V2),]
  xbat=sort(unique(Bathymetry$V1))
  ybat=sort(unique(Bathymetry$V2))
  reshaped=as.matrix(reshape(Bathymetry,idvar="V1",timevar="V2",v.names="V3", direction="wide"))
}
CL1='black'
fn.Fig1=function(DATA,species,numInt) 
{
  id=match(species,names(Data))
  a=Data[Data[,id]>0,]
  CL=rgb(.1,.1,.2,alpha=0.2)
  DATA$LAT=as.numeric(substr(DATA$SLAT,1,5))     #6 minute blocks
  DATA$LONG=as.numeric(substr(DATA$SLONG,1,5))  
  A=aggregate(hrs.trawld~LONG+LAT,DATA,sum)
  Ymax=max(A$hrs.trawld)
  Ymin=min(A$hrs.trawld)
  Breaks=c(0,seq(Ymin,Ymax,length.out=(numInt)))
  b=range(DATA$SLAT)
  ab=range(DATA$SLONG)
  #Colfunc <- colorRampPalette(c("yellow","red"))
  Colfunc <- colorRampPalette(c("grey90","grey10"))
  Couleurs=c("white",Colfunc(numInt-1))
  numberLab=10
  colLeg=(rep(c("black",rep("transparent",numberLab-1)),(numInt+1)/numberLab))
  AA=A
  AA$LAT.cen=AA$LAT-.05
  AA$LONG.cen=AA$LONG+.05 
  AA=AA[order(AA$LAT.cen),]
  lat=unique(AA$LAT.cen)
  Reshaped=as.matrix(reshape(subset(AA,select=c(LONG.cen,LAT.cen,hrs.trawld)),idvar="LONG.cen",timevar="LAT.cen",v.names="hrs.trawld", direction="wide"))	
  Reshaped=Reshaped[order(Reshaped[,1]),]
  lon=Reshaped[,1]
  Reshaped=Reshaped[,-1]	
  LNG=aa
  LAT=seq(round(b[1]),round(b[2]),0.5)
  aa=round(ab[1]):round(ab[2])
  bb=seq(b[1],b[2],length.out=length(aa))
  PLATE=c(.01,.9,.075,.9)
  plotmap(aa,bb,PLATE,"dark grey",ab,b)
  image(lon,lat,z=Reshaped,xlab="",ylab="",yaxt='n',col =Couleurs,breaks=Breaks,add=T)
  box()
  axis(side = 1, at =LNG, labels = LNG, tcl = .5)
  
  if(s==1)color.legend(quantile(ab,probs=.9),quantile(b,probs=.5),quantile(ab,probs=.975),quantile(b,probs=.05),
                       paste(round(Breaks,0),"hrs"),rect.col=Couleurs,gradient="y",col=colLeg,cex=1)
  box()
  axis(2,seq(min(lat),max(lat),.25),-seq(min(lat),max(lat),.25))
  with(a,points(SLONG,SLAT,pch=21,col='white',bg=CL1,cex=1.2))
  legend("topleft",species,cex=1.5,bty="n")
  if(add.depth=="YES")
  {
    contour(xbat, ybat, reshaped[,2:ncol(reshaped)], zlim=c(-1,-300),
            nlevels = 5,labcex=1,col="grey30",add=T)
  }
}
fn.fig("Figure 1",1600,2400)
smart.par(n.plots=length(Species),MAR=c(2,2,.1,1),OMA=c(1,2,.2,.2),MGP=c(2.5,.5,0))
for(s in 1:length(Species))  fn.Fig1(DATA=Data,species=Species[s],numInt=20)
mtext("Longitude (?E)",1,outer=T,cex=1.5,line=-0.25)
mtext("Latitude (?S)",2,las=3,outer=T,cex=1.5,line=0.6)

#inset OZ
library(PBSmapping)
data(worldLLhigh)
par(fig=c(0.575,0.975,.05,0.225), new = T,mgp=c(1,.5,0),las=1)
#par(fig=c(0.575,0.975,.15,0.45), new = T,mgp=c(1,.5,0),las=1)
plotMap(worldLLhigh, xlim=c(113,155),ylim=c(-44.5,-11),plt = c(.1, 1, 0.075, 1),
        col="grey30",tck = 0.025, tckMinor = 0.0125, xlab="",ylab="",axes=F)
box(lwd=2)
polygon(x=c(116,116,120,120),y=c(-18,-21,-21,-18),lwd=1.5,col="grey65")
dev.off()


#Figure 2 Barplot of interactions
fn.2=function(species,what)
{
  id=match(species,names(Data))
  a=Data[Data[,id]>0,]
  a$D=10*round(a$depth/10)
  TabMn=aggregate(a[,id]~StartDate.mn,a,sum)
  TabHr=aggregate(a[,id]~Start.hour,a,sum)
  Tabdpth=aggregate(a[,id]~D,a,sum)
  b=Data
  b$D=10*round(b$depth/10)
  b$dummy=1
  Efrt.Mn=aggregate(hrs.trawld~StartDate.mn,b,sum)
  Efrt.Mn=subset(Efrt.Mn,StartDate.mn%in%sort(TabMn$StartDate.mn))
  Shots.Mn=aggregate(dummy~StartDate.mn,b,sum)
  Shots.Mn=subset(Shots.Mn,StartDate.mn%in%sort(TabMn$StartDate.mn))
  Efrt.Hr=aggregate(hrs.trawld~Start.hour,b,sum)
  Efrt.Hr=subset(Efrt.Hr,Start.hour%in%sort(TabHr$Start.hour))
  Shots.Hr=aggregate(dummy~Start.hour,b,sum)
  Shots.Hr=subset(Shots.Hr,Start.hour%in%sort(TabHr$Start.hour))
  Efrt.dpth=aggregate(hrs.trawld~D,b,sum)
  Efrt.dpth=subset(Efrt.dpth,D%in%sort(Tabdpth$D))
  Shots.dpth=aggregate(dummy~D,b,sum)
  Shots.dpth=subset(Shots.dpth,D%in%sort(Tabdpth$D))
  LISTA=list(list(TabMn,Efrt.Mn,Shots.Mn),list(TabHr,Efrt.Hr,Shots.Hr),list(Tabdpth,Efrt.dpth,Shots.dpth))
  names(LISTA)=c("Month","Hour","Depth (m)")
  third.axis="grey35"
  for(l in 1:length(LISTA))
  {
    barplot(LISTA[[l]][[1]][,2],col="grey95",cex.names=.85,names.arg=LISTA[[l]][[1]][,1],ylim=c(0,max(LISTA[[l]][[1]][,2])))
    par(new = T)
    plot(LISTA[[l]][[2]][,1],LISTA[[l]][[2]][,2],type='l',lwd=2,col=1,axes=F, xlab=NA, ylab=NA,ylim=c(0,max(LISTA[[l]][[2]][,2])))
    if(s==2)axis(side = 4)
    par(new = T)
    plot(LISTA[[l]][[3]][,1],LISTA[[l]][[3]][,2],type='l',lwd=2,col=third.axis,axes=F, xlab=NA, ylab=NA,ylim=c(0,max(LISTA[[l]][[3]][,2])))
    if(s==2)axis(side = 4.25, line=4.5,col=third.axis,col.axis=third.axis)
    box()
    if(l==1)mtext(species,3,line=0.25,cex=1.5)
    mtext(names(LISTA)[l],1,cex=1.25,line=2)
  }
  mtext(side = 4, line = 2, 'Hours trawled',outer=T,las=3,cex=1.35)
  mtext(side = 4, 'Fishing events',outer=T,las=3,col=third.axis,cex=1.35,line=6.5)
  mtext("Interactions",2,line=0,las=3,outer=T,cex=1.35)
}
fn.fig("Figure 2",2100,2400)
par(mfcol=c(3,2),mar=c(2.5,2,2,1),oma=c(1,1.5,1,8),mgp=c(2.5,.5,0),las=1)
for(s in 1:length(Species))  fn.2(species=Species[s],what='catch')
dev.off()

#Barplot of interactions by landings
fn.3=function(species,what)
{
  id=match(species,names(Data))
  a=Data[Data[,id]>0,]
  a$D=10*round(a$depth/10)
  TabMn=aggregate(a[,id]~StartDate.mn,a,sum)
  TabHr=aggregate(a[,id]~Start.hour,a,sum)
  Tabdpth=aggregate(a[,id]~D,a,sum)
  b=Data
  b$D=10*round(b$depth/10)
  b$dummy=1
  Efrt.Mn=aggregate((TotalCatch/1000)~StartDate.mn,b,sum)
  Efrt.Mn=subset(Efrt.Mn,StartDate.mn%in%sort(TabMn$StartDate.mn))
  
  Efrt.Hr=aggregate((TotalCatch/1000)~Start.hour,b,sum)
  Efrt.Hr=subset(Efrt.Hr,Start.hour%in%sort(TabHr$Start.hour))
  Efrt.dpth=aggregate((TotalCatch/1000)~D,b,sum)
  Efrt.dpth=subset(Efrt.dpth,D%in%sort(Tabdpth$D))
  LISTA=list(list(TabMn,Efrt.Mn,Shots.Mn),list(TabHr,Efrt.Hr,Shots.Hr),list(Tabdpth,Efrt.dpth,Shots.dpth))
  names(LISTA)=c("Month","Hour","Depth (m)")
  third.axis="grey35"
  for(l in 1:length(LISTA))
  {
    barplot(LISTA[[l]][[1]][,2],col="grey95",cex.names=.85,names.arg=LISTA[[l]][[1]][,1],ylim=c(0,max(LISTA[[l]][[1]][,2])))
    par(new = T)
    plot(LISTA[[l]][[2]][,1],LISTA[[l]][[2]][,2],type='l',lwd=2,col=1,axes=F, xlab=NA, ylab=NA,ylim=c(0,max(LISTA[[l]][[2]][,2])))
    if(s==2)axis(side = 4)
    box()
    if(l==1)mtext(species,3,line=0.25,cex=1.5)
    mtext(names(LISTA)[l],1,cex=1.25,line=2)
  }
  mtext(side = 4, line = 2, 'Total landings (tons)',outer=T,las=3,cex=1.35)
  mtext("Interactions",2,line=0,las=3,outer=T,cex=1.35)
}
fn.fig("Interactions by landings",2100,2400)
par(mfcol=c(3,2),mar=c(2.5,2,2,1),oma=c(1,1.5,1,4),mgp=c(2.5,.5,0),las=1)
for(s in 1:length(Species))  fn.3(species=Species[s],what='catch')
dev.off()

#Table 1
SawfishGreenALIVE=sum(Data$SawfishGreenALIVE)
SawfishGreenDEAD=sum(Data$SawfishGreenDEAD)
SawfishNarrowALIVE=sum(Data$SawfishNarrowALIVE)
SawfishNarrowDEAD=sum(Data$SawfishNarrowDEAD)

Tot.records=nrow(Data)
Tot.effort=sum(Data$hrs.trawld)
Tab1=data.frame(GreenALIVE=SawfishGreenALIVE,GreenDEAD=SawfishGreenDEAD,
                NarrowALIVE=SawfishNarrowALIVE,hNarrowDEAD=SawfishNarrowDEAD,
                Tot.records=Tot.records,Tot.effort=Tot.effort)
write.csv(Tab1,"Tab1.csv")

#Figure S1. histogram of number caught
library("plotrix")
fn.S1=function(species)
{
  id=match(species,names(Data))
  a=Data[Data[,id]>0,]
  Percent.occur=round(100*nrow(a)/nrow(Data),2)
  TAB=table(Data[,id])
  barplot(TAB, log="y",ylab="",xlab="")
  box()
  legend('topright',paste(Percent.occur,"% of records",sep=""),bty='n',title=species,cex=1.25)
}
fn.fig("S1_Semi_log_histogram_num_caught",2000,2400)
smart.par(n.plots=length(Species),MAR=c(1.75,2.5,.1,1),OMA=c(1.75,2,.5,.1),MGP=c(2,.5,0))
for(s in 1:length(Species))fn.S1(species=Species[s]) 
mtext("Number of indviduals caught",1,outer=T,cex=1.5)
mtext("Number of records",2,line=0.5,las=3,outer=T,cex=1.5)
dev.off()

#Figure S2. Interactions by year
CL1='deepskyblue3'
CL2='forestgreen'
fn.S2=function(DATA,VAR,numInt) 
{
  id=match(Species[1],names(Data))
  a=Data[Data[,id]>0,]
  iid=match(Species[2],names(Data))
  a1=Data[Data[,iid]>0,]
  
  vaR=sort(unique(a[,match(VAR,names(a))]))
  #vaR=sort(unique(Data[,match(VAR,names(Data))]))
  CL=rgb(.1,.1,.2,alpha=0.2)
  DATA$LAT=as.numeric(substr(DATA$SLAT,1,5))     #6 minute blocks
  DATA$LONG=as.numeric(substr(DATA$SLONG,1,5))  
  A=aggregate(hrs.trawld~DATA[,match(VAR,names(DATA))]+LONG+LAT,DATA,sum)
  Ymax=max(A$hrs.trawld)
  Ymin=min(A$hrs.trawld)
  Breaks=c(0,seq(Ymin,Ymax,length.out=(numInt)))
  b=range(DATA$SLAT)
  ab=range(DATA$SLONG)
  Colfunc <- colorRampPalette(c("yellow","red"))
  Couleurs=c("white",Colfunc(numInt-1))
  numberLab=10
  colLeg=(rep(c("black",rep("transparent",numberLab-1)),(numInt+1)/numberLab))
  smart.par(n.plots=length(vaR),MAR=c(2,2,1,1),OMA=c(1.75,2,.5,.1),MGP=c(1,.5,0))
  for(y in 1:length(vaR))
  {
    AA=subset(A,A[,1]==vaR[y])
    AA$LAT.cen=AA$LAT-.05
    AA$LONG.cen=AA$LONG+.05 
    AA=AA[order(AA$LAT.cen),]
    lat=unique(AA$LAT.cen)
    Reshaped=as.matrix(reshape(subset(AA,select=c(LONG.cen,LAT.cen,hrs.trawld)),idvar="LONG.cen",timevar="LAT.cen",v.names="hrs.trawld", direction="wide"))	
    Reshaped=Reshaped[order(Reshaped[,1]),]
    lon=Reshaped[,1]
    Reshaped=Reshaped[,-1]	
    image(lon,lat,z=Reshaped,xlab="",ylab="",col =Couleurs,breaks=Breaks,xlim=ab,ylim=b)
    if(y==1)color.legend(quantile(ab,probs=.9),quantile(b,probs=.5),quantile(ab,probs=.975),quantile(b,probs=.05),
                         paste(round(Breaks,0),"hrs"),rect.col=Couleurs,gradient="y",col=colLeg,cex=.7)
    box()
    aa=a[a[,match(VAR,names(a))]==vaR[y],]
    with(aa,points(SLONG,SLAT,pch=21,col='white',bg=CL1,cex=1.5))   
    aa1=a1[a1[,match(VAR,names(a1))]==vaR[y],]
    with(aa1,points(SLONG,SLAT,pch=21,col='white',bg=CL2,cex=1.5))   
    
    legend('bottom',paste(vaR[y]),bty='n',cex=1.1)
    legend('topleft',c(paste(sum(aa[,id]),Species[1]),paste(sum(aa1[,iid]),Species[2])),
           bty='n',cex=1.1,text.col=c(CL1,CL2))
  }
  mtext("Longitude (?E)",1,outer=T)
  mtext("Latitude (?S)",2,line=0.5,las=3,outer=T)
}
fn.fig("S2_interactions_by_year",2400,2400)
fn.S2(DATA=Data,VAR='StartDate.yr',numInt=20)
dev.off()


#AUC figure. Compare models
source("C:/Matias/Analyses/SOURCE_SCRIPTS/Smart_par.R")
fun.pred.auc=function(moDl,NMS,test)
{
  #pred probability
  if(class(moDl)[1]=="train"){pred.prob=predict(moDl,test,type='prob',n.trees=moDl$bestTune$n.trees)}else
  {
    pred.prob=predict(moDl,test,'response')
    pred.prob=data.frame(no=pred.prob)
    pred.prob$yes=1-pred.prob$no
  }
  #auc
  auc=roc(test$Occurrence, pred.prob[,match('yes',names(pred.prob))])  #.5 is random model; 1 is perfect model
  Preds=cbind(pred.prob)
  names(Preds)=paste(NMS,names(Preds),sep="_")
  return(list(Preds=Preds,auc=auc))
}
best.mdl=function(modl)
{
  Fits=names(modl)[-match(c("test","train","DATA"),names(modl))]
  AUCs=vector('list',length(Fits))
  names(AUCs)=Fits
  for(f in 1:length(Fits))
  {
    AUCs[[f]]=fun.pred.auc(moDl=modl[[match(Fits[f],names(modl))]],NMS=Fits[f],test=modl[[match('test',names(modl))]])   
  }
  
  return(list(AUC=AUCs))
}
for(s in 1:length(Species))Best.MDL[[s]]=best.mdl(modl=Model.out[[s]])

fn.fig("AUC",1600,2400)
smart.par(n.plots=length(Species),MAR=c(4,1,1,1),OMA=rep(1,4),MGP=c(2.5,.7,0))
for(s in 1:length(Species))
{
  plot(Best.MDL[[s]]$AUC[[1]]$auc)
  aucs=rep(NA,length(Best.MDL[[s]]$AUC))
  aucs[1]=Best.MDL[[s]]$AUC[[1]]$auc$auc[1]
  CL=2:length(Best.MDL[[s]]$AUC)
  for(l in CL)
  {
    lines(Best.MDL[[s]]$AUC[[l]]$auc,col=CL[l-1])
    aucs[l]=Best.MDL[[s]]$AUC[[l]]$auc$auc[1]
  }
  nombres=names(Best.MDL[[s]]$AUC)  
  #nombres=read.table(text = nombres, sep = ".", as.is = TRUE)$V2
  legend('bottomright',paste(nombres,round(aucs,3)),lty=1,col=c(1,CL),bty='n',title="AUC",cex=1)
  legend("topleft",Species[s],bty='n')
}
dev.off()


#Marginal effect
fn.marg.eff=function(d,modl)
{
  #Season effect
  Ef.seq=seq(min(round(d$hrs.trawld,1)),round(Eff.quantile[match("90%",names(Eff.quantile))],1),by=.1)
  Eff.quantile=quantile(d$hrs.trawld,probs=seq(0,1,.05))
  nd <- data.frame(Season = factor(rep(c("AuWin","SprSu"),each=length(Ef.seq)),levels=levels(d$Season)), 
                   VESSEL = factor(names(rev(sort(table(d$VESSEL)))[1]),levels=levels(d$VESSEL)),
                   Day_night = factor(names(rev(sort(table(d$Day_night)))[1]),levels=levels(d$Day_night)),
                   long = mean(d$long),
                   TotalCatch = mean(d$TotalCatch),
                   depth = mean(d$depth),
                   hrs.trawld = rep(Ef.seq,2))
  PREDS=predict(modl,nd,type='response',se.fit =T)
  nd$pred=1-PREDS$fit
  nd$SE=PREDS$se.fit
  
  Ymax=2*max(nd$pred+1.96*nd$SE)
  with(subset(nd,Season=='AuWin'),plot(hrs.trawld,pred,ylab="",xlab="",axes=FALSE,pch=19,ylim=c(0,Ymax)))
  with(subset(nd,Season=='AuWin'),segments(hrs.trawld,pred+1.96*SE,hrs.trawld,pred-1.96*SE))
  with(subset(nd,Season=='SprSu'),points(hrs.trawld,pred,pch=19,col="grey60"))
  with(subset(nd,Season=='SprSu'),segments(hrs.trawld,col="grey60",pred+1.96*SE,hrs.trawld,pred-1.96*SE))
  
  box()
  Q=quantile(nd$pred)
  axis(1,nd$hrs.trawld,F,tck=-0.015)
  axis(1,seq(min(nd$hrs.trawld),max(nd$hrs.trawld),.2),F,tck=-0.03)
  if(s==2) axis(1,seq(min(nd$hrs.trawld),max(nd$hrs.trawld),.2),seq(min(nd$hrs.trawld),max(nd$hrs.trawld),.2),tck=-0.03)  
  axis(2,at=c(Q,Ymax*.5,Ymax),labels=c(round(Q,3),0.5,1))
  axis.break(2,1.1*Q[length(Q)])
  legend("topright",Species[s],bty='n',cex=1.25)
  if(s==1) legend("topleft",c("Au-Wi","Sp-Su"),pch=19,col=c("black","grey60"),bty='n')
  
  
  #Average seasonal effect
  nd <- data.frame(Season = factor(c("AuWin","SprSu")), 
                   VESSEL = factor(names(rev(sort(table(d$VESSEL)))[1]),levels=levels(d$VESSEL)),
                   Day_night = factor(names(rev(sort(table(d$Day_night)))[1]),levels=levels(d$Day_night)),
                   long = mean(d$long),
                   TotalCatch = mean(d$TotalCatch),
                   depth = mean(d$depth),
                   hrs.trawld = mean(d$hrs.trawld))
  PREDS.average=predict(modl,nd,type='response',se.fit =T)
  Avg.pred=paste(round(1-PREDS.average$fit,4),round(PREDS.average$se.fit,4),sep="?")
  names(Avg.pred)=c("AuWin","SprSu")
  return(Avg.pred)
}
Avg.prob.season=Model.out
fn.fig("Seasonal_effect",1600,2400)
smart.par(n.plots=length(Species),MAR=c(2,2,.1,1),OMA=c(1,2,.2,.2),MGP=c(2.5,.5,0))
par(las=1)
for(s in 1:length(Species)) Avg.prob.season[[s]]=fn.marg.eff(d=Model.out[[s]]$DATA,modl=FinalModel[[s]]$Modl)
mtext("Hours trawled",1,outer=T,cex=1.5,line=-.25)
mtext("Probability",2,outer=T,line=.5,las=3,cex=1.5)
dev.off()

write.csv(Avg.prob.season,"Avg.prob.season.csv")

#Export Anova table
NMS=colnames(FinalModel$SawfishNarrow$Anova.tab)
TABL.Anova=cbind(FinalModel$SawfishNarrow$Anova.tab,FinalModel$SawfishGreen$Anova.tab)
colnames(TABL.Anova)=c(paste("SawfishNarrow",NMS,sep="_"),
                       paste("SawfishGreen",NMS,sep="_"))

NMS=colnames(FinalModel$SawfishNarrow$COEFS)
TABL.COEF=cbind(FinalModel$SawfishNarrow$COEFS,FinalModel$SawfishGreen$COEFS)
colnames(TABL.COEF)=c(paste("SawfishNarrow",NMS,sep="_"),
                       paste("SawfishGreen",NMS,sep="_"))

write.csv(TABL.Anova,"Anova.csv",row.names=T)
write.csv(TABL.COEF,"COEF.csv",row.names=T)



#######

#NOTE used
NOT.used="NO"
if(NOT.used=="YES")
{
  fn.effort.levels=function(species)
  {
    #select if doing analyses on normalised data or not
    d=Model.d[,c(species,TERMS)]
    
    #remove noise
    id=d[d[,match(species,names(d))]>0,]
    d=subset(d,hrs.trawld<=max(id$hrs.trawld)) 
    
    #set character variables to factors   
    idd=which(!names(d)%in%c(COVS,species))
    for(f in 1:length(idd)) d[,idd[f]]=as.factor(d[,idd[f]])
    d$Occurrence=d[,species]
    d=d[,-match(species,names(d))]
    d$Occurrence=ifelse(d$Occurrence>0,1,0)
    
    #weight   
    d$hrs.trawld.bin=cut(d$hrs.trawld,breaks=seq(0,round(max(d$hrs.trawld)),.5))#30 min bins (the way the fishers operate)
    Tab.weight=as.matrix(table(d$hrs.trawld.bin))
    Tab.weight=data.frame(hrs.trawld.bin=row.names(Tab.weight),weights=1/Tab.weight)
    d=merge(d,Tab.weight,by="hrs.trawld.bin")
    
    #Formula=formula(paste("Occurrence",paste(TERMS,collapse="+"),sep="~"))
    
    
    #Mixed effect binomial
    BIN=lme(Occurrence ~ Season + Day_night + long + depth + TotalCatch + 
              hrs.trawld, data = d, random = ~1| VESSEL)
    
    #BIN=model<- glm(Formula,data =d, family="binomial", maxit=500)
    #ZIB<- glmmadmb(Formula,data=d,family="binomial",zeroInflation=TRUE)
    
    
    return(list(DATA=d,BIN=BIN))
  }
  
  Model.out=vector('list',length(Species))
  names(Model.out)=Species
  Model.show=Model.out
  system.time({for(s in 1:length(Species))Model.out[[s]]=fn.effort.levels(species=Species[s])})
  
  
  
  #compute standard error for predictions
  fn.pred.glmm.SE=function(mod,dat)
  {
    dat$pred=predict(mod, newdata = dat, type = "response")
    Designmat <- model.matrix(eval(eval(mod$call$fixed)[-2]), dat[,-ncol(dat)])
    predvar <- diag(Designmat %*% mod$varFix %*% t(Designmat))
    dat$SE <- sqrt(predvar) 
    # library(AICcmodavg)
    # a=predictSE(mod, nd, se.fit = TRUE)
    # nd$pred=a$fit
    # nd$SE=a$se.fit
    return(dat)
  }
  
  #Show model fits, predictions and term effects
  fn.show=function(MOD)
  {
    
    d=MOD$DATA
    BIN=MOD$BIN
    
    Anova.Table=anova(BIN)
    T.table=summary(BIN)$tTable
    
    #R2=100*r.squaredGLMM(BIN)   #percentage
    

    #Create data for predicting effects
    
    #Season effect
    nd <- data.frame(Season = factor(rep(c("AuWin","SprSu"),each=20),levels=levels(d$Season)), 
                     VESSEL = factor(names(rev(sort(table(d$VESSEL)))[1]),levels=levels(d$VESSEL)),
                     Day_night = factor(names(rev(sort(table(d$Day_night)))[1]),levels=levels(d$Day_night)),
                     long = mean(d$long),
                     TotalCatch = mean(d$TotalCatch),
                     depth = mean(d$depth),
                     hrs.trawld = rep(seq(min(d$hrs.trawld),max(d$hrs.trawld),length.out=20),2))
    nd=fn.pred.glmm.SE(mod=BIN,dat=nd)
    
    with(subset(nd,Season=='AuWin'),plot(hrs.trawld,pred,ylim=c(0,.015)))
    with(subset(nd,Season=='AuWin'),segments(hrs.trawld,pred+1.96*SE,hrs.trawld,pred-1.96*SE))
    with(subset(nd,Season=='SprSu'),points(hrs.trawld,pred,col=2))
    with(subset(nd,Season=='SprSu'),segments(hrs.trawld,col=2,pred+1.96*SE,hrs.trawld,pred-1.96*SE))
    
    #Depth
    nd <- data.frame(Season = factor("AuWin",levels=levels(d$Season)), 
                     VESSEL = factor(names(rev(sort(table(d$VESSEL)))[1]),levels=levels(d$VESSEL)),
                     Day_night = factor(names(rev(sort(table(d$Day_night)))[1]),levels=levels(d$Day_night)),
                     long = mean(d$long),
                     TotalCatch = mean(d$TotalCatch),
                     depth = seq(min(d$depth),max(d$depth),length.out=20),
                     hrs.trawld = mean(d$hrs.trawld))
    nd=fn.pred.glmm.SE(mod=BIN,dat=nd)
    
    with(nd,plot(depth,pred))
    with(nd,segments(depth,pred+1.96*SE,depth,pred-1.96*SE))
    
    
    #Total catch
    nd <- data.frame(Season = factor("AuWin",levels=levels(d$Season)), 
                     VESSEL = factor(names(rev(sort(table(d$VESSEL)))[1]),levels=levels(d$VESSEL)),
                     Day_night = factor(names(rev(sort(table(d$Day_night)))[1]),levels=levels(d$Day_night)),
                     long = mean(d$long),
                     TotalCatch =seq(min(d$TotalCatch),max(d$TotalCatch),length.out=20) ,
                     depth = mean(d$depth) ,
                     hrs.trawld = mean(d$hrs.trawld))
    nd=fn.pred.glmm.SE(mod=BIN,dat=nd)
    
    with(nd,plot(TotalCatch,pred))
    with(nd,segments(TotalCatch,pred+1.96*SE,TotalCatch,pred-1.96*SE))
    
    
    
  }
  Model.show[[s]]=fn.show(MOD=Model.out[[s]])
  
  
  
  #ind=sapply(d,is.character)  #find all character variables
  
  
  #1.   Simulate rare event data to test effect of SMOTE and check if BRT can retrieve data pattern
  #1.   Simulate rare event data to test effect of SMOTE and check if BRT can retrieve data pattern
  #Simulation-testing of the effect of SMOTE
  Do.sim.test="NO"
  #METRIC="Sens"
  #METRIC='Accuracy'
  METRIC='ROC'
  if(Do.sim.test=="YES")
  {
    #1. Original data where prob(event)~f(Var1,Var2)
    set.seed(1234)
    N=5000
    #N=nrow(Data)
    Var1 = rnorm(N,10,2)           # some continuous variables 
    Var2 = rnorm(N,10,2)
    Var3 = rnorm(N,10,2)
    Var4 = rnorm(N,10,2)
    Var5 = rnorm(N,10,2)
    Factor1 <- factor(sample(as.character(1:12),N,replace=TRUE),levels=1:12)   #some factors
    Factor2 <- factor(sample(letters[3:6],N,replace=TRUE))
    Factor3 <- factor(sample(letters[7:11],N,replace=TRUE))
    Factor4 <- factor(sample(letters[12:13],N,replace=TRUE))
    Factor5 <- factor(sample(letters[14:20],N,replace=TRUE))
    df = data.frame(Var1=Var1,Var2=Var2,Var3=Var3,Var4=Var4,Var5=Var5,
                    Factor1=Factor1,Factor2=Factor2,Factor3=Factor3,
                    Factor4=Factor4,Factor5=Factor5)
    
    #df$z=(3-2*df$Var1) *rnorm(N,1,.2)
    #df$z=with(df,ifelse(!Factor1%in%c(5:9),z*1.5,z)) 
    # linear combination with a bias
    #df$z=(1+1*df$Var1-1*df$Var2+1.5*df$Var3-1.5*df$Var4) *runif(N,.4,1.6)  #balanced
    df$z=(-19+1*df$Var1-1*df$Var2+1.5*df$Var3-1.5*df$Var4) *runif(N,.1,1.9)   #imbalanced
    #df$z=(-9+2*df$Var1-1*df$Var2) *runif(N,.6,1.4)
    #df$z=(9+-2*df$Var1)*runif(N,.6,1.4)
    df$pr = 1/(1+exp(-df$z))         # pass through an inv-logit function
    df$Occurrence = rbinom(N,1,df$pr)      # bernoulli response variable
    df$Occurrence=with(df,ifelse(Occurrence==1,'yes','no'))
    df$Occurrence=factor(df$Occurrence,levels=c('yes','no')) 
    100*prop.table(table(df$Occurrence))
    #splom(~df[, 1:10], groups = df$Occurrence)
    
    
    par(mfcol=c(3,2),mar=c(2,3,1,2),mgp=c(1.8,.7,0))
    with(df,plot(Var1,pr,pch=19))
    with(df,plot(Var2,pr,pch=19))
    with(df,plot(Var3,pr,pch=19))
    with(df,plot(Var4,pr,pch=19))
    with(df,plot(Var5,pr,pch=19))
    
    
    par(mfcol=c(5,2),mar=c(2,3,1,2),mgp=c(1.8,.7,0))
    with(df,plot(cut(Var1,10),Occurrence,ylab="",xlab="var1",yaxt='n'))
    with(df,plot(cut(Var2,10),Occurrence,ylab="",xlab="var2"))
    with(df,plot(cut(Var3,10),Occurrence,ylab="",xlab="var3"))
    with(df,plot(cut(Var4,10),Occurrence,ylab="",xlab="var4"))
    with(df,plot(cut(Var4,10),Occurrence,ylab="",xlab="var5"))
    
    with(df,plot(Factor1,Occurrence,ylab="",xlab="Factor1"))
    with(df,plot(Factor2,Occurrence,ylab="",xlab="Factor2"))
    with(df,plot(Factor3,Occurrence,ylab="",xlab="Factor3"))
    with(df,plot(Factor4,Occurrence,ylab="",xlab="Factor4"))
    with(df,plot(Factor5,Occurrence,ylab="",xlab="Factor5"))
    
    
    df1=df[,-match(c("z","pr"),names(df))]
    
    ind <- createDataPartition(df1$Occurrence,times = 1,p = 0.7,list = FALSE)
    train <- df1[ind,]
    test <- df1[-ind,]
    
    #train.smote <- SMOTE(Occurrence~ ., train, perc.over = 100, perc.under=200)
    # prop.table(table(df1$Occurrence))
    # prop.table(table(train$Occurrence))
    # prop.table(table(test$Occurrence))
    # prop.table(table(train.smote$Occurrence))
    
    
    ctrl <- trainControl(method = "repeatedcv",
                         number = 10,     #split data in ten ways to build 10 models
                         repeats = 3,     #repeat 10 fold cross validation 3 times
                         classProbs = TRUE,     #needed to score models using AUC
                         summaryFunction = twoClassSummary,
                         search = "grid")     #find optimal collection of parameters using default grid
    
    fit.glm <- train(Occurrence~.,data=train, method='glm', metric = METRIC,trControl=ctrl)
    #fit.gbm=train(Occurrence~.,data=train,method="gbm", metric = METRIC,verbose=F,trControl=ctrl)  
    
    
    ctrl$sampling='smote'
    fit.glm.smoted <- train(Occurrence~.,data=train, method='glm', metric = METRIC,trControl=ctrl)
    #fit.gbm.smoted=train(Occurrence~.,data=train,method="gbm", metric = METRIC,verbose=F,trControl=ctrl)  
    
    
    #df$z=(-19+1*df$Var1-1*df$Var2+1.5*df$Var3-1.5*df$Var4) *runif(N,.1,1.9)
    exp(coef(fit.glm$finalModel))
    anova(fit.glm$finalModel,test="Chisq")
    anova(fit.glm.smoted$finalModel,test="Chisq")
    
    
    test$pred.glm=predict(fit.glm, test, type="raw")
    #test$pred.gbm=predict(fit.gbm, test, type="raw")
    
    test$pred.glm.smote=predict(fit.glm.smoted, test, type="raw")
    #test$pred.gbm.smote=predict(fit.gbm.smoted, test, type="raw")
    
    
    TAB.glm=with(test,table(Occurrence,pred.glm))
    #with(test,table(Occurrence,pred.gbm))
    
    TAB.glm.smote=with(test,table(Occurrence,pred.glm.smote))
    #with(test,table(Occurrence,pred.gbm.smote))
    
    #Plot 
    TAB.1=with(test,table(Occurrence))
    TAB.2=with(test,table(pred.glm))
    TAB.3=with(test,table(pred.glm.smote))
    
    fn.fig("Simulation",2400,2400)
    b=barplot(GLM.Model.out[[1]]$Tab[,2],names.arg=c("Observed","Binomial","ZIP","Hurdle"),cex.names=1.25)
    text(b[,1],GLM.Model.out[[1]]$Tab[,2]+5,c("",paste(100*GLM.Model.out[[1]]$Prop.corrct,"%")))
    box()
    legend("top",Species[s],cex=1.25,bty='n')
    
    mtext("Number of interactions",2,las=3,line=-.25,outer=T,cex=1.5)
    dev.off()
    
    
    
    #Alternative
    N=30000
    percent.pos.event=1
    odds.AuWin=10
    odds.SprSum=1
    df=data.frame(
      Season=factor(sample(c("AuWin","SprSum"),N,replace=T)),
      VESSEL=factor(sample(c(paste("Ves",1:10,sep=".")),N,replace=T)),
      Day_night=factor(sample(c("Day","Night"),N,replace=T)),
      long=rnorm(N,mean(Model.d$long),sd(Model.d$long)),
      depth=rnorm(N,mean(Model.d$depth),sd(Model.d$depth)),
      TotalCatch=rnorm(N,mean(Model.d$TotalCatch),sd(Model.d$TotalCatch)/5),
      Ann.eff=rnorm(N,mean(Model.d$Ann.eff),sd(Model.d$Ann.eff)/3),
      hrs.trawld=rnorm(N,mean(Model.d$hrs.trawld),sd(Model.d$hrs.trawld)/2)
    )
    df$Occurrence=0
    ind=sample(1:nrow(df),percent.pos.event/100*N,replace=F)
    df$Occurrence[ind]=1
    df$Season[ind]=sample(c("AuWin","SprSum"),length(ind),replace=T,prob=c(odds.AuWin,odds.SprSum))
    #df$Occurrence[ind]=ifelse(df$Season[ind]=="AuWin",sample(0:1,length(ind),replace=T,prob=c(.01,Cond.prob.AuWin)),
    #                         sample(0:1,length(ind),replace=T,prob=c(.5,Cond.prob.SprSum)))
    #df$Occurrence[ind]=ifelse(df$hrs.trawld[ind]<2,0,df$Occurrence[ind])
    100*prop.table(table(df$Occurrence,df$Season))
    
    
    ind <- createDataPartition(df$Occurrence,times = 1,p = 0.7,list = FALSE)
    train <- df[ind,]
    test <- df[-ind,]
    
    train1=train
    train1$Occurrence=ifelse(train1$Occurrence==1,"yes","no")
    train1$Occurrence=factor(train1$Occurrence,levels=c("yes","no"))
    train.smote <- SMOTE(Occurrence~ ., train1, perc.over = 100, perc.under=500)
    train.smote$Occurrence=ifelse(train.smote$Occurrence=="no",0,1)
    prop.table(table(train.smote$Occurrence))
    
    
    fit.glm <- glm(Occurrence~Season+VESSEL+Day_night+long+depth+TotalCatch+Ann.eff+hrs.trawld,data=train, family=binomial)
    
    fit.glm.smoted <- glm(Occurrence~Season+VESSEL+Day_night+long+depth+TotalCatch+Ann.eff+hrs.trawld,data=train.smote, family=binomial)
    
    exp(coef(fit.glm))
    exp(coef(fit.glm.smoted))
    # 
    anova(fit.glm,test="Chisq")
    anova(fit.glm.smoted,test="Chisq")
    
    
    
    ###
    biCatTemp <- 0
    for (i in 3:length(dummy.coef(fit.glm)))
    {
      # Compute mean coef. for each categorical variable and populate "biCatTemp" vector
      biCatTemp[i-1]<- mean(dummy.coef(GLMbi)[[i]])
    }
    # Sum mean coefs. of categorical variables
    biCatFact <- sum(biCatTemp) # This is the constant factor for exp. vars. to use in Eq. 9 (see Maunder and Punt, 2004)
    
    exp(dummy.coef(fit.glm)$Season)/(1+exp(dummy.coef(fit.glm)$Season))
    exp(dummy.coef(fit.glm.smoted)$Season)/(1+exp(dummy.coef(fit.glm.smoted)$Season))
    
    # Compute the sum of coefficients (intercept + year coeffs + constant factors)
    BiSumCoeff <- dummy.coef(GLMbi)[[1]] + dummy.coef(GLMbi)[[2]] + biCatFact
    # Extract the year effects by applying the inverse of the logit link (Eq. 9 in Maunder and Punt, 2004) 
    BiIndex <- exp(BiSumCoeff)/(1+exp(BiSumCoeff))
    
    
    
    ###
    
    
    
    
    Eff.seq=seq(min(test$hrs.trawld),max(test$hrs.trawld),length.out=20)
    Eff.seq=c(Eff.seq,Eff.seq)
    dummy=test[1:length(Eff.seq),]
    dummy[,-match(c("Season","hrs.trawld"),names(dummy))]=dummy[1,-match(c("Season","hrs.trawld"),names(dummy))]
    dummy$hrs.trawld=mean(test$hrs.trawld)
    # dummy$hrs.trawld=Eff.seq
    dummy$Season=rep(c("AuWin","SprSum"),each=length(Eff.seq)/2)
    dummy=subset(dummy,select=-Occurrence)
    
    dummy$pred.glm=predict(fit.glm, dummy, type="response")
    dummy$pred.glm.smote=predict(fit.glm.smoted, dummy, type="response")
    
    plot(dummy$pred.glm.smote[1:20])
    
    test$pred.glm=predict(fit.glm, test, type="response")
    test$pred.glm.smote=predict(fit.glm.smoted, test, type="response")
    
    
    auc.glm=roc(test$Occurrence, test$pred.glm)  #.5 is random model; 1 is perfect model
    auc.glm.smote=roc(test$Occurrence, test$pred.glm.smote)
    
    plot(auc.glm)
    lines(auc.glm.smote,col=2)
    legend("bottomright",paste(c("glm","glm.smote"),round(c(auc.glm$auc[1],auc.glm.smote$auc[1]),2),sep=" "),bty='n',lty=1,col=1:2)
    
    
    
    test$pred.glm=ifelse(test$pred.glm>=.5,1,0)
    test$pred.glm.smote=ifelse(test$pred.glm.smote>=.5,1,0)
    # 
    # table(test$Occurrence)
    # a=subset(test,Occurrence==1)
    # sum(with(a,Occurrence==pred.glm))
    # sum(with(a,Occurrence==pred.glm.smote))
    # 
    # 
    # 
    # TAB.glm=with(test,table(Occurrence,pred.glm))
    # #with(test,table(Occurrence,pred.gbm))
    # 
    # TAB.glm.smote=with(test,table(Occurrence,pred.glm.smote))
    # #with(test,table(Occurrence,pred.gbm.smote))
    # 
    # #Plot 
    TAB.1=with(test,table(Occurrence))
    TAB.2=with(test,table(pred.glm))
    TAB.3=with(test,table(pred.glm.smote))
    
    TAB.1;TAB.2;TAB.3
  }
  
  #2.   fit glms to raw data
  fn.model.glm.all=function(species,normalised,FRML)
  {
    set.seed(1234)
    
    #select if doing analyses on normalised data or not
    if(normalised=="NO")d=Model.d[,c(species,TERMS)]
    if(normalised=="YES")d=Model.d.scaled[,c(species,TERMS)]
    
    #set character variables to factors
    idd=which(!names(d)%in%c(COVS,species))
    for(f in 1:length(idd))d[,idd[f]]=as.factor(d[,idd[f]])
    d$Occurrence=d[,species]
    d=d[,-match(species,names(d))]
    d.bin=d
    d.bin$Occurrence=with(d.bin,ifelse(Occurrence>0,1,0))
    ind <- createDataPartition(d$Occurrence,times = 1,p = 0.7,list = FALSE)
    
    BIN=model<- glm(Occurrence ~ Season+VESSEL+Day_night+long+depth+TotalCatch+offset(hrs.trawld), data = d.bin[ind,], family="binomial", maxit=500)
    ZIP<-zeroinfl(FRML, data = d[ind,], dist = "poisson")
    HURDLE<-hurdle(FRML, data = d[ind,], dist = "poisson")
    
    test.bin=d.bin[-ind,]
    test=d[-ind,]
    BIN.pred=predict(BIN,newdata=test.bin,type='response')
    ZIP.pred=predict(ZIP,newdata=d[-ind,-match("Occurrence",names(d))],type='prob')
    HURDLE.pred=predict(HURDLE,newdata=d[-ind,-match("Occurrence",names(d))],type='prob')
    
    ZIP.pred=rowSums(ZIP.pred[,2:ncol(ZIP.pred)])
    HURDLE.pred=rowSums(HURDLE.pred[,2:ncol(HURDLE.pred)])
    
    BIN.pred=ifelse(BIN.pred>=.5,1,0)
    ZIP.pred=ifelse(ZIP.pred>=.5,1,0)
    HURDLE.pred=ifelse(HURDLE.pred>=.5,1,0)
    
    
    test.bin$Pred=BIN.pred
    Prop.correct.bin=with(subset(test.bin,Occurrence==1),sum(Occurrence==Pred))
    test$Pred.zip=ZIP.pred
    test$Pred.hur=HURDLE.pred
    Prop.correct.zip=with(subset(test,Occurrence==1),sum(Occurrence==Pred.zip))
    Prop.correct.hur=with(subset(test,Occurrence==1),sum(Occurrence==Pred.hur))
    
    Prop.corrct=c(Prop.correct.bin,Prop.correct.zip,Prop.correct.hur)
    
    Tab.obs=as.matrix(t(table(d.bin[-ind,]$Occurrence)))
    Tab.pred.glm=as.matrix(t(table(BIN.pred)))
    Tab.pred.ZIP=as.matrix(t(table(ZIP.pred)))
    Tab.pred.Hur=as.matrix(t(table(HURDLE.pred)))
    NMS=c("0","1")
    id=which(!NMS%in%colnames(Tab.pred.glm))
    if(length(id>0)) Tab.pred.glm=c(Tab.pred.glm,0)
    id=which(!NMS%in%colnames(Tab.pred.ZIP))
    if(length(id>0)) Tab.pred.ZIP=c(Tab.pred.ZIP,0)
    id=which(!NMS%in%colnames(Tab.pred.Hur))
    if(length(id>0)) Tab.pred.Hur=c(Tab.pred.Hur,0)
    
    Tab=rbind(Tab.obs,Tab.pred.glm,Tab.pred.ZIP,Tab.pred.Hur)
    return(list(Prop.corrct=Prop.corrct,Tab=Tab,BIN=BIN, ZIP=ZIP,HURDLE=HURDLE))
  }
  
  GLM.Model.out=vector('list',length(Species))
  names(GLM.Model.out)=Species
  FRML= GLM.Model.out
  FRML[[1]]=as.formula(Occurrence ~ Season+VESSEL+Day_night+depth+hrs.trawld|
                         Season+VESSEL+Day_night+long+depth+hrs.trawld)
  FRML[[2]]=as.formula(Occurrence ~ Season+VESSEL+Day_night+long+hrs.trawld|
                         Season+VESSEL+Day_night+long+depth+TotalCatch+hrs.trawld)
  
  # FRML[[1]]=as.formula(Occurrence ~ Season+VESSEL+Day_night+depth+offset(hrs.trawld)|
  #                        Season+VESSEL+Day_night+long+depth+TotalCatch+offset(hrs.trawld))
  # FRML[[2]]=as.formula(Occurrence ~ Season+VESSEL+Day_night+long+offset(hrs.trawld)|
  #                        Season+VESSEL+Day_night+long+depth+TotalCatch+offset(hrs.trawld))
  system.time({
    for(s in 1:length(Species))GLM.Model.out[[s]]=fn.model.glm.all(species=Species[s],normalised="NO",FRML[[s]]) 
  })
  
  
  #3.   Fit machine learning models to data
  fn.model=function(species,normalised,Optimised)
  {
    set.seed(1234)
    
    #select if doing analyses on normalised data or not
    if(normalised=="NO")d=Model.d[,c(species,TERMS)]
    if(normalised=="YES")d=Model.d.scaled[,c(species,TERMS)]
    
    #set character variables to factors
    idd=which(!names(d)%in%c(COVS,species))
    for(f in 1:length(idd))d[,idd[f]]=as.factor(d[,idd[f]])
    
    d$Occurrence=d[,species]
    d=d[,-match(species,names(d))]
    
    #set response variable to presence/absence
    d$Occurrence=ifelse(d$Occurrence>0,1,0)
    d$Occurrence=factor(ifelse(d$Occurrence==1,"yes","no"),levels=c("yes","no"))    #have presence as first levels
    
    #partition data
    ind <- createDataPartition(d$Occurrence,times = 1,p = 0.7,list = FALSE)
    train <- d[ind,]
    test <- d[-ind,]
    # 100*prop.table(table(d$Occurrence))     #percent yes or no
    # 100*prop.table(table(train$Occurrence))
    # 100*prop.table(table(test$Occurrence))
    
    
    #Oversample rare event before model
    #train.smote <- SMOTE(Occurrence~ ., train, perc.over = 100, perc.under=200)
    #prop.table(table(train.smote$Occurrence))
    
    
    #register cluster to train in parallel (only works for Xgboost)
    cl <-makeCluster(CORES, type = "SOCK") 
    registerDoSNOW(cl)
    
    #Control
    ctrl <- trainControl(method = "repeatedcv",
                         number = 10,     #split data in ten ways to build 10 models
                         repeats = 3,     #repeat 10 fold cross validation 3 times
                         classProbs = TRUE,     #needed to score models using AUC
                         summaryFunction = twoClassSummary,
                         search = "grid")     #find optimal collection of parameters using default grid
    
    #we can compare them using caret's resamples() command. Given that we asked caret to
    # calculate the class probabilities ("classprobs = TRUE") and use 
    # the twoClassFunction summary function ("summaryFunction = twoClassSummary"), 
    # this gives us the sensitivity, specificity, and ROC values for every iteration of 
    # the $finalModel. (For logistic regression, this is just the 10 resamples, but for 
    # the others, it refers to the iterations with the tuning parameters that were 
    # ultimately selected as best.)
    
    #fit original model
    # fit=train(Occurrence~.,data=train,method="gbm",distribution="bernoulli",
    #           metric="ROC",   #will maximise ROC, the area under the curve
    #           verbose=F,trControl=ctrl)
    
    
    #Over and under sampling 
    #ctrl$sampling <- "smote"   
    
    #fit smote model
    #ctrl$seeds <- fit$control$seeds #start from same seed
    # fit.smote=train(Occurrence~.,data=train.smote,method="gbm",distribution="bernoulli",
    #                 metric="ROC",verbose=F,trControl=ctrl)
    
    #8. glm
    fit.glm.smote <- train(Occurrence~.,data=train, method='glm', metric = Optimised,
                           trControl=ctrl)
    train.inflated=train
    train.inflated$Occurrence=ifelse(train.inflated$Occurrence=="no",0,1)
    ZIP<-zeroinfl(FRML[[s]], data = train.inflated, dist = "poisson")
    HURDLE<-hurdle(FRML[[s]], data = train.inflated, dist = "poisson")
    
    
    #9. lda
    fit.lda.smote <- train(Occurrence~.,data=train, method='lda', metric = Optimised,
                           trControl=ctrl)
    
    #10. Support Vector Machine
    fit.SVM.smote <- train(Occurrence~.,data=train, method='svmRadial', metric = Optimised,
                           trControl=ctrl)
    
    #11. Naive Bayes
    fit.NB.smote <- train(Occurrence~.,data=train, method='nb', metric = Optimised,
                          trControl=ctrl)
    
    #12. RPART
    fit.rpart.smote <- train(Occurrence~.,data=train, method='rpart', metric = Optimised,
                             trControl=ctrl) 
    
    #13. Random forest
    fit.rf.smote <- train(Occurrence~.,data=train, method='rf', metric = Optimised,
                          trControl=ctrl)
    
    #14. gbm
    fit.gbm.smote <- train(Occurrence~.,data=train, method='gbm', metric = Optimised,
                           trControl=ctrl)
    
    #Xgboost
    # tune.grid <- expand.grid(eta = c(0.05, 0.075, 0.1),     #eta is learning rate
    #                          nrounds = c(50, 75, 100),      #nrounds is number of trees
    #                          max_depth = 6:8,
    #                          min_child_weight = c(2.0, 2.25, 2.5),
    #                          colsample_bytree = c(0.3, 0.4, 0.5),
    #                          gamma = 0,
    #                          subsample = 1)     #create all permutations combos
    # fit.xgboost <- train(Occurrence~.,data=train,
    #                   method = "xgbTree",
    #                   tuneGrid = tune.grid,
    #                   trControl = ctrl)
    
    
    
    #clear parameters for parallel computation
    stopCluster(cl)
    registerDoSEQ()
    remove(cl)
    
    return(list(test=test,train=train,
                fit.GLM.smote=fit.glm.smote,fit.lda.smote=fit.lda.smote,fit.SVM.smote=fit.SVM.smote,
                fit.NB.smote=fit.NB.smote,fit.RPART.smote=fit.rpart.smote,fit.RF.smote=fit.rf.smote,
                fit.GBM.smote=fit.gbm.smote))
  }
  
  Model.out=vector('list',length(Species))
  names(Model.out)=Species
  Best.MDL=Model.out
  METRIC="Sens"
  #METRIC="ROC"
  system.time({
    for(s in 1:length(Species))Model.out[[s]]=fn.model(species=Species[s],normalised="NO",Optimised=METRIC) 
  })
  
  
  ####
  #Select best algorithm
  source("C:/Matias/Analyses/SOURCE_SCRIPTS/Smart_par.R")
  fun.pred.auc=function(moDl,NMS,test)
  {
    pred.class=predict(moDl,test,type='raw',n.trees=moDl$bestTune$n.trees)  #predict yes or not
    pred.prob=predict(moDl,test,type='prob',n.trees=moDl$bestTune$n.trees)#pred probability
    #auc
    auc=roc(test$Occurrence, pred.prob[,match('yes',names(pred.prob))])  #.5 is random model; 1 is perfect model
    
    Preds=cbind(pred.class,pred.prob)
    names(Preds)=paste(NMS,names(Preds),sep="_")
    
    return(list(Preds=Preds,auc=auc,ConfMat=confusionMatrix(pred.class,test$Occurrence)))
  }
  best.mdl=function(modl)
  {
    Fits=names(modl)[-match(c("test","train"),names(modl))]
    
    #term importance
    VARIMP=lapply(modl[Fits],varImp)
    smart.par(n.plots=length(VARIMP),MAR=c(4,4,1,1),OMA=rep(1,4),MGP=c(2.5,.7,0))
    for(l in 1:length(VARIMP))
    {
      x=VARIMP[[l]][[1]]
      par(las=1)
      barplot(x[,1],main=names(VARIMP)[l],names.arg=rownames(x),horiz =T)
    }
    plot(1:1,ylab="",xlab="",xaxt='n',yaxt='n',col="transparent",main=Species[s])  
    box(col="white")
    
    
    # Table comparison
    results <-resamples(modl[Fits])
    
    
    Table=summary(results)
    
    # boxplot comparison
    #notes:
    #   Sensitivity = proportion of positive catch predicted correctly
    #   Specificity = proportion of negative catch predicted correctly
    bwplot(results)
    
    # Dot-plot comparison
    #dotplot(results)
    
    #Get accuracy and out of sample error
    #not relevant in this case due to high number of zero catches
    # fn.pred.Acc=function(mod)
    # {
    #   Pred <- predict(mod, modl$test, type="raw")
    #   Accuracy <- sum(Pred == modl$test$Occurrence) / length(Pred)
    #   
    #   return(list(Confusion=confusionMatrix(Pred,modl$test$Occurrence, positive = "yes"),
    #         Table=data.frame(Model=names(modl)[i],Accuracy=Accuracy,OutOfSampleError = 1-Accuracy)))
    # }
    # Store=vector('list',length(Fits))
    # for(i in 1:length(Store))Store[[i]]=fn.pred.Acc(mod=modl[[match(Fits[i],names(modl))]])
    # Acc.table=do.call(rbind,Store$Table)
    # Acc.table=Acc.table[order(Acc.table$Accuracy),]
    
    AUCs=vector('list',length(Fits))
    names(AUCs)=Fits
    for(f in 1:length(Fits))
    {
      AUCs[[f]]=fun.pred.auc(moDl=modl[[match(Fits[f],names(modl))]],NMS=Fits[f],test=modl[[match('test',names(modl))]])   
    }
    
    return(list(Table=Table,AUC=AUCs))
  }
  
  for(s in 1:length(Species))Best.MDL[[s]]=best.mdl(modl=Model.out[[s]])
  
  #confusion matrix
  
  #AUC
  smart.par(n.plots=length(Species),MAR=c(4,4,1,1),OMA=rep(1,4),MGP=c(2.5,.7,0))
  for(s in 1:length(Species))
  {
    plot(Best.MDL[[s]]$AUC[[1]]$auc)
    aucs=rep(NA,length(Best.MDL[[s]]$AUC))
    aucs[1]=Best.MDL[[s]]$AUC[[1]]$auc$auc[1]
    CL=2:length(Best.MDL[[s]]$AUC)
    for(l in CL)
    {
      lines(Best.MDL[[s]]$AUC[[l]]$auc,col=CL[l-1])
      aucs[l]=Best.MDL[[s]]$AUC[[l]]$auc$auc[1]
    }
    nombres=names(Best.MDL[[s]]$AUC)  
    nombres=read.table(text = nombres, sep = ".", as.is = TRUE)$V2
    legend('bottomright',paste(nombres,round(aucs,3)),lty=1,col=c(1,CL),bty='n',title="AUC",cex=0.75)
    legend("topleft",Species[s],bty='n')
  }
  
  
  #Show marginal effects
  InvLgt=function(x) exp(x)/(1+exp(x))
  prob.to.odd=function(p) p/(1-p)
  odd.to.prob=function(odds) odds/(1+odds)
  odds.ratio=function(odd1,odd2) odd1/odd2
  odds.ratio_probs=function(prob1,prob2) (prob1/(1-prob1))/(prob2/(1-prob2))
  library(pdp)
  library(doParallel)  #for parallel processing of pdf
  
  
  cl <- makeCluster(CORES) 
  registerDoParallel(cl) # register the parallel backend
  CL=c("black","grey80")
  for(t in 1:length(TERMS))
  {
    var=subset(Model.out[[1]]$train,select=TERMS[t])
    if(is.numeric(var[,1]))plot(min(var[,1]):max(var[,1]),min(var[,1]):max(var[,1]),col='transparent',xaxt='n',ylim=c(0,1),ylab="",xlab=TERMS[t])
    if(is.factor(var[,1]))plot(1:length(unique(var[,1])),1:length(unique(var[,1])),col='transparent',xaxt='n',ylim=c(0,1),ylab="",xlab=TERMS[t])
    for(s in 1:length(Species))
    {
      Train=Model.out[[s]]$train
      modl=Model.out[[s]]$fit.GBM.smote
      var=subset(Train,select=TERMS[t])
      if(is.factor(var[,1])) Pred.Grid=data.frame(factor(levels(var[,1]),levels=levels(var[,1])))
      if(is.numeric(var[,1])) Pred.Grid=data.frame(seq(min(var[,1]),max(var[,1]),length.out=20))
      names(Pred.Grid)=TERMS[t]
      pdp <- partial(modl, pred.var =TERMS[t],pred.grid =Pred.Grid,
                     train=subset(Train,select=-Occurrence),parallel = TRUE)
      d=as.data.frame.matrix(pdp)
      
      
      d$Odds=exp(d$yhat)
      d$Prob=InvLgt(d$yhat)
      d$Odds_ratio=d$Odds/d$Odds[1]
      if(TERMS[t]=="VESSEL") d[,1]=1:nrow(d)
      
      if(is.numeric(var[,1]))points(d[,1],d$Prob,pch=19,col=CL[s],cex=1.25)
      if(is.factor(var[,1])) points(1:nrow(d),d$Prob,pch=19,col=CL[s],cex=1.25)
      
      
      #barplot(d$Prob,xlab=TERMS[t],ylab="",pch=19,names.arg=d[,1],ylim=c(0,1))
    }
    axis(1,d[,1])
    box()
    
  }
  stopCluster(cl)
  
  
  #Predictions from binomial and Zero inflated models
  fn.fig("GLM_prediction_interaction",2400,2400)
  smart.par(n.plots=length(Species),MAR=c(1,2,1,1),OMA=rep(1,4),MGP=c(2.5,.7,0))
  for(s in 1:length(Species))
  {
    b=barplot(GLM.Model.out[[1]]$Tab[,2],names.arg=c("Observed","Binomial","ZIP","Hurdle"),cex.names=1.25)
    text(b[,1],GLM.Model.out[[1]]$Tab[,2]+5,c("",paste(100*GLM.Model.out[[1]]$Prop.corrct,"%")))
    box()
    legend("top",Species[s],cex=1.25,bty='n')
  }
  mtext("Number of interactions",2,las=3,line=-.25,outer=T,cex=1.5)
  dev.off()
  
  
}

