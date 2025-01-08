library(stringr)
library(RMark)
library(dplyr)
library(lme4)
library(ggplot2)
library(MuMIn)
library(emmeans)
library(ggeffects)
source("_rMarkAbund_FUN_v013.013.R")
setwd("C:/Users/jearl/OneDrive - Michigan State University/UGA/Overall Density/Resubmission analysis")
source("C:/Users/jearl/OneDrive - University of Georgia/Baiting All/Redo uptake data/outputfunction.R")
bt<-read.table("./finalcapturedata.txt", header=TRUE, colClasses = c("factor", "factor", "factor", "factor", "character", "numeric"))
area<-read.csv("effarealist5Aug.csv")

#STEP 1- Use model selection to determine most-supported population parameters for each species
raccoon<- subset(bt, Species=="Raccoon")
opossum<- subset(bt, Species=="Opossum")
time.intervals<-c(0,0,0,0,0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,4,
                  0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0)
raccoon2<-raccoon[!(raccoon$Grid=="P1" | raccoon$Grid=="P7"),] #removed these that weren't trapped the entire period

multigroup.models=function () {aa.process=(process.data(raccoon2, model="RDHuggins", 
                                                        time.intervals=time.intervals, groups = c("Sex", "Habitat"))) #Building the model combinations
aa.ddl=make.design.data(aa.process) 
S.const=list(formula=~1) #constant survival
S.time=list(formula=~time) #time-varying survival 
S.habitat=list(formula=~Habitat) #Habitat-varying survival
S.sex=list(formula=~Sex) #Sex-varying survival
p.sex=list(formula=~Sex, share=TRUE) #sex-varying capture probability, recapture set the same
p.habitat=list(formula=~Habitat, share=TRUE) #habitat-varying capture probability, recapture set the same
p.const=list(formula=~1, share=TRUE) #constant capture probability, recapture set the same
p.time=list(formula=~session, share=TRUE) #time-varying capture probability, recapture set the same
p.sessionhab=list(formula=~session+Habitat, share=TRUE) #session and habitat-varying capture probability, recapture set the same
GammaDoublePrime.time=list(formula=~time, share=TRUE) #random emigration
GammaDoublePrime.const=list(formula=~1, share=TRUE) #constant emigration
cml=create.model.list("RDHuggins")
results=mark.wrapper(cml, data=aa.process, ddl=aa.ddl, adjust=TRUE)
return(results)
}

multiresults2=multigroup.models()
multiresults2 #displays model table
multiresults2[[28]] #gets info from top model

#Opossum- repeat of above with fewer models
opossum2<-opossum[!(opossum$Grid=="P1" | opossum$Grid=="P7"),]
multigroup.modelsopp=function () {aa.process=(process.data(opossum2, model="RDHuggins", time.intervals=time.intervals, 
                                                           groups = c("Sex", "Habitat")))
aa.ddl=make.design.data(aa.process)
S.const=list(formula=~1) #constant survival
S.time=list(formula=~time) #time-varying survival
S.sex=list(formula=~Sex) #sex-varying survival
p.sex=list(formula=~Sex, share=TRUE) #sex-varying capture probability, recapture set the same
p.habitat=list(formula=~Habitat, share=TRUE, link="logit") #habitat-varying capture probability, recapture set the same
p.const=list(formula=~1, share=TRUE) #constant capture probability, recapture set the same
p.time=list(formula=~session, share=TRUE) #time-varying capture probability, recapture set the same
p.sexhabitat=list(formula=~Sex+Habitat, share=TRUE, link="logit") #sex and habitat-varying capture probability, recapture set the same
GammaDoublePrime.const=list(formula=~1, share=TRUE, link="logit") #constant emigration
cml=create.model.list("RDHuggins")
results=mark.wrapper(cml, data=aa.process, ddl=aa.ddl, adjust=TRUE)
return(results)
}

multiresultsopp=multigroup.modelsopp()
multiresultsopp
multiresultsopp[[5]]

#STEP 2- Use the population parameters from the top model for each species and apply these back to each grid and sex 
#to derive sex- and grid-specific abundances
btm=process.data(raccoon,model="RDHuggins",groups=c("Sex", "Grid"), 
                 time.intervals = time.intervals) 
dipper=make.design.data(btm)
dipper$p$Habitat<- substring(dipper$p$Grid, 1, 1)
dipper$c$Habitat<- substring(dipper$c$Grid, 1, 1)
dipper$p<-dipper$p %>% mutate(Habitat =
                                                   case_when(
                                                     Habitat=="B"  ~ "Bottomland",
                                                     Habitat=="W"~ "Wetland",
                                                     Habitat=="R"  ~ "Riparian",
                                                     Habitat=="P"~ "Pine"))
dipper$c<-dipper$c %>% mutate(Habitat =
                                        case_when(
                                          Habitat=="B"  ~ "Bottomland",
                                          Habitat=="W"~ "Wetland",
                                          Habitat=="R"  ~ "Riparian",
                                          Habitat=="P"~ "Pine"))

dipper$S<- dipper$S %>% mutate(fix =
                                         case_when(
                                           Sex=="Male"  ~ 0.77622,
                                           Sex=="Female"~ 0.85752))

dipper$p<- dipper$p %>% mutate(fix =
                                         case_when(
                                           Habitat== "Bottomland" & session==1~ 0.0633,
                                           Habitat== "Bottomland" &  session==5 ~0.0707,
                                           Habitat== "Bottomland" & session==9 ~0.1243,
                                           Habitat== "Bottomland" & session==13 ~0.0469,
                                           Habitat== "Bottomland" & session==15 ~0.03166,
                                           Habitat== "Bottomland" & session==16~0.0467,
                                           Habitat== "Bottomland" & session==18 ~0.0212,
                                           Habitat== "Bottomland" & session==19 ~0.1418,
                                             Habitat== "Pine" & session==1~ 0.0490,
                                             Habitat== "Pine" & session==5~ 0.0547,
                                             Habitat== "Pine" & session==9 ~ 0.0975,
                                             Habitat== "Pine" & session==13~ 0.0361,
                                             Habitat== "Pine" & session==15  ~ 0.0242,
                                             Habitat== "Pine" & session==16~ 0.0360,
                                             Habitat== "Pine" & session==18  ~ 0.0162,
                                             Habitat== "Pine" & session==19~ 0.1119,
                                            Habitat== "Riparian" & session==1~ 0.0485,
                                           Habitat== "Riparian" & session==5~ 0.0543,
                                           Habitat== "Riparian" & session==9 ~ 0.0968,
                                           Habitat== "Riparian" & session==13~ 0.0357,
                                           Habitat== "Riparian" & session==15  ~ 0.024,
                                           Habitat== "Riparian" & session==16~ 0.03563,
                                           Habitat== "Riparian" & session==18  ~ 0.016,
                                           Habitat== "Riparian" & session==19~ 0.111,
                                          Habitat== "Wetland" & session==1~ 0.0318,
                                          Habitat== "Wetland" & session==5~ 0.0356,
                                          Habitat== "Wetland" &  session==9 ~ 0.0645,
                                          Habitat== "Wetland" &  session==13~ 0.0233,
                                          Habitat== "Wetland" & session==15  ~ 0.0156,
                                          Habitat== "Wetland" &  session==16~ 0.0232,
                                          Habitat== "Wetland" & session==18  ~ 0.0104,
                                          Habitat== "Wetland" &  session==19~ 0.0742))


dipper$c<- dipper$c %>% mutate(fix =
                                 case_when(
                                   Habitat== "Bottomland" & session==1~ 0.0633,
                                   Habitat== "Bottomland" &  session==5 ~0.0707,
                                   Habitat== "Bottomland" & session==9 ~0.1243,
                                   Habitat== "Bottomland" & session==13 ~0.0469,
                                   Habitat== "Bottomland" & session==15 ~0.03166,
                                   Habitat== "Bottomland" & session==16~0.0467,
                                   Habitat== "Bottomland" & session==18 ~0.0212,
                                   Habitat== "Bottomland" & session==19 ~0.1418,
                                   Habitat== "Pine" & session==1~ 0.0490,
                                   Habitat== "Pine" & session==5~ 0.0547,
                                   Habitat== "Pine" & session==9 ~ 0.0975,
                                   Habitat== "Pine" & session==13~ 0.0361,
                                   Habitat== "Pine" & session==15  ~ 0.0242,
                                   Habitat== "Pine" & session==16~ 0.0360,
                                   Habitat== "Pine" & session==18  ~ 0.0162,
                                   Habitat== "Pine" & session==19~ 0.1119,
                                   Habitat== "Riparian" & session==1~ 0.0485,
                                   Habitat== "Riparian" & session==5~ 0.0543,
                                   Habitat== "Riparian" & session==9 ~ 0.0968,
                                   Habitat== "Riparian" & session==13~ 0.0357,
                                   Habitat== "Riparian" & session==15  ~ 0.024,
                                   Habitat== "Riparian" & session==16~ 0.03563,
                                   Habitat== "Riparian" & session==18  ~ 0.016,
                                   Habitat== "Riparian" & session==19~ 0.111,
                                   Habitat== "Wetland" & session==1~ 0.0318,
                                   Habitat== "Wetland" & session==5~ 0.0356,
                                   Habitat== "Wetland" &  session==9 ~ 0.0645,
                                   Habitat== "Wetland" &  session==13~ 0.0233,
                                   Habitat== "Wetland" & session==15  ~ 0.0156,
                                   Habitat== "Wetland" &  session==16~ 0.0232,
                                   Habitat== "Wetland" & session==18  ~ 0.0104,
                                   Habitat== "Wetland" &  session==19~ 0.0742))
dipper$GammaDoublePrime <-dipper$GammaDoublePrime %>% mutate(fix =
                                                                       case_when(
                                                                         time==1~ 0.1215,
                                                                         time==5~ 0.4361,
                                                                         time==9~ 0.7212,
                                                                         time==13~ 0.863,
                                                                         time==15~0.3442,
                                                                         time==16~0.5808
                                                                        ))
                                                                                                                                              

dipper$GammaPrime <-dipper$GammaPrime %>% mutate(fix =case_when(
  time==5~ 0.4361,
  time==9~ 0.7212,
  time==13~ 0.863,
  time==15~0.3442,
  time==16~0.5808,
  time==18~0.8051
))

racctotal<-mark(btm, dipper)
v<-data.frame(racctotal$group.labels)
v2<-data.frame(Group=v[rep(seq_len(nrow(v)), each = 8), ])
dates<-data.frame(session=c("Spring 2017", "Spring 2018", 
                            "Spring 2019", "Spring 2020",
                            "Fall 2020", "Spring 2021",
                            "Fall 2021", "Spring 2022"))
Seasons<-data.frame(Session=rep(dates$session, times=47))
EstimatesRacc<-bind_cols(racctotal$results$derived, v2, Seasons)
EstimatesRacc$Species<-paste("Raccoon")

##opossum

br2=process.data(opossum,model="RDHuggins", groups=c("Sex", "Grid"), time.intervals = time.intervals)
dipper.ddl2=make.design.data(br2)
dipper.ddl2$S<- dipper.ddl2$S %>% mutate(fix =
                                           case_when(
                                             Sex=="Male"  ~ 0.448,
                                             Sex=="Female"~ 0.448))


dipper.ddl2$p<- dipper.ddl2$p %>% mutate(fix =
                                           case_when(
                                             session==1~ 0.0348,
                                             session==5~ 0.0244,
                                             session==9 ~ 0.0977,
                                             session==13~ 0.0661,
                                             session==15  ~ 0.0736,
                                             session==16~ 0.1023,
                                             session==18  ~ 0.0968,
                                             session==19~ 0.10734))

dipper.ddl2$c<- dipper.ddl2$c%>% mutate(fix =
                                          case_when(
                                            session==1~ 0.0348,
                                            session==5~ 0.0244,
                                            session==9 ~ 0.0977,
                                            session==13~ 0.0661,
                                            session==15  ~ 0.0736,
                                            session==16~ 0.1023,
                                            session==18  ~ 0.0968,
                                            session==19~ 0.10734))

dipper.ddl2$GammaDoublePrime <-dipper.ddl2$GammaDoublePrime %>% mutate(fix =
                                                                         case_when(
                                                                           time==1~ 0.2025,
                                                                           time==5~ 0.2025,
                                                                           time==9  ~ 0.2025,
                                                                           time==13~ 0.2025,
                                                                           time==15  ~ 0.2025,
                                                                           time==16~ 0.2025))
dipper.ddl2$GammaPrime <-dipper.ddl2$GammaPrime %>% mutate(fix =case_when(
  time==5~ 0.2025,
  time==9~ 0.2025,
  time==13  ~ 0.2025,
  time==15~ 0.2025,
  time==16  ~ 0.2025,
  time==18~0.2025))

opptotal<-(mark(br2, dipper.ddl2))
Opp<-data.frame(opptotal$group.labels)
Oppv2<-data.frame(Group=Opp[rep(seq_len(nrow(Opp)), each = 8), ])
dates<-data.frame(session=c("Spring 2017", "Spring 2018", 
                            "Spring 2019", "Spring 2020",
                            "Fall 2020", "Spring 2021",
                            "Fall 2021", "Spring 2022"))
SeasonsOpp<-data.frame(Session=rep(dates$session, times=48))
EstimatesOpp<-bind_cols(opptotal$results$derived, Oppv2, SeasonsOpp)
EstimatesOpp$Species<-paste("Opossum")
Estimates<-bind_rows(EstimatesRacc, EstimatesOpp)
Estimates$grid<- str_sub(Estimates$Group, -2,-1)
Estimates$Year<- str_sub(Estimates$Session, -4,-1)
Estimates<-Estimates %>% 
  mutate(Sex = case_when( 
    grepl("\\Male",Group ) ~ "Male",
    grepl("\\Female", Group) ~ "Female"
  )) 

Estimates<-Estimates %>% 
  mutate(Habitat = case_when( 
    grepl("\\GridB",Group ) ~ "Bottomland",
    grepl("\\GridP", Group) ~ "Pine",
    grepl("\\GridR",Group ) ~ "Riparian",
    grepl("\\GridW", Group) ~ "Wetland"
  )) 

Estimates<-Estimates %>% 
  mutate(Season = case_when( 
    grepl("\\Spring",Session) ~ "Spring",
    grepl("\\Fall", Session) ~ "Fall"
  )) 

#STEP 3- Divide the abundances by the effective trapping areas (based on home ranges) to calculate densities by each desired 
#grouping (by habitat, by grid, by sex grid, etc.)
SEED=666
set.seed(SEED)
TOL=1e-4
mVec=racctotal$results$derived[["N Population Size"]]$estimate
vMat=racctotal$results$derived.vcv[["N Population Size"]]  
nBoot=1000

bootDF=data.frame(t(MASS::mvrnorm(nBoot,mVec,vMat,tol=TOL)),
                  stringsAsFactors = T)

names(bootDF)=fn.genLab(1:ncol(bootDF),"boot")

Raccest<-subset(Estimates, Species=="Raccoon")
TotalRacc<-cbind(Raccest, bootDF)
eff<-read.csv("effarealist5Aug.csv")
areaRacc<- subset(eff, Species=="Raccoon")

for( R in 1:nrow(TotalRacc)){ ## R=25
  TotalRacc$effArea[R] = 
    areaRacc$effArea.km2[
      as.character(areaRacc$grid) == 
        as.character(TotalRacc$grid[R]) &
        as.character(areaRacc$Sex) == 
        as.character(TotalRacc$Sex[R]) &
        as.character(areaRacc$Season) == 
        as.character(TotalRacc$Season[R])]	
}

bootCols=which( grepl("boot",names(TotalRacc)))
varCols=which(! grepl("boot",names(TotalRacc)))
xInvEffArea= matrix(rep(1/TotalRacc$effArea,length(bootCols)),
                    ncol=length(bootCols))
TotalRacc[bootCols]=TotalRacc[bootCols]*xInvEffArea
rm(xInvEffArea)

varCols=which(! grepl("boot",names(TotalRacc)))
COLS= which(substr(names(TotalRacc),1,4)=="boot")
DensGridSexRacc = data.frame(
  grid=TotalRacc$grid,
  Habitat=TotalRacc$Habitat,
  Season=TotalRacc$Season,
  year=TotalRacc$Year,
  sex=TotalRacc$Sex,
  mean=apply(TotalRacc[,COLS],1,mean),
  se=apply(TotalRacc[,COLS],1,sd),
  lci=apply(TotalRacc[,COLS],1,quantile,0.025),
  uci=apply(TotalRacc[,COLS],1,quantile,0.975) )
DensGridSexRacc$Species<-paste("Raccoon")
varCols=which(! grepl("boot",names(TotalRacc)))
tt5<-aggregate(TotalRacc[,-varCols], by=list(Grid=TotalRacc$grid,Habitat=TotalRacc$Habitat, 
                                             Season=TotalRacc$Season, Year=TotalRacc$Year), sum)
COLS= which(substr(names(tt5),1,4)=="boot")
DensGridRacc2 = data.frame(
  Grid=tt5$Grid,
  Habitat=tt5$Habitat,
  Season=tt5$Season,
  Year=tt5$Year,
  mean=apply(tt5[,COLS],1,mean),
  se=apply(tt5[,COLS],1,sd),
  lci=apply(tt5[,COLS],1,quantile,0.025),
  uci=apply(tt5[,COLS],1,quantile,0.975) )
DensGridRacc2$Species<-paste("Raccoon")

varCols=which(! grepl("boot",names(tt5)))
tt6<-aggregate(tt5[,-varCols], by=list(Habitat=tt5$Habitat, 
                                             Season=tt5$Season), mean)
COLS= which(substr(names(tt6),1,4)=="boot")
DensHabSeason = data.frame(
  Habitat=tt6$Habitat,
  Season=tt6$Season,
  mean=apply(tt6[,COLS],1,mean),
  se=apply(tt6[,COLS],1,sd),
  lci=apply(tt6[,COLS],1,quantile,0.025),
  uci=apply(tt6[,COLS],1,quantile,0.975) )
DensHabSeason$Species<-paste("Raccoon")



##Opossum density
fn.genLab=function(labs,str="lab"){ ## labs= c(1:100); str="lab"
  maxChar=ifelse(max(nchar(labs))==1,2,max(nchar(labs)))
  paste0(str,formatC(labs,width=maxChar,flag = "0"))
}

TOL=1e-4
mVec=opptotal$results$derived[["N Population Size"]]$estimate
vMat=opptotal$results$derived.vcv[["N Population Size"]]  
nBoot=1000

bootDF=data.frame(t(MASS::mvrnorm(nBoot,mVec,vMat,tol=TOL)),
                  stringsAsFactors = T)
Oppest<-subset(Estimates, Species=="Opossum")
names(bootDF)=fn.genLab(1:ncol(bootDF),"boot")
TotalOpp<-cbind(Oppest, bootDF)
areaOpp<- subset(eff, Species=="Opossum")

for( R in 1:nrow(TotalOpp)){ ## R=25
  TotalOpp$effArea[R] = 
    areaOpp$effArea.km2[
      as.character(areaOpp$grid) == 
        as.character(TotalOpp$grid[R]) &
        as.character(areaOpp$Sex) == 
        as.character(TotalOpp$Sex[R]) &
        as.character(areaOpp$Season) == 
        as.character(TotalOpp$Season[R])]	
}

bootCols=which( grepl("boot",names(TotalOpp)))
varCols=which(! grepl("boot",names(TotalOpp)))
xInvEffArea= matrix(rep(1/TotalOpp$effArea,length(bootCols)),
                    ncol=length(bootCols))
TotalOpp[bootCols]=TotalOpp[bootCols]*xInvEffArea
rm(xInvEffArea)

varCols=which(! grepl("boot",names(TotalOpp)))
COLS= which(substr(names(TotalOpp),1,4)=="boot")
DensGridSexOpp = data.frame(
  grid=TotalOpp$grid,
  Habitat=TotalOpp$Habitat,
  Season=TotalOpp$Season,
  year=TotalOpp$Year,
  sex=TotalOpp$Sex,
  mean=apply(TotalOpp[,COLS],1,mean),
  se=apply(TotalOpp[,COLS],1,sd),
  lci=apply(TotalOpp[,COLS],1,quantile,0.025),
  uci=apply(TotalOpp[,COLS],1,quantile,0.975) )
DensGridSexOpp$Species<-paste("Opossum")
varCols=which(! grepl("boot",names(TotalOpp)))
tt6<-aggregate(TotalOpp[,-varCols], by=list(Grid=TotalOpp$grid,Habitat=TotalOpp$Habitat, 
                                            Season=TotalOpp$Season, Year=TotalOpp$Year), sum)
COLS= which(substr(names(tt6),1,4)=="boot")
DensGridOpp = data.frame(
  Grid=tt6$Grid,
  Habitat=tt6$Habitat,
  Season=tt6$Season,
  Year=tt6$Year,
  mean=apply(tt6[,COLS],1,mean),
  se=apply(tt6[,COLS],1,sd),
  lci=apply(tt6[,COLS],1,quantile,0.025),
  uci=apply(tt6[,COLS],1,quantile,0.975) )    
DensGridOpp$Species<-paste("Opossum")


STAT<-bind_rows(DensGridRacc2, DensGridOpp)
full2<-bind_rows(DensGridSexRacc, DensGridSexOpp)
full2<-full2[!(full2$grid=="P1" & full2$year>=2020|(full2$grid=="P7" & full2$year<2020)),] #remove the grids in seasons that weren't trapped
STAT<-STAT[!(STAT$Grid=="P1" & STAT$Year>=2020|(STAT$Grid=="P7" & STAT$Year<2020)),]

write.csv(STAT, "GridDensities.csv")
write.csv(full2, "GridSexDensities.csv")

GridDensity<-read.csv("GridDensities.csv")
GridDensitySex<-read.csv("GridSexDensities.csv")

GridDensitySex<-GridDensitySex[!(GridDensitySex$grid=="P1" & GridDensitySex$year>=2020|(GridDensitySex$grid=="P7" & GridDensitySex$year<2020)),]
GridDensity<-GridDensity[!(GridDensity$Grid=="P1" & GridDensity$Year>=2020|(GridDensity$Grid=="P7" & GridDensity$Year<2020)),]
names(GridDensity)[names(GridDensity) == "mean"] <- "Density"


STAT<-read.csv("GridDensities.csv")
full2<-read.csv("GridSexDensities.csv")
full2<-full2[!(full2$grid=="P1" & full2$year>=2020|(full2$grid=="P7" & full2$year<2020)),]
STAT<-STAT[!(STAT$Grid=="P1" & STAT$Year>=2020|(STAT$Grid=="P7" & STAT$Year<2020)),]

#STEP 4- compare the densities derived above
#Comparing raccoon vs opossum
racc<- subset(STAT, Species=="Raccoon")
opp<- subset(STAT, Species=="Opossum")

bothmod<-(lmer(mean~Season*Habitat*Species+(1|Grid:Habitat)+(1|Year), 
               na.action=na.pass, data=STAT, REML=F))
dredge(bothmod)
output(bothmod)
bothfit<-(lmer(mean~Species+Habitat+(1|Grid:Habitat)+(1|Year), 
               na.action=na.pass, data=STAT, REML=T))

emmeans(bothfit, pairwise~Habitat, type="response")
emmip(bothfit, ~Species|Habitat, CIs=TRUE, CIarg=aes(lwd=1, color="black", linetype=1), type="response")

#Raccoon habitat and season

raccfull<-(lmer(mean~Season*Habitat+(1|Grid:Habitat)+(1|Year), 
                na.action=na.pass, data=racc, REML=F))
dredge(raccfull)
raccfit<-(lmer(mean~Habitat+Season+(1|Year), 
               na.action=na.pass, data=racc, REML=T))
r.squaredGLMM(raccfit)
emmeans(raccfit, pairwise~Habitat|Season, type="response")
emmip(raccfit, ~Season, CIs=TRUE, CIarg=aes(lwd=1, color="black", linetype=1), type="response")
ggpredict(raccfit, ci.lvl=0.95, terms=c("Season", "Habitat"))
summary(raccfit)

#opossum habitat and season
oppfull<-(lmer(mean~Season*Habitat+(1|Grid:Habitat)+(1|Year), 
               na.action=na.pass, data=opp, REML=F))
dredge(oppfull)
oppfit<-(lmer(mean~Season*Habitat+(1|Grid:Habitat)+(1|Year), 
              na.action=na.pass, data=opp, REML=T))
emmip(oppfit, ~Season|Habitat, CIs=TRUE, CIarg=aes(lwd=1, color="black", linetype=1), type="response")
emmeans(oppfit, pairwise~Season|Habitat, type="response")
emmeans(oppfit, pairwise~Habitat|Season, type="response")
ggpredict(oppfit, ci.lvl=0.95, terms=c("Habitat"), condition = c(Season="Fall"))
ggpredict(oppfit, ci.lvl=0.95, terms=c("Habitat"), condition = c(Season="Spring"))

ggpredict(oppfit, ci.lvl=0.95, terms=c("Season"), condition = c(Habitat="Wetland"))

#Raccoon sex ratio
raccdens<- subset(full2, Species=="Raccoon")
oppdens<-subset(full2, Species=="Opossum")
raccratio<-(lmer(mean~Season*Habitat*sex+(1|grid:Habitat)+(1|year), 
                 na.action=na.pass, data=raccdens, REML=F))
dredge(raccratio)
output(raccratio)
raccratiofit<-(lmer(mean~Season*Habitat*sex+(1|grid:Habitat)+(1|year), 
                    na.action=na.pass, data=raccdens, REML=T))

emmeans(raccratiofit, pairwise~sex|Season|Habitat, type="response")
emmip(raccratiofit, ~sex|Season|Habitat, CIs=TRUE, CIarg=aes(lwd=1, color="black", linetype=1), type="response")
ggpredict(raccratiofit, ci.lvl=0.95, terms=c("sex"), condition = c(Habitat="Bottomland", Season="Spring"))
ggpredict(raccratiofit, ci.lvl=0.95, terms=c("sex"), condition = c(Habitat="Bottomland", Season="Fall"))
ggpredict(raccratiofit, ci.lvl=0.95, terms=c("sex"), condition = c(Habitat="Riparian", Season="Fall"))

#Opossum sex ratio
oppratio<-(lmer(mean~Season*Habitat*sex+(1|grid:Habitat)+(1|year), 
                na.action=na.pass, data=oppdens, REML=F))
dredge(oppratio)
output(oppratio)
oppratiofit<-(lmer(mean~Season*Habitat+sex+(1|grid:Habitat)+(1|year), 
                na.action=na.pass, data=oppdens, REML=T))
summary(oppratiofit)
ggpredict(oppratiofit, ci.lvl=0.95, terms=c("sex"), condition = c(Habitat="Wetland"))

##AgeRatio
agedata<-read.csv("ageratio.csv")
raccage<- subset(agedata, species=="Raccoon")
oppage<- subset(agedata, species=="Opossum")
raccmodel<-(glmer(AGE~siteType*Season+(1|year), 
              family=binomial(link="logit"), na.action=na.pass, data=raccage))

dredge(raccmodel)
fitmodel<-(glmer(AGE~Season+(1|year), 
                 family=binomial(link="logit"), na.action=na.pass, data=raccage))
emmeans(fitmodel, pairwise~Season, type="response")
summary(fitmodel)
ggpredict(fitmodel, terms=c("Season"))

r.squaredGLMM(fitmodel)


modelopp<-(glmer(AGE~siteType*Season+(1|year), 
                 family=binomial(link="logit"), na.action=na.pass, data=oppage))

dredge(modelopp)
fitmodelopp<-(glmer(AGE~Season+(1|year), 
                    family=binomial(link="logit"), na.action=na.pass, data=oppage))
ggpredict(fitmodelopp, terms=c("Season"))

