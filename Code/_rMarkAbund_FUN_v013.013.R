## _rMarkAbund_FUN_v013.01.R


## changed interpolation of age and sex
## changed p-value calculation for z-test in boot function
## see tag "## Changed 2021.10.02 ##"

## changed the getData function to read summary hoime range data


#####################################################################

fn.getCaptureData=function(FN.Date="__Data/sessionMinMaxDates.csv",
	FN.Capture="__Data/3yr_abundance_capture_data_Lonnie.csv",
      Interpolate=T){

#FN.Date="__Data/sessionMinMaxDates2.csv"
#FN.Capture="__Data/capture_data_adjust.csv"

  ################################################################
  #read in capture database for all three years and all species
  
  dateData=read.csv(FN.Date)
  cols<-c("session","grid")
  dateData[,cols] <- lapply(dateData[,cols],as.factor)
  dateData$startDate = as.Date(as.character(dateData$startDate), "%m/%d/%Y")
  dateData$endDate   = as.Date(as.character(dateData$endDate), "%m/%d/%Y")
  
  
  allData = read.csv(FN.Capture)
  COLS= c("date", "Grid","Session","site.type", "Trap",   
  	"Species","animal_id","Recapture", "Age", "Sex","Testes.Descended",
  	"Swollen.Vulva" , "Teats.Visible", "Lactating", "Offspring" ,             
  	"Weight..kg.")        
  allData=allData[,COLS]; rm(COLS)
  names(allData)= sapply(names(allData), function(x) { 
  	## x = "Grid.Session"
  	y=gsub(".","_",gsub("..","_",x,fixed=T),fixed=T)
        y=strsplit(y,"_")[[1]]
        y[1]= tolower(y[1])
        if(length(y) > 1)
          y[-1] = sapply( y[-1], function(a) ## a=y[2]
  		paste0(toupper(substr(a,1,1)),
  		tolower(substr(a,2,nchar(a)))))
        y=paste(y,collapse="") })
  
  allData$date = as.Date(as.character(allData$date), "%m/%d/%Y")
  allData$year=factor(format(allData$date,"%Y"))
  
  ##-------------------------------------------------------------
  ##-------------------------------------------------------------
  ##-------------------------------------------------------------
  ##-------------------------------------------------------------
  ## manual corrections of some dates...
  ## NEED TO CONFIRM WITH RAW DATA SHEETS
  ##-------------------------------------------------------------
  ##-------------------------------------------------------------
  ##-------------------------------------------------------------
  ##-------------------------------------------------------------

  for(COL in 1:ncol(allData)){ ## COL=2
    if( class(allData[,COL]) == "character")
      allData[,COL] = factor(allData[,COL])
  }
  
  levels(allData$grid)[levels(allData$grid)=="B2-1"]="B2"
  
  ROW = which(allData$year == "2017" & allData$grid == "B4" &
  	allData$date == "2017-01-15")
  ## sum(allData$animalId == allData$animalId[ROW])
  allData$date[ROW] = "2017-01-14"
  
  ROW = which(allData$year == "2017" & allData$grid == "B4" &
  	allData$date == "2017-01-17")
  ## sum(allData$animalId %in% allData$animalId[ROW] & allData$year %in% allData$year[ROW])
  ## allData$grid[ allData$animalId %in% allData$animalId[ROW]]
  allData$date[ROW] = "2017-01-14"
  
  ROW = which(allData$year == "2017" & allData$grid == "P2" &
  	allData$date == "2017-01-18")
  ## sum(allData$animalId %in% allData$animalId[ROW] & allData$year %in% allData$year[ROW])
  ## allData$grid[ allData$animalId %in% allData$animalId[ROW]]
  allData$date[ROW] = "2017-01-17"
  
  ROW = which(allData$year == "2018" & allData$grid == "B6" &
  	allData$date == "2018-02-13")
  ## sum(allData$animalId %in% allData$animalId[ROW] & allData$year %in% allData$year[ROW])
  ## allData$grid[ allData$animalId %in% allData$animalId[ROW]]
  allData$date[ROW] = "2018-02-03"
  
  ROW = which(allData$year == "2018" & allData$grid == "B1" &
  	allData$date == "2018-02-17")
  ## sum(allData$animalId %in% allData$animalId[ROW] & allData$year %in% allData$year[ROW])
  ## allData$grid[ allData$animalId %in% allData$animalId[ROW]]
  allData$date[ROW] = "2018-02-27"
  
  allData$trapDay =allData$endDate=allData$startDate=allData$session=NA
  
  for(R in 1:nrow(allData)){ ## R = 671
    Y=as.numeric(as.character(allData$year[R]))
    G=as.character(allData$grid[R])
    D= allData$date[R]; ## Y; G; D
  
    ROW= which( dateData$year == Y & dateData$grid == G &
  	dateData$startDate <= D & dateData$endDate >= D)
    if(length(ROW) > 1) stop(" length(ROW) > 1" )
    if(length(ROW) == 0) stop(" length(ROW) == 0" )
  
    allData$session[R] = as.character(dateData$session[ROW])
    allData$startDate[R]	= as.character(dateData$startDate[ROW])
    allData$endDate[R]	= as.character(dateData$endDate[ROW])
  }                  
  
    allData$session=factor(  allData$session)
    allData$startDate=as.Date(allData$startDate)
    allData$endDate  =as.Date(allData$endDate)
    allData$trapDay    = as.numeric(allData$date-
  					allData$startDate)+1
  
  grid_nDays=list()
  for(G in levels(allData$grid)){ ## G="B1"
    grid_nDays[[G]]=sapply(levels(allData$year) ,function(y) ## y="2017"
  	as.numeric(diff(range(allData$date[allData$grid==G & allData$year==y]))))
  }
  mean(unlist(grid_nDays))

    #############################################################################
    ## Some age strings converted to dates by Excel; Correct these
    #############################################################################
   
    LEVS=levels(allData$weightKg)
    LEVS[LEVS=="" |LEVS=="NR" | LEVS == "N/A"]=NA
    
    oldNA=sum(is.na(LEVS))
    LEVS=as.numeric(LEVS)
    newNA=sum(is.na(LEVS))
    levels(allData$weightKg)=LEVS
    allData$weightKg=as.numeric(as.character(allData$weightKg)); rm(LEVS)
    
    LEVS=levels(allData$age)
    LEVS[LEVS=="" | LEVS == "N/A"]=NA
    LEVS[grepl("Juv",LEVS,ignore.case=T) | grepl("Year",LEVS,ignore.case=T)]="1.1"
    LEVS[grepl("Adult",LEVS,ignore.case=T) ]="5.1"
    LEVS
    
    for(N in 1:length(LEVS)){ ## N=4
      X=LEVS[N]
      numX=as.numeric(X); numX; X
      if( !is.na(numX)) {
          if(numX > 40000){
        	  X = gsub("0","",format(as.Date(numX, 
			origin = "1899-12-30"),"%m-%d"))
          } 
       }
    
      numX=as.numeric(X); numX
      if( is.na(numX)) {
         X=gsub("+","",gsub("<","",gsub(">","",gsub("~","",X,
    		fixed=T),fixed=T),fixed=T),fixed=T)
         X=mean(as.numeric(strsplit(X,"-",fixed=T)[[1]]))
      }
       
      LEVS[N]=X
    }
    
    oldNA=sum(is.na(LEVS))
    LEVS=as.numeric(LEVS)
    newNA=sum(is.na(LEVS))
    if(oldNA != newNA) fn.FatalError("oldNA != newNA")
    
    newLEVS=rep(NA,length(LEVS))
    newLEVS[LEVS < 2]="Young"; newLEVS[LEVS >= 2] = "Adult"
    levels(allData$age)=newLEVS; rm(LEVS,newLEVS)
    allData$age=factor(as.character(allData$age))

    #############################################################################
    ## getting obvious missing sex information 
    #############################################################################
   
    allData=allData[order(allData$date),]; rownames(allData)=NULL

    STR=levels(allData$sex); STR
    STR[STR=="y"]="Yes";STR[STR=="Y"]="Yes";STR[STR=="N"]="No";
    STR[STR==""]=NA;STR[STR=="N/A"]=NA
    levels(allData$sex)=STR

    STR=levels(allData$testesDescended); STR
    STR[STR=="y"]="Yes";STR[STR=="Y"]="Yes";STR[STR=="N"]="No";
    STR[STR==""]=NA;STR[STR=="N/A"]=NA
    levels(allData$testesDescended)=STR

    STR=levels(allData$swollenVulva); STR
    STR[STR=="y"]="Yes";STR[STR=="Y"]="Yes";STR[STR=="N"]="No";
    STR[STR==""]=NA;STR[STR=="N/A"]=NA
    levels(allData$swollenVulva)=STR

    STR=levels(allData$teatsVisible); STR
    STR[STR=="y"]="Yes";STR[STR=="Y"]="Yes";STR[STR=="N"]="No";
    STR[STR==""]=NA;STR[STR=="N/A"]=NA
    levels(allData$teatsVisible)=STR

    STR=levels(allData$lactating); STR
    STR[STR=="y"]="Yes";STR[STR=="Y"]="Yes";STR[STR=="N"]="No";
    STR[STR==""]=NA;STR[STR=="N/A"]=NA;STR[STR=="NR"]=NA
    levels(allData$lactating)=STR

    for(ID in unique(allData$animalId) ){ ## ID=215
      ROWS= which(allData$animalId == ID)
      Sex = unique(as.character(allData$sex[ROWS]))
      Sex=Sex[!is.na(Sex)]
      allData[ROWS,]

      if(length(Sex) > 1) ## if multiple sexes recorded in error
		## treat as missing
	      Sex=NA
      if(length(Sex) == 1) 
		allData$sex[allData$animalId == ID] = Sex
      if(length(Sex) == 0) {
        
        if(	any(!is.na(as.character(
		allData$testesDescended[ROWS])))) Sex=c(Sex,"M")
        if(	any(
		!is.na(allData$swollenVulva[ROWS]) &
		!is.na(allData$teatsVisible[ROWS]) &
		!is.na(allData$lactating[ROWS])) ) Sex=c(Sex,"F")
	  if(length(Sex) > 1) stop("  check sex")
	  if(length(Sex) == 1)
		allData$sex[allData$animalId == ID] = Sex
      }
    }

  if(Interpolate){

    #############################################################################
    ## interpolate age 
    #############################################################################
    
    allData$age02=NA
    for(ID in unique(allData$animalId)){ ## ID = 14
      ages=as.character(allData$age[allData$animalId==ID])
      dates=allData$date[allData$animalId==ID]
      yrs=as.numeric(format(dates,"%Y"))
    
      if(length(ages)>1 & sum(is.na(ages)) < length(ages)  ){
        if( length(unique(yrs))==1     ){
    	ages[]=unique(ages)[1]
        } else {
    	if(is.na(ages[1])){
            ages[1]= na.omit(ages)[1]
          }
    
          if(ages[1]=="Adult") ages[]="Adult"
    	if(ages[1]=="Young"){
            ages[ yrs==yrs[1] ] = "Young"
            ages[ yrs > yrs[1] ] = "Adult"
          }
        }
      }
      allData$age02[allData$animalId==ID]=ages
    }
    
    allData$age02=factor(allData$age02)
    sum(is.na(allData$age));sum(is.na(allData$age02))
    
    mData=na.omit(allData[,c("age02","weightKg","species")]); rownames(mData)=NULL
    mod1=glm(age02~species*weightKg,data=mData,family=binomial)
    mod2=glm(age02~species*(weightKg+I(weightKg^2)),data=mData,family=binomial)
    mod=mod1
    if(AIC(mod2) < AIC(mod1)) mod=mod2
    
    roc_obj <- pROC::roc(mData$age02,predict(mod,type="response"),quiet=T)
    Cut= as.numeric(pROC::coords(roc_obj, "best", "threshold")[1])
    
    allData$age03=allData$age02
    ROWS= which(is.na(allData$age03))
    pred= predict(mod,data.frame(weightKg=allData$weightKg[ROWS],
    				    species=allData$species[ROWS]),type="response")
    
    lrROWS=ROWS[ !is.na(pred) & pred < Cut]
    grROWS=ROWS[ !is.na(pred) & pred >= Cut]
    allData$age03[lrROWS]=levels(mData$age02)[1]
    allData$age03[grROWS]=levels(mData$age02)[2]
    
    allData$age03
    
    allData$age04=NA
    for(ID in unique(allData$animalId)){ ## ID = 14
      ages=as.character(allData$age03[allData$animalId==ID])
      dates=allData$date[allData$animalId==ID]
      yrs=as.numeric(format(dates,"%Y"))
    
      if(length(ages)>1 & sum(is.na(ages)) < length(ages)  ){
        if( length(unique(yrs))==1     ){
    	ages[]=unique(ages)[1]
        } else {
    	if(is.na(ages[1])){
            ages[1]= na.omit(ages)[1]
          }
    
          if(ages[1]=="Adult") ages[]="Adult"
    	if(ages[1]=="Young"){
            ages[ yrs==yrs[1] ] = "Young"
            ages[ yrs > yrs[1] ] = "Adult"
          }
        }
      }
      allData$age04[allData$animalId==ID]=ages
    }
    
    allData$age04=factor(allData$age04)
    sum(is.na(allData$age03));sum(is.na(allData$age04))
    
    mData=na.omit(allData[,c("age04","grid","species")]); rownames(mData)=NULL
    mod=glm(age04~grid*species,data=mData,family=binomial)
    
    roc_obj <- pROC::roc(mData$age04,predict(mod,type="response"),quiet=T)
    Cut= as.numeric(pROC::coords(roc_obj, "best", "threshold")[1])
    
    allData$age05=allData$age04
    ROWS= which(is.na(allData$age05))
    pred= predict(mod,data.frame(grid=allData$grid[ROWS],
    	species = allData$species[ROWS]),type="response")
    
    lrROWS=ROWS[ !is.na(pred) & pred < Cut]
    grROWS=ROWS[ !is.na(pred) & pred >= Cut]
    allData$age05[lrROWS]=levels(mData$age04)[1]
    allData$age05[grROWS]=levels(mData$age04)[2]
    
    sum(is.na(allData$age04));sum(is.na(allData$age05))
    
    allData$age00=allData$age; 
    allData$age=allData$age05
    allData$age=relevel(allData$age,"Young")
    allData=allData[,-which(grepl("age0",names(allData)))]; names(allData)
    
    #########################################################################
    ## interpolate sex based on age and weight
    #########################################################################
    
    LEVS=levels(allData$sex)
    LEVS[LEVS=="" | LEVS == "N/A"]=NA
    LEVS
    levels(allData$sex)=LEVS; rm(LEVS)
    allData$sex=factor(as.character(allData$sex))
    
    allData$sex02=NA
    for(ID in unique(allData$animalId)){ ## ID = 557
      sexes=as.character(allData$sex[allData$animalId==ID])
      if(length(sexes)>1 & sum(is.na(sexes)) < length(sexes)  ){
        sexes[]=unique(sexes)[1]
      }
      allData$sex02[allData$animalId==ID]=sexes
    }
    
    allData$sex02=factor(allData$sex02)
    sum(is.na(allData$sex));sum(is.na(allData$sex02))
    
    mData=na.omit(allData[,c("sex02","weightKg","age","species")]); rownames(mData)=NULL
    mods=list(
      mod1=glm(sex02~species*weightKg,data=mData,family=binomial),
      mod2=glm(sex02~species*(weightKg+I(weightKg^2)),data=mData,family=binomial),
      mod3=glm(sex02~species*age,data=mData,family=binomial),
      mod4=glm(sex02~species*age*weightKg,data=mData,family=binomial),
      mod5=glm(sex02~species*age*(weightKg+I(weightKg^2)),data=mData,family=binomial,
    	control=glm.control(maxit=1000))
    )
    AICvals=sapply(mods,AIC)
    mod=mods[[which(AICvals==min(AICvals))]]
    
    roc_obj <- pROC::roc(mData$sex02,predict(mod,type="response"),quiet=T)
    Cut= as.numeric(pROC::coords(roc_obj, "best", "threshold")[1])
    
    allData$sex03=allData$sex02
    ROWS= which(is.na(allData$sex03))
    pred= predict(mod,data.frame(weightKg=allData$weightKg[ROWS],
    	age=allData$age[ROWS],
    	species=allData$species[ROWS]),type="response")
    
    lrROWS=ROWS[ !is.na(pred) & pred < Cut]
    grROWS=ROWS[ !is.na(pred) & pred >= Cut]
    allData$sex03[lrROWS]=levels(mData$sex02)[1]
    allData$sex03[grROWS]=levels(mData$sex02)[2]
    
    allData$sex03
    
    allData$sex04=NA
    for(ID in unique(allData$animalId)){ ## ID = 14
      sexes=as.character(allData$sex03[allData$animalId==ID])
      if(length(sexes)>1 & sum(is.na(sexes)) < length(sexes)  ){
        sexes[]=unique(sexes)[1]
      }
      allData$sex04[allData$animalId==ID]=sexes
    }
    
    allData$sex04=factor(allData$sex04)
    sum(is.na(allData$sex03));sum(is.na(allData$sex04))
    
    mData=na.omit(allData[,c("sex04","grid","species")]); rownames(mData)=NULL
    mod=glm(sex04~grid*species,data=mData,family=binomial)
    
    roc_obj <- pROC::roc(mData$sex04,predict(mod,type="response"),quiet=T)
    Cut= as.numeric(pROC::coords(roc_obj, "best", "threshold")[1])
    
    allData$sex05=allData$sex04
    ROWS= which(is.na(allData$sex05))
    pred= predict(mod,data.frame(grid=allData$grid[ROWS],
    	species = allData$species[ROWS]),type="response")
    
    lrROWS=ROWS[ !is.na(pred) & pred < Cut]
    grROWS=ROWS[ !is.na(pred) & pred >= Cut]
    allData$sex05[lrROWS]=levels(mData$sex04)[1]
    allData$sex05[grROWS]=levels(mData$sex04)[2]
    
    sum(is.na(allData$sex04));sum(is.na(allData$sex05))
    
    allData$sex00=allData$sex; 
    allData$sex=allData$sex05
  
    allData=allData[,-which(grepl("sex0",names(allData)))]; names(allData)
  
  } ## end if Interpolate

  sum(is.na(allData$sex))
  allData=drop.levels(allData[!is.na(allData$sex),],reorder=F)
 
  
  allData$ageStr=allData$age
  allData$sexStr=allData$sex
  allData$age=as.numeric(allData$age)-1 ## 0 ="Young"; 1 = "Adult"
  allData$sex=as.numeric(allData$sex)-1 ## 0 ="Female"; 1 = "Male"
  allData
 
}

#####################################################################
if(1==2){
  rawData=sppData;covNames=c("grid","year","sex","age","siteType"); 
  time.units=1; byYear=F
}

time.units=1
fn.convertCapGen<-function(rawData,covNames=NULL,byYear=T,time.units=1)  
{

  names(rawData)
  maxTime=100; if(maxTime < time.units) maxTime=100*time.units

  COLS=c("date","year","trapDay","grid","session","animalId")
  COLS=c(COLS,covNames); 
  COLS=COLS[!duplicated(COLS)]

  rawData=rawData[,COLS]
  rawData$round= as.numeric(substr(as.character(rawData$session),2,2))
  rawData$trapDay=rawData$trapDay+ max(rawData$trapDay)*(rawData$round-1)
  LEVS= seq(min(rawData$trapDay),max(rawData$trapDay))

  if(! byYear){
    rawData$trapDay=rawData$trapDay + (maxTime)*(
		as.numeric(rawData$year)-1)
    LEVS= unlist(lapply( (as.numeric(unique(rawData$year))-1),
	function(x) LEVS+(maxTime*x)) )
  }

  rawData$trapDayStr=factor(rawData$trapDay, levels=LEVS)


  rawData$idVar=paste0(rawData$animalId,"_",rawData$grid)
  if(byYear) rawData$idVar=paste0(rawData$idVar,"_",rawData$year)
  rawData$idVar=factor(rawData$idVar)
  
  chDF=tiVec=nDays=NULL
  for(ID in levels(rawData$idVar)){ ## ID=levels(rawData$idVar)[1]
    xData=rawData[rawData$idVar==ID,]
    rownames(xData)=NULL
    xData$idVar= as.character(xData$idVar)
    CH=table(xData$idVar,xData$trapDayStr)
    TI= diff( as.numeric(colnames(CH)))
    TI[TI <= time.units] = 0
    TI[TI>0]=1
    if(nrow(CH) != 1) stop("check CH")
    CH= paste(as.character(CH),collapse="")
    if(is.null(nDays)) nDays=nchar(CH)
    if(nchar(CH) != nDays) stop("check nDays")

    if(is.null(tiVec)) tiVec=TI
    if(! all(tiVec == TI) ) stop("tiVec")

    xData$ch=CH; 
    chDF=rbind(chDF, xData[1,])
  }
  rownames(chDF)=NULL
  if(!byYear) chDF$year="all.years"
  list(chDF=chDF,tiVec=tiVec,nDays=nDays)

}

############################################################################
############################################################################


if(1==2){
  mod=bestMod;mParamList=paramList;mData=mList$chDF;
  mTime=mList$tiVec;mType=MOD.TYPE;mGrps="grid";
  mName="bestMod"; mEffArea=effArea;nBoot=10

}

fn.refitMarkMods=function(mod,mParamList=paramList,mData=mList$chDF,
	mTime=mList$tiVec,mType=MOD.TYPE,mGrps="grid",mName="bestMod",
	mEffArea=NULL,nBoot=100){

  ## REMEMBER THIS MODEL ASSUMES:
	## 1. original mark model run with grids as groups
      ## 2. 6 grids/habitat type to calculate effective areas
  
  RES.Abund= RES.Dens = list()
  names(mod)
  if(!all(substring(mod$group.labels,1,4)=="grid"))
	fn.FatalError("original model assumes grps are girds")

  ##################################################################
  ##################################################################
  ## grid

  ##----------------------------------------------------------------
  ## boot

  TOL=1e-4
  mVec=mod$results$derived[["N Population Size"]]$estimate
  vMat=mod$results$derived.vcv[["N Population Size"]]  
  if(!lqmm::is.positive.definite(vMat))
    vMat=lqmm::make.positive.definite(vMat, tol=TOL)

  bootDF=data.frame(t(MASS::mvrnorm(nBoot,mVec,vMat,tol=TOL)))
  names(bootDF)=fn.genLab(1:ncol(bootDF),"boot")

  varDF= mod$design.data$f0
  varDF=varDF[,sapply(varDF,function(x) length(unique(x))>1)]
  if(any(substr(names(varDF),1,4)=="boot"))
	stop("col names can not start with 'boot' ")

  summDF = aggregate( bootDF, by=list(grid=varDF$grid),mean)
  COLS= names(summDF)[which(substr(names(summDF),1,4)=="boot")]

  RES.Abund[["gridDF"]] = data.frame(
		grid=summDF$grid,
		mean=apply(summDF[,COLS],1,mean),
		se=apply(summDF[,COLS],1,sd),
		lci=apply(summDF[,COLS],1,quantile,0.025),
		uci=apply(summDF[,COLS],1,quantile,0.975) )

  RES.Dens[["gridDF"]]=NA
  if(!is.null(mEffArea)){

    if(!all(as.character(summDF$grid) %in% 
	  as.character(mEffArea$grid)))
        fn.FatalError("check grids in effARea" )
 
    xEffArea= mEffArea[sapply(as.character(summDF$grid),function(x)
	which(as.character(mEffArea$grid)==x)),]$effArea
    xInvEffArea= 1/(matrix(rep(xEffArea,length(COLS)),ncol=length(COLS)))

    summDF[,COLS]=summDF[,COLS]*xInvEffArea  
    RES.Dens[["gridDF"]]= data.frame(
		grid=summDF$grid,
		mean=apply(summDF[,COLS],1,mean),
		se=apply(summDF[,COLS],1,sd),
		lci=apply(summDF[,COLS],1,quantile,0.025),
		uci=apply(summDF[,COLS],1,quantile,0.975) )
  }


  ##################################################################
  ##################################################################
  ## run new model

  PREFIX= paste0(mod$output,"_habSex")
  grps=c("siteType","sex")
  mod = try(mark(data = mData, model.name = NULL,
	model = mType,groups=grps, 
	time.intervals=mTime,model.parameters=mParamList,
	prefix = PREFIX,output=F,silent=T) )

  if( ( any(class(mod)=="try-error") | is.null(mod$results) |
	any(is.na(mod$results$derived.vcv[["N Population Size"]])) |
	all(mod$results$derived[["N Population Size"]]$se==0) |
	any(mod$results$derived[["N Population Size"]]$estimate> nMax) ))
	stop("modelError:: habSex")

  TOL=1e-4
  mVec=mod$results$derived[["N Population Size"]]$estimate
  vMat=mod$results$derived.vcv[["N Population Size"]]  
  if(!lqmm::is.positive.definite(vMat))
    vMat=lqmm::make.positive.definite(vMat, tol=TOL)

  bootDF=data.frame(t(MASS::mvrnorm(nBoot,mVec,vMat,tol=TOL)))
  names(bootDF)=fn.genLab(1:ncol(bootDF),"boot")

  varDF= mod$design.data$f0
  varDF=varDF[,sapply(varDF,function(x) length(unique(x))>1)]
  levels(varDF$sex)= c("Female","Male")
  if(any(substr(names(varDF),1,4)=="boot"))
	stop("col names can not start with 'boot' ")

  ##################################################################
  ##################################################################
  ## habitatType

  summDF = aggregate( bootDF, 
	by=list(siteType=varDF$siteType,
	        sex=varDF$sex),mean)
  summDF = aggregate( summDF[,-(1:2)], 
	by=list(siteType=summDF$siteType),sum)
  COLS= names(summDF)[which(substr(names(summDF),1,4)=="boot")]

  RES.Abund[["habDF"]] = data.frame(
		siteType=summDF$siteType,
		mean=apply(summDF[,COLS],1,mean),
		se=apply(summDF[,COLS],1,sd),
		lci=apply(summDF[,COLS],1,quantile,0.025),
		uci=apply(summDF[,COLS],1,quantile,0.975) )

  RES.Dens[["habDF"]]=NA
  if(!is.na(mEffArea)){
    
    habEffArea= aggregate(mEffArea$effArea,
		by=list(habType=mEffArea$habType),sum)
    names(habEffArea)[ncol(habEffArea)]="effArea"

    if(!all(toupper(substr(as.character(summDF$siteType),1,1)) %in% 
	  toupper(as.character(mEffArea$habType))))
        fn.FatalError("check habTypes in effARea" )

    levels(habEffArea$habType) = 
	levels(summDF$siteType)[sapply(levels(habEffArea$habType),
	function(x)which(toupper(substr(levels(
		summDF$siteType),1,1))==x))]

    xEffArea= habEffArea[sapply(as.character(summDF$siteType),function(x)
	which(as.character(habEffArea$habType)==x)),]$effArea
    xInvEffArea= 1/(matrix(rep(xEffArea,length(COLS)),ncol=length(COLS)))

    summDF[,COLS]=summDF[,COLS]*xInvEffArea ## 6 grids/habitat type
    RES.Dens[["habDF"]]= data.frame(
		siteType=summDF$siteType,
		mean=apply(summDF[,COLS],1,mean),
		se=apply(summDF[,COLS],1,sd),
		lci=apply(summDF[,COLS],1,quantile,0.025),
		uci=apply(summDF[,COLS],1,quantile,0.975) )
  }

  ##################################################################
  ##################################################################
  ## sex

  summDF = aggregate( bootDF, 
	by=list(siteType=varDF$siteType,
	        sex=varDF$sex),mean)
  summDF = aggregate( summDF[,-(1:2)], 
	by=list(sex=summDF$sex),mean)
  COLS= names(summDF)[which(substr(names(summDF),1,4)=="boot")]

  RES.Abund[["sexDF"]] = data.frame(
		sex=summDF$sex,
		mean=apply(summDF[,COLS],1,mean),
		se=apply(summDF[,COLS],1,sd),
		lci=apply(summDF[,COLS],1,quantile,0.025),
		uci=apply(summDF[,COLS],1,quantile,0.975) )

  RES.Dens[["sexDF"]]=NA
  if(!is.na(mEffArea)){

    sexEffArea= mean(habEffArea$effArea)
    summDF[,COLS]=summDF[,COLS]/ sexEffArea ## 6 grids/habitatType
    RES.Dens[["sexDF"]]= data.frame(
		sex=summDF$sex,
		mean=apply(summDF[,COLS],1,mean),
		se=apply(summDF[,COLS],1,sd),
		lci=apply(summDF[,COLS],1,quantile,0.025),
		uci=apply(summDF[,COLS],1,quantile,0.975) )
  }

  ##################################################################
  ##################################################################
  ## habitat & sex

  summDF = aggregate( bootDF, 
	by=list(siteType=varDF$siteType,
	        sex=varDF$sex),mean)
  COLS= names(summDF)[which(substr(names(summDF),1,4)=="boot")]

  RES.Abund[["habSexDF"]] = data.frame(
		siteType=summDF$siteType,
		sex	  =summDF$sex,
		mean=apply(summDF[,COLS],1,mean),
		se=apply(summDF[,COLS],1,sd),
		lci=apply(summDF[,COLS],1,quantile,0.025),
		uci=apply(summDF[,COLS],1,quantile,0.975) )

  RES.Dens[["habSexDF"]]=NA
  if(!is.na(mEffArea)){

   habSexEffArea=rbind(habEffArea,habEffArea)
    if(!all(as.character(summDF$grid) %in% 
	  as.character(mEffArea$grid)))
        fn.FatalError("check grids in effARea" )

    xEffArea= habEffArea[sapply(as.character(summDF$siteType),function(x)
	which(as.character(habEffArea$habType)==x)),]$effArea
    xInvEffArea= 1/(matrix(rep(xEffArea,length(COLS)),ncol=length(COLS)))

    summDF[,COLS]=summDF[,COLS]*xInvEffArea  
    RES.Dens[["habSexDF"]]= data.frame(
		siteType=summDF$siteType,
		sex	  =summDF$sex,
		mean=apply(summDF[,COLS],1,mean),
		se=apply(summDF[,COLS],1,sd),
		lci=apply(summDF[,COLS],1,quantile,0.025),
		uci=apply(summDF[,COLS],1,quantile,0.975) )
  }

  RES = list(RES.Abund=RES.Abund,RES.Dens=RES.Dens)
  RES

  ## RES.Dens$habSexDF$mean[1:4]+RES.Dens$habSexDF$mean[5:8]
  ## RES.Dens$habDF$mean

}

############################################################################
############################################################################


if(1==2){
  mod=bestMod;mParamList=paramList;mData=mList$chDF;
  mTime=mList$tiVec;mType=MOD.TYPE;mGrps="grid";
  mName="bestMod"; mEffArea=effArea;nBoot=100;.SEED=SEED
  pAdjMethod="bonferroni"

  MU=mod$results$derived[["N Population Size"]]$estimate
  SE=mod$results$derived[["N Population Size"]]$se
  VMAT=mod$results$derived.vcv[["N Population Size"]] 

##  x=t(MASS::mvrnorm(10000,MU,VMAT))
##  mu=apply(x,1,mean);se=apply(x,1,sd)
##  layout(1:2);   plot(MU,mu);plot(SE,se)

}

fn.refitMarkMods2=function(mod,mParamList=paramList,mData=mList$chDF,
	mTime=mList$tiVec,mType=MOD.TYPE,mGrps="grid",mName="bestMod",
	mEffArea=NULL,nBoot=100,pAdjMethod=NA,.SEED=SEED){

  ## REMEMBER THIS MODEL ASSUMES:
	## 1. original mark model run with grids as groups
      ## 2. 6 grids/habitat type to calculate effective areas
      ## 3. pooling of errors is done using parametric bootstrap:
		## http://ms.mcmaster.ca/~bolker/emdbook/book.pdf
		## pg 337
  
  RES.Abund= RES.Dens = list()

  ## names(mod)
  if(!all(substring(mod$group.labels,1,4)=="grid"))
	fn.FatalError("original model assumes grps are girds")

  ##################################################################
  ##################################################################
  ## run new model

  PREFIX= paste0(mod$output,"_gridSex")
  grps=c("grid","sex")

    set.seed(.SEED); TRY=0
    while(TRY <= 10 ){
      TRY=TRY+1; mod=NULL; 
      ## cleanup MARK models
      fVec=list.files(full.names=T,pattern=modNo)
      fVec=fVec[substr(fVec,nchar(fVec)-3,nchar(fVec)) %in%
  	  c(".inp",".out",".res",".vcv",".tmp")]
      unlink(fVec)

        mod = try(mark(data = mData, model.name = NULL,
		model = mType,groups=grps, 
		time.intervals=mTime,model.parameters=mParamList,
		prefix = PREFIX,output=F,silent=T) )

      if(!( any(class(mod)=="try-error") | is.null(mod$results) |
	      any(is.na(mod$results$derived.vcv[["N Population Size"]])) |
		all(mod$results$derived[["N Population Size"]]$se==0) |
		any(mod$results$derived[["N Population Size"]]$estimate > nMax) ))
		break
    }
    TRY

  if( ( any(class(mod)=="try-error") | is.null(mod$results) |
	any(is.na(mod$results$derived.vcv[["N Population Size"]])) |
	all(mod$results$derived[["N Population Size"]]$se==0) |
	any(mod$results$derived[["N Population Size"]]$estimate> nMax) ))
	stop("modelError:: gridSex")

  ##################################################################
  ##################################################################
  ## BOOT

  TOL=1e-4
  mVec=mod$results$derived[["N Population Size"]]$estimate
  vMat=mod$results$derived.vcv[["N Population Size"]]  
  if(!lqmm::is.positive.definite(vMat))
    	vMat=as.matrix(Matrix::nearPD(vMat,corr=F,do2eigen = T,
	conv.norm.type="I",doDykstra=T)$mat)
    	## vMat=lqmm::is.positive.definite(x)

  if(!lqmm::is.positive.definite(vMat))
	stop("	vMat not positive.definite")

  bootDF=data.frame(t(MASS::mvrnorm(nBoot,mVec,vMat,tol=TOL)),
                    stringsAsFactors = T)
  names(bootDF)=fn.genLab(1:ncol(bootDF),"boot")

dim(bootDF); length(mVec); dim(vMat)
bMeans= apply(bootDF,1,mean); plot(mVec,bMeans); lines(c(0,100),c(0,100))
bSD= (apply(bootDF,1,sd)^2); 
	plot(diag(mod$results$derived.vcv[["N Population Size"]]),bSD); lines(c(0,100),c(0,100))

  varDF= mod$design.data$S
  varDF=varDF[,sapply(varDF,function(x) length(unique(x))>1)]
  levels(varDF$sex)= c("Female","Male")
  varDF$habitatType= factor(substr(varDF$grid,1,1))

  orig.summDF.A = aggregate( bootDF, 
	by=list(grid=varDF$grid,sex=varDF$sex),mean)
  orig.summDF.A=cbind(cap.hist=T,
	effArea=NA,
	habitatType=factor(substr(orig.summDF.A$grid,1,1)),
	gridSex=factor(paste0(orig.summDF.A$grid,"_",orig.summDF.A$sex)),
	orig.summDF.A)

  if(any(substr(names(varDF),1,4)=="boot"))
	stop("col names can not start with 'boot' ")
  bootCols=which( grepl("boot",names(orig.summDF.A)))
  varCols=which(! grepl("boot",names(orig.summDF.A)))

  xgs = as.character(sapply(levels(orig.summDF.A$sex),function(x)
	paste0(mEffArea$grid,"_",x)))
  xgs = unique(xgs[which( !xgs %in% levels(orig.summDF.A$gridSex))])

  xDF=orig.summDF.A[1:length(xgs),]
  xDF$cap.hist=F;xDF$gridSex=xgs;
  xDF$grid= sapply(xgs,function(x) strsplit(x,"_")[[1]][1])
  xDF$sex = sapply(xgs,function(x) strsplit(x,"_")[[1]][2])
  xDF$habitatType=substr(xDF$grid,1,1)
  xDF[,bootCols]=0
  
  orig.summDF.A=drop.levels(rbind(orig.summDF.A,xDF),reorder=T)
  orig.summDF.A=orig.summDF.A[order(orig.summDF.A$gridSex),]
  rownames(orig.summDF.A)=NULL

  if(!is.null(mEffArea)){
    if(!all(as.character(orig.summDF.A$grid) %in% 
	  as.character(mEffArea$grid)))
        fn.FatalError("check grids in effARea" )

    for( R in 1:nrow(orig.summDF.A)){ ## R=25
	orig.summDF.A$effArea[R] = 
		mEffArea$effArea.km2[
		as.character(mEffArea$grid) == 
		as.character(orig.summDF.A$grid[R]) &
		as.character(mEffArea$sex) == 
		as.character(orig.summDF.A$sex[R]) ]	
    }
  }

  xInvEffArea= matrix(rep(1/orig.summDF.A$effArea,length(bootCols)),
	ncol=length(bootCols))
  orig.summDF.D=orig.summDF.A
  orig.summDF.D[bootCols]=orig.summDF.D[bootCols]*xInvEffArea
  rm(xInvEffArea)

  RES.Abund[["sexGridBootDF"]]	= orig.summDF.A
  RES.Dens[["sexGridBootDF"]]		= orig.summDF.D

  ##################################################################
  ##################################################################
  ## GRID

    summDF = aggregate( orig.summDF.A[,-varCols], 
		by=list(grid=orig.summDF.A$grid,
			  habitatType=orig.summDF.A$habitatType),sum)
    chDF= aggregate( orig.summDF.A$cap.hist, 
		by=list(grid=orig.summDF.A$grid,
			  habitatType=orig.summDF.A$habitatType),sum) 
   
    summDF
    COLS= which(substr(names(summDF),1,4)=="boot")
    RES.Abund[["gridDF"]] = data.frame(
		grid=summDF$grid,
		sex="All",
		habType=summDF$habitatType,
		ch.df= chDF$x/max(chDF$x),
		mean=apply(summDF[,COLS],1,mean),
		se=apply(summDF[,COLS],1,sd),
		lci=apply(summDF[,COLS],1,quantile,0.025),
		uci=apply(summDF[,COLS],1,quantile,0.975) )

  RES.Dens[["gridDF"]]=NA
  if(!is.null(mEffArea)){
    summDF = aggregate( orig.summDF.D[,-varCols], 
		by=list(grid=orig.summDF.D$grid,
			  habitatType=orig.summDF.D$habitatType),sum)
    chDF= aggregate( orig.summDF.D$cap.hist, 
		by=list(grid=orig.summDF.D$grid,
			  habitatType=orig.summDF.D$habitatType),sum)
    COLS= which(substr(names(summDF),1,4)=="boot")
    RES.Dens[["gridDF"]] = data.frame(
		grid=summDF$grid,
		sex="All",
		habType=summDF$habitatType,
		ch.df= chDF$x/max(chDF$x),
		mean=apply(summDF[,COLS],1,mean),
		se=apply(summDF[,COLS],1,sd),
		lci=apply(summDF[,COLS],1,quantile,0.025),
		uci=apply(summDF[,COLS],1,quantile,0.975) )
  }

  ##################################################################
  ##################################################################
  ## SEX

    summDF = aggregate( orig.summDF.A[,-varCols], 
		by=list(sex=orig.summDF.A$sex),mean)
    chDF= aggregate( orig.summDF.A$cap.hist, 
		by=list(sex=orig.summDF.A$sex),mean)
    COLS= which(substr(names(summDF),1,4)=="boot")
    RES.Abund[["sexDF"]] = data.frame(
		grid="All",
		sex=summDF$sex,
		habType="All",
		ch.df= chDF$x/max(chDF$x),
		mean=apply(summDF[,COLS],1,mean),
		se=apply(summDF[,COLS],1,sd),
		lci=apply(summDF[,COLS],1,quantile,0.025),
		uci=apply(summDF[,COLS],1,quantile,0.975) )

  RES.Dens[["sexDF"]]=NA
  if(!is.null(mEffArea)){
    summDF = aggregate( orig.summDF.D[,-varCols], 
		by=list(sex=orig.summDF.D$sex),mean)
    chDF= aggregate( orig.summDF.D$cap.hist, 
		by=list(sex=orig.summDF.D$sex),mean)
    COLS= which(substr(names(summDF),1,4)=="boot")
    RES.Dens[["sexDF"]] = data.frame(
		grid="All",
		sex=summDF$sex,
		habType="All",
		ch.df= chDF$x/max(chDF$x),
		mean=apply(summDF[,COLS],1,mean),
		se=apply(summDF[,COLS],1,sd),
		lci=apply(summDF[,COLS],1,quantile,0.025),
		uci=apply(summDF[,COLS],1,quantile,0.975) )
  }

  ##################################################################
  ##################################################################
  ## HABITAT

    summDF = aggregate( orig.summDF.A[,-varCols], 
		by=list(grid=orig.summDF.A$grid,
			  habitatType=orig.summDF.A$habitatType),sum)
    chDF= aggregate( orig.summDF.A$cap.hist, 
		by=list(grid=orig.summDF.A$grid,
			  habitatType=orig.summDF.A$habitatType),sum)    

    COLS= which(substr(names(summDF),1,4)=="boot")
    summDF = summDF.H.A = aggregate( summDF[,COLS], 
		by=list(habitatType=summDF$habitatType),mean)
    chDF= aggregate( chDF$x, 
		by=list(habitatType=chDF$habitatType),mean)
  
    COLS= which(substr(names(summDF),1,4)=="boot")
    RES.Abund[["habDF"]] = data.frame(
		grid="All",
		sex="All",
		habType=summDF$habitatType,
		ch.df= chDF$x/max(chDF$x),
		mean=apply(summDF[,COLS],1,mean),
		se=apply(summDF[,COLS],1,sd),
		lci=apply(summDF[,COLS],1,quantile,0.025),
		uci=apply(summDF[,COLS],1,quantile,0.975) )

    
    ##----------------------------------------------------------------
    summDF=summDF[rev(order(apply(summDF[,COLS],1,mean) )),]
    COLS= which(substr(names(summDF),1,4)=="boot")
    summDF.diff=summDF.ratio=summDF.pc=NULL
    for(R1 in 1:(nrow(summDF)-1)){
      for(R2 in (R1+1):(nrow(summDF))){
  	  ## R1=1;R2=2
        NAM= paste0(summDF$habitatType[R1],"-",
			 summDF$habitatType[R2])
        VAL=summDF[R1,COLS]-summDF[R2,COLS]
        summDF.diff=rbind(summDF.diff,cbind( stat = NAM, VAL))

        NAM= paste0("100*(", summDF$habitatType[R1],"-",
			 summDF$habitatType[R2],")/",
			 summDF$habitatType[R2])
        VAL=100*(summDF[R1,COLS]-summDF[R2,COLS])/summDF[R2,COLS]
        summDF.pc=rbind(summDF.pc,cbind( stat = NAM, VAL))

        NAM= paste0(summDF$habitatType[R1],"/",
			 summDF$habitatType[R2])
        VAL=summDF[R1,COLS]/summDF[R2,COLS]
        summDF.ratio=rbind(summDF.ratio,cbind( stat = NAM, VAL))
      }
    }

    summDF=summDF.diff
    COLS= which(substr(names(summDF),1,4)=="boot")
    RES.Abund[["habDiffDF"]] = data.frame(
		grid="All",
		sex="All",
		habType=summDF$stat,
		ch.df= NA,
		mean=apply(summDF[,COLS],1,mean),
		se=apply(summDF[,COLS],1,sd),
		lci=apply(summDF[,COLS],1,quantile,0.025),
		uci=apply(summDF[,COLS],1,quantile,0.975) )

    summDF=summDF.pc
    COLS= which(substr(names(summDF),1,4)=="boot")
    RES.Abund[["habPCDiffDF"]] = data.frame(
		grid="All",
		sex="All",
		habType=summDF$stat,
		ch.df= NA,
		mean=apply(summDF[,COLS],1,mean),
		se=apply(summDF[,COLS],1,sd),
		lci=apply(summDF[,COLS],1,quantile,0.025),
		uci=apply(summDF[,COLS],1,quantile,0.975) )

    summDF=summDF.ratio
    COLS= which(substr(names(summDF),1,4)=="boot")
    RES.Abund[["habRatioDF"]] = data.frame(
		grid="All",
		sex="All",
		habType=summDF$stat,
		ch.df= NA,
		mean=apply(summDF[,COLS],1,mean),
		se=apply(summDF[,COLS],1,sd),
		lci=apply(summDF[,COLS],1,quantile,0.025),
		uci=apply(summDF[,COLS],1,quantile,0.975) )


  RES.Dens[["habDF"]]=RES.Dens[["habRatioDF"]]=
	RES.Dens[["habPCDiffDF"]]=RES.Dens[["habDiffDF"]]=NA
  if(!is.null(mEffArea)){
    summDF = aggregate( orig.summDF.D[,-varCols], 
		by=list(grid=orig.summDF.D$grid,
			  habitatType=orig.summDF.D$habitatType),sum)
    chDF= aggregate( orig.summDF.D$cap.hist, 
		by=list(grid=orig.summDF.D$grid,
			  habitatType=orig.summDF.D$habitatType),sum)    

    COLS= which(substr(names(summDF),1,4)=="boot")
    summDF = summDF.H.D = aggregate( summDF[,COLS], 
		by=list(habitatType=summDF$habitatType),mean)
    chDF= aggregate( chDF$x, 
		by=list(habitatType=chDF$habitatType),mean)
  
    COLS= which(substr(names(summDF),1,4)=="boot")
    RES.Dens[["habDF"]] = data.frame(
		grid="All",
		sex="All",
		habType=summDF$habitatType,
		ch.df= chDF$x/max(chDF$x),
		mean=apply(summDF[,COLS],1,mean),
		se=apply(summDF[,COLS],1,sd),
		lci=apply(summDF[,COLS],1,quantile,0.025),
		uci=apply(summDF[,COLS],1,quantile,0.975) )

    ##----------------------------------------------------------------
    summDF=summDF[rev(order(apply(summDF[,COLS],1,mean) )),]
    COLS= which(substr(names(summDF),1,4)=="boot")
    summDF.diff=summDF.ratio=summDF.pc=NULL
    for(R1 in 1:(nrow(summDF)-1)){
      for(R2 in (R1+1):(nrow(summDF))){
  	  ## R1=1;R2=2
        NAM= paste0(summDF$habitatType[R1],"-",
			 summDF$habitatType[R2])
        VAL=summDF[R1,COLS]-summDF[R2,COLS]
        summDF.diff=rbind(summDF.diff,cbind( stat = NAM, VAL))

        NAM= paste0("100*(", summDF$habitatType[R1],"-",
			 summDF$habitatType[R2],")/",
			 summDF$habitatType[R2])
        VAL=100*(summDF[R1,COLS]-summDF[R2,COLS])/summDF[R2,COLS]
        summDF.pc=rbind(summDF.pc,cbind( stat = NAM, VAL))

        NAM= paste0(summDF$habitatType[R1],"/",
			 summDF$habitatType[R2])
        VAL=summDF[R1,COLS]/summDF[R2,COLS]
        summDF.ratio=rbind(summDF.ratio,cbind( stat = NAM, VAL))
      }
    }

    summDF=summDF.diff
    COLS= which(substr(names(summDF),1,4)=="boot")
    RES.Dens[["habDiffDF"]] = data.frame(
		grid="All",
		sex="All",
		habType=summDF$stat,
		ch.df= NA,
		mean=apply(summDF[,COLS],1,mean),
		se=apply(summDF[,COLS],1,sd),
		lci=apply(summDF[,COLS],1,quantile,0.025),
		uci=apply(summDF[,COLS],1,quantile,0.975) )

    summDF=summDF.pc
    COLS= which(substr(names(summDF),1,4)=="boot")
    RES.Dens[["habPCDiffDF"]] = data.frame(
		grid="All",
		sex="All",
		habType=summDF$stat,
		ch.df= NA,
		mean=apply(summDF[,COLS],1,mean),
		se=apply(summDF[,COLS],1,sd),
		lci=apply(summDF[,COLS],1,quantile,0.025),
		uci=apply(summDF[,COLS],1,quantile,0.975) )

    summDF=summDF.ratio
    COLS= which(substr(names(summDF),1,4)=="boot")
    RES.Dens[["habRatioDF"]] = data.frame(
		grid="All",
		sex="All",
		habType=summDF$stat,
		ch.df= NA,
		mean=apply(summDF[,COLS],1,mean),
		se=apply(summDF[,COLS],1,sd),
		lci=apply(summDF[,COLS],1,quantile,0.025),
		uci=apply(summDF[,COLS],1,quantile,0.975) )
  }

  

  ##################################################################
  ##################################################################
  ## HABITAT & SEX

    summDF = summDF.HS.A = aggregate( orig.summDF.A[,-varCols], 
		by=list(sex=orig.summDF.A$sex,
			  habitatType=orig.summDF.A$habitatType),mean)
    chDF= aggregate( orig.summDF.A$cap.hist, 
		by=list(sex=orig.summDF.A$sex,
			  habitatType=orig.summDF.A$habitatType),mean)    

  
    COLS= which(substr(names(summDF),1,4)=="boot")
    RES.Abund[["sexHabDF"]] = data.frame(
		grid="All",
		sex=summDF$sex,
		habType=summDF$habitatType,
		ch.df= chDF$x/max(chDF$x),
		mean=apply(summDF[,COLS],1,mean),
		se=apply(summDF[,COLS],1,sd),
		lci=apply(summDF[,COLS],1,quantile,0.025),
		uci=apply(summDF[,COLS],1,quantile,0.975) )

  RES.Dens[["sexHabDF"]]=NA
  if(!is.null(mEffArea)){
    summDF = summDF.HS.D = aggregate( orig.summDF.D[,-varCols], 
		by=list(sex=orig.summDF.D$sex,
			  habitatType=orig.summDF.D$habitatType),mean)
    chDF= aggregate( orig.summDF.D$cap.hist, 
		by=list(sex=orig.summDF.D$sex,
			  habitatType=orig.summDF.D$habitatType),mean)    

  
    COLS= which(substr(names(summDF),1,4)=="boot")
    RES.Dens[["sexHabDF"]] = data.frame(
		grid="All",
		sex=summDF$sex,
		habType=summDF$habitatType,
		ch.df= chDF$x/max(chDF$x),
		mean=apply(summDF[,COLS],1,mean),
		se=apply(summDF[,COLS],1,sd),
		lci=apply(summDF[,COLS],1,quantile,0.025),
		uci=apply(summDF[,COLS],1,quantile,0.975) )
  }


  ##############################################################
  ## z.test of differences

  summDF.H.A=cbind(sex=factor("All"),summDF.H.A)
  summDF.H.D=cbind(sex=factor("All"),summDF.H.D)

  t1= expand.grid(
	  h1=levels(summDF.H.A$habitatType),
	  s1=levels(summDF.H.A$sex),
 	  h2=levels(summDF.H.A$habitatType),
	  s2=levels(summDF.H.A$sex),stringsAsFactors = T )	
  
  t1$h1=factor(as.character(t1$h1),c("B","R","W","P"))
  t1$h2=factor(as.character(t1$h2),c("B","R","W","P"))

  t1$lab1=factor(paste0( "(",t1$h1," ",t1$s1,")"))
  t1$lab2=factor(paste0( "(",t1$h2," ",t1$s2,")"))
  ROWS= which( t1$lab1 != t1$lab2 &
		 (t1$h1 == t1$h2 | t1$s1 == t1$s2))
  t1=t1[ROWS,]
  t1=t1[order(t1$h1,t1$h2),]
  ROWS= which(as.numeric(t1$h1) < as.numeric(t1$h2))
  t1=t1[ROWS,]
  t1$contrast= paste0(t1$lab1, " - ", t1$lab2)
  t1$contrast=gsub("(","",gsub(")","",gsub(" All","",
		  gsub("B ","Bottomland ",
		  gsub("P ","Pine ",
		  gsub("R ","Riparian ",
		  gsub("W ","Wetland ",t1$contrast))))),fixed=T),fixed=T)
  t1$p.value=t1$z=t1$se=t1$mean=NA

  resA=resD=t1; bootCols=which(grepl("boot",names(summDF.H.A)))
  for(R in 1:nrow(t1)){ ## R=5
    r1=which(summDF.H.A$habitat==as.character(t1$h1[R]) &
		 summDF.H.A$sex==as.character(t1$s1[R]) )
    r2=which(summDF.H.A$habitat==as.character(t1$h2[R]) &
		 summDF.H.A$sex==as.character(t1$s2[R]) )

    

    m1=mean(as.numeric(summDF.H.A[r1,bootCols]))
    m2=mean(as.numeric(summDF.H.A[r2,bootCols]))
    se1=sd(as.numeric(summDF.H.A[r1,bootCols]))
    se2=sd(as.numeric(summDF.H.A[r2,bootCols]))
    resA$mean[R]= m1-m2; resA$se[R]= sqrt( (se1^2)+(se2^2) ); 
    resA$z[R]=((m1-m2)-0)/sqrt( (se1^2)+(se2^2) )
    resA$p.value[R]= 2*pnorm(-abs(resA$z[R]))

    m1=mean(as.numeric(summDF.H.D[r1,bootCols]))
    m2=mean(as.numeric(summDF.H.D[r2,bootCols]))
    se1=sd(as.numeric(summDF.H.D[r1,bootCols]))
    se2=sd(as.numeric(summDF.H.D[r2,bootCols]))
    resD$mean[R]= m1-m2; resD$se[R]= sqrt( (se1^2)+(se2^2) ); 
    resD$z[R]=((m1-m2)-0)/sqrt( (se1^2)+(se2^2) )
    resD$p.value[R]= 2*pnorm(-abs(resD$z[R]))
  }
  COLS = which(! names(resA) %in% c("h1","h2","s1","s2","lab1","lab2"))
  if(!is.na(pAdjMethod)) resA$p.value= p.adjust(resA$p.value,pAdjMethod)
  if(!is.na(pAdjMethod)) resD$p.value= p.adjust(resD$p.value,pAdjMethod)
  round(resD$p.value,3)


  RES.Abund[[ "contrasts.habitat" ]]=resA[,COLS]
  RES.Dens[[ "contrasts.habitat" ]]=resD[,COLS]
  rm(t1,resA,resD,COLS)

  ##--------------------------------------------------------------

  t1= expand.grid(
	  h1=levels(summDF.HS.A$habitatType),
	  s1=levels(summDF.HS.A$sex),
	  h2=levels(summDF.HS.A$habitatType),
	  s2=levels(summDF.HS.A$sex), stringsAsFactors = T)	
  t1$h1=factor(as.character(t1$h1),c("B","R","W","P"))
  t1$h2=factor(as.character(t1$h2),c("B","R","W","P"))
  t1$s1=factor(as.character(t1$s1),c("Male","Female"))
  t1$s2=factor(as.character(t1$s2),c("Male","Female"))

  t1$lab1=factor(paste0( "(",t1$h1," ",t1$s1,")"))
  t1$lab2=factor(paste0( "(",t1$h2," ",t1$s2,")"))
  ROWS= which( t1$lab1 != t1$lab2 &
		 (t1$h1 == t1$h2))
  t1=t1[ROWS,]
  t1=t1[order(t1$s1,t1$h1),]
  ROWS= which(as.numeric(t1$s1) < as.numeric(t1$s2))
  t1=t1[ROWS,]
  t1$contrast= paste0(t1$lab1, " - ", t1$lab2)
  t1$contrast=
		  gsub("B ","Bottomland ",
		  gsub("P ","Pine ",
		  gsub("R ","Riparian ",
		  gsub("W ","Wetland ",t1$contrast))))
  t1$p.value=t1$z=t1$se=t1$mean=NA

  resA=resD=t1; bootCols=which(grepl("boot",names(summDF.HS.A)))
  for(R in 1:nrow(t1)){ ## R=1
    r1=which(summDF.HS.A$habitat==as.character(t1$h1[R]) &
		 summDF.HS.A$sex==as.character(t1$s1[R]) )
    r2=which(summDF.HS.A$habitat==as.character(t1$h2[R]) &
		 summDF.HS.A$sex==as.character(t1$s2[R]) )

    m1=mean(as.numeric(summDF.HS.A[r1,bootCols]))
    m2=mean(as.numeric(summDF.HS.A[r2,bootCols]))
    se1=sd(as.numeric(summDF.HS.A[r1,bootCols]))
    se2=sd(as.numeric(summDF.HS.A[r2,bootCols]))
    resA$mean[R]= m1-m2; resA$se[R]= sqrt( (se1^2)+(se2^2) ); 
    resA$z[R]=((m1-m2)-0)/sqrt( (se1^2)+(se2^2) )
    resA$p.value[R]= 2*pnorm(-abs(resA$z[R]))

    m1=mean(as.numeric(summDF.HS.D[r1,bootCols]))
    m2=mean(as.numeric(summDF.HS.D[r2,bootCols]))
    se1=sd(as.numeric(summDF.HS.D[r1,bootCols]))
    se2=sd(as.numeric(summDF.HS.D[r2,bootCols]))
    resD$mean[R]= m1-m2; resD$se[R]= sqrt( (se1^2)+(se2^2) ); 
    resD$z[R]=((m1-m2)-0)/sqrt( (se1^2)+(se2^2) )
    resD$p.value[R]= 2*pnorm(-abs(resD$z[R]))

  }
  COLS = which(! names(resA) %in% c("h1","h2","s1","s2","lab1","lab2"))
  if(!is.na(pAdjMethod)) resA$p.value= p.adjust(resA$p.value,pAdjMethod)
  if(!is.na(pAdjMethod)) resD$p.value= p.adjust(resD$p.value,pAdjMethod)

  RES.Abund[[ "contrasts.sex" ]]=resA[,COLS]
  RES.Dens[[ "contrasts.sex" ]]=resD[,COLS]
  rm(t1,resA,resD,COLS)

  ##--------------------------------------------------------------

  RES = list(RES.Abund=RES.Abund,RES.Dens=RES.Dens)
  names(RES[[2]])
  invisible(RES)

  ## RES.Dens$habSexDF$mean[1:4]+RES.Dens$habSexDF$mean[5:8]
  ## RES.Dens$habDF$mean

}



############################################################################
############################################################################
if(1==2){
  nBoot=1000; mod= bestMod; addVars= addVars;  effArea=effArea
}

fn.parametricBoot = function(nBoot, mod,addVars, effArea=NA){
  ## nBoot=10; mod=bestMod; addVars=c("grid","session"); effArea=NA

  TOL=1e-4
  mVec=mod$results$derived[["N Population Size"]]$estimate
  vMat=mod$results$derived.vcv[["N Population Size"]]  
  if(!lqmm::is.positive.definite(vMat))
    vMat=lqmm::make.positive.definite(vMat, tol=TOL)


  bootDF=data.frame(t(MASS::mvrnorm(nBoot,mVec,vMat,tol=TOL)))
  names(bootDF)=fn.genLab(1:ncol(bootDF),"boot")

  varDF= mod$design.data$f0
  varDF=varDF[,sapply(varDF,function(x) length(unique(x))>1)]
  if(any(substr(names(varDF),1,4)=="boot"))
	stop("col names can not start with 'boot' ")

  addVars=addVars[addVars %in% names(varDF)]
  if(length(addVars)==0)stop("add vars not in data")


  ## dim(bootDF);dim(varDF); names(varDF)
  
  varDF$summVar = factor(sapply(1:nrow(varDF),function(x) ## x=1
	paste( sapply(varDF[x,addVars],function(y) as.character(y)),
	collapse="_") ))
  
  BY=lapply(addVars,function(x) varDF[,x]);names(BY)=addVars
  summDF = aggregate( bootDF, by=BY,sum)
  
  summDF$siteType = factor(substr(as.character(summDF$grid),1,1))
  
  COLS= names(summDF)[which(substr(names(summDF),1,4)=="boot")]
  gridDF=aggregate( summDF[,COLS], by=list(grid=summDF$grid),mean)
  habDF=aggregate( summDF[,COLS], by=list(siteType=summDF$siteType),mean)

  RES.Abund= list(
  	boot.abund.gridDF= data.frame(
		grid=gridDF$grid,
		mean=apply(gridDF[,COLS],1,mean),
		se=apply(gridDF[,COLS],1,sd),
		lci=apply(gridDF[,COLS],1,quantile,0.025),
		uci=apply(gridDF[,COLS],1,quantile,0.975) ),

  	boot.abund.habDF= data.frame(
		siteType=habDF$siteType,
		mean=apply(habDF[,COLS],1,mean),
		se=apply(habDF[,COLS],1,sd),
		lci=apply(habDF[,COLS],1,quantile,0.025),
		uci=apply(habDF[,COLS],1,quantile,0.975) )
  )

  RES.Dens = list(
  	boot.dens.gridDF= NA,
  	boot.dens.habDF= NA 
  )

  if(!is.na(effArea)){
    
    gridDF[,COLS]=gridDF[,COLS]/effArea
    habDF[,COLS] =habDF[,COLS]/effArea
    RES.Dens= list(
  	boot.dens.gridDF= data.frame(
		grid=gridDF$grid,
		mean=apply(gridDF[,COLS],1,mean),
		se=apply(gridDF[,COLS],1,sd),
		lci=apply(gridDF[,COLS],1,quantile,0.025),
		uci=apply(gridDF[,COLS],1,quantile,0.975) ),

  	boot.dens.habDF= data.frame(
		siteType=habDF$siteType,
		mean=apply(habDF[,COLS],1,mean),
		se=apply(habDF[,COLS],1,sd),
		lci=apply(habDF[,COLS],1,quantile,0.025),
		uci=apply(habDF[,COLS],1,quantile,0.975) )
    )
  }

  c(RES.Abund,RES.Dens)

}


############################################################################
############################################################################


fn.genLab=function(labs,str="lab"){ ## labs= c(1:100); str="lab"
  maxChar=ifelse(max(nchar(labs))==1,2,max(nchar(labs)))
  paste0(str,formatC(labs,width=maxChar,flag = "0"))
}


############################################################################
############################################################################

fn.initPar=function(){
  par(mar=c(2.5,2.5,1,1),mgp=c(1.25,0,0),tck=0.02,cex=1,mex=1,cex.axis=0.8)
}

############################################################################
############################################################################

fn.line2user <- function(line, side) {
  ## copied from https://stackoverflow.com/questions/29125019/
	## get-margin-line-locations-mgp-in-user-coordinates
  lh <- par('cin')[2] * par('cex') * par('lheight')
  x_off <- diff(grconvertX(0:1, 'inches', 'user'))
  y_off <- diff(grconvertY(0:1, 'inches', 'user'))
  switch(side,
         `1` = par('usr')[3] - line * y_off * lh,
         `2` = par('usr')[1] - line * x_off * lh,
         `3` = par('usr')[4] + line * y_off * lh,
         `4` = par('usr')[2] + line * x_off * lh,
         stop("side must be 1, 2, 3, or 4", call.=FALSE))
}

############################################################################
############################################################################

fn.MakeTransparent=function(color,alpha=0.25){ ## color=seasun(3)[2];alpha=0.25
  if(alpha>1)alpha=alpha/255; color
  for(r in 1:length(color)){ ## r=color[1]
    colorAlpha=as.vector(col2rgb(color[r],T)/255);colorAlpha[4]=alpha
    colorAlpha=rgb(colorAlpha[1],colorAlpha[2],colorAlpha[3],colorAlpha[4])
    color[r]=colorAlpha
  }
  color
}

############################################################################
############################################################################

fn.FatalError<-function(msg="Error"){

  cat("\n","There is an ERRORRRRRRRRR........","\n")
  i<-0
  pb<-txtProgressBar(min=0,max=10,initial=0,char="=",width=25,label="Error",style=3)
  while(i>=0){ ## i<-1
   setTxtProgressBar(pb,i)
   if(i==9){
     cat(msg)
     i<-1
    } else i<-i+1
  }
}


#############################################################################
#############################################################################

fn.MakeCombos<-function(MainVars,Interactions=T){
  AllCombos<-Params<-NULL; PSno<-1  
  for(c in 1:length(MainVars)){##c<-2  
    Combn<-combn(MainVars, c,simplify=T)  
    Combn<-sapply(1:ncol(Combn),function(x)paste(Combn[,x],collapse=":",sep=""))  
    Params<-c(Params,Combn)  
  }
  if(Interactions){
    for(c in 1:length(Params)){##c<-1  
      Combn<-combn(Params, c,simplify=F)  
      AllCombos<-c(AllCombos,Combn)  
    }  
    UseComboVec<-rep(T,length(AllCombos))  
    for(c in 1:length(AllCombos)){##c<-127  
      for(p in 1:length(AllCombos[[c]])){##p<-1  
        Vars<-unlist(strsplit(AllCombos[[c]][p],":"))  
        if(length(Vars)>1){  
          cVars<-NULL  
          for(v in 1:(length(Vars)-1)){##v<-1  
            Combn<-combn(Vars, v,simplify=T)  
            Combn<-sapply(1:ncol(Combn),function(x)paste(Combn[,x],  
      		collapse=":",sep=""))
            cVars<-c(cVars,Combn)  
          } ## end for v  
          if(!sum(sapply(cVars,function(x)(
			length(which(AllCombos[[c]][-p]==x))>0)))==  
    		length(cVars))UseComboVec[c]<-F    
        }## end if  
      }  ## end for p  
    } ##end for c  
    AllCombos<-AllCombos[UseComboVec]  
    AllCombos<-sapply(1:length(AllCombos),function(x)
			paste(AllCombos[[x]],collapse="+",sep=""))
  } else {
	AllCombos=as.vector(sapply(Params,function(x)gsub(":","+",x)))
  }
  AllCombos
}##end function
##MakeCombos(MainVars)

