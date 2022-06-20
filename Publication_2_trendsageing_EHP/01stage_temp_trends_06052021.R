###########################################################################
############## AGE SPECIFIC TIME VAR FINAL  ###############################
########################################################################

# LOAD PACKAGES
Sys.setenv(LANG = "en")
library(data.table); library(dlnm); library(mixmeta) ; library(gnm)
library(tsModel)  ; library(splines) ; library(mgcv)
library(reshape2); library (dplyr); library(lubridate)

# SET DIRECTORY OF DATA INPUT
dir <- "Y:/CCH/Evan/temporal_trends_paper"
dir <- "~/Desktop/PhD/publication1/final code/clean"

# LOAD DATA
munciptable <- readRDS(paste0(dir,"/municipalitylookuptable.rds"))
deathrecords <- readRDS(paste0(dir,"/deathrecordsmuncipality.rds"))

# SAVE ONLY DEATHS BEFORE 2018
deathrecords <- deathrecords[(deathrecords$yy=="2018"),]
deathrecords$NAME <- as.character(deathrecords$NAME)

deathrecordst <- deathrecords %>%
  dplyr::group_by(NAME) %>%
  summarise(deaths=sum(deaths))

# CREATE AGE CATE3GORY
deathrecords$age_cat <- NA
deathrecords$age_cat <- deathrecords$age_death
deathrecords$age_cat <- ifelse(deathrecords$age_cat<=64,"1", ifelse(deathrecords$age_cat<=79,"2", ifelse(deathrecords$age_cat<=150, "3")))

# CREATE LOOKUP TABLE
lookuptable <- deathrecords %>%
  group_by(NAME,comm_resi, assigned_number,Kantonname) %>%
  summarise(KANTONSNUM = mean(KANTONSNUM), deaths=sum(deaths))

lookuptable <- lookuptable[order(lookuptable$assigned_number),]

listcantons <- order(unique(lookuptable$KANTONSNUM))
listcantonsname <- (unique(lookuptable$Kantonname))
#26 cantons

# SPECIFICATION OF THE EXPOSURE FUNCTION
varfun <- "bs"
vardegree <- 2
varper <- c(10, 75,90)

# SPECIFICATION OF THE LAG FUNCTION
lag <- 21
lagnk <- 3

# DEFINE PERIODS
tabperiod <- matrix(NA, nrow=5, ncol=2)
tabperiod[,1] <- seq(1969,2017, by=10)
tabperiod[,2] <- tabperiod[,1] + 9
tabperiod[5,2] <- 2017  ## LAST YEAR IS 2017

# MATRIX OF MEAN/IQR TEAMP PER CANTON & PERIOD
avertmean_sp <- matrix(NA, nrow=5*length(listcantons), ncol=1)
rangetmean_sp <- matrix(NA, nrow=5*length(listcantons), ncol=1)

# TEMPERATURE DISTRIBUTION
tmeanper_list <- list()

## 1) SELECT AGE GROUP
for (a in c(1:3)){
  
  # TELL ME THE ITERATION
  cat("",a)
  
  # SELECT AGE GROUP
  deathrecordsage <- subset(deathrecords, age_cat == a)
  
  # CREATE EMPTY OBJECTS FOR STORING RESULTS - BY AGE CATEGORY
  ncoef <- length(varper) + ifelse(varfun=="bs",vardegree,1)
  coefall <- matrix(NA,(length(listcantons)*nrow(tabperiod)),ncoef+3,
                    dimnames=list(rep(names(listcantons),each=5)))
  
  coefall[,1] <- rep(1:26,  each=5)
  coefall[,2] <- rep(1:5,26)
  coefall[,3] <- rep(a, each=length(listcantons)*5)
  vcovall <- vector("list",length(listcantons)*5)
  names(vcovall) <- rep(listcantonsname, each=5)
  
  redall <- list()
  
  ## 2)  SELECT CANTON
  for(i in seq(length(listcantons))) {
    
    # TELL ME THE ITERATION
    cat("",i)
    
    # SELECT MUNICPIALITIES FOR REGION i
    listmuncip <- lookuptable$assigned_number[lookuptable$KANTONSNUM==i]
    
    # EXTRACT THE MUNCIPALITY DATA AND APPEND THEM
    data <- NULL
    for(j in seq(length(listmuncip))) {
      temp <- data.table(readRDS(paste0(inpath, paste0("/temp.muncip",listmuncip[j]), ".rds")))
      temp <- cbind(data.table(IDmuncip=listmuncip[j]), temp)
      deathrecordssel <- subset(deathrecordsage, assigned_number==listmuncip[j])
      deathrecordssel <- dplyr::select(deathrecordssel,date,deaths)
      deathrecordssel <- deathrecordssel %>% group_by(date) %>%
        summarise(deaths=sum(deaths), .groups = 'drop')
      tempmerged <- merge(temp, deathrecordssel, by="date", all=TRUE)
      data <- rbind(data, tempmerged)
    }
    
    # TRANSFORM NA TO 0
    data <- data %>%
      mutate(deaths = replace(deaths,is.na(deaths),0))
    
    # GENERATE VARIABLES FOR THE ANALYSIS
    data$time <- as.numeric(data$date)
    data$year <- factor(as.numeric(format(data$date, "%Y")))
    data$month <- factor(as.numeric(format(data$date, "%m")))
    data$dow <- factor(factor(wday(data$date)))
    
    # COMPUTE TEMPERATURE PERCENTILES
    predper <- c(seq(0,1,0.1),2:98,seq(99,100,0.1))
    tmeanper_list[[i]] <- quantile(data$aver.temp, predper/100, na.rm=T)
    
    # CREATE THE CB of temperature
    varknots <- quantile(data$aver.temp, varper/100, na.rm=T)
    boundknt <- range(data$aver.temp)
    lagknots <- logknots(lag,lagnk)
    cbtemp <- crossbasis(data$aver.temp, lag=lag,
                         argvar=list(fun=varfun, degree=vardegree, knots=varknots,
                                     Boundary.knots=boundknt),
                         arglag=list(knots=lagknots),
                         group = factor(data$IDmuncip))
    
    # CREATE STRATUM VARIABLE
    data$stratum <- with(data, factor(IDmuncip):year:month:dow)
    
    # ELIMINATE EMPTY STRATA (FOR CORRECT COMPUTATION OF CI IN gnm)
    ind <- tapply(data$deaths, data$stratum, sum)[data$stratum]   
    
    ## ) SELECT DECADE
    for (k in 1:5){
      
      cat("",k)
      
      # SELECT SUBPERIOD
      data.sub <- subset(data, year %in% c(as.numeric(tabperiod[k,1]):as.numeric(tabperiod[k,2])) ) ## dont run , outside
      
      # STORE AVERAGE AND IQR PER SUBPERIOD
      avertmean_sp[((i-1)*5)+k,1] <- mean(data.sub$aver.temp, na.rm=T)
      rangetmean_sp[((i-1)*5)+k,1] <- IQR(data.sub$aver.temp, na.rm=T)
      
      # # CREATE STRATUM VARIABLE
      # data.sub$stratum <- with(data.sub, factor(IDmuncip):year:month:dow)
      # 
      # # ELIMINATE EMPTY STRATA (FOR CORRECT COMPUTATION OF CI IN gnm)
      # ind <- tapply(data.sub$deaths, data.sub$stratum, sum)[data.sub$stratum]
      
      # cbtemp <- crossbasis(data.sub$aver.temp, lag=lag,
      #                      argvar=list(fun=varfun, degree=vardegree, knots=varknots,
      #                                  Boundary.knots=boundknt),
      #                      arglag=list(knots=lagknots),
      #                      group = factor(data.sub$IDmuncip))
      # 
      
      # RUN MODEL THE MODEL
      mod <- gnm(deaths ~ cbtemp, eliminate=stratum,
                 family=quasipoisson(), data=data, na.action="na.exclude",
                 subset=ind>0 & year %in% c(as.numeric(tabperiod[k,1]):as.numeric(tabperiod[k,2]))) 
      
      # REDUCTION TO OVERALL CUMULATIVE
      redall[[((i-1)*5)+k]]<- crossreduce(cbtemp,mod,cen=mean(data.sub$aver.temp,na.rm=T))
      coefall[((i-1)*5)+k,-c(1:3)] <- coef(redall[[((i-1)*5)+k]])
      vcovall[[((i-1)*5)+k]] <- vcov(redall[[((i-1)*5)+k]])
    }
  }
  
  # EXTRACT OBJECTS WITH RESUTLS
  assign(paste0("coefall", a), coefall)
  assign(paste0("vcovall", a), vcovall)   
  assign(paste0("redall", a), redall)   
  
}

# SAVE SPACE 
rm(list=ls()[! ls() %in% c("coefall1", "coefall2", "coefall3",
                           "vcovall1", "vcovall2", "vcovall3", 
                          
                           "avertmean_sp", "rangetmean_sp", "tmeanper_list")])





## save environment
load("~/Desktop/PhD/publication1/final code/clean/01stagedata_allyecb.Rdata")

