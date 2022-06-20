
Sys.setenv(LANG = "en")

library(dlnm) ; library(mixmeta) ; library(splines) ; library(tsModel)
library(mgcv) 
library(data.table); library(reshape2); library(ggplot2)
library(dplyr); library(GGally); library(RColorBrewer); library(cowplot)


#####
#QAIC function
QAIC <- function(model){
  phi <- summary(model)$dispersion
  loglik <-sum(dpois(model$y, model$fitted.values, log=TRUE))
  return(-2*loglik+2*summary(model)$df[3]*phi)
}

load("data_UK_part.Rdata")
fulldataset <- data_fin

plot.ts(data_fin$East$all)
plot.ts(data_fin$East$ERA5)

##########################################
###### UK GCDS ##########################
#########################################
## DEFINITON OF THE PARAMETERS
# SPECIFICATION OF THE EXPOSURE FUNCTION
varfun <- "bs"
vardegree <- 2
varper <- c(10,75,90)

# SPECIFICATION OF THE LAG FUNCTION
lag <- 21
lagnk <- 3

# DEGREE OF FREEDOM FOR SEASONALITY
dfseas <- 8

# DEFINE THE OUTCOME
out <- "all"

# MODEL FORMULA
formula <- all ~ cb + dow + ns(date,df=round(dfseas*length(date)/365.25))

# DEFINE VECTOR WITH NAMES (ORDERED AS THE FULLDATASET!!!!!)
namesreg <- names(fulldataset)

## 1ST STAGE
# CREATE EMPTY OBJECTS TO STORE RESULTS
ncoef <- length(varper) + ifelse(varfun=="bs",vardegree,1)

coefallinsitu <- matrix(NA,length(namesreg),ncoef,dimnames=list(namesreg))
vcovallinsitu <- vector("list",length(namesreg))

coefallhaduk5 <- matrix(NA,length(namesreg),ncoef,dimnames=list(namesreg))
vcovallhaduk5 <- vector("list",length(namesreg))

coefallhaduk5_unw <- matrix(NA,length(namesreg),ncoef,dimnames=list(namesreg))
vcovallhaduk5_unw <- vector("list",length(namesreg))

coefallERA5 <- matrix(NA,length(namesreg),ncoef,dimnames=list(namesreg))
vcovallERA5 <- vector("list",length(namesreg))

coefallERA5_unw <- matrix(NA,length(namesreg),ncoef,dimnames=list(namesreg))
vcovallERA5_unw <- vector("list",length(namesreg))

qaicres <- matrix(NA, ncol=5, nrow=10)
rownames(qaicres) <- namesreg
names(vcovallinsitu) <- names(vcovallhaduk5) <- names(vcovallhaduk5_unw) <- names(vcovallERA5) <- names(vcovallERA5_unw) <- namesreg

redallinsitu <- redallhaduk5 <- redallhaduk5_unw <- redallERA5<- redallERA5_unw <- list()

# TABLE HEAT
tab_heat <- matrix(NA, nrow=10, ncol=5)
rownames(tab_heat) <- namesreg


tab_heat_insitu <- tab_heat_Haduk5 <- tab_heat_Haduk5_unw <- tab_heat_ERA5 <- tab_heat_ERA5_unw <- matrix(NA, nrow=length(fulldataset), ncol=3)
rownames(tab_heat_insitu) <- rownames(tab_heat_Haduk5) <- rownames(tab_heat_Haduk5_unw) <- rownames(tab_heat_ERA5) <- rownames(tab_heat_ERA5_unw) <- rep(names(fulldataset))
tab_heat_insitu[,1] 

# TABLE COLD
tab_cold <- matrix(NA, nrow=10, ncol=5)
rownames(tab_cold) <- namesreg


tab_cold_insitu <- tab_cold_Haduk5 <- tab_cold_Haduk5_unw <- tab_cold_ERA5 <- tab_cold_ERA5_unw <- matrix(NA, nrow=length(fulldataset), ncol=3)
rownames(tab_cold_insitu) <- rownames(tab_cold_Haduk5) <- rownames(tab_cold_Haduk5_unw) <- rownames(tab_cold_ERA5) <- rownames(tab_cold_ERA5_unw) <- rep(names(fulldataset))



colnames(tab_heat) <- colnames(tab_cold) <- c("insitu", "haduk", "haduk_unw", "era5", "era5_unw")

# OBJECTS FOR THE ESTIMATIONS OF ATTRIBUTABLE FRACTIONS
# LOAD FUNCTION
source("attrdl.r")

# CREATE THE VECTORS TO STORE THE TOTAL MORTALITY (ACCOUNTING FOR MISSING)
totdeathall <- rep(NA,length(namesreg))
names(totdeathall) <- namesreg

# CREATE THE MATRIX TO STORE THE ATTRIBUTABLE DEATHS USING THE CITYBLUPS
matsim_insitu <- matsim_haduk <- matsim_haduk_unw <- matsim_era5  <- matsim_era5_unw <- matrix(NA,length(namesreg),3,
                                                                                               dimnames=list(namesreg,c("cold","heat","total")))

# NUMBER OF SIMULATION RUNS FOR COMPUTING EMPIRICAL CI
nsim <- 1000

# CREATE THE ARRAY TO STORE THE CI OF ATTRIBUTABLE DEATHS USING THE CITYBLUPS
arraysim_insitu <- arraysim_haduk <- arraysim_haduk_unw <- arraysim_era5 <- arraysim_era5_unw <- array(NA,dim=c(length(namesreg),3,nsim),dimnames=list(namesreg,
                                                                                                                                                       c("cold","heat","total")))


## create list to store MMTs
mmt_list <- matrix(NA, nrow=10, ncol=5)
# DEFINITION SEQ OF PERCENTILES
predper <- c(seq(0,1,0.1),2:98,seq(99,100,0.1))


# RUN LOOP ACROSS REGIONS
for(i in seq(length(fulldataset))) {
  
  # PRINT
  cat(i,"")
  
  # EXTRACT THE DATA FOR REGION i
  data <- fulldataset[[i]]
  
  #### - Insitu ####
  # DEFINE CROSS-BASIS
  argvar <- list(fun=varfun, degree=vardegree,
                 knots=quantile(data$Insitu,varper/100,na.rm=T), 
                 Bound=range(data$Insitu,na.rm=T))
  #Bound=range(tmeanper,na.rm=T))
  
  arglag <- list(fun="ns",knots=logknots(lag,lagnk))
  cb <- crossbasis(data$Insitu,lag=lag,argvar=argvar,arglag=arglag)
  
  # RUN THE MODEL AND OBTAIN PREDICTIONS
  model <- glm(formula,data,family=quasipoisson,na.action="na.exclude")
  
  # REDUCTION TO OVERALL CUMULATIVE
  pred <- crosspred(cb,model,cen=mean(data$Insitu,na.rm=T))
  plot(pred)
  
  redtemp <- crossreduce(cb,model,cen=mean(data$Insitu,na.rm=T))
  mmt <- redtemp$predvar[which.min(redtemp$fit)[[1]]]
  mmt_list[i,1] <- mmt
  redallinsitu[[i]] <- crossreduce(cb,model,cen=mmt)
  plot(redallinsitu[[1]])
  
  coefallinsitu[i,] <- coef(redallinsitu[[i]])
  vcovallinsitu[[i]] <- vcov(redallinsitu[[i]])
  qaicres[i,1] <- QAIC(model)
  
  # STORE RR
  cutheat <- quantile(data$Insitu,0.99, na.rm=T)
  cutcold <- quantile(data$Insitu,0.01, na.rm=T)
  cut <- crossreduce(cb,model, at=c(cutheat,cutcold), cen=mmt)
  tab_heat[i,1] <- paste0(formatC(cut$RRfit[2],2, format="f"), sep=" (", 
                          formatC(cut$RRlow[2],2, format="f"), sep=" - ",  
                          formatC(cut$RRhigh[2],2, format="f"), sep=")")
  
  tab_cold[i,1] <- paste0(formatC(cut$RRfit[1],2, format="f"), sep=" (", 
                          formatC(cut$RRlow[1],2, format="f"), sep=" - ",  
                          formatC(cut$RRhigh[1],2, format="f"), sep=")")
  
  tab_heat_insitu[i,1:3]<- cbind(cut$RRfit[2], cut$RRlow[2],cut$RRhigh[2])
  tab_cold_insitu[i,1:3] <- cbind(cut$RRfit[1], cut$RRlow[1],cut$RRhigh[1])
  
  
  # ESTIMATION OF THE ATTRIBUTABLE FRACTIONS
  # ESTIMATION OF THE ATTRBUTABLE FRACTIONS
  matsim_insitu[i,"heat"] <- attrdl(data$Insitu,cb,data$all,coef=coefallinsitu[i,],
                                    vcov=vcovallinsitu[[i]],type="an",dir="forw",
                                    range=c((quantile(data$Insitu, 0.75, na.rm=T)),(quantile(data$Insitu, 1, na.rm=T))), cen=mmt)
  
  matsim_insitu[i,"cold"] <- attrdl(data$Insitu,cb,data$all,coef=coefallinsitu[i,],
                                    vcov=vcovallinsitu[[i]],type="an",dir="forw",
                                    range=c((quantile(data$Insitu, 0.0, na.rm=T)),(quantile(data$Insitu, 0.25, na.rm=T))), cen=mmt)
  
  matsim_insitu[i,"total"] <- attrdl(data$Insitu,cb,data$all,coef=coefallinsitu[i,],
                                     vcov=vcovallinsitu[[i]],type="an",dir="forw", cen=mmt)
  
  # COMPUTE EMPIRICAL OCCURRENCES OF THE ATTRIBUTABLE DEATHS
  # USED TO DERIVE CONFIDENCE INTERVALS
  arraysim_insitu[i,"heat",] <- attrdl(data$Insitu,cb,data$all,coef=coefallinsitu[i,],
                                       vcov=vcovallinsitu[[i]],type="an",dir="forw", sim=T,nsim=nsim,
                                       range=c((quantile(data$Insitu, 0.75, na.rm=T)),(quantile(data$Insitu, 1, na.rm=T))), cen=mmt)
  
  arraysim_insitu[i,"cold",] <- attrdl(data$Insitu,cb,data$all,coef=coefallinsitu[i,],
                                       vcov=vcovallinsitu[[i]],type="an",dir="forw", sim=T,nsim=nsim,
                                       range=c((quantile(data$Insitu, 0.0, na.rm=T)),(quantile(data$Insitu, 0.25, na.rm=T))), cen=mmt)
  
  arraysim_insitu[i,"total",] <- attrdl(data$Insitu,cb,data$all,coef=coefallinsitu[i,],
                                        vcov=vcovallinsitu[[i]],type="an",dir="forw", sim=T, nsim=nsim,cen=mmt)
  
  #### - HADUK5 ####
  # DEFINE THE CROSSBASIs
  argvar <- list(fun=varfun, degree=vardegree,
                 knots=quantile(data$Haduk5,varper/100,na.rm=T), 
                 Bound=range(data$Haduk5,na.rm=T))
  #Bound=range(tmeanper,na.rm=T))
  
  arglag <- list(fun="ns",knots=logknots(lag,lagnk))
  cb <- crossbasis(data$Haduk5,lag=lag,argvar=argvar,arglag=arglag)
  
  # RUN THE MODEL AND OBTAIN PREDICTIONS
  model <- glm(formula,data,family=quasipoisson,na.action="na.exclude")
  
  # REDUCTION TO OVERALL CUMULATIVE
  pred <- crosspred(cb,model,cen=mean(data$Haduk5,na.rm=T))
  
  
  redtemp <- crossreduce(cb,model,cen=mean(data$Haduk5,na.rm=T))
  mmt <- redtemp$predvar[which.min(redtemp$fit)[[1]]]
  mmt_list[i,2] <- mmt
  redallhaduk5[[i]] <- crossreduce(cb,model,cen=mmt)
  
  coefallhaduk5[i,] <- coef(redallhaduk5[[i]])
  vcovallhaduk5[[i]] <- vcov(redallhaduk5[[i]])
  qaicres[i,2] <- QAIC(model)
  
  cutheatHaduk5 <- quantile(data$Haduk5,0.99, na.rm=T)
  cutcoldHaduk5 <- quantile(data$Haduk5,0.01, na.rm=T)
  cuthaduk5 <- crossreduce(cb,model, at=c(cutheatHaduk5,cutcoldHaduk5), cen=mmt)
  
  
  tab_heat[i,2] <- paste0(formatC(cuthaduk5$RRfit[2],2, format="f"), sep=" (", 
                          formatC(cuthaduk5$RRlow[2],2, format="f"), sep=" - ",  
                          formatC(cuthaduk5$RRhigh[2],2, format="f"), sep=")")
  
  tab_cold[i,2] <- paste0(formatC(cuthaduk5$RRfit[1],2, format="f"), sep=" (", 
                          formatC(cuthaduk5$RRlow[1],2, format="f"), sep=" - ",  
                          formatC(cuthaduk5$RRhigh[1],2, format="f"), sep=")")
  
  tab_heat_Haduk5[i,1:3] <- cbind(cuthaduk5$RRfit[2], cuthaduk5$RRlow[2],cuthaduk5$RRhigh[2])
  tab_cold_Haduk5[i,1:3] <- cbind(cuthaduk5$RRfit[1], cuthaduk5$RRlow[1],cuthaduk5$RRhigh[1])
  
  # ESTIMATION OF THE ATTRIBUTABLE FRACTIONS
  matsim_haduk[i,"heat"] <- attrdl(data$Haduk5,cb,data$all,coef=coefallhaduk5[i,],
                                   vcov=vcovallhaduk5[[i]],type="an",dir="forw",
                                   range=c((quantile(data$Haduk5, 0.75, na.rm=T)),(quantile(data$Haduk5, 1, na.rm=T))), cen=mmt)
  
  
  matsim_haduk[i,"cold"] <- attrdl(data$Haduk5,cb,data$all,coef=coefallhaduk5[i,],
                                   vcov=vcovallhaduk5[[i]],type="an",dir="forw",
                                   range=c((quantile(data$Haduk5, 0.0, na.rm=T)),(quantile(data$Haduk5, 0.25, na.rm=T))), cen=mmt)
  
  
  matsim_haduk[i,"total"] <- attrdl(data$Haduk5,cb,data$all,coef=coefallhaduk5[i,],
                                    vcov=vcovallhaduk5[[i]],type="an",dir="forw", cen=mmt)
  
  # COMPUTE EMPIRICAL OCCURRENCES OF THE ATTRIBUTABLE DEATHS
  # USED TO DERIVE CONFIDENCE INTERVALS
  arraysim_haduk[i,"heat",] <- attrdl(data$Haduk5,cb,data$all,coef=coefallhaduk5[i,],
                                      vcov=vcovallhaduk5[[i]],type="an",dir="forw", sim=T,nsim=nsim,
                                      range=c((quantile(data$Haduk5, 0.75, na.rm=T)),(quantile(data$Haduk5, 1, na.rm=T))), cen=mmt)
  
  
  arraysim_haduk[i,"cold",] <- attrdl(data$Haduk5,cb,data$all,coef=coefallhaduk5[i,],
                                      vcov=vcovallhaduk5[[i]],type="an",dir="forw", sim=T,nsim=nsim,
                                      range=c((quantile(data$Haduk5, 0.0, na.rm=T)),(quantile(data$Haduk5, 0.25, na.rm=T))), cen=mmt)

  
  arraysim_haduk[i,"total",] <- attrdl(data$Haduk5,cb,data$all,coef=coefallhaduk5[i,],
                                       vcov=vcovallhaduk5[[i]],type="an",dir="forw", sim=T, nsim=nsim,cen=mmt)
  
  
  #### - UNWEIGHTED HADUK5####
  # DEFINE THE CROSSBASIs
  argvar <- list(fun=varfun, degree=vardegree,
                 knots=quantile(data$haduk5_unw,varper/100,na.rm=T), 
                 Bound=range(data$haduk5_unw,na.rm=T))
  #Bound=range(tmeanper,na.rm=T))
  
  arglag <- list(fun="ns",knots=logknots(lag,lagnk))
  cb <- crossbasis(data$haduk5_unw,lag=lag,argvar=argvar,arglag=arglag)
  
  # RUN THE MODEL AND OBTAIN PREDICTIONS
  model <- glm(formula,data,family=quasipoisson,na.action="na.exclude")
  
  # REDUCTION TO OVERALL CUMULATIVE
  pred <- crosspred(cb,model,cen=mean(data$haduk5_unw,na.rm=T))
  redtemp <- crossreduce(cb,model,cen=mean(data$haduk5_unw,na.rm=T))
  mmt <- redtemp$predvar[which.min(redtemp$fit)[[1]]]
  mmt_list[i,3] <- mmt
  redallhaduk5_unw[[i]] <- crossreduce(cb,model,cen=mmt)
  
  coefallhaduk5_unw[i,] <- coef(redallhaduk5_unw[[i]])
  vcovallhaduk5_unw[[i]] <- vcov(redallhaduk5_unw[[i]])
  qaicres[i,3] <- QAIC(model)
  
  cutheatHaduk5_unw <- quantile(data$haduk5_unw,0.99, na.rm=T)
  cutcoldHaduk5_unw <- quantile(data$haduk5_unw,0.01, na.rm=T)
  cuthaduk5_unw <- crossreduce(cb,model, at=c(cutheatHaduk5_unw,cutcoldHaduk5_unw), cen=mmt)
  
  
  tab_heat[i,3] <- paste0(formatC(cuthaduk5_unw$RRfit[2],2, format="f"), sep=" (", 
                          formatC(cuthaduk5_unw$RRlow[2],2, format="f"), sep=" - ",  
                          formatC(cuthaduk5_unw$RRhigh[2],2, format="f"), sep=")")
  
  tab_cold[i,3] <- paste0(formatC(cuthaduk5_unw$RRfit[1],2, format="f"), sep=" (", 
                          formatC(cuthaduk5_unw$RRlow[1],2, format="f"), sep=" - ",  
                          formatC(cuthaduk5_unw$RRhigh[1],2, format="f"), sep=")")
  
  tab_heat_Haduk5_unw[i,1:3] <- cbind(cuthaduk5_unw$RRfit[2], cuthaduk5_unw$RRlow[2],cuthaduk5_unw$RRhigh[2])
  tab_cold_Haduk5_unw[i,1:3] <- cbind(cuthaduk5_unw$RRfit[1], cuthaduk5_unw$RRlow[1],cuthaduk5_unw$RRhigh[1])
  
  # ESTIMATION OF THE ATTRIBUTABLE FRACTIONS
  matsim_haduk_unw[i,"heat"] <- attrdl(data$haduk5_unw,cb,data$all,coef=coefallhaduk5_unw[i,],
                                       vcov=vcovallhaduk5_unw[[i]],type="an",dir="forw",
                                       range=c((quantile(data$haduk5_unw, 0.75, na.rm=T)),(quantile(data$haduk5_unw, 1, na.rm=T))), cen=mmt)
  
  matsim_haduk_unw[i,"cold"] <- attrdl(data$haduk5_unw,cb,data$all,coef=coefallhaduk5_unw[i,],
                                       vcov=vcovallhaduk5_unw[[i]],type="an",dir="forw",
                                       range=c((quantile(data$haduk5_unw, 0.0, na.rm=T)),(quantile(data$haduk5_unw, 0.25, na.rm=T))), cen=mmt)
  
  matsim_haduk_unw[i,"total"] <- attrdl(data$haduk5_unw,cb,data$all,coef=coefallhaduk5_unw[i,],
                                        vcov=vcovallhaduk5_unw[[i]],type="an",dir="forw", cen=mmt)
  
  # COMPUTE EMPIRICAL OCCURRENCES OF THE ATTRIBUTABLE DEATHS
  # USED TO DERIVE CONFIDENCE INTERVALS
  arraysim_haduk_unw[i,"heat",] <- attrdl(data$haduk5_unw,cb,data$all,coef=coefallhaduk5_unw[i,],
                                          vcov=vcovallhaduk5_unw[[i]],type="an",dir="forw", sim=T,nsim=nsim,
                                          range=c((quantile(data$haduk5_unw, 0.75, na.rm=T)),(quantile(data$haduk5_unw, 1, na.rm=T))), cen=mmt)
  
  arraysim_haduk_unw[i,"cold",] <- attrdl(data$haduk5_unw,cb,data$all,coef=coefallhaduk5_unw[i,],
                                          vcov=vcovallhaduk5_unw[[i]],type="an",dir="forw", sim=T,nsim=nsim,
                                          range=c((quantile(data$haduk5_unw, 0.0, na.rm=T)),(quantile(data$haduk5_unw, 0.25, na.rm=T))), cen=mmt)
  
  arraysim_haduk_unw[i,"total",] <- attrdl(data$haduk5_unw,cb,data$all,coef=coefallhaduk5_unw[i,],
                                           vcov=vcovallhaduk5_unw[[i]],type="an",dir="forw", sim=T, nsim=nsim,cen=mmt)
  
  
  #### - ERA5 ####
  # DEFINE THE CROSSBASIS
  argvar <- list(fun=varfun, degree=vardegree,
                 knots=quantile(data$ERA5,varper/100,na.rm=T), 
                 Bound=range(data$ERA5,na.rm=T))
  #Bound=range(tmeanper,na.rm=T))
  
  arglag <- list(fun="ns",knots=logknots(lag,lagnk))
  cb <- crossbasis(data$ERA5,lag=lag,argvar=argvar,arglag=arglag)
  
  # RUN THE MODEL AND OBTAIN PREDICTIONS
  model <- glm(formula,data,family=quasipoisson,na.action="na.exclude")
  
  # REDUCTION TO OVERALL CUMULATIVE
  pred <- crosspred(cb,model,cen=mean(data$ERA5,na.rm=T))
  redtemp <- crossreduce(cb,model,cen=mean(data$ERA5,na.rm=T))
  mmt <- redtemp$predvar[which.min(redtemp$fit)[[1]]]
  mmt_list[i,4] <- mmt
  redallERA5[[i]] <- crossreduce(cb,model,cen=mmt)
  
  coefallERA5[i,] <- coef(redallERA5[[i]])
  vcovallERA5[[i]] <- vcov(redallERA5[[i]])
  qaicres[i,4] <- QAIC(model)
  
  cutheatERA5 <- quantile(data$ERA5,0.99, na.rm=T)
  cutcoldERA5 <- quantile(data$ERA5,0.01, na.rm=T)
  cutERA5 <- crossreduce(cb,model, at=c(cutheatERA5,cutcoldERA5), cen=mmt)
  
  tab_heat[i,4] <- paste0(formatC(cutERA5$RRfit[2],2, format="f"), sep=" (", 
                          formatC(cutERA5$RRlow[2],2, format="f"), sep=" - ",  
                          formatC(cutERA5$RRhigh[2],2, format="f"), sep=")")
  
  tab_cold[i,4] <- paste0(formatC(cutERA5$RRfit[1],2, format="f"), sep=" (", 
                          formatC(cutERA5$RRlow[1],2, format="f"), sep=" - ",  
                          formatC(cutERA5$RRhigh[1],2, format="f"), sep=")")
  
  tab_heat_ERA5[i,1:3] <- cbind(cutERA5$RRfit[2], cutERA5$RRlow[2],cutERA5$RRhigh[2])
  tab_cold_ERA5[i,1:3] <- cbind(cutERA5$RRfit[1], cutERA5$RRlow[1],cutERA5$RRhigh[1])
  
  # ESTIMATION OF THE ATTRIBUTABLE FRACTIONS
  matsim_era5[i,"heat"] <- attrdl(data$ERA5,cb,data$all,coef=coefallERA5[i,],
                                  vcov=vcovallERA5[[i]],type="an",dir="forw",
                                  range=c((quantile(data$ERA5, 0.75, na.rm=T)),(quantile(data$ERA5, 1, na.rm=T))), cen=mmt)
  
  matsim_era5[i,"cold"] <- attrdl(data$ERA5,cb,data$all,coef=coefallERA5[i,],
                                  vcov=vcovallERA5[[i]],type="an",dir="forw",
                                  range=c((quantile(data$ERA5, 0.0, na.rm=T)),(quantile(data$ERA5, 0.25, na.rm=T))), cen=mmt)
  
  matsim_era5[i,"total"] <- attrdl(data$ERA5,cb,data$all,coef=coefallERA5[i,],
                                   vcov=vcovallERA5[[i]],type="an",dir="forw", cen=mmt)
  
  # COMPUTE EMPIRICAL OCCURRENCES OF THE ATTRIBUTABLE DEATHS
  # USED TO DERIVE CONFIDENCE INTERVALS
  arraysim_era5[i,"heat",] <- attrdl(data$ERA5,cb,data$all,coef=coefallERA5[i,],
                                     vcov=vcovallERA5[[i]],type="an",dir="forw", sim=T,nsim=nsim,
                                     range=c((quantile(data$ERA5, 0.75, na.rm=T)),(quantile(data$ERA5, 1, na.rm=T))), cen=mmt)
  
  arraysim_era5[i,"cold",] <- attrdl(data$ERA5,cb,data$all,coef=coefallERA5[i,],
                                     vcov=vcovallERA5[[i]],type="an",dir="forw", sim=T,nsim=nsim,
                                     range=c((quantile(data$ERA5, 0.0, na.rm=T)),(quantile(data$ERA5, 0.25, na.rm=T))), cen=mmt)
  
  arraysim_era5[i,"total",] <- attrdl(data$ERA5,cb,data$all,coef=coefallERA5[i,],
                                      vcov=vcovallERA5[[i]],type="an",dir="forw", sim=T, nsim=nsim,cen=mmt)
  
  
  
  #### - UNWEIGHTED ERA5 ####
  # DEFINE THE CROSSBASIS
  argvar <- list(fun=varfun, degree=vardegree,
                 knots=quantile(data$ERA5_unw,varper/100,na.rm=T), 
                 Bound=range(data$ERA5_unw,na.rm=T))
  
  
  arglag <- list(fun="ns",knots=logknots(lag,lagnk))
  cb <- crossbasis(data$ERA5_unw,lag=lag,argvar=argvar,arglag=arglag)
  
  # RUN THE MODEL AND OBTAIN PREDICTIONS
  model <- glm(formula,data,family=quasipoisson,na.action="na.exclude")
  
  # REDUCTION TO OVERALL CUMULATIVE
  pred <- crosspred(cb,model,cen=mean(data$ERA5_unw,na.rm=T))
  redtemp <- crossreduce(cb,model,cen=mean(data$ERA5_unw,na.rm=T))
  mmt <- redtemp$predvar[which.min(redtemp$fit)[[1]]]
  mmt_list[i,5] <- mmt
  
  redallERA5_unw[[i]] <- crossreduce(cb,model,cen=mmt)
  
  coefallERA5_unw[i,] <- coef(redallERA5_unw[[i]])
  vcovallERA5_unw[[i]] <- vcov(redallERA5_unw[[i]])
  qaicres[i,5] <- QAIC(model)
  
  cutheatERA5_unw <- quantile(data$ERA5_unw,0.99, na.rm=T)
  cutcoldERA5_unw <- quantile(data$ERA5_unw,0.01, na.rm=T)
  cutERA5_unw <- crossreduce(cb,model, at=c(cutheatERA5_unw,cutcoldERA5_unw), cen=mmt)
  
  
  tab_heat[i,5] <- paste0(formatC(cutERA5_unw$RRfit[2],2, format="f"), sep=" (", 
                          formatC(cutERA5_unw$RRlow[2],2, format="f"), sep=" - ",  
                          formatC(cutERA5_unw$RRhigh[2],2, format="f"), sep=")")
  
  tab_cold[i,5] <- paste0(formatC(cutERA5_unw$RRfit[1],2, format="f"), sep=" (", 
                          formatC(cutERA5_unw$RRlow[1],2, format="f"), sep=" - ",  
                          formatC(cutERA5_unw$RRhigh[1],2, format="f"), sep=")")
  
  tab_heat_ERA5_unw[i,1:3] <- cbind(cutERA5_unw$RRfit[2], cutERA5_unw$RRlow[2],cutERA5_unw$RRhigh[2])
  tab_cold_ERA5_unw[i,1:3] <- cbind(cutERA5_unw$RRfit[1], cutERA5_unw$RRlow[1],cutERA5_unw$RRhigh[1])
  
  
  # ESTIMATION OF THE ATTRIBUTABLE FRACTIONS
  matsim_era5_unw[i,"heat"] <- attrdl(data$ERA5_unw,cb,data$all,coef=coefallERA5_unw[i,],
                                      vcov=vcovallERA5_unw[[i]],type="an",dir="forw",
                                      range=c((quantile(data$ERA5_unw, 0.75, na.rm=T)),(quantile(data$ERA5_unw, 1, na.rm=T))), cen=mmt)
  
  matsim_era5_unw[i,"cold"] <- attrdl(data$ERA5_unw,cb,data$all,coef=coefallERA5_unw[i,],
                                      vcov=vcovallERA5_unw[[i]],type="an",dir="forw",
                                      range=c((quantile(data$ERA5_unw, 0.0, na.rm=T)),(quantile(data$ERA5_unw, 0.25, na.rm=T))), cen=mmt)
  
  matsim_era5_unw[i,"total"] <- attrdl(data$ERA5_unw,cb,data$all,coef=coefallERA5_unw[i,],
                                       vcov=vcovallERA5_unw[[i]],type="an",dir="forw", cen=mmt)
  
  # COMPUTE EMPIRICAL OCCURRENCES OF THE ATTRIBUTABLE DEATHS
  # USED TO DERIVE CONFIDENCE INTERVALS
  arraysim_era5_unw[i,"heat",] <- attrdl(data$ERA5_unw,cb,data$all,coef=coefallERA5_unw[i,],
                                         vcov=vcovallERA5_unw[[i]],type="an",dir="forw", sim=T,nsim=nsim,
                                         range=c((quantile(data$ERA5_unw, 0.75, na.rm=T)),(quantile(data$ERA5_unw, 1, na.rm=T))), cen=mmt)
  
  arraysim_era5_unw[i,"cold",] <- attrdl(data$ERA5_unw,cb,data$all,coef=coefallERA5_unw[i,],
                                         vcov=vcovallERA5_unw[[i]],type="an",dir="forw", sim=T,nsim=nsim,
                                         range=c((quantile(data$ERA5_unw, 0.0, na.rm=T)),(quantile(data$ERA5_unw, 0.25, na.rm=T))), cen=mmt)
  
  arraysim_era5_unw[i,"total",] <- attrdl(data$ERA5_unw,cb,data$all,coef=coefallERA5_unw[i,],
                                          vcov=vcovallERA5_unw[[i]],type="an",dir="forw", sim=T, nsim=nsim,cen=mmt)
  
  totdeathall[i] <- sum(data$all,na.rm=T)
  
}




an_insitu <- apply(arraysim_insitu,c(1,2),quantile,0.5, na.rm=T)
an_insitulow <- apply(arraysim_insitu,c(1,2),quantile,0.025, na.rm=T)
an_insituhigh <- apply(arraysim_insitu,c(1,2),quantile,0.975, na.rm=T)
rownames(an_insitu) <- rownames(an_insitulow) <- rownames(an_insituhigh) <- namesreg

an_haduk <-  apply(arraysim_haduk,c(1,2),quantile,0.5, na.rm=T)
an_haduklow <- apply(arraysim_haduk,c(1,2),quantile,0.025, na.rm=T)
an_hadukhigh <- apply(arraysim_haduk,c(1,2),quantile,0.975, na.rm=T)
rownames(an_haduk) <- rownames(an_haduklow) <- rownames(an_hadukhigh) <- namesreg

an_Haduk_unw <- apply(arraysim_haduk_unw,c(1,2),quantile,0.5, na.rm=T)
an_Haduklow_unw <- apply(arraysim_haduk_unw,c(1,2),quantile,0.025, na.rm=T)
an_Hadukhigh_unw <- apply(arraysim_haduk_unw,c(1,2),quantile,0.975, na.rm=T)
rownames(an_Haduk_unw) <- rownames(an_Haduklow_unw) <- rownames(an_Hadukhigh_unw) <- namesreg

an_era5 <- apply(arraysim_era5,c(1,2),quantile,0.5, na.rm=T)
an_era5low <- apply(arraysim_era5,c(1,2),quantile,0.025, na.rm=T)
an_era5high <- apply(arraysim_era5,c(1,2),quantile,0.975, na.rm=T)
rownames(an_era5) <- rownames(an_era5low) <- rownames(an_era5high) <- namesreg

an_era5_unw <- apply(arraysim_era5_unw,c(1,2),quantile,0.5, na.rm=T)
an_era5low_unw <- apply(arraysim_era5_unw,c(1,2),quantile,0.025, na.rm=T)
an_era5high_unw <- apply(arraysim_era5_unw,c(1,2),quantile,0.975, na.rm=T)
rownames(an_era5_unw) <- rownames(an_era5low_unw) <- rownames(an_era5high_unw) <- namesreg


# ATTRIBUTABLE FRACTIONS
af_insitu <- an_insitu/totdeathall*100
af_insitulow <- an_insitulow/totdeathall*100
af_insituhigh <- an_insituhigh/totdeathall*100

af_Haduk5 <- an_haduk/totdeathall*100
af_Haduk5low <- an_haduklow/totdeathall*100
af_Haduk5high <- an_hadukhigh/totdeathall*100

af_Haduk5_unw <- an_Haduk_unw/totdeathall*100
af_Haduk5low_unw <- an_Haduklow_unw/totdeathall*100
af_Hadukhigh5_unw <- an_Hadukhigh_unw/totdeathall*100

af_ERA5 <- an_era5/totdeathall*100
af_ERA5low <- an_era5low/totdeathall*100
af_ERA5high <- an_era5high/totdeathall*100

af_ERA5_unw <- an_era5_unw/totdeathall*100
af_ERA5low_unw <- an_era5low_unw/totdeathall*100
af_ERA5high_unw <- an_era5high_unw/totdeathall*100










