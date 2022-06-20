########################################################
##### 2ND STAGE TEMP TRENDS ############################
########################################################


########################
### 1st stage coef #####
########################


# SET INPUT LOCATION FILE
load("~/Desktop/PhD/publication1/final code/clean/02stage_temp_trends.Rdata")
dir_stage01 <- "~/Desktop/PhD/publication1/final code/clean"
load(paste0(dir_stage01,"/01stagedata_allyecb.Rdata"))

# APPEND THE RESULTS
coefall_age <- rbind(coefall1,coefall2,coefall3)
vcovall_age <- c(vcovall1, vcovall2, vcovall3)

colnames(coefall_age) <- c("canton", "indtime", "age_cat", paste0("coeff", c(1:5)))



########################
### variables #########
#######################
#inpath <- "Y:/CCH/Evan/temporal_trends_paper/data"
#indregion <- read.csv(paste0(inpath,"/regions.csv"), sep=";")
#indregion <- indregion[,2]

#popdens <- readRDS(paste0(inpath,"/densitylong.rds"))
#proportion_age <- readRDS(paste0(inpath,"/proportion_age.RDS"))
#proportion_65 <- as.matrix(proportion_age[,3])
#proportion_80 <- as.matrix(proportion_age[,4])



##########################
### mixmeta  ##########
##########################
mvall1 <- mixmeta(coefall_age[,-c(1:3)] ~ rep(avertmean_sp,3),vcovall_age,control=list(showiter=T),method="ml")

mvall1 <- mixmeta(coefall_age[,-c(1:3)] ~ rep(avertmean_sp,3) + rep(rangetmean_sp,3),vcovall_age,control=list(showiter=T),method="ml")


mvall1 <- mixmeta(coefall_age[,-c(1:3)] ~ rep(avertmean_sp,3) + rep(rangetmean_sp,3),rep(indregion,3),vcovall_age,control=list(showiter=T),method="ml")

mvall1 <- mixmeta(coefall_age[,-c(1:3)] ~ rep(avertmean_sp,3) + rep(rangetmean_sp,3), random = ~1 | rep(indregion,3),vcovall_age,control=list(showiter=T),method="ml")


mvall1 <- mixmeta(coefall_age[,-c(1:3)] ~ rep(avertmean_sp,3) + rep(rangetmean_sp,3) +
                    as.factor(coefall_age[,"indtime"]) +
                    as.factor(coefall_age[,"age_cat"]),vcovall_age,
                  random = ~1 | rep(indregion,3),
                  control=list(showiter=T),method="ml")


mvall2 <- mixmeta(coefall_age[,-c(1:3)] ~ rep(avertmean_sp,3) + 
                   as.factor(coefall_age[,"indtime"]) +
                   as.factor(coefall_age[,"age_cat"]),vcovall_age,
                 random = ~1 | rep(indregion,3),
                 control=list(showiter=T),method="ml")

# WALD TEST
fwald <- function(model,var) {
  ind <- grep(var,names(coef(model)))
  coef <- coef(model)[ind]
  vcov <- vcov(model)[ind,ind]
  waldstat <- coef%*%solve(vcov)%*%coef
  df <- length(coef)
  return(1-pchisq(waldstat,df))
}


summary(mvall)
summary(mvall2)

fwald(mvall,"indtime")
fwald(mvall,"age_cat")
fwald(mvall,"avertmean_sp")


################
###### blup ####
################

blupall.age <- blup(mvall,vcov=T)
blup <- blup(mvall,vcov=T)

################################################################################
# RE-CENTERING


# DEFINE MINIMUM MORTALITY VALUES: EXCLUDE LOW AND VERY HOT TEMPERATURE
for(i in seq(length(dlist))) {
  data <- dlist[[i]]
  predvar <- quantile(data$tmean,1:99/100,na.rm=T)
  # REDEFINE THE FUNCTION USING ALL THE ARGUMENTS (BOUNDARY KNOTS INCLUDED)
  argvar <- list(x=predvar,fun=varfun,
                 knots=quantile(data$tmean,varper/100,na.rm=T),degree=vardegree,
                 Bound=range(data$tmean,na.rm=T))
  bvar <- do.call(onebasis,argvar)
  minperccity[i] <- (1:99)[which.min((bvar%*%blup[[i]]$blup))]
  mintempcity[i] <- quantile(data$tmean,minperccity[i]/100,na.rm=T)
  
  ### ERC PLOTS
  pred <- crosspred(bvar,coef=blup[[i]]$blup,vcov=blup[[i]]$vcov,
                    model.link="log",by=0.1,cen=mintempdistrict[i])
  plot(pred,type="n",ylim=c(0,2.0),yaxt="n",lab=c(6,5,7),xlab="Temp",ylab="RR", main=names(dlist_total)[[i]])
}




rm(list=ls()[! ls() %in% c("mvall", "blupall.age", "coefall_age", "vcovall_age", "avertmean_sp", "indregion")])
load("~/Desktop/PhD/publication1/final code/clean/02stage_temp_trends.Rdata")








