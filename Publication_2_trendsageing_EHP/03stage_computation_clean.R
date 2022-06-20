

# LOAD LIBRARIES
library(mixmeta)
install.packages(c("lubridate", "dlnm"))
library(lubridate)
library(dlnm)
library(MASS)
library(writexl)

# LOAD OBJECTS
#dir <- "Y:/CCH/Evan/temporal_trends_paper/clean"
dir <-  "~/Desktop/PhD/publication1/final code/clean"
# LOAD DATA
munciptable <- readRDS(paste0(dir,"/municipalitylookuptable.rds"))
deathrecords <- readRDS(paste0(dir,"/deathrecordsmuncipality.rds"))


load(paste0(dir,"/agedeaths_19692017.Rdata"))
load(paste0(dir,"/02stage_temp_trends.Rdata"))



# SPECIFICATION OF THE EXPOSURE FUNCTION
varfun <- "bs"
vardegree <- 2
varper <- c(10, 75,90)

# SPECIFICATION OF THE LAG FUNCTION
lag <- 21
lagnk <- 3

# NUMBER OF SIMULATIONS
nsim <- 1000

# DEFINE PERIODS
tabperiod <- matrix(NA, nrow=5, ncol=2)
tabperiod[,1] <- seq(1969,2017, by=10)
tabperiod[,2] <- tabperiod[,1] + 9
tabperiod[5,2] <- 2017  ## LAST YEAR IS 2017

# CREATE MATRIX FOR ATTRIBUTABLE MORTALITY
totdeathall <- array(NA,dim=c(length(dlist)+1,6,4), dimnames=list(c(paste0("cant",seq(1:26)),"total"),
                                      c("1969/78", "1979/88","1989/99","1999/08","2009/2017","overall"),
                                      c("age1","age2","age3","allages")))

# CREATE ARRAY TO STORE SIMULATIONS
arraysim <- array(NA,dim=c(length(dlist),5,3,3,nsim+1),
                        dimnames=list(c(paste0("cant",seq(1:26))),
                                      c("1969/78", "1979/88","1989/99","1999/08","2009/2017"),
                                      c("age1","age2","age3"),  
                                      c("cold", "heat","total"), 
                                      c("est",paste0("sim",seq(nsim)))))

minperccant <- array(NA,dim=c(length(dlist),5,3),
                        dimnames=list(paste0("cant",seq(1:26)),
                                      c("1969/78", "1979/88","1989/99","1999/08","2009/2017"),
                                      c("age1","age2","age3")))

mintempcant <- array(NA,dim=c(length(dlist),5,3),
                        dimnames=list(paste0("cant",seq(1:26)),
                                      c("1969/78", "1979/88","1989/99","1999/08","2009/2017"),
                                      c("age1","age2","age3")))
# LOAD FUNCTION ATTRDL
#source("attrdl.r")

# LOOP 1 - ACROSS CANTONS
for(i in seq(length(dlist))) {
  
  cat("",i)
  
  data <- dlist[[i]]
  data$yy <- year(data$date)
  
  varknots <- quantile(data$aver.temp,varper/100, na.rm=T)
  lagknots <- logknots(lag, 3)
  boundknt <- range(data$aver.temp, na.rm=TRUE)
  
  # LOOP 2 - ACCROSS PERIODS
  for (j in seq(nrow(tabperiod))) {
    
    # SELECT SUBPERIOD
    data.sub <- subset(data, yy %in% c(as.numeric(tabperiod[j,1]):as.numeric(tabperiod[j,2])))
    
    # LOOP 3 - ACROSS AGE CATEGORIES
    for (k in seq(length(unique(coefall_age[,"age_cat"])))){

      # EXTRACT PARAMETERS OF THE AGE-SPECIFIC CATEGORY
      ind <- coefall_age[,"age_cat"]==k & coefall_age[,"indtime"]==j & coefall_age[,"canton"]==i
      blupsel <- blupall.age[ind]
      coef <- blupsel[[1]]$blup
      vcov <- blupsel[[1]]$vcov
      
      # ESTIMATE MMT
      predvar <- quantile(data.sub$aver.temp,75:98/100,na.rm=T)
      argvar <- list(x=predvar,fun=varfun,
                     knots=varknots,
                     Bound=boundknt)
      if(!is.null (vardegree)) argvar$degree <- vardegree
      bvar <- do.call(onebasis,argvar)  

      minperccant[i,j,k] <- (75:98)[which.min(bvar%*%coef)]
      mintempcant[i,j,k] <- quantile(data.sub$aver.temp,minperccant[i,j,k]/100,na.rm=T)
      
      # argvar <- list(fun=varfun,degree=vardegree, knots=varknots, ######################
      #            Bound=boundknt)
      # arglag <- list(knots=logknots(lag,lagnk)) ######################
      # cb <- crossbasis(data.sub$aver.temp,argvar=argvar,arglag=arglag) ######################
      
      # SELECT CENTERING POINT (ACCORDING TO i,j,k)
      cen <- mintempcant[i,j,k]

      # RE-COMPUTE THE  ARGVAR
      argvar <- list(fun=varfun,degree=vardegree,
                     knots=quantile(data$aver.temp,varper/100,na.rm=T),
                     Bound=range(data$aver.temp,na.rm=T))

      # DERIVE THE CENTERED BASIS
      bvar <- do.call(onebasis,c(list(x=data.sub[,"aver.temp"]),argvar))
      cenvec <- do.call(onebasis,c(list(x=cen),argvar))
      bvarcen <- scale(bvar,center=cenvec,scale=F)
  
      # INDICATOR FOR COLD/HEAT DAYS
      indheat <- data.sub[,"aver.temp"]>cen
      indcold <- data.sub[,"aver.temp"]<cen
    
      # DEFINE OUTCOME VARIABLE
      out <- data.sub[,3+k]
      
      #COMPUTE THE DAILY CONTRIBUTIONS OF ATTRIBUTABLE DEATHS
      an <- (1-exp(-bvarcen%*%coef))*out

      # POINT ESTIMATE
      arraysim[i,j,k,"heat",1] <- sum(an[indheat])
      arraysim[i,j,k,"cold",1] <- sum(an[indcold])
      arraysim[i,j,k,"total",1] <- sum(an)
    
      
      # arraysim[i,j,k,"heat",1] <- attrdl(data.sub$aver.temp,cb,out,coef=coef,
      #                                    vcov=vcov,type="an",dir="forw",
      #                                    range=c(cen,100), cen=cen)
      # 
      # arraysim[i,j,k,"cold",1] <- attrdl(data.sub$aver.temp,cb,out,coef=coef,
      #                                    vcov=vcov,type="an",dir="forw",
      #                                    range=c(-100,cen), cen=cen)
      # 
      # arraysim[i,j,k,"total",1] <- attrdl(data.sub$aver.temp,cb,out,coef=coef,
      #                                    vcov=vcov,type="an",dir="forw",
      #                                     range=c(-100,100), cen=cen)
      
      # DEFINE DENOMINATOR
      totdeathall[i,j,k] <- sum(out,na.rm=T) 

      # arraysim[i,j,k,"heat",-1] <- attrdl(data.sub$aver.temp,cb,out,coef=coef,
      #                                    vcov=vcov,type="an",dir="forw", sim=T,nsim=nsim,
      #                                       range=c(cen,100), cen=cen)
      # 
      # arraysim[i,j,k,"cold",-1] <- attrdl(data.sub$aver.temp,cb,out,coef=coef,
      #                                    vcov=vcov,type="an",dir="forw", sim=T,nsim=nsim,
      #                                       range=c(-100,cen), cen=cen)
      # 
      # arraysim[i,j,k,"total",-1] <- attrdl(data.sub$aver.temp,cb,out,coef=coef,
      #                                    vcov=vcov,type="an",dir="forw", sim=T, nsim=nsim,
      #                                 cen=cen)
      
      
      # SAMPLE COEF ASSUMING A MULTIVARIATE NORMAL DISTRIBUTION
      set.seed(13041975+j+k)
      coefsim <- mvrnorm(nsim,coef,vcov)

      # LOOP 4 - ACROSS ITERATIONS
      for(s in seq(nsim)) {

          # COMPUTE THE DAILY CONTRIBUTIONS OF ATTRIBUTABLE DEATHS
          an <- (1-exp(-bvarcen%*%coefsim[s,]))*out

          arraysim[i,j,k,"heat",1+s] <- sum(an[indheat])
          arraysim[i,j,k,"cold",1+s] <- sum(an[indcold])
          arraysim[i,j,k,"total",1+s] <- sum(an)
      }
    }
  }
}
install.packages("reshape2")
library(reshape2)
load(paste0(dir,"/res_attr_tempvar_final.RData"))

dir <-  "~/Desktop/PhD/publication1/final code/clean"
setwd("/Users/evandeschrijver/Desktop/PhD/publication1/final code/clean")
rate_pop_65 <- read.csv("rates_65.csv", sep=";")
rate_pop_65$overall <-  rowSums(rate_pop_65[,2:6])

rate_pop_6579<- read.csv("rates_6579.csv", sep=";")
rate_pop_6579$overall <-  rowSums(rate_pop_6579[,2:6])

rate_pop_80 <- read.csv("rates_80.csv", sep=";")
rate_pop_80$overall <-  rowSums(rate_pop_80[,2:6])


totpoprate <- array(NA,dim=c(length(dlist)+1,6,4), dimnames=list(c(paste0("cant",seq(1:26)),"total"),
                                                                  c("1969/78", "1979/88","1989/98","1999/08","2009/2017","mean"),
                                                                  c("age1","age2","age3","allages")))


totpoprate[,,1] <- data.matrix(rate_pop_65[,2:7])
totpoprate[,,2] <- data.matrix(rate_pop_6579[,2:7])
totpoprate[,,3] <- data.matrix(rate_pop_80[,2:7])
totpoprate[,,4]  <- totpoprate[,,1] + totpoprate[,,2] + totpoprate[,,3] 

# OBTAIN SUMMARIES

std_rate <- ratecanton <- ratecanton5thdec <- ancanton <- afcanton <- ratecanton5th <- array(NA,dim=c(length(dlist)+1,6,4,3,3),
                        dimnames=list(c(paste0("cant",seq(1:26)),"total"),
                                      c("1969/78", "1979/88","1989/98","1999/08","2009/2017","overall"),
                                      c("age1","age2","age3","allages"),  
                                      c("cold", "heat","total"), 
                                      c("est","ci.l","ci.h")))

totdeathall[,6,1] <- rowSums(totdeathall[,-6,1])
totdeathall[,6,2] <- rowSums(totdeathall[,-6,2])
totdeathall[,6,3] <- rowSums(totdeathall[,-6,3])
totdeathall[,,4] <- rowSums(totdeathall[,,-4])
totdeathall[27,,1] <- colSums(totdeathall[-27,,1])
totdeathall[27,,2] <- colSums(totdeathall[-27,,2])
totdeathall[27,,3] <- colSums(totdeathall[-27,,3])
totdeathall[27,,4] <- colSums(totdeathall[-27,,4])
totdeathall[,,4] <- totdeathall[,,1] + totdeathall[,,2] + totdeathall[,,3] 

totpoprate[,,4] <- totpoprate[,,1] + totpoprate[,,2] + totpoprate[,,3]
##### -- AN -- #####
## POINT ESTIMATES

# CANTON SPECIFIC ESTIMATES
ancanton[-27,-6,-4,,1] <- arraysim[,,,,1]

# canton by age 
ancanton[-27,6,-4,,1] <- apply(arraysim[,,,,1], c(1,3:4),sum)

# decade age group
ancanton[27,-6,-4,,1] <- apply(arraysim[,,,,1], 2:4,sum)
ancanton[27,-6,-4,,2] <- apply(apply(arraysim[,,,,-1], 2:5,sum), c(1:3),quantile, 0.025)
ancanton[27,-6,-4,,3] <- apply(apply(arraysim[,,,,-1], 2:5,sum), c(1:3),quantile, 0.975)


(((ancanton[27,-6,4,,1] /9)*10)/7874669)*10000
(((ancanton[27,-6,4,,2]  /9)*10)/7874669)*10000
(((ancanton[27,-6,4,,3] /9)*10)/8100000)*10000


# ESTIMATE TOTAL AGE
ancanton[-27,-6,4,,1] <- apply(arraysim[,,,,1], c(1:2,4),sum)
ancanton[-27,6,4,,1] <- apply(arraysim[,,,,1], c(1,4),sum)

# ESTIMATE TOTAL 
ancanton[27,6,4,,1] <- apply(arraysim[,,,,1], 4,sum)
ancanton[27,6,4,,2] <- apply(apply(arraysim[,,,,-1], c(4,5),sum),1,quantile, 0.025)
ancanton[27,6,4,,3] <- apply(apply(arraysim[,,,,-1], c(4,5),sum),1,quantile, 0.975)


#estimates by decade
ancanton[27,-6,4,,1] <- apply(arraysim[,,,,1], c(2,4),sum)
ancanton[27,-6,4,,2] <- apply(apply(arraysim[,,,,-1], c(2,4:5),sum), c(1:2), quantile,0.025)
ancanton[27,-6,4,,3] <- apply(apply(arraysim[,,,,-1], c(2,4:5),sum), c(1:2), quantile,0.975)


#estimates by age group OVERALL
ancanton[27,6,-4,,1] <- apply(arraysim[,,,,1], c(3,4),sum)
ancanton[27,6,-4,,2] <- apply(apply(arraysim[,,,,-1], c(3,4:5),sum), c(1:2), quantile,0.025)
ancanton[27,6,-4,,3] <- apply(apply(arraysim[,,,,-1], c(3,4:5),sum), c(1:2), quantile,0.975)




## CONFIDENCE INTERVALS
# COMPUTE CI CANTON
ancanton[-27,-6,-4,,2] <- apply(arraysim[,,,,-1],c(1:4), quantile,0.025)
ancanton[-27,-6,-4,,3] <- apply(arraysim[,,,,-1],c(1:4), quantile,0.975)

# ESTIMATE CI TOTAL DECADE
ancanton[-27,6,-4,,2] <- apply(apply(arraysim[,,,,-1], c(1,3:5),sum), c(1:3),quantile, 0.025)
ancanton[-27,6,-4,,3] <- apply(apply(arraysim[,,,,-1], c(1,3:5),sum), c(1:3),quantile, 0.975)



# ESTIMATE CI TOTAL CH



# ESTIMATE TOTAL AGE
ancanton[-27,-6,4,,2] <- apply(apply(arraysim[,,,,-1], c(1:2,4:5),sum), c(1:3),quantile, 0.025) 
ancanton[-27,-6,4,,3] <- apply(apply(arraysim[,,,,-1], c(1:2,4:5),sum), c(1:3),quantile, 0.975) 
ancanton[-27,6,4,,2] <- apply(apply(arraysim[,,,,-1], c(1,4:5),sum), c(1:2),quantile, 0.025) 
ancanton[-27,6,4,,3] <- apply(apply(arraysim[,,,,-1], c(1,4:5),sum), c(1:2),quantile, 0.975) 

(ancanton[-27,-6,4,,3] )




##### -- AF -- #####
## POINT ESTIMATES
# CANTON SPECIFIC ESTIMATES
afcanton[,,,"cold",1] <- (ancanton[,,,"cold",1] / totdeathall)*100
afcanton[,,,"cold",2] <- (ancanton[,,,"cold",2] / totdeathall)*100
afcanton[,,,"cold",3] <- (ancanton[,,,"cold",3] / totdeathall)*100

afcanton[,,,"heat",1] <- (ancanton[,,,"heat",1] / totdeathall)*100
afcanton[,,,"heat",2] <- (ancanton[,,,"heat",2] / totdeathall)*100
afcanton[,,,"heat",3] <- (ancanton[,,,"heat",3] / totdeathall)*100

afcanton[,,,"total",1] <- (ancanton[,,,"total",1] / totdeathall)*100
afcanton[,,,"total",2] <- (ancanton[,,,"total",2] / totdeathall)*100
afcanton[,,,"total",3] <- (ancanton[,,,"total",3] / totdeathall)*100






##### -- rates -- ##### per 100.000
## POINT ESTIMATES
# CANTON SPECIFIC ESTIMATES
ratecanton[,,,"cold",1] <- ((ancanton[,,,"cold",1]/10) / totpoprate)*100000
ratecanton[,,,"cold",2] <- ((ancanton[,,,"cold",2]/10) / totpoprate)*100000
ratecanton[,,,"cold",3] <- ((ancanton[,,,"cold",3]/10) / totpoprate)*100000

ratecanton[,,,"heat",1] <- ((ancanton[,,,"heat",1]/10) / totpoprate)*100000
ratecanton[,,,"heat",2] <- ((ancanton[,,,"heat",2]/10) / totpoprate)*100000
ratecanton[,,,"heat",3] <- ((ancanton[,,,"heat",3]/10) / totpoprate)*100000

ratecanton5th[,,,"cold",1] <- ((ancanton[,,,"cold",1]/9) / totpoprate)*100000
ratecanton5th[,,,"cold",2] <- ((ancanton[,,,"cold",2]/9) / totpoprate)*100000
ratecanton5th[,,,"cold",3] <- ((ancanton[,,,"cold",3]/9) / totpoprate)*100000

ratecanton5th[,,,"heat",1]  <- ((ancanton[,,,"heat",1]/9) / totpoprate)*100000
ratecanton5th[,,,"heat",2]  <- ((ancanton[,,,"heat",2]/9) / totpoprate)*100000
ratecanton5th[,,,"heat",3] <- ((ancanton[,,,"heat",3]/9) / totpoprate)*100000





## AF by decade
#estimates by decade
str(ancanton)
ancanton[27,-6,4,,1] <- apply(arraysim[,,,,1], c(2,4),sum)
ancanton[27,-6,4,,2] <- apply(apply(arraysim[,,,,-1], c(2,4:5),sum), c(1:2), quantile,0.025)
ancanton[27,-6,4,,3] <- apply(apply(arraysim[,,,,-1], c(2,4:5),sum), c(1:2), quantile,0.975)

## by decade 
AF_dec_total <- ancanton[27,1:6,"allages","total","est"]/totdeathall[27,,"allages"]*100
AF_dec_total_low <- ancanton[27,1:6,"allages","total",2]/totdeathall[27,,"allages"]*100
AF_dec_total_high <- ancanton[27,1:6,"allages","total",3]/totdeathall[27,,"allages"]*100

## by age group
AF_dec_total <- ancanton[27,6,1:3,1:3,"est"]/totdeathall[27,6,1:3]*100
AF_dec_total_low <- ancanton[27,6,1:3,1:3,2]/totdeathall[27,6,1:3]*100
AF_dec_total_high <- ancanton[27,6,1:3,1:3,3]/totdeathall[27,6,1:3]*100


## RATE 
str(totpoprate)
## by DECADE
AF_dec_total <- ancanton[27,1:6,"allages","total","est"]/10/totpoprate[27,1:6,4]*100000
AF_dec_total_low <- ancanton[27,1:6,"allages","total",2]/10/totpoprate[27,1:6,4]*100000
AF_dec_total_high <- ancanton[27,1:6,"allages","total",3]/10/totpoprate[27,1:6,4]*100000


## BY AGE GROUP
str(ancanton)
AF_dec_total <- ancanton[27,6,1:3,1:3,"est"]/10/totpoprate[27,6,1:3]*100000
AF_dec_total_low <- ancanton[27,6,1:3,1:3,2]/10/totpoprate[27,6,1:3]*100000
AF_dec_total_high <- ancanton[27,6,1:3,1:3,3]/10/totpoprate[27,6,1:3]*100000







## get total population
rowMeans(u65_pop[,1:5]) #5,911,385
rowMeans(pop_6579[,1:5]) #766,879
rowMeans(pop_80[,1:5]) #240,269

## get age and decade specific population
pop_u65_dec <- u65_pop[27,1:6]
pop_6579_dec <- pop_6579[27,1:6]
pop_80_dec <- pop_80[27,1:6]
tot_pop_dec <- pop_u65_dec + pop_6579_dec + pop_80_dec

totpoprate

#get decade and age specific percentages
perc_u65_dec <- pop_u65_dec/tot_pop_dec
perc_6579_dec <- pop_6579_dec/tot_pop_dec
perc_80_dec <- pop_80_dec/tot_pop_dec



####
std_heat_rate <- t(t(ratecanton[,,1,"heat",1])*perc_u65_dec)+
  t(t(ratecanton[,,2,"heat",1])*perc_6579_dec)+
  t(t(ratecanton[,,3,"heat",1])*perc_80_dec)

std_heat_ratelow <- t(t(ratecanton[,,1,"heat",2])*perc_u65_dec)+
  t(t(ratecanton[,,2,"heat",2])*perc_6579_dec)+
  t(t(ratecanton[,,3,"heat",2])*perc_80_dec)

std_heat_ratehigh <- t(t(ratecanton[,,1,"heat",3])*perc_u65_dec)+
  t(t(ratecanton[,,2,"heat",3])*perc_6579_dec)+
  t(t(ratecanton[,,3,"heat",3])*perc_80_dec)


std_cold_rate <-  t(t(ratecanton[,,1,"cold",1])*perc_u65_dec)+
  t(t(ratecanton[,,2,"cold",1])*perc_6579_dec)+
  t(t(ratecanton[,,3,"cold",1])*perc_80_dec)

std_cold_ratelow <-  t(t(ratecanton[,,1,"cold",2])*perc_u65_dec)+
  t(t(ratecanton[,,2,"cold",2])*perc_6579_dec)+
  t(t(ratecanton[,,3,"cold",2])*perc_80_dec)

std_cold_ratehigh <-  t(t(ratecanton[,,1,"cold",3])*perc_u65_dec)+
  t(t(ratecanton[,,2,"cold",3])*perc_6579_dec)+
  t(t(ratecanton[,,3,"cold",3])*perc_80_dec)

##decade 5
std_heat_rate5th <- t(t(ratecanton5th[,,1,"heat",1])*perc_u65_dec)+
  t(t(ratecanton5th[,,2,"heat",1])*perc_6579_dec)+
  t(t(ratecanton5th[,,3,"heat",1])*perc_80_dec)

std_cold_rate5th <- t(t(ratecanton5th[,,1,"cold",1])*perc_u65_dec)+
  t(t(ratecanton5th[,,2,"cold",1])*perc_6579_dec)+
  t(t(ratecanton5th[,,3,"cold",1])*perc_80_dec)


setwd("/Users/evandeschrijver/Desktop/PhD/publication1/final code/figures_june")

std_heat_rate[,5] <- std_heat_rate5th[,5]
std_cold_rate[,5] <- std_cold_rate5th[,5]
std_cold_rate <- as.data.frame(std_cold_rate)
std_heat_rate <- as.data.frame(std_heat_rate)

sum(std_heat_rate$`1999/08`)

write.csv(std_heat_rate,"std_heat_rate.csv")
write.csv(std_cold_rate,"std_cold_rate.csv")







 stdpop <- array(NA,dim=c(length(dlist)+1,6,4), dimnames=list(c(paste0("cant",seq(1:26)),"total"),
                                                              c("1969/78", "1979/88","1989/98","1999/08","2009/2017","mean"),
                                                              c("age1","age2","age3","allages")))
 

stdpop[,,1] <- t(t(totpoprate[,,4]) *perc_u65_dec)
stdpop[,,2] <- t(t(totpoprate[,,4]) *perc_6579_dec)
stdpop[,,3] <- t(t(totpoprate[,,4]) * perc_80_dec)
stdpop[,,4] <- t(stdpop[,,1]+stdpop[,,2]+stdpop[,,3])



## annual deaths check by age group 
u65_cold_deaths <- 5911385/100000*7
cold_6579_deaths <- 766879/100000*239
cold_80_deaths <- 240269/100000*1275





#### MMT MMP
minperccant <- as.data.frame(minperccant)
write_xlsx(minperccant, paste0(dir,"/results/MMTMMP/MMP_wide.xlsx"))
minperccant_long <- melt(minperccant, id.vars=c("canton"))

colnames(minperccant_long) <- c("canton", "age", "MMP")
minperccant_long$decade <- NA
write_xlsx(minperccant_long,paste0(dir,"/results/MMTMMP/MMP_long.xlsx"))





mintempcant <- as.data.frame(mintempcant)
mintempcant <- tibble::rowid_to_column(mintempcant, "canton")
write_xlsx(mintempcant, paste0(dir,"/results/MMTMMP/MMT_wide.xlsx"))

mintempcant_long <- melt(mintempcant, id.vars=c("canton"))

colnames(mintempcant_long) <- c("canton", "age", "MMP")
mintempcant_long$decade <- NA
write_xlsx(mintempcant_long,paste0(dir,"/results/MMTMMP/MMT_long.xlsx"))



library(tidyr)
library(reshape2)


ancanton <- as.data.frame(ancanton)
write_xlsx(ancanton,paste0(dir,"/results/MMTMMP/ancanton.xlsx"))

ratecanton <- as.data.frame(ratecanton)
write_xlsx(ratecanton,paste0(dir,"/rates/ratecantons.xlsx"))

plot(plots[[1]])
