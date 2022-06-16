

library(ecmwfr); library(RNetCDF); library(rgdal); library(raster); library(ncdf4)
library(datetime); library(maptools); library(lubridate) ; library(tidyverse) ; library(rnaturalearth); library(RColorBrewer); library(viridis); library(maps); library(mapdata); library(ncdf4); library(tictoc); library(exactextractr); library(sf); library(reshape2)


setwd("~/Desktop/PhD/github/Copernicus")
dir <- "~/Desktop/PhD/github/Copernicus"

# 1) Inlcude personal login data for API key
wf_set_key(user ="xxxx", key = "xxxx-xxxx-xxxxx-xxxxxx", service = "cds")

# 2) put in general informaiton of data extraction

# model
model <- "acces_cm2"
# global coordinates of the raster data you want coverage
Nc <- 90
Wc <- -180
Sc <- -90
Ec <- 180

#This can be useful when looping requests
cmip_model <- model
region_name <- "world"
cat("",cmip_model)


# 3)We make the request
request1 <- list(
  temporal_resolution = "daily",
  experiment = "historical",
  level = "single_levels",
  variable = "near_surface_air_temperature",
  model = paste0(cmip_model),
  area = c(Nc,Wc,Sc,Ec), #N,W,S,E
  format = "tgz",
  dataset_short_name = "projections-cmip6",
  target = "historical.tar.gz"
)

# 4) We pull the request (Do not forget to update your own UID number)
download_pull <- function(x) {
  wf_request(user = "xxxxx",
             request = x,  
             transfer = TRUE, 
             path = paste0(dir,cmip_model))
  
}


# 5) CREATE LIST OF REQUESTS
request_list <- list(request1)

#6 ) WE EXTRACT THE DATA DOWNLOAD IT
walk(request_list, download_pull)

#7 ) IF YOU DONWLOADED TAR FILE, THIS IS NECESSARY. IF DOWNLOADED A NETCDF FILE: SKIP THIS
untar(tarfile = "historical.tar.gz", exdir = "GFDL_historical")



CMIP5.proj <- "~/Desktop/PhD/github/Copernicus/temp_data"
#Set the colour for a plot
my.palette <- colorRampPalette(brewer.pal(9,"RdBu"))(100)


#  1 load the shapefile 
locations.sp  <- st_read("shapefile/district.shp")

# 2 load the raseterbirck 
##Get all the nc files within the CLIMM folder
CLIMM.NC.files <- list.files(paste0(CMIP5.proj))

#It has multiple folders because it was too big to download all at once
CLIMM.NC.files


###load the nc file
CLIMM.SSPRCP.NCfilek <- CLIMM.NC.files
CLIMM.SSPRCP.NCfilek.brick <- brick(paste0(CMIP5.proj,"/",CLIMM.SSPRCP.NCfilek), varname = "tas")

summary(CLIMM.SSPRCP.NCfilek.brick$X1980.01.01) ##tas (Temperature) is in Kelvin --> we will deal with that in a later step




#  3. We extract the raster data using the "exact_extract" function.     

# Finally extracting the raster data based on the polygons!

tic()
CLIMM.SSPRCP.NCfilek.daily <- exact_extract(CLIMM.SSPRCP.NCfilek.brick, locations.sp, fun="mean")

# create a common ID number nad name and bind the days to the locations
ID <-1:nrow(locations.sp)
NAME <- locations.sp$NAME
CLIMM.SSPRCP.NCfilek.daily <- cbind(ID,NAME,CLIMM.SSPRCP.NCfilek.daily)
toc()


#plot first 5 days with ID and name of district
head(CLIMM.SSPRCP.NCfilek.daily[,1:7])


#Slight issue, it is not time series format yet! Lets convert it into long format so we can actually work with it.
# 4. We melt the extracted data from wide to long format (which is preferred for time series)        
#Melt the dataframe to get the Date variable and Tmean.C 

CLIMM.SSPRCP.NCfilek.daily.df <- melt(id.vars=c("ID","NAME") ,CLIMM.SSPRCP.NCfilek.daily,variable.name="Date", value.name="Tmean.C")
colnames(CLIMM.SSPRCP.NCfilek.daily.df) <- c("ID","NAME","Date","Tmean.C")
head(CLIMM.SSPRCP.NCfilek.daily.df)


#Now the dates are still an issue. Lets use the gsub function to solve this. Lets also create an indicator for the  year and convert it all from Kelvin to degrees Celsuis (This is quicker  in a dataframe rather than doing it earlier in the raster itself)

# 5. We create year and date variables 
CLIMM.SSPRCP.NCfilek.daily.df$Date <- as.Date(gsub("mean.X", "", as.character(CLIMM.SSPRCP.NCfilek.daily.df$Date)), format=c("%Y.%m.%d"))
head(CLIMM.SSPRCP.NCfilek.daily.df) 

# NOW ALSO CREATE A YEAR VARIABLE
CLIMM.SSPRCP.NCfilek.daily.df$year <- year(CLIMM.SSPRCP.NCfilek.daily.df$Date)


# 6. We convert the data from Kelvin to degrees celsius.                     

# FROM KELIVIN TO TMEAN.C --> I do it in the end because this is quicker than running it over raster
CLIMM.SSPRCP.NCfilek.daily.df$Tmean.C <- CLIMM.SSPRCP.NCfilek.daily.df$Tmean.C - 273.15

head(CLIMM.SSPRCP.NCfilek.daily.df)



#We are finally there! We now have the ID, NAME, Date and temperature in degrees Celsius with an indicator for year in long format! But lets make it Tidy using the 'tidyr'package. These consist of tools to help create tidy data, where each column is a variable, each row is an observation, and each cell contains a single value. 'tidyr' contains tools for changing the shape (pivoting) and hierarchy (nesting and 'unnesting') of a dataset

## create tiblle with nested list of TS per municpiality --> maybe tidyr way of storing the data
CLIMM.SSPRCP.NCfilek.daily.list <- CLIMM.SSPRCP.NCfilek.daily.df %>% 
  select(ID,NAME,Date,Tmean.C,year) %>% 
  nest(TS=c(Date,year,Tmean.C))

#Now we have a list dataset by district in a dataframe! 
head(CLIMM.SSPRCP.NCfilek.daily.list)

#Lets see what that looks like for the first district:
head(CLIMM.SSPRCP.NCfilek.daily.list$TS[[3]])

# plot time series for first district
p <- ggplot(CLIMM.SSPRCP.NCfilek.daily.list$TS[[3]], aes(x=Date, y=Tmean.C)) +
  geom_line() + theme_bw()+
  xlab("")
p



#save the results
saveRDS(CLIMM.SSPRCP.NCfilek.daily.list, paste0(dir,"/GFDL_histroical_1973_1998.rds"))
