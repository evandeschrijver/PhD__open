# packages

library(tidyverse)
library(sf)
library(rnaturalearth)
library(raster)
library(eurostat)
library(rcartocolor)
library(ggspatial)
library(elevatr)
library(cowplot)
library(ragg)
library(eurostat)



#############################
######### CH ELEVATION ######
#############################


# data
setwd("C:/Users/ed19w187/Desktop/PhD/PhD/publication_GCD/FINAL_CODE/dominic")
monitors <- st_read("./spatial/spatial/monitors.shp") # monitor stations
canton.pg <- st_read("./spatial/spatial/CH_canton.shp") %>% # cantons
  st_transform(st_crs(monitors)) %>% 
  st_zm()

# pop density (I tried something, but finally I didn't use it)
pop_dens <- raster("./spatial/spatial/gpw_v4_population_density_rev11_2020_30_sec.tif")
pop_dens <- crop(pop_dens, canton.pg) %>% mask(canton.pg)
pop_dens[pop_dens < 200] <- NA
pop_dens <- as.data.frame(pop_dens, xy = TRUE, na.rm = TRUE)
names(pop_dens)[3] <- "dens.pop"


# background maps

# for inset map
limits_map1 <- ne_countries(scale = 110, returnclass = "sf") %>% 
  st_crop(xmin = -15, xmax = 25, ymin = 33, ymax = 60)

limits_map2 <- get_eurostat_geospatial(resolution = "03", year = "2021", nuts_level = "0") %>% 
  st_crop(xmin = 5.7, xmax = 10.6, ymin = 45.7, ymax = 47.95)
plot(limits_map2)


lakes110 <- ne_download(scale=10, type = 'lakes', category = 'physical')
mapview(lakes110)


#roads <- ne_download(scale= 'large', type = 'roads', category = 'cultural')
#mapview(roads)

#roads <- ne_download(scale= 'large', type = 'airports', category = 'cultural')
#mapview(roads, zcol='type')

# lakes for main map
lakes <- st_read("./spatial/spatial/Large_lakes.shp") %>% 
  st_transform(st_crs(canton.pg)) %>% 
  st_intersection(canton.pg)

lakes_ch <- st_read("./spatial/spatial/lakes_ch.shp") %>% 
  st_transform(st_crs(canton.pg)) %>% 
  st_intersection(canton.pg)

mapview(lakes)

switzerland <- ne_countries(country = 'Switzerland', type='countries', scale='large')
try <- ne_download(country = 'Switzerland', type='countries', scale='large')



# eleveation data for main map
dtm <- get_elev_raster(canton.pg, z = 8)
dtm <- crop(dtm, canton.pg) %>% mask(canton.pg)
dtm <- as.data.frame(dtm, xy = TRUE, na.rm = TRUE)
names(dtm)[3] <- "alt"

# categories of elevation
dtm <- mutate(dtm, alt_cat = cut(alt, c(100, 300, 500, 700, 1000, 1500, 2000, 2500, 3000, 4000, 4500)))

# elevation colours
col_dtm <- c("#556231", "#6d7b46", "#888c5f", "#949257", "#9f9152", "#a1814a", "#936a4a", 
             "#855743", "#75433c", "#5c342f")


# inset map
inset_map <- st_as_sfc(st_bbox(canton.pg))

# projections
rob <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m no_defs"

# map limits
ext <- st_as_sfc(st_bbox(c(xmin = -11, xmax = 20, ymax = 33, ymin = 60), crs = 4326)) %>% 
  st_transform(rob) %>% 
  st_bbox()

# inset
m1 <- ggplot(inset_map) +
  geom_sf(data = limits_map1, fill = "grey90", size = .2, colour = "grey60") +
  geom_sf(fill = NA, colour = "red", linetype = "dashed", size = .4) +
  coord_sf(crs = rob, 
           expand = FALSE, xlim = ext[c(1, 3)], ylim = ext[c(2, 4)]) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.margin = margin(t = 0, r = 0, b = -5, l = -5, unit = "pt"))



# main map
density_ch_p <- ggplot() + 
  geom_sf(data = limits_map2, fill = "grey90", size = .4) +
  geom_tile(data = dtm, aes(x, y, fill = alt_cat)) +
  #geom_tile(data = pop_dens, aes(x, y), alpha = .9, fill = "black") +
  geom_sf(data = canton.pg, 
          col = "white",fill = NA,
          size = .6) + 
  geom_sf(data = lakes, fill = "#D6F1FF", colour = NA) +
  geom_sf(data = lakes_ch, fill = "#D6F1FF", colour = NA) +
  geom_sf(data = monitors, 
          shape = 21, 
          size = 3.1,  
          alpha = .6,
          fill = "#fc4e2a",
          colour = "black") +
  scale_fill_manual(values = col_dtm) +
  coord_sf(expand = FALSE) +
  annotation_scale(style = "ticks", width_hint = .1) +
  guides(fill = guide_colorsteps(barheight = 10, barwidth = .5)) +
  labs(fill = "Elevation (m)", x = "", y = "") + 
  theme_bw() +
  theme(legend.text.align = 1)

# joint inset map + main map and export
agg_png("map1.png", width = 12, height = 8, units = 'in', res = 300)
ggdraw() + 
  draw_plot(density_ch_p) +
  draw_plot(m1, x = 0.04, y = .72, width = 0.2, height = 0.2) # relative position depends on plot height + width
invisible(dev.off())






#############################
######### UK ELEVATION ######
#############################

monitors_ch <- st_read("./spatial/spatial/monitors.shp") 

canton.pg <- st_read("./spatial/spatial/EnglandWales_rgn_BNG.shp") %>% # monitor stations %>% # cantons
  st_transform(st_crs(monitors_ch)) %>% 
  st_zm()


monitors <- st_read("./spatial/spatial/monitors_UK.shp")%>% # monitor stations %>% # cantons
  st_transform(st_crs(monitors_ch)) %>% 
  st_zm()

# background maps

# for inset map
limits_map1 <- ne_countries(scale = 10, returnclass = "sf") %>% 
  st_crop(xmin = -15, xmax = 25, ymin = 33, ymax = 65)

# main map
limits_map2 <- get_eurostat_geospatial(resolution = "03", year = "2021", nuts_level = "0") %>% 
  st_crop(xmin = -6.5, xmax = 2.1, ymin = 49.7, ymax = 56.4)


# eleveation data for main map
dtm <- get_elev_raster(canton.pg, z = 8)
dtm <- crop(dtm, canton.pg) %>% mask(canton.pg)
dtm <- as.data.frame(dtm, xy = TRUE, na.rm = TRUE)
names(dtm)[3] <- "alt"

# categories of elevation
dtm <- mutate(dtm, alt_cat = cut(alt, c(-10,100,200,300,400,500,1000)))

# elevation colours
col_dtm <- c("#556231", "#888c5f","#9f9152", "#936a4a", 
              "#75433c", "#5c342f")


# inset map
inset_map <- st_as_sfc(st_bbox(canton.pg))


# projections
rob <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m no_defs"

# map limits
ext <- st_as_sfc(st_bbox(c(xmin = -11, xmax = 20, ymax = 33, ymin = 60), crs = 4326)) %>% 
  st_transform(rob) %>% 
  st_bbox()

# inset
m1 <- ggplot(inset_map) +
  geom_sf(data = limits_map1, fill = "grey90", size = .2, colour = "grey60") +
  geom_sf(fill = NA, colour = "red", linetype = "dashed", size = .4) +
  coord_sf(crs = rob, 
           expand = FALSE, xlim = ext[c(1, 3)], ylim = ext[c(2, 4)]) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.margin = margin(t = 0, r = 0, b = -5, l = -5, unit = "pt"))



# main map
density_ch_p <- ggplot() + 
  geom_sf(data = limits_map2, fill = "grey90", size = .4) +
  geom_tile(data = dtm, aes(x, y, fill = alt_cat)) +
  #geom_tile(data = pop_dens, aes(x, y), alpha = .9, fill = "black") +
  geom_sf(data = canton.pg, 
          col = "white",fill = NA,
          size = .6) + 
  geom_sf(data = monitors, 
          shape = 21, 
          size = 2,  
          alpha = .6,
          fill = "#fc4e2a",
          colour = "black") +
  scale_fill_manual(values = col_dtm) +
  coord_sf(expand = FALSE) +
  annotation_scale(style = "ticks", width_hint = .1) +
  guides(fill = guide_colorsteps(barheight = 10, barwidth = .5)) +
  labs(fill = "Elevation (m)", x = "", y = "") + 
  theme_bw() +
  theme(legend.text.align = 1)

# joint inset map + main map and export
agg_png("map2.png", width = 7, height = 14, units = 'in', res = 300)
ggdraw() + 
  draw_plot(density_ch_p) +
  draw_plot(m1, x = 0.59, y = .59, width = 0.2, height = 0.2) # relative position depends on plot height + width
invisible(dev.off())



















###################
######### UK ######
###################
setwd("C:/Users/ed19w187/Desktop/PhD/PhD/publication_GCD/FINAL_CODE/dominic")

monitors_ch <- st_read("./spatial/spatial/monitors.shp") 

UK.pg <- st_read("./spatial/spatial/EnglandWales_rgn_BNG.shp") %>% # monitor stations %>% # cantons
  st_transform(st_crs(monitors_ch)) %>% 
  st_zm()

monitors <- st_read("./spatial/spatial/monitors_UK.shp")%>% # monitor stations %>% # cantons
  st_transform(st_crs(monitors_ch)) %>% 
  st_zm()

ERA5UK <- st_read("ERA5_16h.shp")%>% # monitor stations %>% # cantons
  st_transform(st_crs(monitors_ch)) %>% 
  st_zm()

haduk5 <- st_read("./spatial/spatial/polygonraster.shp")
st_crs(haduk5) <- "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=m +no_defs"
haduk5 <- st_transform(haduk5, crs = st_crs(UK.pg))

# crop era5 to switzerland
London <- UK.pg[UK.pg$geo_label %in% "London",]
ERA5UK <- st_transform(ERA5UK, crs=st_crs(UK.pg))
Ldn_era5 <- ERA5UK[London,]


London <- UK.pg[UK.pg$geo_label %in% "London",]
haduk5 <- st_transform(haduk5, crs=st_crs(UK.pg))
Ldn_haduk5 <- haduk5[London,]

monitors <- st_transform(monitors, crs = st_crs(UK.pg))
Ldn_monitors<- monitors[London,]


extent(Ldn_haduk5)
r <- raster(ncol=180, nrow=180)
extent(r) <- extent(Ldn_haduk5)
Ldn_haduk52 <- rasterize(Ldn_haduk5,r,"X1993_01_0")
Ldn_haduk52 <- as.data.frame(Ldn_haduk52, xy = TRUE, na.rm = TRUE)
# categories of elevation
Ldn_haduk52 <- mutate(Ldn_haduk52, layer_cat = cut(layer, c(1,2,3,4,5,6)))


extent(Ldn_era5)
r <- raster(ncol=180, nrow=180)
extent(r) <- extent(Ldn_era5)
Ldn_era5 <- rasterize(Ldn_era5,r,"X1993_0")
# eleveation data for main map
Ldn_era52 <- as.data.frame(Ldn_era5, xy = TRUE, na.rm = TRUE)
# categories of elevation
Ldn_era52$layer <- Ldn_era52$layer-272.32
Ldn_era52 <- mutate(Ldn_era52, layer_cat = cut(layer, c(1,2,3,4,5,6)))
# elevation colours
col_ldn_tabsd <- rev(c("#B2182B", "#D6604D" ,"#F4A582", "#FDDBC7" ,"#D1E5F0"))
col_ldn_era<- rev(c("#B2182B", "#D6604D" ,"#F4A582", "#FDDBC7" ,"#D1E5F0"))


limits_Ldn <- get_eurostat_geospatial(resolution = "03", year = "2021", nuts_level = "1") %>% 
  st_crop(xmin = -1.2, xmax = 0.7, ymin = 51, ymax = 52.2)

maxt <- max(max(Ldn_era52$layer),max(Ldn_haduk52$layer))
# min of color scale
mint <- min(min(Ldn_era52$layer),min(Ldn_haduk52$layer))


UK_LDN <- rbind(limits2[8:9,],limits2[11,])

geom_sf(data = limits_ldn, fill = "grey90", size = .2, colour = "grey40")+
  geom_sf(data = UK_LDN, fill = "grey70", size = .2, colour = "grey40") +
  geom_tile(data = Ldn_era52, aes(x, y, fill = layer_cat)) +



# London  era
Ldn_era_map <- ggplot() + 
  geom_sf(data = limits_Ldn, 
          col = "grey40",fill = "grey80", alpha=0.6,
          size = .6) + 
  geom_tile(data = Ldn_era52, aes(x, y, fill = layer_cat)) +
  geom_sf(data = London, 
          col = "black",fill = NA,
          size = .6) + 
  geom_sf(data = Ldn_monitors, 
          shape = 21, 
          size = 2,  
          alpha = .6,
          fill = "#fc4e2a",
          colour = "black") +
  guides(fill = guide_colorsteps(barheight = 0.5, barwidth = 12, direction="horizontal")) +
  labs(fill = "?C", x = "", y = "") + 
  theme_minimal() +  annotation_scale(style = "ticks", line_width = 0.2) + 
  scale_fill_manual(values = col_ldn_era)+
  theme(legend.text.align = 1,
        legend.position = "bottom",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.ticks.y= element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# London TABSD 
Ldn_tabsd_map <- ggplot() + 
  geom_sf(data = limits_Ldn, 
          col = "grey40",fill = "grey80", alpha=0.6,
          size = .6) + 
  geom_tile(data = Ldn_haduk52, aes(x, y, fill = layer_cat)) +
  geom_sf(data = London, 
          col = "black",fill = NA,
          size = .6) + 
  geom_sf(data = Ldn_monitors, 
          shape = 21, 
          size = 2,  
          alpha = .6,
          fill = "#fc4e2a",
          colour = "black") +
  guides(fill = guide_colorsteps(barheight = 0.5, barwidth = 12, direction="horizontal")) +
  labs(fill = "?C", x = "", y = "") + 
  theme_minimal() +  annotation_scale(style = "ticks", line_width = 0.2) + 
  scale_fill_manual(values = col_ldn_tabsd) + 
  theme(legend.text.align = 1,
        legend.position = "bottom",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.ticks.y= element_blank()) 
#ldn_legend <- get_legend(Ldn_tabsd_map)

Ldn_total_map <- plot_grid(Ldn_tabsd_map,Ldn_era_map,ldn_legend, rel_widths = c(1,1,0.4), ncol=3)




# for inset map
limits_ldn <- ne_countries(scale = 10, returnclass = "sf") %>% 
    st_crop(xmin = -1.2, xmax = 0.7, ymin = 51, ymax = 52.2)


limits2 <- get_eurostat_geospatial(resolution = "03", year = "2021", nuts_level = "1") %>% 
  st_crop(xmin = -7, xmax = 2, ymin = 49.5, ymax = 58)

# inset map
inset_map <- st_as_sfc(st_bbox(limits_Ldn))

# projections
rob <- "+proj=longlat +datum=WGS84 +no_defs"

# map limits
ext <- st_as_sfc(st_bbox(c(xmin = -7, xmax = 2, ymin = 49.5, ymax = 58), crs = 4326)) %>% 
  st_transform(rob) %>% 
  st_bbox()

LDN <- UK.pg[London,]
LDN <- LDN[3,]

# inset
m1 <- ggplot(inset_map) +
  geom_sf(data = limits2, fill = "grey90", size = .2, colour = "grey40") +
  geom_sf(data = UK.pg, fill = "grey80", size = .2, colour = "grey40") +
  geom_sf(data = LDN, fill = "red", size = .2, colour = "grey40", alpha=0.6) +
  geom_sf(fill = NA, colour = "red", linetype = "dashed", size = .4) +
  coord_sf(crs = rob, 
           expand = FALSE, xlim = ext[c(1, 3)], ylim = ext[c(2, 4)]) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))





# joint inset map + main map and export
agg_png("map_ldn.png", width = 7, height = 14, units = 'in', res = 300)
ggdraw() + 
  draw_plot(Ldn_tabsd_map) +
  draw_plot(m1, x = 0.04, y = 0.53, width = 0.40, height = 0.38) # relative position depends on plot height + width
invisible(dev.off())
































































###################
#########CH VALAIS ######
###################

setwd("C:/Users/ed19w187/Desktop/PhD/PhD/publication_GCD/FINAL_CODE/dominic")
monitors <- st_read("./spatial/spatial/monitors.shp") # monitor stations
canton.pg <- st_read("./spatial/spatial/CH_canton.shp") #%>% # cantons
  st_transform(st_crs(monitors)) %>% 
  st_zm()


setwd("C:/Users/ed19w187/Desktop/PhD/PhD/publication_GCD/FINAL_CODE/dominic")
monitors  <- st_read("./spatial/spatial/monitors.shp") # monitor stations
canton.pg <- st_read("./spatial/spatial/CH_canton.shp") %>% # cantons
  st_transform(st_crs(monitors)) %>% 
  st_zm()

setwd("C:/Users/ed19w187/Desktop/PhD/PhD/publication_GCD/shapefiles")
# Tabsd
Tabsd <- st_read("Tabsd_layer.shp")
era5_10 <- st_read("era5_10.shp") %>% # cantons
  st_transform(st_crs(monitors)) %>% 
  st_zm()

# == reformat if needed
# crs of Tabsd to WGS84
Tabsd <- st_transform(Tabsd, crs = st_crs(monitors))
# crs monitors to WGS84
era5_10 <- st_transform(era5_10, crs = st_crs(canton.pg))


# crop era5 to switzerland
Valais <- canton.pg[canton.pg$NAME %in% "Valais",]
st_crs(era5_10) <- crs(canton.pg)
Valais_era5 <- era5_10[Valais,]
plot(Valais_era5)


st_crs(Tabsd) <- crs(canton.pg)
Valais_tabsd <- Tabsd[Valais,]
plot(Valais_tabsd)


Valais_monitors<- monitors[Valais,]
plot(Valais_monitors)


# kelvin to celsius era5
Valais_era5$X2010_0 <- Valais_era5$X2010_0 - 273.15 



extent(Valais_tabsd)
r <- raster(ncol=180, nrow=180)
extent(r) <- extent(Valais_tabsd)
Valais_tabsd <- rasterize(Valais_tabsd,r,"X2010_0")
# eleveation data for main map
Valais_tabsd2 <- as.data.frame(Valais_tabsd, xy = TRUE, na.rm = TRUE)
# categories of elevation
Valais_tabsd2 <- mutate(Valais_tabsd2, layer_cat = cut(layer, c(-18,-15,-12,-9,-6,-3,0,5)))


extent(Valais_era5)
r <- raster(ncol=180, nrow=180)
extent(r) <- extent(Valais_era5)
Valais_era5 <- rasterize(Valais_era5,r,"X2010_0")
# eleveation data for main map
Valais_era52 <- as.data.frame(Valais_era5, xy = TRUE, na.rm = TRUE)
# categories of elevation
Valais_era52 <- mutate(Valais_era52, layer_cat = cut(layer, c(-18,-15,-12,-9,-6,-3,0,5)))
# elevation colours
col_valais_era <- rev(c("#B2182B", "#D6604D" ,"#F4A582", "#FDDBC7" ,"#D1E5F0" ,"#92C5DE", "#4393C3"))




limits_map2 <- get_eurostat_geospatial(resolution = "03", year = "2021", nuts_level = "2") %>% 
  st_crop(xmin = 6.3, xmax = 8.9, ymin = 45.5, ymax = 47)


CH_vs <- rbind(limits_map2[1:4,], limits_map2[8,])

# VALAIS TABSD 
Valais_tabsd_map <- ggplot() + 
  geom_sf(data = limits_map2, fill = "grey90", size = .2, colour = "grey40")+
  geom_sf(data = CH_vs, fill = "grey80", size = .2, colour = "grey40") +
  geom_tile(data = Valais_tabsd2, aes(x, y, fill = layer_cat)) +
  geom_sf(data = Valais, 
          col = "black",fill = NA,
          size = .6) + 
  geom_sf(data = Valais_monitors, 
          shape = 21, 
          size = 3,  
          alpha = .6,
          fill = "#fc4e2a",
          colour = "black") +
  guides(fill = guide_colorsteps(barheight = 15, barwidth = 1,direction="vertical"))+
  labs(fill = "?C", x = "", y = "") + 
  theme_minimal() +  annotation_scale(style = "ticks", line_width = 0.2) + 
  scale_fill_manual(values = col_valais_era) +
  theme(legend.text.align = 1,
        legend.position = "right",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.ticks.y= element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#legend_valais <- get_legend(Valais_tabsd_map)


# VALAIS TABSD 
Valais_era5_map <-ggplot() + 
  geom_sf(data = limits_map2, fill = "grey90", size = .2, colour = "grey40")+
  geom_sf(data = CH_vs, fill = "grey80", size = .2, colour = "grey40") +
  geom_tile(data = Valais_era52, aes(x, y, fill = layer_cat))+
  geom_sf(data = Valais, 
          col = "black",fill = NA,
          size = .6) + 
  geom_sf(data = Valais_monitors, 
          shape = 21, 
          size = 3,  
          alpha = .6,
          fill = "#fc4e2a",
          colour = "black") +
  guides(fill = guide_colorsteps(barheight = 6, barwidth = 2)) +
  labs(fill = "Temperature ?C", x = "", y = "") + 
  theme_minimal() +  annotation_scale(style = "ticks", line_width = 0.2) + 
  scale_fill_manual(values = col_valais_era) +
  theme(legend.text.align = 1,
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.ticks.y= element_blank()) 

valais_total <- plot_grid(Valais_tabsd_map,Valais_era5_map, legend_valais, ncol=3, rel_widths = c(1,1,0.3))




# for inset map
limits_valais <- ne_countries(scale = 10, returnclass = "sf") %>% 
  st_crop(xmin = 6.7, xmax = 8.5, ymin = 45.88, ymax = 46.7)
limits_valais2 <- limits_valais[3,]
limits_valais2 <- canton.pg[Valais,]
limits_valais2 <- limits_valais2[4,]

limits2 <- get_eurostat_geospatial(resolution = "03", year = "2021", nuts_level = "2") %>% 
  st_crop(xmin = 5.8, xmax = 10.6, ymin = 45.8, ymax = 47.9)

# inset map
inset_map <- st_as_sfc(st_bbox(limits_valais))

# projections
rob <- "+proj=longlat +datum=WGS84 +no_defs"


# map limits
ext <- st_as_sfc(st_bbox(c(xmin=5.8, xmax = 10.6, ymin = 45.8, ymax = 47.9), crs = 4326)) %>% 
  st_transform(rob) %>% 
  st_bbox()



# inset
m1 <- ggplot(inset_map) +
  geom_sf(data = limits2, fill = "grey90", size = .2, colour = "grey40") +
  geom_sf(data = canton.pg, fill = "grey80", size=.2, colour="grey40") +
  geom_sf(data = limits_valais2, fill = "red", alpha=0.6) +

  geom_sf(fill = NA, colour = "red", linetype = "dashed", size = .4) +
  coord_sf(crs = rob, 
           expand = FALSE, xlim = ext[c(1, 3)], ylim = ext[c(2, 4)]) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))



# joint inset map + main map and export
agg_png("map_ldn.png", width = 7, height = 14, units = 'in', res = 300)
ggdraw() + 
  draw_plot(Valais_tabsd_map) +
  draw_plot(m1, x = 0.13, y = .64, width = 0.3, height = 0.3) # relative position depends on plot height + width
invisible(dev.off())















############################
######### CH - whole  ######
############################
setwd("C:/Users/ed19w187/Desktop/PhD/PhD/Rfiles/Switzerland/data/density data/popcount_1km_2010/CHdensity")
densityCH <- st_read("weightedtempcell.shp")

setwd("C:/Users/ed19w187/Desktop/PhD/PhD/publication_GCD/FINAL_CODE/dominic")
monitors <- st_read("./spatial/spatial/monitors.shp") # monitor stations
canton.pg <- st_read("./spatial/spatial/CH_canton.shp") %>% # cantons
  st_transform(st_crs(monitors)) %>% 
  st_zm()


setwd("C:/Users/ed19w187/Desktop/PhD/PhD/publication_GCD/FINAL_CODE/dominic")
monitors  <- st_read("./spatial/spatial/monitors.shp") # monitor stations
canton.pg <- st_read("./spatial/spatial/CH_canton.shp") %>% # cantons
  st_transform(st_crs(monitors)) %>% 
  st_zm()

setwd("C:/Users/ed19w187/Desktop/PhD/PhD/publication_GCD/shapefiles")
# Tabsd
Tabsd <- st_read("Tabsd_layer.shp")
era5_10 <- st_read("era5_10.shp") %>% # cantons
  st_transform(st_crs(monitors)) %>% 
  st_zm()

# == reformat if needed
# crs of Tabsd to WGS84
Tabsd <- st_transform(Tabsd, crs = st_crs(monitors))
# crs monitors to WGS84
era5_10 <- st_transform(era5_10, crs = st_crs(canton.pg))
era5_10 <- era5_10[canton.pg,]

densityCH <- st_transform(densityCH, crs = st_crs(monitors))
densityCH <- densityCH[canton.pg,]

# kelvin to celsius era5
era5_10$X2010_0 <- era5_10$X2010_0 - 273.15 


# elevation colours
col_era <- c("#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD")





# ERA5
CH_era5_map <- ggplot() + 
  geom_sf(data = era5_10, aes(fill =  X2010_0)) +
  #geom_tile(data = pop_dens, aes(x, y), alpha = .9, fill = "black") +
  geom_sf(data = canton.pg, 
          col = "white",fill = NA,
          size = .6) + 
  geom_sf(data = monitors, 
          shape = 21, 
          size = 2,  
          alpha = .6,
          fill = "#fc4e2a",
          colour = "black") +
  scale_fill_manual(values = col_era) +
  guides(fill = guide_colorsteps(barheight = 6, barwidth = 2)) +
  labs(fill = "Temperature ?C", x = "", y = "") + 
  theme_minimal() +
  theme(legend.text.align = 1,
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y= element_blank())



# density
density_ch_p <- ggplot() + 
  geom_sf(data = densityCH, aes(fill =  Sum_g_4___)) +
  #geom_tile(data = pop_dens, aes(x, y), alpha = .9, fill = "black") +
  geom_sf(data = canton.pg, 
          col = "black",fill = NA,
          size = 1) + 
  geom_sf(data = monitors, 
          shape = 21, 
          size = 2,  
          alpha = .6,
          fill = "#fc4e2a",
          colour = "black") +
  scale_fill_viridis_c(breaks = round(exp(c(4,5,6,7,8,9))), trans = "log")+
  guides(fill = guide_colorsteps(barheight = 6, barwidth = 2)) +
  labs(fill = "Population density", x = "", y = "") + 
  theme_minimal() +
  theme(legend.text.align = 1,
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y= element_blank())













#############################
######### CH DENSITY ######
#############################
setwd("C:/Users/ed19w187/Desktop/PhD/PhD/Rfiles/Switzerland/data/density data/popcount_1km_2010/CHdensity")
densityCH <- st_read("weightedtempcell.shp")


# data
setwd("C:/Users/ed19w187/Desktop/PhD/PhD/publication_GCD/FINAL_CODE/dominic")
monitors <- st_read("./spatial/spatial/monitors.shp") # monitor stations
canton.pg <- st_read("./spatial/spatial/CH_canton.shp") %>% # cantons
  st_transform(st_crs(monitors)) %>% 
  st_zm()


densityCH <- st_transform(densityCH, crs=st_crs(canton.pg))
densityCH <- densityCH[canton.pg,]



extent(densityCH)
r <- raster(ncol=180, nrow=180)
extent(r) <- extent(densityCH)
densityCH <- rasterize(densityCH,r,"Sum_g_4___")
# eleveation data for main map
densityCH2 <- as.data.frame(densityCH, xy = TRUE, na.rm = TRUE)
# categories of elevation
densityCH2 <- mutate(densityCH2, layer_cat = cut(layer, c(-5,50,200,400,1000,2000,4000,8000,12000,30000)))
# elevation colours




# pop density (I tried something, but finally I didn't use it)
pop_dens <- raster("./spatial/spatial/gpw_v4_population_density_rev11_2020_30_sec.tif")
pop_dens <- crop(pop_dens, canton.pg) %>% mask(canton.pg)
pop_dens <- as.data.frame(pop_dens, xy = TRUE, na.rm = TRUE)
names(pop_dens)[3] <- "dens.pop"
pop_dens <- mutate(pop_dens, layer_cat = cut(dens.pop, c(-5,10,20,50,150,400,1100,3000,8000,30000)))






library(scales)
q_colors =  10
v_colors =  viridis(q_colors, option= "D")
collist <- c(v_colors)

col_valais_dens<- c("#440154FF", "#482878FF", "#31688EFF", "#26828EFF", "#1F9E89FF", "#35B779FF", "#6DCD59FF", "#B4DE2CFF", "#FDE725FF")


setwd("C:/Users/ed19w187/Desktop/PhD/PhD/publication_GCD/FINAL_CODE/dominic")

# background maps

# for inset map
limits_map1 <- ne_countries(scale = 10, returnclass = "sf") %>% 
  st_crop(xmin = -15, xmax = 25, ymin = 33, ymax = 60)

# main map
limits_map2 <- get_eurostat_geospatial(resolution = "03", year = "2021", nuts_level = "0") %>% 
  st_crop(xmin = 5.7, xmax = 10.6, ymin = 45.7, ymax = 48)

# lakes for main map
lakes <- st_read("./spatial/spatial/Large_lakes.shp") %>% 
  st_transform(st_crs(canton.pg)) %>% 
  st_intersection(canton.pg)

lakes_ch <- st_read("./spatial/spatial/lakes_ch.shp") %>% 
  st_transform(st_crs(canton.pg)) %>% 
  st_intersection(canton.pg)


# inset map
inset_map <- st_as_sfc(st_bbox(canton.pg))
plot(inset_map)
# projections
rob <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m no_defs"

# map limits
ext <- st_as_sfc(st_bbox(c(xmin = -11, xmax = 20, ymax = 33, ymin = 60), crs = 4326)) %>% 
  st_transform(rob) %>% 
  st_bbox()

limits_map1 <- mutate(limits_map1, dummy = ifelse(name_nl == "Zwitserland", "yes", "no"))
try <- subset(limits_map1, name_nl="Zwitserland")
limits_mapd <- limits_map1[canton.pg,]
limits_mapd <- limits_mapd[5,]


try <- limits_map1$dummy
col_dummy <- c("grey90","grey40")


# inset
m1 <- ggplot(inset_map) +
  geom_sf(data = limits_map1, fill = "grey90", size = .2, colour = "grey40") +
  geom_sf(data = limits_mapd, fill = "grey70", size=0.2, colour="grey40") +
  geom_sf(fill = NA, colour = "red", linetype = "dashed", size = .4) +
  coord_sf(crs = rob, 
           expand = FALSE, xlim = ext[c(1, 3)], ylim = ext[c(2, 4)]) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.margin = margin(t = 0, r = 0, b = -5, l = -5, unit = "pt")) 



# main map
density_ch_p <- ggplot() + 
  geom_sf(data = limits_map2, fill = "grey90", size = .4) +
 # geom_tile(data = densityCH2, aes(x, y, fill = layer_cat)) +
  geom_tile(data = pop_dens, aes(x, y,fill=dens.pop)) +
  geom_sf(data = canton.pg, 
          col = "white",fill = NA,
          size = .6) + 
  geom_sf(data = lakes, fill = "#D6F1FF", colour = NA) +
  geom_sf(data = lakes_ch, fill = "#D6F1FF", colour = NA) +
  geom_sf(data = monitors, 
          shape = 21, 
          size = 3.1,  
          alpha = .6,
          fill = "#fc4e2a",
          colour = "black") +
  #scale_fill_manual(values = col_valais_dens) +
  scale_fill_viridis_c(breaks = round(exp(c(2,4,5,6,7,8,9,10))), trans = "log") +
  annotation_scale(style = "ticks", line_width = 0.2, 
                        height = unit(0.5, "cm")) + 
  
  coord_sf(expand = FALSE) + # to put the scale into the map!!!!! 
  guides(fill = guide_colorsteps(barheight = 0.5, barwidth = 20,direction="horizontal",  label.position = "top")) +
  labs(fill = "Population density", x = "", y = "") + 
  theme_minimal() +
  theme(legend.text.align = 1,
        legend.position = "bottom",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.ticks.y= element_blank()) 

# joint inset map + main map and export
agg_png("map1.png", width = 12, height = 8, units = 'in', res = 300)
ggdraw() + 
  draw_plot(density_ch_p) +
  draw_plot(m1, x = 0.1, y = .755, width = 0.21, height = 0.21) # relative position depends on plot height + width
invisible(dev.off())






























#############################
######### UK DENSITY ######
#############################
setwd("C:/Users/ed19w187/Desktop/PhD/PhD/Rfiles/Switzerland/data/density data/popcount_1km_2010/CHdensity")




###################
######### UK ######
###################
setwd("C:/Users/ed19w187/Desktop/PhD/PhD/publication_GCD/FINAL_CODE/dominic")

monitors_ch <- st_read("./spatial/spatial/monitors.shp") 

UK.pg <- st_read("./spatial/spatial/EnglandWales_rgn_BNG.shp") %>% # monitor stations %>% # cantons
  st_transform(st_crs(monitors_ch)) %>% 
  st_zm()

monitors <- st_read("./spatial/spatial/monitors_UK.shp")%>% # monitor stations %>% # cantons
  st_transform(st_crs(monitors_ch)) %>% 
  st_zm()

ERA5UK <- st_read("ERA5_16h.shp")%>% # monitor stations %>% # cantons
  st_transform(st_crs(monitors_ch)) %>% 
  st_zm()

haduk5 <- st_read("./spatial/spatial/polygonraster.shp")
st_crs(haduk5) <- "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=m +no_defs"
haduk5 <- st_transform(haduk5, crs = st_crs(UK.pg))

# crop era5 to switzerland
London <- UK.pg[UK.pg$geo_label %in% "London",]
ERA5UK <- st_transform(ERA5UK, crs=st_crs(UK.pg))
Ldn_era5 <- ERA5UK[London,]



# pop density (I tried something, but finally I didn't use it)
pop_dens <- raster("./spatial/spatial/gpw_v4_population_density_rev11_2020_30_sec.tif")
pop_dens <- crop(pop_dens, UK.pg) %>% mask(UK.pg)
pop_dens <- as.data.frame(pop_dens, xy = TRUE, na.rm = TRUE)
names(pop_dens)[3] <- "dens.pop"
pop_dens <- mutate(pop_dens, layer_cat = cut(dens.pop, c(-5,10,20,50,150,400,1100,3000,8000,30000)))


# for inset map
limits_map1 <- ne_countries(scale = 10, returnclass = "sf") %>% 
  st_crop(xmin = -15, xmax = 25, ymin = 33, ymax = 65)

# main map
limits_map2 <- get_eurostat_geospatial(resolution = "03", year = "2021", nuts_level = "0") %>% 
  st_crop(xmin = -6.5, xmax = 2.1, ymin = 49.7, ymax = 56.4)





# main map
limits_map2 <- get_eurostat_geospatial(resolution = "03", year = "2021", nuts_level = "0") %>% 
  st_crop(xmin = -7, xmax = 5.5, ymin = 47.5 , ymax = 60)






# inset map
inset_map <- st_as_sfc(st_bbox(UK.pg))

# projections
rob <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m no_defs"

# map limits
ext <- st_as_sfc(st_bbox(c(xmin = -11, xmax = 20, ymax = 33, ymin = 60), crs = 4326)) %>% 
  st_transform(rob) %>% 
  st_bbox()



limits_mapd <- limits_map1[UK.pg,]


# inset
m1 <- ggplot(inset_map) +
  geom_sf(data = limits_map1, fill = "grey90", size = .2, colour = "grey40") +
  geom_sf(data = UK.pg, fill = "grey80", size = .2,colour="grey40") +
  geom_sf(fill = NA, colour = "red", linetype = "dashed", size = .4) +
  coord_sf(crs = rob, 
           expand = FALSE, xlim = ext[c(1, 3)], ylim = ext[c(2, 4)]) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.margin = margin(t = 0, r = 0, b = -5, l = -5, unit = "pt"))



# main map
density_uk_p <- ggplot() + 
  geom_sf(data = limits_map2, fill = "grey90", size = .4) +
  geom_tile(data = pop_dens, aes(x, y,fill=dens.pop)) +
  geom_sf(data = canton.pg, 
          col = "white",fill = NA,
          size = .6) + 
  geom_sf(data = monitors, 
          shape = 21, 
          size = 2.1,  
          alpha = .6,
          fill = "#fc4e2a",
          colour = "black") +
  #scale_fill_manual(values = col_valais_dens) +
  scale_fill_viridis_c(breaks = round(exp(c(2,4,5,6,7,8,9,10))), trans = "log") +
  coord_sf(expand = FALSE) +
  annotation_scale(style = "ticks", width_hint = .1) +
  guides(fill = guide_colorsteps(barheight = 0.5, barwidth = 20,direction="horizontal",  label.position = "top")) +
  labs(fill = "Population density", x = "", y = "") + 
  theme_minimal() +
  theme(legend.text.align = 1,
        legend.position = "bottom",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.ticks.y= element_blank())






# joint inset map + main map and export
agg_png("map1.png", width = 12, height = 8, units = 'in', res = 300)
ggdraw() + 
  draw_plot(density_uk_p) +
  draw_plot(m1, x = 0.57, y = 0.68, width = 0.22, height = 0.22) # relative position depends on plot height + width

ggdraw() + 
  draw_plot(density_uk_p) +
  draw_plot(m1, x = 0.55, y = 0.70, width = 0.22, height = 0.22) # relative position depends on plot height + width

invisible(dev.off())






















