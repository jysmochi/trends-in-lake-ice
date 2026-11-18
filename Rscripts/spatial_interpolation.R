
rm(list = ls())

#> load packages
library(tidyverse)
library(sf)
library(stars)
library(rnaturalearth)
library(rnaturalearthdata)
library(usmap)
library(spatstat)
library(geodata)
library(Metrics)
library(leaflet)
library(paletteer)
library(maps)

#> read in data
ds <- read_csv("https://raw.githubusercontent.com/mingjys/trends-in-lake-ice/refs/heads/main/data/corrs20.csv")

ds <- ds %>%
  rename(lakename = lakename20,
         latitude = latitude20,
         longitude = longitude20,
         country = country20)

ds <- drop_na(ds)

fish <- function(x) {
  z <- 1/2 * log( (1+x)/(1-x) )
  return(z)
}

ds <- ds %>%
  arrange(latitude, -latitude) %>%
  filter(duplicated(latitude) == FALSE) %>%
  mutate(fi.corr = fish(corr))

################## Inverse Distance Weighting ##################

#> split ds
ds.fi <- ds %>%
  filter(country == 'FINLAND')
ds.nord <- ds %>%
  filter(country %in% c('FINLAND', 'SWEDEN'))
ds.usa <- ds %>%
  filter(country == 'UNITED STATES')

#> obtain map geography
# finland_adm <- gadm('finland', level=1)
# sweden_adm <- gadm('sweden', level=1)
# usa_adm <- gadm('usa', level=1)

#> set observation windows
finland_obs_window <- owin(xrange=c(20.55121, 31.5828), 
                           yrange=c(59.72097, 70.09145))
nord_obs_window <- owin(xrange=c(10.96139, 31.5828),
                        yrange=c(55.33625, 70.09145))
usa_obs_window <- owin(xrange=c(-93, -71),
                        yrange=c(40, 50))

#> create point pattern
ppp_finland <- ppp(ds.fi$longitude, ds.fi$latitude, 
                   marks=ds.fi$fi.corr, window=finland_obs_window)
ppp_nord <- ppp(ds.nord$longitude, ds.nord$latitude,
                marks=ds.nord$fi.corr, window=nord_obs_window)
ppp_usa <- ppp(ds.usa$longitude, ds.usa$latitude,
               marks=ds.usa$fi.corr, window = usa_obs_window)

#> run inverse distance weighting
idw_finland <- idw(ppp_finland, power=0.05, at="pixels")
idw_nord <- idw(ppp_nord, power = 0.05, at="pixels")
idw_usa <- idw(ppp_usa, power = 0.05, at='pixels')

#> plot idw objects
plot(idw_finland, col = heat.colors(20))
plot(idw_nord, col = heat.colors(20))
plot(idw_usa, col = heat.colors(20))

idw_finland_points <- idw(ppp_finland, power = 0.05, at="points")
plot(idw_finland_points, col=heat.colors(64))

#> find optimal power
powers <- seq(0.001, 10, 0.01)
mse_result <- NULL
for(power in powers) {
  cv_idw <- idw(ppp_finland, power=power, at="points")
  mse_result <- c(mse_result,
                  Metrics::mse(ppp_finland$marks, cv_idw))
}
optimal_power_fi <- powers[which.min(mse_result)]
mse_result <- NULL
for(power in powers) {
  cv_idw <- idw(ppp_nord, power=power, at="points")
  mse_result <- c(mse_result,
                  Metrics::mse(ppp_nord$marks, cv_idw))
}
optimal_power_nord <- powers[which.min(mse_result)]
mse_result <- NULL
for(power in powers) {
  cv_idw <- idw(ppp_usa, power=power, at="points")
  mse_result <- c(mse_result,
                  Metrics::mse(ppp_usa$marks, cv_idw))
}
optimal_power_usa <- powers[which.min(mse_result)]

#> run idw on nordic and us regions
#> note that the optimal power does not work for the nordic region
idw_nord <- idw(ppp_nord, power = optimal_power_usa, at="pixels", se=FALSE)
idw_usa <- idw(ppp_usa, power = optimal_power_usa, at='pixels', se=FALSE)

plot(idw_nord, col = heat.colors(20))
plot(idw_usa, col = heat.colors(20))

nord <- ne_countries(country = c("finland", "sweden"), type = "countries", scale = "medium")

usa2 <- ne_countries(country = "United States of America", type = "countries", scale = "large")

#> plots
idw_nord_ds <- as.data.frame(idw_nord) %>%
  rename("Correlation" = value)
ggplot() + 
  geom_raster(data = idw_nord_ds, aes(x=x, y=y, fill=Correlation)) + 
  scale_fill_gradient2(
    low="blue4",
    mid='purple',
    high='#ffffff',
    midpoint=-0.4
  ) + 
  geom_sf(data = nord, fill = NA, color = "black", size = 0.3) + 
  coord_sf() + 
  theme_minimal() + 
  labs(x = "Longitude", 
       y = "Latitude",
       title = "Nordic Fisher Z-Transformed Correlation IDW")

idw_usa_ds <- as.data.frame(idw_usa) %>%
  rename("Correlation" = value)
ggplot() + 
  geom_raster(data = idw_usa_ds, aes(x=x, y=y, fill=Correlation)) + 
  scale_fill_gradient2(
    low="blue4",
    mid='purple',
    high='#ffffff',
    midpoint=-0.4
  ) + 
  geom_sf(data = usa2, fill = NA, color = "black", size = 0.3) + 
  coord_sf(xlim = c(-92,-72),
           ylim = c(40.5,50)) + 
  theme_minimal()+ 
  labs(x = "Longitude", 
       y = "Latitude",
       title = "USA Fisher Z-Transformed Correlation IDW")

#> run idw on nordic and us regions >>> get se
#> note that the optimal power does not work for the nordic region
idw_nord_se <- idw(ppp_nord, power = optimal_power_usa, at="pixels", se=TRUE)
idw_usa_se <- idw(ppp_usa, power = optimal_power_usa, at='pixels', se=TRUE)

#> plots
idw_nord_se_ds <- as.data.frame(idw_nord_se$SE) %>%
  rename("SE" = value)
ggplot() + 
  geom_raster(data = idw_nord_se_ds, aes(x=x, y=y, fill=SE)) + 
  scale_fill_gradient2(
    low="white",
    mid='red',
    high='red4',
    midpoint=0.15
  ) + 
  geom_sf(data = nord, fill = NA, color = "black", size = 0.3) + 
  coord_sf() + 
  theme_minimal() + 
  labs(x = "Longitude", 
       y = "Latitude",
       title = "Nordic Fisher Z-Transformed Correlation IDW Standard Error")
idw_usa_se_ds <- as.data.frame(idw_usa_se$SE) %>%
  rename("SE" = value)
ggplot() + 
  geom_raster(data = idw_usa_se_ds, aes(x=x, y=y, fill=SE)) + 
  scale_fill_gradient2(
    low="white",
    mid='red',
    high='red4',
    midpoint=0.15
  ) + 
  geom_sf(data = usa2, fill = NA, color = "black", size = 0.3) + 
  coord_sf(xlim = c(-92,-72),
           ylim = c(40.5,50)) + 
  theme_minimal()+ 
  labs(x = "Longitude", 
       y = "Latitude",
       title = "USA Fisher Z-Transformed Correlation IDW Standard Error")
