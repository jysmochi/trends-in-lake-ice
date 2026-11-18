rm(list = ls())

library(naniar)
library(ggpubr)
library(gridExtra)
library(car)
library(VGAM)
library(funtimes)
library(Rssa)
library(forecast)
library(class)
library(ncdf4)
library(raster)
library(sp)
library(tidyverse)
library(sf)
library(this.path)

#################################################################################################################

# set up and functions

#>>># NOTE: to be able to download the brightness temperature data 
#>>>#       from the github repo, you must first save this file locally. 
#>>>#       the code below downloads the brightness temperature data 
#>>>#       to the location of this file. 

#> set dir to source file location
file_path <- this.path::here()
setwd(file_path)

#> read data
ice = read.csv("https://raw.githubusercontent.com/mingjys/trends-in-lake-ice/refs/heads/main/data/liag_ice.csv", header=TRUE)
download.file(url="https://github.com/mingjys/trends-in-lake-ice/raw/refs/heads/main/data/uat4_tb_v04r00_avrg_chtlt_s197812_e202312_c20240104.nc", destfile="uat4_tb_v04r00_avrg_chtlt_s197812_e202312_c20240104.nc", mode="wb")
nc_data = nc_open("uat4_tb_v04r00_avrg_chtlt_s197812_e202312_c20240104.nc")
phys_char <- read.csv("https://raw.githubusercontent.com/mingjys/trends-in-lake-ice/refs/heads/main/data/liag_info.csv", header=TRUE)

long = ncvar_get(nc_data, "longitude")
lat = ncvar_get(nc_data, "latitude")
tb = ncvar_get(nc_data, "brightness_temperature")

nc_close(nc_data)

#> filter and clean data
ice2 = ice %>%
  dplyr::select(c(iceon_year, iceoff_year, duration, latitude, longitude, lakename, lakeorriver, country, lakecode)) %>% 
  filter(iceon_year >= 1978,
         iceoff_year >= 1978,
         duration > 0,
         lakeorriver == " L") %>%
  mutate(lakename = trimws(lakename),
         lakeorriver = trimws(lakeorriver),
         country = trimws(country),
         lakecode = trimws(lakecode))

lake_info <- phys_char %>%
  mutate(lakename = trimws(lakename),
         lakecode = trimws(lakecode),
         continent = trimws(continent),
         country = trimws(country),
         state = trimws(state),
         contributor = trimws(contributor),
         inlet_streams = trimws(inlet_streams),
         power_plant_discharge = trimws(power_plant_discharge)) %>%
  dplyr::select(-comments)

count_ice_35 = ice2 %>%
  group_by(lakecode) %>%
  dplyr::summarise(n = n()) %>%
  filter(n > 35)

more_lakes_35 = count_ice_35$lakecode

count_ice_30 = ice2 %>%
  group_by(lakecode) %>%
  dplyr::summarise(n = n()) %>%
  filter(n > 30)

more_lakes_30 = count_ice_30$lakecode

not35_30 <- setdiff(more_lakes_30, more_lakes_35)

count_ice_20 = ice2 %>%
  group_by(lakecode) %>%
  dplyr::summarise(n = n()) %>%
  filter(n > 20)

more_lakes_20 = count_ice_20$lakecode

not35_20 <- setdiff(more_lakes_20, more_lakes_35)

#> input: lakecode
#> output: filtered ice2 data to only include data from lakecode
f_lake_data = function(code) {
  ds_lake = ice2 %>%
    filter(lakecode == code) %>%
    arrange(iceoff_year)
  return(ds_lake)
}

#> input: lakecode
#> output: logit transform and run ssa on icd from lakecode
f_icessa = function(code) {
  spec = f_lake_data(code)
  ice = spec$duration
  ice_ts = ts(ice, start = 1, frequency=1)
  a = min(ice) * 0.9
  b = max(ice) * 1.1
  icelog = log((ice_ts-a)/(b-ice_ts))
  ssaout = ssa(icelog, L = ceiling(length(ice_ts)/2))
  return(ssaout)
}

#> input: lakecode
#> output: print output ssa plots for lakecode data
f_ssa_plots = function(code) {
  specific_lake = f_lake_data(code)
  name = specific_lake$lakename[1]
  code = specific_lake$lakecode[1]
  spec_lake_ssa = f_icessa(code)
  
  print(plot(spec_lake_ssa, type = "values", main = paste0(name, " ", code, " Eigenvalues")))
  print(plot(spec_lake_ssa, type = "vectors", main = paste0(name, " ", code, " Eigenvectors")))
  print(plot(spec_lake_ssa, type = "paired", main = paste0(name, " ", code, " Paired Eigenvectors")))
  print(plot(wcor(spec_lake_ssa), main = paste0(name, " ", code, " W-Correlation Matrix")))  
}

#> input: lakecode
#> output: print output ssa plots in individual pdfs
f_ssa_plots_ind = function(code) {
  specific_lake = f_lake_data(code)
  name = specific_lake$lakename[1]
  code = specific_lake$lakecode[1]
  spec_lake_ssa = f_icessa(code)
  
  pdf(paste0(name, " values.pdf"))
  print(plot(spec_lake_ssa, type = "values", main = paste0(name, " ", code, " Singular Values")))
  dev.off()
  pdf(paste0(name, " vectors.pdf"))
  print(plot(spec_lake_ssa, type = "vectors", main = paste0(name, " ", code, " Singular Vectors")))
  dev.off()
  pdf(paste0(name, " paired.pdf"))
  print(plot(spec_lake_ssa, type = "paired", main = paste0(name, " ", code, " Paired Singular Vectors")))
  dev.off()
  pdf(paste0(name, " wcor.pdf"))
  print(plot(wcor(spec_lake_ssa), main = paste0(name, " ", code, " W-Correlation Matrix")))  
  dev.off()
}

#> example usage and testing
# mylake <- f_lake_data("JK48")
# name <- mylake$lakename[1]
# code <- mylake$lakecode[1]
# mylake_ssa <- f_icessa(code)
# 
# pdf("ahtarinjarvi_eigenvalues.pdf")
# print(plot(mylake_ssa, type = "values", main = paste0(name, " ", code, " Eigenvalues")))
# dev.off()
# pdf("ahtarinjarvi_eigenvectors.pdf")
# print(plot(mylake_ssa, type = "vectors", main = paste0(name, " ", code, " Eigenvectors")))
# dev.off()
# pdf("ahtarinjarvi_paired_eigenvectors.pdf")
# print(plot(mylake_ssa, type = "paired", main = paste0(name, " ", code, " Paired Eigenvectors")))
# dev.off()
# pdf("ahtarinjarvi_wcor.pdf")
# print(plot(wcor(mylake_ssa), main = paste0(name, " ", code, " W-Correlation Matrix")))
# dev.off()

#> input: lakecode, x
#> does: apply f_icessa, reconstruct data, take residuals, log transform
#> output: (x=1) empirical residuals, (x=2) reconstruct object, (x=3) logit transformed residuals, (x=4) bounds
f_iceres = function(code, x) {
  log_ssa = f_icessa(code)
  lake_row = ssa_lakes %>%
    filter(lakecode == code)
  group = lake_row$groups[[1]]
  icerec = reconstruct(log_ssa, groups = group)
  iceresiduals = residuals(icerec)
  ires = iceresiduals
  sqrtabsires = sqrt(abs(ires))
  a2 = min(sqrtabsires) * 0.9 # these two values
  b2 = max(sqrtabsires) * 1.1
  bounds = c(a2, b2)
  res_trans = log((sqrtabsires - a2)/(b2 - sqrtabsires))
  if (x == 1) {
    return(ires)
  } else if (x == 2) {
    return(icerec)
  } else if (x == 3) {
    return(res_trans)
  } else if (x == 4) {
    return(bounds)
  }
}

#> input: lakecode
#> returns reconstructed series plots for lakecode data
f_reconplot = function(code) {
  ssa_ice = f_icessa(code)
  recon = f_iceres(code, 2)
  name = f_lake_data(code)[1,6]
  print(plot(recon, main = paste0(name, " ", code, ": Reconstructed Series")))
}

#> input: lakecode
#> pdf plot of reconstructed series for lakecode data
f_reconpdf <- function(code) {
  ssa_ice <- f_icessa(code)
  recon <- f_iceres(code, 2)
  name <- f_lake_data(code)[1,6]
  pdf(paste0(name, " recon.pdf"))
  print(plot(recon, main = paste0(name, " ", code, ": Reconstructed Series")))
  dev.off()
}

#> input: lakecode
#> output: df matrix of trend series
yeardur = function(code) {
  alldata = f_lake_data(code)
  recon_obj = f_iceres(code, 2)
  years = alldata$iceoff_year
  years = sort(years)
  country_vec = alldata$country
  code_vec = alldata$lakecode
  name_vec = alldata$lakename
  # df = data.frame(matrix(nrow = length(recon_obj$F1), ncol = 3))
  df = data.frame(year = years, recon = recon_obj$F1, country = country_vec, lakename = name_vec, lakecode = code_vec)
  return(df)
}

#> OVERWRITTEN
#> input: lakecode
#> output: df matrix of trend series, INCLUDES RESIDUALS
yeardur2 = function(code) {
  alldata = f_lake_data(code)
  recon_obj = f_iceres(code, 2)
  recon_obj_res = resssatheta(code)
  years = alldata$iceoff_year
  years = sort(years)
  country_vec = alldata$country
  code_vec = alldata$lakecode
  name_vec = alldata$lakename
  # df = data.frame(matrix(nrow = length(recon_obj$F1), ncol = 3))
  df = data.frame(year = years, recon = recon_obj$F1, reconres = recon_obj_res$F1, country = country_vec, lakename = name_vec, lakecode = code_vec)
  return(df)
}

#> input: lakecode
#> output: df matrix of trend series, INCLUDES RESIDUALS
yeardur2 = function(code) {
  alldata = f_lake_data(code)
  recon_obj = f_iceres(code, 2)
  years = alldata$iceoff_year
  years = sort(years)
  country_vec = alldata$country
  code_vec = alldata$lakecode
  name_vec = alldata$lakename
  # df = data.frame(matrix(nrow = length(recon_obj$F1), ncol = 3))
  df = data.frame(year = years, recon = recon_obj$F1, reconres = recon_obj[length(recon_obj)-1][[1]], country = country_vec, lakename = name_vec, lakecode = code_vec)
  return(df)
}

#> input: reconstructed trend df (shown later)
#> output: reconstructed trend df with missing years added back in
addNA = function(df) {
  min = min(df$year)
  max = max(df$year)
  country = df$country[1]
  name = df$lakename[1]
  code = df$lakecode[1]
  nadf = data.frame(matrix(data = NA, nrow = max-min+1, ncol = 2))
  nadf[,1] = min:max
  for (i in min:max) {
    for (j in 1:nrow(df)) {
      if (df[j,1] == i) {
        nadf[i-min+1,2] = df[j,2]
      }
    }
  }
  n = nrow(nadf)
  country_vec = rep(country, n)
  name_vec = rep(name, n)
  code_vec = rep(code, n)
  nadf = cbind(nadf, country_vec, name_vec, code_vec)
  nadf = nadf %>%
    rename(year = X1,
           recon = X2,
           country = country_vec,
           lakename = name_vec,
           lakecode = code_vec)
  return(nadf)
}


#> input: reconstructed trend df (shown later) (INCLUDES RESIDUALS)
#> output: reconstructed trend df with missing years added back in
addNA2 = function(df) {
  min = min(df$year)
  max = max(df$year)
  country = df$country[1]
  name = df$lakename[1]
  code = df$lakecode[1]
  nadf = data.frame(matrix(data = NA, nrow = max-min+1, ncol = 2))
  nadf[,1] = min:max
  for (i in min:max) {
    for (j in 1:nrow(df)) {
      if (df[j,1] == i) {
        nadf[i-min+1,2] = df[j,2]
        nadf[i-min+1,3] = df[j,3]
      }
    }
  }
  n = nrow(nadf)
  country_vec = rep(country, n)
  name_vec = rep(name, n)
  code_vec = rep(code, n)
  nadf = cbind(nadf, country_vec, name_vec, code_vec)
  nadf = nadf %>%
    rename(year = X1,
           recon = X2,
           reconres = V3,
           country = country_vec,
           lakename = name_vec,
           lakecode = code_vec)
  return(nadf)
}

# retrieve coords: output (name, country, lat, long)
#> input: lakecode
#> output: name, country, lat, long
coords = function(code) {
  data = f_lake_data(code)
  spec_coord = c()
  spec_coord[1] = data$lakename[1]
  spec_coord[2] = data$country[1]
  spec_coord[3] = data$latitude[1]
  spec_coord[4] = data$longitude[1]
  return(spec_coord)
}

#> input: lakecode
#> output: coordinates of lake: c(lat, long)
justcoords = function(code) {
  data = f_lake_data(code)
  spec_coord = c()
  spec_coord[1] = data$latitude[1]
  spec_coord[2] = data$longitude[1]
  return(spec_coord)
}

#> input: lakecode
#> output: lat index for tb (brightness temperature) grid
latindex = function(code) {
  lati = justcoords(code)[1]
  for (i in 1:length(lat)) {
    lower = lat[i] - 1.25
    upper = lat[i] + 1.25
    if (lati >= lower & lati < upper) {
      return(i)
    }
  }
}

#> input: lakecode
#> output: long index for tb grid
longindex = function(code) {
  longi = justcoords(code)[2]
  for (i in 1:length(long)) {
    lower = long[i] - 1.25
    upper = long[i] + 1.25
    if (longi >= lower & longi < upper) {
      return(i)
    }
  }
}

#> input: lakecode
#> output: tibble with monthly tb values at lake grid location
f_laketemp = function(code) {
  latc = latindex(code)
  longc = longindex(code)
  yeartemp = c()
  yearlytemp = as.data.frame(matrix(ncol = 12))
  colnames(yearlytemp) = c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December")
  for (i in 1:45) {
    for (j in 1:12) {
      dim = 12 * (i-1) + j + 1
      yeartemp[j] = tb[longc, latc, dim]
    }
    yearlytemp[i,] = yeartemp
    yeartemp = c()
  }
  df = tibble(
    year = 1978:2022,
    yearlytemp
  )
  return(df)
}

#> input: lakecode
#> output: append tb data to icd data for lakecode. calculate simple 
#> summary statistics for tb data (avg, sd, min, max, range)
f_specdattemp2 = function(code) {
  temps = f_laketemp(code)
  dat = ice2 %>%
    filter(lakecode == code) %>%
    arrange(iceoff_year) %>%
    dplyr::select(-c(iceon_year, lakeorriver, longitude)) %>%
    left_join(temps, by = join_by(iceoff_year == year)) %>%
    mutate(avg_bt = NA,
           sd_bt = NA,
           min_bt = NA,
           max_bt = NA)
  for (j in 1:nrow(dat)) {
    dat$avg_bt[j] = mean(as.numeric(dat[j, 7:18]))
    dat$sd_bt[j] = sd(dat[j, 7:18])
    dat$min_bt[j] = min(dat[j, 7:18])
    dat$max_bt[j] = max(dat[j, 7:18])
  }
  dat = dat %>%
    mutate(range_bt = max_bt - min_bt)
  return(dat)
}

#> input: lakecode
#> does: apply ssa to logit transformed residuals
resssa2 = function(code) {
  restrans = f_iceres(code, 3)
  ssares = ssa(restrans, L = ceiling(length(restrans)/2))
  return(ssares)
}

#> input: lakecode
#> output: print ssa plots for logit transformed residuals
resssa2plot = function(code) {
  specific_lake = f_lake_data(code)
  name = specific_lake$lakename[1]
  restrans = f_iceres(code, 3)
  ssares = ssa(restrans, L = ceiling(length(restrans)/2))
  print(plot(ssares, main = paste0(name, " Transformed Residual Eigenvalues")))
  print(plot(ssares, type = "vectors", main = paste0(name, " Transformed Residual Eigenvectors")))
  print(plot(ssares, type = "paired", main = paste0(name, " Transformed Residual Paired Eigenvectors")))
  print(plot(wcor(ssares), main = paste0(name, " Transformed Residual W-Correlation Matrix")))
}

#> input: lakecode
#> output: individual pdfs of ssa plots for logit transformed residuals
resssa2pdf = function(code) {
  specific_lake = f_lake_data(code)
  name = specific_lake$lakename[1]
  restrans = f_iceres(code, 3)
  ssares = ssa(restrans, L = ceiling(length(restrans)/2))
  pdf(paste0(name, " resssa values.pdf"))
  print(plot(ssares, main = paste0(name, " Transformed Residual Singular Values")))
  dev.off()
  pdf(paste0(name, " resssa vectors.pdf"))
  print(plot(ssares, type = "vectors", main = paste0(name, " Transformed Residual Singular Vectors")))
  dev.off()
  pdf(paste0(name, " resssa paired.pdf"))
  print(plot(ssares, type = "paired", main = paste0(name, " Transformed Residual Paired Singular Vectors")))
  dev.off()
  pdf(paste0(name, " resssa wcor.pdf"))
  print(plot(wcor(ssares), main = paste0(name, " Transformed Residual W-Correlation Matrix")))
  dev.off()
}

#> input: lakecode
#> output: reconstruct object for logit transformed residuals
resssatheta = function(code) {
  ssares = resssa2(code)
  lake_row = ssa_lakes %>%
    filter(lakecode == code)
  group = lake_row$resgroups[[1]]
  recon = reconstruct(ssares, groups = group)
  return(recon)
}

#> input: lakecode
#> output: pdf of reconstruct object for logit transformed residuals
resssathetapdf = function(code) {
  specific_lake = f_lake_data(code)
  name = specific_lake$lakename[1]
  ssares = resssa2(code)
  lake_row = ssa_lakes %>%
    filter(lakecode == code)
  group = lake_row$resgroups[[1]]
  recon = reconstruct(ssares, groups = group)
  pdf(paste0(name, " resrecon.pdf"))
  print(plot(recon, main = paste0(name, " ", code, ": Residuals Reconstructed Series")))
  dev.off()
}

#> input: lakecode
#> output: undo logit transform
expitsq = function(code) {
  df <- addNA2(yeardur2(code))
  bds <- f_iceres(code, 4)
  l = bds[1]
  u = bds[2]
  for (i in 1:nrow(df)) {
    # test
    # df[i,2] = 1.4826*((u*exp(df[i,2]) + l)/(1+exp(df[i,2])))^2
    df[i,3] = 1.4826*((u*exp(abs(df[i,3])) + l)/(1+exp(abs(df[i,3]))))^2
  }
  return(df)
}

#> input: lakecode
#> output: t-statistics
tstat = function(code) {
  dat <- expitsq(code)
  # dat <- dat %>%
  #   mutate(recon = recon - mean(recon, na.rm=TRUE))
  # first_recon <- dat[which.min(dat$year),]$recon
  # dat <- dat %>%
  #   mutate(recon = recon - first_recon)
  all_recon <- dat %>%
    filter(year == 1991)
  all_recon <- all_recon$recon
  dat <- dat %>%
    mutate(recon = recon - all_recon)
  # df = data.frame(matrix(data = NA, nrow = nrow(mu), ncol = 2))
  # df[,1] = mu[,1]
  # for (i in 1:nrow(mu)) {
  #   if (is.na(mu[i,2]) == FALSE) {
  #     df[i,2] = mu[i,2]/sigma[i,2]
  #   }
  # }
  tstats <- c()
  for (i in 1:nrow(dat)) {
    if (is.na(dat$recon[i]) == FALSE) {
      # dat$recon[i] - (first lvl)
      tstats <- c(tstats, (dat$recon[i]) / dat$reconres[i])
    } else {
      tstats <- c(tstats, NA)
    }
  }
  dat <- cbind(dat, tstats)
  return(dat)
}

#> input: lakecode
#> output: monthly tb for grid containing lakecode
gettempemp = function(code) {
  latc = latindex(code)
  longc = longindex(code)
  temps = c()
  yeartemp = c()
  for (i in 1:45) {
    for (j in 1:12) {
      dim = 12 * (i-1) + j + 1
      yeartemp[j] = tb[longc, latc, dim]
    }
    temps = append(temps, yeartemp)
    yeartemp = c()
  }
  return(temps)
}

#> input: lakecode
#> output: yearly max tb for grid containing lakecode
getmaxtemp = function(code) {
  latc = latindex(code)
  longc = longindex(code)
  temps = c()
  yeartemp = c()
  for (i in 1:45) {
    for (j in 1:12) {
      dim = 12 * (i-1) + j + 1
      yeartemp[j] = tb[longc, latc, dim]
    }
    temps[i] = max(yeartemp)
    yeartemp = c()
  }
  return(temps)
}

#> input: lakecode
#> output: year and ice coverage duration (icd) of lakecode
durs = function(code) {
  data = f_lake_data(code)
  year = data$iceoff_year
  #year = sort(year)
  dur = data$duration
  durations = c()
  for (i in 1:length(year)) {
    if (year[i] > 1978) {
      durations = append(durations, dur[i])
    }
  }
  n = length(durations)
  m = 45-n
  durations = append(durations, rep(NA,m))
  return(durations)
}

#> input: lakecode
#> output: yearly average tb for grid containing lakecode
gettemp = function(code) {
  latc = latindex(code)
  longc = longindex(code)
  temps = c()
  yeartemp = c()
  for (i in 1:45) {
    for (j in 1:12) {
      dim = 12 * (i-1) + j + 1
      yeartemp[j] = tb[longc, latc, dim]
    }
    temps[i] = mean(yeartemp)
    yeartemp = c()
  }
  return(temps)
}

#> input: vector of lakecodes
#> output: plot of average tb vs icd
scattertemp = function(codes) {
  plot(gettemp(codes[1]), durs(codes[1]), main="Brightness Temperature vs Ice Cover Duration",
       xlab = "Brightness Temperature", ylab = "Ice Cover Duration", xlim = c(255.5,269.5),ylim=c(0,250))
  # legend(265,250, legend=c("SA < 20", "20 < SA < 40", "40 < SA < 100", "SA > 100"), 
  #        col=c("blue", "green", "orange","red"),lty=1, cex=0.8)
  for (i in 2:length(codes)) {
    # sa = areas[indices[i]]
    # if (sa < 20) {
    #   color = "blue"
    # } else if (sa >= 20 & sa < 40) {
    #   color = "green"
    # } else if (sa >= 40 & sa < 100) {
    #   color = "orange"
    # } else if (sa >= 100) {
    #   color = "red"
    # }    
    points(gettemp(codes[i]), durs(codes[i]))
  }
}

#> input: vector of lakecodes
#> output: ds of durations and avg tb
tbicd = function(codes) {
  n = length(codes)
  vec1 = gettemp(codes[1])
  vec2 = durs(codes[1])
  ds = cbind(vec1, vec2)
  for (i in 2:n) {
    tempvec = gettemp(codes[i])
    dursvec = durs(codes[i])
    ds_pass = cbind(tempvec, dursvec)
    ds = rbind(ds, ds_pass)
  }
  return(ds)
}

#> input: vector of lakecodes
#> output: df of avg tb, duration, latitude
tbicdlat = function(codes) {
  n = length(codes)
  vec1 = gettemp(codes[1])
  vec2 = durs(codes[1])
  len_vec3 = length(vec1)
  lake_data = coords(codes[1])
  lat_val = lake_data[3]
  vec3 = rep(lat_val, len_vec3)
  ds = cbind(vec1, vec2, vec3)
  for (i in 2:n) {
    tempvec = gettemp(codes[i])
    dursvec = durs(codes[i])
    len_lat_vec = length(tempvec)
    lake_data = coords(codes[i])
    lat_val = lake_data[3]
    lat_vec = rep(lat_val, len_lat_vec)
    ds_more = cbind(tempvec, dursvec, lat_vec)
    ds = rbind(ds, ds_more)
  }
  return(ds)
}

#> input: vector of lakecodes, name of plot
#> output: plot of avg tb vs duration with color as surface area
scattertempsize = function(codes, name) {
  plot(gettemp(codes[1]), durs(codes[1]), main=name, col="green", 
       xlab = "Brightness Temperature", ylab = "Ice Cover Duration", 
       xlim = c(255.5,269.5),ylim=c(0,250))
  legend(265,250, legend=c("SArea < 20", "20 < SArea < 40", "40 < SArea < 100", "SArea > 100"),
         col=c("blue", "green", "orange","red"),lty=1, cex=0.8)
  for (i in 2:length(codes)) {
    dscode <- areas %>%
      filter(lakecode == codes[i])
    sarea <- dscode$surface_area
    if (is.na(sarea) == TRUE) {
      #> do nothing
    } else {
      if (sarea < 20) {
        color = "blue"
      } else if (sarea >= 20 & sarea < 40) {
        color = "green"
      } else if (sarea >= 40 & sarea < 100) {
        color = "orange"
      } else if (sarea >= 100) {
        color = "red"
      }
      points(gettemp(codes[i]), durs(codes[i]), col=color)
    }
  }
}

#> input: vector of lakecodes, name of plot
#> output: plot of max tb vs duration with color as surface area
scattertempsizemax = function(codes, name) {
  plot(getmaxtemp(codes[1]), durs(codes[1]), main=name, col = "green",
       xlab = "Brightness Temperature", ylab = "Ice Cover Duration", 
       xlim = c(255.5,273),ylim=c(0,250))
  legend(265,250, legend=c("SA < 20", "20 < SA < 40", "40 < SA < 100", "SA > 100"),
         col=c("blue", "green", "orange","red"),lty=1, cex=0.8)
  for (i in 2:length(codes)) {
    dscode <- areas %>%
      filter(lakecode == codes[i])
    sarea <- dscode$surface_area
    if (is.na(sarea) == TRUE) {
      #> do nothing
    } else {
      if (sarea < 20) {
        color = "blue"
      } else if (sarea >= 20 & sarea < 40) {
        color = "green"
      } else if (sarea >= 40 & sarea < 100) {
        color = "orange"
      } else if (sarea >= 100) {
        color = "red"
      }
      points(getmaxtemp(codes[i]), durs(codes[i]), col=color)
    }
  }
}

#> input: vector of lakecodes, name of plot
#> output: plot of avg tb vs duration with color as surface area
scattertempelev = function(codes, name) {
  plot(gettemp(codes[1]), durs(codes[1]), main=name, col="yellow3", 
       xlab = "Brightness Temperature", ylab = "Ice Cover Duration", 
       xlim = c(255.5,269.5),ylim=c(0,250))
  legend(265,250, legend=c("Elevation < 100", "100 < Elevation < 200", "200 < Elevation < 300", "300 < Elevation < 400", "Elevation > 400"),
         col=c("blue", "green", "yellow3", "orange","red"),lty=1, cex=0.8)
  for (i in 2:length(codes)) {
    dscode <- elev %>%
      filter(lakecode == codes[i])
    thiselev <- dscode$elevation
    if (is.na(thiselev) == TRUE) {
      #> do nothing
    } else {
      if (thiselev < 100) {
        color = "blue"
      } else if (thiselev >= 100 & thiselev < 200) {
        color = "green"
      }else if (thiselev >= 200 & thiselev < 300) {
        color = "yellow3"
      } else if (thiselev >= 300 & thiselev < 400) {
        color = "orange"
      } else if (thiselev >= 400) {
        color = "red"
      }
      points(gettemp(codes[i]), durs(codes[i]), col=color)
    }
  }
}

#> input: endpoints
#> output: linear intercept
intercept = function(a, b) {
  if (a == b) {
    return(1)
  }
  b = b - 263
  a = a - 263
  slope = b - a
  x = -a / slope
  return(x)
}

#> function deprecated
# getyeartb = function(code,i) {
#   latc = latindex(code)
#   longc = longindex(code)
#   yeartemp = c()
#   for (j in 1:12) {
#     dim = 12 * (i-1) + j + 1
#     yeartemp[j] = tb[longc, latc, dim]
#   }
#   return(yeartemp)
# }

#> input: lakecode
#> output: df of centered duration above 263 tb
getpass3 = function(code) {
  latc = latindex(code)
  longc = longindex(code)
  yeartemp = c()
  dur = c()
  for (i in 1:45) {
    month = c()
    months = c()
    ints = c()
    for (j in 1:12) {
      dim = 12 * (i-1) + j + 1
      yeartemp[j] = tb[longc, latc, dim]
    }
    for (k in 2:12) {
      if (yeartemp[k] > 263 & yeartemp[k-1] > 263) {
        month[k] = NA
      }
      if (yeartemp[k] < 263 & yeartemp[k-1] < 263) {
        month[k] = NA
      }
      if (yeartemp[k] < 263 & yeartemp[k-1] > 263) {
        month[k] = k
      }
      if (yeartemp[k] > 263 & yeartemp[k-1] < 263) {
        month[k] = k
      }
    }
    months = na.omit(month)
    for (m in 1:length(months)) {
      ints[m] = intercept(yeartemp[months[m]-1], yeartemp[months[m]])
      ints[m] = ints[m] + months[m] - 1
    }
    n = length(ints)
    for (p in 1:n) {
      if (n == 1) {
        dur[i] = ints[1]
      }
      if (n == 2) {
        dur[i] = ints[1] + 12 - ints[2]
      }
      if (n == 3) {
        dur[i] = ints[1] + ints[3] - ints[2]
      }
      if (n == 4) {
        dur[i] = ints[1] + ints[3] - ints[2] + 12 - ints[4]
      }
      if (n == 5) {
        dur[i] = ints[1] + ints[3] - ints[2] + ints[5] - ints[4]
      }
      if (n == 6) {
        dur[i] = ints[1] + ints[3] - ints[2] + ints[5] - ints[4] + 12 - ints[6]
      }
    }
  }
  dur20 = dur[1:20]
  meandur20 = mean(dur20)
  meandur2045 = rep(meandur20,45)
  centraldur = dur - meandur2045
  return(centraldur)
}

#> deprecated
# freezetempdur = function(monthup, monthdn) {
#   durf = c()
#   for (i in 2:45) {
#     durf = monthup[i] + monthdn[i-1]
#   }
#   return(durf)
# }


#> input: vector of lakecodes
#> output: temps and durations
tbicd = function(codes) {
  n = length(codes)
  vec1 = gettemp(codes[1])
  vec2 = durs(codes[1])
  ds = cbind(vec1, vec2)
  for (i in 2:n) {
    tempvec = gettemp(codes[i])
    dursvec = durs(codes[i])
    ds_pass = cbind(tempvec, dursvec)
    ds = rbind(ds, ds_pass)
  }
  return(ds)
}

#> input: vector of lakecodes
#> output: temps, durations, and latitudes
tbicdlat = function(codes) {
  n = length(codes)
  vec1 = gettemp(codes[1])
  vec2 = durs(codes[1])
  lakename <- f_lake_data(codes[1])$lakename[1]
  len_vec3 = length(vec1)
  lake_data = coords(codes[1])
  lat_val = lake_data[3]
  vec3 = rep(lat_val, len_vec3)
  vec4 <- rep(lakename, len_vec3)
  ds = cbind(vec1, vec2, vec3, vec4)
  for (i in 2:n) {
    tempvec = gettemp(codes[i])
    dursvec = durs(codes[i])
    lakename <- f_lake_data(codes[i])$lakename[1]
    len_lat_vec = length(tempvec)
    lake_data = coords(codes[i])
    lat_val = lake_data[3]
    lat_vec = rep(lat_val, len_lat_vec)
    name_vec <- rep(lakename, len_lat_vec)
    ds_more = cbind(tempvec, dursvec, lat_vec, name_vec)
    ds = rbind(ds, ds_more)
  }
  return(ds)
}

#> input: vector of lakecodes
#> output: dataframe of lakecodes and corresponding first periods
f_period <- function(codes) {
  df <- data.frame(lakecode=codes, period1=rep(0,length(codes)))
  for (i in codes) {
    ssa_obj <- f_icessa(i)
    # print("ssa done")
    groups_obj <- ssa_lakes$groups[which(ssa_lakes$lakecode==i)][[1]]
    # print("groups done")
    trendstat <- ssa_lakes %>%
      filter(lakecode == i) %>%
      dplyr::select(notrend)
    # print("trend done")
    per_obj <- parestimate(ssa_obj, groups=groups_obj, method="esprit")
    # print("par done")
    if (trendstat == FALSE) {
      df[which(df$lakecode==i),'period1'] <- per_obj$F2$periods[1]
    } else {
      df[which(df$lakecode==i),'period1'] <- per_obj$F1$periods[1]
    }
  }
  return(df)
}

#> input: vector of lakecodes
#> output: dataframe of lakecodes and corresponding minimum period
f_minperiod <- function(codes) {
  df <- data.frame(lakecode=codes, minperiod=rep(0,length(codes)))
  fvec <- c('F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8')
  for (i in codes) {
    ssa_obj <- f_icessa(i)
    # print("ssa done")
    groups_obj <- ssa_lakes$groups[which(ssa_lakes$lakecode==i)][[1]]
    # print("groups done")
    trendstat <- ssa_lakes %>%
      filter(lakecode == i) %>%
      dplyr::select(notrend)
    # print('trendstat done')
    per_obj <- parestimate(ssa_obj, groups=groups_obj, method="esprit")
    # print(per_obj)
    # return(per_obj)
    # print('period done')
    # print(length(per_obj))
    per_vec <- c()
    for (j in 1:length(per_obj)) {
      per_vec <- c(per_vec, per_obj[[fvec[j]]][['periods']][1])
    }
    # print(per_vec)
    # print(!is.infinite(per_vec))
    per_vec <- per_vec[!is.infinite(per_vec)]
    df[which(df$lakecode==i),'minperiod'] <- min(per_vec)
  }
  return(df)
}

#> input: vector of lakecodes
#> output: df of lakecodes and all periods
f_allperiod <- function(codes) {
  df <- data.frame(lakecode=character(), period=numeric())
  fvec <- c('F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8')
  for (i in codes) {
    ssa_obj <- f_icessa(i)
    # print("ssa done")
    groups_obj <- ssa_lakes$groups[which(ssa_lakes$lakecode==i)][[1]]
    # print("groups done")
    # trendstat <- ssa_lakes %>%
    #   filter(lakecode == i) %>%
    #   dplyr::select(notrend)
    # print('trendstat done')
    per_obj <- parestimate(ssa_obj, groups=groups_obj, method="esprit")
    # print(per_obj)
    # return(per_obj)
    # print('period done')
    # print(length(per_obj))
    per_vec <- c()
    for (j in 1:length(per_obj)) {
      per_vec <- c(per_vec, per_obj[[fvec[j]]][['periods']][1])
    }
    # print(per_vec)
    # print(!is.infinite(per_vec))
    per_vec <- per_vec[!is.infinite(per_vec)]
    df_spec <- data.frame(lakecode=rep(i,length(per_vec)),period=per_vec)
    df <- rbind(df, df_spec)
  }
  return(df)
}

#################################################################################################################

#> ssa not suitable for these lakes, removed from data set
lake_removal = c("JJM18", "JJM3", "JJM4", "JJM7", "JJM8", "KMS14", "PJD7")

ssa_lakes = data.frame(more_lakes_35)

ssa_lakes = ssa_lakes %>%
  rename(lakecode = more_lakes_35) %>%
  filter(!lakecode %in% lake_removal) 

#> singular value grouping
groups = list(
  # mendota DMR1
  list(c(1), c(2,3), c(4,5), c(6,7)),
  # monona DMR2
  list(c(1,2), c(4,5)),
  # orsasjon GW240
  list(c(1), c(3,4)),
  # nackten GW341
  list(c(1,4), c(2,3), c(5,6)),
  # jukkasjarvi GW592
  list(c(1,8), c(2,3), c(4,5), c(6,7)),
  # allequash JJM1
  list(c(1), c(3,4)),
  # escanaba JJM15
  list(c(1,4,8), c(2,3), c(6,7)),
  # nebish JJM17
  list(c(1,5), c(2,3)),
  # muskellunge JJM2
  list(c(1), c(2,3), c(4,5), c(6,7)),
  # shell JJM22
  list(c(5), c(1,2), c(3,4)), 
  # spruce JJM23
  list(c(1), c(2,3)),
  # big green JJM24
  list(c(1,2), c(3,4), c(5,6), c(7,8), c(9,10)), # no trend
  # geneva JJM27
  list(c(3,6), c(1,2), c(4,5)),
  # superior JJM33
  list(c(2,3,6), c(4,5), c(7,8)),
  # sparkling JJM6
  list(c(1), c(2,3), c(4,5), c(6,7)),
  # oulujarvi JK01
  list(c(1), c(4,5)),
  # kallavesi JK02
  list(c(1,2), c(3,4)),
  # nasijarvi JK03
  list(c(1), c(2,3)),
  # vesijarvi JK04
  list(c(1,2), c(3,4), c(6,7)),
  # paijanne JK05
  list(c(2), c(6,7)),
  # muurasjarvi JK06
  list(c(1,2), c(3,4), c(5,6)),
  # kalmarinjarvi JK07
  list(c(1,2), c(4,5)),
  # summasjarvi JK08
  list(c(1,2), c(3,4), c(6,7)),
  # pielavesi savia JK09
  list(c(1,4), c(2,3), c(5,6)),
  # hankavesi JK10
  list(c(1,2), c(3,4), c(5,6)),
  # lappajarvi JK12
  list(c(1,2), c(3,4), c(5,6), c(7,8)),
  # kukkia JK14
  list(c(1,2), c(3,4)),
  # ala rieveli JK17
  list(c(1,2), c(4,5)),
  # vesijarvi JK18
  list(c(1,2), c(4,5)),
  # jaasjarvi JK19
  list(c(1,2), c(3,4)),
  # iisvesi JK21
  list(c(1,2), c(3,4), c(5,6), c(7,8), c(9,10)),
  # paajarvi JK22
  list(c(1,2), c(3,4), c(7,8)),
  # kivijarvi JK23
  list(c(1), c(2,3), c(4,5)),
  # saanijarvi JK24
  list(c(1,2), c(4,5)),
  # haukivesi JK25
  list(c(1), c(2,3), c(5,6), c(7,8)),
  # inari nellim JK26
  list(c(3,4), c(1,2)),
  # kilipsjarvi JK28
  list(c(1), c(2,3)),
  # lentua JK31
  list(c(1,2), c(3,4), c(5,6)),
  # lestijarvi JK32
  list(c(1), c(2,3), c(4,5)),
  # mutusjarvi JK34
  list(c(1), c(2,3), c(5,6)),
  # oijarvi JK36
  list(c(1), c(7,8)),
  # ounasjarvi JK37
  list(c(1), c(2,3), c(5,6)),
  # pielinen JK40
  list(c(1,2), c(3,4), c(5,6)),
  # rehja JK42
  list(c(1,2), c(3,4), c(5,6)),
  # visuvesi JK45
  list(c(1,4), c(2,3)),
  # palovesi JK46
  list(c(1,2), c(3,4)),
  # simpelejarvi JK47
  list(c(1,2), c(4,5), c(6,7)),
  # ahtarinjarvi JK48
  list(c(1,2), c(3,4)),
  # kuivajarvi JK49
  list(c(1,2), c(3,4)),
  # george KMS10
  list(c(1), c(2,3), c(4,5), c(6,7)),
  # cazenovia LR2
  list(c(1,2)), # no trend
  # gull MICH03
  list(c(3,6), c(1,2), c(4,5), c(7,8)), # check
  # fair MICH06
  list(c(3,4), c(1,2), c(7,8)),
  # baikal NG1
  list(c(1,2), c(3,4), c(5,6), c(8,9)),
  # mueggelesee RA1
  list(c(1), c(2,3))
)

#> lakes have no trend component
lake_code_notrend = c("JJM24", "LR2")

#> add groups to dataframe
newcol = rep(NA,nrow(ssa_lakes))
ssa_lakes = cbind(ssa_lakes, data.frame(newcol))
ssa_lakes = ssa_lakes %>%
  rename(groups = newcol)


for (i in 1:nrow(ssa_lakes)) {
  ssa_lakes$groups[i] = groups[i]
}

ssa_lakes = ssa_lakes %>%
  mutate(notrend = lakecode %in% lake_code_notrend)

#> reconstruct plots for ssa
# pdf("ssa_recon_35.pdf")
for (i in ssa_lakes$lakecode) {
  f_reconplot(i)
}
# dev.off()

#> reformat data
allrecon = tibble(year = integer(), recon = double(), country = character(), lakename = character(), lakecode = character())
for (i in 1:nrow(ssa_lakes)) {
  if (ssa_lakes$notrend[i] == FALSE) {
    newdat = addNA(yeardur(ssa_lakes$lakecode[i]))
    allrecon = allrecon %>%
      bind_rows(newdat)
  }
}

#> plot of reconstructed trends: nordic
# pdf("finland_sweden_trend_plot.pdf")
allrecon %>%
  filter(country %in% c("SWEDEN", "FINLAND")) %>%
  ggplot(aes(x = year, y = recon)) +
  geom_line(aes(group = lakecode, color = country)) +
  scale_color_manual(values = c("#56B4E9", "#56B4E9")) +
  theme_bw() + 
  labs(# title = "Reconstructed Lake Ice Cover Duration Trends",
       title = "Finland and Sweden",
       y = "Reconstructed Ice Coverage Duration",
       x = "Year") + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")
# dev.off()

#> plot of reconstructed trends: !nordic
# pdf("usa_russia_germany_trend_plot.pdf")
allrecon %>%
  filter(country %in% c("UNITED STATES", "RUSSIA", "GERMANY")) %>%
  ggplot(aes(x = year, y = recon)) +
  geom_line(aes(group = lakecode, color = country)) + 
  scale_color_manual(values = c("#E69F00", "#CC79A7", "#000000")) +
  theme_bw() + 
  labs(title = "Germany (in orange), Russia (in purple), and North America (in black)",
       y = "Reconstructed Ice Coverage Duration",
       x = "Year") + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")
# dev.off()

#> append tb data to dataframe
ds_final = f_specdattemp2(ssa_lakes$lakecode[1])
for (i in 2:nrow(ssa_lakes)) {
  ds_final = rbind(ds_final, f_specdattemp2(ssa_lakes$lakecode[i]))
}
ds_final <- ds_final %>%
  left_join(ssa_lakes)

# saveRDS(ds_final, file = "tb_icd_groups.rds")

ds_nogroup <- ds_final %>%
  dplyr::select(-c(groups, notrend))

# write_csv(ds_nogroup, file="tb_icd_nogroup.csv")

# pdf("resssa2plot.pdf")
# for (i in ssa_lakes$lakecode) {
#   resssa2plot(i)
# }
# dev.off()

#> singular value groups for residual ssa
resgroups = list(
  # mendota DMR1
  list(c(1), c(2,3), c(4,5)),
  # monona DMR2
  list(c(1,2), c(3,4)),
  # orsasjon GW240
  list(c(1), c(2,3)),
  # nackten GW341
  list(c(1,6), c(2,3), c(4,5)),
  # jukkasjarvi GW592
  list(c(3,4), c(1,2)),
  # allequash JJM1
  list(c(1,2), c(3,4)),
  # escanaba JJM15
  list(c(5), c(1,2), c(3,4)),
  # nebish JJM17
  list(c(4), c(1,2)),
  # muskellunge JJM2
  list(c(1,2), c(3,4)),
  # shell JJM22
  list(c(7), c(1,2), c(3,4), c(5,6)), 
  # spruce JJM23
  list(c(7), c(1,2), c(3,4), c(5,6)),
  # big green JJM24
  list(c(1), c(2,3), c(4,5)),
  # geneva JJM27
  list(c(5,6,7), c(1,2), c(3,4)),
  # superior JJM33
  list(c(1,6), c(2,3), c(4,5)),
  # sparkling JJM6
  list(c(3,4), c(1,2), c(5,6)),
  # oulujarvi JK01
  list(c(3), c(1,2), c(4,5), c(6,7)),
  # kallavesi JK02
  list(c(3,4), c(1,2)),
  # nasijarvi JK03
  list(c(1,2), c(3,4)),
  # vesijarvi JK04
  list(c(1,4), c(2,3)),
  # paijanne JK05
  list(c(1), c(2,3), c(4,5)),
  # muurasjarvi JK06
  list(c(1), c(2,3), c(4,5)),
  # kalmarinjarvi JK07
  list(c(1,2,3), c(4,5)),
  # summasjarvi JK08
  list(c(1), c(2,3), c(6,7)),
  # pielavesi savia JK09
  list(c(3,4), c(1,2), c(5,6)),
  # hankavesi JK10
  list(c(1), c(2,3), c(4,5)),
  # lappajarvi JK12
  list(c(1), c(2,3)),
  # kukkia JK14
  list(c(1), c(2,3)),
  # ala rieveli JK17
  list(c(1), c(2,3)),
  # vesijarvi JK18
  list(c(6), c(1,2), c(3,4)),
  # jaasjarvi JK19
  list(c(5), c(1,2), c(3,4)),
  # iisvesi JK21
  list(c(1), c(2,3), c(4,5)),
  # paajarvi JK22
  list(c(1), c(2,3)),
  # kivijarvi JK23
  list(c(3), c(1,2)),
  # saanijarvi JK24
  list(c(3,4), c(1,2)),
  # haukivesi JK25
  list(c(1), c(2,3), c(4,5), c(6,7)),
  # inari nellim JK26
  list(c(1,2), c(3,4)),
  # kilipsjarvi JK28
  list(c(1), c(2,3)),
  # lentua JK31
  list(c(1), c(2,3)),
  # lestijarvi JK32
  list(c(7), c(1,2), c(3,4)),
  # mutusjarvi JK34
  list(c(4), c(1,2)),
  # oijarvi JK36
  list(c(1), c(2,3)),
  # ounasjarvi JK37
  list(c(1), c(2,3), c(4,5)),
  # pielinen JK40
  list(c(1), c(2,3)),
  # rehja JK42
  list(c(1), c(4,5), c(6,7)),
  # visuvesi JK45
  list(c(1), c(2,3), c(4,5)),
  # palovesi JK46
  list(c(1), c(2,3), c(4,5), c(6,7)),
  # simpelejarvi JK47
  list(c(1,2), c(3,4), c(5,6)),
  # ahtarinjarvi JK48
  list(c(2), c(3,4), c(5,6)),
  # kuivajarvi JK49
  list(c(1), c(2,3), c(4,5)),
  # george KMS10
  list(c(3,4), c(1,2), c(5,6)),
  # cazenovia LR2
  list(c(1), c(2,3), c(4,5), c(6,7)), # no trend
  # gull MICH03
  list(c(1,2,3,4), c(5,6)), # check
  # fair MICH06
  list(c(3), c(1,2)),
  # baikal NG1
  list(c(1), c(2,3)),
  # mueggelesee RA1
  list(c(1), c(2,3))
)

#> calculating correlation coefficient for each region 
# ds_final %>%
#   group_by(country) %>%
#   summarise(
#     corr = cor(duration, avg_bt)
#   )

#> add res groups to dataframe
newcol = rep(NA,nrow(ssa_lakes))
ssa_lakes = cbind(ssa_lakes, data.frame(newcol))
ssa_lakes = ssa_lakes %>%
  rename(resgroups = newcol)

for (i in 1:nrow(ssa_lakes)) {
  ssa_lakes$resgroups[i] = resgroups[i]
}

#> t stats for lakes
ds_tstats = tibble(year = integer(), recon = double(), 
                   reconres = double(), country = character(), 
                   lakename = character(), lakecode = character(),
                   tstats = double())
for (i in 1:nrow(ssa_lakes)) {
  if (ssa_lakes$notrend[i] == FALSE) {
    newdat = tstat(ssa_lakes$lakecode[i])
    ds_tstats = ds_tstats %>%
      bind_rows(newdat)
  }
}

#> nordic plot of t stats
ds_tstats %>%
  filter(
    country %in% c("FINLAND", "SWEDEN")
  ) %>%
  # filter(
  #   ! (lakecode %in% c("GW592", "JK09", "JK21"))
  # ) %>%
  ggplot(aes(x = year, y = tstats)) + 
  geom_line(aes(group = lakecode)) + 
  labs(title = "Finland + Sweden t-stats (+outliers)")

#> non nordic plot of t stats
ds_tstats %>%
  filter(
    country %in% c("UNITED STATES", "RUSSIA", "GERMANY")
  ) %>%
filter(
  ! (lakecode %in% c("GW592", "JK09", "JK21", "MICH03", "DMR1"))
) %>%
  ggplot(aes(x = year, y = tstats)) + 
  geom_line(aes(group = lakecode)) + 
  labs(title = "US + Russia + Germany t-stats (-outliers)")

ds_tstats %>%
  filter(tstats < -10)

ds_tstats %>%
  filter(tstats > 10)

#> calculating correlation coefficient for each lake 
corrs <- ds_final %>%
  group_by(lakecode) %>%
  summarise(
    corr = cor(duration, avg_bt)
  )

ds_final %>%
  group_by(country) %>%
  summarise(
    corr = cor(duration, avg_bt)
  )

lakename <- rep(0, nrow(corrs)) -> latitude -> longitude -> country
for (i in 1:nrow(corrs)) {
  dat <- f_lake_data(corrs$lakecode[i])
  lakename[i] <- dat$lakename[1]
  latitude[i] <- dat$latitude[1]
  longitude[i] <- dat$longitude[1]
  country[i] <- dat$country[1]
}

corrs <- cbind(corrs, lakename, latitude, longitude, country)

# write_csv(corrs, file="corrs.csv")

areas <- data.frame(corrs$lakecode) %>%
  rename(lakecode = corrs.lakecode)

sarea_info <- lake_info %>%
  dplyr::select(c(lakecode, surface_area))

areas <- areas %>%
  left_join(sarea_info) %>%
  mutate(surface_area = ifelse(surface_area > 0, surface_area, NA))

elev <- data.frame(corrs$lakecode) %>%
  rename(lakecode = corrs.lakecode)

elev_info <- lake_info %>%
  dplyr::select(c(lakecode, elevation))

elev <- elev %>%
  left_join(elev_info) %>%
  mutate(elevation = ifelse(elevation > 0, elevation, NA))

# pdf("scattertempsize_new.pdf")
# scattertempsize(areas$lakecode, "Ice Cover Duration vs Brightness Temperature and Surface Area")
# dev.off()

# pdf("scattertempelev.pdf")
# scattertempelev(areas$lakecode, "Ice Cover Duration vs Brightness Temperature and Elevation")
# dev.off()

tb_icd_ds = tbicd(areas$lakecode)
tb_icd_ds = na.omit(tb_icd_ds)
tb_icd_ds = scale(tb_icd_ds)

# kmeans clustering
km.out = kmeans(tb_icd_ds, centers = 2, nstart = 20)
cluster_id = factor(km.out$cluster)
tb_icd_ds_km = data.frame(cbind(tb_icd_ds, cluster_id))
kmeansplot = ggplot(tb_icd_ds_km, aes(vec1, vec2, color = cluster_id)) + 
  geom_point(shape = 1) + 
  theme_classic() + 
  xlab("Brightness Temperature") + 
  ylab("Ice Coverage Duration") + 
  ggtitle("K-Means Clustering: k = 2") + 
  theme(legend.position = "none")

# hierarchical clustering
dist_mat = dist(tb_icd_ds, method = "euclidean")
# average
hclust_avg = hclust(dist_mat, method = "average")
cut_avg = cutree(hclust_avg, k = 2)
tb_icd_ds_havg = data.frame(cbind(tb_icd_ds, cut_avg))
hclustavg_plot = ggplot(tb_icd_ds_havg, aes(vec1, vec2, color = cut_avg)) + 
  geom_point(shape = 1) + 
  theme_classic() + 
  xlab("Brightness Temperature") + 
  ylab("Ice Coverage Duration") +  
  ggtitle("Hierarchical Clustering: Method = Average") + 
  theme(legend.position = "none")
# complete
hclust_com = hclust(dist_mat, method = "complete")
cut_com = cutree(hclust_com, k = 2)
tb_icd_ds_hcom = data.frame(cbind(tb_icd_ds, cut_com))
hclustcom_plot = ggplot(tb_icd_ds_hcom, aes(vec1, vec2, color = cut_com)) + 
  geom_point(shape = 1) + 
  theme_classic() + 
  xlab("Brightness Temperature") + 
  ylab("Ice Coverage Duration") + 
  ggtitle("hierarchical clustering: method = complete") + 
  theme(legend.position = "none")
# single
hclust_sin = hclust(dist_mat, method = "single")
cut_sin = cutree(hclust_sin, k = 2)
tb_icd_ds_hsin = data.frame(cbind(tb_icd_ds, cut_sin))
hclustsin_plot = ggplot(tb_icd_ds_hsin, aes(vec1, vec2, color = cut_sin)) + 
  geom_point(shape = 1) + 
  theme_classic() + 
  xlab("Brightness Temperature") + 
  ylab("Ice Coverage Duration") + 
  ggtitle("hierarchical clustering: method = single") + 
  theme(legend.position = "none")
# centroid
hclust_cen = hclust(dist_mat, method = "centroid")
cut_cen = cutree(hclust_cen, k = 2)
tb_icd_ds_hcen = data.frame(cbind(tb_icd_ds, cut_cen))
hclustcen_plot = ggplot(tb_icd_ds_hcen, aes(vec1, vec2, color = cut_cen)) + 
  geom_point(shape = 1) + 
  theme_classic() + 
  xlab("Brightness Temperature") + 
  ylab("Ice Coverage Duration") + 
  ggtitle("hierarchical clustering: method = centroid") + 
  theme(legend.position = "none")

# pdf("kmeansnew.pdf")
# kmeansplot
# dev.off()
# 
# pdf("hclustavgnew.pdf")
# hclustavg_plot
# dev.off()

tb_icd_lat_ds = tbicdlat(areas$lakecode)

tb_icd_lat_ds = na.omit(tb_icd_lat_ds)
tb_icd_lat_ds = data.frame(tb_icd_lat_ds)
tb_icd_lat_ds = tb_icd_lat_ds %>%
  mutate(latitude = factor(ifelse(vec3 >= 60, "above 60", "below 60")),
         vec1 = as.numeric(vec1),
         vec2 = as.numeric(vec2))
lat_plot = ggplot(tb_icd_lat_ds, aes(x = vec1, y = vec2, color = latitude)) + 
  scale_color_manual(values = c("#000000", "cyan3")) +
  geom_point(alpha = 1, size = 2, shape=1) + 
  theme_classic() + 
  xlab("Brightness Temperature") + ylab("Ice Coverage Duration") + 
  ggtitle("Ice Coverage Duration vs Brightness Temperature: Latitude")

# pdf("icd_tb_latnew.pdf")
# lat_plot
# dev.off()

lake.codes <- corrs$lakecode
lake.codes <- lake.codes[-which(lake.codes=="LR2")]

lake.tstats <- tstat(lake.codes[1])
for (i in 2:length(lake.codes)) {
  lake.tstats <- rbind(lake.tstats, tstat(lake.codes[i]))
}


#> plot of t-statistics
plot(NULL, xlim=c(1977, 2021), ylim=c(-5, 5), 
     main="Lake Ice Cover Duration t-statistics",
     xlab="Year",
     ylab="t-stat")
for (i in 1:length(lake.codes)) {
  this.code <- lake.codes[i]
  this.lake <- lake.tstats %>%
    filter(lakecode == this.code)
  this.country <- this.lake$country[1]
  if (this.country %in% c("SWEDEN", "FINLAND")) {
    line.color <- "#56B4E9"
  } else if (this.country == "UNITED STATES") {
    line.color <- "#E69F00"
  } else if (this.country == "RUSSIA") {
    line.color <- "#CC79A7"
  } else if (this.country == "GERMANY") {
    line.color <- "#000000"
  }
  lines(x=this.lake$year, y=this.lake$tstats, col=line.color)
}
abline(h=0, lty=2, col="darkgrey")
abline(h=1.96, lty=2)
abline(h=-1.96, lty=2)
legend("bottomleft",
       legend=c("Nordic", "United States", "Russia", "Germany"),
       col = c("#56B4E9", "#E69F00", "#CC79A7", "#000000"),
       pch=19,
       bty="n",
       inset=c(0.01,0.01))
# caption: explain jukkasjarvi has large negative outlier

#> creating >20 years dataset
ds_20 = f_specdattemp2(more_lakes_20[1])
for (i in 2:length(more_lakes_20)) {
  ds_20 = rbind(ds_20, f_specdattemp2(more_lakes_20[i]))
}

corrs20 <- ds_20 %>%
  group_by(lakecode) %>%
  summarise(
    corr = cor(duration, avg_bt)
  )

lakename20 <- rep(0, nrow(corrs20)) -> latitude20 -> longitude20 -> country20
for (i in 1:nrow(corrs20)) {
  dat <- f_lake_data(corrs20$lakecode[i])
  lakename20[i] <- dat$lakename[1]
  latitude20[i] <- dat$latitude[1]
  longitude20[i] <- dat$longitude[1]
  country20[i] <- dat$country[1]
}

corrs20 <- cbind(corrs20, lakename20, latitude20, longitude20, country20)

# write_csv(corrs20, file="corrs20.csv")

ds_20 %>%
  filter(country %in% c("UNITED STATES", "SWEDEN", "FINLAND")) %>%
  filter(!(lakecode %in% c("JGL02", "JLG01", "MINN42", "MINN34", "MINN2"))) %>%
  ggplot(aes(x = avg_bt, y = duration)) + 
  geom_point(shape = "O") + 
  theme_classic() + 
  ggtitle("Average Brightness Temperature vs Ice Cover Duration") + 
  xlab("Average Brightness Temperature (K)") +
  ylab("Ice Cover Duration (Days)")

#> trend plot of just NA and Nord
trend.df <- data.frame(year=integer(),
                       recon=double(),
                       country=character(),
                       lakename=character(),
                       lakecode=character())
for (i in lake.codes) {
  if (f_lake_data(i)[1,"country"] %in% c("UNITED STATES", "SWEDEN", "FINLAND")) {
    trend.df <- rbind(trend.df, addNA(yeardur(i)))
  }
}
trend.df %>%
  mutate(region = ifelse(country=="UNITED STATES", country, "NORDIC")) %>%
  ggplot(aes(x=year,y=recon,group=lakecode,color=region)) + 
  geom_line() + 
  scale_color_manual(values=c("#56B4E9", "#E69F00")) + 
  theme_classic() + 
  labs(title = "Lake Ice Cover Duration Trend",
       x = "Year",
       y = "Reconstructed Duration") + 
  theme(legend.position = "inside", legend.position.inside = c(0.11,0.1))


#> trend plot of all lakes
trend.df <- data.frame(year=integer(),
                       recon=double(),
                       country=character(),
                       lakename=character(),
                       lakecode=character())
for (i in lake.codes) {
  trend.df <- rbind(trend.df, addNA(yeardur(i)))
}
trend.df %>%
  mutate(region = ifelse(country %in% c("SWEDEN", "FINLAND"), "NORDIC", country)) %>%
  ggplot(aes(x=year,y=recon,group=lakecode,color=region)) + 
  geom_line() + 
  scale_color_manual(values=c("#000000", "#56B4E9", "#CC79A7", "#E69F00")) + 
  theme_classic() + 
  labs(title = "Lake Ice Cover Duration Trend",
       x = "Year",
       y = "Reconstructed Duration") + 
  theme(legend.position = "inside", legend.position.inside = c(0.1,0.12))

#> plot of first periods
firstpers <- f_period(lake.codes)
firstpers_plot <- firstpers %>%
  ggplot(aes(x=period1)) + 
  geom_histogram(breaks=seq(2,21,1), color="black") + 
  labs(x="Length of First Period",
       y="Count") + 
  theme_classic() + 
  coord_cartesian(ylim=c(0,13))

#> plot of minimum periods
minpers <- f_minperiod(lake.codes)
minpers <- minpers %>% 
  left_join(corrs, by="lakecode") %>%
  dplyr::select(c("lakecode", "minperiod", "country")) %>%
  mutate(region = ifelse(country %in% c("SWEDEN", "FINLAND"), "NORDIC", country)) %>%
  dplyr::select(-"country")
minper_plot <- minpers %>%
  ggplot(aes(x=minperiod, fill=region)) + 
  geom_histogram(breaks=seq(2,21,1), color="black") + 
  scale_fill_manual(values=c("#000000", "#56B4E9", "#CC79A7", "#E69F00")) +
  labs(title="Minimum Ice Cover Duration Period Length",
       x="Length of Minimum Period",
       y="Count") + 
  theme_classic() +
  coord_cartesian(ylim=c(0,25)) + 
  theme(legend.position = "inside", legend.position.inside = c(0.85,0.85))

# pdf("minperiod_hist.pdf")
# minper_plot
# dev.off()

#> plot of all periods
allpers <- f_allperiod(lake.codes)
allpers <- allpers %>%
  left_join(corrs, by="lakecode") %>%
  dplyr::select(c("lakecode", "period", "country")) %>%
  mutate(region = ifelse(country %in% c("SWEDEN", "FINLAND"), "NORDIC", country)) %>%
  dplyr::select(-"country")
allpers_plot <- allpers %>%
  ggplot(aes(x=period, fill=region)) + 
  geom_histogram(breaks=seq(2,60,1), color="black") + 
  scale_fill_manual(values=c("#000000", "#56B4E9", "#CC79A7", "#E69F00")) +
  labs(title="Ice Cover Duration Period Length: 0 to 60 Years",
       x="Length of Period",
       y="Count") + 
  theme_classic() +
  # coord_cartesian(ylim=c(0,25)) + 
  theme(legend.position = "inside", legend.position.inside = c(0.85,0.85)) +   
  scale_x_continuous(breaks=seq(0,60,10)) + 
  scale_y_continuous(breaks=seq(0,30,5))

# pdf("some_period_hist.pdf")
# allpers_plot
# dev.off()

#> summary statistics about the periods
# allpers %>%
#   summarise(
#     avg = mean(period),
#     sdev = sd(period),
#     med = median(period),
#     iqr = IQR(period)
#   )

#> export pdf plots of ssa output
# for (i in 1:nrow(ssa_lakes)) {
#   f_ssa_plots_ind(ssa_lakes$lakecode[i])
# }

#> export pdf plots of reconstructed series
# for (i in 1:nrow(ssa_lakes)) {
#   f_reconpdf(ssa_lakes$lakecode[i])
# }

#> export pdf plots of residuals ssa output
# for (i in 1:nrow(ssa_lakes)) {
#   resssa2pdf(ssa_lakes$lakecode[i])
# }

#> export pdf plots of residuals reconstructed series
# for (i in 1:nrow(ssa_lakes)) {
#   resssathetapdf(ssa_lakes$lakecode[i])
# }

#> table for Supplementary Materials
ssa_lakes_copy <- ssa_lakes
# ssa_lakes_copy <- ssa_lakes_copy %>% 
#   mutate(lakename=0)
for (i in 1:nrow(ssa_lakes_copy)) {
  ssa_lakes_copy$groups[i] <- gsub("c", paste(as.character(ssa_lakes$groups[i][[1]]), collapse=', '), replacement="")
  ssa_lakes_copy$resgroups[i] <- gsub("c", paste(as.character(ssa_lakes$resgroups[i][[1]]), collapse=', '), replacement="")
  ssa_lakes_copy$lakename[i] <- f_lake_data(ssa_lakes_copy$lakecode[i])$lakename[1]
}
ssa_lakes_copy <- cbind(ssa_lakes_copy, latitude=rep(0, nrow(ssa_lakes_copy)))
ssa_lakes_copy <- cbind(ssa_lakes_copy, longitude=rep(0, nrow(ssa_lakes_copy)))
for (i in 1:nrow(ssa_lakes_copy)) {
  thiscoord <- justcoords(ssa_lakes_copy$lakecode[i])
  ssa_lakes_copy$latitude[i] <- thiscoord[1]
  ssa_lakes_copy$longitude[i] <- thiscoord[2]
}
ssa_lakes_copy <- ssa_lakes_copy %>%
  left_join(areas, by='lakecode')
ssa_lakes_copy <- ssa_lakes_copy %>%
  dplyr::select(lakecode, lakename, latitude, longitude, surface_area, notrend, groups, resgroups) %>%
  mutate(notrend = ifelse(notrend==FALSE, '0', '1')) %>%
  rename(
    "Code" = lakecode,
    "Name" = lakename,
    "Latitude" = latitude,
    "Longitude" = longitude,
    "Surface Area" = surface_area,
    "No Trend" = notrend,
    "Groups" = groups,
    "Residual Groups" = resgroups
  )







