.libPaths(c(.libPaths(), "C:/Users/mem446/R/win-library/3.3", "C:/Program Files/R/R-3.3.1/library"))

library('raster')
library('rgeos')
library('rgdal')
library('plyr')
library('dplyr')
library('parallel')
library('foreach')
library('doParallel')

source("C:/Users/mem446/Desktop/CRTV_HWA/HWA_spread_wrap_CRTV5.R")
mydir <- "C:/Users/mem446/Desktop/CRTV_HWA"
timestep <- readLines(paste0(mydir, "/lockfile"),n = 1)
timestep <- as.numeric(timestep)
if(timestep>0){
  HWA_spread_par(timestep, mydir)

  input_inf <- paste0(mydir, "/output/HWA/inf-level-", timestep, ".img")
  input_lu <- paste0(mydir, "/land-use-input/land-use-", timestep, ".img")
  output <- paste0(mydir, "/land-use-input/land-use-", timestep+5,".img")
  input_infR <- raster(input_inf)#; names(input_infR) <- "inf_level"
  input_luR <- raster(input_lu)
  # AFT <- readOGR(dsn = mydir, layer = "parcels-lev2") #raster("testAFT.img")
  # 
  # inf_lev <- extract(input_infR, AFT, fun=max, sp=TRUE)
  # inf_lev$lu_cut <- as.factor(as.numeric(paste0(1,inf_lev$cut, inf_lev$inf_level)))
  # 
  # level <- sort(unique(inf_lev$lu_cut), index.return = FALSE)
  # level_df <- as.data.frame(level)
  # i=0
  # for(n in level){
  #   i = i+1
  #   inf_lev$ranks[inf_lev$lu_cut == n] <- i
  #   level_df$rank[level_df[,1]==n] <- i
  # }
  # level_df <- level_df[,c(2,1)]
  # 
  # lu_ranks <- rasterize(as(inf_lev, "SpatialPolygons"), input_luR, field=inf_lev@data[,"ranks"])
  # lu_cutr <- reclassify(lu_ranks, level_df)
  # lu <- overlay(lu_cutr, input_luR, fun=function(x,y) {ifelse(is.na(x), y, x)})
  
  lu_cutr <- input_infR+10
  s <- stack(lu_cutr, input_luR)
  lu <- calc(s, fun=function(x) {ifelse(is.na(x[1]), x[2], x[1])})
  lu[is.na(lu)] <- 0
  writeRaster(lu,output,overwrite = T,format="HFA", datatype="INT2S")
}
