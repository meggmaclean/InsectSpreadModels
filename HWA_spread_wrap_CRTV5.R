
HWA_spread_par <- function(timestep, mydir){
  cl <- makeCluster(3)
  registerDoParallel(cl)
  
  start_year <- 2010+timestep
  years = c(4, 3, 2, 1, 0)
  
###################################
## General Inputs
###################################
tsuga_in <- raster(paste0(mydir, "/output/agbiomass/tsugcana/AGBiomass_",timestep-5,".img"))# tsuga agbiomass raster at timestep

ecoregion_in <- raster(paste0(mydir,"/CTRV_eco_100_feb_2018.img")) # ecoregion map to grab temps
ecoregion <- getValues(ecoregion_in)
eco_lookup <- read.table(paste0(mydir, "/CTRV_march2017_current_rcp85_pnet_EcoregionParameters2.txt"), skip = 1, header = TRUE)

PRJ<-projection(ecoregion_in)
EXT<-extent(ecoregion_in)

remove(ecoregion_in)

###################################
## Timestep Inputs
###################################

hwa_in <- raster(paste0(mydir, "/output/HWA/hwa-pop-", timestep-5,".img"))# hwa population raster at timestep
hwa_pop <- getValues(hwa_in)

projection(tsuga_in)<-PRJ
extent(tsuga_in)<-EXT
projection(hwa_in)<-PRJ
extent(hwa_in)<-EXT

compareRaster(tsuga_in, hwa_in, stopiffalse = TRUE, showwarning = TRUE)

###################################
## Outputs
###################################
hwa_out <- paste0(mydir, "/output/HWA/hwa-pop-", timestep, ".img")
inf_out <- paste0(mydir, "/output/HWA/hwa-inf-lev-", timestep, ".img")

###################################
## Functions
###################################

# Infestated cells advance on a yearly basis
N_adv <- function(N.0, K){
  N_t <- 8*N.0 # from McAvoy et al. 2017 <- number not constrained well
  ifelse(N_t>=K, K, N_t)
}

# Winter temp mortality
Mw <- function(WT){
  1/(1+exp(-(-5.7290+(-0.3311*WT)))) # from McAvoy et al. 2017
}

# Neighbors for 5 years is 5x a 1yr matrix
who3 <- matrix(1, ncol=1499, nrow = 1499) # up to 75 km away from source
who3[750,750] <- 0
who4 <- matrix(1, ncol=29, nrow = 29) # up to 1.5 km away from source
who4[15,15] <- 0

neigh <- function(num, who){
  raster::adjacent(tsuga_in, num, directions = who, pairs = FALSE, target = Uninfected_cells, include = FALSE)
}

inf_l <- function(r1, r2){
  ifelse(r2>0, ifelse(r1>100, ifelse(r1>0.4*r2, ifelse(r1>0.9*r2, 3, 2), 1), 0), 0)
}

###################################
## Current Locations and initializing rasters
###################################
hwa_curr <- hwa_in # hwa_in needs to have the same raster size/cells etc as tsuga_in
hwa_cells <- which(hwa_pop>100)
tsuga_agb <- getValues(tsuga_in)
tsuga_cells <- which(tsuga_agb>0)
Uninfected_cells <- setdiff(tsuga_cells, hwa_cells)

K_in <- tsuga_in*(139680*30/13.1563) # compute carry capacity
K_val <- getValues(K_in)

compareRaster(K_in, hwa_curr, stopiffalse = TRUE, showwarning = TRUE)
remove(hwa_in, tsuga_agb, tsuga_cells)

###################################
## Disperse!
###################################
#@ for loop converted to run in parallel with foreach loop
hwa_par <- foreach(i = hwa_cells, .combine = 'rbind') %dopar% {
  K <- K_val[i]
  N.1 <- hwa_pop[i]

## Dynamics
  eco_i <- paste0("eco", ecoregion[i])
  eco_file <- as.character(eco_lookup$climateFileName[match(eco_i, eco_lookup$EcoregionParameters)])
  eco_vals <- read.table(paste0(mydir, "/", eco_file), header = TRUE)
  for(y in years){
    start_yeary = start_year - y
    
  # Minimum Winter Temp
    WT_dec <- eco_vals$Tmin[eco_vals$Year == start_yeary-1 & eco_vals$Month == 12]
    WT_jan <- eco_vals$Tmin[eco_vals$Year == start_yeary & eco_vals$Month == 1]
    WT_feb <- eco_vals$Tmin[eco_vals$Year == start_yeary & eco_vals$Month == 2]
    WT_march <- eco_vals$Tmin[eco_vals$Year == start_yeary & eco_vals$Month == 3]

  # What is the yearly spring HWA population
    WT <- min(c(WT_march, WT_feb, WT_jan, WT_dec))
    Wmort <- Mw(WT)
    remove(WT_dec, WT_jan, WT_feb, WT_march, WT)

    N.max <- N_adv(N.1, K)
    N.1 <- N.max - Wmort*N.max; remove(Wmort, N.max) # summer mortality?
    N.1 <- ifelse(N.1>0, N.1, 0)
  }
  remove(eco_vals, eco_file)
  
  #@ created two column matrix to store cell number and cell value of any newly infested cells - this will be returned by the foreach loop
  hwa_return <- matrix(c(i, N.1), nrow = 1, ncol = 2)
  
##### @@@@@
  if(N.1 > 1000) {
# Where my neighbors at?!
    neighbors4 <- neigh(i, who4) 
    
    if(N.1/K > 0.04){
      neighbors3 <- neigh(i, who3)
      neighbors3s <- setdiff(neighbors3, neighbors4); remove(neighbors3)
      sam_n3s <- round((length(neighbors3s)*0.09)*(N.1/K)) # up to 9% of total uninfested cells randomly chosen for infestation
      neighbors3s <- sample(neighbors3s, sam_n3s, replace = FALSE); remove(sam_n3s)
      N.n1 <- round((N.1*0.000000001)+101) # from Fitzpatrick et al. 2012
      N.n1 <- ifelse(N.n1>0, N.n1, 0)

      for(l in neighbors3s){
        hwa_return <- rbind(hwa_return, c(l, N.n1))
      }
      remove(neighbors3s, N.n1)
    }

    for(m in neighbors4){
      N.n0 <- rbinom(1,1,0.9) # 90% chance of being infested
      N.n1 <- round(N.n0*((N.1/4800000)+101)) # from Fitzpatrick et al. 2012
      N.n1 <- ifelse(N.n1>0, N.n1, 0)
      hwa_return <- rbind(hwa_return, c(m, N.n1))
      remove(N.n0, N.n1)
    }
    remove(neighbors4)
  }
  
  #@ return hwa_return as the result of the foreach loop
  return(hwa_return)
}

##### @@@
#@ summarise values in hwa_par
hwa_par <- ddply(.data = data.frame(hwa_par), .(X1), summarize, sum(X2))
stopCluster(cl)

#@ update hwa_curr raster with new infestation values from values from hwa_par
for (row in 1:nrow(hwa_par)) {
  cell_num <- hwa_par[row,1]
  hwa_new <- hwa_curr[cell_num]
  hwa_new <- ifelse(is.na(hwa_new),0, hwa_new) + hwa_par[row,2]
  K_par <- K_val[cell_num]
  hwa_curr[cell_num] <- ifelse(hwa_new > K_par, K_par, hwa_new)
  remove(cell_num, hwa_new, K_par)
}

###################################
## Plots and exports
###################################

# Export current HWA pop img
writeRaster(hwa_curr, hwa_out, format = "HFA", overwrite = TRUE)

# Calculate and export infestation level img
hwa_inf <- hwa_curr; hwa_inf <- calc(hwa_curr, fun = function(x){x[x>0] <- 0; return(x)})
hwa_inf <- overlay(hwa_curr, K_in, fun = inf_l)

writeRaster(hwa_inf, inf_out, format = "HFA", overwrite = TRUE, datatype="INT2S")

}
