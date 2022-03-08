#############################################################################################################
# upload packages
#############################################################################################################

library(dplyr)
library(rlist)
library(tibble)

#############################################################################################################
# usefull functions
#############################################################################################################

# function to convert list into df
convert_list_df <- function(df_list){
  df <- list()
  for (i in 1:length(df_list)) {
    df <- rbind(df, df_list[[i]])
  }
  df <- data.frame(df)
  return(df)
}

#############################################################################################################
#############################################################################################################
# comparison between the modelled and measured absorption coefficient (Babs)
#############################################################################################################
#############################################################################################################

#############################################################################################################
# import the calculated Mie coefficients
#############################################################################################################

# set wd
setwd("Z:/Clarissa/Data_Optical_calculation/R_Code/Code_final/ICEdust_samples/Maeli2/datasets/Python_output/original/Mie_coef")

# create a list using the selected pattern
filenames <- list.files(pattern="*.csv", full.names=TRUE)

# import the selected files which correspond to the Mie coefficients calculated at different wavelengths (WL)
Datafile <- list()

for (i in 1:length(filenames)){
  Datafile[[i]] <- read.csv(filenames[[i]], sep = ",")
  Datafile[[i]]$Time  <- as.POSIXct(Datafile[[i]]$Time ,  format =  "%Y-%m-%d %H:%M:%S") # transform as POSIXct
  colnames(Datafile[[i]])[3] <- "date"
}

# subset the data from 30 min to 2.5h after the dust injection peak
for (i in 1:length(Datafile)){
  Datafile[[i]] <- Datafile[[i]][which(Datafile[[i]]$date >= "2019-01-21 11:31:00" & Datafile[[i]]$date <= "2019-01-21 13:31:00"),]
}

#############################################################################################################
# get Babs
#############################################################################################################

# set new wd
setwd("Z:/Clarissa/Data_Optical_calculation/R_Code/Code_final/ICEdust_samples/Maeli2/Nep_Aet_Mie/original/")

Babs_calculated <- Datafile

# select Xnk, date, and Babs columns in the dfs
# X = dynamic shape factor
# n = real part of the complex refractive index 
# k = imaginary part of the complex refractive index
# Babs is calculated for a number of X-n-k scenarios
# X, n, and k parameters were used to correct size distribution measurements from GRIMM and SMPS
# Babs is calculated based on the corrected size distribution measurements

for (i in 1:length(Datafile)){
  Babs_calculated[[i]] <- Babs_calculated[[i]][, c(2,3,6)]
}

# change the name of column 3 (Babs) for each df in the list
# each df in the list correspond to Bsca calculated at a certain WL
# each df contains Babs calculated for different X-n-k scenarios

WL <- c(370, 450, 470, 520, 550, 590, 660, 700, 880, 950)

for (i in 1:length(Babs_calculated)) {
  colnames(Babs_calculated[[i]])[3] <- sprintf("Babs_mod_%d", WL[i])
}

# remove data at 450, 550, and 700 nm as unnecessary for the comparison
rm_WL <- c("Babs_mod_450", "Babs_mod_550", "Babs_mod_700")

rm_data <- list()
for (i in 1:length(Babs_calculated)) {
  if ((colnames(Babs_calculated[[i]])[3] %in% rm_WL) == F){
    rm_data[[i]] <- Babs_calculated[[i]]}
  else {
    rm_data[[i]] <- NULL
  }
}

# remove Null terms from the list
Babs_calculated <- rm_data[lengths(rm_data) != 0]

#############################################################################################################
# split the data based on the Xnk value
#############################################################################################################

# create a function to split the data based on the Xnk value
split_data_byXnk <- function(df) {
  
  # find the Xnk range
  Xnk_range <- unique(df$Xnk)
  
  # create a list with the calculated Babs for each scenario
  WL_scenario <- list()
  for (i in 1:length(Xnk_range)) {
    WL_scenario[[i]] <- df[which(df$Xnk %in% Xnk_range[i] == TRUE),]
    WL_scenario[[i]]$Xnk <- as.character(WL_scenario[[i]]$Xnk) # convert from 1 level factor to charachter
  }
  
  return(WL_scenario) 
  
}

# apply the function to all dfs in the list
Babs_calculated_scenarios <- list()

for (i in 1:length(Babs_calculated)) {
  Babs_calculated_scenarios[[i]] <- split_data_byXnk(Babs_calculated[[i]])
}

#############################################################################################################
# import the measured Babs
#############################################################################################################

# import data
Babs_data <- read.csv("Z:/Clarissa/Data_Optical_calculation/R_Code/Code_final/ICEdust_samples/Maeli2/Nep_Aet_Mie/original/data/Maeli2_Babs_12min_original.csv")

Babs_data$date  <- as.POSIXct(Babs_data$date ,  format =  "%Y-%m-%d %H:%M:%S") # transform as POSIXct

# subset the data from 30 min to 2.5h after the dust injection peak
Babs_data <- Babs_data[which(Babs_data$date >= "2019-01-21 11:31:00" & Babs_data$date <= "2019-01-21 13:31:00"),]

# get mean Babs measured values only
Babs <- Babs_data[,1:8]

# get error Babs measured values only
Babs_error <- Babs_data[,-c(2:8, 16:22)]

# give column names
Aet_WL <- c(370, 470, 520, 590, 660, 880, 950)

for (i in 2:length(Babs)) {
  colnames(Babs)[i] <- sprintf("Babs_meas_%d", Aet_WL[i-1])
}

#############################################################################################################
# ensure calculated and measured data have the same number of observations
#############################################################################################################

Babs$date  <- as.POSIXct(Babs$date ,  format =  "%Y-%m-%d %H:%M:%S") # transform as POSIXct

for (i in 1:length(Babs_calculated_scenarios)) {
  for (k in 1:length(Babs_calculated_scenarios[[i]])) {
    
    Babs_calculated_scenarios[[i]][[k]]$date <- as.POSIXct(Babs_calculated_scenarios[[i]][[k]]$date ,  format =  "%Y-%m-%d %H:%M:%S") # transform as POSIXct
    
  }
}

# join Babs calculated and measured dfs
Babs_comparison <- Babs_calculated_scenarios

for (i in 1:length(Babs_calculated_scenarios)) {
  for (k in 1:length(Babs_calculated_scenarios[[i]])) {
    
    Babs_comparison[[i]][[k]] <- full_join(Babs_calculated_scenarios[[i]][[k]], Babs[,c(1,i+1)])
    
  }
}

# remove NAs
for (i in 1:length(Babs_comparison)) {
  for (k in 1:length(Babs_comparison[[i]])) {
    Babs_comparison[[i]][[k]] <- na.omit(Babs_comparison[[i]][[k]])
  }
}

#############################################################################################################
# add the X, n,  k columns to the dfs 
#############################################################################################################

add_Xnk <- function(df) {
  
  # convert Xnk column as charachter
  df$Xnk <- as.character(df$Xnk)
  # add X, n, k columns to df
  df$X <- NA
  df$n <- NA
  df$k <- NA
  
  for (i in 1:nrow(df)) {
    df$X[i] <- unlist(strsplit(df$Xnk[i]," "))[1]
    df$n[i] <- unlist(strsplit(df$Xnk[i]," "))[2]
    df$k[i] <- unlist(strsplit(df$Xnk[i]," "))[3]
  }
  
  return(df)
  
}

for (i in 1:length(Babs_comparison)) {
  for (k in 1:length(Babs_comparison[[i]])) {
    Babs_comparison[[i]][[k]] <- add_Xnk(Babs_comparison[[i]][[k]])
  }
}

# convert the X-n-k values as numeric
for (i in 1:length(Babs_comparison)) {
  for (k in 1:length(Babs_comparison[[i]])) {
    Babs_comparison[[i]][[k]]$X <- as.numeric(as.character(Babs_comparison[[i]][[k]]$X))
    Babs_comparison[[i]][[k]]$n <- as.numeric(as.character(Babs_comparison[[i]][[k]]$n))
    Babs_comparison[[i]][[k]]$k <- as.numeric(as.character(Babs_comparison[[i]][[k]]$k))
  }
}

#############################################################################################################
# check the correlation and calculate the statistics
#############################################################################################################

# prepare the data
for (i in 1:length(Babs_comparison)) {
  for (k in 1:length(Babs_comparison[[i]])) {
    colnames(Babs_comparison[[i]][[k]])[3] <- "model"
    colnames(Babs_comparison[[i]][[k]])[4] <- "obs"
  }
}

# remove the combinations with X = 1
rm_X1 <- list()
for (i in 1:length(Babs_comparison)) {
  rm_X1[[i]] <- list()
}

for (i in 1:length(Babs_comparison)){
  for (k in 1:length(Babs_comparison[[i]])){
    if (Babs_comparison[[i]][[k]]$X[1] != 1){
      rm_X1[[i]][[k]] <- Babs_comparison[[i]][[k]]
    }
  }
}

# remove Null terms from the list
for (i in 1:length(rm_X1)){
  rm_X1[[i]] <- rm_X1[[i]][lengths(rm_X1[[i]]) != 0]
}


#############################################################################################################
# check the correlation

library(stats)
library(Metrics)


test_corr <- rm_X1
for (i in 1:length(rm_X1)) {
  for (k in 1:length(rm_X1[[i]])) {
    test_corr[[i]][[k]] <- summary(lm(model ~ obs, rm_X1[[i]][[k]]))
  }
}

# get the r.squared value
test_R2 <-  rm_X1
for (i in 1:length(rm_X1)) {
  for (k in 1:length(rm_X1[[i]])) {
    test_R2[[i]][[k]] <- round(test_corr[[i]][[k]][["r.squared"]],1)
  }
}

# calculate RMSE
test_RMSE <-  rm_X1
for (i in 1:length(rm_X1)) {
  for (k in 1:length(rm_X1[[i]])) {
    test_RMSE[[i]][[k]] <- rmse(rm_X1[[i]][[k]]$obs, rm_X1[[i]][[k]]$model)/sd(rm_X1[[i]][[k]]$obs)*100
  }
}

#############################################################################################################
# for each WL, create a single df containing a summary of the statistic for each X-n-k scenario
#############################################################################################################

# Test for 1 WL
df_WL = Aet_WL[1]
df_list = rm_X1[[1]]
df_list_RMSE = test_RMSE[[1]]
df_list_R2 = test_R2[[1]]

stat_comparison <- data.frame(matrix(ncol = 7, nrow = length(df_list))) # the length of tested scenario
colnames(stat_comparison) <- c("WL", "Xnk", "X", "n", "k", "RMSE", "R2")

for (i in 1:length(df_list)) {
  stat_comparison$WL[i] <- df_WL
  stat_comparison$Xnk[i] <- df_list[[i]]$Xnk[1]
  stat_comparison$X[i] <- df_list[[i]]$X[1]
  stat_comparison$n[i] <- df_list[[i]]$n[1]
  stat_comparison$k[i] <- df_list[[i]]$k[1]
  stat_comparison$RMSE[i] <- df_list_RMSE[[i]]
  stat_comparison$R2[i] <- df_list_R2[[i]]
}

#############################################################################################################
# function

build_stat_df <- function(df_WL, df_list, df_list_RMSE, df_list_R2) {
  stat_comparison <- data.frame(matrix(ncol = 7, nrow = length(df_list))) # the length of tested scenario
  colnames(stat_comparison) <- c("WL", "Xnk", "X", "n", "k", "RMSE", "R2")
  
  for (i in 1:length(df_list)) {
    stat_comparison$WL[i] <- df_WL
    stat_comparison$Xnk[i] <- df_list[[i]]$Xnk[1]
    stat_comparison$X[i] <- df_list[[i]]$X[1]
    stat_comparison$n[i] <- df_list[[i]]$n[1]
    stat_comparison$k[i] <- df_list[[i]]$k[1]
    stat_comparison$RMSE[i] <- df_list_RMSE[[i]]
    stat_comparison$R2[i] <- df_list_R2[[i]]
  }
  
  return(stat_comparison)
}

#############################################################################################################
# apply the function to all df

test <- build_stat_df(Aet_WL[1], rm_X1[[1]], test_RMSE[[1]], test_R2[[1]])

identical(test, stat_comparison)

stat_comparison <- list()
for (i in 1:length(rm_X1)) {
  stat_comparison[[i]] <- build_stat_df(Aet_WL[i], rm_X1[[i]], test_RMSE[[i]], test_R2[[i]])
}

Babs_OG_results <- stat_comparison

#############################################################################################################
# Repeat the calculations after subtracting the error to the measured Babs
#############################################################################################################

# remove terms from R environment except for:
rm(list=setdiff(ls(), c("Bsca_OG_results", "Bsca_pSD1_results","Bsca_mSD1_results",
                        "Babs_OG_results", "Babs_pSD1_results","Babs_mSD1_results",
                        "convert_list_df", "add_Xnk", "build_stat_df", "split_data_byXnk", 
                        "Babs_calculated_scenarios", "Babs_data")))


# get mean Babs measured values only
Babs <- Babs_data[,1:8]

# get error Babs measured values only
Babs_error <- Babs_data[,-c(2:8, 16:22)]

# subtract the error to Babs
for (i in 2:length(Babs)) {
  Babs[,i] <- Babs[,i] - Babs_error[,i]
}

# give column names
Aet_WL <- c(370, 470, 520, 590, 660, 880, 950)

for (i in 2:length(Babs)) {
  colnames(Babs)[i] <- sprintf("Babs_meas_%d", Aet_WL[i-1])
}

#############################################################################################################
# ensure calculated and measured data have the same number of observations
#############################################################################################################

Babs$date  <- as.POSIXct(Babs$date ,  format =  "%Y-%m-%d %H:%M:%S") # transform as POSIXct

# join Babs calculated and measured dfs
Babs_comparison <- Babs_calculated_scenarios

for (i in 1:length(Babs_calculated_scenarios)) {
  for (k in 1:length(Babs_calculated_scenarios[[i]])) {
    Babs_comparison[[i]][[k]] <- full_join(Babs_calculated_scenarios[[i]][[k]], Babs[,c(1,i+1)])
  }
}

# remove Nas
for (i in 1:length(Babs_comparison)) {
  for (k in 1:length(Babs_comparison[[i]])) {
    Babs_comparison[[i]][[k]] <- na.omit(Babs_comparison[[i]][[k]])
  }
}

#############################################################################################################
# add the X, n,  k columns to the dfs 
#############################################################################################################

for (i in 1:length(Babs_comparison)) {
  for (k in 1:length(Babs_comparison[[i]])) {
    Babs_comparison[[i]][[k]] <- add_Xnk(Babs_comparison[[i]][[k]])
  }
}

# convert the X-n-k values as numeric
for (i in 1:length(Babs_comparison)) {
  for (k in 1:length(Babs_comparison[[i]])) {
    Babs_comparison[[i]][[k]]$X <- as.numeric(as.character(Babs_comparison[[i]][[k]]$X))
    Babs_comparison[[i]][[k]]$n <- as.numeric(as.character(Babs_comparison[[i]][[k]]$n))
    Babs_comparison[[i]][[k]]$k <- as.numeric(as.character(Babs_comparison[[i]][[k]]$k))
  }
}

#############################################################################################################
# check the correlation and calculate the statistics
#############################################################################################################

# prepare the data
for (i in 1:length(Babs_comparison)) {
  for (k in 1:length(Babs_comparison[[i]])) {
    colnames(Babs_comparison[[i]][[k]])[3] <- "model"
    colnames(Babs_comparison[[i]][[k]])[4] <- "obs"
  }
}

# remove the combinations with X = 1
rm_X1 <- list()
for (i in 1:length(Babs_comparison)) {
  rm_X1[[i]] <- list()
}

for (i in 1:length(Babs_comparison)){
  for (k in 1:length(Babs_comparison[[i]])){
    if (Babs_comparison[[i]][[k]]$X[1] != 1){
      rm_X1[[i]][[k]] <- Babs_comparison[[i]][[k]]
    }
  }
}

# remove Null terms from the list
for (i in 1:length(rm_X1)){
  rm_X1[[i]] <- rm_X1[[i]][lengths(rm_X1[[i]]) != 0]
}


#############################################################################################################
# check the correlation

test_corr <- rm_X1
for (i in 1:length(rm_X1)) {
  for (k in 1:length(rm_X1[[i]])) {
    test_corr[[i]][[k]] <- summary(lm(model ~ obs, rm_X1[[i]][[k]]))
  }
}

# get the r.squared value
test_R2 <-  rm_X1
for (i in 1:length(rm_X1)) {
  for (k in 1:length(rm_X1[[i]])) {
    test_R2[[i]][[k]] <- round(test_corr[[i]][[k]][["r.squared"]],1)
  }
}

# calculate RMSE
test_RMSE <-  rm_X1
for (i in 1:length(rm_X1)) {
  for (k in 1:length(rm_X1[[i]])) {
    test_RMSE[[i]][[k]] <- rmse(rm_X1[[i]][[k]]$obs, rm_X1[[i]][[k]]$model)/sd(rm_X1[[i]][[k]]$obs)*100
  }
}

#############################################################################################################
# for each WL, create a single df containing a summary of the statistic for each X-n-k scenario
#############################################################################################################

stat_comparison <- list()
for (i in 1:length(rm_X1)) {
  stat_comparison[[i]] <- build_stat_df(Aet_WL[i], rm_X1[[i]], test_RMSE[[i]], test_R2[[i]])
}

Babs_mSD1_results <- stat_comparison

#############################################################################################################
# Repeat the calculations after adding the error to the measured Babs
#############################################################################################################

# remove terms from R environment except for:
rm(list=setdiff(ls(), c("Bsca_OG_results", "Bsca_pSD1_results","Bsca_mSD1_results",
                        "Babs_OG_results", "Babs_pSD1_results","Babs_mSD1_results",
                        "convert_list_df", "add_Xnk", "build_stat_df", "split_data_byXnk", 
                        "Babs_calculated_scenarios", "Babs_data")))

# get mean Babs measured values only
Babs <- Babs_data[,1:8]

# get error Babs measured values only
Babs_error <- Babs_data[,-c(2:8, 16:22)]

# add the error to Babs
for (i in 2:length(Babs)) {
  Babs[,i] <- Babs[,i] + Babs_error[,i]
}

# give column names
Aet_WL <- c(370, 470, 520, 590, 660, 880, 950)

for (i in 2:length(Babs)) {
  colnames(Babs)[i] <- sprintf("Babs_meas_%d", Aet_WL[i-1])
}

#############################################################################################################
# ensure calculated and measured data have the same number of observations
#############################################################################################################

Babs$date  <- as.POSIXct(Babs$date ,  format =  "%Y-%m-%d %H:%M:%S") # transform as POSIXct

# join Babs calculated and measured dfs
Babs_comparison <- Babs_calculated_scenarios

for (i in 1:length(Babs_calculated_scenarios)) {
  for (k in 1:length(Babs_calculated_scenarios[[i]])) {
    
    Babs_comparison[[i]][[k]] <- full_join(Babs_calculated_scenarios[[i]][[k]], Babs[,c(1,i+1)])
    
  }
}

# remove Nas
for (i in 1:length(Babs_comparison)) {
  for (k in 1:length(Babs_comparison[[i]])) {
    Babs_comparison[[i]][[k]] <- na.omit(Babs_comparison[[i]][[k]])
  }
}

#############################################################################################################
# add the X, n,  k columns to the dfs 
#############################################################################################################

for (i in 1:length(Babs_comparison)) {
  for (k in 1:length(Babs_comparison[[i]])) {
    Babs_comparison[[i]][[k]] <- add_Xnk(Babs_comparison[[i]][[k]])
  }
}

# convert the X-n-k values as numeric
for (i in 1:length(Babs_comparison)) {
  for (k in 1:length(Babs_comparison[[i]])) {
    Babs_comparison[[i]][[k]]$X <- as.numeric(as.character(Babs_comparison[[i]][[k]]$X))
    Babs_comparison[[i]][[k]]$n <- as.numeric(as.character(Babs_comparison[[i]][[k]]$n))
    Babs_comparison[[i]][[k]]$k <- as.numeric(as.character(Babs_comparison[[i]][[k]]$k))
  }
}

#############################################################################################################
# check the correlation and calculate the statistics
#############################################################################################################

# prepare the data
for (i in 1:length(Babs_comparison)) {
  for (k in 1:length(Babs_comparison[[i]])) {
    colnames(Babs_comparison[[i]][[k]])[3] <- "model"
    colnames(Babs_comparison[[i]][[k]])[4] <- "obs"
  }
}

# remove the combinations with X = 1
rm_X1 <- list()
for (i in 1:length(Babs_comparison)) {
  rm_X1[[i]] <- list()
}

for (i in 1:length(Babs_comparison)){
  for (k in 1:length(Babs_comparison[[i]])){
    if (Babs_comparison[[i]][[k]]$X[1] != 1){
      rm_X1[[i]][[k]] <- Babs_comparison[[i]][[k]]
    }
  }
}

# remove Null terms from the list
for (i in 1:length(rm_X1)){
  rm_X1[[i]] <- rm_X1[[i]][lengths(rm_X1[[i]]) != 0]
}


#############################################################################################################
# check the correlation

test_corr <- rm_X1
for (i in 1:length(rm_X1)) {
  for (k in 1:length(rm_X1[[i]])) {
    test_corr[[i]][[k]] <- summary(lm(model ~ obs, rm_X1[[i]][[k]]))
  }
}

# get the r.squared value
test_R2 <-  rm_X1
for (i in 1:length(rm_X1)) {
  for (k in 1:length(rm_X1[[i]])) {
    test_R2[[i]][[k]] <- round(test_corr[[i]][[k]][["r.squared"]],1)
  }
}

# calculate RMSE
test_RMSE <-  rm_X1
for (i in 1:length(rm_X1)) {
  for (k in 1:length(rm_X1[[i]])) {
    test_RMSE[[i]][[k]] <- rmse(rm_X1[[i]][[k]]$obs, rm_X1[[i]][[k]]$model)/sd(rm_X1[[i]][[k]]$obs)*100
  }
}

#############################################################################################################
# for each WL, create a single df containing a summary of the statistic for each X-n-k scenario
#############################################################################################################

stat_comparison <- list()
for (i in 1:length(rm_X1)) {
  stat_comparison[[i]] <- build_stat_df(Aet_WL[i], rm_X1[[i]], test_RMSE[[i]], test_R2[[i]])
}

Babs_pSD1_results <- stat_comparison

# remove terms from R environment except for:
rm(list=setdiff(ls(), c("Bsca_OG_results", "Bsca_pSD1_results","Bsca_mSD1_results",
                        "Babs_OG_results", "Babs_pSD1_results","Babs_mSD1_results")))

