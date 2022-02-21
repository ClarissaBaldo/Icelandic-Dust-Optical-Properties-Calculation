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
# import calculated Mie coefficients
#############################################################################################################

# set wd
setwd("Z:/Clarissa/Data_Optical_calculation/R_Code/Code_final/ICEdust_samples/Maeli2/datasets/Python_output/original/Mie_coef")

# create a list using a selected pattern
filenames <- list.files(pattern="*.csv", full.names=TRUE)

# import selected files
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
# get Bext and Bsca calculated based on size distribution data
#############################################################################################################

# set new wd
setwd("Z:/Clarissa/Data_Optical_calculation/R_Code/Code_final/ICEdust_samples/Maeli2/Nep_Aet_Mie/original/")

Mie_coeff_calculated <- Datafile

# select Xnk, date, Bext and Bsca columns in dfs
for (i in 1:length(Datafile)){
  Mie_coeff_calculated[[i]] <- Mie_coeff_calculated[[i]][, c(2:5)]
}

# add WL column
WL <- c(370, 450, 470, 520, 550, 590, 660, 700, 880, 950)

for (i in 1:length(Mie_coeff_calculated)) {
  Mie_coeff_calculated[[i]] <- add_column(Mie_coeff_calculated[[i]], "WL" = WL[[i]], .before = 1)
}

# remove data at 450, 550, and 700 nm
rm_WL <- c(450, 550, 700)

rm_data <- list()
for (i in 1:length(Mie_coeff_calculated)) {
  if ((Mie_coeff_calculated[[i]][1,1] %in% rm_WL)==F){
    rm_data[[i]] <- Mie_coeff_calculated[[i]]}
  else {
    rm_data[[i]] <- NULL
  }
}

# remove Null terms from the list
Mie_coeff_calculated <- rm_data[lengths(rm_data) != 0]

#############################################################################################################
# ensure that measured and modelled Mie coefficients have the same time length
#############################################################################################################

# import measured SSA
SSA_data <- read.csv("Z:/Clarissa/Data_Optical_calculation/R_Code/Code_final/ICEdust_samples/Maeli2/Nep_Aet_Mie/original/data/Maeli2_SSA_12min_original.csv")

SSA_data$date  <- as.POSIXct(SSA_data$date ,  format =  "%Y-%m-%d %H:%M:%S") # transform as POSIXct

# subset the data from 30 min to 2.5h after the dust injection peak
SSA_data <- SSA_data[which(SSA_data$date >= "2019-01-21 11:31:00" & SSA_data$date <= "2019-01-21 13:31:00"),]

# measured and modelled Mie coefficients must have the same time length
for (i in 1:length(Mie_coeff_calculated)) {
  Mie_coeff_calculated[[i]] <- Mie_coeff_calculated[[i]][which(Mie_coeff_calculated[[i]]$date %in% SSA_data$date),]
}

#############################################################################################################
# split data based on Xnk
#############################################################################################################

# create a function to split data based on Xnk
split_data_byXnk <- function(df) {
  
  # find the Xnk range
  Xnk_range <- unique(df$Xnk)
  
  # create a list with calculated Bext/Bsca for each scenario
  WL_scenario <- list()
  for (i in 1:length(Xnk_range)) {
    WL_scenario[[i]] <- df[which(df$Xnk %in% Xnk_range[i] == TRUE),]
    WL_scenario[[i]]$Xnk <- as.character(WL_scenario[[i]]$Xnk) # convert from 1 level factor to charachter
  }
  
  return(WL_scenario) 
  
}

# apply the function to all dfs in the list
Mie_coeff_calculated_scenarios <- list()

for (i in 1:length(Mie_coeff_calculated)) {
  Mie_coeff_calculated_scenarios[[i]] <- split_data_byXnk(Mie_coeff_calculated[[i]])
}

#############################################################################################################
# calculate experiment-averaged SSA modelled
#############################################################################################################

# packages needed
library(lmodel2)
library(dplyr)

# lmodel2(y ~ x)
# both x and y have errors, but these are assumed to be constant and not input
# this is why we use SMA, but you must set x for the obs you trust more (x for the so-called "ground-truth")

# to calculate the experiment-averaged SSA modelled, perform SMA regression Bsca vs Bext
test <- lmodel2(Bsca ~ Bext, data = Mie_coeff_calculated_scenarios[[1]][[1]])

SSA_SMA <- Mie_coeff_calculated_scenarios
for (i in 1:length(Mie_coeff_calculated_scenarios)) {
  for (j in 1:length(Mie_coeff_calculated_scenarios[[i]])) {
    SSA_SMA[[i]][[j]] <- lmodel2(Bsca ~ Bext, data = Mie_coeff_calculated_scenarios[[i]][[j]])
  }
}

###########################################################################################################
# build a df with the experiment-averaged SSA modelled results from SMA regression
###########################################################################################################

# calculate sd of experiment-averaged SSA modelled from SMA regression using the estimated confidence intervals
calculate_SMA_sd <- Mie_coeff_calculated_scenarios
for (i in 1:length(SSA_SMA)) {
  for (j in 1:length(SSA_SMA[[i]])) {
  calculate_SMA_sd[[i]][[j]] <- data.frame(test1 = abs(SSA_SMA[[i]][[j]][["confidence.intervals"]][["2.5%-Slope"]][3] - SSA_SMA[[i]][[j]][["regression.results"]][["Slope"]][3]),
                                      test2 = abs(SSA_SMA[[i]][[j]][["confidence.intervals"]][["97.5%-Slope"]][3] - SSA_SMA[[i]][[j]][["regression.results"]][["Slope"]][3]))
  }
}

# select the values with the largest sd
SSA_SMA_sd <- calculate_SMA_sd
for (i in 1:length(calculate_SMA_sd)) {
  for (j in 1:length(calculate_SMA_sd[[i]])) {
    if (calculate_SMA_sd[[i]][[j]]$test1 > calculate_SMA_sd[[i]][[j]]$test2) {
      SSA_SMA_sd[[i]][[j]] <- calculate_SMA_sd[[i]][[j]]$test1
    } else  {
      SSA_SMA_sd[[i]][[j]] <- calculate_SMA_sd[[i]][[j]]$test2
    }
  }
}

# get the slope of SMA regression line
SSA_SMA_slope <- Mie_coeff_calculated_scenarios
for (i in 1:length(SSA_SMA)) {
  for (j in 1:length(SSA_SMA[[i]])) {
  SSA_SMA_slope[[i]][[j]] <- SSA_SMA[[i]][[j]][["regression.results"]][["Slope"]][3]
  }
}

# add WL column
Aet_WL <- c(370, 470, 520, 590, 660, 880, 950)

# output experiment-averaged SSA modelled from SMA regression
SSA_SMA_data <- Mie_coeff_calculated_scenarios

for (i in 1:length(SSA_SMA_sd)) {
  for (j in 1:length(SSA_SMA_sd[[i]])) {
    
    SSA_SMA_data[[i]][[j]] <- data.frame("WL" = Aet_WL[i],
                                         "Xnk" = Mie_coeff_calculated_scenarios[[i]][[j]][1,2],
                                         "SSA" = SSA_SMA_slope[[i]][[j]], 
                                         "SSA_sd" = SSA_SMA_sd[[i]][[j]],
                                         "SSA_rsd" = SSA_SMA_sd[[i]][[j]]/SSA_SMA_slope[[i]][[j]]*100)
  }
}

# combine the data by WL
SSA_Xnk_WL <- list()
for (i in 1:length(SSA_SMA_data)) {
  SSA_Xnk_WL[[i]] <- list.rbind(SSA_SMA_data[[i]])
}

#############################################################################################################
# add X, n, k columns to dfs
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

for (i in 1:length(SSA_Xnk_WL)) {
  SSA_Xnk_WL[[i]] <- add_Xnk(SSA_Xnk_WL[[i]])
}

# convert X-n-k values as numeric
for (i in 1:length(SSA_Xnk_WL)) {
    SSA_Xnk_WL[[i]]$X <- as.numeric(as.character(SSA_Xnk_WL[[i]]$X))
    SSA_Xnk_WL[[i]]$n <- as.numeric(as.character(SSA_Xnk_WL[[i]]$n))
    SSA_Xnk_WL[[i]]$k <- as.numeric(as.character(SSA_Xnk_WL[[i]]$k))
}

# remove data X=1
for (i in 1:length(SSA_Xnk_WL)) {
  SSA_Xnk_WL[[i]] <- SSA_Xnk_WL[[i]][which(SSA_Xnk_WL[[i]]$X != 1),]
}

# remove terms from R environment except for:
rm(list=setdiff(ls(), c("Bsca_OG_results", "Bsca_pSD1_results","Bsca_mSD1_results",
                        "Babs_OG_results", "Babs_pSD1_results","Babs_mSD1_results",
                        "SSA_Xnk_WL")))

