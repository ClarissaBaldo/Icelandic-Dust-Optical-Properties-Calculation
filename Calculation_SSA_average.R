#############################################################################################################
# set working directory and upload packages
#############################################################################################################

# get the working directory
setwd("Z:/Clarissa/Data_Optical_calculation/R_Code/Code_final/ICEdust_samples/Maeli2/SSA_average/original")

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
# import measured Mie coefficients 
#############################################################################################################

# import processed data as 12 min avergae
Maeli2_Babs_data <- read.csv("Z:/Clarissa/Data_Optical_calculation/R_Code/Code_final/ICEdust_samples/Maeli2/Nep_Aet_Mie/original/data/Maeli2_Babs_12min_original.csv")
Maeli2_Bsca_data <- read.csv("Z:/Clarissa/Data_Optical_calculation/R_Code/Code_final/ICEdust_samples/Maeli2/Nep_Aet_Mie/original/data/Maeli2_Bsca_fitted_12min_original.csv")

Maeli2_Babs_data$date  <- as.POSIXct(Maeli2_Babs_data$date  ,  format =  "%Y-%m-%d %H:%M:%S")
Maeli2_Bsca_data$date  <- as.POSIXct(Maeli2_Bsca_data$date  ,  format =  "%Y-%m-%d %H:%M:%S")

Maeli2_Babs_data <- Maeli2_Babs_data[which(Maeli2_Babs_data$date >= "2019-01-21 11:31:00" & Maeli2_Babs_data$date <= "2019-01-21 13:31:00"),]
Maeli2_Bsca_data <- Maeli2_Bsca_data[which(Maeli2_Bsca_data$date >= "2019-01-21 11:31:00" & Maeli2_Bsca_data$date <= "2019-01-21 13:31:00"),]


# split data and sd values
Maeli2_Babs <- Maeli2_Babs_data[,1:8]
Maeli2_Babs_error <- Maeli2_Babs_data[,c(1,9:15)]

Maeli2_Bsca <- Maeli2_Bsca_data[,1:8]
Maeli2_Bsca_error <- Maeli2_Bsca_data[,c(1,9:15)]

# prepare data for fitting
data_fitting <- list()

for (i in 2:length(Maeli2_Babs)) {
  data_fitting[[i]] <- data.frame("date" = Maeli2_Babs[,1], 
                                  "Babs" = Maeli2_Babs[,i], 
                                  "Bsca" = Maeli2_Bsca[,i], 
                                  "Babs_sd" = Maeli2_Babs_error[,i], 
                                  "Bsca_sd" = Maeli2_Bsca_error[,i])
}

# drop null terms from the list
data_fitting <- data_fitting[lengths(data_fitting) != 0]

# drop all rows with na values
for (i in 1:length(data_fitting)) {
  data_fitting[[i]] <- na.omit(data_fitting[[i]])
}

###########################################################################################################
###########################################################################################################
# Reduced Major Axis Regression (SMA)
###########################################################################################################
###########################################################################################################

# packages needed
library(lmodel2)
library(dplyr)

# perform SMA regression for comparison with deming results
data_SMA <- data_fitting

# lmodel2(y ~ x)
# both x and y have errors, but these are assumed to be constant and not input
# this is why we use SMA, but you must set x for the obs you trust more (x for the so-called "ground-truth")

# perform SMA regression
SMA_fit <- list()
for (i in 1:length(data_SMA)) {
  SMA_fit [[i]] <- lmodel2(Bsca ~ Babs, data = data_SMA[[i]])
}

# save the main results of the SMA regression
SSA_SMA <- list()
for (i in 1:length(SMA_fit)) {
  SSA_SMA[[i]] <- cbind(SMA_fit[[i]][["regression.results"]][,-c(4:5)], SMA_fit[[i]][["confidence.intervals"]][,-1], SMA_fit[[i]][["rsquare"]])
}

write.xlsx(SSA_SMA, sheetName = Aet_WL,
           "Z:/Clarissa/Data_Optical_calculation/R_Code/Code_final/ICEdust_samples/Maeli2/SSA_average/original/Maeli2_SMA_original.xlsx", row.names = F)

###########################################################################################################
# build a df with the exp-average SSA results from SMA regression
###########################################################################################################

# calculate sd of exp-average SSA from SMA regression using the estimated confidence intervals
calculate_SMA_sd <- list()

for (i in 1:length(SSA_SMA)) {
  calculate_SMA_sd[[i]] <- data.frame(test1 = abs(SSA_SMA[[i]]$`2.5%-Slope`[3] - SSA_SMA[[i]]$Slope[3]),
                                  test2 = abs(SSA_SMA[[i]]$`97.5%-Slope`[3] - SSA_SMA[[i]]$Slope[3]))
}

# select the values with the largest sd
SSA_SMA_sd <- c()
for (i in 1:length(calculate_SMA_sd)) {
  if (calculate_SMA_sd[[i]]$test1 > calculate_SMA_sd[[i]]$test2) {
    SSA_SMA_sd[i] <- calculate_SMA_sd[[i]]$test1
  } else  {
    SSA_SMA_sd[i] <- calculate_SMA_sd[[i]]$test2
  }
}

# get the slope of the SMA regression
SSA_SMA_slope <- c()
for (i in 1:length(SSA_SMA)) {
  SSA_SMA_slope[i] <- SSA_SMA[[i]]$Slope[3]
}

# output exp-average SSA from SMA regression
SSA_SMA_data <- data.frame("WL" = Aet_WL,
                           "m" = SSA_SMA_slope, 
                           "m_sd" = SSA_SMA_sd,
                           "m_rsd" = SSA_SMA_sd/SSA_SMA_slope*100)

SSA_average_SMA <- SSA_SMA_data
SSA_average_SMA$SSA <- (1 + 1/SSA_average_SMA$m)^(-1)
SSA_average_SMA$SSA_sd <- SSA_average_SMA$SSA*SSA_average_SMA$m_sd/(SSA_average_SMA$m)
SSA_average_SMA$SSA_rsd <- SSA_average_SMA$SSA_sd/SSA_average_SMA$SSA*100


# save the results
write.csv(SSA_average_SMA, "Z:/Clarissa/Data_Optical_calculation/R_Code/Code_final/ICEdust_samples/Maeli2/SSA_average/original/Maeli2_SSA_exp_average_SMA_original.csv", row.names = F)

