#############################################################################################################
# set working directory and upload packages
#############################################################################################################

setwd("Z:/Clarissa/Data_Optical_calculation/R_Code/Code_final/ICEdust_samples/Maeli2/Nep_Aet_Mie/original")

library(dplyr)
library(rlist)
library(tibble)
library(ggplot2)

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
# Nephelometer (Nep) data - Scattering coefficients (Bsca)
#############################################################################################################
#############################################################################################################

# import Nep data (1-sec resolution)
Nep <- read.csv("Z:/Clarissa/Data_Optical_calculation/R_Code/Code_final/ICEdust_samples/Maeli2/datasets/tidy_datasets/Maeli2_Nep_data_tidy.csv")

# import correction factors for each wavelength (12-min resolution)
Trunc_correction_700nm <- read.csv("Z:/Clarissa/Data_Optical_calculation/R_Code/Code_final/ICEdust_samples/Maeli2/datasets/Python_output/original/Truncation_corrections/Maeli2_truncation_stat_700nm_original.csv")
Trunc_correction_550nm <- read.csv("Z:/Clarissa/Data_Optical_calculation/R_Code/Code_final/ICEdust_samples/Maeli2/datasets/Python_output/original/Truncation_corrections/Maeli2_truncation_stat_550nm_original.csv")
Trunc_correction_450nm <- read.csv("Z:/Clarissa/Data_Optical_calculation/R_Code/Code_final/ICEdust_samples/Maeli2/datasets/Python_output/original/Truncation_corrections/Maeli2_truncation_stat_450nm_original.csv")

Trunc_correction_700nm$date  <- as.POSIXct(Trunc_correction_700nm$date  ,  format =  "%Y-%m-%d %H:%M:%S")
Trunc_correction_550nm$date   <- as.POSIXct(Trunc_correction_550nm$date  ,  format =  "%Y-%m-%d %H:%M:%S")
Trunc_correction_450nm$date   <- as.POSIXct(Trunc_correction_450nm$date  ,  format =  "%Y-%m-%d %H:%M:%S")

############################################################################################################
# calculate 12-min average 
############################################################################################################

# function to calculate 12-min average by sec
timeAverage_12min <- function(df) {
  # create a complete sequence of date-time 
  date_seq <-  seq(as.POSIXct(df$date[1]), as.POSIXct(df$date[nrow(df)]), by = "sec")
  date_seq <- as.data.frame(date_seq)
  colnames(date_seq) <- "date"
  
  # merge the complete sequence of date-time and the original dataset 
  # the missing date-time will be shown now as NA
  df <- full_join(date_seq, df, by = NULL)
  
  # to calculate 12-min average split the data in #df = dim(df)[1]/720 of 720 rows
  
  # first, calculate the nrow to delete to average constant time intervals
  if ((dim(df)[1]- round(dim(df)[1]/720)*720) < 0) {
    df <-df[1:(dim(df)[1]- (dim(df)[1]- (round(dim(df)[1]/720)-1)*720)),]
  } else {df <-df[1:(dim(df)[1]- (dim(df)[1]- round(dim(df)[1]/720)*720)),]} 
  
  # split the dfs in constant time intervals
  df_split <- split(df,rep(1:round(dim(df)[1]/720),each=720))
  
  # calculate the average
  df_TA <- data.frame(matrix(ncol = ncol(df), nrow = length(df_split)))
  colnames(df_TA) <- colnames(df)
  for (i in 1:length(df_split)) {
    for (j in 2:length(df_split[[i]])) {
      df_TA[i,j] <- mean(df_split[[i]][,j], na.rm = T)
    }
  }
  
  # get the time sequence at 12 min
  df_TA$date <- seq(as.POSIXct(df_split[[1]][1,1]), as.POSIXct(df_split[[length(df_split)]][1,1]), by = "12 min")
  
  return(df_TA)
}

# ensure the date-time column is as "date"
colnames(Nep)[1] <- "date"

# convert Nep$date as POSIXct
Nep$date <- as.POSIXct(Nep$date,  format =  "%Y-%m-%d %H:%M:%S")

# work on the data starting as the model results
Nep_subset <- subset(Nep, Nep$date >= Trunc_correction_700nm$date[1])

# verify the presence of duplicates in the row sequence
Nep_subset <- Nep_subset[-which(duplicated(Nep_subset$date) == T),]

# calculate 12-min average
Nep_TA <- timeAverage_12min(Nep_subset)

# ensure measured and model data have the same time interval
Nep_TA <- subset(Nep_TA,Nep_TA$date <= Trunc_correction_700nm$date[nrow(Trunc_correction_700nm)])

############################################################################################################
# calculate 12-min standard deviation (SD) - step1
############################################################################################################

# function to calculate 12-min SD by sec
timeSD_12min <- function(df) {
  # create a complete sequence of date-time 
  date_seq <-  seq(as.POSIXct(df$date[1]), as.POSIXct(df$date[nrow(df)]), by = "sec")
  date_seq <- as.data.frame(date_seq)
  colnames(date_seq) <- "date"
  
  # merge the complete sequence of date-time and the original dataset 
  # the missing date-time will be shown now as NA
  df <- full_join(date_seq, df, by = NULL)
  
  # to calculate 12-min SD split the data in #df = dim(df)[1]/720 of 720 rows
  
  # first, calculate the nrow to delete in order to obtain constant time intervals
  if ((dim(df)[1]- round(dim(df)[1]/720)*720) < 0) {
    df <-df[1:(dim(df)[1]- (dim(df)[1]- (round(dim(df)[1]/720)-1)*720)),]
  } else {df <-df[1:(dim(df)[1]- (dim(df)[1]- round(dim(df)[1]/720)*720)),]} 
  
  # split the dfs in constant time intervals
  df_split <- split(df,rep(1:round(dim(df)[1]/720),each=720))
  
  # calculate SD
  df_sd <- data.frame(matrix(ncol = ncol(df), nrow = length(df_split)))
  colnames(df_sd) <- colnames(df)
  for (i in 1:length(df_split)) {
    for (j in 2:length(df_split[[i]])) {
      df_sd[i,j] <- sd(df_split[[i]][,j], na.rm = T)
    }
  }
  
  # get the time-averaged sequence
  df_sd$date <- seq(as.POSIXct(df_split[[1]][1,1]), as.POSIXct(df_split[[length(df_split)]][1,1]), by = "12 min")
  
  return(df_sd)
}

# calculate 12-min average SD
Nep_TA_error_step1 <- timeSD_12min(Nep_subset)

# ensure measured and model data have the same time interval
Nep_TA_error_step1 <- subset(Nep_TA_error_step1, Nep_TA_error_step1$date <= Trunc_correction_700nm$date[nrow(Trunc_correction_700nm)])

############################################################################################################
# add to the error further 5% uncertainty due to gas calibration - step 2
############################################################################################################

# calculate gas calibration uncertainty
gas_cal_error <- Nep_TA

for (i in 2:length(Nep_TA)) {
  gas_cal_error[,i] <- (Nep_TA[,i]/100)*5
}

# add the gas calib error to the uncertainty in step1
Nep_TA_error_step2 <- Nep_TA_error_step1

for (i in 2:length(Nep_TA_error_step1)) {
  Nep_TA_error_step2[,i] <- sqrt((Nep_TA_error_step1[,i]^2) + (gas_cal_error[,i]^2))
}

############################################################################################################
# correct Nep data for truncation 
############################################################################################################

# select the columns with Bsca values 
Nep_TA_corrected <- Nep_TA[,c(1:4)]

Nep_TA_corrected$Scatt_blue <- Nep_TA$Scatt_blue*Trunc_correction_450nm$mean
Nep_TA_corrected$Scatt_green <- Nep_TA$Scatt_green*Trunc_correction_550nm$mean
Nep_TA_corrected$Scatt_red <- Nep_TA$Scatt_red*Trunc_correction_700nm$mean

############################################################################################################
# add the uncertainty related to the truncation correction - step3
############################################################################################################

# select the columns with the error on Bsca values 
Nep_TA_error_step2 <- Nep_TA_error_step2[,1:4]

# add the uncertainty related to the truncation correction
# C = A*B, SD(C) = SD(A*B) = C*SQRT(dA/A^2 + dB/B^2) 
Nep_TA_error_step3 <- Nep_TA_error_step2

# 450 nm
Nep_TA_error_step3$Scatt_blue <- Nep_TA_corrected$Scatt_blue*sqrt(((Nep_TA_error_step2$Scatt_blue/Nep_TA$Scatt_blue)^2) + ((Trunc_correction_450nm$sd/Trunc_correction_450nm$mean)^2))

# 550 nm
Nep_TA_error_step3$Scatt_green <- Nep_TA_corrected$Scatt_green*sqrt(((Nep_TA_error_step2$Scatt_green/Nep_TA$Scatt_green)^2) + ((Trunc_correction_550nm$sd/Trunc_correction_550nm$mean)^2))

# 700 nm
Nep_TA_error_step3$Scatt_red <- Nep_TA_corrected$Scatt_red*sqrt(((Nep_TA_error_step2$Scatt_red/Nep_TA$Scatt_red)^2) + ((Trunc_correction_700nm$sd/Trunc_correction_700nm$mean)^2))

############################################################################################################
# predict Bsca at the aethalometer wavelengths (370, 470, 520, 590, 660, 880, 950 nm)
############################################################################################################

# fit the data
# Bsca = A*WL^(-AE)

# Nep wavelength (WL)
x <-  c(450,550,700)

# for each time get Bsca at the 3 WLs
y <- list()

for (i in 1:nrow(Nep_TA_corrected)) {
  y[[i]] <- c(Nep_TA_corrected[i,2], Nep_TA_corrected[i,3], Nep_TA_corrected[i,4])
}

# Weight variable in regression, is a measure of how important an observation is to your model due to different reasons 
# (eg. may be in terms of reliability of measurement or inverse of variance estimate). 
# therefore, some observations may be more important/ weigh higher than others.

# calculate the inverse variance = 1/SD^2
# weights related to Bsca_corrected ("y")
Weights_y <- list()

for (i in 1:nrow(Nep_TA_corrected)) {
  Weights_y[[i]] <- c((1/Nep_TA_error_step3[i,2])^2, (1/Nep_TA_error_step3[i,3])^2, (1/Nep_TA_error_step3[i,4])^2)
}

# build a dataframe x-y
xy <- list()

for (i in 1:length(y)) {
  xy[[i]] <- data.frame(cbind(x,y[[i]]))
  colnames(xy[[i]]) <- c("x", "y")
}

# define the starting a and b value for the power-law fitting
test_ab <- list()
a <- list()
b <- list()

for (i in 1:length(xy)) {
  if (("TRUE" %in% is.na(xy[[i]]$y))|("TRUE" %in% is.na(Weights_y[[i]]))){ # if any term is NA
    test_ab[[i]] <- NA
    a[[i]] <- NA
    b[[i]] <- NA
  }else{
    test_ab[[i]] <- lm(log(y) ~ log(x), data = xy[[i]])
    a[[i]] <- test_ab[[i]][["coefficients"]][1]
    b[[i]] <- test_ab[[i]][["coefficients"]][2]
  }
}

# fit the data using the power law
model_powerLaw <- list()

for (i in 1:length(xy)) {
  if (("TRUE" %in% is.na(xy[[i]]$y))|("TRUE" %in% is.na(Weights_y[[i]]))){ # if any term is NA
    model_powerLaw[[i]] <- NA
  }else{  
    model_powerLaw[[i]] <- nls(y ~ a*x^(b), start = list(a=exp(a[[i]]), b=b[[i]]), data = xy[[i]], control=nls.control(maxiter=1000), weights = Weights_y[[i]])
  }
}

#############################################################################################################
# extrapolate Bsca at 370, 470, 520, 590, 660, 880, 950 nm
#############################################################################################################

# the new WL class must be as dataframe
# Nep WL were also added to estimate the uncertainty related to the predicted values
Aet_WL <- c(370, 470, 520, 590, 660, 880, 950, 450, 550, 700)
new_WL <- list()

for (i in 1:length(Aet_WL)) {
  new_WL[[i]] <- data.frame(x = Aet_WL[i])    
}

# build a function to calculate Bsca at selected WL for all the date-time
predict_Bsca <- function(new_WL) {
  # for each row (date-time) is applied a specific power-law fitting
  # predict Bsca at different time for the selected WL
  Bsca <- list()
  for (i in 1:length(model_powerLaw)) {
    if (("TRUE" %in% is.na(xy[[i]]$y))|("TRUE" %in% is.na(Weights_y[[i]]))){ # if any term is NA
      Bsca[[i]] <- NA
    }else{
      Bsca[[i]] <- predict(model_powerLaw[[i]], new_WL)
    }
  }
  
  # bind the results together
  Bsca_fitted <- unlist(Bsca)
  Bsca_fitted <- as.data.frame(Bsca_fitted)
  return(Bsca_fitted)
}

# find the results at 370, 470, 520, 590, 660, 880, 950 nm
Bsca_fitted_list <- list()
for (i in 1:length(new_WL)) {
  Bsca_fitted_list[[i]] <- predict_Bsca(new_WL[[i]])
  colnames(Bsca_fitted_list[[i]]) <- sprintf("Bsca_%d", Aet_WL[i])
}

# bind the results together
Bsca_fitted_df <- list.cbind(Bsca_fitted_list)

# add the date column 
Bsca_fitted_df <- cbind(Nep_TA$date, Bsca_fitted_df)
Bsca_fitted_df <- data.frame(Bsca_fitted_df)
colnames(Bsca_fitted_df)[1] <- "date"

#############################################################################################################
# calculate the error on Bsca predicted - step4
#############################################################################################################

# Bsca measured vs Bsca calculated at 450, 550 and 700 nm

# for each time get Bsca measured at 450, 550 and 700 nm (Y)
y <- list()

for (i in 1:nrow(Nep_TA_corrected)) {
  y[[i]] <- c(Nep_TA_corrected[i,2], Nep_TA_corrected[i,3], Nep_TA_corrected[i,4])
}


# for each time get Bsca predicted at 450, 550 and 700 nm (X)
x <- list()

for (i in 1:nrow(Bsca_fitted_df)) {
  x[[i]] <- c(Bsca_fitted_df$Bsca_450[i], Bsca_fitted_df$Bsca_550[i], Bsca_fitted_df$Bsca_700[i])
}

# build a dataframe x-y
xy <- list()

for (i in 1:length(y)) {
  xy[[i]] <- data.frame(cbind(x[[i]],y[[i]]))
  colnames(xy[[i]]) <- c("x", "y")
}

# perform the linear regression between Bsca measured and calculated at 450, 550 and 700 nm
# weights are based on the uncertainty on Bsca measured 

Bsca_fitted_error <- list()

for (i in 1:length(xy)) {
  if (("TRUE" %in% is.na(xy[[i]]$y))|("TRUE" %in% is.na(Weights_y[[i]]))){ # if any term is NA
    Bsca_fitted_error[[i]] <- NA
  }else{
    Bsca_fitted_error[[i]] <- lm(y ~ x, data = xy[[i]],  weights = Weights_y[[i]])
  }
}

# calculate the root mean square error (RMSE) on Bsca predicted 
# which corresponds to the standard deviation on Bsca predicted

Bsca_fitted_RMSE <- list()

for (i in 1:length(Bsca_fitted_error)) {
  if ("TRUE" %in% is.na(xy[[i]]$y)){ # if any term is NA
    Bsca_fitted_RMSE[[i]] <- NA
  }else{
    Bsca_fitted_RMSE[[i]] <- sqrt(mean(Bsca_fitted_error[[i]]$residuals^2))
  }
}

Bsca_fitted_RMSE_df <- unlist(Bsca_fitted_RMSE)
Bsca_fitted_RMSE_df <- data.frame("date" = Bsca_fitted_df$date, "RMSE" = Bsca_fitted_RMSE_df)

# remove the columns of Bsca at 450, 550 and 700 nm
Bsca_fitted_df <- Bsca_fitted_df[,-c(9:11)]

#############################################################################################################

# calculate the average systematic error on Bsca corrected at 450, 550 and 700 nm

# calculate Bsca corrected rsd at 450, 550 and 700 nm
Nep_TA_corrected_rsd <- Nep_TA_corrected
for (i in 2:length(Nep_TA_corrected_rsd)) {
  Nep_TA_corrected_rsd[,i] <- Nep_TA_error_step3[,i]/Nep_TA_corrected[,i]*100
}

# calculate the mean of Bsca corrected rsd at 450, 550 and 700 nm
Nep_TA_corrected_rsd$mean_rsd <- rowMeans(Nep_TA_corrected_rsd[,2:4])

# calculate the systematic error on Bsca fitted based on the mean of Bsca corrected rsd at 450, 550 and 700 nm
Bsca_fitted_systematic_error_df <- Bsca_fitted_df

for (i in 2:length(Bsca_fitted_systematic_error_df)) {
  Bsca_fitted_systematic_error_df[,i] <- Bsca_fitted_df[,i]/100*Nep_TA_corrected_rsd$mean_rsd
}

#############################################################################################################

# calculate the total error on Bsca fitted as the sum of the systematic error and RMSE
Nep_TA_error_step4 <- Bsca_fitted_systematic_error_df
for (i in 2:length(Bsca_fitted_systematic_error_df)) {
  Nep_TA_error_step4[,i] <- sqrt(Bsca_fitted_systematic_error_df[,i]^2 + Bsca_fitted_RMSE_df[,2]^2)
}

#############################################################################################################
#############################################################################################################
# Aethalometer (Aet) data - Absorption coefficients (Babs)
#############################################################################################################
#############################################################################################################

# remove terms from R environment except for:
rm(list=setdiff(ls(), c("convert_list_df", 
                        "Bsca_fitted_df", "Nep_TA_error_step4", 
                        "Trunc_correction_450nm", "Trunc_correction_550nm", "Trunc_correction_700nm")))

# import Aet data (2-min resolution)
Aet <- read.csv("Z:/Clarissa/Data_Optical_calculation/R_Code/Code_final/ICEdust_samples/Maeli2/Aet_corrections/original/Maeli2_Aet_Babs_2min_original.csv")

############################################################################################################
# calculate 12-min average
############################################################################################################

# function to calculate 12 min average by min
timeAverage_12min <- function(df) {
  # create a complete sequence of date-time 
  date_seq <-  seq(as.POSIXct(df$date[1]), as.POSIXct(df$date[nrow(df)]), by = "min")
  date_seq <- as.data.frame(date_seq)
  colnames(date_seq) <- "date"
  
  # merge the complete sequence of date-time and the original dataset 
  # the missing date-time will be shown now as NA
  df <- full_join(date_seq, df, by = NULL)
  
  # to calculate 12-min average split the data in #df = dim(df)[1]/12 of 12 rows
  
  # first, calculate the nrow to delete to average constant time intervals
  if ((dim(df)[1]- round(dim(df)[1]/12)*12) < 0) {
    df <-df[1:(dim(df)[1]- (dim(df)[1]- (round(dim(df)[1]/12)-1)*12)),]
  } else {df <-df[1:(dim(df)[1]- (dim(df)[1]- round(dim(df)[1]/12)*12)),]} 
  
  # split the df in constant time intervals
  df_split <- split(df,rep(1:round(dim(df)[1]/12),each=12))
  
  # calculate the average
  df_TA <- data.frame(matrix(ncol = ncol(df), nrow = length(df_split)))
  colnames(df_TA) <- colnames(df)
  for (i in 1:length(df_split)) {
    for (j in 2:length(df_split[[i]])) {
      df_TA[i,j] <- mean(df_split[[i]][,j], na.rm = T)
    }
  }
  
  # get the time sequence at 12 min
  df_TA$date <- seq(as.POSIXct(df_split[[1]][1,1]), as.POSIXct(df_split[[length(df_split)]][1,1]), by = "12 min")
  
  return(df_TA)
}

# ensure the date-time column is as "date"
colnames(Aet)[1] <- "date"

# convert as POSIXct
Aet$date <- as.POSIXct(Aet$date,  format =  "%Y-%m-%d %H:%M:%S")

# work on the data starting as the model results
Aet_subset <- subset(Aet, Aet$date >= Trunc_correction_700nm$date[1])

# verify the presence of duplicates in the row sequence
Aet_subset <- Aet_subset[which(duplicated(Aet_subset$date) == F),]

# calculate 12-min average
Aet_data <- Aet_subset[,1:8]
Aet_TA <- timeAverage_12min(Aet_data)

############################################################################################################
# calculate 12-min average error
############################################################################################################

# function to calculate SD over 12-min intervals considering the original data uncertainty
timeSD_12min <- function(df) {
  # create a complete sequence of date-time 
  date_seq <-  seq(as.POSIXct(df$date[1]), as.POSIXct(df$date[nrow(df)]), by = "min")
  date_seq <- as.data.frame(date_seq)
  colnames(date_seq) <- "date"
  
  # merge the complete sequence of date-time and the original dataset 
  # the missing date-time will be shown now as NA
  df <- full_join(date_seq, df, by = NULL)
  
  # to calculate 12-min SD split the data in #df = dim(df)[1]/12 of 12 rows
  
  # first, calculate the nrow to delete in order to obtain constant time intervals
  if ((dim(df)[1]- round(dim(df)[1]/12)*12) < 0) {
    df <-df[1:(dim(df)[1]- (dim(df)[1]- (round(dim(df)[1]/12)-1)*12)),]
  } else {df <-df[1:(dim(df)[1]- (dim(df)[1]- round(dim(df)[1]/12)*12)),]} 
  
  # split the dfs in constant time intervals
  df_split <- split(df,rep(1:round(dim(df)[1]/12),each=12))
  
  # calculate SD
  df_error <- data.frame(matrix(ncol = ncol(df), nrow = length(df_split)))
  colnames(df_error) <- colnames(df)
  for (i in 1:length(df_split)) {
    for (j in 2:length(df_split[[i]])) {
      df_error[i,j] <- sd(df_split[[i]][,j], na.rm = T)
    }
  }
  
  # get the time sequence at 12 min
  df_error$date <- seq(as.POSIXct(df_split[[1]][1,1]), as.POSIXct(df_split[[length(df_split)]][1,1]), by = "12 min")
  
  return(df_error)
}

# calculate error over 12-min intervals
Aet_TA_error <- timeSD_12min(Aet_data)

############################################################################################################
# calculate the systematic error derived from Aet corrections for each 12-min intervals
############################################################################################################

# get rsd of 2-min Babs data, which is related to Aet corrections
Aet_error_rsd <- Aet_subset[,-c(2:15)]

# calculate the average rsd over 12-min intervals
Aet_TA_Sys_error_rsd <- timeAverage_12min(Aet_error_rsd)

# retrieve the average systemaic error
Aet_TA_Sys_error <- Aet_TA_Sys_error_rsd
colnames(Aet_TA_Sys_error) <- colnames(Aet_TA_error)

for (i in 2:length(Aet_TA)) {
  Aet_TA_Sys_error[,i] <- Aet_TA[,i]/100*Aet_TA_Sys_error_rsd[,i]
}

# sum the SD from 12-min time average and the 12-min averaged systematic error 
Aet_TA_error_total <- Aet_TA_error

for (i in 2:length(Aet_TA_error)) {
  Aet_TA_error_total[,i] <- sqrt(Aet_TA_error[,i]^2 + Aet_TA_Sys_error[,i]^2)
}

# ensure Bsca and Babs have the same time intervals
Bsca_fitted_df <- Bsca_fitted_df[which(Bsca_fitted_df$date %in% Aet_TA$date),]
Nep_TA_error_step4 <- Nep_TA_error_step4[which(Nep_TA_error_step4$date %in% Aet_TA$date),]

# remove the following from the environment
rm("Aet_TA_Sys_error")
rm("Aet_TA_Sys_error_rsd")
rm("Aet_TA_error")
rm("Aet_error_rsd")

#############################################################################################################
#############################################################################################################
# calculate the single scattering albedo (SSA) and extinction coefficient (Bext) at 370, 470, 520, 590, 660, 880, 950 nm
#############################################################################################################
#############################################################################################################

# calculate SSA
SSA <- Aet_TA

for (i in 2:length(Aet_TA)) {
  SSA[,i] <- Bsca_fitted_df[,i]/(Bsca_fitted_df[,i]+Aet_TA[,i])
}

# give column names
Aet_WL <- c(370, 470, 520, 590, 660, 880, 950)

for (i in 2:length(SSA)) {
  colnames(SSA)[i] <- sprintf("SSA_%d", Aet_WL[i-1])
}

# calculate the error related to SSA

# first calculate Bext and the error related to Bext, Bext = Bsca + Babs

# calculate Bext
Bext <- SSA
colnames(Bext) <- c("date", Aet_WL)

for (i in 2:length(Aet_TA)) {
  Bext[,i] <- Aet_TA[,i] + Bsca_fitted_df[,i]
}

# calculate Bext error (Bsca (predicted) error + Babs error)
Bext_error <- Bext

for (i in 2:length(Aet_TA)) {
  Bext_error[,i] <- sqrt(Aet_TA_error_total[,i]^2 + Nep_TA_error_step4[,i]^2)
}

# calculate the error related to SSA (Bext error * Bsca (predicted) error)
SSA_error <- SSA
for (i in 2:length(SSA)) {
  SSA_error[,i] <- SSA[,i]*sqrt((Bext_error[,i]/Bext[,i])^2 + (Nep_TA_error_step4[,i]/Bsca_fitted_df[,i])^2)
}

#####################################################################################################
# before saving out the results
#####################################################################################################

# calculate rsd and tidy up dfs before saving the output

# tidy up Aet data
Aet_TA_rsd <- Aet_TA_error_total

for (i in 2:length(Aet_TA_error_total)) {
  Aet_TA_rsd[,i] <- Aet_TA_error_total[,i]/Aet_TA[,i]*100
  colnames(Aet_TA_rsd)[i] <- sprintf("rsd_%d", Aet_WL[i-1])
  colnames(Aet_TA_error_total)[i] <- sprintf("sd_%d", Aet_WL[i-1])
}

Babs <- cbind(Aet_TA, Aet_TA_error_total[,-1], Aet_TA_rsd[,-1])

# tidy up Nep data
Nep_TA_rsd_step4 <- Bsca_fitted_df

for (i in 2:length(Bsca_fitted_df)) {
  Nep_TA_rsd_step4[,i] <- Nep_TA_error_step4[,i]/Bsca_fitted_df[,i]*100
  colnames(Nep_TA_rsd_step4)[i] <- sprintf("rsd_%d", Aet_WL[i-1])
  colnames(Nep_TA_error_step4)[i] <- sprintf("sd_%d", Aet_WL[i-1])
}

Bsca <- cbind(Bsca_fitted_df, Nep_TA_error_step4[,-1], Nep_TA_rsd_step4[,-1])

# tidy up SSA data
SSA_rsd <- SSA_error

for (i in 2:length(SSA_error)) {
  SSA_rsd[,i] <- SSA_error[,i]/SSA[,i]*100
  colnames(SSA_rsd)[i] <- sprintf("rsd_%d", Aet_WL[i-1])
  colnames(SSA_error)[i] <- sprintf("sd_%d", Aet_WL[i-1])
}

SSA <- cbind(SSA, SSA_error[,-1], SSA_rsd[,-1])

# tidy up Bext data
Bext_rsd <- Bext_error

for (i in 2:length(Bext_error)) {
  Bext_rsd[,i] <- Bext_error[,i]/Bext[,i]*100
  colnames(Bext_rsd)[i] <- sprintf("rsd_%d", Aet_WL[i-1])
  colnames(Bext_error)[i] <- sprintf("sd_%d", Aet_WL[i-1])
  colnames(Bext)[i] <- sprintf("Bext_%d", Aet_WL[i-1])
}

Bext <- cbind(Bext, Bext_error[,-1], Bext_rsd[,-1])

#####################################################################################################
# plot data from 30 min to 2.5 hours after the dust injection peak
#####################################################################################################

# define the time beaks on the x axis
breaks_tidy <- seq(from=as.POSIXct("2019-01-21 11:31","%Y-%m-%d %H:%M", tz="UTC"),
                   to=as.POSIXct("2019-01-21 13:31","%Y-%m-%d %H:%M", tz="UTC"),
                   by="30 min")

# subset Babs data
Babs <- Babs[which(Babs$date >= "2019-01-21 11:31:00" & Babs$date <= "2019-01-21 13:31:00"),]

# plot Babs
ggplot()+
  
  geom_line(data = na.omit(Babs), aes(x= date, y= Babs_370, color = 'df1'), linetype = "dashed")+
  geom_line(data = na.omit(Babs), aes(x= date, y= Babs_470, color = 'df2'), linetype = "dashed")+
  geom_line(data = na.omit(Babs), aes(x= date, y= Babs_520, color = 'df3'), linetype = "dashed")+
  geom_line(data = na.omit(Babs), aes(x= date, y= Babs_590, color = 'df4'), linetype = "dashed")+
  geom_line(data = na.omit(Babs), aes(x= date, y= Babs_660, color = 'df5'), linetype = "dashed")+
  geom_line(data = na.omit(Babs), aes(x= date, y= Babs_880, color = 'df6'), linetype = "dashed")+
  geom_line(data = na.omit(Babs), aes(x= date, y= Babs_950, color = 'df7'), linetype = "dashed")+
  
  geom_point(data = Babs, aes(x= date, y= Babs_370, color = 'df1', shape = 'df1'),  size=3, stroke = 2)+
  geom_point(data = Babs, aes(x= date, y= Babs_470, color = 'df2', shape = 'df2'),  size=3, stroke = 2)+
  geom_point(data = Babs, aes(x= date, y= Babs_520, color = 'df3', shape = 'df3'),  size=3, stroke = 2)+
  geom_point(data = Babs, aes(x= date, y= Babs_590, color = 'df4', shape = 'df4'),  size=3, stroke = 2)+
  geom_point(data = Babs, aes(x= date, y= Babs_660, color = 'df5', shape = 'df5'),  size=3, stroke = 2)+
  geom_point(data = Babs, aes(x= date, y= Babs_880, color = 'df6', shape = 'df6'),  size=3, stroke = 2)+
  geom_point(data = Babs, aes(x= date, y= Babs_950, color = 'df7', shape = 'df7'),  size=3, stroke = 2)+
  
  geom_errorbar(data = na.omit(Babs), aes(x = date, y = Babs_370, ymin = Babs_370 - sd_370, ymax= Babs_370 + sd_370, color = 'df1'), width = 100) +
  geom_errorbar(data = na.omit(Babs), aes(x = date, y = Babs_470, ymin = Babs_470 - sd_470, ymax= Babs_470 + sd_470, color = 'df2'), width = 100) +
  geom_errorbar(data = na.omit(Babs), aes(x = date, y = Babs_520, ymin = Babs_520 - sd_520, ymax= Babs_520 + sd_520, color = 'df3'), width = 100) +
  geom_errorbar(data = na.omit(Babs), aes(x = date, y = Babs_590, ymin = Babs_590 - sd_590, ymax= Babs_590 + sd_590, color = 'df4'), width = 100) +
  geom_errorbar(data = na.omit(Babs), aes(x = date, y = Babs_660, ymin = Babs_660 - sd_660, ymax= Babs_660 + sd_660, color = 'df5'), width = 100) +
  geom_errorbar(data = na.omit(Babs), aes(x = date, y = Babs_880, ymin = Babs_880 - sd_880, ymax= Babs_880 + sd_880, color = 'df6'), width = 100) +
  geom_errorbar(data = na.omit(Babs), aes(x = date, y = Babs_950, ymin = Babs_950 - sd_950, ymax= Babs_950 + sd_950, color = 'df7'), width = 100) +
  
  scale_colour_manual(name = bquote(lambda~(nm)),
                      labels = c("370 nm", "470 nm", "520 nm",
                                 "590 nm", "660 nm", "880 nm", "950 nm"),
                      values = c("df1" = "#00bfff",
                                 "df2" = "#189ed8",
                                 "df3" = "#1c7fb1",
                                 "df4" = "#19618b",
                                 "df5" = "#134467",
                                 "df6" = "#0a2a44",
                                 "df7" = "#001224"))+
  
  scale_shape_manual(name = bquote(lambda~(nm)),
                     labels = c("370 nm", "470 nm", "520 nm",
                                "590 nm", "660 nm", "880 nm", "950 nm"),
                     values = c("df1" = 19,
                                "df2" = 19,
                                "df3" = 19,
                                "df4" = 19,
                                "df5" = 19,
                                "df6" = 19,
                                "df7" = 19))+
  
  scale_x_datetime(breaks = breaks_tidy, labels = seq(30, 150, 30))+
  scale_y_continuous(labels = function(x) format(x, scientific = F), limits = c(0,160))+
  theme_bw()+
  ylab(bquote(beta[abs]~(Mm^-1)))+
  xlab("Time (min)")+
  ggtitle(bquote(bold("Absorption coefficient - "*beta[abs](lambda))))+
  
  theme(plot.title = element_text(size = 18, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text.align = 0,
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        text = element_text(size=18),
        axis.text = element_text(size=18, colour = "black"),
        legend.text = element_text(size=18),
        legend.title = element_text(size=18),
        legend.direction = "vertical",
        legend.position = c(0.8, 0.70),
        legend.key.size = unit(0.5, "cm"))

# subset Bsca data
Bsca <- Bsca[which(Bsca$date >= "2019-01-21 11:31:00" & Bsca$date <= "2019-01-21 13:31:00"),]

# plot Bsca
ggplot()+
  
  geom_line(data = na.omit(Bsca), aes(x= date, y= Bsca_370, color = 'df1'), linetype = "dashed")+
  geom_line(data = na.omit(Bsca), aes(x= date, y= Bsca_470, color = 'df2'), linetype = "dashed")+
  geom_line(data = na.omit(Bsca), aes(x= date, y= Bsca_520, color = 'df3'), linetype = "dashed")+
  geom_line(data = na.omit(Bsca), aes(x= date, y= Bsca_590, color = 'df4'), linetype = "dashed")+
  geom_line(data = na.omit(Bsca), aes(x= date, y= Bsca_660, color = 'df5'), linetype = "dashed")+
  geom_line(data = na.omit(Bsca), aes(x= date, y= Bsca_880, color = 'df6'), linetype = "dashed")+
  geom_line(data = na.omit(Bsca), aes(x= date, y= Bsca_950, color = 'df7'), linetype = "dashed")+
  
  geom_point(data = Bsca, aes(x= date, y= Bsca_370, color = 'df1', shape = 'df1'),  size=3, stroke = 2)+
  geom_point(data = Bsca, aes(x= date, y= Bsca_470, color = 'df2', shape = 'df2'),  size=3, stroke = 2)+
  geom_point(data = Bsca, aes(x= date, y= Bsca_520, color = 'df3', shape = 'df3'),  size=3, stroke = 2)+
  geom_point(data = Bsca, aes(x= date, y= Bsca_590, color = 'df4', shape = 'df4'),  size=3, stroke = 2)+
  geom_point(data = Bsca, aes(x= date, y= Bsca_660, color = 'df5', shape = 'df5'),  size=3, stroke = 2)+
  geom_point(data = Bsca, aes(x= date, y= Bsca_880, color = 'df6', shape = 'df6'),  size=3, stroke = 2)+
  geom_point(data = Bsca, aes(x= date, y= Bsca_950, color = 'df7', shape = 'df7'),  size=3, stroke = 2)+
  
  geom_errorbar(data = na.omit(Bsca), aes(x = date, y = Bsca_370, ymin = Bsca_370 - sd_370, ymax= Bsca_370 + sd_370, color = 'df1'), width = 100) +
  geom_errorbar(data = na.omit(Bsca), aes(x = date, y = Bsca_470, ymin = Bsca_470 - sd_470, ymax= Bsca_470 + sd_470, color = 'df2'), width = 100) +
  geom_errorbar(data = na.omit(Bsca), aes(x = date, y = Bsca_520, ymin = Bsca_520 - sd_520, ymax= Bsca_520 + sd_520, color = 'df3'), width = 100) +
  geom_errorbar(data = na.omit(Bsca), aes(x = date, y = Bsca_590, ymin = Bsca_590 - sd_590, ymax= Bsca_590 + sd_590, color = 'df4'), width = 100) +
  geom_errorbar(data = na.omit(Bsca), aes(x = date, y = Bsca_660, ymin = Bsca_660 - sd_660, ymax= Bsca_660 + sd_660, color = 'df5'), width = 100) +
  geom_errorbar(data = na.omit(Bsca), aes(x = date, y = Bsca_880, ymin = Bsca_880 - sd_880, ymax= Bsca_880 + sd_880, color = 'df6'), width = 100) +
  geom_errorbar(data = na.omit(Bsca), aes(x = date, y = Bsca_950, ymin = Bsca_950 - sd_950, ymax= Bsca_950 + sd_950, color = 'df7'), width = 100) +
  
  
  scale_colour_manual(name = bquote(lambda~(nm)),
                      labels = c("370 nm", "470 nm", "520 nm",
                                 "590 nm", "660 nm", "880 nm", "950 nm"),
                      values = c("df1" = "#00bfff",
                                 "df2" = "#189ed8",
                                 "df3" = "#1c7fb1",
                                 "df4" = "#19618b",
                                 "df5" = "#134467",
                                 "df6" = "#0a2a44",
                                 "df7" = "#001224"))+
  
  scale_shape_manual(name = bquote(lambda~(nm)),
                     labels = c("370 nm", "470 nm", "520 nm",
                                "590 nm", "660 nm", "880 nm", "950 nm"),
                     values = c("df1" = 19,
                                "df2" = 19,
                                "df3" = 19,
                                "df4" = 19,
                                "df5" = 19,
                                "df6" = 19,
                                "df7" = 19))+
  
  scale_x_datetime(breaks = breaks_tidy, labels = seq(30, 150, 30))+
  scale_y_continuous(labels = function(x) format(x, scientific = F), limits = c(0,1200))+
  theme_bw()+
  ylab(bquote(beta[sca]~(Mm^-1)))+
  xlab("Time (min)")+
  ggtitle(bquote(bold("Scattering coefficient - "*beta[sca](lambda))))+
  
  theme(plot.title = element_text(size = 18, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text.align = 0,
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        text = element_text(size=18),
        axis.text = element_text(size=18, colour = "black"),
        legend.text = element_text(size=18),
        legend.title = element_text(size=18),
        legend.direction = "vertical",
        legend.position = c(0.8, 0.70),
        legend.key.size = unit(0.5, "cm"))

# subset Bext data
Bext <- Bext[which(Bext$date >= "2019-01-21 11:31:00" & Bext$date <= "2019-01-21 13:31:00"),]

# plot Bext
ggplot()+
  
  geom_line(data = na.omit(Bext), aes(x= date, y= Bext_370, color = 'df1'), linetype = "dashed")+
  geom_line(data = na.omit(Bext), aes(x= date, y= Bext_470, color = 'df2'), linetype = "dashed")+
  geom_line(data = na.omit(Bext), aes(x= date, y= Bext_520, color = 'df3'), linetype = "dashed")+
  geom_line(data = na.omit(Bext), aes(x= date, y= Bext_590, color = 'df4'), linetype = "dashed")+
  geom_line(data = na.omit(Bext), aes(x= date, y= Bext_660, color = 'df5'), linetype = "dashed")+
  geom_line(data = na.omit(Bext), aes(x= date, y= Bext_880, color = 'df6'), linetype = "dashed")+
  geom_line(data = na.omit(Bext), aes(x= date, y= Bext_950, color = 'df7'), linetype = "dashed")+
  
  geom_point(data = Bext, aes(x= date, y= Bext_370, color = 'df1', shape = 'df1'),  size=3, stroke = 2)+
  geom_point(data = Bext, aes(x= date, y= Bext_470, color = 'df2', shape = 'df2'),  size=3, stroke = 2)+
  geom_point(data = Bext, aes(x= date, y= Bext_520, color = 'df3', shape = 'df3'),  size=3, stroke = 2)+
  geom_point(data = Bext, aes(x= date, y= Bext_590, color = 'df4', shape = 'df4'),  size=3, stroke = 2)+
  geom_point(data = Bext, aes(x= date, y= Bext_660, color = 'df5', shape = 'df5'),  size=3, stroke = 2)+
  geom_point(data = Bext, aes(x= date, y= Bext_880, color = 'df6', shape = 'df6'),  size=3, stroke = 2)+
  geom_point(data = Bext, aes(x= date, y= Bext_950, color = 'df7', shape = 'df7'),  size=3, stroke = 2)+
  
  geom_errorbar(data = na.omit(Bext), aes(x = date, y = Bext_370, ymin = Bext_370 - sd_370, ymax= Bext_370 + sd_370, color = 'df1'), width = 100) +
  geom_errorbar(data = na.omit(Bext), aes(x = date, y = Bext_470, ymin = Bext_470 - sd_470, ymax= Bext_470 + sd_470, color = 'df2'), width = 100) +
  geom_errorbar(data = na.omit(Bext), aes(x = date, y = Bext_520, ymin = Bext_520 - sd_520, ymax= Bext_520 + sd_520, color = 'df3'), width = 100) +
  geom_errorbar(data = na.omit(Bext), aes(x = date, y = Bext_590, ymin = Bext_590 - sd_590, ymax= Bext_590 + sd_590, color = 'df4'), width = 100) +
  geom_errorbar(data = na.omit(Bext), aes(x = date, y = Bext_660, ymin = Bext_660 - sd_660, ymax= Bext_660 + sd_660, color = 'df5'), width = 100) +
  geom_errorbar(data = na.omit(Bext), aes(x = date, y = Bext_880, ymin = Bext_880 - sd_880, ymax= Bext_880 + sd_880, color = 'df6'), width = 100) +
  geom_errorbar(data = na.omit(Bext), aes(x = date, y = Bext_950, ymin = Bext_950 - sd_950, ymax= Bext_950 + sd_950, color = 'df7'), width = 100) +
  
  
  scale_colour_manual(name = bquote(lambda~(nm)),
                      labels = c("370 nm", "470 nm", "520 nm",
                                 "590 nm", "660 nm", "880 nm", "950 nm"),
                      values = c("df1" = "#00bfff",
                                 "df2" = "#189ed8",
                                 "df3" = "#1c7fb1",
                                 "df4" = "#19618b",
                                 "df5" = "#134467",
                                 "df6" = "#0a2a44",
                                 "df7" = "#001224"))+
  
  scale_shape_manual(name = bquote(lambda~(nm)),
                     labels = c("370 nm", "470 nm", "520 nm",
                                "590 nm", "660 nm", "880 nm", "950 nm"),
                     values = c("df1" = 19,
                                "df2" = 19,
                                "df3" = 19,
                                "df4" = 19,
                                "df5" = 19,
                                "df6" = 19,
                                "df7" = 19))+
  
  scale_x_datetime(breaks = breaks_tidy, labels = seq(30, 150, 30))+
  scale_y_continuous(labels = function(x) format(x, scientific = F), limits = c(0,1200))+
  theme_bw()+
  ylab(bquote(beta[ext]~(Mm^-1)))+
  xlab("Time (min)")+
  ggtitle(bquote(bold("Extinction coefficient - "*beta[ext](lambda))))+
  
  theme(plot.title = element_text(size = 18, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text.align = 0,
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        text = element_text(size=18),
        axis.text = element_text(size=18, colour = "black"),
        legend.text = element_text(size=18),
        legend.title = element_text(size=18),
        legend.direction = "vertical",
        legend.position = c(0.8, 0.70),
        legend.key.size = unit(0.5, "cm"))

# subset SSA data
SSA <-SSA[which(SSA$date >= "2019-01-21 11:31:00" & SSA$date <= "2019-01-21 13:31:00"),]

# plot SSA
ggplot()+
  
  geom_line(data = na.omit(SSA), aes(x= date, y= SSA_370, color = 'df1'), linetype = "dashed")+
  geom_line(data = na.omit(SSA), aes(x= date, y= SSA_470, color = 'df2'), linetype = "dashed")+
  geom_line(data = na.omit(SSA), aes(x= date, y= SSA_520, color = 'df3'), linetype = "dashed")+
  geom_line(data = na.omit(SSA), aes(x= date, y= SSA_590, color = 'df4'), linetype = "dashed")+
  geom_line(data = na.omit(SSA), aes(x= date, y= SSA_660, color = 'df5'), linetype = "dashed")+
  geom_line(data = na.omit(SSA), aes(x= date, y= SSA_880, color = 'df6'), linetype = "dashed")+
  geom_line(data = na.omit(SSA), aes(x= date, y= SSA_950, color = 'df7'), linetype = "dashed")+
  
  geom_point(data = SSA, aes(x= date, y= SSA_370, color = 'df1', shape = 'df1'),  size=3, stroke = 2)+
  geom_point(data = SSA, aes(x= date, y= SSA_470, color = 'df2', shape = 'df2'),  size=3, stroke = 2)+
  geom_point(data = SSA, aes(x= date, y= SSA_520, color = 'df3', shape = 'df3'),  size=3, stroke = 2)+
  geom_point(data = SSA, aes(x= date, y= SSA_590, color = 'df4', shape = 'df4'),  size=3, stroke = 2)+
  geom_point(data = SSA, aes(x= date, y= SSA_660, color = 'df5', shape = 'df5'),  size=3, stroke = 2)+
  geom_point(data = SSA, aes(x= date, y= SSA_880, color = 'df6', shape = 'df6'),  size=3, stroke = 2)+
  geom_point(data = SSA, aes(x= date, y= SSA_950, color = 'df7', shape = 'df7'),  size=3, stroke = 2)+
  
  geom_errorbar(data = na.omit(SSA), aes(x = date, y = SSA_370, ymin = SSA_370 - sd_370, ymax= SSA_370 + sd_370, color = 'df1'), width = 100) +
  geom_errorbar(data = na.omit(SSA), aes(x = date, y = SSA_470, ymin = SSA_470 - sd_470, ymax= SSA_470 + sd_470, color = 'df2'), width = 100) +
  geom_errorbar(data = na.omit(SSA), aes(x = date, y = SSA_520, ymin = SSA_520 - sd_520, ymax= SSA_520 + sd_520, color = 'df3'), width = 100) +
  geom_errorbar(data = na.omit(SSA), aes(x = date, y = SSA_590, ymin = SSA_590 - sd_590, ymax= SSA_590 + sd_590, color = 'df4'), width = 100) +
  geom_errorbar(data = na.omit(SSA), aes(x = date, y = SSA_660, ymin = SSA_660 - sd_660, ymax= SSA_660 + sd_660, color = 'df5'), width = 100) +
  geom_errorbar(data = na.omit(SSA), aes(x = date, y = SSA_880, ymin = SSA_880 - sd_880, ymax= SSA_880 + sd_880, color = 'df6'), width = 100) +
  geom_errorbar(data = na.omit(SSA), aes(x = date, y = SSA_950, ymin = SSA_950 - sd_950, ymax= SSA_950 + sd_950, color = 'df7'), width = 100) +
  
  scale_colour_manual(name = bquote(lambda~(nm)),
                      labels = c("370 nm", "470 nm", "520 nm",
                                 "590 nm", "660 nm", "880 nm", "950 nm"),
                      values = c("df1" = "#00bfff",
                                 "df2" = "#189ed8",
                                 "df3" = "#1c7fb1",
                                 "df4" = "#19618b",
                                 "df5" = "#134467",
                                 "df6" = "#0a2a44",
                                 "df7" = "#001224"))+
  
  scale_shape_manual(name = bquote(lambda~(nm)),
                     labels = c("370 nm", "470 nm", "520 nm",
                                "590 nm", "660 nm", "880 nm", "950 nm"),
                     values = c("df1" = 19,
                                "df2" = 19,
                                "df3" = 19,
                                "df4" = 19,
                                "df5" = 19,
                                "df6" = 19,
                                "df7" = 19))+
  
  scale_x_datetime(breaks = breaks_tidy, labels = seq(30, 150, 30))+
  scale_y_continuous(labels = function(x) format(x, scientific = F), limits = c(0,1.3))+
  theme_bw()+
  ylab(bquote(SSA))+
  xlab("Time (min)")+
  ggtitle(bquote(bold("Single scattering albedo - "*SSA(lambda))))+
  
  theme(plot.title = element_text(size = 18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text.align = 0,
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        text = element_text(size=18),
        axis.text = element_text(size=18, colour = "black"),
        legend.text = element_text(size=18),
        legend.title = element_text(size=18),
        legend.direction = "vertical",
        legend.position = c(0.8, 0.3),
        legend.key.size = unit(0.5, "cm"))
