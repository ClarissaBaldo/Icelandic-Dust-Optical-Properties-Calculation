#############################################################################################################
# set working directory and upload packages
#############################################################################################################
# get the working directory
setwd("Z:/Clarissa/Data_Optical_calculation/R_Code/Code_final/ICEdust_samples/Maeli2/Nep_Aet_Mie/original")

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
# Nep measured Mie coefficients 
#############################################################################################################
#############################################################################################################

# import Nep data
Nep <- read.csv("Z:/Clarissa/Data_Optical_calculation/R_Code/Code_final/ICEdust_samples/Maeli2/datasets/tidy_datasets/Maeli2_Nep_data_tidy.csv")

# import correction factor for each wavelength
Trunc_correction_700nm <- read.csv("Z:/Clarissa/Data_Optical_calculation/R_Code/Code_final/ICEdust_samples/Maeli2/datasets/Python_output/original/Truncation_corrections/Maeli2_truncation_stat_700nm_original.csv")
Trunc_correction_550nm <- read.csv("Z:/Clarissa/Data_Optical_calculation/R_Code/Code_final/ICEdust_samples/Maeli2/datasets/Python_output/original/Truncation_corrections/Maeli2_truncation_stat_550nm_original.csv")
Trunc_correction_450nm <- read.csv("Z:/Clarissa/Data_Optical_calculation/R_Code/Code_final/ICEdust_samples/Maeli2/datasets/Python_output/original/Truncation_corrections/Maeli2_truncation_stat_450nm_original.csv")

Trunc_correction_700nm$date  <- as.POSIXct(Trunc_correction_700nm$date  ,  format =  "%Y-%m-%d %H:%M:%S")
Trunc_correction_550nm$date   <- as.POSIXct(Trunc_correction_550nm$date  ,  format =  "%Y-%m-%d %H:%M:%S")
Trunc_correction_450nm$date   <- as.POSIXct(Trunc_correction_450nm$date  ,  format =  "%Y-%m-%d %H:%M:%S")

############################################################################################################
# calculate 12-min average 
############################################################################################################

# function to calculate 12 min average by sec
timeAverage_12min <- function(df) {
  # create a complete sequence of date-time 
  date_seq <-  seq(as.POSIXct(df$date[1]), as.POSIXct(df$date[nrow(df)]), by = "sec")
  date_seq <- as.data.frame(date_seq)
  colnames(date_seq) <- "date"
  
  # merge the complete sequence of date-time and the original dataset 
  # the missing date-time will be shown now as NA
  df <- full_join(date_seq, df, by = NULL)
  
  # to calculate 12 minute average split the data in #df = dim(df)[1]/720 of 720 rows
  
  # first, calculate the nrow to delete to average constant time intervals
  if ((dim(df)[1]- round(dim(df)[1]/720)*720) < 0) {
    df <-df[1:(dim(df)[1]- (dim(df)[1]- (round(dim(df)[1]/720)-1)*720)),]
  } else {df <-df[1:(dim(df)[1]- (dim(df)[1]- round(dim(df)[1]/720)*720)),]} 
  
  # split the df in constant time intervals
  df_split <- split(df,rep(1:round(dim(df)[1]/720),each=720))
  
  # calculate the average
  df_TA <- data.frame(matrix(ncol = ncol(df), nrow = length(df_split)))
  colnames(df_TA) <- colnames(df)
  for (i in 1:length(df_split)) {
    for (j in 2:length(df_split[[i]])) {
      df_TA[i,j] <- mean(df_split[[i]][,j], na.rm = T)
    }
  }
  
  # get the time-averaged sequence
  df_TA$date <- seq(as.POSIXct(df_split[[1]][1,1]), as.POSIXct(df_split[[length(df_split)]][1,1]), by = "12 min")
  
  return(df_TA)
}

# ensure the date-time column is as "date"
colnames(Nep)[1] <- "date"

# convert Nep$date as POSIXct
Nep$date <- as.POSIXct(Nep$date,  format =  "%Y-%m-%d %H:%M:%S")

# work on the data starting as in the model results
Nep_subset <- subset(Nep, Nep$date >= Trunc_correction_700nm$date[1])

# verify the presence of duplicates in the row sequence
Nep_subset <- Nep_subset[-which(duplicated(Nep_subset$date) == T),]

# calculatee 12 min average
Nep_TA <- timeAverage_12min(Nep_subset)

# ensure measured and model data have the same time interval
Nep_TA <- subset(Nep_TA,Nep_TA$date <= Trunc_correction_700nm$date[nrow(Trunc_correction_700nm)])

############################################################################################################
# calculate 12-min STDV - step1
############################################################################################################

# function to calculate 12 min sd by sec
timeSTDV_12min <- function(df) {
  # create a complete sequence of date-time 
  date_seq <-  seq(as.POSIXct(df$date[1]), as.POSIXct(df$date[nrow(df)]), by = "sec")
  date_seq <- as.data.frame(date_seq)
  colnames(date_seq) <- "date"
  
  # merge the complete sequence of date-time and the original dataset 
  # the missing date-time will be shown now as NA
  df <- full_join(date_seq, df, by = NULL)
  
  # to calculate 12 minute average split the data in #df = dim(df)[1]/720 of 720 rows
  
  # first, calculate the nrow to delete to average constant time intervals
  if ((dim(df)[1]- round(dim(df)[1]/720)*720) < 0) {
    df <-df[1:(dim(df)[1]- (dim(df)[1]- (round(dim(df)[1]/720)-1)*720)),]
  } else {df <-df[1:(dim(df)[1]- (dim(df)[1]- round(dim(df)[1]/720)*720)),]} 
  
  # split the df in constant time intervals
  df_split <- split(df,rep(1:round(dim(df)[1]/720),each=720))
  
  # calculate sd
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

# calculate 12 min average STDV
Nep_TA_error_step1 <- timeSTDV_12min(Nep_subset)

# ensure measured and model data have the same time interval
Nep_TA_error_step1 <- subset(Nep_TA_error_step1, Nep_TA_error_step1$date <= Trunc_correction_700nm$date[nrow(Trunc_correction_700nm)])

############################################################################################################
# add to STDV further 5% uncertainty due to gas calibration - step 2
############################################################################################################

# calculate gas calibration STDV
gas_cal_error <- Nep_TA

for (i in 2:length(Nep_TA)) {
  
  gas_cal_error[,i] <- (Nep_TA[,i]/100)*5
}

# add the gas calib error to the 12 min STDV error
Nep_TA_error_step2 <- Nep_TA_error_step1

for (i in 2:length(Nep_TA_error_step1)) {
  
  Nep_TA_error_step2[,i] <- sqrt((Nep_TA_error_step1[,i]^2) + (gas_cal_error[,i]^2))
}

############################################################################################################
# plot raw data
############################################################################################################

Cairo(file = "Maeli2_Nep_original.png", 
      bg="white",
      type="png",
      units="in", 
      width=7*3, 
      height=6*3, 
      pointsize=12*3, 
      dpi = 72*2)

par(mfrow= c(1,1))

plot(Nep_TA$date, Nep_TA$Scatt_blue,
     main = bquote(bold("Maeli2, Nephelometer,"~beta[sca]~"raw data")),
     xlab = "Time",
     ylab = bquote(~beta[sca]~(Mm^-1)),
     ylim = c(0,2000), col="blue", pch = 20)
points(Nep_TA$date, Nep_TA$Scatt_green, col="green", pch = 20)
points(Nep_TA$date, Nep_TA$Scatt_red, col="red", pch = 20)
arrows(x0=Nep_TA$date, y0=Nep_TA$Scatt_blue - Nep_TA_error_step2$Scatt_blue, 
       x1=Nep_TA$date, y1=Nep_TA$Scatt_blue + Nep_TA_error_step2$Scatt_blue, code=3, angle=90, length=0.1, col="blue", lwd=1)
arrows(x0=Nep_TA$date, y0=Nep_TA$Scatt_green - Nep_TA_error_step2$Scatt_green, 
       x1=Nep_TA$date, y1=Nep_TA$Scatt_green + Nep_TA_error_step2$Scatt_green, code=3, angle=90, length=0.1, col="green", lwd=1)
arrows(x0=Nep_TA$date, y0=Nep_TA$Scatt_red - Nep_TA_error_step2$Scatt_red, 
       x1=Nep_TA$date, y1=Nep_TA$Scatt_red + Nep_TA_error_step2$Scatt_red, code=3, angle=90, length=0.1, col="red", lwd=1)

legend("topright", inset =.02, box.lty = 0,
       legend = c("450 nm", "550 nm", "700 nm"),
       pch = c(20,20,20),
       col = c("blue","green", "red"))

############################################################################################################
# correct data for truncation 
############################################################################################################

# work only on the measured scattering 
Nep_TA_corrected <- Nep_TA[,c(1:4)]

Nep_TA_corrected$Scatt_blue <- Nep_TA$Scatt_blue*Trunc_correction_450nm$mean
Nep_TA_corrected$Scatt_green <- Nep_TA$Scatt_green*Trunc_correction_550nm$mean
Nep_TA_corrected$Scatt_red <- Nep_TA$Scatt_red*Trunc_correction_700nm$mean

############################################################################################################
# add the uncertainty related to the truncation correction - step3
############################################################################################################

# work only on the measured scattering 
Nep_TA_error_step2 <- Nep_TA_error_step2[,1:4]

# add the uncertainty related to the truncation correction
# C = A*B, STDV(C) = STDV(A*B) = C*SQRT(dA/A^2 + dB/B^2) 
Nep_TA_error_step3 <- Nep_TA_error_step2

# 450 nm
Nep_TA_error_step3$Scatt_blue <- Nep_TA_corrected$Scatt_blue*sqrt(((Nep_TA_error_step2$Scatt_blue/Nep_TA$Scatt_blue)^2) + ((Trunc_correction_450nm$sd/Trunc_correction_450nm$mean)^2))

# 550 nm
Nep_TA_error_step3$Scatt_green <- Nep_TA_corrected$Scatt_green*sqrt(((Nep_TA_error_step2$Scatt_green/Nep_TA$Scatt_green)^2) + ((Trunc_correction_550nm$sd/Trunc_correction_550nm$mean)^2))

# 700 nm
Nep_TA_error_step3$Scatt_red <- Nep_TA_corrected$Scatt_red*sqrt(((Nep_TA_error_step2$Scatt_red/Nep_TA$Scatt_red)^2) + ((Trunc_correction_700nm$sd/Trunc_correction_700nm$mean)^2))

############################################################################################################
# plot data corrected
############################################################################################################

Cairo(file = "Maeli2_Nep_corrected_original.png", 
      bg="white",
      type="png",
      units="in", 
      width=7*3, 
      height=6*3, 
      pointsize=12*3, 
      dpi = 72*2)

par(mfrow= c(1,1))

plot(Nep_TA_corrected$date, Nep_TA_corrected$Scatt_blue,
     main = bquote(bold("Maeli2, Nephelometer,"~beta[sca]~"corrected for truncation")),
     xlab = "Time",
     ylab = bquote(~beta[sca]~(Mm^-1)),
     ylim = c(0,2000), col="blue", pch = 20)
points(Nep_TA_corrected$date, Nep_TA_corrected$Scatt_green, col="green", pch = 20)
points(Nep_TA_corrected$date, Nep_TA_corrected$Scatt_red, col="red", pch = 20)
arrows(x0=Nep_TA_corrected$date, y0=Nep_TA_corrected$Scatt_blue - Nep_TA_error_step3$Scatt_blue, 
       x1=Nep_TA_corrected$date, y1=Nep_TA_corrected$Scatt_blue + Nep_TA_error_step3$Scatt_blue, code=3, angle=90, length=0.1, col="blue", lwd=1)
arrows(x0=Nep_TA_corrected$date, y0=Nep_TA_corrected$Scatt_green - Nep_TA_error_step3$Scatt_green, 
       x1=Nep_TA_corrected$date, y1=Nep_TA_corrected$Scatt_green + Nep_TA_error_step3$Scatt_green, code=3, angle=90, length=0.1, col="green", lwd=1)
arrows(x0=Nep_TA_corrected$date, y0=Nep_TA_corrected$Scatt_red - Nep_TA_error_step3$Scatt_red, 
       x1=Nep_TA_corrected$date, y1=Nep_TA_corrected$Scatt_red + Nep_TA_error_step3$Scatt_red, code=3, angle=90, length=0.1, col="red", lwd=1)

legend("topright", inset =.02, box.lty = 0,
       legend = c("450 nm", "550 nm", "700 nm"),
       pch = c(20,20,20),
       col = c("blue","green", "red"))

############################################################################################################
# predict Bsca at the Aethalometer wavelengths 370, 470, 520, 590, 660, 880, 950 nm
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
# Therefore, some observations may be more important/ weigh higher than others.

# calculate the inverse variance = 1/STDV^2
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
# extrapolate Bsca values at the Aet WL
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
  # for each row (date-time) is applied a specific powerLaw fit
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

# find the results for all the aethalometer WL
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

# calculate the root mean square error on Bsca predicted 
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

############################################################################################################
# plot fitted data
############################################################################################################

par(mfrow= c(1,1))

plot(Bsca_fitted_df$date, Bsca_fitted_df$Bsca_370,
     main = bquote(bold("Maeli2, Nephelometer,"~beta[sca]~"predicted")),
     xlab = "Time",
     ylab = bquote(beta[sca]~(Mm^-1)),
     ylim = c(0,2000), col="violet", pch = 20)
points(Bsca_fitted_df$date, Bsca_fitted_df$Bsca_470, col="blue", pch = 20)
points(Bsca_fitted_df$date, Bsca_fitted_df$Bsca_520, col="green", pch = 20)
points(Bsca_fitted_df$date, Bsca_fitted_df$Bsca_590, col="yellow", pch = 20)
points(Bsca_fitted_df$date, Bsca_fitted_df$Bsca_660, col="hotpink", pch = 20)
points(Bsca_fitted_df$date, Bsca_fitted_df$Bsca_880, col="red", pch = 20)
points(Bsca_fitted_df$date, Bsca_fitted_df$Bsca_950, col="red4", pch = 20)

arrows(x0=Bsca_fitted_df$date, y0=Bsca_fitted_df$Bsca_370 - Nep_TA_error_step4$Bsca_370, 
       x1=Bsca_fitted_df$date, y1=Bsca_fitted_df$Bsca_370 + Nep_TA_error_step4$Bsca_370, code=3, angle=90, length=0.1, col="violet", lwd=1)
arrows(x0=Bsca_fitted_df$date, y0=Bsca_fitted_df$Bsca_470 - Nep_TA_error_step4$Bsca_470, 
       x1=Bsca_fitted_df$date, y1=Bsca_fitted_df$Bsca_470 + Nep_TA_error_step4$Bsca_470, code=3, angle=90, length=0.1, col="blue", lwd=1)
arrows(x0=Bsca_fitted_df$date, y0=Bsca_fitted_df$Bsca_520 - Nep_TA_error_step4$Bsca_520, 
       x1=Bsca_fitted_df$date, y1=Bsca_fitted_df$Bsca_520 + Nep_TA_error_step4$Bsca_520, code=3, angle=90, length=0.1, col="green", lwd=1)
arrows(x0=Bsca_fitted_df$date, y0=Bsca_fitted_df$Bsca_590 - Nep_TA_error_step4$Bsca_590, 
       x1=Bsca_fitted_df$date, y1=Bsca_fitted_df$Bsca_590 + Nep_TA_error_step4$Bsca_590, code=3, angle=90, length=0.1, col="yellow", lwd=1)
arrows(x0=Bsca_fitted_df$date, y0=Bsca_fitted_df$Bsca_660 - Nep_TA_error_step4$Bsca_660, 
       x1=Bsca_fitted_df$date, y1=Bsca_fitted_df$Bsca_660 + Nep_TA_error_step4$Bsca_660, code=3, angle=90, length=0.1, col="hotpink", lwd=1)
arrows(x0=Bsca_fitted_df$date, y0=Bsca_fitted_df$Bsca_880 - Nep_TA_error_step4$Bsca_880, 
       x1=Bsca_fitted_df$date, y1=Bsca_fitted_df$Bsca_880 + Nep_TA_error_step4$Bsca_880, code=3, angle=90, length=0.1, col="red", lwd=1)
arrows(x0=Bsca_fitted_df$date, y0=Bsca_fitted_df$Bsca_950 - Nep_TA_error_step4$Bsca_950, 
       x1=Bsca_fitted_df$date, y1=Bsca_fitted_df$Bsca_950 + Nep_TA_error_step4$Bsca_950, code=3, angle=90, length=0.1, col="red4", lwd=1)

legend("topright", inset =.02, box.lty = 0,
       legend = c("370 nm", "470 nm", "520 nm",
                  "590 nm", "660 nm", "880 nm",
                  "950 nm"),
       pch = c(20,20,20,20,20,20,20),
       col = c("violet","blue","green","yellow", "hotpink", "red", "red4"))

#############################################################################################################
#############################################################################################################
# Aet measured Mie coefficients 
#############################################################################################################
#############################################################################################################
# remove terms from R environment except for:
rm(list=setdiff(ls(), c("convert_list_df", 
                        "Bsca_fitted_df", "Nep_TA_error_step4", 
                        "Trunc_correction_450nm", "Trunc_correction_550nm", "Trunc_correction_700nm")))

# import data
Aet <- read.csv("Z:/Clarissa/Data_Optical_calculation/R_Code/Code_final/ICEdust_samples/Maeli2/Aet_corrections/original/Maeli2_Aet_Babs_2min_original.csv")

############################################################################################################
# calculate 12-min average
############################################################################################################

# function to calculate the time average at different time resoluiton, by min (here 12min)
timeAverage_12min <- function(df) {
  # create a complete sequence of date-time 
  date_seq <-  seq(as.POSIXct(df$date[1]), as.POSIXct(df$date[nrow(df)]), by = "min")
  date_seq <- as.data.frame(date_seq)
  colnames(date_seq) <- "date"
  
  # merge the complete sequence of date-time and the original dataset 
  # the missing date-time will be shown now as NA
  df <- full_join(date_seq, df, by = NULL)
  
  # to calculate 12 minute average split the data in #df = dim(df)[1]/12 of 12 rows
  
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
  
  # get the time-averaged sequence
  df_TA$date <- seq(as.POSIXct(df_split[[1]][1,1]), as.POSIXct(df_split[[length(df_split)]][1,1]), by = "12 min")
  
  return(df_TA)
}

# ensure the date-time column is as "date"
colnames(Aet)[1] <- "date"

# convert as POSIXct
Aet$date <- as.POSIXct(Aet$date,  format =  "%Y-%m-%d %H:%M:%S")

# work on the data starting as the model values (use time intervals in any df in calculated Bsca, Bsca list)
Aet_subset <- subset(Aet, Aet$date >= Trunc_correction_700nm$date[1])

# verify the presence of duplicates in the row sequence
Aet_subset <- Aet_subset[which(duplicated(Aet_subset$date) == F),]

# calculate 12min average
Aet_data <- Aet_subset[,1:8]
Aet_TA <- timeAverage_12min(Aet_data)

############################################################################################################
# calculate 12-min average error
############################################################################################################

# function to calculate STDV over 12 min intervals considering the original data uncertainty
timeSTDV_12min <- function(df) {
  # create a complete sequence of date-time 
  date_seq <-  seq(as.POSIXct(df$date[1]), as.POSIXct(df$date[nrow(df)]), by = "min")
  date_seq <- as.data.frame(date_seq)
  colnames(date_seq) <- "date"
  
  # merge the complete sequence of date-time and the original dataset 
  # the missing date-time will be shown now as NA
  df <- full_join(date_seq, df, by = NULL)
  
  # to calculate 12 minute average split the data in #df = dim(df)[1]/12 of 12 rows
  
  # first, calculate the nrow to delete to average constant time intervals
  if ((dim(df)[1]- round(dim(df)[1]/12)*12) < 0) {
    df <-df[1:(dim(df)[1]- (dim(df)[1]- (round(dim(df)[1]/12)-1)*12)),]
  } else {df <-df[1:(dim(df)[1]- (dim(df)[1]- round(dim(df)[1]/12)*12)),]} 
  
  # split the df in constant time intervals
  df_split <- split(df,rep(1:round(dim(df)[1]/12),each=12))
  
  # calculate the average
  df_error <- data.frame(matrix(ncol = ncol(df), nrow = length(df_split)))
  colnames(df_error) <- colnames(df)
  for (i in 1:length(df_split)) {
    for (j in 2:length(df_split[[i]])) {
      df_error[i,j] <- sd(df_split[[i]][,j], na.rm = T)
    }
  }
  
  # get the time-averaged sequence
  df_error$date <- seq(as.POSIXct(df_split[[1]][1,1]), as.POSIXct(df_split[[length(df_split)]][1,1]), by = "12 min")
  
  return(df_error)
}

# calculate error over 12 min intervals
Aet_TA_error <- timeSTDV_12min(Aet_data)

############################################################################################################
# calculate the systematic error derived from the Aet corrections for each 12-min intervals
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

# rm component of Aet error from the environment
rm("Aet_TA_Sys_error")
rm("Aet_TA_Sys_error_rsd")
rm("Aet_TA_error")
rm("Aet_error_rsd")

#############################################################################################################
# plot Babs 
#############################################################################################################

par(mfrow= c(1,1))

plot(Aet_TA$date, Aet_TA$Babs_370,
     xlab = "Time",
     ylab = bquote(beta[abs]~(Mm^-1)),
     ylim = c(0, 400),
     main = bquote(bold("Maeli2, absorption coefficient -"~beta[abs](lambda))),
     col="violet", pch = 20)

points(Aet_TA$date, Aet_TA$Babs_470, col="blue", pch = 20)
points(Aet_TA$date, Aet_TA$Babs_520, col="green", pch = 20)
points(Aet_TA$date, Aet_TA$Babs_590, col="yellow", pch = 20)
points(Aet_TA$date, Aet_TA$Babs_660, col="hotpink", pch = 20)
points(Aet_TA$date, Aet_TA$Babs_880, col="red", pch = 20)
points(Aet_TA$date, Aet_TA$Babs_950, col="red4", pch = 20)

arrows(x0=Aet_TA$date, y0=Aet_TA$Babs_370 - Aet_TA_error_total$Babs_370, 
       x1=Aet_TA$date, y1=Aet_TA$Babs_370 + Aet_TA_error_total$Babs_370, code=3, angle=90, length=0.1, col="violet", lwd=1)
arrows(x0=Aet_TA$date, y0=Aet_TA$Babs_470 - Aet_TA_error_total$Babs_470, 
       x1=Aet_TA$date, y1=Aet_TA$Babs_470 + Aet_TA_error_total$Babs_470, code=3, angle=90, length=0.1, col="blue", lwd=1)
arrows(x0=Aet_TA$date, y0=Aet_TA$Babs_520 - Aet_TA_error_total$Babs_520, 
       x1=Aet_TA$date, y1=Aet_TA$Babs_520 + Aet_TA_error_total$Babs_520, code=3, angle=90, length=0.1, col="green", lwd=1)
arrows(x0=Aet_TA$date, y0=Aet_TA$Babs_590 - Aet_TA_error_total$Babs_590, 
       x1=Aet_TA$date, y1=Aet_TA$Babs_590 + Aet_TA_error_total$Babs_590, code=3, angle=90, length=0.1, col="yellow", lwd=1)
arrows(x0=Aet_TA$date, y0=Aet_TA$Babs_660 - Aet_TA_error_total$Babs_660, 
       x1=Aet_TA$date, y1=Aet_TA$Babs_660 + Aet_TA_error_total$Babs_660, code=3, angle=90, length=0.1, col="hotpink", lwd=1)
arrows(x0=Aet_TA$date, y0=Aet_TA$Babs_880 - Aet_TA_error_total$Babs_880, 
       x1=Aet_TA$date, y1=Aet_TA$Babs_880 + Aet_TA_error_total$Babs_880, code=3, angle=90, length=0.1, col="red", lwd=1)
arrows(x0=Aet_TA$date, y0=Aet_TA$Babs_950 - Aet_TA_error_total$Babs_950, 
       x1=Aet_TA$date, y1=Aet_TA$Babs_950 + Aet_TA_error_total$Babs_950, code=3, angle=90, length=0.1, col="red4", lwd=1)


legend("topright", inset =.02, box.lty = 0,
       legend = c("370 nm", "470 nm", "520 nm",
                  "590 nm", "660 nm", "880 nm", "950 nm"),
       pch = c(20,20,20,20,20,20,20),
       col = c("violet","blue","green","yellow","hotpink", "red", "red4"))

dev.off()

#############################################################################################################
# calculate SSA
#############################################################################################################

# for each WL calculate SSA
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

# calculate Bext error (Bsca (predicted) error + Babs* error)
Bext_error <- Bext

for (i in 2:length(Aet_TA)) {
  Bext_error[,i] <- sqrt(Aet_TA_error_total[,i]^2 + Nep_TA_error_step4[,i]^2)
}

# calculate the error related to SSA (Bext error * Bsca (predicted) error)
SSA_error <- SSA
for (i in 2:length(SSA)) {
  SSA_error[,i] <- SSA[,i]*sqrt((Bext_error[,i]/Bext[,i])^2 + (Nep_TA_error_step4[,i]/Bsca_fitted_df[,i])^2)
}

#############################################################################################################
# plot SSA
#############################################################################################################

par(mfrow= c(1,1))

plot(SSA$date, SSA$SSA_370,
     ylab = "SSA",
     xlab = "Time",
     ylim = c(0,1.5),
     main = bquote(bold("Maeli2, single scattering albedo - SSA"*(lambda))),
     col="violet", pch = 20)
points(SSA$date, SSA$SSA_470, col="blue", pch = 20)
points(SSA$date, SSA$SSA_520, col="green", pch = 20)
points(SSA$date, SSA$SSA_590, col="yellow", pch = 20)
points(SSA$date, SSA$SSA_660, col="hotpink", pch = 20)
points(SSA$date, SSA$SSA_880, col="red", pch = 20)
points(SSA$date, SSA$SSA_950, col="red4", pch = 20)

arrows(x0=SSA$date, y0=SSA$SSA_370 - SSA_error$SSA_370, 
       x1=SSA$date, y1=SSA$SSA_370 + SSA_error$SSA_370, code=3, angle=90, length=0.1, col="violet", lwd=1)
arrows(x0=SSA$date, y0=SSA$SSA_470 - SSA_error$SSA_470, 
       x1=SSA$date, y1=SSA$SSA_470 + SSA_error$SSA_470, code=3, angle=90, length=0.1, col="blue", lwd=1)
arrows(x0=SSA$date, y0=SSA$SSA_520 - SSA_error$SSA_520, 
       x1=SSA$date, y1=SSA$SSA_520 + SSA_error$SSA_520, code=3, angle=90, length=0.1, col="green", lwd=1)
arrows(x0=SSA$date, y0=SSA$SSA_590 - SSA_error$SSA_590, 
       x1=SSA$date, y1=SSA$SSA_590 + SSA_error$SSA_590, code=3, angle=90, length=0.1, col="yellow", lwd=1)
arrows(x0=SSA$date, y0=SSA$SSA_660 - SSA_error$SSA_660, 
       x1=SSA$date, y1=SSA$SSA_660 + SSA_error$SSA_660, code=3, angle=90, length=0.1, col="hotpink", lwd=1)
arrows(x0=SSA$date, y0=SSA$SSA_880 - SSA_error$SSA_880, 
       x1=SSA$date, y1=SSA$SSA_880 + SSA_error$SSA_880, code=3, angle=90, length=0.1, col="red", lwd=1)
arrows(x0=SSA$date, y0=SSA$SSA_950 - SSA_error$SSA_950, 
       x1=SSA$date, y1=SSA$SSA_950 + SSA_error$SSA_950, code=3, angle=90, length=0.1, col="red4", lwd=1)

#############################################################################################################
#############################################################################################################
# save the output
#############################################################################################################
#############################################################################################################

# calculate rsd and tidy up dfs before saving the output as .csv

# tidy Aet_TA
Aet_TA_rsd <- Aet_TA_error_total

for (i in 2:length(Aet_TA_error_total)) {
  Aet_TA_rsd[,i] <- Aet_TA_error_total[,i]/Aet_TA[,i]*100
  colnames(Aet_TA_rsd)[i] <- sprintf("rsd_%d", Aet_WL[i-1])
  colnames(Aet_TA_error_total)[i] <- sprintf("sd_%d", Aet_WL[i-1])
}

Babs <- cbind(Aet_TA, Aet_TA_error_total[,-1], Aet_TA_rsd[,-1])

# tidy Nep data
Nep_TA_rsd_step4 <- Bsca_fitted_df

for (i in 2:length(Bsca_fitted_df)) {
  Nep_TA_rsd_step4[,i] <- Nep_TA_error_step4[,i]/Bsca_fitted_df[,i]*100
  colnames(Nep_TA_rsd_step4)[i] <- sprintf("rsd_%d", Aet_WL[i-1])
  colnames(Nep_TA_error_step4)[i] <- sprintf("sd_%d", Aet_WL[i-1])
}

Bsca <- cbind(Bsca_fitted_df, Nep_TA_error_step4[,-1], Nep_TA_rsd_step4[,-1])

# tidy SSA
SSA_rsd <- SSA_error

for (i in 2:length(SSA_error)) {
  SSA_rsd[,i] <- SSA_error[,i]/SSA[,i]*100
  colnames(SSA_rsd)[i] <- sprintf("rsd_%d", Aet_WL[i-1])
  colnames(SSA_error)[i] <- sprintf("sd_%d", Aet_WL[i-1])
}

SSA <- cbind(SSA, SSA_error[,-1], SSA_rsd[,-1])

# tidy Bext
Bext_rsd <- Bext_error

for (i in 2:length(Bext_error)) {
  Bext_rsd[,i] <- Bext_error[,i]/Bext[,i]*100
  colnames(Bext_rsd)[i] <- sprintf("rsd_%d", Aet_WL[i-1])
  colnames(Bext_error)[i] <- sprintf("sd_%d", Aet_WL[i-1])
  colnames(Bext)[i] <- sprintf("Bext_%d", Aet_WL[i-1])
}

Bext <- cbind(Bext, Bext_error[,-1], Bext_rsd[,-1])

# save the output
write.csv(Babs, "Z:/Clarissa/Data_Optical_calculation/R_Code/Code_final/ICEdust_samples/Maeli2/Nep_Aet_Mie/original/data/Maeli2_Babs_12min_original.csv", row.names = F)
write.csv(Bsca, "Z:/Clarissa/Data_Optical_calculation/R_Code/Code_final/ICEdust_samples/Maeli2/Nep_Aet_Mie/original/data/Maeli2_Bsca_fitted_12min_original.csv", row.names = F)
write.csv(SSA, "Z:/Clarissa/Data_Optical_calculation/R_Code/Code_final/ICEdust_samples/Maeli2/Nep_Aet_Mie/original/data/Maeli2_SSA_12min_original.csv", row.names = F)
write.csv(Bext, "Z:/Clarissa/Data_Optical_calculation/R_Code/Code_final/ICEdust_samples/Maeli2/Nep_Aet_Mie/original/data/Maeli2_Bext_12min_original.csv", row.names = F)
