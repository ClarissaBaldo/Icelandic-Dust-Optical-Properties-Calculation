#############################################################################################################
# set working directory and upload packages
#############################################################################################################

setwd("Z:/Clarissa/Data_Optical_calculation/R_Code/Code_final/ICEdust_samples/Maeli2/Aet_corrections/original")
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

timeAverage_2min <- function(df) {
  
  # create a complete sequence of date-time 
  date_seq <-  seq(as.POSIXct(df$date[1]), as.POSIXct(df$date[nrow(df)]), by = "sec")
  date_seq <- as.data.frame(date_seq)
  colnames(date_seq) <- "date"
  
  # merge the complete sequence of date-time and the original dataset 
  # the missing date-time will be shown now as NA
  df <- full_join(date_seq, df, by = NULL)
  
  # to calculate 2 minute average split the data in #df = dim(df)[1]/120 of 120 rows
  
  # first, calculate the nrow to delete to average constant time intervals
  if ((dim(df)[1]- round(dim(df)[1]/120)*120) < 0) {
    df <-df[1:(dim(df)[1]- (dim(df)[1]- (round(dim(df)[1]/120)-1)*120)),]
  } else {df <-df[1:(dim(df)[1]- (dim(df)[1]- round(dim(df)[1]/120)*120)),]} 
  
  # split the df in constant time intervals
  df_split <- split(df,rep(1:round(dim(df)[1]/120),each=120))
  
  # calculate the average
  df_TA <- data.frame(matrix(ncol = ncol(df), nrow = length(df_split)))
  colnames(df_TA) <- colnames(df)
  for (i in 1:length(df_split)) {
    for (j in 2:length(df_split[[i]])) {
      df_TA[i,j] <- mean(df_split[[i]][,j], na.rm = T)
    }
  }
  
  # get the time-averaged sequence
  df_TA$date <- seq(as.POSIXct(df_split[[1]][1,1]), as.POSIXct(df_split[[length(df_split)]][1,1]), by = "2 min")
  
  return(df_TA)
}

timeAverage_2min_error <- function(df) {
  
  # create a complete sequence of date-time 
  date_seq <-  seq(as.POSIXct(df$date[1]), as.POSIXct(df$date[nrow(df)]), by = "sec")
  date_seq <- as.data.frame(date_seq)
  colnames(date_seq) <- "date"
  
  # merge the complete sequence of date-time and the original dataset 
  # the missing date-time will be shown now as NA
  df <- full_join(date_seq, df, by = NULL)
  
  # to calculate 2 minute average split the data in #df = dim(df)[1]/120 of 120 rows
  
  # first, calculate the nrow to delete to average constant time intervals
  if ((dim(df)[1]- round(dim(df)[1]/120)*120) < 0) {
    df <-df[1:(dim(df)[1]- (dim(df)[1]- (round(dim(df)[1]/120)-1)*120)),]
  } else {df <-df[1:(dim(df)[1]- (dim(df)[1]- round(dim(df)[1]/120)*120)),]} 
  
  # split the df in constant time intervals
  df_split <- split(df,rep(1:round(dim(df)[1]/120),each=120))
  
  # calculate the average
  df_sd <- data.frame(matrix(ncol = ncol(df), nrow = length(df_split)))
  colnames(df_sd) <- colnames(df)
  for (i in 1:length(df_split)) {
    for (j in 2:length(df_split[[i]])) {
      df_sd[i,j] <- sd(df_split[[i]][,j], na.rm = T)
    }
  }
  
  # get the time-averaged sequence
  df_sd$date <- seq(as.POSIXct(df_split[[1]][1,1]), as.POSIXct(df_split[[length(df_split)]][1,1]), by = "2 min")
  
  return(df_sd)
}

#############################################################################################################
#############################################################################################################
# Aethalometer 
#############################################################################################################
#############################################################################################################

# import data
Aet <- read.csv("Z:/Clarissa/Data_Optical_calculation/R_Code/Code_final/ICEdust_samples/Maeli2/datasets/tidy_datasets/Maeli2_Aet_data_tidy.csv", header = T)

Aet$date <- as.POSIXct(Aet$date,  format =  "%Y-%m-%d %H:%M:%S")
# Round to minutes
Aet$date <- strptime(Aet$date,  format = "%Y-%m-%d %H:%M")
Aet$date <- as.POSIXct(Aet$date,  format =  "%Y-%m-%d %H:%M:%S")

# work on the Aet starting from the samplig interval
Aet <- subset(Aet, Aet$date >= "2019-01-21 11:07:00")

# ensure that the time serie of the Aet data is complete and every 2 min
df = Aet
date_seq <-  seq(as.POSIXct(df$date[1]), as.POSIXct(df$date[nrow(df)]), by = "2 min")
date_seq <- as.data.frame(date_seq)
colnames(date_seq) <- "date"

Aet <- left_join(date_seq, Aet)

#############################################################################################################
# plot attenuation original data
#############################################################################################################

par(mfrow= c(1,1))

plot(Aet$date, Aet$ATTN.370,
     main = bquote(bold("Maeli2, attenuation - ATTN"*(lambda))),
     xlab = "Time",
     ylab = "ATTN (%)", pch = 20)
points(Aet$date, Aet$ATTN.470, col = "blue", pch = 20)
points(Aet$date, Aet$ATTN.520, col = "green", pch = 20)
points(Aet$date, Aet$ATTN.590, col = "red", pch = 20)
points(Aet$date, Aet$ATTN.660, col = "purple", pch = 20)
points(Aet$date, Aet$ATTN.880, col = "grey", pch = 20)
points(Aet$date, Aet$ATTN.950, col = "yellow", pch = 20)

legend("topright", inset =.02, box.lty = 0,
       legend = c("370 nm", "470 nm", "520 nm", "590 nm", "660 nm", "880 nm","950 nm"),
       pch = c(20,20,20,20,20,20,20),
       col = c("black","blue","green", "red", "purple", "grey", "yellow"))



#############################################################################################################
# calculate the attenuation coefficient: BATTN = [dATTN/dt]*A/V
#############################################################################################################

# from the original data consider only the column with attenuation data 
Aet_ATTN <- Aet[,c(1,9:15)]

# transform the ATTN% into a ratio
for (i in 2:length(Aet_ATTN)) {
  Aet_ATTN[,i] <-  Aet_ATTN[,i]/100
}

# when ATTN is > 100%, the filter in the instrument is replaced with a new one
# so, when correcting the data, different filters should be treated separately

#############################################################################################################
# calculate dATTN/dt

# verify the dt interval is constant = 2 min
dt <- c()
for (i in 1:nrow(Aet_ATTN)) {
  dt[i] <- (difftime(Aet_ATTN$date[i+1], Aet_ATTN$date[i], units = c("min"))==2)
}

dt 

# now, calculate dATTN/dt
ATTN_dt <- data.frame(matrix(ncol = length(Aet_ATTN), nrow = nrow(Aet_ATTN)-1))

for (i in 2:length(Aet_ATTN)) {
  for (j in 1:(nrow(Aet_ATTN)-1)) {
    ATTN_dt[j,i] <- (Aet_ATTN[j+1,i] -  Aet_ATTN[j,i])/2     
  }
}

colnames(ATTN_dt) <- colnames(Aet_ATTN)
ATTN_dt$date <- Aet_ATTN$date[-nrow(Aet_ATTN)]

#############################################################################################################
# calculate BATTN = [dATTN/dt]*A/V

# A = 0.5 cm^2 = 0.00005 m^2
# V = 0.002 m^3/min 
# the final BATTN should be in m^-1
# (remember to convert from m^-1 into Mm^-1 multiply by 10^6)

B_ATTN <- ATTN_dt

for (i in 2:length(ATTN_dt)) {
        B_ATTN[i] <- (ATTN_dt[i]*0.00005/0.002)
}

B_ATTN <- data.frame(B_ATTN)
colnames(B_ATTN) <- c("date",
                      "BATTN_370", 
                      "BATTN_470",
                      "BATTN_520",
                      "BATTN_590",
                      "BATTN_660",
                      "BATTN_880",
                      "BATTN_950")

#############################################################################################################
# calculate the measurement average ATTN
#############################################################################################################

# build a function to calculate measurement average for each time interval
calculate_measurement_average <- function(df) {
  # create an empty df
  # calculating the averge for each time interval, the nrow of the output will be equal to nrow(df) -1
  output_mean <- data.frame(matrix(ncol = length(df), nrow = nrow(df)-1))
  
  for (i in 2:length(df)) {
    for (j in 1:(nrow(df[-1,]))) {
      output_mean[j,i] <- mean(c(df[j,i], df[j+1,i]), na.rm = T)
    }
  }
  
  colnames(output_mean) <- colnames(df)
  output_mean$date <- df$date[-nrow(df)]
  return(output_mean) 
}

# calculate the measurement average for ATTN
ATTN_mean <- calculate_measurement_average(Aet_ATTN)

#############################################################################################################
# calculate the error on ATTN

# build a function to calculate the error related to ATTN measurement average for each time interval
calculate_measurement_error <- function(df) {
  # create an empty df
  # calculating the averge for each time interval, the nrow of the output will be equal to nrow(df) -1
  output_sd <- data.frame(matrix(ncol = length(df), nrow = nrow(df)-1))
  
  for (i in 2:length(df)) {
    for (j in 1:(nrow(df[-1,]))) {
      output_sd[j,i] <- sd(c(df[j,i], df[j+1,i]), na.rm = T)
    }
  }
  
  colnames(output_sd) <- colnames(df)
  output_sd$date <- df$date[-nrow(df)]
  return(output_sd) 
}

# apply the function to calculate the error on ATTN
ATTN_error <- calculate_measurement_error(Aet_ATTN)

#############################################################################################################
# bind ATTN and BATTN dataframes

# ensure the date column from the B_ATTN dataframe is removed
Aet_TA <- data.frame(cbind(ATTN_mean, B_ATTN[,-1]))
Aet_TA <- as.data.frame(Aet_TA)

# give new colnames
colnames(Aet_TA) <- c("date", 
                         "ATTN_370", 
                         "ATTN_470",
                         "ATTN_520",
                         "ATTN_590",
                         "ATTN_660",
                         "ATTN_880",
                         "ATTN_950",
                         "BATTN_370", 
                         "BATTN_470",
                         "BATTN_520",
                         "BATTN_590",
                         "BATTN_660",
                         "BATTN_880",
                         "BATTN_950")

# round date to seconds 
Aet_TA$date <- strptime(Aet_TA$date,  format = "%Y-%m-%d %H:%M")
Aet_TA$date <- as.POSIXct(Aet_TA$date,  format =  "%Y-%m-%d %H:%M:%S")

# reset the index
row.names(Aet_TA) <- c(1:nrow(Aet_TA))

# select the values starting from the bigginning of filter sampling
Aet_TA <- subset(Aet_TA, Aet_TA$date <= "2019-01-21  14:00:00")

#############################################################################################################
# plot BATTN at different wavelengths
#############################################################################################################

par(mfrow= c(1,1))

plot(Aet_TA$date, Aet_TA$BATTN_370,
     main = bquote(bold("Maeli2, attenuation coefficient -"~beta[ATTN](lambda))),
     xlab = "Time",
     ylab = bquote(beta[ATTN]~(m^-1)),
     ylim = c(0,0.001), pch = 20)
points(Aet_TA$date, Aet_TA$BATTN_470, col = "blue", pch = 20)
points(Aet_TA$date, Aet_TA$BATTN_520, col = "green", pch = 20)
points(Aet_TA$date, Aet_TA$BATTN_590, col = "red", pch = 20)
points(Aet_TA$date, Aet_TA$BATTN_660, col = "purple", pch = 20)
points(Aet_TA$date, Aet_TA$BATTN_880, col = "grey", pch = 20)
points(Aet_TA$date, Aet_TA$BATTN_950, col = "yellow", pch = 20)
legend("topright", inset =.02, box.lty = 0,
       legend = c("370 nm", "470 nm", "520 nm", "590 nm", "660 nm", "880 nm","950 nm"),
       pch = c(20,20,20,20,20,20,20),
       col = c("black","blue","green", "red", "purple", "grey", "yellow"))



#############################################################################################################
#############################################################################################################
# Nephelometer
#############################################################################################################
#############################################################################################################

# remove terms from R environment except for:
rm(list=setdiff(ls(), c("Aet_TA", "ATTN_error",
                        "timeAverage_2min", "timeAverage_2min_error", "convert_list_df", "calculate_measurement_average")))

# import data
Nep <- read.csv("Z:/Clarissa/Data_Optical_calculation/R_Code/Code_final/ICEdust_samples/Maeli2/datasets/tidy_datasets/Maeli2_Nep_data_tidy.csv")

# ensure the date-time column is as "date"
colnames(Nep)[1] <- "date"

# set the date as.POSIXct
Nep$date <- as.POSIXct(Nep$date, format ="%Y-%m-%d %H:%M:%S")

# work on the data starting as for the Aet data
Nep_subset <- subset(Nep, Nep$date >= Aet_TA$date[1])

# get only date and the columns with the relevant measurements (scattering)
Nep_subset  <- Nep_subset[,c(1:4)]

# verify the presence of duplicates in the row sequence
Nep_subset <- Nep_subset[which(duplicated(Nep_subset$date) == F),]

# get 2-min average
Nep_TA <- timeAverage_2min(Nep_subset)

# get 2-min average error - step1
Nep_TA_error_step1 <- timeAverage_2min_error(Nep_subset)

# calculate the measurement average within intervals for Nep to match Aet measurements
Nep_TA <- calculate_measurement_average(Nep_TA)

# calculate the error on the measurement average - step2
Nep_TA_error_step2 <- Nep_TA

for (i in 2:length(Nep_TA_error_step1)) {
  for (j in 1:nrow(Nep_TA_error_step1)-1) {
    Nep_TA_error_step2[j,i] <- 1/2*(Nep_TA_error_step1[j,i] + Nep_TA_error_step1[j+1,i]) # just average the two SD similar to a linear interpolation
  }
}

# get the same time interval as Aet
Nep_TA <- subset(Nep_TA, Nep_TA$date <= Aet_TA$date[nrow(Aet_TA)])

Nep_TA_error_step2 <- subset(Nep_TA_error_step2, Nep_TA_error_step2$date <= Aet_TA$date[nrow(Aet_TA)])

# convert the data as m-1
# (remember to convert from m^-1 into Mm^-1 multiply by 10^6)
for (i in 2:length(Nep_TA)) {
  Nep_TA[,i] <- Nep_TA[,i]/10^6
  
}

for (i in 2:length(Nep_TA_error_step2)) {
  Nep_TA_error_step2[,i] <- Nep_TA_error_step2[,i]/10^6
  
}

# add to STDV further 5% uncertainty due to gas calibration - step 3

# calculate gas calibration STDV
gas_cal_error <- Nep_TA

for (i in 2:length(Nep_TA)) {
  
  gas_cal_error[,i] <- (Nep_TA[,i]/100)*5
}

# add the gas calib error to the 12 min STDV error
Nep_TA_error_step3 <- Nep_TA_error_step2

for (i in 2:length(Nep_TA_error_step2)) {
  Nep_TA_error_step3[,i] <- sqrt((Nep_TA_error_step2[,i]^2) + (gas_cal_error[,i]^2))
}

############################################################################################################
# plot data
############################################################################################################

plot(Nep_TA$date, Nep_TA$Scatt_blue,
     main = bquote(bold("Maeli2, Nephelometer, scattering coefficient -"~beta[sca](lambda))),
     xlab = "Time",
     ylab = bquote(beta[sca]~(m^-1)),
     ylim = c(0,0.002), 
     col="blue", pch = 20)
points(Nep_TA$date, Nep_TA$Scatt_green, col="green", pch = 20)
points(Nep_TA$date, Nep_TA$Scatt_red, col="red", pch = 20)

arrows(x0=Nep_TA$date, y0=Nep_TA$Scatt_blue - Nep_TA_error_step3$Scatt_blue, 
       x1=Nep_TA$date, y1=Nep_TA$Scatt_blue + Nep_TA_error_step3$Scatt_blue, code=3, angle=90, length=0.1, col="blue", lwd=1)
arrows(x0=Nep_TA$date, y0=Nep_TA$Scatt_green - Nep_TA_error_step3$Scatt_green, 
       x1=Nep_TA$date, y1=Nep_TA$Scatt_green + Nep_TA_error_step3$Scatt_green, code=3, angle=90, length=0.1, col="green", lwd=1)
arrows(x0=Nep_TA$date, y0=Nep_TA$Scatt_red - Nep_TA_error_step3$Scatt_red, 
       x1=Nep_TA$date, y1=Nep_TA$Scatt_red + Nep_TA_error_step3$Scatt_red, code=3, angle=90, length=0.1, col="red", lwd=1)

legend("topright", inset =.02, box.lty = 0,
       legend = c("450 nm", "550 nm", "700 nm"),
       pch = c(20,20,20),
       col = c("blue","green", "red"))



############################################################################################################
# correct data for truncation 
############################################################################################################

# import correction factor for each wavelength - 12 min average
Trunc_correction_700nm <- read.csv("Z:/Clarissa/Data_Optical_calculation/R_Code/Code_final/ICEdust_samples/Maeli2/datasets/Python_output/original/Truncation_corrections/Maeli2_truncation_stat_700nm_original.csv")
Trunc_correction_550nm <- read.csv("Z:/Clarissa/Data_Optical_calculation/R_Code/Code_final/ICEdust_samples/Maeli2/datasets/Python_output/original/Truncation_corrections/Maeli2_truncation_stat_550nm_original.csv")
Trunc_correction_450nm <- read.csv("Z:/Clarissa/Data_Optical_calculation/R_Code/Code_final/ICEdust_samples/Maeli2/datasets/Python_output/original/Truncation_corrections/Maeli2_truncation_stat_450nm_original.csv")

Trunc_correction_700nm$date  <- as.POSIXct(Trunc_correction_700nm$date  ,  format =  "%Y-%m-%d %H:%M:%S")
Trunc_correction_550nm$date   <- as.POSIXct(Trunc_correction_550nm$date  ,  format =  "%Y-%m-%d %H:%M:%S")
Trunc_correction_450nm$date   <- as.POSIXct(Trunc_correction_450nm$date  ,  format =  "%Y-%m-%d %H:%M:%S")

# increase the time resolution of the truncation factors (every 2 min)
# interpolate the data using a linear model, the ratio decreases in a quite linear way
ratio_700nm <- predict(lm(mean ~ date , data = Trunc_correction_700nm), data.frame(date = Nep_TA$date))
ratio_550nm <- predict(lm(mean ~ date , data = Trunc_correction_550nm), data.frame(date = Nep_TA$date))
ratio_450nm <- predict(lm(mean ~ date , data = Trunc_correction_450nm), data.frame(date = Nep_TA$date))

# correct data for truncation 
Nep_TA_corrected <- Nep_TA

Nep_TA_corrected$Scatt_blue <- Nep_TA$Scatt_blue*ratio_450nm
Nep_TA_corrected$Scatt_green <- Nep_TA$Scatt_green*ratio_550nm
Nep_TA_corrected$Scatt_red <- Nep_TA$Scatt_red*ratio_700nm

Nep_TA_corrected$date  <- as.POSIXct(Nep_TA_corrected$date  ,  format =  "%Y-%m-%d %H:%M:%S")

############################################################################################################
# add the uncertainty related to the truncation correction - step4
############################################################################################################

# interpolate also truncation sd using a linear model
ratio_700nm_error <- predict(lm(sd ~ date , data = Trunc_correction_700nm), data.frame(date = Nep_TA$date))
ratio_550nm_error <- predict(lm(sd ~ date , data = Trunc_correction_550nm), data.frame(date = Nep_TA$date))
ratio_450nm_error <- predict(lm(sd ~ date , data = Trunc_correction_450nm), data.frame(date = Nep_TA$date))

# add the uncertainty related to the truncation correction
# C = A*B, STDV(C) = STDV(A*B) = C*SQRT(dA/A^2 + dB/B^2) 
Nep_TA_error_step4 <- Nep_TA_error_step3

# 450 nm
Nep_TA_error_step4$Scatt_blue <- Nep_TA_corrected$Scatt_blue*sqrt(((Nep_TA_error_step3$Scatt_blue/Nep_TA$Scatt_blue)^2) + ((ratio_450nm_error/ratio_450nm)^2))

# 550 nm
Nep_TA_error_step4$Scatt_green <- Nep_TA_corrected$Scatt_green*sqrt(((Nep_TA_error_step3$Scatt_green/Nep_TA$Scatt_green)^2) + ((ratio_550nm_error/ratio_550nm)^2))

# 700 nm
Nep_TA_error_step4$Scatt_red <- Nep_TA_corrected$Scatt_red*sqrt(((Nep_TA_error_step3$Scatt_red/Nep_TA$Scatt_red)^2) + ((ratio_700nm_error/ratio_700nm)^2))

############################################################################################################
# plot data corrected
############################################################################################################

plot(Nep_TA$date, Nep_TA_corrected$Scatt_blue,
     main = bquote(bold("Maeli2, Nephelometer,"~beta[sca]~"corrected for truncation")),
     xlab = "Time",
     ylab = bquote(~beta[sca]~(m^-1)),
     ylim = c(0,0.002), col="blue", pch = 20)
points(Nep_TA$date, Nep_TA_corrected$Scatt_green, col="green", pch = 20)
points(Nep_TA$date, Nep_TA_corrected$Scatt_red, col="red", pch = 20)
arrows(x0=Nep_TA$date, y0=Nep_TA_corrected$Scatt_blue - Nep_TA_error_step4$Scatt_blue, 
       x1=Nep_TA$date, y1=Nep_TA_corrected$Scatt_blue + Nep_TA_error_step4$Scatt_blue, code=3, angle=90, length=0.1, col="blue", lwd=1)
arrows(x0=Nep_TA$date, y0=Nep_TA_corrected$Scatt_green - Nep_TA_error_step4$Scatt_green, 
       x1=Nep_TA$date, y1=Nep_TA_corrected$Scatt_green + Nep_TA_error_step4$Scatt_green, code=3, angle=90, length=0.1, col="green", lwd=1)
arrows(x0=Nep_TA$date, y0=Nep_TA_corrected$Scatt_red - Nep_TA_error_step4$Scatt_red, 
       x1=Nep_TA$date, y1=Nep_TA_corrected$Scatt_red + Nep_TA_error_step4$Scatt_red, code=3, angle=90, length=0.1, col="red", lwd=1)

legend("topright", inset =.02, box.lty = 0,
       legend = c("450 nm", "550 nm", "700 nm"),
       pch = c(20,20,20),
       col = c("blue","green", "red"))



#############################################################################################################
#############################################################################################################
#############################################################################################################
# Estimate the correction parameters for BATTN
#############################################################################################################
#############################################################################################################
#############################################################################################################

# remove terms from R environment except for:
rm(list=setdiff(ls(), c("Aet_TA", "ATTN_error", 
                        "Nep_TA", "Nep_TA_error_step3",
                        "Nep_TA_corrected", "Nep_TA_error_step4", 
                        "convert_list_df")))

#############################################################################################################
#############################################################################################################
# Scattering effect correction - SE_correction
#############################################################################################################
#############################################################################################################

# Method from Collaud Coen et al. 2010 and Di Biagio et al. 2019

# SE_correction(WL) =A^(d-1)*c*WL^(-AE*(d-1))

# WL = wavelength
# AE = Angstrom exponent
# Bsca(WL) = A*WL^AE (powerLaw fit)

# c = 3.29*10^-4
# d = 0.564

#############################################################################################################
# Bsca(WL) = A*WL^AE
#############################################################################################################

# Nep wavelength (WL)
x <-  c(450,550,700)

# for each time get Bsca corrected for truncation at the 3 WLs
y <- list()

for (i in 1:nrow(Nep_TA_corrected)) {
  y[[i]] <- c(Nep_TA_corrected[i,2], Nep_TA_corrected[i,3], Nep_TA_corrected[i,4])
}

# Weight variable in regression, is a measure of how important an observation is to your model due to different reasons 
# (eg. may be in terms of reliability of measurement or inverse of variance estimate). 
# Therefore, some observations may be more important/ weigh higher than others.

# calculate the inverse variance = 1/STDV^2
# this wil be used as weights related to Bsca_corrected (Y)
Weights_y <- list()

for (i in 1:nrow(Nep_TA_error_step4)) {
  Weights_y[[i]] <- c((1/Nep_TA_error_step4[i,2])^2, (1/Nep_TA_error_step4[i,3])^2, (1/Nep_TA_error_step4[i,4])^2)
}

# build a dataframe x-y
xy <- list()

for (i in 1:length(y)) {
  xy[[i]] <- data.frame(cbind(x,y[[i]]))
  colnames(xy[[i]]) <- c("x", "y")
}

# define the starting a and b values for the powerLaw fitting
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

# fit the data using powerLaw fitting
model_powerLaw <- list()

for (i in 1:length(xy)) {
  if (("TRUE" %in% is.na(xy[[i]]$y))|("TRUE" %in% is.na(Weights_y[[i]]))){ # if any term is NA
    model_powerLaw[[i]] <- NA
  }else{  
    model_powerLaw[[i]] <- nls(y ~ a*x^(b), start = list(a=exp(a[[i]]), b=b[[i]]), data = xy[[i]], control=nls.control(maxiter=1000), weights = Weights_y[[i]])
  }
}

# get the a and b parameters
# a = A
# b = AE

A <- list()
AE <- list()

for (i in 1:length(model_powerLaw)) {
  if (("TRUE" %in% is.na(xy[[i]]$y))|("TRUE" %in% is.na(Weights_y[[i]]))){ # if any term is NA
    A[[i]] <- NA
    AE[[i]] <- NA
  }else{
    A[[i]] <-  summary(model_powerLaw[[i]])[["coefficients"]][1]
    AE[[i]] <-  summary(model_powerLaw[[i]])[["coefficients"]][2]
  }
}

#############################################################################################################
# get the error on A and AE parameters

# a = A
# b = AE

# save the error values
A_error <- list()
AE_error <- list()

for (i in 1:length(model_powerLaw)) {
  if (("TRUE" %in% is.na(xy[[i]]$y))|("TRUE" %in% is.na(Weights_y[[i]]))){ # if any term is NA
    A_error[[i]] <- NA
    AE_error[[i]] <- NA
  }else{
    A_error[[i]] <-  summary(model_powerLaw[[i]])[["coefficients"]][3]
    AE_error[[i]] <-  summary(model_powerLaw[[i]])[["coefficients"]][4]
  }
}

# unlist the results
AE_error <- unlist(AE_error)
A_error <- unlist(A_error)

#############################################################################################################
# check the fitted lines
#############################################################################################################

par(mfrow= c(2,2))

plot(xy[[1]]$x, xy[[1]]$y,
     ylab = bquote(beta[sca]~(m^-1)),
     xlab = bquote(lambda~(nm)), pch = 20)
lines(xy[[1]]$x, predict(model_powerLaw[[1]]))
mtext(bquote("Power-law fit of"~beta[sca](lambda)~"vs"~lambda~", at 11:07"), side = 3, cex = 0.75) 
text(550, 0.00168, sprintf("Y = %sX^%s", 
                       formatC(A[[1]], format = "e", digits = 2), 
                       formatC(AE[[1]], format = "f", digits = 2)),
                       col = "red",
                       font = 2)

plot(xy[[16]]$x, xy[[16]]$y,
     ylab = bquote(beta[sca]~(m^-1)),
     xlab = bquote(lambda~(nm)), pch = 20)
lines(xy[[16]]$x, predict(model_powerLaw[[16]]))
mtext(bquote("Power-law fit of"~beta[sca](lambda)~"vs"~lambda~", at 11:37"), side = 3, cex = 0.75)  
text(550, 0.00085, sprintf("Y = %sX^%s", 
                       formatC(A[[16]], format = "e", digits = 2), 
                       formatC(AE[[16]], format = "f", digits = 2)),
                       col = "red",
                       font = 2)


plot(xy[[31]]$x, xy[[31]]$y,
     ylab = bquote(beta[sca]~(m^-1)),
     xlab = bquote(lambda~(nm)), pch = 20)
lines(xy[[31]]$x, predict(model_powerLaw[[31]]))
mtext(bquote("Power-law fit of"~beta[sca](lambda)~"vs"~lambda~", at 12:07"), side = 3, cex = 0.75)  
text(550, 0.000510, sprintf("Y = %sX^%s", 
                       formatC(A[[31]], format = "e", digits = 2), 
                       formatC(AE[[31]], format = "f", digits = 2)),
                       col = "red",
                       font = 2)


plot(xy[[46]]$x, xy[[46]]$y,
     ylab = bquote(beta[sca]~(m^-1)),
     xlab = bquote(lambda~(nm)), pch = 20)
lines(xy[[46]]$x, predict(model_powerLaw[[46]]))
mtext(bquote("Power-law fit of"~beta[sca](lambda)~"vs"~lambda~", at 12:37"), side = 3, cex = 0.75)  
text(550, 0.000325, sprintf("Y = %sX^%s", 
                       formatC(A[[46]], format = "e", digits = 2), 
                       formatC(AE[[46]], format = "f", digits = 2)),
                       col = "red",
                       font = 2)




#############################################################################################################
# extrapolate Bsca values at the Aet WL
#############################################################################################################

# the new WL class must be as dataframe
# Nep WL were also added to estimate the uncertainty related to the predicted Bsca values
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
Bsca_fitted_df <- cbind(Nep_TA_corrected$date, Bsca_fitted_df)
Bsca_fitted_df <- data.frame(Bsca_fitted_df)
colnames(Bsca_fitted_df)[1] <- "date"

#############################################################################################################
# calculate the error on Bsca predicted - step 5
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
# weights are based on the uncertainty on Bsca measured (corrected values)

Bsca_fitted_error <- list()

for (i in 1:length(xy)) {
  if (("TRUE" %in% is.na(xy[[i]]$y))|("TRUE" %in% is.na(Weights_y[[i]]))){ # if any term is NA
    Bsca_fitted_error[[i]] <- NA
  }else{
    Bsca_fitted_error[[i]] <- lm(y ~ x, data = xy[[i]],  weights = Weights_y[[i]])
  }
}

# calculate the root mean squared error (RMSE) on Bsca predicted 
# which corresponds to the standard deviation on Bsca predicted

Bsca_fitted_RMSE <- list()

for (i in 1:length(Bsca_fitted_error)) {
  if (("TRUE" %in% is.na(xy[[i]]$y))|("TRUE" %in% is.na(Weights_y[[i]]))){ # if any term is NA
    Bsca_fitted_RMSE[[i]] <- NA
  }else{
    Bsca_fitted_RMSE[[i]] <- sqrt(mean(Bsca_fitted_error[[i]]$residuals^2))
  }
}

Bsca_fitted_RMSE_df <- unlist(Bsca_fitted_RMSE)
Bsca_fitted_RMSE_df <- data.frame("date" = Nep_TA_corrected$date, "RMSE" = Bsca_fitted_RMSE_df)

# remove the columns of Bsca at 450, 550 and 700 nm
Bsca_fitted_df <- Bsca_fitted_df[,-c(9:11)]

#############################################################################################################
# calculate the average systematic error on Bsca corrected at 450, 550 and 700 nm

# calculate Bsca corrected rsd at 450, 550 and 700 nm
Nep_TA_corrected_rsd <- Nep_TA_corrected
for (i in 2:length(Nep_TA_corrected_rsd)) {
  Nep_TA_corrected_rsd[,i] <- Nep_TA_error_step4[,i]/Nep_TA_corrected[,i]*100
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
Nep_TA_error_step5 <- Bsca_fitted_systematic_error_df
for (i in 2:length(Bsca_fitted_systematic_error_df)) {
  Nep_TA_error_step5[,i] <- sqrt(Bsca_fitted_systematic_error_df[,i]^2 + Bsca_fitted_RMSE_df[,2]^2)
}

############################################################################################################
# plot fitted data
############################################################################################################

par(mfrow= c(1,1))

plot(Bsca_fitted_df$date, Bsca_fitted_df$Bsca_370,
     main = bquote(bold("Maeli2, Nephelometer,"~beta[sca]~"predicted values")),
     xlab = "Time",
     ylab = bquote(beta[sca]~(m^-1)),
     ylim = c(0,0.002), col="violet", pch = 20)
points(Bsca_fitted_df$date, Bsca_fitted_df$Bsca_470, col="blue", pch = 20)
points(Bsca_fitted_df$date, Bsca_fitted_df$Bsca_520, col="green", pch = 20)
points(Bsca_fitted_df$date, Bsca_fitted_df$Bsca_590, col="yellow", pch = 20)
points(Bsca_fitted_df$date, Bsca_fitted_df$Bsca_660, col="hotpink", pch = 20)
points(Bsca_fitted_df$date, Bsca_fitted_df$Bsca_880, col="red", pch = 20)
points(Bsca_fitted_df$date, Bsca_fitted_df$Bsca_950, col="red4", pch = 20)

arrows(x0=Bsca_fitted_df$date, y0=Bsca_fitted_df$Bsca_370 - Nep_TA_error_step5$Bsca_370, 
       x1=Bsca_fitted_df$date, y1=Bsca_fitted_df$Bsca_370 + Nep_TA_error_step5$Bsca_370, code=3, angle=90, length=0.1, col="violet", lwd=1)
arrows(x0=Bsca_fitted_df$date, y0=Bsca_fitted_df$Bsca_470 - Nep_TA_error_step5$Bsca_470, 
       x1=Bsca_fitted_df$date, y1=Bsca_fitted_df$Bsca_470 + Nep_TA_error_step5$Bsca_470, code=3, angle=90, length=0.1, col="blue", lwd=1)
arrows(x0=Bsca_fitted_df$date, y0=Bsca_fitted_df$Bsca_520 - Nep_TA_error_step5$Bsca_520, 
       x1=Bsca_fitted_df$date, y1=Bsca_fitted_df$Bsca_520 + Nep_TA_error_step5$Bsca_520, code=3, angle=90, length=0.1, col="green", lwd=1)
arrows(x0=Bsca_fitted_df$date, y0=Bsca_fitted_df$Bsca_590 - Nep_TA_error_step5$Bsca_590, 
       x1=Bsca_fitted_df$date, y1=Bsca_fitted_df$Bsca_590 + Nep_TA_error_step5$Bsca_590, code=3, angle=90, length=0.1, col="yellow", lwd=1)
arrows(x0=Bsca_fitted_df$date, y0=Bsca_fitted_df$Bsca_660 - Nep_TA_error_step5$Bsca_660, 
       x1=Bsca_fitted_df$date, y1=Bsca_fitted_df$Bsca_660 + Nep_TA_error_step5$Bsca_660, code=3, angle=90, length=0.1, col="hotpink", lwd=1)
arrows(x0=Bsca_fitted_df$date, y0=Bsca_fitted_df$Bsca_880 - Nep_TA_error_step5$Bsca_880, 
       x1=Bsca_fitted_df$date, y1=Bsca_fitted_df$Bsca_880 + Nep_TA_error_step5$Bsca_880, code=3, angle=90, length=0.1, col="red", lwd=1)
arrows(x0=Bsca_fitted_df$date, y0=Bsca_fitted_df$Bsca_950 - Nep_TA_error_step5$Bsca_950, 
       x1=Bsca_fitted_df$date, y1=Bsca_fitted_df$Bsca_950 + Nep_TA_error_step5$Bsca_950, code=3, angle=90, length=0.1, col="red4", lwd=1)

legend("topright", inset =.02, box.lty = 0,
       legend = c("370 nm", "470 nm", "520 nm",
                  "590 nm", "660 nm", "880 nm",
                  "950 nm"),
       pch = c(20,20,20,20,20,20,20),
       col = c("violet","blue","green","yellow", "hotpink", "red", "red4"))


#############################################################################################################
# calculate SE_correction at the Aet WL  
#############################################################################################################

# remove terms from R environment except for:
rm(list=setdiff(ls(), c("Aet_TA", "ATTN_error", 
                        "Nep_TA", "Nep_TA_error_step3",
                        "Nep_TA_corrected", "Nep_TA_error_step4", 
                        "A", "AE", "A_error", "AE_error",
                        "Bsca_fitted_df", "Nep_TA_error_step5",
                        "convert_list_df")))

# SE_correction(WL) =A^(d-1)*c*WL^(-AE*(d-1))
# c = 3.29*10^-4
# d = 0.564

# convert AE from list to df
AE_df <- unlist(AE)
AE_df <- data.frame("date" = Nep_TA_corrected$date, "AE" = AE_df)

# convert A from list to df
A_df <- unlist(A)
A_df <- data.frame("date" = Nep_TA_corrected$date, "A" = A_df)

# calculate the SE_correction at all the Aet WL
Aet_WL <- c(370, 470, 520, 590, 660, 880, 950)
# AE is saved already as -AE

SE_correction_list <- list()

for (i in 2:length(Bsca_fitted_df)) {
  SE_correction_list[[i]] <- ((3.29*10^-4)*(A_df[,2]^(0.564-1)))*(Aet_WL[i-1]^(AE_df[,2]*(0.564-1)))
  SE_correction_list[[i]] <- data.frame(SE_correction_list[[i]])
  colnames(SE_correction_list[[i]]) <- sprintf("SE_corr_%d", Aet_WL[i-1])
}

# drop empty terms in the list
SE_correction_list <- SE_correction_list[lengths(SE_correction_list) != 0]

# bind the results together
SE_correction_df <- list.cbind(SE_correction_list)

# add the date column 
SE_correction_df <- cbind(Nep_TA_corrected$date,SE_correction_df)
SE_correction_df <- data.frame(SE_correction_df)
colnames(SE_correction_df)[1] <- "date"

#############################################################################################################
# calculate SE_correction error
#############################################################################################################

# SE_correction(WL) =A^(d-1)*c*WL^(-AE*(d-1))
# c = 3.29*10^-4
# d = 0.564

A_df$A_error <- A_error
AE_df$AE_error <- AE_error

# calculate dSE 
# A^(d-1) error = (d-1)*dA/A
# WL^(-AE*(d-1)) error = ln(WL)*(d-1)*dAE
# dSE = SE*sqrt((c*(d-1)*dA/A)^2 + (ln(WL)*(d-1)*dAE)^2)
SE_correction_df_error <- SE_correction_df
for (i in 2:length(SE_correction_df)) {
  SE_correction_df_error[,i] <- SE_correction_df[,i]*sqrt(((0.564-1)*(A_df$A_error/A_df$A))^2 + (log(Aet_WL[i-1])*(0.564-1)*AE_df$AE_error)^2)
}

#############################################################################################################
# plot SE_correction at the Aet WL
#############################################################################################################

par(mfrow= c(1,1))

plot(SE_correction_df$date, SE_correction_df$SE_corr_370,
     ylab = bquote(alpha),
     ylim = c(0,0.04),
     xlab = "Time",
     main = bquote(bold("Maeli2, scattering effect correction -"~alpha(lambda))),
     col="violet", pch = 20)

points(SE_correction_df$date, SE_correction_df$SE_corr_470, col="blue", pch = 20)
points(SE_correction_df$date, SE_correction_df$SE_corr_520, col="green", pch = 20)
points(SE_correction_df$date, SE_correction_df$SE_corr_590, col="yellow", pch = 20)
points(SE_correction_df$date, SE_correction_df$SE_corr_660, col="hotpink", pch = 20)
points(SE_correction_df$date, SE_correction_df$SE_corr_880, col="red", pch = 20)
points(SE_correction_df$date, SE_correction_df$SE_corr_950, col="red4", pch = 20)

arrows(x0=SE_correction_df$date, y0=SE_correction_df$SE_corr_370 - SE_correction_df_error$SE_corr_370, 
       x1=SE_correction_df$date, y1=SE_correction_df$SE_corr_370 + SE_correction_df_error$SE_corr_370, code=3, angle=90, length=0.1, col="violet", lwd=1)
arrows(x0=SE_correction_df$date, y0=SE_correction_df$SE_corr_470 - SE_correction_df_error$SE_corr_470, 
       x1=SE_correction_df$date, y1=SE_correction_df$SE_corr_470 + SE_correction_df_error$SE_corr_470, code=3, angle=90, length=0.1, col="blue", lwd=1)
arrows(x0=SE_correction_df$date, y0=SE_correction_df$SE_corr_520 - SE_correction_df_error$SE_corr_520, 
       x1=SE_correction_df$date, y1=SE_correction_df$SE_corr_520 + SE_correction_df_error$SE_corr_520, code=3, angle=90, length=0.1, col="green", lwd=1)
arrows(x0=SE_correction_df$date, y0=SE_correction_df$SE_corr_590 - SE_correction_df_error$SE_corr_590, 
       x1=SE_correction_df$date, y1=SE_correction_df$SE_corr_590 + SE_correction_df_error$SE_corr_590, code=3, angle=90, length=0.1, col="yellow", lwd=1)
arrows(x0=SE_correction_df$date, y0=SE_correction_df$SE_corr_660 - SE_correction_df_error$SE_corr_660, 
       x1=SE_correction_df$date, y1=SE_correction_df$SE_corr_660 + SE_correction_df_error$SE_corr_660, code=3, angle=90, length=0.1, col="hotpink", lwd=1)
arrows(x0=SE_correction_df$date, y0=SE_correction_df$SE_corr_880 - SE_correction_df_error$SE_corr_880, 
       x1=SE_correction_df$date, y1=SE_correction_df$SE_corr_880 + SE_correction_df_error$SE_corr_880, code=3, angle=90, length=0.1, col="red", lwd=1)
arrows(x0=SE_correction_df$date, y0=SE_correction_df$SE_corr_950 - SE_correction_df_error$SE_corr_950, 
       x1=SE_correction_df$date, y1=SE_correction_df$SE_corr_950 + SE_correction_df_error$SE_corr_950, code=3, angle=90, length=0.1, col="red4", lwd=1)


legend("topleft", inset =.02, box.lty = 0,
       legend = c("370 nm", "470 nm", "520 nm",
                  "590 nm", "660 nm", "880 nm", "950 nm"),
       pch = c(20,20,20,20,20,20,20),
       col = c("violet","blue","green","yellow","hotpink", "red", "red4"))

#############################################################################################################
#############################################################################################################
# Loading effect correction - LE_correction
#############################################################################################################
#############################################################################################################

# Method from C2010 and Di Biagio et al. 2019

# LE_correction = ((1/f(WL)) - 1) * (ATT(WL)/50%) + 1

# both methods require the term f(WL) which is function of the single scattering albedo (SSA):

# f(WL)  = a(1-SSA(WL))+1

# where a = 0.74 from C2010

#############################################################################################################
# calculate SSA at all Aet WL
#############################################################################################################

# remove terms from R environment except for:
rm(list=setdiff(ls(), c("Aet_TA", "ATTN_error", 
                        "Nep_TA", "Nep_TA_error_step3",
                        "Nep_TA_corrected", "Nep_TA_error_step4", 
                        "A_df", "AE_df",
                        "Bsca_fitted_df", "Nep_TA_error_step5",
                        "convert_list_df",
                        "SE_correction_df", "SE_correction_df_error",
                        "Aet_WL")))


# SSA(WL) <- Bsca(WL)/Bext(WL) or Bsca(WL)/(Bsca(WL)+Babs(WL))

# Bsca is from the Nephelometer
# Bext is from te CAPS_Red (it did not work)
# Babs is from Aet (we can use BATTN - SE_corr, the scattering correction)

# first apply the SE correction to the Aet data
# Babs(minus SE_corr) = BATTN(WL) - SE_corr(WL)*Bsca(WL)

# only select relevant columns from Nep and Aet data 
BATTN <- Aet_TA[c(1, 9:15)]

# calculate Babs(minus SE_corr) = BATTN(WL) - SE_corr(WL)*Bsca(WL) at the Aet WL
Babs_minus_SE_list <- list()

for (i in 2:length(Bsca_fitted_df)) {
  Babs_minus_SE_list[[i]] <- BATTN[,i] - (SE_correction_df[,i]*Bsca_fitted_df[,i])
  Babs_minus_SE_list[[i]] <- data.frame(Babs_minus_SE_list[[i]])
  colnames(Babs_minus_SE_list[[i]]) <- sprintf("Babs_minus_SE_%d", Aet_WL[i-1])
}

# drop empty terms in the list
Babs_minus_SE_list <- Babs_minus_SE_list[lengths(Babs_minus_SE_list) != 0]

# bind the results together
Babs_minus_SE_df <- list.cbind(Babs_minus_SE_list)

# add the date column 
Babs_minus_SE_df <- cbind(Nep_TA_corrected$date, Babs_minus_SE_df)
Babs_minus_SE_df <- data.frame(Babs_minus_SE_df)
colnames(Babs_minus_SE_df)[1] <- "date"

# calculate the error related to Babs_minus (Se correction error * Bsca (predicted) error)
Babs_minus_SE_df_error <- Babs_minus_SE_df 
for (i in 2:length(Babs_minus_SE_df)) {
  Babs_minus_SE_df_error[,i] <- Babs_minus_SE_df[,i]*sqrt((SE_correction_df_error[,i]/SE_correction_df[,i])^2 + (Nep_TA_error_step5[,i]/Bsca_fitted_df[,i])^2)
}

#############################################################################################################
# plot Babs_minus_SE at all WL
#############################################################################################################

par(mfrow= c(1,1))

plot(Babs_minus_SE_df$date, Babs_minus_SE_df$Babs_minus_SE_370,
     ylab = bquote(beta[abs]*"*"~(m^-1)),
     xlab = "Time",
     ylim = c(0,0.0012),
     main = bquote(bold("Maeli2,"~beta[abs](lambda)*"*")),
     col="violet", pch = 20)
points(Babs_minus_SE_df$date, Babs_minus_SE_df$Babs_minus_SE_470, col="blue", pch = 20)
points(Babs_minus_SE_df$date, Babs_minus_SE_df$Babs_minus_SE_520, col="green", pch = 20)
points(Babs_minus_SE_df$date, Babs_minus_SE_df$Babs_minus_SE_590, col="yellow", pch = 20)
points(Babs_minus_SE_df$date, Babs_minus_SE_df$Babs_minus_SE_660, col="hotpink", pch = 20)
points(Babs_minus_SE_df$date, Babs_minus_SE_df$Babs_minus_SE_880, col="red", pch = 20)
points(Babs_minus_SE_df$date, Babs_minus_SE_df$Babs_minus_SE_950, col="red4", pch = 20)

arrows(x0=Babs_minus_SE_df$date, y0=Babs_minus_SE_df$Babs_minus_SE_370 - Babs_minus_SE_df_error$Babs_minus_SE_370, 
       x1=Babs_minus_SE_df$date, y1=Babs_minus_SE_df$Babs_minus_SE_370 + Babs_minus_SE_df_error$Babs_minus_SE_370, code=3, angle=90, length=0.1, col="violet", lwd=1)
arrows(x0=Babs_minus_SE_df$date, y0=Babs_minus_SE_df$Babs_minus_SE_470 - Babs_minus_SE_df_error$Babs_minus_SE_470, 
       x1=Babs_minus_SE_df$date, y1=Babs_minus_SE_df$Babs_minus_SE_470 + Babs_minus_SE_df_error$Babs_minus_SE_470, code=3, angle=90, length=0.1, col="blue", lwd=1)
arrows(x0=Babs_minus_SE_df$date, y0=Babs_minus_SE_df$Babs_minus_SE_520 - Babs_minus_SE_df_error$Babs_minus_SE_520, 
       x1=Babs_minus_SE_df$date, y1=Babs_minus_SE_df$Babs_minus_SE_520 + Babs_minus_SE_df_error$Babs_minus_SE_520, code=3, angle=90, length=0.1, col="green", lwd=1)
arrows(x0=Babs_minus_SE_df$date, y0=Babs_minus_SE_df$Babs_minus_SE_590 - Babs_minus_SE_df_error$Babs_minus_SE_590, 
       x1=Babs_minus_SE_df$date, y1=Babs_minus_SE_df$Babs_minus_SE_590 + Babs_minus_SE_df_error$Babs_minus_SE_590, code=3, angle=90, length=0.1, col="yellow", lwd=1)
arrows(x0=Babs_minus_SE_df$date, y0=Babs_minus_SE_df$Babs_minus_SE_660 - Babs_minus_SE_df_error$Babs_minus_SE_660, 
       x1=Babs_minus_SE_df$date, y1=Babs_minus_SE_df$Babs_minus_SE_660 + Babs_minus_SE_df_error$Babs_minus_SE_660, code=3, angle=90, length=0.1, col="hotpink", lwd=1)
arrows(x0=Babs_minus_SE_df$date, y0=Babs_minus_SE_df$Babs_minus_SE_880 - Babs_minus_SE_df_error$Babs_minus_SE_880, 
       x1=Babs_minus_SE_df$date, y1=Babs_minus_SE_df$Babs_minus_SE_880 + Babs_minus_SE_df_error$Babs_minus_SE_880, code=3, angle=90, length=0.1, col="red", lwd=1)
arrows(x0=Babs_minus_SE_df$date, y0=Babs_minus_SE_df$Babs_minus_SE_950 - Babs_minus_SE_df_error$Babs_minus_SE_950, 
       x1=Babs_minus_SE_df$date, y1=Babs_minus_SE_df$Babs_minus_SE_950 + Babs_minus_SE_df_error$Babs_minus_SE_950, code=3, angle=90, length=0.1, col="red4", lwd=1)


legend("topright", inset =.02, box.lty = 0,
       legend = c("370 nm", "470 nm", "520 nm",
                  "590 nm", "660 nm", "880 nm", "950 nm"),
       pch = c(20,20,20,20,20,20,20),
       col = c("violet","blue","green","yellow","hotpink", "red", "red4"))

#############################################################################################################
# calculate SSA
#############################################################################################################

# remove terms from R environment except for:

rm(list=setdiff(ls(), c("Aet_TA", "ATTN_error", 
                        "Nep_TA", "Nep_TA_error_step3",
                        "Nep_TA_corrected", "Nep_TA_error_step4", 
                        "A_df", "AE_df",
                        "Bsca_fitted_df", "Nep_TA_error_step5",
                        "convert_list_df",
                        "SE_correction_df", "SE_correction_df_error",
                        "Aet_WL", "BATTN",
                        "Babs_minus_SE_df", "Babs_minus_SE_df_error")))


# calculate SSA at the Aet WL
SSA_list <- list()

for (i in 2:length(Bsca_fitted_df)) {
  SSA_list[[i]] <- Bsca_fitted_df[,i]/(Babs_minus_SE_df[,i] + Bsca_fitted_df[,i])
  SSA_list[[i]] <- data.frame(SSA_list[[i]])
  colnames(SSA_list[[i]]) <- sprintf("SSA_%d", Aet_WL[i-1])
}

# drop empty terms in the list
SSA_list <- SSA_list[lengths(SSA_list) != 0]

# bind the results together
SSA_df <- list.cbind(SSA_list)

# add the date column 
SSA_df <- cbind(Nep_TA_corrected$date, SSA_df)
SSA_df <- data.frame(SSA_df)
colnames(SSA_df)[1] <- "date"


# calculate the error related to SSA

# first calculate Bext and the error related to Bext, Bext = Bsca + Babs

# calculate Bext
Bext_df <- SSA_df
colnames(Bext_df) <- c("date", Aet_WL)

for (i in 2:length(Babs_minus_SE_df)) {
  Bext_df[,i] <- Babs_minus_SE_df[,i] + Bsca_fitted_df[,i]
}

# calculate Bext error (Bsca (predicted) error + Babs* error)
Bext_df_error <- Bext_df

for (i in 2:length(Babs_minus_SE_df)) {
  Bext_df_error[,i] <- sqrt((Babs_minus_SE_df_error[,i])^2 + (Nep_TA_error_step5[,i])^2)
}

# calculate the error related to SSA (Bext error * Bsca (predicted) error)
SSA_df_error <- SSA_df 
for (i in 2:length(SSA_df)) {
  SSA_df_error[,i] <- SSA_df[,i]*sqrt((Bext_df_error[,i]/Bext_df[,i])^2 + (Nep_TA_error_step5[,i]/Bsca_fitted_df[,i])^2)
}

#############################################################################################################
# plot SSA at all WL
#############################################################################################################

par(mfrow= c(1,1))

plot(SSA_df$date, SSA_df$SSA_370,
     ylab = "SSA",
     xlab = "Time",
     ylim = c(0,1.5),
     main = bquote(bold("Maeli2, single scattering albedo - SSA"*(lambda))),
     col="violet", pch = 20)
points(SSA_df$date, SSA_df$SSA_470, col="blue", pch = 20)
points(SSA_df$date, SSA_df$SSA_520, col="green", pch = 20)
points(SSA_df$date, SSA_df$SSA_590, col="yellow", pch = 20)
points(SSA_df$date, SSA_df$SSA_660, col="hotpink", pch = 20)
points(SSA_df$date, SSA_df$SSA_880, col="red", pch = 20)
points(SSA_df$date, SSA_df$SSA_950, col="red4", pch = 20)

arrows(x0=SSA_df$date, y0=SSA_df$SSA_370 - SSA_df_error$SSA_370, 
       x1=SSA_df$date, y1=SSA_df$SSA_370 + SSA_df_error$SSA_370, code=3, angle=90, length=0.1, col="violet", lwd=1)
arrows(x0=SSA_df$date, y0=SSA_df$SSA_470 - SSA_df_error$SSA_470, 
       x1=SSA_df$date, y1=SSA_df$SSA_470 + SSA_df_error$SSA_470, code=3, angle=90, length=0.1, col="blue", lwd=1)
arrows(x0=SSA_df$date, y0=SSA_df$SSA_520 - SSA_df_error$SSA_520, 
       x1=SSA_df$date, y1=SSA_df$SSA_520 + SSA_df_error$SSA_520, code=3, angle=90, length=0.1, col="green", lwd=1)
arrows(x0=SSA_df$date, y0=SSA_df$SSA_590 - SSA_df_error$SSA_590, 
       x1=SSA_df$date, y1=SSA_df$SSA_590 + SSA_df_error$SSA_590, code=3, angle=90, length=0.1, col="yellow", lwd=1)
arrows(x0=SSA_df$date, y0=SSA_df$SSA_660 - SSA_df_error$SSA_660, 
       x1=SSA_df$date, y1=SSA_df$SSA_660 + SSA_df_error$SSA_660, code=3, angle=90, length=0.1, col="hotpink", lwd=1)
arrows(x0=SSA_df$date, y0=SSA_df$SSA_880 - SSA_df_error$SSA_880, 
       x1=SSA_df$date, y1=SSA_df$SSA_880 + SSA_df_error$SSA_880, code=3, angle=90, length=0.1, col="red", lwd=1)
arrows(x0=SSA_df$date, y0=SSA_df$SSA_950 - SSA_df_error$SSA_950, 
       x1=SSA_df$date, y1=SSA_df$SSA_950 + SSA_df_error$SSA_950, code=3, angle=90, length=0.1, col="red4", lwd=1)

legend("bottomright", inset =.02, box.lty = 0,
       legend = c("370 nm", "470 nm", "520 nm",
                  "590 nm", "660 nm", "880 nm", "950 nm"),
       pch = c(20,20,20,20,20,20,20),
       col = c("violet","blue","green","yellow","hotpink", "red", "red4"))

#############################################################################################################
# calculate f(WL) at the Aet WL
#############################################################################################################

# f(WL)  = a(1-SSA(WL))+1

a_parameter <- 0.74

f_WL_list <- list()

for (i in 2:length(SSA_df)) {
  f_WL_list[[i]] <- a_parameter*(1 - SSA_df[,i]) + 1  
  f_WL_list[[i]] <- data.frame(f_WL_list[[i]])
  colnames(f_WL_list[[i]]) <- sprintf("f_WL_%d", Aet_WL[i-1])
}

# drop empty terms in the list
f_WL_list <- f_WL_list[lengths(f_WL_list) != 0]

# bind the results together
f_WL_df <- list.cbind(f_WL_list)

# add the date column 
f_WL_df <- cbind(Nep_TA_corrected$date,f_WL_df)
f_WL_df <- data.frame(f_WL_df)
colnames(f_WL_df)[1] <- "date"

# calculate the error related to f_WL (a_parameter*SSA error)
f_WL_df_error <- f_WL_df 
for (i in 2:length(SSA_df_error)) {
  f_WL_df_error[,i] <- SSA_df_error[,i]*a_parameter
}

#############################################################################################################
# plot f_WL at the Aet WL
#############################################################################################################

par(mfrow= c(1,1))

plot(f_WL_df$date, f_WL_df$f_WL_370,
     ylab = bquote(f),
     xlab = "Time",
     main = bquote(bold("Maeli2,"~f(lambda))),
     col="violet", pch = 20, ylim = c(0,2))
points(f_WL_df$date, f_WL_df$f_WL_470, col="blue", pch = 20)
points(f_WL_df$date, f_WL_df$f_WL_520, col="green", pch = 20)
points(f_WL_df$date, f_WL_df$f_WL_590, col="yellow", pch = 20)
points(f_WL_df$date, f_WL_df$f_WL_660, col="hotpink", pch = 20)
points(f_WL_df$date, f_WL_df$f_WL_880, col="red", pch = 20)
points(f_WL_df$date, f_WL_df$f_WL_950, col="red4", pch = 20)

arrows(x0=f_WL_df$date, y0=f_WL_df$f_WL_370 - f_WL_df_error$f_WL_370, 
       x1=f_WL_df$date, y1=f_WL_df$f_WL_370 + f_WL_df_error$f_WL_370, code=3, angle=90, length=0.1, col="violet", lwd=1)
arrows(x0=f_WL_df$date, y0=f_WL_df$f_WL_470 - f_WL_df_error$f_WL_470, 
       x1=f_WL_df$date, y1=f_WL_df$f_WL_470 + f_WL_df_error$f_WL_470, code=3, angle=90, length=0.1, col="blue", lwd=1)
arrows(x0=f_WL_df$date, y0=f_WL_df$f_WL_520 - f_WL_df_error$f_WL_520, 
       x1=f_WL_df$date, y1=f_WL_df$f_WL_520 + f_WL_df_error$f_WL_520, code=3, angle=90, length=0.1, col="green", lwd=1)
arrows(x0=f_WL_df$date, y0=f_WL_df$f_WL_590 - f_WL_df_error$f_WL_590, 
       x1=f_WL_df$date, y1=f_WL_df$f_WL_590 + f_WL_df_error$f_WL_590, code=3, angle=90, length=0.1, col="yellow", lwd=1)
arrows(x0=f_WL_df$date, y0=f_WL_df$f_WL_660 - f_WL_df_error$f_WL_660, 
       x1=f_WL_df$date, y1=f_WL_df$f_WL_660 + f_WL_df_error$f_WL_660, code=3, angle=90, length=0.1, col="hotpink", lwd=1)
arrows(x0=f_WL_df$date, y0=f_WL_df$f_WL_880 - f_WL_df_error$f_WL_880, 
       x1=f_WL_df$date, y1=f_WL_df$f_WL_880 + f_WL_df_error$f_WL_880, code=3, angle=90, length=0.1, col="red", lwd=1)
arrows(x0=f_WL_df$date, y0=f_WL_df$f_WL_950 - f_WL_df_error$f_WL_950, 
       x1=f_WL_df$date, y1=f_WL_df$f_WL_950 + f_WL_df_error$f_WL_950, code=3, angle=90, length=0.1, col="red4", lwd=1)


legend("bottomright", inset =.02, box.lty = 0,
       legend = c("370 nm", "470 nm", "520 nm",
                  "590 nm", "660 nm", "880 nm", "950 nm"),
       pch = c(20,20,20,20,20,20,20),
       col = c("violet","blue","green","yellow","hotpink", "red", "red4"))


#############################################################################################################
# calculate the LE_correction at the Aet WL
#############################################################################################################
# remove terms from R environment except for:
rm(list=setdiff(ls(), c("Aet_TA", "ATTN_error", 
                        "Nep_TA", "Nep_TA_error_step3",
                        "Nep_TA_corrected", "Nep_TA_error_step4", 
                        "A_df", "AE_df",
                        "Bsca_fitted_df", "Nep_TA_error_step5",
                        "convert_list_df",
                        "SE_correction_df", "SE_correction_df_error",
                        "Aet_WL", 
                        "Babs_minus_SE_df", "Babs_minus_SE_df_error",
                        "SSA_df", "SSA_df_error",
                        "f_WL_df", "f_WL_df_error",
                        "Bext_df", "Bext_df_error")))

# C2010 Method = ((1/f(WL)) - 1) * (ATT(WL)/50%) + 1

# get the ATTN data
ATTN <- Aet_TA[,c(1:8)]

# subset also ATTN error data
ATTN_error <- subset(ATTN_error, ATTN_error$date <= ATTN$date[nrow(ATTN)])

LE_correction_list <- list()

for (i in 2:length(f_WL_df)) {
  LE_correction_list[[i]] <- ((1/f_WL_df[,i]) - 1) * ((ATTN[,i]*100)/50) + 1 
  LE_correction_list[[i]] <- data.frame(LE_correction_list[[i]])
  colnames(LE_correction_list[[i]]) <- sprintf("LE_corr_%d", Aet_WL[i-1])
}

# drop empty terms in the list
LE_correction_list <- LE_correction_list[lengths(LE_correction_list) != 0]

# bind the results together
LE_correction_df <- list.cbind(LE_correction_list)

# add the date column 
LE_correction_df <- cbind(Nep_TA_corrected$date,LE_correction_df)
LE_correction_df <- data.frame(LE_correction_df)
colnames(LE_correction_df)[1] <- "date"

# calculate the error related to LE_correction (f_WL error*ATTN error)
LE_correction_df_error <- LE_correction_df 
for (i in 2:length(f_WL_df)) {
  LE_correction_df_error[,i] <- LE_correction_df[,i]*sqrt((f_WL_df_error[,i]/f_WL_df[,i])^2 + (ATTN_error[,i]/ATTN[,i])^2)
}

#############################################################################################################
# plot LE_corr at the Aet WL
#############################################################################################################


par(mfrow= c(1,1))

plot(LE_correction_df$date, LE_correction_df$LE_corr_370,
     ylab = bquote(R),
     xlab = "Time",
     ylim = c(0,1.6),
     main = bquote(bold("Maeli2, loading effect correction -"~R(lambda))),
     col="violet", pch = 20)

points(LE_correction_df$date, LE_correction_df$LE_corr_470, col="blue", pch = 20)
points(LE_correction_df$date, LE_correction_df$LE_corr_520, col="green", pch = 20)
points(LE_correction_df$date, LE_correction_df$LE_corr_590, col="yellow", pch = 20)
points(LE_correction_df$date, LE_correction_df$LE_corr_660, col="hotpink", pch = 20)
points(LE_correction_df$date, LE_correction_df$LE_corr_880, col="red", pch = 20)
points(LE_correction_df$date, LE_correction_df$LE_corr_950, col="red4", pch = 20)

arrows(x0=LE_correction_df$date, y0=LE_correction_df$LE_corr_370 - LE_correction_df_error$LE_corr_370, 
       x1=LE_correction_df$date, y1=LE_correction_df$LE_corr_370 + LE_correction_df_error$LE_corr_370, code=3, angle=90, length=0.1, col="violet", lwd=1)
arrows(x0=LE_correction_df$date, y0=LE_correction_df$LE_corr_470 - LE_correction_df_error$LE_corr_470, 
       x1=LE_correction_df$date, y1=LE_correction_df$LE_corr_470 + LE_correction_df_error$LE_corr_470, code=3, angle=90, length=0.1, col="blue", lwd=1)
arrows(x0=LE_correction_df$date, y0=LE_correction_df$LE_corr_520 - LE_correction_df_error$LE_corr_520, 
       x1=LE_correction_df$date, y1=LE_correction_df$LE_corr_520 + LE_correction_df_error$LE_corr_520, code=3, angle=90, length=0.1, col="green", lwd=1)
arrows(x0=LE_correction_df$date, y0=LE_correction_df$LE_corr_590 - LE_correction_df_error$LE_corr_590, 
       x1=LE_correction_df$date, y1=LE_correction_df$LE_corr_590 + LE_correction_df_error$LE_corr_590, code=3, angle=90, length=0.1, col="yellow", lwd=1)
arrows(x0=LE_correction_df$date, y0=LE_correction_df$LE_corr_660 - LE_correction_df_error$LE_corr_660, 
       x1=LE_correction_df$date, y1=LE_correction_df$LE_corr_660 + LE_correction_df_error$LE_corr_660, code=3, angle=90, length=0.1, col="hotpink", lwd=1)
arrows(x0=LE_correction_df$date, y0=LE_correction_df$LE_corr_880 - LE_correction_df_error$LE_corr_880, 
       x1=LE_correction_df$date, y1=LE_correction_df$LE_corr_880 + LE_correction_df_error$LE_corr_880, code=3, angle=90, length=0.1, col="red", lwd=1)
arrows(x0=LE_correction_df$date, y0=LE_correction_df$LE_corr_950 - LE_correction_df_error$LE_corr_950, 
       x1=LE_correction_df$date, y1=LE_correction_df$LE_corr_950 + LE_correction_df_error$LE_corr_950, code=3, angle=90, length=0.1, col="red4", lwd=1)


legend("bottomright", inset =.02, box.lty = 0,
       legend = c("370 nm", "470 nm", "520 nm",
                  "590 nm", "660 nm", "880 nm", "950 nm"),
       pch = c(20,20,20,20,20,20,20),
       col = c("violet","blue","green","yellow","hotpink", "red", "red4"))

#############################################################################################################
#############################################################################################################
#############################################################################################################
# Calculate Babs for all the Aet WL
#############################################################################################################
#############################################################################################################
#############################################################################################################

# remove terms from R environment except for:
rm(list=setdiff(ls(), c("Aet_TA", "ATTN_error", 
                        "Nep_TA", "Nep_TA_error_step3",
                        "Nep_TA_corrected", "Nep_TA_error_step4", 
                        "A_df", "AE_df",
                        "Bsca_fitted_df", "Nep_TA_error_step5",
                        "convert_list_df",
                        "SE_correction_df", "SE_correction_df_error",
                        "Aet_WL", 
                        "Babs_minus_SE_df", "Babs_minus_SE_df_error",
                        "SSA_df", "SSA_df_error",
                        "f_WL_df", "f_WL_df_error",
                        "Bext_df", "Bext_df_error",
                        "LE_correction_df", "LE_correction_df_error")))


# Babs = (BATTN(WL) - SE_correction(WL)*Bsca(WL))/LE_correction(WL)*Cref(WL)

# extrapolate the Cref at all the Aet WL using linear regression
# x = WL, y = Cref (Di Biagio et al. 2017 for Sahel/Niger)
xy <- data.frame(x = c(450, 660), y = c(4.18,3.83))
model_linear <- lm(y ~ x, data = xy)
new_x <- data.frame(x = c(370, 470, 520, 590, 660, 880, 950))
C_ref <- predict(model_linear, new_x)
C_ref_df <- data.frame(WL = c(370, 470, 520, 590, 660, 880, 950), C_ref = C_ref )

# calculate Cref error (10% Di Biagio et al., 2019)
C_ref_df_error <- C_ref_df
C_ref_df_error$C_ref <- C_ref_df$C_ref*0.1

# get the BATTN values
BATTN <- Aet_TA[,c(1,9:15)]

# calculate Babs
Babs_list <- list()

for (i in 2:length(BATTN)) {
  Babs_list[[i]] <- (BATTN[,i] - (SE_correction_df[,i]*Bsca_fitted_df[,i]))/(LE_correction_df[,i]*C_ref_df[i-1,2])
  Babs_list[[i]] <- data.frame(Babs_list[[i]])
  colnames(Babs_list[[i]]) <- sprintf("Babs_%d", Aet_WL[i-1])
}

# drop empty terms in the list
Babs_list <- Babs_list[lengths(Babs_list) != 0]

# bind the results together
Babs_df <- list.cbind(Babs_list)

# add the date column 
Babs_df <- cbind(Nep_TA_corrected$date, Babs_df)
Babs_df <- data.frame(Babs_df)
colnames(Babs_df)[1] <- "date"

# calculate Babs error (SE_corr error)*(Bsca (predicted) error)*(LE_corr error)*(Cref error)
Babs_df_error <- Babs_df

for (i in 2:length(Babs_df)) {
  Babs_df_error[,i] <- Babs_df_error[,i]*sqrt((SE_correction_df_error[,i]/SE_correction_df[,i])^2 + 
                                              (Nep_TA_error_step5[,i]/Bsca_fitted_df[,i])^2 + 
                                              (LE_correction_df_error[,i]/LE_correction_df[,i])^2 +
                                              (C_ref_df_error[i-1,2]/C_ref_df[i-1,2])^2 )
}

#############################################################################################################
# plot Babs at the Aet WL
#############################################################################################################

par(mfrow= c(1,1))

plot(Babs_df$date, Babs_df$Babs_370,
     xlab = "Time",
     ylab = bquote(beta[abs]~(m^-1)),
     ylim = c(0, 0.0005),
     main = bquote(bold("Maeli2, absorption coefficient -"~beta[abs](lambda))),
     col="violet", pch = 20)

points(Babs_df$date, Babs_df$Babs_470, col="blue", pch = 20)
points(Babs_df$date, Babs_df$Babs_520, col="green", pch = 20)
points(Babs_df$date, Babs_df$Babs_590, col="yellow", pch = 20)
points(Babs_df$date, Babs_df$Babs_660, col="hotpink", pch = 20)
points(Babs_df$date, Babs_df$Babs_880, col="red", pch = 20)
points(Babs_df$date, Babs_df$Babs_950, col="red4", pch = 20)

arrows(x0=Babs_df$date, y0=Babs_df$Babs_370 - Babs_df_error$Babs_370, 
       x1=Babs_df$date, y1=Babs_df$Babs_370 + Babs_df_error$Babs_370, code=3, angle=90, length=0.1, col="violet", lwd=1)
arrows(x0=Babs_df$date, y0=Babs_df$Babs_470 - Babs_df_error$Babs_470, 
       x1=Babs_df$date, y1=Babs_df$Babs_470 + Babs_df_error$Babs_470, code=3, angle=90, length=0.1, col="blue", lwd=1)
arrows(x0=Babs_df$date, y0=Babs_df$Babs_520 - Babs_df_error$Babs_520, 
       x1=Babs_df$date, y1=Babs_df$Babs_520 + Babs_df_error$Babs_520, code=3, angle=90, length=0.1, col="green", lwd=1)
arrows(x0=Babs_df$date, y0=Babs_df$Babs_590 - Babs_df_error$Babs_590, 
       x1=Babs_df$date, y1=Babs_df$Babs_590 + Babs_df_error$Babs_590, code=3, angle=90, length=0.1, col="yellow", lwd=1)
arrows(x0=Babs_df$date, y0=Babs_df$Babs_660 - Babs_df_error$Babs_660, 
       x1=Babs_df$date, y1=Babs_df$Babs_660 + Babs_df_error$Babs_660, code=3, angle=90, length=0.1, col="hotpink", lwd=1)
arrows(x0=Babs_df$date, y0=Babs_df$Babs_880 - Babs_df_error$Babs_880, 
       x1=Babs_df$date, y1=Babs_df$Babs_880 + Babs_df_error$Babs_880, code=3, angle=90, length=0.1, col="red", lwd=1)
arrows(x0=Babs_df$date, y0=Babs_df$Babs_950 - Babs_df_error$Babs_950, 
       x1=Babs_df$date, y1=Babs_df$Babs_950 + Babs_df_error$Babs_950, code=3, angle=90, length=0.1, col="red4", lwd=1)

legend("topright", inset =.02, box.lty = 0,
       legend = c("370 nm", "470 nm", "520 nm",
                  "590 nm", "660 nm", "880 nm", "950 nm"),
       pch = c(20,20,20,20,20,20,20),
       col = c("violet","blue","green","yellow","hotpink", "red", "red4"))

#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################
# Iterative calculation - for Babs and SSA
#############################################################################################################
#############################################################################################################
# re-calculate the loading effect correction 
# SSA is now estimated from Babs_df and Bsca_fitted_df
#############################################################################################################
#############################################################################################################

# Method from C2010 and Di Biagio et al. 2019

# LE_correction = ((1/f(WL)) - 1) * (ATT(WL)/50%) + 1

# both methods require the term f(WL) which is function of the single scattering albedo (SSA):

# f(WL)  = a(1-SSA(WL))+1

# where a is the linear fit between BATTN(WL)(X) and ATTN(WL)(y)

#############################################################################################################
# re-calculate SSA at all WL
#############################################################################################################

SSA_df_new <- SSA_df

for (i in 2:length(SSA_df_new)) {
  SSA_df_new[,i] <- Bsca_fitted_df[,i]/(Babs_df[,i] + Bsca_fitted_df[,i])
}

# re-calculate the error related to SSA

# first re-calculate Bext and the error related to Bext, Bext = Bsca + Babs

# re-calculate Bext
Bext_df_new <- SSA_df_new
colnames(Bext_df_new) <- c("date", Aet_WL)

for (i in 2:length(Babs_df)) {
  Bext_df_new[,i] <- Babs_df[,i] + Bsca_fitted_df[,i]
}

# re-calculate Bext error (Bsca (predicted) error + Babs error)
Bext_df_error_new <- Bext_df_new

for (i in 2:length(Babs_df)) {
  Bext_df_error_new[,i] <- sqrt((Babs_df_error[,i])^2 + (Nep_TA_error_step5[,i])^2)
}

# calculate the error related to SSA (Bext error * Bsca (predicted) error)
SSA_df_error_new <- SSA_df_new 
for (i in 2:length(SSA_df_new)) {
  SSA_df_error_new[,i] <- SSA_df_new[,i]*sqrt((Bext_df_error_new[,i]/Bext_df_new[,i])^2 + (Nep_TA_error_step5[,i]/Bsca_fitted_df[,i])^2)
}

#############################################################################################################
# plot SSA at all WL
#############################################################################################################

par(mfrow= c(1,1))

plot(SSA_df_new$date, SSA_df_new$SSA_370,
     ylab = "SSA",
     xlab = "Time",
     ylim = c(0,1.5),
     main = bquote(bold("Maeli2, single scattering albedo - SSA"*(lambda))),
     col="violet", pch = 20)
points(SSA_df_new$date, SSA_df_new$SSA_470, col="blue", pch = 20)
points(SSA_df_new$date, SSA_df_new$SSA_520, col="green", pch = 20)
points(SSA_df_new$date, SSA_df_new$SSA_590, col="yellow", pch = 20)
points(SSA_df_new$date, SSA_df_new$SSA_660, col="hotpink", pch = 20)
points(SSA_df_new$date, SSA_df_new$SSA_880, col="red", pch = 20)
points(SSA_df_new$date, SSA_df_new$SSA_950, col="red4", pch = 20)

arrows(x0=SSA_df_new$date, y0=SSA_df_new$SSA_370 - SSA_df_error_new$SSA_370, 
       x1=SSA_df_new$date, y1=SSA_df_new$SSA_370 + SSA_df_error_new$SSA_370, code=3, angle=90, length=0.1, col="violet", lwd=1)
arrows(x0=SSA_df_new$date, y0=SSA_df_new$SSA_470 - SSA_df_error_new$SSA_470, 
       x1=SSA_df_new$date, y1=SSA_df_new$SSA_470 + SSA_df_error_new$SSA_470, code=3, angle=90, length=0.1, col="blue", lwd=1)
arrows(x0=SSA_df_new$date, y0=SSA_df_new$SSA_520 - SSA_df_error_new$SSA_520, 
       x1=SSA_df_new$date, y1=SSA_df_new$SSA_520 + SSA_df_error_new$SSA_520, code=3, angle=90, length=0.1, col="green", lwd=1)
arrows(x0=SSA_df_new$date, y0=SSA_df_new$SSA_590 - SSA_df_error_new$SSA_590, 
       x1=SSA_df_new$date, y1=SSA_df_new$SSA_590 + SSA_df_error_new$SSA_590, code=3, angle=90, length=0.1, col="yellow", lwd=1)
arrows(x0=SSA_df_new$date, y0=SSA_df_new$SSA_660 - SSA_df_error_new$SSA_660, 
       x1=SSA_df_new$date, y1=SSA_df_new$SSA_660 + SSA_df_error_new$SSA_660, code=3, angle=90, length=0.1, col="hotpink", lwd=1)
arrows(x0=SSA_df_new$date, y0=SSA_df_new$SSA_880 - SSA_df_error_new$SSA_880, 
       x1=SSA_df_new$date, y1=SSA_df_new$SSA_880 + SSA_df_error_new$SSA_880, code=3, angle=90, length=0.1, col="red", lwd=1)
arrows(x0=SSA_df_new$date, y0=SSA_df_new$SSA_950 - SSA_df_error_new$SSA_950, 
       x1=SSA_df_new$date, y1=SSA_df_new$SSA_950 + SSA_df_error_new$SSA_950, code=3, angle=90, length=0.1, col="red4", lwd=1)

legend("bottomright", inset =.02, box.lty = 0,
       legend = c("370 nm", "470 nm", "520 nm",
                  "590 nm", "660 nm", "880 nm", "950 nm"),
       pch = c(20,20,20,20,20,20,20),
       col = c("violet","blue","green","yellow","hotpink", "red", "red4"))

#############################################################################################################
# re-calculate f(WL) at the Aet WL
#############################################################################################################

# re-calculate f(WL)  = a(1-SSA(WL))+1

a_parameter <- 0.74

f_WL_df_new <- f_WL_df

for (i in 2:length(SSA_df)) {
  f_WL_df_new[,i] <- a_parameter*(1 - SSA_df_new[,i]) + 1  
}

# re-calculate the error related to f_WL (a_parameter*SSA error)
f_WL_df_error_new <- f_WL_df_new 
for (i in 2:length(SSA_df_error_new)) {
  f_WL_df_error_new[,i] <- SSA_df_error_new[,i]*a_parameter
}

#############################################################################################################
# plot f_WL at the Aet WL
#############################################################################################################

par(mfrow= c(1,1))

plot(f_WL_df_new$date, f_WL_df_new$f_WL_370,
     ylab = bquote(f),
     xlab = "Time",
     ylim = c(0,1.3),
     main = bquote(bold("Maeli2,"~f(lambda))),
     col="violet", pch = 20)

points(f_WL_df_new$date, f_WL_df_new$f_WL_470, col="blue", pch = 20)
points(f_WL_df_new$date, f_WL_df_new$f_WL_520, col="green", pch = 20)
points(f_WL_df_new$date, f_WL_df_new$f_WL_590, col="yellow", pch = 20)
points(f_WL_df_new$date, f_WL_df_new$f_WL_660, col="hotpink", pch = 20)
points(f_WL_df_new$date, f_WL_df_new$f_WL_880, col="red", pch = 20)
points(f_WL_df_new$date, f_WL_df_new$f_WL_950, col="red4", pch = 20)

arrows(x0=f_WL_df_new$date, y0=f_WL_df_new$f_WL_370 - f_WL_df_error_new$f_WL_370, 
       x1=f_WL_df_new$date, y1=f_WL_df_new$f_WL_370 + f_WL_df_error_new$f_WL_370, code=3, angle=90, length=0.1, col="violet", lwd=1)
arrows(x0=f_WL_df_new$date, y0=f_WL_df_new$f_WL_470 - f_WL_df_error_new$f_WL_470, 
       x1=f_WL_df_new$date, y1=f_WL_df_new$f_WL_470 + f_WL_df_error_new$f_WL_470, code=3, angle=90, length=0.1, col="blue", lwd=1)
arrows(x0=f_WL_df_new$date, y0=f_WL_df_new$f_WL_520 - f_WL_df_error_new$f_WL_520, 
       x1=f_WL_df_new$date, y1=f_WL_df_new$f_WL_520 + f_WL_df_error_new$f_WL_520, code=3, angle=90, length=0.1, col="green", lwd=1)
arrows(x0=f_WL_df_new$date, y0=f_WL_df_new$f_WL_590 - f_WL_df_error_new$f_WL_590, 
       x1=f_WL_df_new$date, y1=f_WL_df_new$f_WL_590 + f_WL_df_error_new$f_WL_590, code=3, angle=90, length=0.1, col="yellow", lwd=1)
arrows(x0=f_WL_df_new$date, y0=f_WL_df_new$f_WL_660 - f_WL_df_error_new$f_WL_660, 
       x1=f_WL_df_new$date, y1=f_WL_df_new$f_WL_660 + f_WL_df_error_new$f_WL_660, code=3, angle=90, length=0.1, col="hotpink", lwd=1)
arrows(x0=f_WL_df_new$date, y0=f_WL_df_new$f_WL_880 - f_WL_df_error_new$f_WL_880, 
       x1=f_WL_df_new$date, y1=f_WL_df_new$f_WL_880 + f_WL_df_error_new$f_WL_880, code=3, angle=90, length=0.1, col="red", lwd=1)
arrows(x0=f_WL_df_new$date, y0=f_WL_df_new$f_WL_950 - f_WL_df_error_new$f_WL_950, 
       x1=f_WL_df_new$date, y1=f_WL_df_new$f_WL_950 + f_WL_df_error_new$f_WL_950, code=3, angle=90, length=0.1, col="red4", lwd=1)

legend("bottomright", inset =.02, box.lty = 0,
       legend = c("370 nm", "470 nm", "520 nm",
                  "590 nm", "660 nm", "880 nm", "950 nm"),
       pch = c(20,20,20,20,20,20,20),
       col = c("violet","blue","green","yellow","hotpink", "red", "red4"))

#############################################################################################################
# re-calculate the LE_correction at the Aet WL
#############################################################################################################

# remove terms from R environment except for:
rm(list=setdiff(ls(), c("Aet_TA", "ATTN_error", 
                        "Nep_TA", "Nep_TA_error_step3",
                        "Nep_TA_corrected", "Nep_TA_error_step4", 
                        "A_df", "AE_df",
                        "Bsca_fitted_df", "Nep_TA_error_step5",
                        "convert_list_df",
                        "SE_correction_df", "SE_correction_df_error",
                        "Aet_WL", 
                        "SSA_df_new", "SSA_df_error_new",
                        "f_WL_df_new", "f_WL_df_error_new",
                        "Bext_df_new", "Bext_df_error_new",
                        "LE_correction_df", "LE_correction_df_error",
                        "Babs_df", "Babs_df_error",
                        "C_ref_df", "C_ref_df_error")))


# C2010 Method = ((1/f(WL)) - 1) * (ATT(WL)/50%) + 1

# get the ATTN data
ATTN <- Aet_TA[,c(1:8)]

# re-calculate LE_correction
LE_correction_df_new <- LE_correction_df

for (i in 2:length(LE_correction_df_new)) {
  LE_correction_df_new[,i] <- ((1/f_WL_df_new[,i]) - 1) * ((ATTN[,i]*100)/50) + 1 
}

# re-calculate the error related to LE_correction (f_WL error*ATTN error)
LE_correction_df_error_new <- LE_correction_df_new 
for (i in 2:length(f_WL_df_new)) {
  LE_correction_df_error_new[,i] <- LE_correction_df_new[,i]*sqrt((f_WL_df_error_new[,i]/f_WL_df_new[,i])^2 + (ATTN_error[,i]/ATTN[,i])^2)
}

#############################################################################################################
# plot LE_corr at the Aet WL
#############################################################################################################

par(mfrow= c(1,1))

plot(LE_correction_df_new$date, LE_correction_df_new$LE_corr_370,
     ylab = bquote(R),
     xlab = "Time",
     ylim = c(0,1.6),
     main = bquote(bold("Maeli2, loading effect correction -"~R(lambda))),
     col="violet", pch = 20)

points(LE_correction_df_new$date, LE_correction_df_new$LE_corr_470, col="blue", pch = 20)
points(LE_correction_df_new$date, LE_correction_df_new$LE_corr_520, col="green", pch = 20)
points(LE_correction_df_new$date, LE_correction_df_new$LE_corr_590, col="yellow", pch = 20)
points(LE_correction_df_new$date, LE_correction_df_new$LE_corr_660, col="hotpink", pch = 20)
points(LE_correction_df_new$date, LE_correction_df_new$LE_corr_880, col="red", pch = 20)
points(LE_correction_df_new$date, LE_correction_df_new$LE_corr_950, col="red4", pch = 20)

arrows(x0=LE_correction_df_new$date, y0=LE_correction_df_new$LE_corr_370 - LE_correction_df_error_new$LE_corr_370, 
       x1=LE_correction_df_new$date, y1=LE_correction_df_new$LE_corr_370 + LE_correction_df_error_new$LE_corr_370, code=3, angle=90, length=0.1, col="violet", lwd=1)
arrows(x0=LE_correction_df_new$date, y0=LE_correction_df_new$LE_corr_470 - LE_correction_df_error_new$LE_corr_470, 
       x1=LE_correction_df_new$date, y1=LE_correction_df_new$LE_corr_470 + LE_correction_df_error_new$LE_corr_470, code=3, angle=90, length=0.1, col="blue", lwd=1)
arrows(x0=LE_correction_df_new$date, y0=LE_correction_df_new$LE_corr_520 - LE_correction_df_error_new$LE_corr_520, 
       x1=LE_correction_df_new$date, y1=LE_correction_df_new$LE_corr_520 + LE_correction_df_error_new$LE_corr_520, code=3, angle=90, length=0.1, col="green", lwd=1)
arrows(x0=LE_correction_df_new$date, y0=LE_correction_df_new$LE_corr_590 - LE_correction_df_error_new$LE_corr_590, 
       x1=LE_correction_df_new$date, y1=LE_correction_df_new$LE_corr_590 + LE_correction_df_error_new$LE_corr_590, code=3, angle=90, length=0.1, col="yellow", lwd=1)
arrows(x0=LE_correction_df_new$date, y0=LE_correction_df_new$LE_corr_660 - LE_correction_df_error_new$LE_corr_660, 
       x1=LE_correction_df_new$date, y1=LE_correction_df_new$LE_corr_660 + LE_correction_df_error_new$LE_corr_660, code=3, angle=90, length=0.1, col="hotpink", lwd=1)
arrows(x0=LE_correction_df_new$date, y0=LE_correction_df_new$LE_corr_880 - LE_correction_df_error_new$LE_corr_880, 
       x1=LE_correction_df_new$date, y1=LE_correction_df_new$LE_corr_880 + LE_correction_df_error_new$LE_corr_880, code=3, angle=90, length=0.1, col="red", lwd=1)
arrows(x0=LE_correction_df_new$date, y0=LE_correction_df_new$LE_corr_950 - LE_correction_df_error_new$LE_corr_950, 
       x1=LE_correction_df_new$date, y1=LE_correction_df_new$LE_corr_950 + LE_correction_df_error_new$LE_corr_950, code=3, angle=90, length=0.1, col="red4", lwd=1)

legend("bottomright", inset =.02, box.lty = 0,
       legend = c("370 nm", "470 nm", "520 nm",
                  "590 nm", "660 nm", "880 nm", "950 nm"),
       pch = c(20,20,20,20,20,20,20),
       col = c("violet","blue","green","yellow","hotpink", "red", "red4"))

#############################################################################################################
#############################################################################################################
# re-calculate Babs for all the Aet WL
#############################################################################################################
#############################################################################################################

# remove terms from R environment except for:
rm(list=setdiff(ls(), c("Aet_TA", "ATTN_error", 
                        "Nep_TA", "Nep_TA_error_step3",
                        "Nep_TA_corrected", "Nep_TA_error_step4", 
                        "A_df", "AE_df",
                        "Bsca_fitted_df", "Nep_TA_error_step5",
                        "convert_list_df",
                        "SE_correction_df", "SE_correction_df_error",
                        "Aet_WL",
                        "SSA_df_new", "SSA_df_error_new",
                        "f_WL_df_new", "f_WL_df_error_new",
                        "Bext_df_new", "Bext_df_error_new",
                        "LE_correction_df_new", "LE_correction_df_error_new",
                        "Babs_df", "Babs_df_error",
                        "C_ref_df", "C_ref_df_error")))

# Babs = (BATTN(WL) - SE_correction(WL)*Bsca(WL))/LE_correction(WL)*Cref(WL)

# get the BATTN values
BATTN <- Aet_TA[,c(1,9:15)]

# calculate Babs
Babs_df_new <- Babs_df

for (i in 2:length(Babs_df_new)) {
  Babs_df_new[,i] <- (BATTN[,i] - (SE_correction_df[,i]*Bsca_fitted_df[,i]))/(LE_correction_df_new[,i]*C_ref_df[i-1,2])
}

# re-calculate Babs error (SE_corr error)*(Bsca (predicted) error)*(LE_corr error)*(Cref error)
Babs_df_error_new <- Babs_df_new

for (i in 2:length(Babs_df)) {
  Babs_df_error_new[,i] <- Babs_df_error_new[,i]*sqrt((SE_correction_df_error[,i]/SE_correction_df[,i])^2 + 
                                                (Nep_TA_error_step5[,i]/Bsca_fitted_df[,i])^2 + 
                                                (LE_correction_df_error_new[,i]/LE_correction_df_new[,i])^2 +
                                                (C_ref_df_error[i-1,2]/C_ref_df[i-1,2])^2 )
}

#############################################################################################################
# plot Babs at the Aet WL
#############################################################################################################

par(mfrow= c(1,1))

plot(Babs_df_new$date, Babs_df_new$Babs_370,
     xlab = "Time",
     ylab = bquote(beta[abs]~(m^-1)),
     ylim = c(0, 0.0005),
     main = bquote(bold("Maeli2, absorption coefficient -"~beta[abs](lambda))),
     col="violet", pch = 20)

points(Babs_df_new$date, Babs_df_new$Babs_470, col="blue", pch = 20)
points(Babs_df_new$date, Babs_df_new$Babs_520, col="green", pch = 20)
points(Babs_df_new$date, Babs_df_new$Babs_590, col="yellow", pch = 20)
points(Babs_df_new$date, Babs_df_new$Babs_660, col="hotpink", pch = 20)
points(Babs_df_new$date, Babs_df_new$Babs_880, col="red", pch = 20)
points(Babs_df_new$date, Babs_df_new$Babs_950, col="red4", pch = 20)

arrows(x0=Babs_df_new$date, y0=Babs_df_new$Babs_370 - Babs_df_error_new$Babs_370, 
       x1=Babs_df_new$date, y1=Babs_df_new$Babs_370 + Babs_df_error_new$Babs_370, code=3, angle=90, length=0.1, col="violet", lwd=1)
arrows(x0=Babs_df_new$date, y0=Babs_df_new$Babs_470 - Babs_df_error_new$Babs_470, 
       x1=Babs_df_new$date, y1=Babs_df_new$Babs_470 + Babs_df_error_new$Babs_470, code=3, angle=90, length=0.1, col="blue", lwd=1)
arrows(x0=Babs_df_new$date, y0=Babs_df_new$Babs_520 - Babs_df_error_new$Babs_520, 
       x1=Babs_df_new$date, y1=Babs_df_new$Babs_520 + Babs_df_error_new$Babs_520, code=3, angle=90, length=0.1, col="green", lwd=1)
arrows(x0=Babs_df_new$date, y0=Babs_df_new$Babs_590 - Babs_df_error_new$Babs_590, 
       x1=Babs_df_new$date, y1=Babs_df_new$Babs_590 + Babs_df_error_new$Babs_590, code=3, angle=90, length=0.1, col="yellow", lwd=1)
arrows(x0=Babs_df_new$date, y0=Babs_df_new$Babs_660 - Babs_df_error_new$Babs_660, 
       x1=Babs_df_new$date, y1=Babs_df_new$Babs_660 + Babs_df_error_new$Babs_660, code=3, angle=90, length=0.1, col="hotpink", lwd=1)
arrows(x0=Babs_df_new$date, y0=Babs_df_new$Babs_880 - Babs_df_error_new$Babs_880, 
       x1=Babs_df_new$date, y1=Babs_df_new$Babs_880 + Babs_df_error_new$Babs_880, code=3, angle=90, length=0.1, col="red", lwd=1)
arrows(x0=Babs_df_new$date, y0=Babs_df_new$Babs_950 - Babs_df_error_new$Babs_950, 
       x1=Babs_df_new$date, y1=Babs_df_new$Babs_950 + Babs_df_error_new$Babs_950, code=3, angle=90, length=0.1, col="red4", lwd=1)

legend("topright", inset =.02, box.lty = 0,
       legend = c("370 nm", "470 nm", "520 nm",
                  "590 nm", "660 nm", "880 nm", "950 nm"),
       pch = c(20,20,20,20,20,20,20),
       col = c("violet","blue","green","yellow","hotpink", "red", "red4"))

#############################################################################################################
#############################################################################################################
# save all the parameters used for Aet correction
#############################################################################################################
#############################################################################################################

# calculate rsd and tidy up dfs before saving the output as .csv

# tidy A, AE data
A_df$rsd <- A_df$A_error/A_df$A*100
colnames(A_df) <- c("date", "A", "sd", "rsd")
AE_df$rsd <- AE_df$AE_error/AE_df$AE*100
colnames(AE_df) <- c("date", "AE", "sd", "rsd")

# tidy SE_correction data
SE_correction_df_rsd <- SE_correction_df_error

for (i in 2:length(SE_correction_df_error)) {
  SE_correction_df_rsd[,i] <- SE_correction_df_error[,i]/SE_correction_df[,i]*100
  colnames(SE_correction_df_rsd)[i] <- sprintf("rsd_%d", Aet_WL[i-1])
  colnames(SE_correction_df_error)[i] <- sprintf("sd_%d", Aet_WL[i-1])
}

SE_correction <- cbind(SE_correction_df, SE_correction_df_error[,-1], SE_correction_df_rsd[,-1])

# tidy LE_correction data
LE_correction_df_rsd_new <- LE_correction_df_error_new

for (i in 2:length(LE_correction_df_error_new)) {
  LE_correction_df_rsd_new[,i] <- LE_correction_df_error_new[,i]/LE_correction_df_new[,i]*100
  colnames(LE_correction_df_rsd_new)[i] <- sprintf("rsd_%d", Aet_WL[i-1])
  colnames(LE_correction_df_error_new)[i] <- sprintf("sd_%d", Aet_WL[i-1])
}

LE_correction <- cbind(LE_correction_df_new, LE_correction_df_error_new[,-1], LE_correction_df_rsd_new[,-1])

# tidy f_WL
f_WL_df_rsd_new <- f_WL_df_error_new

for (i in 2:length(f_WL_df_error_new)) {
  f_WL_df_rsd_new[,i] <- f_WL_df_error_new[,i]/f_WL_df_new[,i]*100
  colnames(f_WL_df_rsd_new)[i] <- sprintf("rsd_%d", Aet_WL[i-1])
  colnames(f_WL_df_error_new)[i] <- sprintf("sd_%d", Aet_WL[i-1])
}

f_WL <- cbind(f_WL_df_new, f_WL_df_error_new[,-1], f_WL_df_rsd_new[,-1])

# save the output
write.csv(A_df, "Z:/Clarissa/Data_Optical_calculation/R_Code/Code_final/ICEdust_samples/Maeli2/Aet_corrections/original/Maeli2_Aet_corr_parameter_A_original.csv", row.names = F)
write.csv(AE_df, "Z:/Clarissa/Data_Optical_calculation/R_Code/Code_final/ICEdust_samples/Maeli2/Aet_corrections/original/Maeli2_Aet_corr_parameter_AE_original.csv", row.names = F)
write.csv(SE_correction, "Z:/Clarissa/Data_Optical_calculation/R_Code/Code_final/ICEdust_samples/Maeli2/Aet_corrections/original/Maeli2_Aet_corr_parameter_SE_correction_original.csv", row.names = F)
write.csv(LE_correction, "Z:/Clarissa/Data_Optical_calculation/R_Code/Code_final/ICEdust_samples/Maeli2/Aet_corrections/original/Maeli2_Aet_corr_parameter_LE_correction_original.csv", row.names = F)
write.csv(f_WL, "Z:/Clarissa/Data_Optical_calculation/R_Code/Code_final/ICEdust_samples/Maeli2/Aet_corrections/original/Maeli2_Aet_corr_parameter_f_WL_original.csv", row.names = F)

# tidy Babs
# convert Babs as Mm-1
for (i in 2:length(Babs_df_new)) {
  Babs_df_new[,i] <- Babs_df_new[,i]*10^6
}

for (i in 2:length(Babs_df_error_new)) {
  Babs_df_error_new[,i] <- Babs_df_error_new[,i]*10^6
}

Babs_df_rsd_new <- Babs_df_new

for (i in 2:length(Babs_df_error_new)) {
  Babs_df_rsd_new[,i] <- Babs_df_error_new[,i]/Babs_df_new[,i]*100
  colnames(Babs_df_rsd_new)[i] <- sprintf("rsd_%d", Aet_WL[i-1])
  colnames(Babs_df_error_new)[i] <- sprintf("sd_%d", Aet_WL[i-1])
}

Babs <- cbind(Babs_df_new, Babs_df_error_new[,-1], Babs_df_rsd_new[,-1])

# save the output
write.csv(Babs, "Z:/Clarissa/Data_Optical_calculation/R_Code/Code_final/ICEdust_samples/Maeli2/Aet_corrections/original/Maeli2_Aet_Babs_2min_original.csv", row.names = F)

