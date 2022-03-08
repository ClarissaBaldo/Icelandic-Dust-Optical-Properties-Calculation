#############################################################################################################
# set working directory and upload packages
#############################################################################################################

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

# import the absorption coefficient (Babs) and the scattering coefficient (Bsca) at 12-min resolution
Maeli2_Babs_data <- read.csv("Z:/Clarissa/Data_Optical_calculation/R_Code/Code_final/ICEdust_samples/Maeli2/Nep_Aet_Mie/original/data/Maeli2_Babs_12min_original.csv")
Maeli2_Bsca_data <- read.csv("Z:/Clarissa/Data_Optical_calculation/R_Code/Code_final/ICEdust_samples/Maeli2/Nep_Aet_Mie/original/data/Maeli2_Bsca_fitted_12min_original.csv")

# subset the data from 30 min to 2.5 hours after the dust injection peak
Maeli2_Babs_data$date  <- as.POSIXct(Maeli2_Babs_data$date  ,  format =  "%Y-%m-%d %H:%M:%S")
Maeli2_Bsca_data$date  <- as.POSIXct(Maeli2_Bsca_data$date  ,  format =  "%Y-%m-%d %H:%M:%S")

Maeli2_Babs_data <- Maeli2_Babs_data[which(Maeli2_Babs_data$date >= "2019-01-21 11:31:00" & Maeli2_Babs_data$date <= "2019-01-21 13:31:00"),]
Maeli2_Bsca_data <- Maeli2_Bsca_data[which(Maeli2_Bsca_data$date >= "2019-01-21 11:31:00" & Maeli2_Bsca_data$date <= "2019-01-21 13:31:00"),]

# split data and uncertainty (error/sd)
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

# drop all rows with NA values
for (i in 1:length(data_fitting)) {
  data_fitting[[i]] <- na.omit(data_fitting[[i]])
}

###########################################################################################################
###########################################################################################################
# calculate the average single scattering albedo (SSA) over the duration of the experiment
###########################################################################################################
###########################################################################################################

# exp-average SSA = (1 + 1/mRMA)^-1 
# Moosmuller et al. (2012) (10.1029/2011jd016909) and Di Biagio et al. (2019) (10.5194/acp-19-15503-2019)
# mRMA = slope of the linear regression between Bsca and Babs

###########################################################################################################
# Reduced Major Axis Regression (RMA)
###########################################################################################################

# packages needed
library(lmodel2)
library(dplyr)

# perform RMA regression
data_RMA <- data_fitting

# lmodel2(y ~ x)
# both x and y have errors, but these are assumed to be constant and not input
# this is why we use RMA, but you must set x for the obs you trust more (x for the so-called "ground-truth")

# perform RMA regression
RMA_fit <- list()
for (i in 1:length(data_RMA)) {
  RMA_fit [[i]] <- lmodel2(Bsca ~ Babs, data = data_RMA[[i]])
}

###########################################################################################################
# save the results of the RMA regression
###########################################################################################################

# save the results
SSA_RMA <- list()
for (i in 1:length(RMA_fit)) {
  SSA_RMA[[i]] <- cbind(RMA_fit[[i]][["regression.results"]][,-c(4:5)], RMA_fit[[i]][["confidence.intervals"]][,-1], RMA_fit[[i]][["rsquare"]])
}

# calculate the uncertainty on mRMA using the estimated confidence intervals
calculate_RMA_sd <- list()

for (i in 1:length(SSA_RMA)) {
  calculate_RMA_sd[[i]] <- data.frame(test1 = abs(SSA_RMA[[i]]$`2.5%-Slope`[3] - SSA_RMA[[i]]$Slope[3]),
                                      test2 = abs(SSA_RMA[[i]]$`97.5%-Slope`[3] - SSA_RMA[[i]]$Slope[3]))
}

# select the largest uncertainty
SSA_RMA_sd <- c()
for (i in 1:length(calculate_RMA_sd)) {
  if (calculate_RMA_sd[[i]]$test1 > calculate_RMA_sd[[i]]$test2) {
    SSA_RMA_sd[i] <- calculate_RMA_sd[[i]]$test1
  } else  {
    SSA_RMA_sd[i] <- calculate_RMA_sd[[i]]$test2
  }
}

# get mRMA for different wavelengths (WL)
SSA_RMA_slope <- c()
for (i in 1:length(SSA_RMA)) {
  SSA_RMA_slope[i] <- SSA_RMA[[i]]$Slope[3]
}

# to give column names
Aet_WL <- c(370, 470, 520, 590, 660, 880, 950)

# create a df with mRMA and related uncertainty
SSA_RMA_data <- data.frame("WL" = Aet_WL,
                           "m" = SSA_RMA_slope, 
                           "m_sd" = SSA_RMA_sd,
                           "m_rsd" = SSA_RMA_sd/SSA_RMA_slope*100)

# calculate the average SSA over the duration of the experiment and related uncertainty
# exp-average SSA = (1 + 1/mRMA)^-1
# save the results and related uncertainty
SSA_average_RMA <- SSA_RMA_data
SSA_average_RMA$SSA <- (1 + 1/SSA_average_RMA$m)^(-1)
SSA_average_RMA$SSA_sd <- SSA_average_RMA$SSA*SSA_average_RMA$m_sd/(SSA_average_RMA$m)
SSA_average_RMA$SSA_rsd <- SSA_average_RMA$SSA_sd/SSA_average_RMA$SSA*100

# plot average SSA
library(ggplot2)
ggplot()+
  
  geom_point(data = SSA_average_RMA, aes(x = WL , y = SSA), color = "black", shape = 19, size=3, stroke = 2)+
  geom_line(data = SSA_average_RMA, aes(x = WL , y = SSA), color = "black", linetype = "dashed")+
  geom_errorbar(data = SSA_average_RMA, aes(x = WL, y = SSA, ymin = SSA - SSA_sd, ymax= SSA + SSA_sd), color = "black", width = 10) +
  
  ggtitle(bquote(bold("Experiment-averaged single scattering albedo - "*SSA[avg](lambda))))+
  
  xlab(expression(lambda~(nm)))+
  ylab(expression("SSA"))+
  scale_y_continuous(labels = function(x) format(x, scientific = F), limits =  c(0.6,1.05))+
  
  theme_bw() +
  theme(plot.title = element_text(size = 18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text.align = 0,
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        text = element_text(size=18),
        axis.text = element_text(size=18, colour = "black"),
        legend.text = element_text(size=18),
        legend.title = element_blank(),
        legend.direction = "vertical",
        legend.position = c(0.8, 0.3),
        legend.key.size = unit(0.5, "cm"))


