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

######################################################################################################
######################################################################################################
######################################################################################################
# GRIMM (GRM) corrections
######################################################################################################
######################################################################################################
######################################################################################################

######################################################################################################
# import the files for correction from optical diameter (Dop) to geometric diameter (Dg)
######################################################################################################

setwd("Z:/Clarissa/Data_Optical_Calculation/SIZE")

# n = 1.57-1.63, for k use 0.001 step resolution from 0 to 0.02
k_sequence <- seq(0.000000, 0.02, by = 0.001)

# create file_name seq to use as pattern
file_names <- list()
for (i in 1:length(k_sequence)) {
  file_names[[i]] <- sprintf("_%s0", k_sequence[i])
}

# for k = 0 use pattern as follows
file_names[[1]] <- "_0.000000i"

# create a list using the selected pattern
files_list <- list()

for (i in 1:length(file_names)){
  files_list[[i]] <- list.files(pattern = file_names[[i]]) 
}

# flatten the nested list
files <- unlist(files_list)

# remove terms from R environment except for:
rm(list=setdiff(ls(), c("files","convert_list_df")))


# import the selected correction files
Datafile <- list()

for (i in 1:length(files)){
  Datafile[[i]] <- read.csv(files[[i]], skip = 3, header = F, sep = "")
  colnames(Datafile[[i]]) <- c("DiameterCal", "DiameterCor")
}

# save the correction files, n = 1.57-1.63, k =0.000-0.02
GRM_Dg <- Datafile

######################################################################################################
# replace those diameters giving negative dlogDg with the neighbours average values
######################################################################################################

# keep only the column with the corrected Dg
for (i in 1:length(GRM_Dg)) {
  GRM_Dg[[i]] <- GRM_Dg[[i]][,2]  
}

# replace problematic diameters
for (i in 1:length(GRM_Dg)) {
  for (j in 2:(length(GRM_Dg[[i]])-1)) {
    if (log10(GRM_Dg[[i]][j]/GRM_Dg[[i]][j-1]) < 0.01) {
      GRM_Dg[[i]][j] <- (GRM_Dg[[i]][j+1] + GRM_Dg[[i]][j-1])/2
    } 
    else {} 
  }
}

######################################################################################################
# calculate the new midpoints from Dg
######################################################################################################

# Calculate the geometric midpoint diameters (midpoint Dg) = (D1*D2)^1/2

# create a list of empty array
midpoint_GRM <- list()

for (i in 1:length(GRM_Dg)){
  midpoint_GRM[[i]] <- array()
}

# there are in total 31 new Dg
# for the last Dg, D2 does not exist so the forloop is up to k-1 terms
# round to 6 decimals
for (i in 1:length(GRM_Dg)){
  for (k in 1:30){
    midpoint_GRM[[i]][k] <-(GRM_Dg[[i]][k+1]*GRM_Dg[[i]][k])^(1/2)
    midpoint_GRM[[i]][k] <- round(midpoint_GRM[[i]][k],6)
  }
}

######################################################################################################
# calculate the dlogDg
######################################################################################################

# create a list of empty array
GRM_dlogDg <- list()

for (i in 1:length(GRM_Dg)){
  GRM_dlogDg[[i]] <- array()
}

# there are in total 31 new Dg
# for the last Dg, D2 does not exist so the forloop is up to k-1 terms
for (i in 1:length(GRM_Dg)){
  for (k in 1:30){
    GRM_dlogDg[[i]][k] <- log10(GRM_Dg[[i]][k+1]/GRM_Dg[[i]][k])
  }
}

######################################################################################################
# Import GRM data and tidy up the dataset
######################################################################################################

GRM_data <- read.csv("Z:/Clarissa/Data_Optical_calculation/CLA_DATA/GRIMM/Maeli2_GRM.csv", header = T)

# Convert the date as POSIXct
GRM_data$Date <- as.POSIXct(GRM_data$Date,  format =  "%Y-%m-%d %H:%M:%S")

# Round to minutes
GRM_data$Date <- strptime(GRM_data$Date,  format = "%Y-%m-%d %H:%M")
GRM_data$Date <- as.POSIXct(GRM_data$Date,  format =  "%Y-%m-%d %H:%M:%S")

colnames(GRM_data)[1] <- "date"

######################################################################################################
# GRM - TimeAverage (TA) calculation
######################################################################################################
# GRM and SMPS must have the same time
# I want to transfrom both in 12 min average

# check the date-time columnames is "date"

# select the date since the sampling has started
# Note for Maeli2 the starting sampling time is 11:07, 
GRM_data_subset <- subset(GRM_data, GRM_data$date >= 	"2019-01-21 11:07:00")

# verify the presence of duplicates in the row sequence
GRM_data_subset <- GRM_data_subset[which(duplicated(GRM_data_subset$date) == F),]

#############################################################################################################
# function to calculate the time average at different time resoluiton (here 12min)

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

#############################################################################################################
# apply the function to calculate TA
GRM_data_TA <- timeAverage_12min(GRM_data_subset)

row.names(GRM_data_TA) <- c(1:nrow(GRM_data_TA))

######################################################################################################
# conversion from optical diameter to geometric diameter
######################################################################################################

# replace the new midpoint Dg in the original df (GRM_data_TA)
GRM_data_conv <- list()

for (i in 1:length(midpoint_GRM)){
  GRM_data_conv[[i]] <- GRM_data_TA
  colnames(GRM_data_conv[[i]])[-1] <- midpoint_GRM[[i]]
}

# dN does not change, just replace the optical diameters with the new midpoint Dg

######################################################################################################
# calculate dN/dlogDg for GRIMM data
######################################################################################################

# for each scenario (df in GRM_data_conv), divide the bin-sizes (midpoints-columns) by the corresponding calculated dlogDg value
GRM_data_dNdlogDg <- GRM_data_conv

for (i in 1:length(GRM_data_conv)){
  for (j in 2:length(GRM_data_conv[[i]])) {
    GRM_data_dNdlogDg[[i]][,j] <- GRM_data_conv[[i]][,j]/GRM_dlogDg[[i]][j]
  }
}

# dN does not change, just replace the optical diameters with the new midpoint Dg
######################################################################################################
######################################################################################################
######################################################################################################
# SMPS corrections
######################################################################################################
######################################################################################################
######################################################################################################

######################################################################################################
# Import SMPS data and tidy up the dataset
######################################################################################################

SMPS_data <- read.csv("Z:/Clarissa/Data_Optical_calculation/CLA_DATA/SMPS/SMPS_Maeli2.csv", skip = 1, header = T)

# Combine Date and StartTime
# Check and convert the format of the date first
SMPS_data$Date <- as.POSIXct(SMPS_data$Date,format = "%m/%d/%Y")

# Combine date and time
SMPS_data$DateTime <- as.POSIXct(paste(SMPS_data$Date, SMPS_data$Start.Time), format="%Y-%m-%d %H:%M:%S")

# Round to minutes
SMPS_data$DateTime <- strptime(SMPS_data$DateTime,  format = "%Y-%m-%d %H:%M")
SMPS_data$DateTime <- as.POSIXct(SMPS_data$DateTime,  format =  "%Y-%m-%d %H:%M:%S")

# Put the DateTime column at the front and delete useless columns
SMPS_data <- SMPS_data[, c(111, 1:110)]
SMPS_data <- SMPS_data[,-(2:3)]
colnames(SMPS_data)[1] <- "date"

# remove the last column with total dN 
SMPS_data <- SMPS_data[-length(SMPS_data)]

######################################################################################################
# SMPS - TimeAverage calculation
######################################################################################################
# GRIMM and SMPS must have the same time
# I want to transfrom both in 12 min average

# check the date-time columnames is "date"

# select the date since the sampling has started
SMPS_data_subset <- subset(SMPS_data, SMPS_data$date >= 	"2019-01-21 11:07:00")

# verify the presence of duplicates in the row sequence
SMPS_data_subset <- SMPS_data_subset[which(duplicated(SMPS_data_subset$date) == F),]

# apply the function to calculate TA
SMPS_data_TA <- timeAverage_12min(SMPS_data_subset)

row.names(SMPS_data_TA) <- c(1:nrow(SMPS_data_TA))

# GRM and SMPS data must have the same time interval
SMPS_data_TA <- subset(SMPS_data_TA, SMPS_data_TA$date >= GRM_data_TA[1,1] & SMPS_data_TA$date <= GRM_data_TA[nrow(GRM_data_TA),1])

######################################################################################################
# get the mobility diameter (Dm)
######################################################################################################

# for Dm = 881.7 nm there are no measurement records
# remove the last column corresponding to Dm = 881.7 nm
SMPS_data_TA <- SMPS_data_TA[-length(SMPS_data_TA)]

# SMPS Dm in nm
SMPS_Dm <- colnames(SMPS_data_TA)[-1]

for (i in 1:length(SMPS_Dm)) {
  SMPS_Dm[i] <- unlist(strsplit(SMPS_Dm[i], "X"))[2]
}

SMPS_Dm <- as.numeric(as.character(SMPS_Dm))

######################################################################################################
# conversion from mobility diameter (Dm) to geometric diameter (Dg)
######################################################################################################

# transform Dm into Dg
# Dg = Dm/dyn_shape_factor * Cun_Dg/Cun_Dm

# the dynamic shape factor is our variable

# set the dynamic shape factor from 1.6 to 2
dyn_shape_factor <- c(1,seq(1.6, 2, by = 0.1))

# create a list to convert Dm into Dg according to the different dynamic shape factors
SMPS_Dm_list <- list()

for (i in 1:length(dyn_shape_factor)) {
  SMPS_Dm_list[[i]] <- SMPS_Dm
}

######################################################################################################
# calculate the Cunningham slip correction factor (Cun)
######################################################################################################
# the Cunningham slip correction factor is calculated assuming a mean free path of 66 nm (air condition) 

# Dp (particle diameter) > 100 nm, 1 + (66/Dp)*2.52
# Dp < 100 nm, 1 + (66/Dp)*(2.34 + 1.05*exp(-0.39*Dp/66))

# define a function to calculate the Cun correction factor
convert_test <- function(x){
  convert_out <- x
  for (i in 1:length(x)){
    for (j in 1:length(x[[i]])){
      if (x[[i]][j] > 100) {
        
        convert_out[[i]][j] <- 1 + ((66/x[[i]][j])*2.52)
        
      } else {
        
        convert_out[[i]][j] <- 1 + (66/x[[i]][j])*(2.34 + 1.05*exp(-0.39*x[[i]][j]/66))
      }
    }
  }
  return(convert_out)
}

# calculate the correction factor for Dm (Cun_Dm)
Cun_Dm <- convert_test(SMPS_Dm_list)

# transform Dm into Dg
# the equation is "Dg = Dm/dyn_shape_factor * Cun_Dg/Cun_Dm"
# this means for each Dm, find the optimal Dg which best meets the condition above
# that means we need to find the root for "Dm/dyn_shape_factor * Cun_Dg/Cun_Dm - Dg = 0"

# create the structure for the list of the SMPS corrected gemetric diameters, but the values being NA
SMPS_Dg <- SMPS_Dm_list

for (i in 1:length(SMPS_Dg)){
  for(j in 1:length(SMPS_Dg[[i]])){
    SMPS_Dg[[i]][j] <- NA
  }
}

# calculate Dg
# the loop reads each single Dm, calculate the function and only save the estimated Dg values
for (i in 1:length(SMPS_Dm_list)){
  for (j in 1:length(SMPS_Dm_list[[i]])){
    x <- SMPS_Dg[[i]][[j]]
    f_test <- function(x) (SMPS_Dm_list[[i]][j]/dyn_shape_factor[[i]]*(convert_test(x)/Cun_Dm[[i]][j]) - x)
    SMPS_Dg[[i]][j] <- print(uniroot(f_test, lower=0, upper=1e+09, extendInt = "yes")$root)
  }
}

# check the output
# for dynamic shape factor = 1
# the result should be equal to the original Dm values
round(SMPS_Dg[[1]] - SMPS_Dm_list[[1]])

# feed any resulting x in the function, the output should be close to 0
i = 2
j = 71

f_test <- function(x) (SMPS_Dm_list[[i]][j]/dyn_shape_factor[[i]]*(convert_test(x)/Cun_Dm[[i]][j])- x)

x = SMPS_Dg[[i]][j]
f_test(x) 

# convert Dg in um
for (i in 1:length(SMPS_Dg)) {
  SMPS_Dg[[i]] <- SMPS_Dg[[i]]*0.001
}

######################################################################################################
# replace the new midpoint Dg in the original dataset (SMPS_data)
######################################################################################################

SMPS_data_conv <- list()

# dN does not change, just replace Dm with Dg in um
for (i in 1:length(SMPS_Dg)) {
  SMPS_data_conv[[i]] <- SMPS_data_TA
  colnames(SMPS_data_conv[[i]])[-1] <- round(SMPS_Dg[[i]],6)
}

######################################################################################################
# calculate dN/dlogDg for SMPS data
######################################################################################################

# for each scenario (df in SMPS_data_conv), divide the bin-sizes (midpoints-columns) by the corresponding calculated dlogDg value
SMPS_data_dNdlogDg <- SMPS_data_conv

for (i in 1:length(SMPS_data_conv)){
  for (j in 2:length(SMPS_data_conv[[i]])) {
    SMPS_data_dNdlogDg[[i]][,j] <- SMPS_data_conv[[i]][,j]*64
  }
}

######################################################################################################
######################################################################################################
######################################################################################################
# Merging the datasets
######################################################################################################
######################################################################################################
######################################################################################################

# remove terms from R environment except for:
rm(list=setdiff(ls(), c("GRM_data_conv", "SMPS_data_conv", "files", "dyn_shape_factor", "GRM_data_TA", "SMPS_data_TA", 
                        "midpoint_GRM", "SMPS_Dg", "convert_list_df", "GRM_dlogDg", "SMPS_data_dNdlogDg", "GRM_data_dNdlogDg")))

######################################################################################################
# get the details for the n-k scenarios
######################################################################################################

# extrapolate the n and k values from the "files" list containing GRM corrections 
GRM_nk <- strsplit(files,"_") # divide the string at the delimiter "_"
n <- list()
k <- list()

for (i in 1:length(GRM_nk)) {
  n[[i]] <- GRM_nk[[i]][5] # select #element in the cut string
  k[[i]] <- strsplit(GRM_nk[[i]][6], "i")# divide the string at the "i"
}

for (i in 1:length(k)) {
  k[[i]] <- k[[i]][[1]][1] # select #element in the cut string
}

# the list of n-k scenarios is copied # times = # X scenarios
nk_scenario <- list()
for (i in 1:length(SMPS_data_conv)) {
  nk_scenario[[i]] <- list()
}

for (i in 1:length(SMPS_data_conv)) {
  for (j in 1:length(n)) {
    nk_scenario[[i]][[j]] <- matrix(data = 1, nrow = nrow(GRM_data_TA), ncol = 2)
    nk_scenario[[i]][[j]][,1] <- n[[j]]
    nk_scenario[[i]][[j]][,2] <- k[[j]]
    nk_scenario[[i]][[j]] <- data.frame(nk_scenario[[i]][[j]])
    colnames(nk_scenario[[i]][[j]]) <- c("n", "k")
  }
}

# check the output
identical(nk_scenario[[1]][[147]], nk_scenario[[2]][[147]])

######################################################################################################
# prepare GRM data for merging
######################################################################################################
# Each n-k scenario is copied # times = # X scenarios
GRM_copy <- list()
for (i in 1:length(SMPS_data_dNdlogDg)) {
  GRM_copy[[i]] <- list()
}

for (i in 1:length(SMPS_data_dNdlogDg)) {
  for (j in 1:length(GRM_data_dNdlogDg)) {
    GRM_copy[[i]][[j]] <- GRM_data_dNdlogDg[[j]]
  }
}

# the same for the GRM diameters
GRM_Dg_copy <- GRM_copy

for (i in 1:length(SMPS_data_dNdlogDg)) {
  for (j in 1:length(GRM_data_dNdlogDg)) {
    GRM_Dg_copy[[i]][[j]] <- midpoint_GRM[[j]]
  }
}

# the same for the dlogDg values
GRM_dlogDg_copy <- GRM_copy

for (i in 1:length(SMPS_data_dNdlogDg)) {
  for (j in 1:length(GRM_data_dNdlogDg)) {
    GRM_dlogDg_copy[[i]][[j]] <- data.frame(cbind(midpoint_GRM[[j]],GRM_dlogDg[[j]]))
    colnames(GRM_dlogDg_copy[[i]][[j]]) <- c("Dg", "dlogDg")
  }
}

# check the output
identical(GRM_copy[[1]][[147]], GRM_copy[[2]][[147]])
# check the output
identical(GRM_Dg_copy[[1]][[147]], GRM_Dg_copy[[5]][[147]])
# check the output
identical(GRM_dlogDg_copy[[1]][[147]][,2], GRM_dlogDg_copy[[5]][[147]][,2])

#####################################################################################################
# prepare SMPS data for merging
######################################################################################################
# Each X scenario is copied # times = # n-k scenarios
SMPS_copy <- list()
for (i in 1:length(SMPS_data_dNdlogDg)) {
  SMPS_copy[[i]] <- list()
}

for (i in 1:length(SMPS_data_dNdlogDg)) {
  for (j in 1:length(GRM_data_dNdlogDg)) {
    SMPS_copy[[i]][[j]] <- SMPS_data_dNdlogDg[[i]]
  }
}

# add a column with the X values for the different scenarios
for (i in 1:length(SMPS_copy)) {
  for (j in 1:length(SMPS_copy[[i]])) {
    SMPS_copy[[i]][[j]] <- add_column(SMPS_copy[[i]][[j]], X = dyn_shape_factor[i],.before = "date")
  }
}

# check the output
identical(SMPS_copy[[1]][[1]], SMPS_copy[[1]][[147]])
identical(SMPS_copy[[1]][[1]], SMPS_copy[[5]][[1]])

# the same for the SMPS diameters
SMPS_Dg_copy <- SMPS_copy

for (i in 1:length(SMPS_copy)) {
  for (j in 1:length(SMPS_copy[[i]])) {
    SMPS_Dg_copy[[i]][[j]] <- SMPS_Dg[[i]]
  }
}

# check the output
identical(SMPS_Dg_copy[[1]][[1]], SMPS_Dg_copy[[1]][[147]])
identical(SMPS_Dg_copy[[1]][[1]], SMPS_Dg_copy[[2]][[1]])

######################################################################################################
# merge the data
######################################################################################################

# merge the scenario info, SMPS and GRM data
# we already selected the SMPS data for Dm >=  0.850 um
# the final size distribution will be equal to the SMPS dataset for Dm >=  0.850 um and GRIMM_Dg >= SMPS_Dg 

# build a function to merge single X_SMPS - GRM - nk scenario dfs
merge_SMPS_GRM <- function(X_SMPS, GRM, nk_scenario) {
  merge_df <- list()
  final_df <- list()
  for (i in 1:length(GRM)) {
    merge_df <- GRM[,which(colnames(GRM) >= colnames(X_SMPS)[length(X_SMPS)])]
    final_df <- data.frame(cbind(nk_scenario, X_SMPS, merge_df[,-1]))
  }
  return(final_df)
}

# apply the function to all SMPS - GRM - nk scenario dfs
merged_data <- SMPS_copy
for (i in 1:length(SMPS_copy)) {
  for (j in 1:length(SMPS_copy[[i]])) {
    merged_data[[i]][[j]] <- merge_SMPS_GRM(SMPS_copy[[i]][[j]], GRM_copy[[i]][[j]], nk_scenario[[i]][[j]])
  }
}

# now merge the diameters
Dg_merged <- merged_data

for (i in 1:length(merged_data)) {
  for (j in 1:length(merged_data[[i]])) {
    GRM_Dg_copy[[i]][[j]] <-  GRM_Dg_copy[[i]][[j]][which(GRM_Dg_copy[[i]][[j]] > max(SMPS_Dg_copy[[i]][[j]]))]
    Dg_merged[[i]][[j]] <- c(SMPS_Dg_copy[[i]][[j]], GRM_Dg_copy[[i]][[j]])
  }
}

# check the output
i = 1
j = 147

# the number of columns in the merged dataframe must be equal to the the length of the merged Dg
test <- merged_data[[i]][[j]][-c(1:4)]
length(test) == length(Dg_merged[[i]][[j]])

# each colnames should correspond to the values in merged_Dg  
k = 6
as.numeric(unlist(strsplit(colnames(test)[k], "X"))[2]) == round(Dg_merged[[i]][[j]][k],6)

######################################################################################################
# Transpose column and rows, as for each time you want to show the size distribution
######################################################################################################

# build a function to transpose df
transpose_df <- function(merged_df){
  TS_merged_df <- data.frame(t(merged_df[,-c(1:4)]))
  # Assign column and row names
  colnames(TS_merged_df) <- GRM_data_TA[,1]
  TS_merged_df  <- rownames_to_column(TS_merged_df, var="Dg")
  
  return(TS_merged_df) 
}

# apply the function to all the dfs in the list merged_data (merged SD)
TS_merged_data <- merged_data

for (i in 1:length(merged_data)) {
  for (j in 1:length(merged_data[[i]])) {
    TS_merged_data[[i]][[j]] <- transpose_df(merged_data[[i]][[j]])
  }
}

# check the output
i = 4
j = 17

# check the dimensions of the transposed df, they should be consistent with the original df
length(TS_merged_data[[i]][[j]]) - nrow(merged_data[[i]][[j]]) == 1
length(merged_data[[i]][[j]]) -4 == nrow(TS_merged_data[[i]][[j]])

# in the TS merged df the Dg values are as character
# this is because they are extrapolated from the colnames
# converting this to numeric can lead to errors
# Thus, replace the Dg columns with the values in the Dg_merged list

for (i in 1:length(Dg_merged)) {
  for (j in 1:length(Dg_merged[[i]])) {
    TS_merged_data[[i]][[j]]$Dg <- Dg_merged[[i]][[j]]
  }
}

######################################################################################################
######################################################################################################
######################################################################################################
# interpolate the data - test Calib, X= 1, n=1.59, k=0
######################################################################################################
######################################################################################################
######################################################################################################

# remove terms from R environment except for:
rm(list=setdiff(ls(), c("GRM_data_conv", "SMPS_data_conv", "files", "dyn_shape_factor", "GRM_data_TA", "SMPS_data_TA",
                        "merged_data", "TS_merged_data", "SMPS_Dg_copy", "GRM_Dg_copy", "convert_list_df", "GRM_dlogDg", 
                        "SMPS_data_dNdlogDg", "GRM_data_dNdlogDg", "GRM_dlogDg_copy")))

# identify the outliers based on the the fitted smooth spline

######################################################################################################
# Interpolate the merged Dg serie 
######################################################################################################

# Interpolate using a smooth spline (3-grade polynom)

# set the 0 values as NA
TS_merged_data_interpolate <- TS_merged_data[[1]][[3]]
TS_merged_data_interpolate[TS_merged_data_interpolate[] == 0] <- NA

# preavious tests showed that interpolating the data using a cubic spline
# can produce several peaks in the upper part of the size distribution
# thus, before interpolating the data we need to transform the dN/dlogDg values
# convert the data as ln(dN/dlogDg)

for (i in 2:length(TS_merged_data_interpolate )) {
  TS_merged_data_interpolate[,i] <- log(TS_merged_data_interpolate[,i])
}


# x is Dg before interpolation
x <- data.frame(TS_merged_data_interpolate[,1])

# y is the number concentration of particles at different times
y <- list()

for (i in 2:length(TS_merged_data_interpolate)){
  y[[i]] <- data.frame(TS_merged_data_interpolate[,i])
}

# drop the NULL terms in the list
y <- y[lengths(y) != 0]

# the smoothspline function does not support missing values
# so remove the missing values from each column
# create a list containing the Dg-dN pair 
# the NA values has been removed so each object in the nested list has different number of rows

xy <- list()

for (i in 1:length(y)) {
  xy[[i]] <- cbind(x[,1],y[[i]][,1])
  xy[[i]] <- na.omit(xy[[i]])
  xy[[i]] <- data.frame(xy[[i]])
  colnames(xy[[i]]) <- c("x", "y")
}

# generate the spline function
# all the observation data are true
# spline.res = smooth.spline(xy$x, xy$y, spar = 0, all.knots = TRUE)

# apply the spline function to each column in the original dataframe
# which correspond to the different dfs in the list
spline.res <- list()

for (i in 1:length(xy)) {
  spline.res[[i]] <- smooth.spline(xy[[i]]$x, xy[[i]]$y, spar = 1, all.knots = F)
}

######################################################################################################
# Build a new size distribution dlogDg = 1/64 & Interpolate the merged Dg serie - test
######################################################################################################

######################################################################################################
# Create a new diameter serie for interpolation

# D1 and D2 are two consecutive diameters
# dlogDg = log10(D2/D1) = 1/64 for SMPS
# midpoint^2 =  D1*D2

# D2 = D1 * 10^(1/64)
# midpoint^2 =  D1*D1 * 10^(1/64)
# midpoint^2 =  D1^2 * 10^(1/64)
# D1^2 = midpoint^2/10^(1/64)
# D2^2 = midpoint^2*10^(1/64)

# calculate the start and end Dg of the size distribution 
# use the first and last midpoint of the serie
midpoint_start <-  TS_merged_data[[1]][[3]][1,1]
D1 <- sqrt((midpoint_start^2)/10^(1/64))
midpoint_end <-  TS_merged_data[[1]][[3]][nrow(TS_merged_data[[1]][[3]]),1]
Dend <- sqrt((midpoint_end^2)*10^(1/64))

# generate the new series of Dg for the whole merged size distribution dlogDg = 1/64

# create an empty list
D1Dend_interpolation <- list()

# assign the starting point of the sequence/list as D1
D1Dend_interpolation[[1]] <- D1

# D2 = D1*10^(1/64); Di = D[[i-1]]*10^(1/64)
# calculate the new element based on the preavious one
# here we set the upper limit to be Dend
# if (D1Dend_interpolation[[i]] >= Dend) {break}

for (i in 2:10^6) {
  D1Dend_interpolation[[i]] <- (D1Dend_interpolation[[i-1]])*10^(1/64)
  if (D1Dend_interpolation[[i]] >= Dend) {
    break }
}

# calculate the new mid point diameters
# midpoint^2 =  D1*D2
midpoint_interpolation <- list()

for (i in 1:(length(D1Dend_interpolation)-1)) { # exclude the last term of the list
  midpoint_interpolation[[i]] <- sqrt(D1Dend_interpolation[[i]]*D1Dend_interpolation[[i+1]])
}

# convert the list into a sequence of numbers
midpoint_interpolation <- convert_list_df(midpoint_interpolation)
midpoint_interpolation <- unlist(midpoint_interpolation[,1])

######################################################################################################
# use a smooth line to interpolate the new Dg serie
######################################################################################################  

# apply the function spline.res to the new Dg 
spline_interpolation <- list()

for (i in 1:length(spline.res)){
  spline_interpolation[[i]] <- data.frame(predict(spline.res[[i]], midpoint_interpolation))
}

# the interpolated data will have different numbers of rows
# use the function full join to merge the data
# nrow is the maximum number of interpolated diameters (midpoint_interpolation)
# ncol is the maximum number of column in the original transposed dataframe

# create a new empty dataframe
Interpolation_output <- matrix(nrow = length(midpoint_interpolation), ncol = length(TS_merged_data[[1]][[3]]))

# the first column in the dataframe corresponds to the new midpoint serie
Interpolation_output[,1] <- midpoint_interpolation

# the full join function will return a dataframe with three columns
# the first column correponds to the new midpoint diameter 
# the second column and the third column are the data from two different joint time series
# you can only join two terms at once
# full join all the terms in the list, the first time interval is constant
# while the second term change
# in the full join output keep only the third column which is the variable [,3]

for (i in 1:length(spline_interpolation)) {
  Interpolation_output[,i+1] <- (full_join(spline_interpolation[[i]], spline_interpolation[[i]], by = "x"))[,3]
}

colnames(Interpolation_output) <- colnames(TS_merged_data[[1]][[3]])

Interpolation_output <- data.frame(Interpolation_output)

# check the output
i = 1
identical(spline_interpolation[[i]][,2], Interpolation_output[,i+1])
identical(Interpolation_output[,1], midpoint_interpolation)

######################################################################################################
# only consider the interpolated value starting/ending from/to valid measurements
######################################################################################################
# Identify the first row that is not NA for each column
# larger volumes > 10 um have dN = 0
# when interpolating, those values become negative or produce spikes
# this is not ok when normalization the data
# these problematic data have to be replaced as NA in the interpolation outputs

# Identify the first and last row that is not NA for each column in the original merged size distributions 

# prepare the data for corrections
correct_Interpolation_output <- list()

for (i in 2:length(Interpolation_output)) {
  correct_Interpolation_output[[i]] <- data.frame(Interpolation_output[,c(1,i)])
}

# from the merged Dg serie before interpolation identify the first and last Dg corresponding to valid measurements
valid_meas <- list()

for (i in 2:length(TS_merged_data_interpolate)) {
  valid_meas[[i]] <- data.frame(TS_merged_data_interpolate[,c(1,i)])
}

# remove the NA values
for (i in 2:length(valid_meas)) {
  valid_meas[[i]] <- na.omit(valid_meas[[i]])
}

# in the interpolation output 
# for each date time
# replace the value outside the Dg interval corresponding to valid measurements
for (i in 2:length(correct_Interpolation_output)) {
  for (j in 1:nrow(correct_Interpolation_output[[i]])) {
    if(correct_Interpolation_output[[i]]$Dg[j] <= valid_meas[[i]]$Dg[1] | correct_Interpolation_output[[i]]$Dg[j] >= valid_meas[[i]]$Dg[nrow(valid_meas[[i]])]){
      correct_Interpolation_output[[i]][j,2] <- NA }
    else {}
  }
}

# save the output 
Interpolation_output_NA <- Interpolation_output

for (i in 2:length(Interpolation_output_NA)) {
  Interpolation_output_NA[,i] <- correct_Interpolation_output[[i]][,2]
}

# now reconvert the values in dN/dlogDg
for (i in 2:length(Interpolation_output_NA)) {
  Interpolation_output_NA[,i] <- exp(Interpolation_output_NA[,i])
}

# also ensure all the negative values are as NA
for (i in 2:length(Interpolation_output_NA)) {
  for (j in 1:nrow(Interpolation_output_NA)) {
    if ((Interpolation_output_NA[j,i] < 0) == TRUE | is.na(Interpolation_output_NA[j,i]) == TRUE) {
      Interpolation_output_NA[j,i] = NA} 
    else{}
  }
}

######################################################################################################
# plot the data - interpolation
######################################################################################################

par(mfrow= c(2,2))

plot(TS_merged_data[[1]][[3]][,1],TS_merged_data[[1]][[3]][,2], log = "xy",
     xaxt = "n",
     yaxt = "n",
     xlab = expression(paste("Diameter (",mu,"m)")),
     ylab = expression(paste("dN/dlogDg (",cm^{-3}*")")),
     main= "dN/dlogDg at 11:07",
     xlim = c(0.01, 100),
     ylim = c(0.0005, 1500), pch = 20)
points(Interpolation_output_NA$Dg, Interpolation_output_NA[,2], col = "green", pch = 20)
mtext("n = 1.59, k = 0, X = 1", cex = 0.75, side=3,line=0, col = "red")
axis(1, at = c(0.1,1,10))
axis(2, at = c(0.1,1,10,100,1000), labels = c(0.1,1,10,100,1000))
legend("bottomleft", inset =.02, box.lty = 0,
       legend = c("Merged data", "Interpolated values"),
       pch = c(20,20,20),
       col = c("black", "green"), cex=0.75)


plot(TS_merged_data[[1]][[3]][,1],TS_merged_data[[1]][[3]][,4], log = "xy",
     xaxt = "n",
     yaxt = "n",
     xlab = expression(paste("Diameter (",mu,"m)")),
     ylab = expression(paste("dN/dlogDg (",cm^{-3}*")")),
     main= "dN/dlogDg at 11:31",
     xlim = c(0.01, 100),
     ylim = c(0.0005, 1500), pch = 20)
points(Interpolation_output_NA$Dg, Interpolation_output_NA[,4], col = "green", pch = 20)
mtext("n = 1.59, k = 0, X = 1", cex = 0.75, side=3,line=0, col = "red")
axis(1, at = c(0.1,1,10))
axis(2, at = c(0.1,1,10,100,1000), labels = c(0.1,1,10,100,1000))
legend("bottomleft", inset =.02, box.lty = 0,
       legend = c("Merged data", "Interpolated values"),
       pch = c(20,20,20),
       col = c("black", "green"), cex=0.75)


plot(TS_merged_data[[1]][[3]][,1],TS_merged_data[[1]][[3]][,6], log = "xy",
     xaxt = "n",
     yaxt = "n",
     xlab = expression(paste("Diameter (",mu,"m)")),
     ylab = expression(paste("dN/dlogDg (",cm^{-3}*")")),
     main= "dN/dlogDg at 11:55",
     xlim = c(0.01, 100),
     ylim = c(0.0005, 1500), pch = 20)
points(Interpolation_output_NA$Dg, Interpolation_output_NA[,6], col = "green", pch = 20)
mtext("n = 1.59, k = 0, X = 1", cex = 0.75, side=3,line=0, col = "red")
axis(1, at = c(0.1,1,10))
axis(2, at = c(0.1,1,10,100,1000), labels = c(0.1,1,10,100,1000))
legend("bottomleft", inset =.02, box.lty = 0,
       legend = c("Merged data", "Interpolated values"),
       pch = c(20,20,20),
       col = c("black", "green"), cex=0.75)


plot(TS_merged_data[[1]][[3]][,1],TS_merged_data[[1]][[3]][,9], log = "xy",
     xaxt = "n",
     yaxt = "n",
     xlab = expression(paste("Diameter (",mu,"m)")),
     ylab = expression(paste("dN/dlogDg (",cm^{-3}*")")),
     main= "dN/dlogDg at 12:31",
     xlim = c(0.01, 100),
     ylim = c(0.0005, 1500), pch = 20)
points(Interpolation_output_NA$Dg, Interpolation_output_NA[,9], col = "green", pch = 20)
mtext("n = 1.59, k = 0, X = 1", cex = 0.75, side=3,line=0, col = "red")
axis(1, at = c(0.1,1,10))
axis(2, at = c(0.1,1,10,100,1000), labels = c(0.1,1,10,100,1000))
legend("bottomleft", inset =.02, box.lty = 0,
       legend = c("Merged data", "Interpolated values"),
       pch = c(20,20,20),
       col = c("black", "green"), cex=0.75)


######################################################################################################
######################################################################################################
######################################################################################################
# build a function to interpolate the data for all df
######################################################################################################
######################################################################################################
######################################################################################################

interpolate_data <- function(TS_merged_df){
  # Interpolate the merged Dg serie 
  # set the 0 values as NA
  TS_merged_data_interpolate <- data.frame(TS_merged_df)
  TS_merged_data_interpolate[TS_merged_data_interpolate[] == 0] <- NA
  
  # convert the data as ln(dN/dlogDg)
  for (i in 2:length(TS_merged_data_interpolate)) {
    TS_merged_data_interpolate[,i] <- log(TS_merged_data_interpolate[,i])
  }
  
  # x is Dg before interpolation
  x <- data.frame(TS_merged_data_interpolate[,1])
  
  # y is the number concentration of particles at different times
  y <- list()
  
  for (i in 2:length(TS_merged_data_interpolate)){
    y[[i]] <- data.frame(TS_merged_data_interpolate[,i])
  }
  
  # drop the NULL terms in the list
  y <- y[lengths(y) != 0]
  
  # the smoothspline function does not support missing values
  # so remove the missing values from each column
  # create a list containing the diameter-dN pair 
  # the NA values has been removed so each object in the nested list has different number of rows
  
  xy <- list()
  
  for (i in 1:length(y)) {
    xy[[i]] <- cbind(x[,1],y[[i]][,1])
    xy[[i]] <- na.omit(xy[[i]])
    xy[[i]] <- data.frame(xy[[i]])
    colnames(xy[[i]]) <- c("x", "y")
  }
  
  # generate the spline function
  # all the observation data are true
  # spline.res = smooth.spline(xy$x, xy$y, spar = 0, all.knots = F)
  
  # apply the spline function to each column in the original dataframe
  # which correspond to the different dfs in the list
  spline.res <- list()
  
  for (i in 1:length(xy)) {
    spline.res[[i]] <- smooth.spline(xy[[i]]$x, xy[[i]]$y, spar = 1, all.knots = F)
  }
  
  ######################################################################################################
  # Create a new diameter serie for interpolation
  
  # D1 and D2 are two consecutive diameters
  # dlogDg = log10(D2/D1) = 1/64 for SMPS
  # midpoint^2 =  D1*D2
  
  # D2 = D1 * 10^(1/64)
  # midpoint^2 =  D1*D1 * 10^(1/64)
  # midpoint^2 =  D1^2 * 10^(1/64)
  # D1^2 = midpoint^2/10^(1/64)
  # D2^2 = midpoint^2*10^(1/64)
  
  # calculate the start and end Dg of the size distribution 
  # use the first and last midpoint of the serie
  midpoint_start <-  TS_merged_df[1,1]
  D1 <- sqrt((midpoint_start^2)/10^(1/64))
  midpoint_end <-  TS_merged_df[nrow(TS_merged_df),1]
  Dend <- sqrt((midpoint_end^2)*10^(1/64))
  
  # generate the new series of Dg for the whole merged size distribution dlogDg = 1/64
  
  # create an empty list
  D1Dend_interpolation <- list()
  
  # assign the starting point of the sequence/list as D1
  D1Dend_interpolation[[1]] <- D1
  
  # D2 = D1*10^(1/64); Di = D[[i-1]]*10^(1/64)
  # calculate the new element based on the preavious one
  # here we set the upper limit to be Dend
  # if (D1Dend_interpolation[[i]] >= Dend) {break}
  
  for (i in 2:10^6) {
    D1Dend_interpolation[[i]] <- (D1Dend_interpolation[[i-1]])*10^(1/64)
    if (D1Dend_interpolation[[i]] >= Dend) {
      break }
  }
  
  # calculate the new mid point diameters
  # midpoint^2 =  D1*D2
  midpoint_interpolation <- list()
  
  for (i in 1:(length(D1Dend_interpolation)-1)) { # exclude the last term of the list
    midpoint_interpolation[[i]] <- sqrt(D1Dend_interpolation[[i]]*D1Dend_interpolation[[i+1]])
  }
  
  # convert the list into a sequence of numbers
  midpoint_interpolation <- convert_list_df(midpoint_interpolation)
  midpoint_interpolation <- unlist(midpoint_interpolation[,1])
  
  ######################################################################################################
  # use a smooth line to interpolate the new Dg serie
   
  # apply the function spline.res to the new Dg 
  spline_interpolation <- list()
  
  for (i in 1:length(spline.res)){
    spline_interpolation[[i]] <- data.frame(predict(spline.res[[i]], midpoint_interpolation))
  }
  
  # the interpolated data will have different numbers of rows
  # use the function full join to merge the data
  # nrow is the maximum number of interpolated diameters (midpoint_interpolation)
  # ncol is the maximum number of column in the original transposed dataframe
  
  # create a new empty dataframe
  Interpolation_output <- matrix(nrow = length(midpoint_interpolation), ncol = length(TS_merged_df))
  
  # the first column in the dataframe corresponds to the new midpoint serie
  Interpolation_output[,1] <- midpoint_interpolation
  
  # the full join function will return a dataframe with three columns
  # the first column correponds to the new midpoint diameter 
  # the second column and the third column are the data from two different joint time series
  # you can only join two terms at once
  # full join all the terms in the list, the first time interval is constant
  # while the second term change
  # in the full join output keep only the third column which is the variable [,3]
  
  for (i in 1:length(spline_interpolation)) {
    Interpolation_output[,i+1] <- (full_join(spline_interpolation[[i]], spline_interpolation[[i]], by = "x"))[,3]
  }
  
  colnames(Interpolation_output) <- colnames(TS_merged_df)
  
  Interpolation_output <- data.frame(Interpolation_output)
  
  # only consider the interpolated value starting/ending from/to valid measurements
  
  # Identify the first and last row that is not NA for each column in the original merged size distributions 
  
  # prepare the data for corrections
  correct_Interpolation_output <- list()
  
  for (i in 2:length(Interpolation_output)) {
    correct_Interpolation_output[[i]] <- data.frame(Interpolation_output[,c(1,i)])
  }
  
  # from the merged serie before interpolation identify the first and last Dg corresponding to valid measurements
  valid_meas <- list()
  
  for (i in 2:length(TS_merged_data_interpolate)) {
    valid_meas[[i]] <- data.frame(TS_merged_data_interpolate[,c(1,i)])
  }
  
  # remove the NA values
  for (i in 2:length(valid_meas)) {
    valid_meas[[i]] <- na.omit(valid_meas[[i]])
  }
  
  # in the interpolation output 
  # for each date time
  # replace the value outside the Dg interval corresponding to valid measurements
  for (i in 2:length(correct_Interpolation_output)) {
    for (j in 1:nrow(correct_Interpolation_output[[i]])) {
      if(correct_Interpolation_output[[i]]$Dg[j] <= valid_meas[[i]]$Dg[1] | correct_Interpolation_output[[i]]$Dg[j] >= valid_meas[[i]]$Dg[nrow(valid_meas[[i]])]){
        correct_Interpolation_output[[i]][j,2] <- NA }
      else {}
    }
  }
  
  # save the output 
  Interpolation_output_NA <- Interpolation_output
  
  for (i in 2:length(Interpolation_output_NA)) {
    Interpolation_output_NA[,i] <- correct_Interpolation_output[[i]][,2]
  }
  
  # reconvert the data as dN/dlogDg
  for (i in 2:length(Interpolation_output_NA)) {
    Interpolation_output_NA[,i] <- exp(Interpolation_output_NA[,i])
  }
  
  # also ensure all the negative values are as NA
  for (i in 2:length(Interpolation_output_NA)) {
    for (j in 1:nrow(Interpolation_output_NA)) {
      if ((Interpolation_output_NA[j,i] < 0) == TRUE | is.na(Interpolation_output_NA[j,i]) == TRUE) {
        Interpolation_output_NA[j,i] = NA} 
      else{}
    }
  }
  
  return(Interpolation_output_NA)
}

test <- interpolate_data(TS_merged_data[[1]][[3]])

identical(test, Interpolation_output_NA)

# drop useless files from the environment
rm(Interpolation_output_NA)

# apply the function to all dfs
interpolated_df <- TS_merged_data

for (i in 1:length(TS_merged_data)) {
  for (j in 1:length(TS_merged_data[[i]])) {
    interpolated_df[[i]][[j]] <- interpolate_data(TS_merged_data[[i]][[j]])
  }
}

######################################################################################################
######################################################################################################
######################################################################################################
# Normalisation
######################################################################################################
######################################################################################################
######################################################################################################

######################################################################################################
# separate GRIMM and SMPS data
######################################################################################################
# to calculate dN total before interpolation we need to split TS_merged_data
# SMPS and GRIMM have different dlogDg
# SMPS dlogDg is 1/64
# GRIMM dlogDg varies depending on the scenario

# get GRIMM merged data before interpolation
GRM_BI <- TS_merged_data

for (i in 1:length(TS_merged_data)) {
  for (j in 1:length(TS_merged_data[[i]])) {
    GRM_BI[[i]][[j]] <- TS_merged_data[[i]][[j]][which(TS_merged_data[[i]][[j]]$Dg >= GRM_Dg_copy[[i]][[j]][1]),]
  }
}

# get SMPS merged data before interpolation
SMPS_BI <- TS_merged_data

for (i in 1:length(TS_merged_data)) {
  for (j in 1:length(TS_merged_data[[i]])) {
    SMPS_BI[[i]][[j]] <- TS_merged_data[[i]][[j]][which(TS_merged_data[[i]][[j]]$Dg < GRM_Dg_copy[[i]][[j]][1]),]
  }
}

# check that the total number of row norm_SMPS - norm_GRM is equal to the number of rom in TS_merged_data
i = 3
j = 89
nrow(TS_merged_data[[i]][[j]]) == nrow(SMPS_BI[[i]][[j]]) + nrow(GRM_BI[[i]][[j]])

######################################################################################################
# build a function to convert dN/dlogDg into dN
######################################################################################################

# first subset the GRM_dlogDg_copy, keep only the dlogDg of the GRIMM diameters used in the merged dataset

for (i in 1:length(GRM_dlogDg_copy)) {
  for (j in 1:length(GRM_dlogDg_copy[[i]])) {
    GRM_dlogDg_copy[[i]][[j]] <- GRM_dlogDg_copy[[i]][[j]][which(GRM_dlogDg_copy[[i]][[j]]$Dg >= GRM_Dg_copy[[i]][[j]][1]),]
  }
}

######################################################################################################
# Example for GRIMM
# for one scenario get dN/dlogDg concentration and the corresponding dlogDg values
test_GRIMM <- GRM_BI[[1]][[7]]
dlogDg_GRIMM <- GRM_dlogDg_copy[[1]][[7]][,2]# dlogDg is at the second column 

# get dN/dlogDg concentration and multiply each bin for the corresponding dlogDg value
# for each row (Dg) apply a different dlogDg
test_GRIMM_dN <- test_GRIMM
for (i in 2:length(test_GRIMM)) {
  for (k in 1:nrow(test_GRIMM)) {
    test_GRIMM_dN[k,i] <- test_GRIMM[k,i]*dlogDg_GRIMM[k]
  }
}

# define the function for GRM
convert_dNdlogDg_into_dN_GRM <- function(df_dNdlogDg,dlogDg) {
  df_dN <- df_dNdlogDg
  for (i in 2:length(df_dNdlogDg)) {
    for (k in 1:nrow(df_dNdlogDg)) {
      df_dN[k,i] <- df_dNdlogDg[k,i]*dlogDg[k]
    }
  }
  return(df_dN)
}

# check the results
test <- convert_dNdlogDg_into_dN_GRM(GRM_BI[[1]][[7]],GRM_dlogDg_copy[[1]][[7]][,2])
identical(test_GRIMM_dN, test)

######################################################################################################
# Example for SMPS
# for one scenario get dN/dlogDg concentration and the corresponding dlogDg values
test_SMPS <- SMPS_BI[[1]][[7]]

# get dN/dlogDg concentration and multiply each bin for the corresponding dlogDg value
# for each row (Dg) apply a different dlogDg
test_SMPS_dN <- test_SMPS
for (i in 2:length(test_SMPS)) {
  test_SMPS_dN[,i] <- test_SMPS[,i]/64
}

# define the function for SMPS
convert_dNdlogDg_into_dN <- function(df_dNdlogDg) {
  df_dN <- df_dNdlogDg
  for (i in 2:length(df_dNdlogDg)) {
    df_dN[,i] <- df_dNdlogDg[,i]/64
  }
  return(df_dN)
}

# check the results
test <- convert_dNdlogDg_into_dN(SMPS_BI[[1]][[7]])
identical(test_SMPS_dN, test)

######################################################################################################
# the same function for SMPS can be used to convert the interpolated data in dN
# for the interpolated data dlogDg is 1/64
SD_AI <- interpolated_df
######################################################################################################

# apply the functions to all dataframes
for (i in 1:length(GRM_BI)) {
  for (j in 1:length(GRM_BI[[i]])) {
    GRM_BI[[i]][[j]] <- convert_dNdlogDg_into_dN_GRM(GRM_BI[[i]][[j]], GRM_dlogDg_copy[[i]][[j]][,2])
    SMPS_BI[[i]][[j]] <- convert_dNdlogDg_into_dN(SMPS_BI[[i]][[j]])
    SD_AI[[i]][[j]] <- convert_dNdlogDg_into_dN(SD_AI[[i]][[j]])
    
  }
}

# check the results
identical(test_GRIMM_dN, GRM_BI[[1]][[7]])
identical(test_SMPS_dN, SMPS_BI[[1]][[7]])

######################################################################################################
# calculate the total dN
######################################################################################################

# build a function to apply to single dfs to calculate the sum of values in each column
# and save the output as a numeric serie
calculate_dN_total <- function(df) {
  # calculate the sum of the terms in each column
  dN_sum <- list()
  for (i in 2:length(df)) {
    dN_sum[[i]] <- sum(df[,i], na.rm = T)
  }
  # drop the NULL terms in the list
  dN_sum <- dN_sum[lengths(dN_sum) != 0]
  
  # convert the list into a numeric serie
  num_serie <- convert_list_df(dN_sum)
  num_serie <- unlist(num_serie)
  
  return(num_serie)
}

# build a large list of dfs to save the output
dN_total <- list()
for (i in 1:length(TS_merged_data)) {
  for (j in 1:length(TS_merged_data[[i]])) {
    dN_total[[i]] <- list()
  }
}

for (i in 1:length(TS_merged_data)) {
  for (j in 1:length(TS_merged_data[[i]])) {
    # build a matrix to save dN total calculation
    dN_total[[i]][[j]] <- matrix(ncol = 4, nrow = nrow(GRM_data_TA))
    dN_total[[i]][[j]] <- data.frame(dN_total[[i]][[j]])
    colnames(dN_total[[i]][[j]]) <- c("date", "dN_BI_SMPS", "dN_BI_GRM", "dN_AI")
  }
}

# calculate the total dN, before interpolation, when replacing the outliers, and after interpolation
for (i in 1:length(dN_total)) {
  for (j in 1:length(dN_total[[i]])) {
    dN_total[[i]][[j]][,1] <- GRM_data_TA$date
    dN_total[[i]][[j]][,2] <- calculate_dN_total(SMPS_BI[[i]][[j]])
    dN_total[[i]][[j]][,3] <- calculate_dN_total(GRM_BI[[i]][[j]])
    dN_total[[i]][[j]][,4] <- calculate_dN_total(SD_AI[[i]][[j]])
  }
}

# to normalise the data use the ratio between dN total before interpolation
# and dN total calculated after interpolating the data

for (i in 1:length(dN_total)) {
  for (j in 1:length(dN_total[[i]])) {
    dN_total[[i]][[j]] <- add_column(dN_total[[i]][[j]],
                                     Norm_ratio = (dN_total[[i]][[j]][,2] + dN_total[[i]][[j]][,3])/dN_total[[i]][[j]][,4],
                                     .after = "dN_AI") 
  } 
}

# build a function to normalise the data
normalise_data <- function(interp_df, dN){
  
  data_norm <- interp_df[,-1]
  for (i in 1:length(data_norm)) {
    for (j in 1:nrow(data_norm)) {
      data_norm[j,i] <- data_norm[j,i]*dN$Norm_ratio[i]
    }
  }
  
  data_norm <- cbind(interp_df[,1], data_norm)
  return(data_norm)
}

# test the function to normalise the data
i = 1
j = 1
test <- normalise_data(interpolated_df[[i]][[j]], dN_total[[i]][[j]])

# apply the function to normalise the data
normalised_df <- interpolated_df

for (i in 1:length(interpolated_df)) {
  for (j in 1:length(interpolated_df[[i]])) {
    normalised_df[[i]][[j]] <- normalise_data(interpolated_df[[i]][[j]], dN_total[[i]][[j]])
  }
}

# test dN total after normalisation is consistent with dN before interpolation corrected for outliers
i = 4
j = 147

test <- convert_dNdlogDg_into_dN(normalised_df[[i]][[j]])
test_dN <- calculate_dN_total(test)
test_dN - (dN_total[[i]][[j]]$dN_BI_SMPS + dN_total[[i]][[j]]$dN_BI_GRM)

######################################################################################################
# plot the data - normalisation
######################################################################################################

par(mfrow= c(2,2))

i = 1
j = 3


plot(interpolated_df[[i]][[j]][,1], interpolated_df[[i]][[j]][,2],
     col = "green",
     xaxt = "n",
     yaxt = "n",
     xlab = expression(paste("Diameter (",mu,"m)")),
     ylab = expression(paste("dN/dlogDg (",cm^{-3}*")")),
     xlim = c(0.01, 35),
     ylim = c(0.0005, 1500),
     main= "dN/dlogDg at 11:07",
     log = "xy", pch = 20)
axis(1, at = c(0.1,1,10), labels = c(0.1,1,10))
axis(2, at = c(0.1,1,10,100,1000), labels = c(0.1,1,10,100,1000))
points(TS_merged_data[[i]][[j]][,1], TS_merged_data[[i]][[j]][,2],
       col = "black", pch = 20)
points(normalised_df[[i]][[j]][,1], normalised_df[[i]][[j]][,2],
       col = "hotpink", pch = 20)
legend("bottomleft", inset =.02, box.lty = 0,
       legend = c("Raw data", "After Interpolation", "Normalised values"),
       pch = c(20,20,20),
       col = c("black", "green", "hotpink" ) ,cex=0.75)
mtext("n = 1.59, k = 0, X = 1", cex = 0.75, side=3,line=0, col = "red")


plot(interpolated_df[[i]][[j]][,1], interpolated_df[[i]][[j]][,4],
     col = "green",
     xaxt = "n",
     yaxt = "n",
     xlab = expression(paste("Diameter (",mu,"m)")),
     ylab = expression(paste("dN/dlogDg (",cm^{-3}*")")),
     xlim = c(0.01, 35),
     ylim = c(0.0005, 1500),
     main= "dN/dlogDg at 11:31",
     log = "xy", pch = 20)
axis(1, at = c(0.1,1,10), labels = c(0.1,1,10))
axis(2, at = c(0.1,1,10,100,1000), labels = c(0.1,1,10,100,1000))
points(TS_merged_data[[i]][[j]][,1], TS_merged_data[[i]][[j]][,4],
       col = "black", pch = 20)
points(normalised_df[[i]][[j]][,1], normalised_df[[i]][[j]][,4],
       col = "hotpink", pch = 20)
legend("bottomleft", inset =.02, box.lty = 0,
       legend = c("Raw data", "After Interpolation", "Normalised values"),
       pch = c(20,20,20),
       col = c("black", "green", "hotpink" ) ,cex=0.75)
mtext("n = 1.59, k = 0, X = 1", cex = 0.75, side=3,line=0, col = "red")


plot(interpolated_df[[i]][[j]][,1], interpolated_df[[i]][[j]][,6],
     col = "green",
     xaxt = "n",
     yaxt = "n",
     xlab = expression(paste("Diameter (",mu,"m)")),
     ylab = expression(paste("dN/dlogDg (",cm^{-3}*")")),
     xlim = c(0.01, 35),
     ylim = c(0.0005, 1500),
     main= "dN/dlogDg at 11:55",
     log = "xy", pch = 20)
axis(1, at = c(0.1,1,10), labels = c(0.1,1,10))
axis(2, at = c(0.1,1,10,100,1000), labels = c(0.1,1,10,100,1000))
points(TS_merged_data[[i]][[j]][,1], TS_merged_data[[i]][[j]][,6],
       col = "black", pch = 20)
points(normalised_df[[i]][[j]][,1], normalised_df[[i]][[j]][,6],
       col = "hotpink", pch = 20)
legend("bottomleft", inset =.02, box.lty = 0,
       legend = c("Raw data", "After Interpolation", "Normalised values"),
       pch = c(20,20,20),
       col = c("black", "green", "hotpink" ) ,cex=0.75)
mtext("n = 1.59, k = 0, X = 1", cex = 0.75, side=3,line=0, col = "red")


plot(interpolated_df[[i]][[j]][,1], interpolated_df[[i]][[j]][,9],
     col = "green",
     xaxt = "n",
     yaxt = "n",
     xlab = expression(paste("Diameter (",mu,"m)")),
     ylab = expression(paste("dN/dlogDg (",cm^{-3}*")")),
     xlim = c(0.01, 35),
     ylim = c(0.0005, 1500),
     main= "dN/dlogDg at 12:31",
     log = "xy", pch = 20)
axis(1, at = c(0.1,1,10), labels = c(0.1,1,10))
axis(2, at = c(0.1,1,10,100,1000), labels = c(0.1,1,10,100,1000))
points(TS_merged_data[[i]][[j]][,1], TS_merged_data[[i]][[j]][,9],
       col = "black", pch = 20)
points(normalised_df[[i]][[j]][,1], normalised_df[[i]][[j]][,9],
       col = "hotpink", pch = 20)
legend("bottomleft", inset =.02, box.lty = 0,
       legend = c("Raw data", "After Interpolation", "Normalised values"),
       pch = c(20,20,20),
       col = c("black", "green", "hotpink" ) ,cex=0.75)
mtext("n = 1.59, k = 0, X = 1", cex = 0.75, side=3,line=0, col = "red")


######################################################################################################
######################################################################################################
######################################################################################################
# Corrrection for particle loss
######################################################################################################
######################################################################################################
######################################################################################################

# remove useless object from the environment
rm(test, test_GRIMM, test_GRIMM_dN, test_SMPS, test_SMPS_dN, SMPS_BI, GRM_BI, SD_AI)

# Import and interpolate the particle loss table
setwd("Z:/Clarissa/Data_Optical_Calculation/Loss_Calculation/PLC_all/Loss_data/")

files_PLC <- list.files(pattern = ".csv")
PLC <- list()

for (i in 1:length(files_PLC)){
  PLC[[i]] <- read.csv(files_PLC[[i]], header = T)
}

# replace te 0 values (100% particle loss) with NA 
for (i in 1:length(PLC)) {
  PLC[[i]][PLC[[i]][] == 0] <- NA
}

# Keep only the shape factors from 1.6 to 2
test <- PLC[c(1,7:11)]

identical(test[[3]], PLC[[8]])

# now save the output as PLC
PLC <- test

#####################################################################################################
# prepare PLC data for calculation
######################################################################################################

# PLC is calculated based on the different X scenario 
# Thus, each PLC files is copied # times = # n-k scenarios

Loss_table <- TS_merged_data

for (i in 1:length(TS_merged_data)) {
  for (j in 1:length(TS_merged_data[[i]])) {
    Loss_table[[i]][[j]] <- PLC[[i]]
  }
}

# check the output 
identical(Loss_table[[1]][[1]],Loss_table[[1]][[147]])
identical(Loss_table[[2]][[1]],Loss_table[[2]][[147]])
identical(Loss_table[[2]][[1]],Loss_table[[1]][[1]])

#####################################################################################################
# PLC interpolation - test
######################################################################################################
# PLC data must be interpolated at the midpoint Dg series used for interpolating the merged SD

# interpolate the original data, method = linear
Loss <- Loss_table[[1]][[1]]
Norm_data <- normalised_df[[1]][[1]]

Loss_interpol <- list()
for (i in 2:length(Loss)) {
  Loss_interpol[[i]] <- approx(Loss[,1], Loss[,i], Norm_data[,1], method="linear", ties = mean)
}

# save the output
Loss_data  <- list()
# save the output
for (i in 2:length(Loss_interpol)) {
  Loss_data[[i]] <- Loss_interpol[[i]][[2]]
}

# the first term in the list is the midpoint Dg
Loss_data[[1]] <- Norm_data[,1]

# bind together the element in the list
Loss_table_interpol <- do.call(cbind, Loss_data)

colnames(Loss_table_interpol) <- colnames(Loss)

Loss_table_interpol <- data.frame(Loss_table_interpol)

######################################################################################################
# plot the data - Loss table
######################################################################################################

par(mfrow= c(1,1))

plot(Loss$dp,Loss$Total.Loss.Aethalometer, 
     col = "black",
     xaxt = "n",
     yaxt = "n",
     xlim = c(0.01,40),
     ylim = c(0.1,100),
     xlab = expression(paste("Diameter (",mu,"m)")),
     ylab = "Loss (%)",
     main = "Summary - Loss %",
     log = "xy", pch = 20)
axis(1, at = c(0.1,1,10), labels = c(0.1,1,10))
axis(2, at = c(0.1,1,10,100), labels = c(0.1,1,10,100))
points(Loss$dp,Loss$Total.Loss.Caps.Blue, col = "blue", pch = 20)
points(Loss$dp,Loss$Total.Loss.Caps.Red, col = "red", pch = 20)
points(Loss$dp,Loss$Total.Loss.Nephelometer, col = "yellow", pch = 20)
points(Loss$dp,Loss$Total.Loss.Grimm, col = "purple", pch = 20)
points(Loss$dp,Loss$Total.Loss.SkyGrimm, col = "green", pch = 20)
points(Loss$dp,Loss$Total.Loss.Filters, col = "pink", pch = 20)
points(Loss$dp,Loss$Total.Loss.SMPS, col = "orange", pch = 20)

legend("bottomright", inset =.02, box.lty = 0,
       legend = c("Aethalometer", "Caps Blue", "Caps Red", "Nephelometer", "Grimm", "SkyGrimm", "Filters", "SMPS"),
       pch = c(20,20,20,20,20,20,20),
       col = c("black", "blue", "red", "yellow", "purple", "green", "pink", "orange"), cex=1)


plot(Loss_table_interpol$dp,Loss_table_interpol$Total.Loss.Aethalometer, 
     col = "black",
     xaxt = "n",
     yaxt = "n",
     xlim = c(0.01,40),
     ylim = c(0.1,100),
     xlab = expression(paste("Diameter (",mu,"m)")),
     ylab = "Loss (%)",
     main = "Summary - Loss % (Interpolated values)",
     log = "xy", pch = 20)
axis(1, at = c(0.1,1,10), labels = c(0.1,1,10))
axis(2, at = c(0.1,1,10,100), labels = c(0.1,1,10,100))
points(Loss_table_interpol$dp,Loss_table_interpol$Total.Loss.Caps.Blue, col = "blue", pch = 20)
points(Loss_table_interpol$dp,Loss_table_interpol$Total.Loss.Caps.Red, col = "red", pch = 20)
points(Loss_table_interpol$dp,Loss_table_interpol$Total.Loss.Nephelometer, col = "yellow", pch = 20)
points(Loss_table_interpol$dp,Loss_table_interpol$Total.Loss.Grimm, col = "purple", pch = 20)
points(Loss_table_interpol$dp,Loss_table_interpol$Total.Loss.SkyGrimm, col = "green", pch = 20)
points(Loss_table_interpol$dp,Loss_table_interpol$Total.Loss.Filters, col = "pink", pch = 20)
points(Loss_table_interpol$dp,Loss_table_interpol$Total.Loss.SMPS, col = "orange", pch = 20)

legend("bottomright", inset =.02, box.lty = 0,
       legend = c("Aethalometer", "Caps Blue", "Caps Red", "Nephelometer", "Grimm", "SkyGrimm", "Filters", "SMPS"),
       pch = c(20,20,20,20,20,20,20),
       col = c("black", "blue", "red", "yellow", "purple", "green", "pink", "orange"), cex=1)

#####################################################################################################
# build a function to interpolate all df
######################################################################################################

Loss_interpolation_function <- function(Loss, Norm_data) {
  
  Loss_interpol <- list()
  for (i in 2:length(Loss)) {
    Loss_interpol[[i]] <- approx(Loss[,1], Loss[,i], Norm_data[,1], method="linear", ties = mean)
  }
  
  # save the output
  Loss_data  <- list()
  # save the output
  for (i in 2:length(Loss_interpol)) {
    Loss_data[[i]] <- Loss_interpol[[i]][[2]]
  }
  
  # the first term in the list is the midpoint Dg
  Loss_data[[1]] <- Norm_data[,1]
  
  # bind together the element in the list
  Loss_table_interpol <- do.call(cbind, Loss_data)
  
  colnames(Loss_table_interpol) <- colnames(Loss)
  
  Loss_table_interpol <- data.frame(Loss_table_interpol)
  
  return(Loss_table_interpol)
}

# apply the function to all df
Loss_df <- normalised_df

for (i in 1:length(normalised_df)) {
  for (j in 1:length(normalised_df[[i]]) ) {
    Loss_df[[i]][[j]] <- Loss_interpolation_function(Loss_table[[i]][[j]], normalised_df[[i]][[j]])    
  }
}

# check the output 
identical(Loss_table_interpol,Loss_df[[1]][[1]])

######################################################################################################
# merge Loss_table_interpolated and the normalization_output
######################################################################################################

# remove useless object from the environment
rm(Loss, Loss_data, Loss_interpol, Norm_data, Loss_table_interpol)

# prepare the Loss column, as you will need to divide the dN values by 1 - Loss/100
for (i in 1:length(Loss_df)) {
  for (j in 1:length(Loss_df[[i]])) {
    for (k in 2:length(Loss_df[[i]][[j]])){
      Loss_df[[i]][[j]][,k] <- 1 - Loss_df[[i]][[j]][,k]/100
    }
  } 
}

# bind loss and data normalised dfs
for (i in 1:length(Loss_df)) {
  for (j in 1:length(Loss_df[[i]])) {
    Loss_df[[i]][[j]] <- cbind(Loss_df[[i]][[j]], normalised_df[[i]][[j]][,-1]) #remove the diameter column since you have two after binding the datasets
    colnames(Loss_df[[i]][[j]])[1] <- "Dg"
    Loss_df[[i]][[j]] <- data.frame(Loss_df[[i]][[j]])
  } 
}

######################################################################################################
######################################################################################################
# Calculate the real size distribution in CESAM
######################################################################################################
######################################################################################################

# Calculate the real size distribution in CESAM after correcting for GRM and SMPS particle loss
# split the table in two as part of the data must be corrected for the SMPS particle Loss
# and part for the GRM particle loss

# get the SMPS Dg values before merging the data
SMPS_Dg <- SMPS_Dg_copy

# get the SMPS data 
SD_CESAM_SMPS <- Loss_df

for (i in 1:length(SD_CESAM_SMPS)) {
  for (j in 1:length(SD_CESAM_SMPS[[i]])) {
    SD_CESAM_SMPS[[i]][[j]] <- subset(SD_CESAM_SMPS[[i]][[j]], SD_CESAM_SMPS[[i]][[j]][,1] <= max(SMPS_Dg[[i]][[j]]))
  }
}

# get the GRIMM data 
SD_CESAM_GRM <- Loss_df

for (i in 1:length(SD_CESAM_GRM)) {
  for (j in 1:length(SD_CESAM_GRM[[i]])) {
    SD_CESAM_GRM[[i]][[j]] <- subset(SD_CESAM_GRM[[i]][[j]], SD_CESAM_GRM[[i]][[j]][,1] > max(SMPS_Dg[[i]][[j]]))
  }
}

# correct the SMPS df for particle loss
SD_CESAM_SMPS_corr <- SD_CESAM_SMPS
colnumber <- which(colnames(SD_CESAM_SMPS[[1]][[1]]) == "Total.Loss.SMPS")
for (i in 1:length(SD_CESAM_SMPS_corr)) {
  for (j in 1:length(SD_CESAM_SMPS_corr[[i]])) {
    for (k in 10:length(SD_CESAM_SMPS_corr[[i]][[j]])) {
      SD_CESAM_SMPS_corr[[i]][[j]][,k] <- SD_CESAM_SMPS_corr[[i]][[j]][,k]/SD_CESAM_SMPS_corr[[i]][[j]][,colnumber]
    } 
  }
}

# correct the GRM df for particle loss
SD_CESAM_GRM_corr <- SD_CESAM_GRM
colnumber <- which(colnames(SD_CESAM_GRM[[1]][[1]]) == "Total.Loss.Grimm")
for (i in 1:length(SD_CESAM_GRM_corr)) {
  for (j in 1:length(SD_CESAM_GRM_corr[[i]])) {
    for (k in 10:length(SD_CESAM_GRM_corr[[i]][[j]])) {
      SD_CESAM_GRM_corr[[i]][[j]][,k] <- SD_CESAM_GRM_corr[[i]][[j]][,k]/SD_CESAM_GRM_corr[[i]][[j]][,colnumber]
    } 
  }
}

# now re-combine the two datasets
SD_CESAM <- SD_CESAM_SMPS_corr

for (i in 1:length(SD_CESAM_SMPS_corr)) {
  for (j in 1:length(SD_CESAM_SMPS_corr[[i]])) {
    SD_CESAM[[i]][[j]] <- rbind(SD_CESAM_SMPS_corr[[i]][[j]], SD_CESAM_GRM_corr[[i]][[j]])
    SD_CESAM[[i]][[j]] <- data.frame(SD_CESAM[[i]][[j]])
  }
}

# check the output
i = 3
j = 4
ncol(SD_CESAM[[i]][[j]]) == ncol(Loss_df[[i]][[j]])
nrow(SD_CESAM[[i]][[j]]) == nrow(Loss_df[[i]][[j]])
SD_CESAM[[i]][[j]]$Dg == Loss_df[[i]][[j]]$Dg

######################################################################################################
# plot the data - CESAM SD
######################################################################################################

par(mfrow= c(2,2))

i = 1
j = 3

plot(normalised_df[[i]][[j]][,1], normalised_df[[i]][[j]][,2],
     col = "hotpink",
     xaxt = "n",
     yaxt = "n",
     xlab = expression(paste("Diameter (",mu,"m)")),
     ylab = expression(paste("dN/dlogDg (",cm^{-3}*")")),
     xlim = c(0.01, 35),
     ylim = c(0.0005, 1500),
     main= "dN/dlogDg at 11:07",
     log = "xy", pch = 20)
axis(1, at = c(0.1,1,10), labels = c(0.1,1,10))
axis(2, at = c(0.1,1,10,100,1000), labels = c(0.1,1,10,100,1000))
points(TS_merged_data[[i]][[j]][,1], TS_merged_data[[i]][[j]][,2],
       col = "black", pch = 20)
points(interpolated_df[[i]][[j]][,1], interpolated_df[[i]][[j]][,2],
       col = "green", pch = 20)
points(normalised_df[[i]][[j]][,1], normalised_df[[i]][[j]][,2],
       col = "hotpink", pch = 20)# repeat so it can overlap the interpolated values
points(SD_CESAM[[i]][[j]][,1], SD_CESAM[[i]][[j]][,10],
       col = "blue", pch = 20)
legend("bottomleft", inset =.02, box.lty = 0,
       legend = c("Raw data", "After Interpolation", "Normalised values", "After GRIMM-SMPS loss correction"),
       pch = c(20,20,20,20),
       col = c("black", "green", "hotpink", "blue") ,cex=0.6)
mtext("n = 1.59, k = 0, X = 1", cex = 0.75, side=3,line=0, col = "red")


plot(normalised_df[[i]][[j]][,1], normalised_df[[i]][[j]][,4],
     col = "hotpink",
     xaxt = "n",
     yaxt = "n",
     xlab = expression(paste("Diameter (",mu,"m)")),
     ylab = expression(paste("dN/dlogDg (",cm^{-3}*")")),
     xlim = c(0.01, 35),
     ylim = c(0.0005, 1500),
     main= "dN/dlogDg at 11:31",
     log = "xy", pch = 20)
axis(1, at = c(0.1,1,10), labels = c(0.1,1,10))
axis(2, at = c(0.1,1,10,100,1000), labels = c(0.1,1,10,100,1000))
points(TS_merged_data[[i]][[j]][,1], TS_merged_data[[i]][[j]][,4],
       col = "black", pch = 20)
points(interpolated_df[[i]][[j]][,1], interpolated_df[[i]][[j]][,4],
       col = "green", pch = 20)
points(normalised_df[[i]][[j]][,1], normalised_df[[i]][[j]][,4],
       col = "hotpink", pch = 20)# repeat so it can overlap the interpolated values
points(SD_CESAM[[i]][[j]][,1], SD_CESAM[[i]][[j]][,12],
       col = "blue", pch = 20)
legend("bottomleft", inset =.02, box.lty = 0,
       legend = c("Raw data", "After Interpolation", "Normalised values", "After GRIMM-SMPS loss correction"),
       pch = c(20,20,20,20),
       col = c("black", "green", "hotpink", "blue") ,cex=0.6)
mtext("n = 1.59, k = 0, X = 1", cex = 0.75, side=3,line=0, col = "red")


plot(normalised_df[[i]][[j]][,1], normalised_df[[i]][[j]][,6],
     col = "hotpink",
     xaxt = "n",
     yaxt = "n",
     xlab = expression(paste("Diameter (",mu,"m)")),
     ylab = expression(paste("dN/dlogDg (",cm^{-3}*")")),
     xlim = c(0.01, 35),
     ylim = c(0.0005, 1500),
     main= "dN/dlogDg at 11:55",
     log = "xy", pch = 20)
axis(1, at = c(0.1,1,10), labels = c(0.1,1,10))
axis(2, at = c(0.1,1,10,100,1000), labels = c(0.1,1,10,100,1000))
points(TS_merged_data[[i]][[j]][,1], TS_merged_data[[i]][[j]][,6],
       col = "black", pch = 20)
points(interpolated_df[[i]][[j]][,1], interpolated_df[[i]][[j]][,6],
       col = "green", pch = 20)
points(normalised_df[[i]][[j]][,1], normalised_df[[i]][[j]][,6],
       col = "hotpink", pch = 20)# repeat so it can overlap the interpolated values
points(SD_CESAM[[i]][[j]][,1], SD_CESAM[[i]][[j]][,14],
       col = "blue", pch = 20)
legend("bottomleft", inset =.02, box.lty = 0,
       legend = c("Raw data", "After Interpolation", "Normalised values", "After GRIMM-SMPS loss correction"),
       pch = c(20,20,20,20),
       col = c("black", "green", "hotpink", "blue") ,cex=0.6)
mtext("n = 1.59, k = 0, X = 1", cex = 0.75, side=3,line=0, col = "red")

plot(normalised_df[[i]][[j]][,1], normalised_df[[i]][[j]][,9],
     col = "hotpink",
     xaxt = "n",
     yaxt = "n",
     xlab = expression(paste("Diameter (",mu,"m)")),
     ylab = expression(paste("dN/dlogDg (",cm^{-3}*")")),
     xlim = c(0.01, 35),
     ylim = c(0.0005, 1500),
     main= "dN/dlogDg at 12:31",
     log = "xy", pch = 20)
axis(1, at = c(0.1,1,10), labels = c(0.1,1,10))
axis(2, at = c(0.1,1,10,100,1000), labels = c(0.1,1,10,100,1000))
points(TS_merged_data[[i]][[j]][,1], TS_merged_data[[i]][[j]][,9],
       col = "black", pch = 20)
points(interpolated_df[[i]][[j]][,1], interpolated_df[[i]][[j]][,9],
       col = "green", pch = 20)
points(normalised_df[[i]][[j]][,1], normalised_df[[i]][[j]][,9],
       col = "hotpink", pch = 20)
points(SD_CESAM[[i]][[j]][,1], SD_CESAM[[i]][[j]][,17],
       col = "blue", pch = 20)
legend("bottomleft", inset =.02, box.lty = 0,
       legend = c("Raw data", "After Interpolation", "Normalised values", "After GRIMM-SMPS loss correction"),
       pch = c(20,20,20,20),
       col = c("black", "green", "hotpink", "blue") ,cex=0.6)
mtext("n = 1.59, k = 0, X = 1", cex = 0.75, side=3,line=0, col = "red")


######################################################################################################
######################################################################################################
# Aethalometer size distribution
######################################################################################################
######################################################################################################

SD_Aethalometer <- SD_CESAM
colnumber <- which(colnames(SD_CESAM[[1]][[1]]) == "Total.Loss.Aethalometer")

for (i in 1:length(SD_CESAM)) {
  for (j in 1:length(SD_CESAM[[i]])) {
    for (k in 10:length(SD_CESAM[[i]][[j]])) {
      SD_Aethalometer[[i]][[j]][,k] <- SD_CESAM[[i]][[j]][,k]*SD_CESAM[[i]][[j]][,colnumber]
    }
  }
}

# check the output
i = 4
j = 18

identical(SD_Aethalometer[[i]][[j]][,1], SD_CESAM[[i]][[j]][,1])

######################################################################################################
# plot the data - Aethalometer SD
######################################################################################################

par(mfrow= c(2,2))

i = 1
j = 3

plot(interpolated_df[[i]][[j]][,1], interpolated_df[[i]][[j]][,2],
     col = "green",
     xaxt = "n",
     yaxt = "n",
     xlab = expression(paste("Diameter (",mu,"m)")),
     ylab = expression(paste("dN/dlogDg (",cm^{-3}*")")),
     xlim = c(0.01, 35),
     ylim = c(0.0005, 1500),
     main= "dN/dlogDg at 11:07",
     log = "xy", pch = 20)
axis(1, at = c(0.1,1,10), labels = c(0.1,1,10))
axis(2, at = c(0.1,1,10,100,1000), labels = c(0.1,1,10,100,1000))
points(TS_merged_data[[i]][[j]][,1], TS_merged_data[[i]][[j]][,2],
       col = "black", pch = 20)
points(normalised_df[[i]][[j]][,1], normalised_df[[i]][[j]][,2],
       col = "hotpink", pch = 20)
points(SD_CESAM[[i]][[j]][,1], SD_CESAM[[i]][[j]][,10],
       col = "blue", pch = 20)
points(SD_Aethalometer[[i]][[j]][,1], SD_Aethalometer[[i]][[j]][,10],
       col = "gold", pch = 20)
legend("bottomleft", inset =.02, box.lty = 0,
       legend = c("Raw data", "After Interpolation", "Normalised values", "After GRIMM-SMPS loss correction", "Aethalometer SD"),
       pch = c(20,20,20,20,20),
       col = c("black", "green", "hotpink", "blue", "gold") ,cex=0.6)
mtext("n = 1.59, k = 0, X = 1", cex = 0.75, side=3,line=0, col = "red")


plot(interpolated_df[[i]][[j]][,1], interpolated_df[[i]][[j]][,4],
     col = "green",
     xaxt = "n",
     yaxt = "n",
     xlab = expression(paste("Diameter (",mu,"m)")),
     ylab = expression(paste("dN/dlogDg (",cm^{-3}*")")),
     xlim = c(0.01, 35),
     ylim = c(0.0005, 1500),
     main= "dN/dlogDg at 11:31",
     log = "xy", pch = 20)
axis(1, at = c(0.1,1,10), labels = c(0.1,1,10))
axis(2, at = c(0.1,1,10,100,1000), labels = c(0.1,1,10,100,1000))
points(TS_merged_data[[i]][[j]][,1], TS_merged_data[[i]][[j]][,4],
       col = "black", pch = 20)
points(normalised_df[[i]][[j]][,1], normalised_df[[i]][[j]][,4],
       col = "hotpink", pch = 20)
points(SD_CESAM[[i]][[j]][,1], SD_CESAM[[i]][[j]][,12],
       col = "blue", pch = 20)
points(SD_Aethalometer[[i]][[j]][,1], SD_Aethalometer[[i]][[j]][,12],
       col = "gold", pch = 20)
legend("bottomleft", inset =.02, box.lty = 0,
       legend = c("Raw data", "After Interpolation", "Normalised values", "After GRIMM-SMPS loss correction", "Aethalometer SD"),
       pch = c(20,20,20,20,20),
       col = c("black", "green", "hotpink", "blue", "gold") ,cex=0.6)
mtext("n = 1.59, k = 0, X = 1", cex = 0.75, side=3,line=0, col = "red")


plot(interpolated_df[[i]][[j]][,1], interpolated_df[[i]][[j]][,6],
     col = "green",
     xaxt = "n",
     yaxt = "n",
     xlab = expression(paste("Diameter (",mu,"m)")),
     ylab = expression(paste("dN/dlogDg (",cm^{-3}*")")),
     xlim = c(0.01, 35),
     ylim = c(0.0005, 1500),
     main= "dN/dlogDg at 11:55",
     log = "xy", pch = 20)
axis(1, at = c(0.1,1,10), labels = c(0.1,1,10))
axis(2, at = c(0.1,1,10,100,1000), labels = c(0.1,1,10,100,1000))
points(TS_merged_data[[i]][[j]][,1], TS_merged_data[[i]][[j]][,6],
       col = "black", pch = 20)
points(normalised_df[[i]][[j]][,1], normalised_df[[i]][[j]][,6],
       col = "hotpink", pch = 20)
points(SD_CESAM[[i]][[j]][,1], SD_CESAM[[i]][[j]][,14],
       col = "blue", pch = 20)
points(SD_Aethalometer[[i]][[j]][,1], SD_Aethalometer[[i]][[j]][,14],
       col = "gold", pch = 20)
legend("bottomleft", inset =.02, box.lty = 0,
       legend = c("Raw data", "After Interpolation", "Normalised values", "After GRIMM-SMPS loss correction", "Aethalometer SD"),
       pch = c(20,20,20,20,20),
       col = c("black", "green", "hotpink", "blue", "gold") ,cex=0.6)
mtext("n = 1.59, k = 0, X = 1", cex = 0.75, side=3,line=0, col = "red")


plot(interpolated_df[[i]][[j]][,1], interpolated_df[[i]][[j]][,9],
     col = "green",
     xaxt = "n",
     yaxt = "n",
     xlab = expression(paste("Diameter (",mu,"m)")),
     ylab = expression(paste("dN/dlogDg (",cm^{-3}*")")),
     xlim = c(0.01, 35),
     ylim = c(0.0005, 1500),
     main= "dN/dlogDg at 12:31",
     log = "xy", pch = 20)
axis(1, at = c(0.1,1,10), labels = c(0.1,1,10))
axis(2, at = c(0.1,1,10,100,1000), labels = c(0.1,1,10,100,1000))
points(TS_merged_data[[i]][[j]][,1], TS_merged_data[[i]][[j]][,9],
       col = "black", pch = 20)
points(normalised_df[[i]][[j]][,1], normalised_df[[i]][[j]][,9],
       col = "hotpink", pch = 20)
points(SD_CESAM[[i]][[j]][,1], SD_CESAM[[i]][[j]][,17],
       col = "blue", pch = 20)
points(SD_Aethalometer[[i]][[j]][,1], SD_Aethalometer[[i]][[j]][,17],
       col = "gold", pch = 20)
legend("bottomleft", inset =.02, box.lty = 0,
       legend = c("Raw data", "After Interpolation", "Normalised values", "After GRIMM-SMPS loss correction", "Aethalometer SD"),
       pch = c(20,20,20,20,20),
       col = c("black", "green", "hotpink", "blue", "gold") ,cex=0.6)
mtext("n = 1.59, k = 0, X = 1", cex = 0.75, side=3,line=0, col = "red")


######################################################################################################
######################################################################################################
# Nephelometer size distribution
######################################################################################################
######################################################################################################

SD_Nephelometer <- SD_CESAM
colnumber <- which(colnames(SD_CESAM[[1]][[1]]) == "Total.Loss.Nephelometer")

for (i in 1:length(SD_CESAM)) {
  for (j in 1:length(SD_CESAM[[i]])) {
    for (k in 10:length(SD_CESAM[[i]][[j]])) {
      SD_Nephelometer[[i]][[j]][,k] <- SD_CESAM[[i]][[j]][,k]*SD_CESAM[[i]][[j]][,colnumber]
    }
  }
}


# check the output
i = 2
j = 147

identical(SD_Nephelometer[[i]][[j]][,1], SD_CESAM[[i]][[j]][,1])

######################################################################################################
# plot the data - Nephelometer SD
######################################################################################################

par(mfrow= c(2,2))

i = 1
j = 3

plot(interpolated_df[[i]][[j]][,1], interpolated_df[[i]][[j]][,2],
     col = "green",
     xaxt = "n",
     yaxt = "n",
     xlab = expression(paste("Diameter (",mu,"m)")),
     ylab = expression(paste("dN/dlogDg (",cm^{-3}*")")),
     xlim = c(0.01, 35),
     ylim = c(0.0005, 1500),
     main= "dN/dlogDg at 11:07",
     log = "xy", pch = 20)
axis(1, at = c(0.1,1,10), labels = c(0.1,1,10))
axis(2, at = c(0.1,1,10,100,1000), labels = c(0.1,1,10,100,1000))
points(TS_merged_data[[i]][[j]][,1], TS_merged_data[[i]][[j]][,2],
       col = "black", pch = 20)
points(normalised_df[[i]][[j]][,1], normalised_df[[i]][[j]][,2],
       col = "hotpink", pch = 20)
points(SD_CESAM[[i]][[j]][,1], SD_CESAM[[i]][[j]][,10],
       col = "blue", pch = 20)
points(SD_Nephelometer[[i]][[j]][,1], SD_Nephelometer[[i]][[j]][,10],
       col = "gold", pch = 20)
legend("bottomleft", inset =.02, box.lty = 0,
       legend = c("Raw data", "After Interpolation", "Normalised values", "After GRIMM-SMPS loss correction", "Nephelometer SD"),
       pch = c(20,20,20,20,20),
       col = c("black", "green", "hotpink", "blue", "gold") ,cex=0.6)
mtext("n = 1.59, k = 0, X = 1", cex = 0.75, side=3,line=0, col = "red")


plot(interpolated_df[[i]][[j]][,1], interpolated_df[[i]][[j]][,4],
     col = "green",
     xaxt = "n",
     yaxt = "n",
     xlab = expression(paste("Diameter (",mu,"m)")),
     ylab = expression(paste("dN/dlogDg (",cm^{-3}*")")),
     xlim = c(0.01, 35),
     ylim = c(0.0005, 1500),
     main= "dN/dlogDg at 11:31",
     log = "xy", pch = 20)
axis(1, at = c(0.1,1,10), labels = c(0.1,1,10))
axis(2, at = c(0.1,1,10,100,1000), labels = c(0.1,1,10,100,1000))
points(TS_merged_data[[i]][[j]][,1], TS_merged_data[[i]][[j]][,4],
       col = "black", pch = 20)
points(normalised_df[[i]][[j]][,1], normalised_df[[i]][[j]][,4],
       col = "hotpink", pch = 20)
points(SD_CESAM[[i]][[j]][,1], SD_CESAM[[i]][[j]][,12],
       col = "blue", pch = 20)
points(SD_Nephelometer[[i]][[j]][,1], SD_Nephelometer[[i]][[j]][,12],
       col = "gold", pch = 20)
legend("bottomleft", inset =.02, box.lty = 0,
       legend = c("Raw data", "After Interpolation", "Normalised values", "After GRIMM-SMPS loss correction", "Nephelometer SD"),
       pch = c(20,20,20,20,20),
       col = c("black", "green", "hotpink", "blue", "gold") ,cex=0.6)
mtext("n = 1.59, k = 0, X = 1", cex = 0.75, side=3,line=0, col = "red")


plot(interpolated_df[[i]][[j]][,1], interpolated_df[[i]][[j]][,6],
     col = "green",
     xaxt = "n",
     yaxt = "n",
     xlab = expression(paste("Diameter (",mu,"m)")),
     ylab = expression(paste("dN/dlogDg (",cm^{-3}*")")),
     xlim = c(0.01, 35),
     ylim = c(0.0005, 1500),
     main= "dN/dlogDg at 11:55",
     log = "xy", pch = 20)
axis(1, at = c(0.1,1,10), labels = c(0.1,1,10))
axis(2, at = c(0.1,1,10,100,1000), labels = c(0.1,1,10,100,1000))
points(TS_merged_data[[i]][[j]][,1], TS_merged_data[[i]][[j]][,6],
       col = "black", pch = 20)
points(normalised_df[[i]][[j]][,1], normalised_df[[i]][[j]][,6],
       col = "hotpink", pch = 20)
points(SD_CESAM[[i]][[j]][,1], SD_CESAM[[i]][[j]][,14],
       col = "blue", pch = 20)
points(SD_Nephelometer[[i]][[j]][,1], SD_Nephelometer[[i]][[j]][,14],
       col = "gold", pch = 20)
legend("bottomleft", inset =.02, box.lty = 0,
       legend = c("Raw data", "After Interpolation", "Normalised values", "After GRIMM-SMPS loss correction", "Nephelometer SD"),
       pch = c(20,20,20,20,20),
       col = c("black", "green", "hotpink", "blue", "gold") ,cex=0.6)
mtext("n = 1.59, k = 0, X = 1", cex = 0.75, side=3,line=0, col = "red")


plot(interpolated_df[[i]][[j]][,1], interpolated_df[[i]][[j]][,9],
     col = "green",
     xaxt = "n",
     yaxt = "n",
     xlab = expression(paste("Diameter (",mu,"m)")),
     ylab = expression(paste("dN/dlogDg (",cm^{-3}*")")),
     xlim = c(0.01, 35),
     ylim = c(0.0005, 1500),
     main= "dN/dlogDg at 12:31",
     log = "xy", pch = 20)
axis(1, at = c(0.1,1,10), labels = c(0.1,1,10))
axis(2, at = c(0.1,1,10,100,1000), labels = c(0.1,1,10,100,1000))
points(TS_merged_data[[i]][[j]][,1], TS_merged_data[[i]][[j]][,9],
       col = "black", pch = 20)
points(normalised_df[[i]][[j]][,1], normalised_df[[i]][[j]][,9],
       col = "hotpink", pch = 20)
points(SD_CESAM[[i]][[j]][,1], SD_CESAM[[i]][[j]][,17],
       col = "blue", pch = 20)
points(SD_Nephelometer[[i]][[j]][,1], SD_Nephelometer[[i]][[j]][,17],
       col = "gold", pch = 20)
legend("bottomleft", inset =.02, box.lty = 0,
       legend = c("Raw data", "After Interpolation", "Normalised values", "After GRIMM-SMPS loss correction", "Nephelometer SD"),
       pch = c(20,20,20,20,20),
       col = c("black", "green", "hotpink", "blue", "gold") ,cex=0.6)
mtext("n = 1.59, k = 0, X = 1", cex = 0.75, side=3,line=0, col = "red")


######################################################################################################
######################################################################################################
######################################################################################################
# save the output
######################################################################################################
######################################################################################################
######################################################################################################

# For each instrument, add columns to help identify scenarios

add_Xnk_scenario <- function(Instrument_SD, merged_df){
  Instrument_SD <- add_column(Instrument_SD,X=merged_df$X[1],.before = "Dg")
  Instrument_SD <- add_column(Instrument_SD,n=merged_df$n[1],.after = "X")
  Instrument_SD <- add_column(Instrument_SD,k=merged_df$k[1],.after = "n")
  Instrument_SD <- add_column(Instrument_SD,Xnk=paste(merged_data[[i]][[j]]$X[1], as.numeric(as.character(merged_data[[i]][[j]]$n[1])), 
                                                      as.numeric(as.character(merged_data[[i]][[j]]$k[1]))),.after = "k")
  
  return(Instrument_SD)
}

# test the function for 1 instrument
i = 2
j = 3
test <- add_Xnk_scenario(SD_Aethalometer[[i]][[j]], merged_data[[i]][[j]])

######################################################################################################
# CESAM
######################################################################################################
# save CESAM data

# add scenario details to all df, instrument = Aethalometer
for (i in 1:length(SD_CESAM)) {
  for (j in 1:length(SD_CESAM[[i]])) {
    SD_CESAM[[i]][[j]] <- add_Xnk_scenario(SD_CESAM[[i]][[j]], merged_data[[i]][[j]])
    SD_CESAM[[i]][[j]] <- SD_CESAM[[i]][[j]][,-c(6:13)]
  }
}

# flatten the nested list
SD_CESAM_list <- do.call(c, SD_CESAM)

# check the output
identical(SD_CESAM[[3]][[1]], SD_CESAM_list[[295]])

# bind all the scenario in a large file
SD_CESAM_df <- convert_list_df(SD_CESAM_list)

# check the output

# flatten the nested list
SD_CESAM_list <- do.call(c, SD_CESAM)

# check the output
identical(SD_CESAM[[3]][[1]], SD_CESAM_list[[295]])

# bind all the scenario in a large file
SD_CESAM_df <- convert_list_df(SD_CESAM_list)

# save the output
test <- SD_CESAM_df

for (i in 5:length(test)) {
  for (j in 1:nrow(test)) {
    test[j,i] <- round(test[j,i], digits = 6) 
  }
}

# save the output
write.csv(test, "Z:/Clarissa/Data_Optical_calculation/R_Code/Code_final/ICEdust_samples/Maeli2/datasets/SD_output/original/Maeli2_SD_CESAM_original.csv", row.names = F)

######################################################################################################
# Aethalometer
######################################################################################################

# add scenario details to all df, instrument = Aethalometer
for (i in 1:length(SD_Aethalometer)) {
  for (j in 1:length(SD_Aethalometer[[i]])) {
    SD_Aethalometer[[i]][[j]] <- add_Xnk_scenario(SD_Aethalometer[[i]][[j]], merged_data[[i]][[j]])
    SD_Aethalometer[[i]][[j]] <- SD_Aethalometer[[i]][[j]][,-c(6:13)]
  }
}

# flatten the nested list
SD_Aethalometer_list <- do.call(c, SD_Aethalometer)

# check the output
identical(SD_Aethalometer[[3]][[1]], SD_Aethalometer_list[[295]])

# bind all the scenario in a large file
SD_Aethalometer_df <- convert_list_df(SD_Aethalometer_list)

# check the output

# check the output
identical(SD_CESAM_df[,5], SD_Aethalometer_df[,5])

# save the output
test <- SD_Aethalometer_df

for (i in 5:length(test)) {
  for (j in 1:nrow(test)) {
    test[j,i] <- round(test[j,i], digits = 6) 
  }
}

# save the output
write.csv(test, "Z:/Clarissa/Data_Optical_calculation/R_Code/Code_final/ICEdust_samples/Maeli2/datasets/SD_output/original/Maeli2_SD_Aethalometer_original.csv", row.names = F)

######################################################################################################
# Nephelometer
######################################################################################################

# add scenario details to all df, instrument = Nephelometer
for (i in 1:length(SD_Nephelometer)) {
  for (j in 1:length(SD_Nephelometer[[i]])) {
    SD_Nephelometer[[i]][[j]] <- add_Xnk_scenario(SD_Nephelometer[[i]][[j]], merged_data[[i]][[j]])
    SD_Nephelometer[[i]][[j]] <- SD_Nephelometer[[i]][[j]][,-c(6:13)]
  }
}

# flatten the nested list
SD_Nephelometer_list <- do.call(c, SD_Nephelometer)

# check the output
identical(SD_Nephelometer[[2]][[1]], SD_Nephelometer_list[[148]])

# bind all the scenario in a large file
SD_Nephelometer_df <- convert_list_df(SD_Nephelometer_list)

# check the output
identical(SD_CESAM_df[,5], SD_Nephelometer_df[,5])
identical(SD_Aethalometer_df[,5], SD_Nephelometer_df[,5])

# save the output
test <- SD_Nephelometer_df

for (i in 5:length(test)) {
  for (j in 1:nrow(test)) {
    test[j,i] <- round(test[j,i], digits = 6) 
  }
}

# save the output
write.csv(test, "Z:/Clarissa/Data_Optical_calculation/R_Code/Code_final/ICEdust_samples/Maeli2/datasets/SD_output/original/Maeli2_SD_Nephelometer_original.csv", row.names = F)

