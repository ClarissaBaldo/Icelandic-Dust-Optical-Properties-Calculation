#############################################################################################################
# after running the codes parts 1-3
# calculate the mean, standard deviation (SD), and relative standard deviation (RSD)
# of X, n, k, and SSA selected based on the RMSE values 
# obtained from the comparison between Bsca and Babs modelled and measured
#############################################################################################################

library(matrixStats)

i=1

# Test for 1 WL and 1 time
Bsca_OG_df = Bsca_OG_results[[i]]
Bsca_mSD1_df = Bsca_mSD1_results[[i]]
Bsca_pSD1_df = Bsca_pSD1_results[[i]]
Babs_OG_df = Babs_OG_results[[i]]
Babs_mSD1_df = Babs_mSD1_results[[i]]
Babs_pSD1_df = Babs_pSD1_results[[i]]

SSA_df = SSA_Xnk_WL[[i]]

# get the quantile range from 0.1 to 0.5 by 0.1
range_quantile <- seq(0.1,0.5,0.1)

# subset the data using different quantiles for RMSE
Bsca_OG_subset_list <- list()
Bsca_mSD1_subset_list <- list()
Bsca_pSD1_subset_list <- list()
Babs_OG_subset_list <- list()
Babs_mSD1_subset_list <- list()
Babs_pSD1_subset_list <- list()

for (i in 1:length(range_quantile)) {
  Bsca_OG_subset_list[[i]] <- Bsca_OG_df[which((Bsca_OG_df$RMSE < quantile(Bsca_OG_df$RMSE, probs = range_quantile[i])) & (Bsca_OG_df$R2 >=0.7)),]
  Bsca_mSD1_subset_list[[i]] <- Bsca_mSD1_df[which((Bsca_mSD1_df$RMSE < quantile(Bsca_mSD1_df$RMSE, probs = range_quantile[i])) & (Bsca_mSD1_df$R2 >=0.7)),]
  Bsca_pSD1_subset_list[[i]] <- Bsca_pSD1_df[which((Bsca_pSD1_df$RMSE < quantile(Bsca_pSD1_df$RMSE, probs = range_quantile[i])) & (Bsca_pSD1_df$R2 >=0.7)),]
  Babs_OG_subset_list[[i]] <- Babs_OG_df[which((Babs_OG_df$RMSE < quantile(Babs_OG_df$RMSE, probs = range_quantile[i])) & (Babs_OG_df$R2 >=0.7)),]
  Babs_mSD1_subset_list[[i]] <- Babs_mSD1_df[which((Babs_mSD1_df$RMSE < quantile(Babs_mSD1_df$RMSE, probs = range_quantile[i])) & (Babs_mSD1_df$R2 >=0.7)),]
  Babs_pSD1_subset_list[[i]] <- Babs_pSD1_df[which((Babs_pSD1_df$RMSE < quantile(Babs_pSD1_df$RMSE, probs = range_quantile[i])) & (Babs_pSD1_df$R2 >=0.7)),]
}

# intersect Bsca and Babs results (OG and plus/minus error)
Xnk_comparison <- list()

for (i in 1:length(range_quantile)) {
  # get row_numbers
  Xnk_comparison[[i]] <- Bsca_OG_subset_list[[i]][which((Bsca_OG_subset_list[[i]]$Xnk %in% Bsca_mSD1_subset_list[[i]]$Xnk) &
                                                          (Bsca_OG_subset_list[[i]]$Xnk %in% Bsca_pSD1_subset_list[[i]]$Xnk) &
                                                          (Bsca_OG_subset_list[[i]]$Xnk %in% Babs_OG_subset_list[[i]]$Xnk) &
                                                          (Bsca_OG_subset_list[[i]]$Xnk %in% Babs_mSD1_subset_list[[i]]$Xnk) &
                                                          (Bsca_OG_subset_list[[i]]$Xnk %in% Babs_pSD1_subset_list[[i]]$Xnk)),]
}

# if nrow(Xnk_comparison) == 0, there are no common results 
# in the list, replace df with nrow == 0 with NA

for (i in 1:length(Xnk_comparison)) {
  if (nrow(Xnk_comparison[[i]]) == 0) {
    
    Xnk_comparison[[i]] <- NA
    
  }
}

# remove NAs from the list
Xnk_comparison <- Xnk_comparison[!is.na(Xnk_comparison)]

# at the top of the list is the term with the lowest uncertainty
Xnk_comparison_results <- Xnk_comparison[[1]]

# to determine SSA model subset SSA results according to Xnk_comparison
SSA_model_results <- SSA_df[which(SSA_df$Xnk %in% Xnk_comparison_results$Xnk),]

# calculate the average values for k, n, X, SSA and RMSE

k_WL <- data.frame(WL = Xnk_comparison_results$WL[1], k_mean = mean(Xnk_comparison_results$k, na.rm = T), 
                   k_sd = sd(Xnk_comparison_results$k, na.rm = T), k_rsd = sd(Xnk_comparison_results$k, na.rm = T)/mean(Xnk_comparison_results$k, na.rm = T)*100)

n_WL <- data.frame(n_mean = mean(Xnk_comparison_results$n, na.rm = T), 
                   n_sd = sd(Xnk_comparison_results$n, na.rm = T), n_rsd = sd(Xnk_comparison_results$n, na.rm = T)/mean(Xnk_comparison_results$n, na.rm = T)*100)

X_WL <- data.frame(X_mean = mean(Xnk_comparison_results$X, na.rm = T), 
                   X_sd = sd(Xnk_comparison_results$X, na.rm = T), X_rsd = sd(Xnk_comparison_results$X, na.rm = T)/mean(Xnk_comparison_results$X, na.rm = T)*100)

SSA_WL <- data.frame(SSA_mean = mean(SSA_model_results$SSA, na.rm = T), 
                     SSA_sd = sd(SSA_model_results$SSA, na.rm = T), SSA_rsd = sd(SSA_model_results$SSA, na.rm = T)/mean(SSA_model_results$SSA, na.rm = T)*100)

RMSE_WL <- data.frame(RMSE_mean = mean(Xnk_comparison_results$RMSE, na.rm = T), 
                      RMSE_sd = sd(Xnk_comparison_results$RMSE, na.rm = T), RMSE_rsd = sd(Xnk_comparison_results$RMSE, na.rm = T)/mean(Xnk_comparison_results$RMSE, na.rm = T)*100)

test1 <- data.frame(cbind(k_WL, n_WL, X_WL, SSA_WL, RMSE_WL))  

#############################################################################################################
# build a function to perform calculation for all dfs
#############################################################################################################

# remove terms from R environment except for:
rm(list=setdiff(ls(), c("Babs_OG_results", "Babs_pSD1_results","Babs_mSD1_results",
                        "Bsca_OG_results", "Bsca_pSD1_results","Bsca_mSD1_results",
                        "test1", "SSA_Xnk_WL")))


get_best_Xnk <- function(Bsca_OG_df, Bsca_mSD1_df, Bsca_pSD1_df, 
                         Babs_OG_df, Babs_mSD1_df, Babs_pSD1_df,
                         SSA_df) {
  
  # get the quantile range from 0.1 to 0.5 by 0.1
  range_quantile <- seq(0.1,0.5,0.1)
  
  # subset the data using different quantiles for RMSE
  Bsca_OG_subset_list <- list()
  Bsca_mSD1_subset_list <- list()
  Bsca_pSD1_subset_list <- list()
  Babs_OG_subset_list <- list()
  Babs_mSD1_subset_list <- list()
  Babs_pSD1_subset_list <- list()
  
  for (i in 1:length(range_quantile)) {
    Bsca_OG_subset_list[[i]] <- Bsca_OG_df[which((Bsca_OG_df$RMSE < quantile(Bsca_OG_df$RMSE, probs = range_quantile[i])) & (Bsca_OG_df$R2 >=0.7)),]
    Bsca_mSD1_subset_list[[i]] <- Bsca_mSD1_df[which((Bsca_mSD1_df$RMSE < quantile(Bsca_mSD1_df$RMSE, probs = range_quantile[i])) & (Bsca_mSD1_df$R2 >=0.7)),]
    Bsca_pSD1_subset_list[[i]] <- Bsca_pSD1_df[which((Bsca_pSD1_df$RMSE < quantile(Bsca_pSD1_df$RMSE, probs = range_quantile[i])) & (Bsca_pSD1_df$R2 >=0.7)),]
    Babs_OG_subset_list[[i]] <- Babs_OG_df[which((Babs_OG_df$RMSE < quantile(Babs_OG_df$RMSE, probs = range_quantile[i])) & (Babs_OG_df$R2 >=0.7)),]
    Babs_mSD1_subset_list[[i]] <- Babs_mSD1_df[which((Babs_mSD1_df$RMSE < quantile(Babs_mSD1_df$RMSE, probs = range_quantile[i])) & (Babs_mSD1_df$R2 >=0.7)),]
    Babs_pSD1_subset_list[[i]] <- Babs_pSD1_df[which((Babs_pSD1_df$RMSE < quantile(Babs_pSD1_df$RMSE, probs = range_quantile[i])) & (Babs_pSD1_df$R2 >=0.7)),]
  }
  
  # intersect Bsca and Babs results (OG and plus/minus error)
  Xnk_comparison <- list()
  
  for (i in 1:length(range_quantile)) {
    # get row_numbers
    Xnk_comparison[[i]] <- Bsca_OG_subset_list[[i]][which((Bsca_OG_subset_list[[i]]$Xnk %in% Bsca_mSD1_subset_list[[i]]$Xnk) &
                                                            (Bsca_OG_subset_list[[i]]$Xnk %in% Bsca_pSD1_subset_list[[i]]$Xnk) &
                                                            (Bsca_OG_subset_list[[i]]$Xnk %in% Babs_OG_subset_list[[i]]$Xnk) &
                                                            (Bsca_OG_subset_list[[i]]$Xnk %in% Babs_mSD1_subset_list[[i]]$Xnk) &
                                                            (Bsca_OG_subset_list[[i]]$Xnk %in% Babs_pSD1_subset_list[[i]]$Xnk)),]
  }
  
  # if nrow(Xnk_comparison) == 0, there are no common results 
  # in the list, replace df with nrow == 0 with NA
  
  for (i in 1:length(Xnk_comparison)) {
    if (nrow(Xnk_comparison[[i]]) == 0) {
      
      Xnk_comparison[[i]] <- NA
      
    }
  }
  
  # remove NAs from the list
  Xnk_comparison <- Xnk_comparison[!is.na(Xnk_comparison)]
  
  # at the top of the list is the term with the lowest uncertainty
  Xnk_comparison_results <- Xnk_comparison[[1]]
  
  # to determine SSA model subset SSA results according to Xnk_comparison
  SSA_model_results <- SSA_df[which(SSA_df$Xnk %in% Xnk_comparison_results$Xnk),]
  
  # calculate the average values for k, n, X, SSA and RMSE
  
  k_WL <- data.frame(WL = Xnk_comparison_results$WL[1], k_mean = mean(Xnk_comparison_results$k, na.rm = T), 
                     k_sd = sd(Xnk_comparison_results$k, na.rm = T), k_rsd = sd(Xnk_comparison_results$k, na.rm = T)/mean(Xnk_comparison_results$k, na.rm = T)*100)
  
  n_WL <- data.frame(n_mean = mean(Xnk_comparison_results$n, na.rm = T), 
                     n_sd = sd(Xnk_comparison_results$n, na.rm = T), n_rsd = sd(Xnk_comparison_results$n, na.rm = T)/mean(Xnk_comparison_results$n, na.rm = T)*100)
  
  X_WL <- data.frame(X_mean = mean(Xnk_comparison_results$X, na.rm = T), 
                     X_sd = sd(Xnk_comparison_results$X, na.rm = T), X_rsd = sd(Xnk_comparison_results$X, na.rm = T)/mean(Xnk_comparison_results$X, na.rm = T)*100)
  
  SSA_WL <- data.frame(SSA_mean = mean(SSA_model_results$SSA, na.rm = T), 
                       SSA_sd = sd(SSA_model_results$SSA, na.rm = T), SSA_rsd = sd(SSA_model_results$SSA, na.rm = T)/mean(SSA_model_results$SSA, na.rm = T)*100)
  
  RMSE_WL <- data.frame(RMSE_mean = mean(Xnk_comparison_results$RMSE, na.rm = T), 
                        RMSE_sd = sd(Xnk_comparison_results$RMSE, na.rm = T), RMSE_rsd = sd(Xnk_comparison_results$RMSE, na.rm = T)/mean(Xnk_comparison_results$RMSE, na.rm = T)*100)
  
  summary_results <- data.frame(cbind(k_WL, n_WL, X_WL, SSA_WL, RMSE_WL))  
  
  return(summary_results)
  
}

# check the function
i=1

test2 <- get_best_Xnk(Bsca_OG_results[[i]],Bsca_mSD1_results[[i]],Bsca_pSD1_results[[i]],
                      Babs_OG_results[[i]],Babs_mSD1_results[[i]],Babs_pSD1_results[[i]],
                      SSA_Xnk_WL[[i]])

identical(test1,test2)

# apply the function to all dfs
Xnk_results_WL <- list()

for (i in 1:length(Bsca_OG_results)) {
  Xnk_results_WL[[i]] <- get_best_Xnk(Bsca_OG_results[[i]],Bsca_mSD1_results[[i]],Bsca_pSD1_results[[i]],
                                      Babs_OG_results[[i]],Babs_mSD1_results[[i]],Babs_pSD1_results[[i]],
                                      SSA_Xnk_WL[[i]])
}

# check the results
i=1

identical(test1, Xnk_results_WL[[i]])

#############################################################################################################
# output the data
Xnk_final <- list.rbind(Xnk_results_WL)

write.csv(Xnk_final, "Z:/Clarissa/Data_Optical_calculation/R_Code/Code_final/ICEdust_samples/Maeli2/Nep_Aet_Mie/original/data/Maeli2_Xnk_average_original.csv", row.names = F)
