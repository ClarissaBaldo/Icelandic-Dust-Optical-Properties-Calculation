# import packages needed
import os
import numpy as np
import pandas as pd
import PyMieScatt as mie
import math

################################################################################################################################
# import the dataset for all scenarios
Nep_all = pd.read_csv('Maeli2_SD_average_original.csv')

# multiply dN/dLog by dLog (1/64) to obtain dN, thus divide the measurements by 64
for i in range(5, Nep_all.shape[1]):
    Nep_all.iloc[:,i] = Nep_all.iloc[:,i]/64

# convert the unit to nanometers
Nep_all['Dg'] = Nep_all['Dg']*1000

# list all scenarios by "Xnk"
scenarios = pd.unique(Nep_all['Xnk'])

print("number of scenarios:",len(scenarios))
print("first scenario:",scenarios[0])
print("last scenario:",scenarios[-1])
print("all scenarios:",scenarios,sep = "\n")

# group data by scenarios
Nep_list = []

for i in range(len(scenarios)):
    Nep_list.append(Nep_all[Nep_all['Xnk'] == scenarios[i]])   

# reset index
for i in range(len(scenarios)):
    Nep_list[i] = Nep_list[i].reset_index(drop=True)
    
################################################################################################################################  
# Correction for truncation #

# define a function to calculate the increment d(theta)  
def calculate_d_theta(df):
    d_theta =[]
    for i in range(df.shape[0]-1):
        d_theta.append(df['theta'][i+1] - df['theta'][i])
    return d_theta

# define a function to calculate the average SU(theta)*sin(theta)
def calculate_SU_sin_theta_average(df):
    SU_sin_theta_average =[]
    for i in range(df.shape[0]-1):
        SU_sin_theta_average.append((df['SU_sin_theta'][i+1] + df['SU_sin_theta'][i])/2)
    return SU_sin_theta_average

#####################################################################################
# define the function to calculate Bsca for theta from 7-170: "calculate_Bsca_7170"

def calculate_Bsca_7170(INPUT_DATA,INPUT_WAVELENGTH):
    
    # define m
    m = complex(INPUT_DATA['n'][0], INPUT_DATA['k'][0])
    
    # calculate the scattering intensities using the angular function:
    output_angular_function = []
    for i in range(INPUT_DATA.shape[0]):
        output_angular_function.append(mie.ScatteringFunction(m,INPUT_WAVELENGTH,INPUT_DATA['Dg'][i]))
    
    # from the angular function output, save the results for theta and SU = SU(theta)*2  
    theta = []
    SU = []
    
    for i in range(len(output_angular_function)):
            theta.append(output_angular_function[i][0])
            SU.append(output_angular_function[i][3]*2)
    
    # calculate sin(theta) and SU*sin(theta)
    sin_theta = []
    SU_sin_theta = []
    
    for i in range(len(theta)):
        sin_theta.append([np.round(math.sin(data),8) for data in theta[i]])    
    
    SU_sin_theta = np.array(SU)*np.array(sin_theta)

    # build a list of dataframes with the output
    # each dataframe correspond to a single Dp
    Output_angular_function = []
    
    for i in range(len(theta)):
        Output_angular_function.append(pd.DataFrame({'theta': theta[i],
                                                     'SU':SU[i],
                                                     'sin_theta': sin_theta[i],
                                                     'SU_sin_theta':SU_sin_theta[i]}))
    # subset the data
    theta_7170 = [df[(df['theta'] >= (np.pi/180)*7) & (df['theta'] <= (np.pi/180)*170)] for df in Output_angular_function]
    
    # reset the index in theta_7170
    theta_7170 = [df.reset_index(drop=True) for df in theta_7170]

    # calculate d(theta) for all the dataframe in the list Output_angular_function, for theta from 7-170
    d_theta_7170 = [calculate_d_theta(df) for df in theta_7170]

    # calculate the average SU*sin(theta) for all the dataframe in the list Output_angular_function, for (theta) from 7-170
    SU_sin_theta_average_7170 = [calculate_SU_sin_theta_average(df) for df in theta_7170]
   
    # calculate SU*sin(theta)(average theta(i+1), theta(i))*d(theta), for theta from 7-170
    SU_sin_theta_dtheta_7170 = []

    for i in range(len(SU_sin_theta_average_7170)):
        SU_sin_theta_dtheta_7170.append(np.array(d_theta_7170[i])*np.array(SU_sin_theta_average_7170[i]))

    # calculate the  sum of SU(theta)*2*sin(theta)*d(theta), for (theta) from 7-170
    sum_SU_sin_theta_dtheta_7170 = []

    for i in range(len(SU_sin_theta_dtheta_7170)):
        sum_SU_sin_theta_dtheta_7170.append(SU_sin_theta_dtheta_7170[i].sum())

    # Calculate Qsca(Dp) = (wavelength^2)/(pi^2*Dp^2)*sum_SU_sin_theta_dtheta_7170(Dp)
    Q_7170 = []

    # the number of rows changes! use INPUT_DATA.shape[0]
    for i in range(INPUT_DATA.shape[0]):
        Q_7170.append(sum_SU_sin_theta_dtheta_7170[i]*((INPUT_WAVELENGTH**2)/((math.pi**2)*(INPUT_DATA['Dg'][i]**2))))
    
    # For each date-time, calculate Bsca for theta from 7-170
    # save a copy of the size distribution of the first scenario in the Nep_list
    Bsca_7170 = INPUT_DATA.copy()

    # add the calculated Qsca as the 5th column in the dataframe Bsca_7170
    Bsca_7170.insert(5,"Qsca",Q_7170)

    # now calculate dN(Dp)*Qsca(Dp)*((pi*(Dp)^2)/4), for the columns in Bsca_7170, 6-Bsca_7170.shape[1]
    for i in range(6,Bsca_7170.shape[1]):
        Bsca_7170.iloc[:,i] = Bsca_7170.iloc[:,i]*Bsca_7170.iloc[:,5]*((math.pi*Bsca_7170.iloc[:,4]**2)/4)
    
    # For each date time, calculate Bsca as the integral of Qsca, which is the sum of dN(Dp)*Qsca(Dp)*((pi*(Dp)^2)/4) contributions
    # Bsca(Mm^-1) = 10^(-6) * sum of dN(Dp)*Qsca(Dp)*((pi*(Dp)^2)/4) contributions
    Bsca_7170_time_serie = []

    for i in range(6,Bsca_7170.shape[1]):
        Bsca_7170_time_serie.append(Bsca_7170.iloc[:,i].sum()*(10**(-6)))
    
    # generate the column for Time steps
    Time = (pd.DataFrame(columns=['NULL'],index=pd.date_range('2019-01-21 11:07:00', '2019-01-21 13:55:00',freq='12T'))
                                         .index.strftime('%Y-%m-%dT%H:%M:%S').tolist())

    # generate a dataframe with the Bsca time series for theta from 7-170:
    Bsca_7170_time_serie = pd.DataFrame({'Time': Time,
                                     'Bsca_7170':Bsca_7170_time_serie})

    return Bsca_7170_time_serie

#####################################################################################
# apply the function for all scenarios and save the results
Bsca_7170_450nm = [calculate_Bsca_7170(data, 450) for data in Nep_list]

#####################################################################################
# combine the "Xnk" marker and the results
for i in range(len(scenarios)):
    Bsca_7170_450nm[i].insert(0,'Xnk', scenarios[i])
    
# combine the results as a single dataframe
df_450nm = pd.concat(Bsca_7170_450nm)
df_450nm

# save the output
df_450nm.to_csv('Maeli2_Nep_mie_7170_450nm_original.csv')