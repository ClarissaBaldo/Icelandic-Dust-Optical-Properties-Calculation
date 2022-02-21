#####################################################################################
# Python codes for Clarissa's optical calcuations
# online user guide: https://pymiescatt.readthedocs.io/en/latest/

# import packages needed
import os
import numpy as np
import pandas as pd
import PyMieScatt as mie

#####################################################################################
#####################################################################################

# move to working directory
os.chdir("Z:\\Clarissa\\Data_Optical_calculation\\R_Code\\Code_final\\ICEdust_samples\\Maeli2\\datasets\\SD_output\\original")

# import the dataset for all scenarios
Aet_all = pd.read_csv('Maeli2_SD_average_original.csv')

# multiply dN/dLog by dLog (1/64) to obtain dN, thus divide the measurements by 64
for i in range(5, Aet_all.shape[1]):
  Aet_all.iloc[:,i] = Aet_all.iloc[:,i]/64

# convert the unit to nanometers
Aet_all['Dg'] = Aet_all['Dg']*1000

# list all scenarios by "Xnk"
scenarios = pd.unique(Aet_all['Xnk'])

# group data by scenarios
Aet_list = []

for i in range(len(scenarios)):
  Aet_list.append(Aet_all[Aet_all['Xnk'] == scenarios[i]]) 

# reset index
for i in range(len(scenarios)):
  Aet_list[i] = Aet_list[i].reset_index(drop=True)

#####################################################################################
#####################################################################################
# "perform_mie_calculation"

def perform_mie_calculation(INPUT_DATA,INPUT_WAVELENGTH):
  """This is to prepare the data using the size distribution from a certain instrument, to calcualte the Mie coefficients, 
       and save relevant results. The INPUT_DATA is a data frame for a single scenario that needs to be processed.
       Use a for loop later to apply this function to all scenarios"""

# for each datetime column, we use a small dataframe to hold the data
df = []

for i in range(5,INPUT_DATA.shape[1]):     # df.shape[1] is the ncol, df.shape[0] is the nrow
  df.append(INPUT_DATA.iloc[:,[3,4,i]])  # the first datetime: df[0], pick column 3,4 and i (in reality, 4,5...)  
# last datetime: df[60]

# remove "NA" for each datetime point
for i in range(len(df)):
  df[i] = df[i].dropna()

# with "NA" removed, we can calculate the mie coefficients for each datetime point using "Mie_SD"
# Mie_SD(m, wavelength, sizeDistributionDiameterBins, sizeDistribution[, nMedium=1.0, asDict=False])
# m and wavelength are fixed for this dataset
# "sizeDistributionDiameterBins" = "Dg"

# set the wavelength   
wavelength = INPUT_WAVELENGTH

# set the m value
m = complex(INPUT_DATA['n'][0],INPUT_DATA['k'][0])

# now start to calculate mie coefficients for each time step 
mie_results = []

for i in range(len(df)):
  mie_results.append(mie.Mie_SD(m,wavelength,df[i]['Dg'],df[i].iloc[:,2],asDict=True))

# extract relevant results from mie coefficients        
Bext = []
Bsca = []
Babs = []

for i in range(len(mie_results)):
  Bext.append(mie_results[i]['Bext'])
Bsca.append(mie_results[i]['Bsca'])
Babs.append(mie_results[i]['Babs'])  

# genereate the column for Time steps
Time = (pd.DataFrame(columns=['NULL'],index=pd.date_range('2019-01-21 11:07:00', '2019-01-21 13:55:00',freq='12T'))
        .index.strftime('%Y-%m-%dT%H:%M:%S').tolist())

Time = [data.replace("T"," ") for data in Time]

# combine relevant results at all time steps
Instrument_output_data = {'Time': Time,
  'Bext': Bext, 
  'Bsca': Bsca,
  'Babs':Babs}
Instrument_output_data = pd.DataFrame(data=Instrument_output_data)
return Instrument_output_data

#####################################################################################
# apply the function for all scenarios and save the results

# data from all scenarios are held at: Aet_list
Aet_output_450nm = [perform_mie_calculation(data, 450) for data in Aet_list]

# combine the "Xnk" marker and the results
for i in range(len(scenarios)):
  Aet_output_450nm[i].insert(0,'Xnk', scenarios[i]) 

# combine the results as a single dataframe
df_450nm = pd.concat(Aet_output_450nm)

# save the output
# move to working directory
os.chdir("Z:\\Clarissa\\Data_Optical_calculation\\R_Code\\Code_final\\ICEdust_samples\\Maeli2\\datasets\\Python_output\\original\\Mie_coef")
df_450nm.to_csv('Maeli2_Aet_mie_450nm_original.csv')
