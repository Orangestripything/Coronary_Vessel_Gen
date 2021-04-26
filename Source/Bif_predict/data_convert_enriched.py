#SCRIPT DATA_CONVERT_ENRICHED concatinates all enriched segment connection matrices and converts them into appropriate form.

#Import modules
import os
import pandas as pd
import numpy as np


def load_data(file):
#LOAD_DATA loads selected file as a pandas dataframe for processing

  #read the selected file into a pandas dataframe
  df = pd.read_csv(rawdata_path + file,header=None,names = rawcolumn)
  
  #shifts index by +1 
  df.index = df.index + 1

  return df


def normalise_df (df):
#NORMALISE_DF converts connection points to x/y coordinates relative to terminal node/bifurcation point.

  #define parent x and y columns
  parentxcolumns = ['P2x','GPx']
  parentycolumns = ['P2y','GPy']

  #define parent x and y columns
  xcolumns = ['Px','Bifx','D1x','D2x']
  ycolumns = ['Py','Bify','D1y','D2y']

  # ^make sure Px and Py are at the end otherwise it will be minusing 0 from other columns 

  #normalise x coordinate of parent connection points
  for column in parentxcolumns:
      df[column] = df[column] - df['Px']

  #normalise y coordinate of parent connection points
  for column in parentycolumns:
      df[column] = df[column] - df['Py']

  #normalise x coordinate of points
  for column in xcolumns:
      df[column] = df[column] - df['D2x']

  #normalise y coordinate of points
  for column in ycolumns:
      df[column] = df[column] - df['D2y']

  return df


#MAIN FUNCTION
#change working directory to 06_Data
os.chdir('../../06_Data')

#select path for rawdata and output csv
rawdata_path = (os.getcwd() + "/RawData/Enriched/3000seg/Segenriched/") #path to load
write_path = (os.getcwd() + "/ProcessedData/WP2T1/3000_seg_enriched_dataset.csv") #path to write to
files = os.listdir(rawdata_path)  #get list of files in rawdata_path directory

#Display number of files to process
print('There are', len(files),'files')

#Set feature columns
rawcolumn = ['Px','Py','D1x','D1y','D2x','D2y','Radius','circleradius','totseg','P2x','P2y','GPx','GPy','OrigSeg','l1','l2','l3','beta','alpha','Bifx','Bify']

#Initialise training dataset dataframe and counter
traindf = pd.DataFrame(columns=rawcolumn)
n = 1

#Loop through all files in directory
for file in files:
    df = load_data(file)  #load file as a pandas dataframe
    traindf = traindf.append(df,ignore_index=True)  #append dataframe to traindf
    print('%d of %d files to process...' %(n,len(files))) #display progress
    n += 1  #increase increment

#Convert connection points to x/y coordinates relative to terminal node/bifurcation point.
traindf = normalise_df(traindf)

#Display number of observations in training dataset
print('Dataset = %d' %(len(traindf)))

#Save pandas dataframe to write path as csv file
traindf.to_csv(write_path)

print('Completed!')