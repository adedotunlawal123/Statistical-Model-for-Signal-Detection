#!/usr/bin/env python3


#A data analysis script to generate effective sound speed from atmospheric profiles.
#This script requires that the atmospheric profiles from which we want to the extract effective sound speed are in the same directory as the script.
#This script is for atmospheric profiles where Z0 is available.
#This script generates another script, dd.py, which you should run to get the final.dat. final.dat is to be used as data matrix in transloss_stat.py

import glob
import os

az = 22.5 # The angle of propagation for getting effective sound speed.
commands=open("dd.py", 'w')
print("""#!/usr/bin/env python3
import glob
import pandas as pd
import numpy as np
import math
f = []""", file=commands)

filelist1 = glob.iglob("*dat")
for datafile in sorted(filelist1): # sort the list above so we can have the columns in order
    print('a = pd.read_csv("{input}")'.format( input=datafile ), file=commands ) 
    print('buu = a.drop(a.index[[0,1,2,3,4,5,6]])', file=commands) # Drop the 2nd to 8th row the csv file, although the 1st to 7th in the dataframe above. Note that the first row became the column names for the dataframe above. 
    print('buu.to_csv("{output}")'.format( output=datafile ), file=commands) # convert this buu dataframe above to csv
    
filelist2 = glob.iglob("*dat")
for datafile in sorted(filelist2): # sort filelist2 so we can keep the order we need to produce a timeseries
    print('cc = pd.read_csv("{input_1}", names=[1,2,3,4,5,6,7,8,9,10,11,12], delimiter = " ")'.format( input_1 = datafile ), file=commands) # note that the column names that are even numbers are Nan
    print('qw = cc.drop(cc.index[0])', file=commands)
    print('ss = qw.dropna(axis="columns", how ="all")', file=commands) # drop all Nan
    print('del ss[11]', file=commands) # delete the column named 11 (it has some commas in the entries; hence making it difficult to convert the dataframe entries to float)
    print('ss = ss.astype(float)', file=commands) # convert all the dataframe entries to float
    print('ss[3] = np.sqrt((1.4*287.0*ss[3])) + ss[5]*math.sin(({}*(math.pi/180))) + ss[7]*math.cos(({}*(math.pi/180)))'.format(az,az), file=commands) # change the temperature column to effective sound speed, note the angle is hardcoded.
    print('ss.drop(ss.columns[[0, 2, 3, 4]], axis = 1, inplace = True)', file=commands) # drop the all columns except the one we need by column index and not column names
    we = datafile # assign a name to the name of the individual files
    wee = we[5:13] # pick out the number (date and time) of individual files to be used to name the resulting file after several cleaning.
    print('ss_{} = ss'.format( wee ), file=commands) # rename the resulting dataframe using the number(date and time).
    print('f.append(ss_{})'.format(wee), file=commands) # append the individual dataframe to f, a list, above.
    

print('final = pd.concat(f, axis =1)', file=commands) # concatenate all the dataframes in f and name it final
print('final.to_csv("final.dat", index = False)', file=commands)
commands.close()
