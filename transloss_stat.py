#!/usr/bin/env python3

import pandas as pd
import numpy as np
from eofs.standard import Eof #https://ajdawson.github.io/eofs/latest/userguide/index.html
#from eofs.examples import example_data_path
from numpy import savetxt
from scipy.stats import norm
from scipy import interpolate
import matplotlib.pyplot as plt
from KDEpy import FFTKDE  #https://kdepy.readthedocs.io/en/latest/
import seaborn as sns
from numpy.random import random
from random import randint
import glob
import os
from scipy.stats import skew
from scipy.stats import kurtosis
import subprocess

sns.set()

az = 22.5 #Propagation direction in degrees
rnge = [20,40,60,80] # Range in Km

# A function to get a data matrix from historical data
def get_datamatrix():
    #convert the final.dat to pandas dataframe
    final = pd.read_csv('final.dat')
    #Transpose the dataframe 'final'
    final_T = final.T
    #convert the dataframe to a numpy array to serve as the data matrix
    a = final_T.values
    a_mean_vector = np.mean(a, axis = 0)
    #center the data 
    a = np.subtract(a, a_mean_vector)
    return a, a_mean_vector
    

a, a_mean_vector = get_datamatrix()

# Create an instance of the Eof solver class
solver = Eof(a, center = False)

# A function to get optimal truncation location for truncated SVD or number of principal components
def get_no_of_true():
    #find the svd for the data matrix 'a'
    U, S, Vh = np.linalg.svd(a)
    ymed=np.median(S)
    shape = a.shape
    no_rows=shape[0]
    no_columns=shape[1]
    B=no_rows/no_columns
    w=0.56*B**3 - 0.95*B**2 + 1.82*B + 1.43
    T=ymed*w
    #Number of singular values greater than the optimal truncation location.
    no_of_true=np.count_nonzero((S > T))
    return no_of_true  

no_of_true = get_no_of_true()


# A function to extract the altitude and density (rho) from atmospheric profile. Note the atmospheric profile must be named 'Z_RHO'
def get_Z_and_RHO():
    Z_RHO = []
    col_names = ['Z(km)', 'T(degk)','U(m/s)','V(m/s)', 'RHO(mbar)','P(kg/m3)']
    examp = pd.read_csv('Z_RHO', skiprows = 8, header = None, sep = ' ')
    examp.dropna(axis = 1, inplace = True)
    examp.columns = col_names
    Z = examp['Z(km)']
    RHO = examp['RHO(mbar)']
    Z = np.array(list(Z))
    Z_RHO.append(Z)
    RHO = np.array(list(RHO))
    Z_RHO.append(RHO)
    return Z_RHO

Z_RHO = get_Z_and_RHO()


# A function to prepend the formatted headers to profile (already in csv format), so epape can use it.
def prepend_multiple_lines(file_name, list_of_lines):
    """Insert given list of strings as a new lines at the beginning of a file"""
    # define name of temporary dummy file
    dummy_file = file_name + '.bak'
    # open given original file in read mode and dummy file in write mode
    with open(file_name, 'r') as read_obj, open(dummy_file, 'w') as write_obj:
        # Iterate over the given list of strings and write them to dummy file as lines
        for line in list_of_lines:
            write_obj.write(line + '\n')
        # Read lines from original file one by one and append them to the dummy file
        for line in read_obj:
            write_obj.write(line)
        # remove original file
        os.remove(file_name)
        # Rename dummy file as the original file
        os.rename(dummy_file, file_name)


# A function to get a sample effective sound speed
def get_sample_Ceff():
    eff_soundspeeds = []
    ceff = 0
    for i in range(no_of_true):
        z = Z_RHO[0]
        #Get coeffiecients of eof
        pc=solver.pcs(npcs=(i+1))
        # Select the ith coeffiecient of eof
        pc=pc[:,i] #This is in a row form
        #Get their KDE to use the probability densities from KDE to generate a cumulative distribution function, cdf
        x, y = FFTKDE(kernel='gaussian', bw='ISJ').fit(pc).evaluate()
        #Get samples via inverse transform sampling method
        cdf_y = np.cumsum(y) #cumulative distribution function, cdf 
        cdf_y = cdf_y/cdf_y.max() #takes care of normalizing cdf to 1.0
        inverse_cdf = interpolate.interp1d(cdf_y,x,fill_value="extrapolate") #gets inverse function of cdf from cdf_y
        uniform_samples = random(int(100)) #generate values in [0,1) having uniform samples (uniform distribution)
        pc_samp = inverse_cdf(uniform_samples) # sampled pc
        #Generate a random number generator to randomly sample from pc_samp
        r = randint(0,99)
        #randomly sample from pc_samp in a try except block
        try:
            rand_pc = pc_samp[r]
        except IndexError:
            r1 = randint(0, 100)
            rand_pc = pc_samp[r1]
        #Generate correspnding EOF
        eof = solver.eofs(neofs=(i+1))
        eof=eof[i]
        #Do a basis expansion with the ith eof
        ith_basis_expnsn = rand_pc * eof
        #Get effective sound speed from summation of basis function 
        #(summation occurs by accumulating after each iteration)
        ceff += ith_basis_expnsn
    ceff += a_mean_vector
    eff_soundspeeds.append(ceff)
    return eff_soundspeeds

eff_soundspeeds = get_sample_Ceff()
   
# A function to develop atmospheric profiles for ePape     
def get_profile():
    #create an atmospheric profile for the sample ceff by first creating a dataframe
    profile = pd.DataFrame({'Z':Z_RHO[0], 'RHO':Z_RHO[1], 'CEFF':eff_soundspeeds[0]}, columns = None, index = np.arange(1,91))
    #convert the dataframe to csv
    profile.to_csv('profile', sep = ' ', header = False, index = False)
    
    
# A function to prepend formatted headers to atmospheric profiles for ePape   
def prepend_atmos_spec():
    #A list containing the formatted headers. This is an argument in the function for prepending
    list_of_lines = ['#% 0, Z0, m, 117.3', '#% 1, Z, km', '#% 2, RHO, kg/m3', '#% 3, CEFF, m/s']
    #Calling the function for prepending
    #prepend_multiple_lines("profile".format(v), list_of_lines)
    prepend_multiple_lines("profile", list_of_lines)
    
# A function to run ePape    
def run_epape(x):
    #Run epape to generate the transfer functions
    epape_cmd = ["ePape", "--singleprop", "--starter self", "--atmosfile profile", "--freq 5", "--azimuth {}".format(az), "--maxrange_km {}".format(x)]
    epape_run = subprocess.Popen(epape_cmd)
    epape_run.wait()
    #This function returns a file called tloss_1d.pe
    
    
# A function for feature engineering. Transmission loss feature is developed from ePape result    
def get_sorted_tranloss(y):
    #create a dataraframe from the tloss_1d.pe file generated by epape
    tf = pd.read_csv('tloss_1d.pe', header = None, sep = ' ')
    #create the signal attenuation, magnitude of transfer functions in a fourth column 
    tf[4] = np.sqrt((tf[2])**2 + (tf[3])**2)
    #signal attenuation as a numpy array
    sig_att = (tf[4]).values
    #change signal attenuation to transmission loss
    trans_loss = 20*(np.log10(sig_att))
    trans_loss_transpose = trans_loss.T
    list_of_trans_loss = list(trans_loss)
    #Add the transmission loss to a list to store all generated transmission loss
    y = y + list_of_trans_loss
    sorted_trans_loss4kde = sorted(y)
    return sorted_trans_loss4kde
    
    
# A function to plot KDE of transmission loss    
def plot_kde(x_trans, y_prob, x_dis):
    plt.figure()
    plt.plot(x_trans, y_prob)
    plt.xlabel('transmission_loss_(db)')
    plt.ylabel('probability density')
    plt.figtext(0.7, 0.8, "az={}deg".format(az))
    plt.figtext(0.7, 0.7, "range={}km".format(x_dis))
    
    
# A function to save the KDE plot   
def save_fig(x):
    plt.savefig('transloss_kde_{}km'.format(x))
    
# A function to get transmission loss KDE for a range    
def get_transloss_kde_for_a_range(x):
    trans_loss_4_kde = []
    for v in range(1,16):
        get_sample_Ceff()
        get_profile()
        prepend_atmos_spec()
        run_epape(x)
        sortddd_transloss = get_sorted_tranloss(trans_loss_4_kde)
    trans_loss_x, prob_den = FFTKDE(kernel='gaussian', bw='ISJ').fit(sortddd_transloss).evaluate()
    plot_kde(trans_loss_x, prob_den, x)
    save_fig(x)
    

# A function to get transmission loss KDE for ranges in rnge list
def get_transloss_kde_for_given_range():
    for dis in rnge:
        get_transloss_kde_for_a_range(dis)
        

#Run the program        
get_transloss_kde_for_given_range()
