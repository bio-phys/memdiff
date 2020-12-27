# System output and input
import sys
import os.path
import pickle

# Mathematical operations
import numpy as np
import scipy as sp
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from matplotlib import colors

# Functions for membrane diffusion
import memdiff.oseen as mdv
import memdiff.immersed as mdc

# Natural constants
kB     = 1.3806485*10**-23  # Boltzmann constant [J/K]

# unit factors
cf     = 10**-5             # dc [nm^2/ns] * cf = dc [cm^2/s]
m2     = 10**4              # dc [m^2/s] * m2 = dc [cm^2/s]  



####################################
#   ITERATIVE CALCULATION OF D_0   #
####################################


def chi2(params,data_sets):
 
    eta_m  = params[0]
    d0     = params[1:] # list with length = number of data sets
    
    # Initialize chi2
    chi2 = 0.0
    # Loop through all data sets
    for j,ds in enumerate(data_sets):
        
        # Initialize the estimated D
        d_est  = np.empty(len(ds.dc))
        # Loop through all data points ... 
        for i in range(len(d_est)):
            # ... and estimate D  
            if ds.approx: 
                # ... approximately.
                d_est[i] = d0[j] + mdv.deltaD_approximation(ds.temperature,ds.eta_f,eta_m,ds.edge[i],ds.height[i])
            else: 
                # ... with the full theory.
                if ds.mono:
                    d_est[i] = d0[j] + mdv.deltaD_mono(ds.temperature,ds.eta_f,eta_m,ds.b,ds.edge[i],ds.height[i],100)
                else:
                    d_est[i] = d0[j] + mdv.deltaD(ds.temperature,ds.eta_f,eta_m,ds.edge[i],ds.height[i],100)
                    
        # Add the chi2 contribution for the current data set
        chi2 += np.sum( np.square((ds.dc - d_est)/ds.err) )
        
    return chi2


def dinf_iter( data_sets, eta_m_guess = 3e-11, d0_guess = 6.0, epsilon = 0.01e-11 ):
    
    # If only one data set is given, put it in a list
    if isinstance(data_sets, list): 
        arguments = data_sets
    else:
        arguments = [data_sets]
        
    # build the list of initial guesses
    init_guess = [eta_m_guess]
    for i in range(len(arguments)):
        init_guess.append(d0_guess)
    
    # Minimize chi^2
    res = sp.optimize.minimize( chi2, init_guess, args = arguments, method = 'Nelder-Mead' )
    # ... and obtain optimized eta_m and D0
    eta_m_opt = res.x[0]
    d0_opt    = res.x[1:]
    
    return eta_m_opt,d0_opt



############################################
#   Definition of the data set structure   #
############################################


class Dataset:
    """
    A data set and all its properties
    """
    
    def __init__(self, name, data_file, temperature = 300, mono = False, label = None, eta_f = 10.2e-04, eta_m_range = [2e-11,7e-11], epsilon = 0.01e-11, b = 2.9e6, approx = False ):
        
        self.name        = name
        self.data        = np.transpose(np.loadtxt( data_file ))
        self.temperature = temperature
        self.mono        = mono
        self.label       = label
        
        self.eta_f       = eta_f
        self.eta_m_range = eta_m_range
        self.epsilon     = epsilon
        self.b           = b
        self.approx      = approx
        
        self.edge        = self.data[0]
        self.height      = (self.data[1]-4.5)/2.0
        self.dc          = self.data[2]*cf
        self.err         = self.data[3]*cf
        
        # Flags to show which analysis has already been performed
        self.calculated_d_inf_iteratively   = False
        
        
    def find_d_inf(self):

        # Calculate D_0 for a range of eta_m at a fixed value of eta_f
        self.eta_m_opt,self.d_inf_opt = dinf_iter( self )
        
        # Calculate SD length from optimized parameters
        self.l_sd_opt  = 0.5*self.eta_m_opt/self.eta_f
        
        # Calculate the single values  of D_0 for each data point again
        if self.approx:
            self.d_inf_points = np.array([mdv.d0_approximation(dpbc,self.temperature,self.eta_f,self.eta_m_opt,self.edge[i],self.height[i]) for i, dpbc in enumerate(self.dc)])
        else:
            if self.mono:
                self.d_inf_points = np.array([mdv.d0_mono(dpbc,self.temperature,self.eta_f,self.eta_m_opt,self.b,self.edge[i],self.height[i],40) for i, dpbc in enumerate(self.dc)])
            else:
                self.d_inf_points = np.array([mdv.d0(dpbc,self.temperature,self.eta_f,self.eta_m_opt,self.edge[i],self.height[i],40) for i, dpbc in enumerate(self.dc)])
        
        # Set the flag for iterative calculation
        self.calculated_d_inf_iteratively = True
        self.calculated_d_inf_by_deviations = False
        
        
    
    # PRINTING AND PLOTTING ROUTINES   

    def overview(self,approx=False):

        if not self.calculated_d_inf_iteratively:
            print("Starting analysis.")
            self.find_d_inf()
            print("Finished analysis.")

        print("D_inf  = %3.3f 10^-7 cm^2/s"  %(self.d_inf_opt*1e7 ))
        print("eta_m  = %3.3f 10^-11 Pa*s*m"%(self.eta_m_opt*1e11))
        print("L_SD   = %3.3f nm"           %(self.l_sd_opt *1e9 ))
    
    
    def newvalues_plot(self,xaxis='edge',save=False,approx=False):

        if not self.calculated_d_inf_iteratively:
            print("Starting analysis.")
            self.find_d_inf()
            print("Finished analysis.")

        if xaxis == 'height':
            xax = self.height
            xlb = r'$H\;[\mathrm{nm}]$'
        else:
            xax = self.edge
            xlb = r'$L_{x,y}\;[\mathrm{nm}]$'

        # Make a plot
        fig, ax = plt.subplots(1, 1, figsize=plt.figaspect(.75)*1.0) #, dpi=300)
        plotsc = 1e7 # scale for the plot

        # Simulation values
        ax.errorbar( xax, self.dc*plotsc, yerr=self.err*plotsc, fmt='o', ls='', label=r'$D_\mathrm{PBC}$', color='C1' )

        # D_inf from the approximation
        ax.errorbar( xax, self.d_inf_points*plotsc, yerr=self.err*plotsc, fmt='o', ls='', label=r'$D_0$', color='C0' )
        ax.axhline( y=np.mean( self.d_inf_points*plotsc ), xmin=0, xmax=50, ls=':', color='C0' )

        # Annotations and Labels
        ax.set(xlabel=xlb)
        ax.set(ylabel=r'$D\;[10^{-7}\;\mathrm{cm}^2/\mathrm{s}]$')
        ax.xaxis.label.set_size(20)
        ax.yaxis.label.set_size(20)
        ax.tick_params(labelsize=14)
        ax.legend(loc='best',numpoints=1,fontsize=20)

        fig.tight_layout()

        if save:
            fig.savefig(self.name+"-dpbc-dinf"+".pdf", format='pdf', dpi=300)
        
