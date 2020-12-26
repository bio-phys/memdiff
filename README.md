![flow chart](images/MemDiff-Picture.png)

Lateral diffusion cefficients in membrane simulations with periodic boundary conditions are subject to substantial hydrodynamic finite-size effects.
You can use the implementation of the correction formulas here to obtain the unperturbed value from the values calculated in the simulation.

## Installation

If you want to keep MemDiff separate from your other Python code, create a new conda environment.

    conda create --name memdiff python=3.7 pip numpy scipy matplotlib jupyter
    conda activate memdiff
 
Clone the repository and install MemDiff via pip.

    git clone https://github.com/bio-phys/memdiff.git
    cd memdiff
    pip install -e . 


## Usage

The package contains three modules:
 - _oseen_: implementation of the Oseen corrections for transmembrane and monotopic components (1,2) as well as the flat-box approximation (1).
 - _immersed_: implementation of correction formulas from the immersed-boundary method (3).
 - _analysis_: a framework to analyse data sets of diffusion coefficients obtained from different box geometries (2).
 
Importing memdiff to your Python code will make all of them available

    import memdiff

To calculate, for example, the corrected diffusion coefficient according to the numerical Oseen correction (1,2), call

    memdiff.oseen.d0(dpbc,T,eta_f,eta_m,l,h,imax)
    
with the following parameters:

    dpbc:  uncorrected diffusion coefficient in cm^2/s, 
    T:     temperature in K
    eta_f: solvent viscosity in Pa*s
    eta_m: membrane surface viscosity in Pa*s*m
    l:     width of the simulation box in nm
    h:     height of the simulation box in nm
    imax:  maximal index of k-space vectors to take into account. 

The higher imax, the more precise and the slower the calculation. 20 is usually a good compromise.
    
There are analogous functions for the approximation formula and for the monotopic correction. See the chart below to decide which one to use.

An example Jupyter notebook is provided that uses data from a simulation of a carbon nanotube porin in a POPC/DOPC membrane to show the usage of the more advanced analysis functions. It shows how to fit the diffusion coefficient and the membrane viscosity from a series of simulations at different values of the box width.

## Which correction formula should I use?

![flow chart](images/membrane-diffusion-flowchart.png)

A PDF version of the flow chart is available ![here](images/membrane-diffusion-flowchart.pdf).

Details and tests are described in our papers (1,2).

## Literature
 - (1) M. Vögele and G. Hummer: Divergent Diffusion Coefficients in Simulations of Fluids and Lipid Membranes; J. Phys. Chem. B, 2016, 120 (33), pp 8722–8732; DOI: 10.1021/acs.jpcb.6b05102
 - (2) M. Vögele, J. Köfinger, and G. Hummer: Hydrodynamics of Diffusion in Lipid Membrane Simulations, Phys. Rev. Lett. 2018, 120, 268104. DOI: 10.1103/PhysRevLett.120.268104 (preprint available at arXiv:1803.04714)
 - (3) B. A. Camley, M. G. Lerner, R. W. Pastor., F. L. Brown: Strong influence of periodic boundary conditions on lateral diffusion in lipid bilayer membranes. J Chem Phys. 2015 Dec 28;143(24):243113. doi: 10.1063/1.4932980.
