# memdiff
Implementation of corrections for diffusion coefficients in membrane simulations.

## Usage
The package contains three files:
 - memdiff.py: implementation of correction formulas from the point-particle method (1,2) as well as the flat-box approximation (1).
 - memdiff_immersed_boundary.py: implementation of correction formulas from the immersed-boundary method (3).
 - analysis_diffusion.py: a framework to analyse data sets of diffusion coefficients obtained from different box geometries (2).
 
## Literature
 - (1) M. Vögele and G. Hummer: Divergent Diffusion Coefficients in Simulations of Fluids and Lipid Membranes; J. Phys. Chem. B, 2016, 120 (33), pp 8722–8732; DOI: 10.1021/acs.jpcb.6b05102
 - (2) M. Vögele, J. Köfinger, and G. Hummer: Hydrodynamics of Diffusion in Lipid Membrane Simulations, in preparation (February 2018)
 - (3) B. A. Camley, M. G. Lerner, R. W. Pastor., F. L. Brown: Strong influence of periodic boundary conditions on lateral diffusion in lipid bilayer membranes. J Chem Phys. 2015 Dec 28;143(24):243113. doi: 10.1063/1.4932980.
