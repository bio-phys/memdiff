# -------------------------------------------
#  Implementation of the Oseen correction
# -------------------------------------------
#
#  References:
# -------------
#
#  M. Vögele and G. Hummer: 
#  Divergent Diffusion Coefficients in Simulations of Fluids and Lipid Membranes; 
#  J. Phys. Chem. B, 2016, 120 (33), pp 87228732; DOI: 10.1021/acs.jpcb.6b05102
# 
#  M. Vögele, J. Köfinger, and G. Hummer: Hydrodynamics of Diffusion in Lipid Membrane Simulations;
#  Phys. Rev. Lett. 2018, 120, 268104. DOI: 10.1103/PhysRevLett.120.268104 
#  (preprint available at arXiv:1803.04714)
#

import numpy as np
import scipy as sp
import scipy.integrate as integrate

kB    = 1.3806485*10**-23  # Boltzmann constant [J/K]
m2    = 10**4              # dc [m^2/s] * m2 = dc [cm^2/s]


# -- Functions for transmembrane inclusions --


def _integral(sigma,lsd):
    u = 1/(np.sqrt(2)*sigma*lsd)
    i = np.pi*sp.special.erfi(u) - sp.special.expi(u**2)
    i *= np.exp(-u**2)/(4*np.pi)
    return i


def _summand(length,height,ix,iy,sigma,lsd):
    # vectors in Fourier space
    kx   = 2*np.pi*ix/length
    ky   = 2*np.pi*iy/length
    k2   = kx**2+ky**2
    k    = np.sqrt(k2)    
    # Calculate the summand
    s    = 1/(k2 + k*np.tanh(k*height)/lsd) + (np.exp(-k2/(2*sigma**2))-1)/(k2+k/lsd)
    return s


def deltaT(length,height,imax,lsd):
    """ 
    Calculates Delta_T as a function of H and L.
    
    Args:
        length: box width [nm]
        height: box height [nm]
        imax:   maximal index of k-space vectors over which to calculate the sum
        lsd:    Saffman-Delbrück length [nm]
        
    Returns:
        Delta_T
    """
    # Values for two succeeding values of sigma are calculated to extrapolate
    deltaT_sigma  = []
    nsigma = np.array([5,6])
    sigma  = nsigma*2*np.pi/length
    # Loop over both values of sigma
    for sig in sigma:
        # Calculate the contribution of the sum
        deltaT_sum = 0.0
        for ix in range(-imax,imax+1):
            for iy in range(-imax,imax+1):
                ii = ix**2 + iy**2
                if not ii == 0 and ii <= imax**2:
                    deltaT_sum += _summand(length,height,ix,iy,sig,lsd)
        deltaT_sum /= length**2
        # Calculate the contribution of the integral
        deltaT_int = _integral(sig,lsd)
        # Calculate deltaT for the respective sigma
        deltaT_sigma.append(deltaT_sum - deltaT_int)
    # Extrapolate to the limit of sigma -> infinity
    deltaT_extrapolated  = nsigma[1]**2*deltaT_sigma[1] - nsigma[0]**2*deltaT_sigma[0]
    deltaT_extrapolated /= (nsigma[1]**2 - nsigma[0]**2)
    return deltaT_extrapolated/2.


def deltaD(T,eta_f,eta_m,l,h,imax):
    """
    Correction for a membrane diffusion coefficient in a periodic box,
    calculated from the numerical Oseen method.
    
    Args:
        T:     temperature [K]
        eta_f: solvent viscosity [Pa*s]
        eta_m: membrane surface viscosity [Pa*s*m]
        l:     box width [nm]
        h:     box height [nm]
        imax:  largest index of k-space vectors over which to calculate the sum
        
    Returns:
        deltaD: diffusion correction [cm^2/s]
    """
    lsd = eta_m/(2*eta_f)*1e09
    dT  = deltaT(l,h,imax,lsd)
    dD  = m2*dT*kB*T/eta_m
    return dD


def d0(dpbc,T,eta_f,eta_m,l,h,imax):
    """
    Corrected membrane diffusion coefficient in a periodic box, 
    calculated from the numerical Oseen method.
    
    Args:
        dpbc:  uncorrected diffusion coefficient [cm^2/s]
        T:     temperature [K]
        eta_f: solvent viscosity [Pa*s]
        eta_m: membrane surface viscosity [Pa*s*m]
        l:     box width [nm]
        h:     box height [nm]
        imax:  largest index of k-space vectors over which to calculate the sum
        
    Returns:
        d0:    corrected diffusion coefficient [cm^2/s]
    """
    d_corr = dpbc - deltaD(T,eta_f,eta_m,l,h,imax) 
    return d_corr


# -- Functions for the approximation --
#
#    The approximate correction was derived for transmembrane components, 
#    but is applicable also for monotopic ones in the appropriate regime.
#    See Vögele et al., Phys. Rev. Lett. 2018, for details.
#


def deltaD_approximation(T,eta_f,eta_m,l,h):
    """
    Correction for a membrane diffusion coefficient in a periodic box,
    calculated from the approximation fomula by Vögele and Hummer.
    
    Args:
        T:     temperature [K]
        eta_f: solvent viscosity [Pa*s]
        eta_m: membrane surface viscosity [Pa*s*m]
        l:     box width [nm]
        h:     box height [nm]
        
    Returns:
        deltaD: diffusion correction [cm^2/s]
    """
    kbt  = kB*T
    lsd  = 1.0e09*eta_m/(2*eta_f)
    xi   = np.log( 1 + (np.pi/2)*h/lsd ) + np.e -1
    corr = kbt/(4*np.pi/m2*eta_m) * ( np.log(l/lsd) - xi )/(1+h/lsd)
    return corr


def d0_approximation(dpbc,T,eta_f,eta_m,l,h):
    """
    Corrected membrane diffusion coefficient in a periodic box, 
    calculated from the approximation fomula by Vögele and Hummer.
    
    Args:
        dpbc:  uncorrected diffusion coefficient [cm^2/s]
        T:     temperature [K]
        eta_f: solvent viscosity [Pa*s]
        eta_m: membrane surface viscosity [Pa*s*m]
        l:     box width [nm]
        h:     box height [nm]
        
    Returns:
        d0:    corrected diffusion coefficient [cm^2/s]
    """
    d0  = dpbc - deltaD_approximation(T,eta_f,eta_m,l,h)
    return d0


# -- Functions for monotopic inclusions --
#
#    See Vögele et al., Phys. Rev. Lett. 2018, for details.
#

def _csch(x):
    if x > 710:
        res = 0.0
    else:
        res = 1./np.sinh(x)
    return res


def _integrand_mono(k,lsd,bstar,sigma):
    k2  = k**2
    aaa = 0.5*k2 + 0.5*k/lsd + bstar
    itg = k*np.exp(-0.5*k2/sigma**2) * aaa/(aaa**2 - bstar**2)
    return itg


def _integral_mono(sigma,lsd,bstar):
    I = integrate.quad(_integrand_mono, 0.0, np.inf, args=(lsd,bstar,sigma))
    return I[0]/(2*np.pi)


def _summand_mono(length,height,ix,iy,sigma,lsd,bstar):
    # vectors in Fourier space
    kx   = 2*np.pi*ix/length
    ky   = 2*np.pi*iy/length
    k2   = kx**2+ky**2
    k    = np.sqrt(k2)
    # help functions
    aa   = 0.5*k2 + 0.5*k/np.tanh(2*height*k)/lsd + bstar
    aaa  = 0.5*k2 + 0.5*k/lsd + bstar
    bb   = 0.5*k*_csch(2*height*k)/lsd + bstar
    # First summand
    s    = aa/(aa**2-bb**2)
    # Second summand
    s   += (np.exp(-0.5*k2/sigma**2)-1)*aaa/(aaa**2-bstar**2)
    return s


def deltaT_mono(length,height,imax,lsd,bstar):
    """ 
    Calculates Delta_T/eta_m for the monotopic correction as a function of H and L.
    
    Args:
        length: box width [nm]
        height: box height [nm]
        imax:   maximal index of k-space vectors over which to calculate the sum
        lsd:    Saffman-Delbrück length [nm]
        bstar:  normalized inter-leaflet friction parameter
        
    Returns:
        Delta_T/eta_m for the monotopic correction.
    """
    # Values for two succeeding values of sigma are calculated to extrapolate
    deltaT_sigma  = []
    nsigma = np.array([5,6])
    sigma  = nsigma*2*np.pi/length
    # Loop over both values of sigma
    for sig in sigma:
        # Calculate the contribution of the sum
        deltaT_sum = 0.0
        for ix in range(-imax,imax+1):
            for iy in range(-imax,imax+1):
                ii = ix**2 + iy**2
                if not ii == 0 and ii <= imax**2:
                    deltaT_sum += _summand_mono(length,height,ix,iy,sig,lsd,bstar)
        deltaT_sum /= length**2
        # Calculate the contribution of the integral
        deltaT_int = _integral_mono(sig,lsd,bstar)
        # Calculate deltaT for the respective sigma
        deltaT_sigma.append(deltaT_sum - deltaT_int)
    # Extrapolate to the limit of sigma -> infinity
    deltaT_extrapolated  = nsigma[1]**2*deltaT_sigma[1] - nsigma[0]**2*deltaT_sigma[0]
    deltaT_extrapolated /= (nsigma[1]**2 - nsigma[0]**2)
    return deltaT_extrapolated/2.


def deltaD_mono(T,eta_f,eta_m,b,l,h,imax):
    """
    Correction for a membrane diffusion coefficient in a periodic box,
    assuming a monotopic membrane inclusion.
    
    Args:
        T:     temperature [K]
        eta_f: solvent viscosity [Pa*s]
        eta_m: membrane surface viscosity [Pa*s*m]
        b:     inter-leaflet friction [Pa*s/m]
        l:     box width [nm]
        h:     box height [nm]
        imax:  largest index of k-space vectors over which to calculate the sum
        
    Returns:
        deltaD: diffusion correction [cm^2/s]
    """
    lsd   = eta_m/(2*eta_f)*1e09
    bstar = 1.0e-18*b/eta_m
    dT    = deltaT_mono(l,h,imax,lsd,bstar)
    return m2*dT*kB*T/eta_m


def d0_mono(dpbc,T,eta_f,eta_m,b,l,h,imax):
    """
    Corrected membrane diffusion coefficient in a periodic box, 
    assuming a monotopic membrane inclusion.
    
    Args:
        dpbc:  uncorrected diffusion coefficient [cm^2/s]
        T:     temperature [K]
        eta_f: solvent viscosity [Pa*s]
        eta_m: membrane surface viscosity [Pa*s*m]
        b:     inter-leaflet friction [Pa*s/m]
        l:     box width [nm]
        h:     box height [nm]
        imax:  largest index of k-space vectors over which to calculate the sum
        
    Returns:
        d0:    corrected diffusion coefficient [cm^2/s]
    """
    d_corr = dpbc - deltaD_mono(T,eta_f,eta_m,b,l,h,imax) 
    return d_corr

