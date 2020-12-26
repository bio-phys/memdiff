import numpy as np
import scipy as sp
import scipy.integrate as integrate


# Implementation of the model by Voegele and Hummer
# If you use this, please cite:


kB    = 1.3806485*10**-23  # Boltzmann constant [J/K]
m2    = 10**4              # dc [m^2/s] * m2 = dc [cm^2/s]



############################################
## Functions for transmembrane inclusions ##
############################################


def integral(sigma,lsd):
    
    u = 1/(np.sqrt(2)*sigma*lsd)
    i = np.pi*sp.special.erfi(u) - sp.special.expi(u**2)
    i *= np.exp(-u**2)/(4*np.pi)
    
    return i



def summand(length,height,ix,iy,sigma,lsd):
    
    # vectors in Fourier space
    kx   = 2*np.pi*ix/length
    ky   = 2*np.pi*iy/length
    k2   = kx**2+ky**2
    k    = np.sqrt(k2)
    
    # Calculate the summand
    s    = 1/(k2 + k*np.tanh(k*height)/lsd) + (np.exp(-k2/(2*sigma**2))-1)/(k2+k/lsd)

    return s



def deltaT(length,height,imax,lsd):
    """ Returns Delta_T/eta_m as a function of H and L
        Parameters:
        length: L [nm]
        height: H [nm]
        imax:   maximal index of k-space vectors over which to calculate the sum
    """
    
    # Values for two succeeding values of sigma are calculated to extrapolate
    deltaT_sigma  = []
    nsigma = np.array([5,6])
    sigma  = nsigma*2*np.pi/length
    
    # Loop over both values of sigma
    for sig in sigma:

        # Calculate the contribution of the sum
        deltaT_sum = 0.0
        for ix in xrange(-imax,imax+1):
            for iy in xrange(-imax,imax+1):
                ii = ix**2 + iy**2
                if not ii == 0 and ii <= imax**2:
                    deltaT_sum += summand(length,height,ix,iy,sig,lsd)
        deltaT_sum /= length**2
        # Calculate the contribution of the integral
        deltaT_int = integral(sig,lsd)
        # Calculate deltaT for the respective sigma
        deltaT_sigma.append(deltaT_sum - deltaT_int)
    
    # Extrapolate to the limit of sigma -> infinity
    deltaT_extrapolated  = nsigma[1]**2*deltaT_sigma[1] - nsigma[0]**2*deltaT_sigma[0]
    deltaT_extrapolated /= (nsigma[1]**2 - nsigma[0]**2)
    
    return deltaT_extrapolated/2.



def deltaD(T,eta_f,eta_m,l,h,imax):
    """Correction for a membrane diffusion coefficient in a periodic box (in cm^2/s)"""
    
    lsd = eta_m/(2*eta_f)*1e09
    dT  = deltaT(l,h,imax,lsd)
    
    return m2*dT*kB*T/eta_m



def d0(dpbc,T,eta_f,eta_m,l,h,imax):
    """Corrected membrane diffusion coefficient in a periodic box (input and output in cm^2/s, dimensions in nm)"""
    
    d_corr = dpbc - deltaD(T,eta_f,eta_m,l,h,imax) 
    
    return d_corr



############################################
## Functions for approximations           ##
## (derived for transmembrane components, ##
## applicable also for monotopic ones     ##
## in the appropriate regime, see paper)  ##
############################################


def deltaD_approximation(T,eta_f,eta_m,l,h):
    """Correction from the approximation fomula by Voegele and Hummer"""
    
    kbt  = kB*T
    lsd  = 1.0e09*eta_m/(2*eta_f)
    xi   = np.log( 1 + (np.pi/2)*h/lsd ) + np.e -1
    corr = kbt/(4*np.pi/m2*eta_m) * ( np.log(l/lsd) - xi )/(1+h/lsd)
    
    return corr


def d0_approximation(dpbc,T,eta_f,eta_m,l,h):
    """Infinite-system diffusion coefficient corrected from the finite-size value by the approximation fomula by Voegele and Hummer"""
    
    d0  = dpbc - deltaD_approximation(T,eta_f,eta_m,l,h)
    
    return d0



########################################
## Functions for monotopic inclusions ##
########################################


def csch(x):
    
    if x > 710:
        res = 0.0
    else:
        res = 1./np.sinh(x)
        
    return res


def integrand_mono(k,lsd,bstar,sigma):
    
    k2  = k**2
    aaa = 0.5*k2 + 0.5*k/lsd + bstar
    itg = k*np.exp(-0.5*k2/sigma**2) * aaa/(aaa**2 - bstar**2)
    
    return itg


def integral_mono(sigma,lsd,bstar):
    
    I = integrate.quad(integrand_mono, 0.0, np.inf, args=(lsd,bstar,sigma))
    
    return I[0]/(2*np.pi)


def summand_mono(length,height,ix,iy,sigma,lsd,bstar):
    
    # vectors in Fourier space
    kx   = 2*np.pi*ix/length
    ky   = 2*np.pi*iy/length
    k2   = kx**2+ky**2
    k    = np.sqrt(k2)
    
    # help functions
    aa   = 0.5*k2 + 0.5*k/np.tanh(2*height*k)/lsd + bstar
    aaa  = 0.5*k2 + 0.5*k/lsd + bstar
    bb   = 0.5*k*csch(2*height*k)/lsd + bstar
    
    # First summand
    s    = aa/(aa**2-bb**2)
    
    # Second summand
    s   += (np.exp(-0.5*k2/sigma**2)-1)*aaa/(aaa**2-bstar**2)

    return s



def deltaT_mono(length,height,imax,lsd,bstar):
    """ Returns Delta_T/eta_m as a function of H and L
        Parameters:
        length: L [nm]
        height: H [nm]
        imax:   maximal index of k-space vectors over which to calculate the sum
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
                    deltaT_sum += summand_mono(length,height,ix,iy,sig,lsd,bstar)
        deltaT_sum /= length**2
        # Calculate the contribution of the integral
        deltaT_int = integral_mono(sig,lsd,bstar)
        # Calculate deltaT for the respective sigma
        deltaT_sigma.append(deltaT_sum - deltaT_int)
    
    # Extrapolate to the limit of sigma -> infinity
    deltaT_extrapolated  = nsigma[1]**2*deltaT_sigma[1] - nsigma[0]**2*deltaT_sigma[0]
    deltaT_extrapolated /= (nsigma[1]**2 - nsigma[0]**2)
    
    return deltaT_extrapolated/2.



def deltaD_mono(T,eta_f,eta_m,b,l,h,imax):
    """Correction for a membrane diffusion coefficient in a periodic box (in cm^2/s)"""
    
    lsd   = eta_m/(2*eta_f)*1e09
    bstar = 1.0e-18*b/eta_m
    dT    = deltaT_mono(l,h,imax,lsd,bstar)
    
    return m2*dT*kB*T/eta_m



def d0_mono(dpbc,T,eta_f,eta_m,b,l,h,imax):
    """Corrected membrane diffusion coefficient in a periodic box (input and output in cm^2/s, dimensions in nm)"""
    
    d_corr = dpbc - deltaD_mono(T,eta_f,eta_m,b,l,h,imax) 
    
    return d_corr

