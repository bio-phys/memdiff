import numpy as np
import scipy as sp
import scipy.integrate as integrate


# Implementation of the model by Camley et al.
# usually beta = 0.828494


kB    = 1.3806485*10**-23  # Boltzmann constant [J/K]
m2    = 10**4              # dc [m^2/s] * m2 = dc [cm^2/s]



############################################
## Functions for transmembrane inclusions ##
############################################


def summand(length,height,lsd,br2,ix,iy):

    pre = 2*np.pi/length
    i2  = ix**2 + iy**2
    k2  = i2*pre**2
    k   = pre*np.sqrt(i2)
    
    s   = np.exp(-k2*br2/2)
    s  /= (k2 + k*np.tanh(k*height)/lsd)
    
    return s


def d_pbc(T,eta_f,eta_m,rad,beta,length,height,imax):
    """Diffusion coefficient in a membrane in a periodic box"""
    
    lsd   = eta_m/(2*eta_f)*1.0e9
    br2   = (beta*rad)**2
    
    dpbc  = 0
    for i in xrange(-imax,imax+1):
        for j in xrange(-imax,imax+1):
            if i*i + j*j > 0:
                dpbc += summand(length,height,lsd,br2,i,j)
    
    dpbc /= eta_m*length**2
    
    m2    = 10**4          # dc [m^2/s] * m2 = dc [cm^2/s]
    dpbc *= m2*kB*T/2.0
    
    return dpbc



def d_inf(T,eta_f,eta_m,rad,beta):
    """Diffusion coefficient in an infinite-membrane system"""
    
    lsd   = eta_m/(2*eta_f)*1.0e9
    brl   = beta*rad/(lsd*np.sqrt(2))
    
    dinf  = np.pi*sp.special.erfi(brl) - sp.special.expi(brl**2)
    dinf *= np.exp(-brl**2)/(4*np.pi)
    dinf /= eta_m
    
    m2    = 10**4          # dc [m^2/s] * m2 = dc [cm^2/s]
    dinf *= m2*kB*T/2.0
    
    return dinf


def deltaD(T,eta_f,eta_m,rad,beta,length,height,imax):
    """Difference of finite-size and infinite-size diffusion coefficient [nm^2/2]"""
    
    dpbc = d_pbc(T,eta_f,eta_m,rad,beta,length,height,imax)
    dinf = d_inf(T,eta_f,eta_m,rad,beta)
    
    return dpbc-dinf



########################################
## Functions for monotopic inclusions ##
########################################


def csch(x):
    
    if x > 710:
        res = 0.0
    else:
        res = 1./np.sinh(x)
        
    return res


def summand_mono(length,height,eta_f,eta_m,b,br2,ix,iy):

    lsd   = eta_m/(2*eta_f)*1.0e9
    
    pre = 2*np.pi/length
    i2  = ix**2 + iy**2
    k2  = i2*pre**2
    k   = pre*np.sqrt(i2)
    
    eta_mono = 0.5*eta_m
    aa  = 0.5*k2 + 0.5*k/np.tanh(2*height*k)/lsd + 1.0e-18*b/eta_m
    bb  = 1.0e-18*b/eta_m + 0.5*k*csch(2*height*k)/lsd
    
    s   = np.exp(-k2*br2/2)
    s  *= aa/(aa**2-bb**2)
    
    return s


def d_pbc_mono(T,eta_f,eta_m,b,rad,beta,length,height,imax):
    """Diffusion coefficient in a membrane in a periodic box"""
    
    lsd   = eta_m/(2*eta_f)*1.0e9
    br2   = (beta*rad)**2
    
    dpbc  = 0
    for i in range(-imax,imax+1):
        for j in range(-imax,imax+1):
            if i*i + j*j > 0:
                dpbc += summand_mono(length,height,eta_f,eta_m,b,br2,i,j)
    
    dpbc /= length**2*eta_m
    
    m2    = 10**4          # dc [m^2/s] * m2 = dc [cm^2/s]
    dpbc *= m2*kB*T/2.0
    
    return dpbc


def integrand_mono(k,eta_f,eta_m,b,br2):
    
    k2  = k**2
    lsd = eta_m/(2*eta_f)*1.0e9
    aaa = 0.5*k2 + 0.5*k/lsd + 1.0e-18*b/eta_m
    itg = k*np.exp(-0.5*k2*br2) * aaa/(aaa**2 - (1.0e-18*b/eta_m)**2)
    
    return itg
    

def d_inf_mono(T,eta_f,eta_m,b,rad,beta):
    """Diffusion coefficient in an infinite-membrane system"""
        
    eta_mono = 0.5*eta_m
    br2      = (beta*rad)**2
    I = integrate.quad(integrand_mono, 0.0, np.inf, args=(eta_f,eta_m,b,br2))
    dinf = I[0]/(2*np.pi)
    
    m2    = 10**4          # dc [m^2/s] * m2 = dc [cm^2/s]
    dinf *= m2*kB*T/(2.0*eta_m)
    
    return dinf


def deltaD_mono(T,eta_f,eta_m,b,rad,beta,length,height,imax):
    """Difference of finite-size and infinite-size diffusion coefficient [nm^2/2]"""
    
    dpbc = d_pbc_mono(T,eta_f,eta_m,b,rad,beta,length,height,imax)
    dinf = d_inf_mono(T,eta_f,eta_m,b,rad,beta)
    
    return dpbc-dinf

