#---this one analyze the profile for using sage:
import numpy as np
var('r','r0','A','B')


def s2profile(x):
    """
    hyperbolic secant square function
    """
    res = x/(np.exp(x)+np.exp(-x))**2
    return res

def gprofile(x):
    """
    Fit the binned distribution to a 1D gaussian profile with a constant
    """
    res = x*np.exp(-0.5*(x)**2)
    return res

def mprofile(x):
    """
    Fit the light distribution to a Moffat profile
    """
    res = x*(1+x**2)**(-beta)
    return res


integral(s2profile(x),x,0,infinity)

forget()
var('x','beta')
assume(beta != 1)
integral(mprofile(x),x)
integral(gprofile(x),x,0,infinity).n()
integrate(gprofile(x),x,0,infinity)
