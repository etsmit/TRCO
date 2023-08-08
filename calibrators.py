"""
Coefficients to compute radio SED of calibrator sources.
"""

import numpy as np

from astropy import units as u

from groundhog.fluxscales import bscales


def poly_pb(nu, coefs):
    """
    Equation
    """
    
    nu = nu.to(u.GHz).value
    
    # Calibrator flux density in Jy.
    snu = np.power(10., np.polyval(coefs, np.log10(nu)))*u.Jy
    
    return snu


def poly_ott(nu, coefs):
    """
    From Table 5 in Ott et al. 1994.
    """
    
    nu = nu.to(u.MHz).value
    lognu = np.log10(nu)

    snu = np.power(10., ( coefs[0] + coefs[1]*lognu + coefs[2]*np.power(lognu, 2.) ))*u.Jy
    
    return snu
    

def compute_sed(freq, scale, source, units='Jy'):
    """
    """
    
    coefs = cal_coefs[scale][source]
    #nu = freq.to(u.GHz).value
    
    # Calibrator flux density in Jy.
    snu = cal_coefs[scale]['method'](freq, coefs)
    #snu = np.power(10., np.polyval(coefs, np.log10(nu)))*u.Jy
    
    if 'K' in units:
        conv = bscales.jy2k(freq)
        snu *= conv
    
    return snu.to(units)


calibrators = ['3C48', '3C123', '3C138', '3C147', '3C196', '3C286', '3C295', '3C348', '3C353', '3C380']


cal_coefs = {'Perley-Butler 2017':{'method' : poly_pb,
                                   '3C48' : [0.04980, -0.1914, -0.7553,  1.3253],
                                   '3C123': [0.00900, -0.0248, -0.1035, -0.7884,  1.8017],
                                   '3C138': [0.02230, -0.0102, -0.1552, -0.4981,  1.0088],
                                   '3C147': [0.02890, -0.0464,  0.0640, -0.2007, -0.6961, 1.4516],
                                   '3C196': [0.02010, -0.0200, -0.1534, -0.8530,  1.2872],
                                   '3C286': [0.03570, -0.1798, -0.4507,  1.2481],
                                   '3C295': [0.03990, -0.0347, -0.2780, -0.7658,  1.4701],
                                   '3C348': [0.00000, -0.0951, -1.0247,  1.8298],
                                   '3C353': [-0.0732, -0.0998, -0.6938,  1.8627],
                                   '3C380': [-0.1566, -0.1794,  0.0976,  0.0947, -0.7909, 1.2320],
                                   }, # Perley and Butler 2017
             'Ott 1994' : {'method': poly_ott,
                           '3C286': [0.956,  0.584, -0.1644],
                           '3C353': [3.148, -0.157, -0.0911],
                           }, # Ott et al. 1994
            } 
