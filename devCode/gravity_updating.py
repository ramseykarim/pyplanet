import numpy as np
import matplotlib.pyplot as plt
import scipy.special as scisp
import math
import atmosphere

H = {'psi':0, 'lat':1, 'r':2, 'n':3, 't':4}

def calcGeoid(planet, r, phi):
    latstep = sign(phi)*0.05
    pclat = np.arange(0.0,phi+latstep,latstep)
    for ii,v in enumerate(pclat):
        radv = v*np.pi/180.0
        vw = np.interp(v,planet.vwlat,planet.vwdat)/1000.0
        omega = planet.omega_m + vw/(r*np.cos(radv))
        g, geoid = gravity(v, r, planet.RJ, planet.Req, planet.Rpol, planet.GM, omega, planet.Jn)
        dlat = latstep*np.pi/180.0
        dr_vec = r*dlat*np.array([geoid[H['t']][0], geoid[H['t']][1]])
        r_vec = geoid[H['r']]+dr_vec
        r = math.sqrt(r_vec[0]**2 + r_vec[1]**2)
    return geoid
def gravity(pclat, r, RJ, Req, Rpol, GM, omega, Jn):
    g_static = GM/r**2
    lat = pclat*np.pi/180.0
    dphi = 0.00001
    Sr = 0.0
    Sp = 0.0
    sp = np.sin(lat)
    sp1 = np.sin(lat + dphi)
    sp0 = np.sin(lat - dphi)
    for i in range(len(Jn)):
        # g(r)
        Sr+= (i+1.0)*Jn[i]*pow( RJ/r, i )*scisp.legendre(i)(sp)
        # g(phi)
        dP = (scisp.legendre(i)(sp1) - scisp.legendre(i)(sp))/dphi
        dP+= (scisp.legendre(i)(sp) - scisp.legendre(i)(sp0))/dphi
        dP*= 0.5
        Sp+= Jn[i]*pow( RJ/r, i )*dP
    # g(r)
    gr = g_static*(1.0 - Sr) - (2.0/3.0)*(omega**2.0)*r*(1.0 - scisp.legendre(2)(sp))
    # g(phi)
    dP = (3.0*sp*np.sqrt(1.0-sp**2))
    gp = (1.0/3.0)*(omega**2.0)*r*dP + g_static*Sp
    gt = math.sqrt(gr**2 + gp**2)
    g = [gr,gp,gt]
    # geoid
    psi = math.atan2(gp, gr)
    norm = [math.cos(lat+psi), math.sin(lat+psi)]
    tang = [-math.sin(lat + psi), math.cos(lat + psi)]
    r_vec = [r*math.cos(lat), r*math.sin(lat)]
    geoid = [psi, lat, r_vec, norm, tang]
    return g, geoid

