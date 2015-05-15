import string
import matplotlib.mlab as mlab
from scipy.interpolate import interp1d
import numpy as np
import utils
import atmosphere
import properties
import os.path

def regrid(atm,regridType=None,Pmin=None,Pmax=None):
    """This puts atm and cloud on the same grid used later for calculations.  It has three different regimes:
        1 - interpolating within given points:  P,T,z are assumed to follow given adiabat and constituents are linearly interpolated
            NOTE FOR NOW IT JUST DOES THE LINEAR INTERPOLATION
        2 - extrapolating outward:  project last slope out
        3 - interpolating inward:  uses a dry adiabat
       It has three types:
        1 - pressure grid points in file
            format: 'filename' (string)
                where
                    'filename' is a file with the pressures in bars (increasing pressure)
        2 - fixed z step NOT IMPLEMENTED YET
            format:  'step unit'
                where
                    'step' step size (float)
                    'unit' km only right now (string)
        3 - fixed number of steps in logP
            format:  'number'
                    where
                        'number' numbers of steps within range (int)
       Pmin/Pmax are optional for type 2/3 (defaults are min/max in gas) and not used in type 1."""

    #set regridType/regrid
    if regridType==None:
        regridType = atm.config.regridType
    if string.lower(regridType) == 'none' or regridType == None:
        print 'No regridding.  Note that there is a risk that not everything is on the same grid...\n'
        return 0
    if string.lower(regridType) == 'auto' or string.lower(regridType) == 'default':
        regridType = '1000'
    regrid = regridType.split()

    #set default Pmin/Pmax
    if Pmin==None:
        Pmin = min(atm.gas[atm.config.C['P']])
    if Pmax==None:
        Pmax = max(atm.gas[atm.config.C['P']])

    #set Pgrid or zstep(not yet)
    if len(regrid) == 1:
        try:
            regrid_numsteps = int(regrid[0])
            regridType = 3
            Pgrid = np.logspace(np.log10(Pmin),np.log10(Pmax),regrid_numsteps)
            print 'Regridding on %d steps' % (regrid_numsteps)
        except ValueError:
            regrid_file = os.path.join(atm.config.path,regrid[0])
            regridType = 1
            Pgrid = _procF_(regrid_file)
            Pmin = min(Pgrid)
            Pmax = max(Pgrid)
    else:
        regrid_stepsize = float(regrid[0])
        regrid_unit = regrid[1]
        regridType = 2
        Pgrid = None
    if Pgrid is None:
        print 'not implemented yet - no regrid'
        return 0
    
    ### Copy over for new array size
    nAtm = len(Pgrid)
    nGas = atm.gas.shape[0]
    gas = np.zeros((nGas,nAtm))
    nCloud = atm.cloud.shape[0]
    cloud = np.zeros((nCloud,nAtm))

    ### Interpolate gas onto the grid
    berr = False
    interpType = 'linear'
    fillval = -999.9
    Pinput = atm.gas[atm.config.C['P']]
    gas[atm.config.C['P']] = Pgrid
    for yvar in atm.config.C:
        if yvar == 'P':
            continue
        fv = interp1d(Pinput,atm.gas[atm.config.C[yvar]],kind=interpType,fill_value=fillval,bounds_error=berr)
        gas[atm.config.C[yvar]] = fv(Pgrid)
    ### Extrapolate gas if needed
    if fillval in gas:
        gpar = (atm.layerProperty[atm.config.LP['g']][-1],atm.layerProperty[atm.config.LP['R']][-1],atm.layerProperty[atm.config.LP['P']][-1])
        gas = extrapolate(gas,fillval,atm,gpar)
    atm.gas = gas

    ### Interpolate cloud onto the grid - fillval=0.0 extrapolates since we assume no clouds outside range and
    ###      don't care about other stuff then either
    fillval = 0.0
    Pinput = atm.cloud[atm.config.Cl['P']]
    cloud[atm.config.Cl['P']] = Pgrid
    for yvar in atm.config.Cl:
        if yvar == 'P':
            continue
        fv = interp1d(Pinput,atm.cloud[atm.config.Cl[yvar]],kind=interpType,fill_value=fillval,bounds_error=berr)
        cloud[atm.config.Cl[yvar]] = fv(Pgrid)
    atm.cloud = cloud

    #renormalize such that zDeep = 0.0
    zDeep = atm.gas[atm.config.C['Z']][-1]
    if abs(zDeep) > 1E-6:
        for i in range(len(atm.gas[atm.config.C['Z']])):
            atm.gas[atm.config.C['Z']][i]-=zDeep
        print 'Deepest levels set from %.1f to 0' % (zDeep)

    #put in DZ
    dz = np.abs(np.diff(gas[atm.config.C['Z']]))*1.0E5  #convert from km to cm (so no unit checking!!!)
    gas[atm.config.C['DZ']] = np.append(np.array([0.0]),dz)

    atm.computeProp(False)
    print 'Regrid:  %d levels' % (nAtm)
    return 1

def _procF_(filename):
    try:
        reg = open(filename,'r')
    except IOError:
        print "file '"+filename+"' not found - no regrid"
        return 0
    print 'Regridding on file '+filename
    Pgrid = []
    for line in reg:
        Pgrid.append(float(line))
    reg.close()
    Pgrid = np.array(Pgrid)
    if not np.all(np.diff(Pgrid) > 0.0):
        print 'Error in regrid:  P not increasing - flipping around and hoping for the best'
        Pgrid = np.fliplr(Pgrid)
    return Pgrid

def extrapolate(gas,fillval,atm,gpar):
    """First extrapolate in, then the rest of the fillvals get extrapolated out"""
    Pmin = min(atm.gas[atm.config.C['P']])
    Pmax = max(atm.gas[atm.config.C['P']])
    #extrapolate constituents in as fixed mixing ratios
    for yvar in atm.config.C:
        if yvar=='Z' or yvar=='P' or yvar=='T' or yvar=='DZ':
            continue
        val = atm.gas[atm.config.C[yvar]][-1]
        for i,P in enumerate(gas[atm.config.C['P']]):
            if P>Pmax:
                gas[atm.config.C[yvar]][i]=val
    #extrapolate T and z as dry adiabat in hydrostatic equilibrium
    for i,P in enumerate(gas[atm.config.C['P']]):
        if P>Pmax:
            dP = P-gas[atm.config.C['P']][i-1]
            P = gas[atm.config.C['P']][i-1]
            T = gas[atm.config.C['T']][i-1]
            cp = 0.0
            for key in properties.specific_heat:
                if key in atm.config.C:
                    cp+=properties.specific_heat[key]*gas[atm.config.C[key]][i]
            cp*=properties.R #since the catalogued values are Cp/R
            dT = (properties.R*T)/(cp*P)*dP
            gas[atm.config.C['T']][i] = T+dT
            amu = 0.0
            for key in properties.amu:
                if key in atm.config.C:
                    amu+=properties.amu[key]*gas[atm.config.C[key]][i]
            g = gpar[0] + 2.0*properties.R*T*np.log(P/gpar[2])/(gpar[1]*amu)/1000.0
            H = properties.R*T/(amu*g)/1000.0
            dz = H*dP/P
            gas[atm.config.C['Z']][i] = gas[atm.config.C['Z']][i-1]-dz
    return gas            

    #extrapolate out along last slope (move return gas above if you think you want this)
    for yvar in atm.config.C:
        if yvar == 'P':
            continue
        gas[atm.config.C[yvar]] = extrapolateOut(gas[atm.config.C['P']],gas[atm.config.C[yvar]],fillval)
def extrapolateOut(x,y,fillval):
    b = mlab.find(y!=fillval)
    if y[0] == fillval:
        slope = (y[b[0]+1] - y[b[0]])/(x[b[0]+1] - x[b[0]])
        intercept = y[b[0]] - slope*x[b[0]]
        i=0
        while y[i] == fillval:
            y[i] = slope*x[i] + intercept
            i+=1
    if y[-1] == fillval:
        slope = (y[b[-1]] - y[b[-1]-1])/(x[b[-1]] - x[b[-1]-1])
        intercept = y[b[-1]] - slope*x[b[-1]]
        i=-1
        while y[i] == fillval:
            y[i] = slope*x[i] + intercept
            i-=1
    return y
    
