import math
import string
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from scipy.interpolate import interp1d
import numpy as np
import utils
import sys
import os
import os.path
import atmosphere as atm

def regrid(atm,regridType=None,Pmin=None,Pmax=None):
    """This puts atm and cloud on the same grid used later for calculations
        regridType is a 2,3 or 4 element string with the following regridding options:
             1:  variable to regrid on ('P' or 'z')
             2:  log or linear regridding ('log' or 'lin') OR filename containing points
             3:  value (either a float or an integer for number of levels) 
             4:  ['unit' - if present assumes that the grid is on values not number of steps]
             e.g. P log 100            ==> 100 layers evenly spaced in log
                  z lin 1 km (default) ==> layers as needed, evenly spaced by 1 km
                  z zpoint.dat         ==> use values in zpoint.dat
        Pmin/Pmax are optional - defaults are min/max in gas.
        Assumes that the gas range encompasses cloud range"""

    ### Determine regridType and parse string
    if regridType==None:
        regridType = atm.regridType
    if string.lower(regridType) == 'none' or regridType == None:
        print 'No regridding.  Note that there is a risk that not everything is on the same grid...\n'
        return 0
    if string.lower(regridType) == 'auto' or string.lower(regridType) == 'default':
        regridType = 'z lin 1 km'
    regrid = regridType.split()
    xvar = string.upper(regrid[0])
    print '---regrid:  '+regridType

    ### Set Pmin/Pmax if needed
    if Pmin==None:
        Pmin = min(atm.gas[atm.C['P']])
    if Pmax==None:
        Pmax = max(atm.gas[atm.C['P']])
    zmin = np.interp(Pmax,atm.gas[atm.C['P']],atm.gas[atm.C['Z']])
    zmax = np.interp(Pmin,atm.gas[atm.C['P']],atm.gas[atm.C['Z']])
    print "Regrid pressure limits:  ",Pmin,Pmax
    print "...corresponding z limits:  ",zmin,zmax
        
    if xvar == 'Z':### flip the array so that z is increasing so that interp works
        atm.gas = np.fliplr(atm.gas)
        atm.cloud = np.fliplr(atm.cloud)
        atm.layerProperty = np.fliplr(atm.layerProperty)
    if not np.all(np.diff(atm.gas[atm.C[xvar]]) > 0.0) or not np.all(np.diff(atm.cloud[atm.Cl[xvar]]) > 0.0):
        print 'Error in regrid:  abscissa not increasing - no regrid.'
        return 0.0
        
    ###xs is the desired abscissa - there are multiple returns embedded if error in regridType
    xs = []
    if len(regrid) == 2:  #read abscissa points from file
        loglin = 'LIN'
        filename = regrid[1]
        filename = os.path.join(atm.path,regrid[1])
        try:
            reg = open(filename,'r')
        except IOError:
            print "file '"+filename+"' not found - no regrid"
            return 0
        print '...interpolating on '+xvar+' from file '+filename
        xs = []
        for line in reg:
            xs.append(float(line))
        reg.close()
        if xs[0] < atm.gas[atm.C['P']][0]:
            print 'lower regrid bound error - will extrapolate'
        if xs[-1] > atm.gas[atm.C['P']][-1]:
            print 'upper regrid bound error - will extrapolate'
    elif len(regrid) == 3:   #compute absicissa and npts
        loglin = string.upper(regrid[1])
        npts = int(regrid[2])
        xs = [0] #BUT REALLY NEED TO COMPUTE AT npts
        print regridType+' not yet supported - no regrid'
        return 0
    elif len(regrid) == 4:
        loglin = string.upper(regrid[1])
        step = float(regrid[2])
        unit = regrid[3]
        if unit in utils.Units:
            step = step*utils.Units[unit]/utils.Units[utils.processingAtmLayerUnit]
        else:
            print "Error in regrid '"+regridType+"' - no regrid"
            return 0
        print '...interpolating on '+xvar+' at '+str(step)+' '+utils.processingAtmLayerUnit+' ('+loglin+')'
        if loglin == 'LIN':
            regridding = True
            v = atm.gas[atm.C[xvar]][0]
            while regridding:
                if xvar=='P':
                    P = v
                else:
                    P = np.interp(v,atm.gas[atm.C['Z']],atm.gas[atm.C['P']])
                    z = v
                if xvar=='P' and P>Pmax:
                    regridding = False
                elif xvar=='Z' and z>zmax:
                    regridding = False
                elif P<Pmin or P>Pmax:
                    pass
                else:
                    xs.append(v)
                v+=step
        elif loglin == 'LOG':
            print regridType+" not yet supported - no regrid"
            print "BUT CAN I BETTER INTEGRATE LOG VERSION INTO LIN VERSION JUST ABOVE?"
            return 0
        else:
            print "Error in regrid '"+regridType+"' - no regrid"
            return 0
    else:
        print "Error in regrid '"+regridType+"' - no regrid"
        return 0
    nAtm = len(xs)
    xs = np.array(xs)

    ### Copy over for new array size
    gas = atm.gas
    nPcs = atm.gas.shape[0]
    atm.gas = np.resize(atm.gas,(nPcs,nAtm))
    cloud = atm.cloud
    nPcs = atm.cloud.shape[0]
    atm.cloud = np.resize(atm.cloud,(nPcs,nAtm))
    layerProperty = atm.layerProperty
    nPcs = atm.layerProperty.shape[0]
    atm.layerProperty = np.resize(atm.layerProperty,(nPcs,nAtm))

    if loglin == 'LOG':
        print 'LOG INTERPOLATION NOT YET SUPPORTED - NO REGRID'
        return 0
    
    berr = False
    interpType = 'linear'
    #interpType = 3
    fillval = -999.9
    print 'interpType = '+interpType
    
    ### Interpolate gas onto the grid
    for yvar in atm.C:
        if yvar == xvar:
            continue
        fv = interp1d(gas[atm.C[xvar]],gas[atm.C[yvar]],kind=interpType,fill_value=fillval,bounds_error=berr)
        atm.gas[atm.C[yvar]] = fv(xs)
        if fillval in atm.gas[atm.C[yvar]]:
            print 'Extrapolating gas '+yvar
            atm.gas[atm.C[yvar]] = extrapolate(xs,atm.gas[atm.C[yvar]],fillval)
        #atm.gas[atm.C[yvar]] = np.interp(xs,gas[atm.C[xvar]],gas[atm.C[yvar]])
    atm.gas[atm.C[xvar]] = xs
    #return 1
    
    ### Interpolate cloud onto the grid
    for yvar in atm.Cl:
        if yvar == xvar:
            continue
        fv = interp1d(cloud[atm.Cl[xvar]],cloud[atm.Cl[yvar]],kind=interpType,fill_value=fillval,bounds_error=berr)
        atm.cloud[atm.Cl[yvar]] = fv(xs)
        if fillval in atm.cloud[atm.Cl[yvar]]:
            print 'Extrapolating cloud '+yvar
            atm.cloud[atm.Cl[yvar]] = extrapolate(xs,atm.cloud[atm.Cl[yvar]],fillval)
        #atm.cloud[atm.Cl[yvar]] = np.interp(xs,cloud[atm.Cl[xvar]],cloud[atm.Cl[yvar]])
    atm.cloud[atm.Cl[xvar]] = xs

    ### Interpolate layerProperties onto the grid
    for yvar in atm.LP:
        if yvar == xvar:
            continue
        fv = interp1d(layerProperty[atm.LP[xvar]],layerProperty[atm.LP[yvar]],kind=interpType,fill_value=fillval,bounds_error=berr)
        atm.layerProperty[atm.LP[yvar]] = fv(xs)
        if fillval in atm.layerProperty[atm.LP[yvar]]:
            print 'Extrapolating layer property '+yvar
            atm.layerProperty[atm.LP[yvar]] = extrapolate(xs,atm.layerProperty[atm.LP[yvar]],fillval)
        #atm.layerProperty[atm.LP[yvar]] = np.interp(xs,layerProperty[atm.LP[xvar]],layerProperty[atm.LP[yvar]])
    atm.layerProperty[atm.LP[xvar]] = xs

    if xvar == 'Z':   ### flip array back
        atm.gas = np.fliplr(atm.gas)
        atm.cloud = np.fliplr(atm.cloud)
        atm.layerProperty = np.fliplr(atm.layerProperty)
    
    return 1


def extrapolate(x,y,fillval):
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
    
