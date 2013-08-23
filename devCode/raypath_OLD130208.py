###This tests the refractive bending math.
import math
import matplotlib.pyplot as plt
import atmosphere
import scipy.special as scisp

raypathdir = {'egress':-1, 'ingress':1, 'tangent':0}
_X = 0
_Y = 1

def compute_ds(b_vec, atm, verbose=False, plot=True):
    """Computes the path length through the atmosphere given:
            b = impact parameter (fractional distance to outer edge
            atmProp = atm.layerProperty """

    b = b_vec[0]
    phi = b_vec[1]  # don't use this yet
    if verbose:
        print 'b = %f (%f)' % (b,b*rl[0])
    rl = atm.layerProperty[atm.LP['R']] # radius of layers
    nr = atm.layerProperty[atm.LP['N']] # refractive index of layers

    ### initialize terms
    rp = [ [rl[0]*math.sqrt(1.0-b*b)], [rl[0]*b] ]                  # this is r_0
    s = [ [-1.0], [0.0] ]                                           # this is s_0
    tinc = [math.asin(b)]                                           # incident angle_0
    nratio = nr[0]/nr[1]                                            # refractive index ratio_0
    ttran = [math.asin(nratio*b)]                                   # transmitted angle_0
    norm = [ [math.cos(tinc[0])],[math.sin(tinc[0])] ]              # n_0
    ds = []
    dsLayer = []
    ds_store = 0.0

    # loop while in atmosphere
    i = 0                      # loops over the path
    layer = 0                  # keeps track which physical layer you are in
    inAtmosphere = True
    direction = 'ingress'
    justPastTangent = False
    while inAtmosphere:
        if verbose:
            print '------------------'

        rnow = rl[layer]
        rnext= rl[layer+raypathdir[direction]]
        sdir = float(raypathdir[direction])
        
        ### update values
        s[_X].append(nratio*s[_X][i] + nratio*math.cos(tinc[i])*norm[_X][i] - math.cos(ttran[i])*norm[_X][i])   # s_X_(i+1)
        s[_Y].append(nratio*s[_Y][i] + nratio*math.cos(tinc[i])*norm[_Y][i] - math.cos(ttran[i])*norm[_Y][i])   # s_Y_(i+1)
        if direction == 'tangent':
            ds_tan = 2.0*rnow*math.sin(phi)
            ds.append(ds_tan)
            dsLayer.append(ds_tan)
            print '\nThis now handles ingress, but not quite tangent and egress or total internal reflection'
        elif justPastTangent:
            ds.append(ds_store)
            justPastTangent = False
        else:
            try:
                ds.append( sdir*( rnow*math.cos(ttran[i]) - math.sqrt(rnext**2.0 - (rnow*math.sin(ttran[i]))**2.0) ) )     # ds_i
            except ValueError:
                ds.append(ds_store)
                inAtmosphere = False
                break
        if direction == 'ingress':
            dsLayer.append(ds[i])
            ds_store = ds[i]
        rp[_X].append(rp[_X][i] + ds[i]*s[_X][i+1])
        rp[_Y].append(rp[_Y][i] + ds[i]*s[_Y][i+1])
        
        norm[_X].append(sdir*rp[_X][i+1]/rnext)
        norm[_Y].append(sdir*rp[_Y][i+1]/rnext)

        tinc.append( math.acos( -1.0*(norm[_X][i+1]*s[_X][i+1] + norm[_Y][i+1]*s[_Y][i+1]) ) )

        if verbose:
            mags = math.sqrt(s[_X][i+1]**2.0 + s[_Y][i+1]**2)
            magn = math.sqrt(norm[_X][i+1]**2.0 + norm[_Y][i+1]**2)
            print '\tstep %d:  layer %d %s ' % (i,layer,direction)
            print '\ts:  (%f,%f) -> (%f,%f)  |%f|  ds: %f' % (s[_X][i],s[_Y][i],s[_X][i+1],s[_Y][i+1],mags,ds[i])
            print '\tnorm:  (%f,%f) -> (%f,%f) |%f|' % (norm[_X][i],norm[_Y][i],norm[_X][i+1],norm[_Y][i+1],magn)
            print '\trp: (%f,%f) -> (%f,%f)' % (rp[_X][i],rp[_Y][i],rp[_X][i+1],rp[_Y][i+1])
            print '\t|rp|: %f -> %f' % (math.sqrt(rp[_X][i]**2 + rp[_Y][i]**2),math.sqrt(rp[_X][i+1]**2+rp[_Y][i+1]**2))
            print '\tinc, tran, tinc:  %f -> %f -> %f' % (180.0*tinc[i]/math.pi,180.0*ttran[i]/math.pi,180.0*tinc[i+1]/math.pi)
            print '\tr:  %f -> %f' % (rnow,rnext)

        # prep for next loop and check for limits/tangent
        if direction == 'tangent':
            direction= 'egress'
            layer += 1
            justPastTangent = True
        i+=1
        layer+=raypathdir[direction]
        # get refractive index ratio
        try:
            nratio = nr[layer]/nr[layer+raypathdir[direction]]
        except IndexError:
            inAtmosphere = False
            break
        try:
            ttran.append(math.asin(nratio*math.sin(tinc[i])))
        except ValueError:
            print 'total internal reflection'
            print 'fix this later...'
            ttran.append(0.0)
            inAtmosphere = False
        if verbose:
            print '\t\ttran -> %f' % (180.0*ttran[i]/math.pi)

        rnow = rl[layer]                        # These two get duplicated at the start of the loop to handle the tangent layer
        rnext= rl[layer+raypathdir[direction]]  #   

        # check for tangent layer
        normAngle = math.atan2(norm[_Y][i],norm[_X][i])
        phi = math.atan2(-s[_X][i],s[_Y][i]) - normAngle
        rtanlyr = rnow*math.cos(phi)
        if direction == 'ingress' and rtanlyr > rnext:
            direction = 'tangent'
            print 'I think this finds it too early...'

        # check if still in atmosphere
        if rnext > rl[0] or rnext < rl[-1]:
            inAtmosphere = False

        ### end loop ###

    # need to round out the 'kept' ds's (need one extra for 'normal' end or pad out for early exit)
    for i in range(len(rl) - len(dsLayer)):
            dsLayer.append(ds_store)

    if plot:
        if plot==True:
            plot = range(len(rl))
            lt = '-x'
        elif plot=='limits':
            plot = [0,len(rl)-1]
            lt = '-'
        plt.figure('raypath')
        plt.plot(rp[_X],rp[_Y],lt)
        v = plt.axis()
        circle = []
        for i,r in enumerate(rl):
            circle.append([[],[]])
            for th in range(361):
                angle = (th-30.0)*math.pi/180.0
                circle[i][0].append(r*math.cos(angle))
                circle[i][1].append(r*math.sin(angle))
            if i in plot: #i==1 or i==len(rv)-1: #i>0:
                plt.plot(circle[i][0],circle[i][1],color='black')
        plt.axis('image')
        #plt.axis([v[0],v[1],v[2]*0.99,v[3]*1.001])

    return dsLayer

def computeG(b_vec, atm):
    
    ratio = 1
    return ratio

###Test functions
def refractTest(layers):
    """Refurns refractive index at layers"""
    n = [1.0]
    n = []
    for r in layers:
        v = 1.0 + (3000.0/r)**2
        #v = 1.0
        n.append(v)
    return n
def computeEdges(z): # not used anymore
    edge = [ (3.0*z[0] - z[1])/2.0 ]
    for i in range(len(z)-1):
        edge.append((z[i]+z[i+1])/2.0)
    edge.append( (3.0*z[-1] - z[-2])/2.0 )
    return edge
def layersTest(rmin=100.0,rmax=20000.0,nlyr=100):
    """Returns layer"""
    dr = (rmax-rmin)/nlyr
    mid = []
    for i in range(nlyr):
        r = rmax - i*dr
        mid.append(r)
    return mid
def testPath(b=0.5,rmin=12000.0,rmax=20000.0,nlyr=50,verbose=False,plot=True):
    #make layers
    mid = layersTest(rmin=rmin,rmax=rmax,nlyr=nlyr)
    n = refractTest(mid)
    ds = compute_ds(b,mid,n,verbose=verbose,plot=plot)
    plt.figure('ds')
    plt.plot(mid,ds)
    return mid,ds
