###This tests the refractive bending math.
import math
import matplotlib.pyplot as plt

raypath = {'egress':-1, 'ingress':1, 'tangent':0}

def compute_dS(b,rv,n,verbose=False,plot=True):
    """Computes the path length through the atmosphere given:
            b = impact parameter (fractional distance to outer edge
            rv = locations of layers.  Next layer is rv+ds
            n = refractive index of layers"""
    b = b*rv[0]
    ##### initialize outside the atmosphere (see pg 92 of notebook)
    # initialize x,y position of ray outside atmosphere
    p = [ [1.001*rv[0]], [b] ]                                      # this is p[0,1][0]
    #  ... now at the atmosphere edge
    p[0].append(math.sqrt(rv[0]**2 - b**2))                         # this is p[0][1]
    p[1].append(b)                                                  # this is p[1][1]

    # initialize s-vector, norm-vector (ds and incAngle)
    s = [ [-1.0], [0.0] ]                                           # this is s[0,1][0]
    normAngle = math.atan2(p[1][1],p[0][1])                         # this is the normal angle at the atmosphere edge
    norm = [ [math.cos(normAngle)], [math.sin(normAngle)] ]         # this is norm[0,1][0]
    incAngle=[normAngle]                                            # this is also an arbitrary value for outside the atmosphere
    ds = []
    dsLayer = []
    ds_store = 0.0
    earlyExit = False

    if verbose:
        print 'b = %f (%f)' % (b,b/rv[0])

    # loop while in atmosphere
    i = 1                      # loops over the path
    layer = 1                  # keeps track which physical layer you are in
    inAtmosphere = True
    direction = 'ingress'
    while inAtmosphere:
        if verbose:
            print '------------------'
        # get refractive index ratio
        try:
            nratio = n[layer-raypath[direction]]/n[layer]
        except:
            inAtmosphere = False
            break

        # get incident angle into i'th layer
        cosThetaInc = -1.0*raypath[direction]*(s[0][i-1]*norm[0][i-1] + s[1][i-1]*norm[1][i-1])   
        try:
            sinThetaInc = math.sqrt(1.0 - cosThetaInc**2)
        except:
            print 'Error in sinThetaInc.  cosThetaInc = ',cosThetaInc
            print 'But we are going to keep going'
            sinThetaInc = 0.0
            cosThetaInc = 1.0
        incAngle.append(math.acos(cosThetaInc))                      # this is incAngle[i]

        # get transmitted angle in i'th layer
        sinThetaTran = nratio*sinThetaInc
        try:
            cosThetaTran = math.sqrt(1.0 - sinThetaTran**2)
        except:
            print 'Error in cosThetaTran.  sinThetaTran = ',sinThetaTran
            print 'But we are going to keep going'
            cosThetaTran = 0.0
            sinThetaTran = 1.0

        # update s-vector in i'th layer
        A = (nratio*cosThetaInc - cosThetaTran)*raypath[direction]
        s[0].append(nratio*s[0][i-1] + A*norm[0][i-1])               # this is s[0][i]
        s[1].append(nratio*s[1][i-1] + A*norm[1][i-1])               # this is s[1][i]
        s_mag = s[0][i]**2 + s[1][i]**2

        r = rv[layer]
        try:
            r_next = rv[layer+raypath[direction]]
        except:
            r_next = r

        # check for tangent layer
        normAngle = math.atan2(norm[1][i-1],norm[0][i-1])
        phi = math.atan2(-s[0][i],s[1][i]) - normAngle
        r_tanlyr = r*math.cos(phi)
        if direction == 'ingress' and r_tanlyr > r_next:
            direction = 'tangent'

        if verbose:
            print 'layer ',layer,direction
            print "x, y, angle, n:  ",p[0][i],p[1][i],incAngle[i]*180.0/math.pi,n[layer]
            if direction != 'egress':
                print 'phi,r,r_next, r_tanlyr  ',180.0*phi/math.pi,r,r_next,r_tanlyr

        # calculate path length in layer
        if direction == 'tangent':
            ds_tan = 2.0*r*math.sin(phi)
            ds.append(ds_tan)
            direction = 'egress'
        else:
            bq = p[0][i]*s[0][i] + p[1][i]*s[1][i]
            try:
                ds.append(-bq - raypath[direction]*math.sqrt(bq**2 + r_next**2 - r**2))      # this is ds[i]
            except:
                print 'Error in (%d, %d).  Appending ds=0.0' % (layer,i)
                print 'This is probably because you are going out through the other side of the atmosphere'
                print 'bq:  ',bq
                print 'r_next^2: ',r_next**2
                print 'r^2:  ',r**2
                print 'taking sqrt of: ',bq**2 + r_next**2 - r**2
                ds.append(0.0)
                inAtmosphere = False
                earlyExit = True

        # write to layer ds that gets returned
        if direction == 'ingress' or direction == 'tangent':
            dsLayer.append(ds[i-1])
            ds_store = ds[i-1]

        # update position
        p[0].append(p[0][i] + ds[i-1]*s[0][i])                         # this is p[0][i+1]
        p[1].append(p[1][i] + ds[i-1]*s[1][i])                         # this is p[1][i+1]

        # update norm-vector
        normAngle = math.atan2(p[1][i+1],p[0][i+1])
        norm[0].append(math.cos(normAngle))                          # this is norm[0][i]
        norm[1].append(math.sin(normAngle))                          # this is norm[1][i]
        norm_mag = norm[0][i]**2 + norm[1][i]**2

        if verbose:
            print 'ds, normAngle:  ',ds[i],normAngle*180.0/math.pi

        # check if still in atmosphere 
        if r_next > rv[0] or r_next < rv[-1]:
            inAtmosphere = False

        #increment counters
        i+=1
        layer+=raypath[direction]

    incAngle.append(incAngle[i-1])  # append to make same number of elements as p    
    if earlyExit:
        earlyExit = False
        for i in range(len(rv) - len(dsLayer)):
            dsLayer.append(0.0)
    else:
        dsLayer.append(ds_store)

    if plot:
        if plot==True:
            plot = range(len(rv))
        elif plot=='limits':
            plot = [1,len(rv)-1]
        plt.figure('raypath')
        plt.plot(p[0],p[1],'-x')
        v = plt.axis()
        circle = []
        for i,r in enumerate(rv):
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


###Test functions
def refractTest(layers):
    """Refurns refractive index at layers"""
    n = [1.0]
    n = []
    for r in layers:
        v = 1.0 + (3500.0/r)**2
        #v = 1.0
        n.append(v)
    #n.append(v)  # add an extra for good measure
    return n
def computeEdges(z):
    edge = [ (3.0*z[0] - z[1])/2.0 ]
    for i in range(len(z)-1):
        edge.append((z[i]+z[i+1])/2.0)
    edge.append( (3.0*z[-1] - z[-2])/2.0 )
    return edge
def layersTest(rmin=100.0,rmax=20000.0,nlyr=100, infinity=1.001):
    """Returns layer mid-points and edges.
       Note that layer 0 is free space"""
    dr = (rmax-rmin)/nlyr
    #print 'rmax,rmin,dr,nlyr: ',rmax,rmin,dr,nlyr
    # set the boundaries
    edges = [infinity*rmax]   # this is the edge at "infinity"
    for i in range(nlyr):
        r = rmax - i*dr
        edges.append(r)
    # set the midpoints (ill-defined for the free space layer 0)
    mid = []
    for i in range(nlyr-1):
        r = (edges[i]+edges[i+1])/2.0
        mid.append(r)
    mid.append(r)           # this is the "interior" (set to the same value as previous mid)
    return mid, edges
def testPath(rmin=100.0,rmax=20000.0,nlyr=100,infinity=1.001,verbose=False,plot=True):
    #make layers
    mid, edge = layersTest(rmin=rmin,rmax=rmax,nlyr=nlyr,infinity=infinity)
    n = refractTest(mid)
    b = 0.22
    ds = compute_dS(b,mid,n,verbose=verbose,plot=plot)
    return ds
