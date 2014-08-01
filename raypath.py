###Implements full 3-D refractive path with geoid gravity
import math
import matplotlib.pyplot as plt
import atmosphere
import numpy as np
import sys
USE_MAYAVI = False
if USE_MAYAVI:
    import mayavi.mlab as mlab
import shape

raypathdir = {'egress':-1, 'ingress':1, 'tangent':0}
# Vectors
_X = 0
_Y = 1
_Z = 2
xHat = np.array([1.0,0.0,0.0])
yHat = np.array([0.0,1.0,0.0])
zHat = np.array([0.0,0.0,1.0])
plotExist = False

class Ray:
    def __init__(self,ds=None,layer4ds=None,r4ds=None,P4ds=None,doppler=None,tip=None,rotate=None,rNorm=None):
        self.ds = ds
        self.layer4ds = layer4ds
        self.r4ds = r4ds
        self.P4ds = P4ds
        self.rNorm = rNorm
        self.tip = tip
        self.rotate = rotate
        self.doppler = doppler
    def update(self,ds=None,layer4ds=None,r4ds=None,P4ds=None,doppler=None,tip=None,rotate=None,rNorm=None):
        if ds != None:
            self.ds = ds
        if layer4ds != None:
            self.layer4ds = layer4ds
        if r4ds != None:
            self.r4ds = r4ds
        if P4ds != None:
            self.P4ds = P4ds
        if rNorm != None:
            self.rNorm = rNorm
        if tip != None:
            self.tip = tip
        if rotate != None:
            self.rotate = rotate
        if doppler != None:
            self.doppler = doppler
def __computeAspect__(Q,f=1.0):
    """Convert the orientation vector [posAng,lat_planetographic] to the rotation angles"""
    tip = -Q[0]*np.pi/180.0   #'tip' to north
    rotate = -math.atan( math.tan(Q[1]*np.pi/180.0) * (1.0-f)**2)   # 'rotate' the sub-earth planetographic latitude
    return tip, rotate   # planetocentric coordinates
def __rotate2planet__(rotate,tip,b):
    """separate out to keep use consistent!"""
    # out_vec = shape.rotZ(tip,shape.rotX(rotate,b))  # the first way, which is seemingly incorrect
    out_vec = shape.rotX(rotate,shape.rotZ(tip,b))
    return out_vec
def __rotate2obs__(rotate,tip,b):
    """This should be opposite to __rotate2planet__..."""
    out_vec = shape.rotZ(tip,shape.rotX(rotate,b))
    return out_vec

def __findEdge__(atm,b,rNorm,tip,rotate,gtype,printdot=True):
    tmp = (b[0]**2 + b[1]**2)
    try:
        zQ_Trial = np.arange(math.sqrt(1.0-tmp)*1.01,0.0,-0.005)
    except ValueError:
        return None, None
    pclat_zQ_Trial = []
    r_zQ_Trial = []
    r_pclat_zQ_Trial = []
    hitPlanet = False
    geoid = shape.Shape(gtype)
    for zQ in zQ_Trial:
        if printdot:
            print '.',
            sys.stdout.flush()
        # Q
        b_vec = np.array([b[0],b[1],zQ])
        # --> P
        b_vec = __rotate2planet__(rotate,tip,b_vec)
        r1 = np.linalg.norm(b_vec) * rNorm
        r_zQ_Trial.append( r1 )
        # get planetocentric latitude/longitude
        pclat = (180.0/math.pi)*math.asin( np.dot(b_vec,yHat)/np.linalg.norm(b_vec) )
        delta_lng = (180.0/math.pi)*math.atan2( np.dot(b_vec,xHat), np.dot(b_vec,zHat) )
        pclat_zQ_Trial.append(pclat)
        r2 = geoid.calcShape(atm,rNorm,pclat,delta_lng)
        r_pclat_zQ_Trial.append( r2 )
        if r1 <= r2:
            hitPlanet = True
            break
    if not hitPlanet:
        return None, None
    # interpolate to value
    xx = np.flipud(np.array(r_zQ_Trial)-np.array(r_pclat_zQ_Trial))
    yy = np.flipud(np.array(zQ_Trial[0:len(r_zQ_Trial)]))
    #-debug-#plt.figure('testing')
    #-debug-#plt.plot(yy,xx)
    try:
        zQ = np.interp(0.0,xx,yy)
        b = np.array([b[0],b[1],zQ])
        edge = rNorm * __rotate2planet__(rotate,tip,b)
    except ValueError:
        edge = None

    del zQ_Trial, pclat_zQ_Trial, r_zQ_Trial, r_pclat_zQ_Trial, geoid, xx, yy
    return edge, b  # same vector but in P and Q coordinate systems

def compute_ds(atm, b, orientation=None, gtype=None, verbose=False, plot=True):
    """Computes the path length through the atmosphere given:
            b = impact parameter (fractional distance to outer edge at that latitude in observer's coordinates)
            orientation = position angle of the planet [0]='tip', [1]='subearth latitude' """
    if gtype==None:
        gtype = atm.config.gtype
    if gtype == 'gravity':
        verbose = True
    if orientation == None:
        orientation = atm.config.orientation
    path = Ray()
    req = atm.layerProperty[atm.config.LP['R']]   # radius of layers along equator
    rNorm = req[0]
    nr = atm.layerProperty[atm.config.LP['N']]    # refractive index of layers
    P = atm.gas[atm.config.C['P']]
    if (b[0]**2 + b[1]**2) > 1.0:
        return path

    if atm.config.limb == 'sec':
        mu = math.sqrt(1.0 - b[0]**2 - b[1]**2)

    f = 1.0 - atm.config.Rpol/atm.config.Req
    tip, rotate = __computeAspect__(orientation,f)
    print 'intersection:  (%.3f, %.3f)    ' % (b[0],b[1]),
    print 'aspect:  (%.4f,  %.4f)' % (tip*180.0/np.pi,rotate*180.0/np.pi)

    print 'Finding atmospheric edge',
    edge, b = __findEdge__(atm,b,rNorm,tip,rotate,gtype)
    if edge == None:
        return path
    r1 = np.linalg.norm(edge)
    
    # get planetocentric latitude/longitude
    pclat = (180.0/math.pi)*math.asin( np.dot(edge,yHat)/np.linalg.norm(edge) )
    delta_lng = (180.0/math.pi)*math.atan2( np.dot(edge,xHat), np.dot(edge,zHat) )
    geoid = shape.Shape(gtype)
    r2 = geoid.calcShape(atm,rNorm,pclat,delta_lng)
    print ' within %.2f m' % (abs(r1-r2)*100.0)
    geoid.printShort()  #short version of print geoid
    if plot:
        plotStuff(atm=atm,r=rNorm,b=b,gtype=gtype,delta_lng=delta_lng,geoid=geoid,tip=tip,rotate=rotate)

    # initialize - everything below is in planetocentric coordinates, so need to rotate s-vector
    start = np.array([0.0, 0.0, -1.0])
    start = __rotate2planet__(rotate,tip,start)
    s = [start]
    n = [geoid.n]
    r = [geoid.r]
    t_inc = [math.acos( -np.dot(s[-1],n[-1]) )]         # incident angle_0
    nratio = nr[0]/nr[1]                                # refractive index ratio_0
    t_tran = [math.asin(nratio*math.sin(t_inc[-1]))]    # transmitted angle_0

    # return array initialization, use a value to keep indexing numbering for loop
    ds = [0.0]
    layer4ds = [0.0]
    r4ds = [0.0]  # not needed, but for information
    P4ds = [0.0]  # not needed, but for information
    doppler = []

    # loop while in atmosphere
    i = 0                      # loops over the path
    layer = 0                  # keeps track which physical layer you are in
    inAtmosphere = True
    direction = 'ingress'
    while inAtmosphere:
        if verbose:
            print '------------------'
            print '\tstep %d:  layer %d %s ' % (i,layer,direction)
            print '\ts = [%.4f, %.4f, %.4f],  ds = %.4f' % (s[-1][_X],s[-1][_Y],s[-1][_Z],ds[-1])
            geoid.print_a_step()
            print '\tt_inc, tran:  %.8f -> %.8f' % (180.0*t_inc[-1]/math.pi,180.0*t_tran[-1]/math.pi)
        # update s-vector
        s.append( nratio*s[i] + raypathdir[direction]*( nratio*math.cos(t_inc[i])*n[i] - math.cos(t_tran[i])*n[i] ) )
        rNowMag = geoid.rmag
        #---#rScale = rNowMag/req[layer]
        #---#rNextMag= req[layer+raypathdir[direction]] * rScale
        rNextMag = geoid.calcShape(atm,req[layer+raypathdir[direction]],pclat,delta_lng)
        rdots = np.dot( r[i],s[i+1] )

        vw = np.interp(pclat,atm.config.vwlat,atm.config.vwdat)/1000.0
        dopp = 1.0 - (atm.config.omega_m*rNowMag*math.cos(pclat*math.pi/180.0) + vw)*math.sin(delta_lng*math.pi/180.0)/3.0E5

        # get ds, checking for exit, errors and warnings...
        try:
            dsp = -rdots + math.sqrt(rdots**2.0 + rNextMag**2.0 - rNowMag**2.0)
            dsm = -rdots - math.sqrt(rdots**2.0 + rNextMag**2.0 - rNowMag**2.0)
            if direction == 'ingress':
                ds_step = dsm
            elif direction == 'egress':
                ds_step = dsp
        except ValueError:  #tangent layer
            if direction == 'ingress':
                print 'In tangent layer  ',
                direction = 'tangent'
                ds_step = -2.0*rdots
                if ds_step > rNorm:
                    inAtmosphere = False
                    print 'Error:  tangent ds too large'
                    break
                rNextMag = rNowMag
            else:
                inAtmosphere = False
                break
        if ds_step < 0.0:
            if direction != 'egress':
                print 'Error:  ds < 0  (%s:  r.s=%f, ds=%f, [%f,%f])' % (direction,rdots,ds_step,dsp,dsm)
            inAtmosphere = False
            break

        if atm.config.limb == 'sec':  #overwrite ds_step with secant version
            ds_step = abs(rNextMag - rNowMag)/mu

        # append to return arrays
        ds.append(ds_step)
        layer4ds.append(layer)
        r4ds.append(rNowMag)
        P4ds.append(atm.gas[atm.config.C['P']][layer])
        doppler.append(dopp)

        # get next step, double-check r value (and compute n, etc)
        rnext = r[i] + ds[i+1]*s[i+1]
        pclat = (180.0/math.pi)*math.asin( np.dot(rnext,yHat)/np.linalg.norm(rnext) )
        delta_lng = (180.0/math.pi)*math.atan2( np.dot(rnext,xHat), np.dot(rnext,zHat) )
        r2 = geoid.calcShape(atm,req[layer+raypathdir[direction]],pclat,delta_lng)
        if abs(r2-rNextMag) > 10.0:
            print 'Warning:  %f != %f' % (r2,rNextMag)
        r.append( rnext )  # which is also geoid.r, or should be 
        n.append( geoid.n )

        # get new incident angle
        layer+=raypathdir[direction]
        if direction=='tangent':
            direction = 'egress'
            t_inc.append(np.pi/2.0)
        else:
            try:
                t_inc.append( math.acos( -raypathdir[direction]*np.dot( s[i+1],n[i+1] ) ) )
            except ValueError:
                print 't_inc ValueError |s_(i+1)| = %f, |n_(i+1)| = %f - set to previous' % ( np.linalg.norm(s[i+1]), np.linalg.norm(n[i+1]) )
                t_inc.append( t_inc[i] )
                inAtmosphere = False
                break
        i+=1

        # get refractive index ratio and check for exit
        try:
            nratio = nr[layer]/nr[layer+raypathdir[direction]]
            try:
                t_tmp = math.asin(nratio*math.sin(t_inc[-1]))
            except ValueError:
                t_tmp = nratio*np.pi/2.0
            t_tran.append(t_tmp)
        except IndexError:
            inAtmosphere = False
        ### end loop ###

    # Get rid of the first entry, which was just used to make indexing in loop consistent
    del ds[0], layer4ds[0], r4ds[0], P4ds[0]
    path.update(ds=ds,layer4ds=layer4ds,r4ds=r4ds,P4ds=P4ds,doppler=doppler,tip=tip,rotate=rotate,rNorm=rNorm)
    if plot:
        plotStuff(r=np.array(r),ray=path)
    del s, r, n, ds, layer4ds, r4ds, P4ds,geoid, req, nr
    return path

def plotStuff(atm=None,r=None,b=None,gtype=None,delta_lng=None,geoid=None,ray=None,tip=None,rotate=None):
    global plotExist
    if ray == None:
        if not plotExist:
            print 'Generating plots'
            plotExist = True
            if USE_MAYAVI:  ###For now just block everything out if not using MAYAVI
                _r = np.zeros( (2,3) )
                edge, b_tmp = __findEdge__(atm,[-0.1,0.0],1.005*r,tip,rotate,gtype)
                _r[0] = edge
                edge, b_tmp = __findEdge__(atm,[0.1,0.0],1.005*r,tip,rotate,gtype)
                _r[1] = edge
                mlab.plot3d(_r[:,_X],_r[:,_Y],_r[:,_Z], color=(0,0,1),tube_radius=250.0)
                edge, b_tmp = __findEdge__(atm,[0.0,-0.05],1.005*r,tip,rotate,gtype)
                _r[0] = edge
                edge, b_tmp = __findEdge__(atm,[0.0,0.05],1.005*r,tip,rotate,gtype)
                _r[1] = edge
                mlab.plot3d(_r[:,_X],_r[:,_Y],_r[:,_Z], color=(0.96,1.0,0.04),tube_radius=250.0)
                #shape.plotMethods(atm,r,delta_lng=delta_lng,gtypes=[gtype])
                #-#shape.plotMethods(atm,r,delta_lng=0.0,gtypes=[gtype],color='g')
                #-#shape.plotMethods(atm,r,delta_lng=90.0,gtypes=[gtype],color='k')
                #-#shape.plotMethods(atm,r,delta_lng=180.0,gtypes=[gtype],color='k')
                #-#shape.plotMethods(atm,r,delta_lng=270.0,gtypes=[gtype],color='k')
                # plot equator, axis and northern 'halo'
                _x = []; _y = []; _z = []
                theta = np.arange(0.0,2.0*np.pi,0.001)
                for t in theta:
                    _x.append(r*np.cos(t))
                    _z.append(r*np.sin(t))
                    _y.append(0.0)
                mlab.plot3d(_x,_y,_z,color=(0,0,0),opacity=0.5,tube_radius=250)
                _x = []; _y = []; _z = []
                for t in theta:
                    _x.append((r/10.0)*np.cos(t))
                    _z.append((r/10.0)*np.sin(t))
                    _y.append(r*1.2)
                mlab.plot3d(_x,_y,_z,color=(0,0,0),opacity=0.5,tube_radius=250)
                _x = [0.0,0.0]; _y = [-1.5*r,1.5*r]; _z = [0.0,0.0]
                mlab.plot3d(_x,_y,_z,color=(0,0,0),opacity=0.5,tube_radius=250)
                del _x, _y, _z, _r
        plt.figure('observer')
        plt.plot(r*b[0],r*b[1],'.',color='k')
    else:
        if USE_MAYAVI:
            mlab.plot3d(r[:,_X],r[:,_Y],r[:,_Z],color=(1,0,0),tube_radius=150)
        plt.figure('raypath-r')
        plt.plot(ray.r4ds, ray.ds)
        plt.figure('raypath-P')
        plt.semilogy(ray.ds,ray.P4ds)
        plt.axis(ymin = ray.P4ds[-1],ymax=ray.P4ds[0])
        #-debug-#plt.plot(ray.r4ds, ray.layer4ds)

###Test functions
def refractTest(layers):
    """Returns fake refractive index at layers"""
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
