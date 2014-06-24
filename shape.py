import numpy as np
import math
import scipy.special as scisp
import matplotlib.pyplot as plt
import matplotlib.colors as clr
USE_MAYAVI = False
if USE_MAYAVI:
    import mayavi.mlab as mlab
import atmosphere
import string
_X = 0
_Y = 1
_Z = 2

class Shape:
    """Calculates radius and normal to planet.  Types are:
            gravity, reference, sphere, ellipse"""
    def __init__(self, gtype = 'reference'):
        self.gtype = gtype
        self.gamma = 0.0
        self.pclat = 0.0
        self.pglat = 0.0
        self.delta_lng = 0.0
        self.omega = 0.0
        self.r = np.zeros(3)
        self.rmag = 0.0
        self.n = np.zeros(3)
        self.t = np.zeros(3)
        self.g = []
        self.gmag = 0.0
        self.g_static = 0.0
        self.default_latstep = 0.05
        self.referenceGeoid = False
        self.referenceRadius = 0.0
        self.referenceDecimate = 10
    def __str__(self):
        s = '-------------------------Geoid---------------------------\n'
        s+= 'type:  %s, gamma = %.2f deg, lat = %.2f deg, lng = %.2f\n' % (self.gtype,self.gamma*180.0/np.pi, self.pclat*180.0/np.pi, self.delta_lng*180.0/np.pi)
        s+= 'r = [%.2f,  %.2f,  %.2f] km, |r| = %.2f km\n' % (self.r[_X], self.r[_Y], self.r[_Z],self.rmag)
        s+= 'n = [%.4f,  %.4f,  %.4f]\n' % (self.n[_X], self.n[_Y], self.n[_Z])
        s+= 't = [%.4f,  %.4f,  %.4f]\n' % (self.t[_X], self.t[_Y], self.t[_Z])
        s+= 'g = [%.6f,  %.6f], |g| = %.6f\n' % (self.g[0], self.g[1], self.gmag)
        s+= '---------------------------------------------------------\n'
        return s
    def printShort(self):
        s = '%s:  gamma = %.2f, lat = %.2f, lng = %.2f, |r|=%.2f km\n' % (self.gtype,self.gamma*180.0/np.pi, self.pclat*180.0/np.pi, self.delta_lng*180.0/np.pi,self.rmag)
        print s
    def print_a_step(self):
        s = '\ttype:  %s, gamma = %.2f deg, lat = %.2f deg, lng = %.2f  |g| = %.2f\n' % (self.gtype,self.gamma*180.0/np.pi, self.pclat*180.0/np.pi, self.delta_lng*180.0/np.pi, self.gmag*1000.0)
        s+= '\tr = [%.2f,  %.2f,  %.2f] km, |r| = %.2f km\n' % (self.r[_X], self.r[_Y], self.r[_Z],self.rmag)
        s+= '\tn = [%.4f,  %.4f,  %.4f]' % (self.n[_X], self.n[_Y], self.n[_Z])
        print s
        return s

    ###------------------------------------General handling function----------------------------------------
    def calcShape(self, planet, r, pclat=90.0, delta_lng=0.0, direction=1.0, gtype=None, plot3d=False, latstep='default', color='k'):
        if gtype == None:
            gtype = self.gtype
        if gtype == 'gravity':
            r = self.calcGeoid(planet, r, pclat=pclat, delta_lng=delta_lng, plot3d=plot3d, latstep=latstep)
        elif gtype == 'reference':
            r = self.calcFromReference(planet, r, pclat=pclat, delta_lng=delta_lng, plot3d=plot3d, latstep=latstep)
        elif gtype == 'sphere':
            r = self.calcEllipse(planet, r, pclat=pclat, delta_lng=delta_lng, gtype=gtype, plot3d=plot3d, color=color)
        elif gtype == 'ellipse':
            r = self.calcEllipse(planet, r, pclat=pclat, delta_lng=delta_lng, gtype=gtype, plot3d=plot3d, color=color)
        else:
            print gtype+' not a valid planet shape'
            r = None
        return r
    
    ###---------------------Below are the specific shape handlers----------------------
    ###'reference' - fits a geoid at the reference pressure and scales
    def calcFromReference(self, planet, r, pclat=90.0, delta_lng=0.0, direction=1.0, plot3d=False, latstep='default',color='k'):
        print 'Need to fix plot etc and add direction!!!!'
        if type(self.referenceGeoid) == bool:
            print 'Calculating reference geoid at P_ref=%f' % (planet.config.p_ref)
            self.referenceRadius=planet.config.Req
            self.calcGeoid(planet,planet.config.Req,-90,0.0,plot3d=False,latstep=latstep)
            tmp = []
            for i in range(len(self.referenceList)):
                if not i%self.referenceDecimate:
                    tmp.append(self.referenceList[i])
            self.referenceGeoid = np.flipud(np.array(tmp))
            self.calcGeoid(planet,planet.config.Req,90,0.0,plot3d=False,latstep=latstep)
            tmp = []
            for i in range(len(self.referenceList)):
                if not i%self.referenceDecimate:
                    tmp.append(self.referenceList[i])
            del tmp[0] 
            self.referenceGeoid = np.vstack( (self.referenceGeoid,np.array(tmp)) )
        rlat  = np.interp(pclat,self.referenceGeoid[:,0],self.referenceGeoid[:,1]) * (r/self.referenceRadius)
        gamma = np.interp(pclat,self.referenceGeoid[:,0],self.referenceGeoid[:,2])

        lat = pclat*np.pi/180.0
        lng = delta_lng*np.pi/180.0
        norm = np.array([0.0,   math.sin(lat+gamma),   math.cos(lat+gamma)])
        norm = rotY(lng,norm)
        tang = np.array([0.0,   math.cos(lat+gamma),  -math.sin(lat+gamma)])
        tang = rotY(lng,tang)
        r_vec = np.array([0.0, rlat*math.sin(lat),       rlat*math.cos(lat)])
        r_vec = rotY(lng,r_vec)
        self.gamma = gamma
        self.pclat = lat
        self.pglat = lat + gamma
        self.delta_lng = lng
        self.r = r_vec
        self.rmag = np.linalg.norm(r_vec)
        self.n = norm
        self.t = tang
        self.g = [0.0,0.0]
        self.gmag = 0.0

        if plot3d:
            _x = []; _y = []; _z = []; _p = []; _r = []; _g = []
            for i,vlat in enumerate(self.referenceGeoid[:,0]):
                lat = (np.pi/180.0)*vlat
                r_vec = self.referenceGeoid[i,1] * (r/self.referenceRadius)*np.array([0.0, math.sin(lat), math.cos(lat)])
                r_vec = rotY(lng,r_vec)
                _x.append(r_vec[_X])
                _y.append(r_vec[_Y])
                _z.append(r_vec[_Z])
                #-debug-#_p.append(vlat)
                #-debug-#_r.append(np.linalg.norm(r_vec))
                #-debug-#_g.append(self.referenceGeoid[i,2])
            print ' need to change over to mlab... '
            #-debug-#fp = open('testReference.out','w')
            #-debug-#for i in range(len(_x)):
            #-debug-#    s = '%f\t%f\t%f\t%f\t%f\t%f\n' % (_x[i],_y[i],_z[i],_p[i],_r[i],_g[i])
            #-debug-#    fp.write(s)
            #-debug-#fp.close()
            del _x, _y, _z, _p, _r, _g
        return self.rmag

    ###'gravity' - does the full thing, but is very time-consuming
    def calcGeoid(self, planet, r, pclat=90.0, delta_lng=0.0, direction=1.0, plot3d=False, latstep='default',color='k'):
        """Starts at equatorial radius and moves north or south to compute geoid at pclat"""
        print 'Need to fix plot etc and add direction!!!!'
        if plot3d:
            _x = []; _y = []; _z = []; _p = []; _r = []; _g = []
        if type(latstep) == str:
            latstep = self.default_latstep
        if pclat == 0.0:
            nsp = 1.0
        else:
            nsp = np.sign(pclat)
        latstep = nsp*latstep
        pclatSteps = np.arange(0.0,pclat+latstep,latstep)

        if self.gtype == 'reference':
            self.referenceList = []

        GM = np.interp(r,planet.layerProperty[planet.config.LP['R']],planet.layerProperty[planet.config.LP['GM']])
        for latv in pclatSteps:
            radlatv = latv*np.pi/180.0
            vw = np.interp(latv,planet.config.vwlat,planet.config.vwdat)/1000.0
            self.omega = planet.config.omega_m + vw/(r*np.cos(radlatv))
            self.__gravity__(latv, delta_lng, r, GM, self.omega, planet.config.Jn, planet.config.RJ)
            if plot3d:
                _x.append(self.r[_X])
                _y.append(self.r[_Y])
                _z.append(self.r[_Z])
                #-debug-#_p.append(latv)
                #-debug-#_r.append(self.rmag)
                #-debug-#_g.append(self.gamma)
            dlat = latstep*np.pi/180.0
            dr_vec = r*dlat*self.t
            r_vec = self.r+dr_vec
            r = np.linalg.norm(r_vec)
            if self.gtype == 'reference':
                self.referenceList.append([latv,r,self.gamma])
        if plot3d:
            print 'need to change over to mlab...'
            #-debug-#for i in range(len(_x)):
            #-debug-#    fp = open('testGeoid.out','w')
            #-debug-#    for i in range(len(_x)):
            #-debug-#        s = '%f\t%f\t%f\t%f\t%f\t%f\n' % (_x[i],_y[i],_z[i],_p[i],_r[i],_g[i])
            #-debug-#        fp.write(s)
            #-debug-#    fp.close()
            del _x, _y, _z, _p, _r, _g
        del pclatSteps
        return self.rmag
    def __gravity__(self, pclat, delta_lng, r, GM, omega, Jn, RJ):
        self.g_static = GM/r**2
        lat = pclat*np.pi/180.0
        lng = delta_lng*np.pi/180.0
        if lat == 0.0:
            nsl = 1.0
        else:
            nsl = np.sign(lat)
        dphi = nsl*0.00001
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
        gr = self.g_static*(1.0 - Sr) - (2.0/3.0)*(omega**2.0)*r*(1.0 - scisp.legendre(2)(sp))
        # g(phi)
        dP = (3.0*sp*np.sqrt(1.0-sp**2))
        gp = (1.0/3.0)*(omega**2.0)*r*dP + self.g_static*Sp
        gt = math.sqrt(gr**2 + gp**2)
        # geoid
        gamma = math.atan2(gp, gr)
        norm = np.array([0.0,   math.sin(lat+gamma),   math.cos(lat+gamma)])
        norm = rotY(lng,norm)
        tang = np.array([0.0,   math.cos(lat+gamma),  -math.sin(lat+gamma)])
        tang = rotY(lng,tang)
        r_vec = np.array([0.0, r*math.sin(lat),       r*math.cos(lat)])
        r_vec = rotY(lng,r_vec)
        self.gamma = gamma
        self.pclat = lat
        self.pglat = lat + gamma
        self.delta_lng = lng
        self.r = r_vec
        self.rmag = np.linalg.norm(r_vec)
        self.n = norm
        self.t = tang
        self.g = [gr,gp]
        self.gmag = gt

    ### 'ellipse' or 'circle' - simply does the equivalent ellipse or circle
    def calcEllipse(self, planet, r, pclat, delta_lng, direction=1.0, gtype='ellipse', plot3d=False, color='k'):
        a = r
        if gtype=='ellipse':
            b = (planet.config.Rpol/planet.config.Req)*r
        else:
            b = r
        lat = pclat*np.pi/180.0
        lng = delta_lng*np.pi/180.0
        if lat == 0.0:
            nsl = 1.0
            lat=1.0E-6   ###WHY????!!!!  This [largely] fixes the ValueError on gamma below
        else:
            nsl = np.sign(lat)

        norm = np.array([0.0,  a*math.sin(lat), b*math.cos(lat)])
        norm = norm/np.linalg.norm(norm)
        norm = rotY(lng,norm)
        tang = np.array([0.0, -b*math.cos(lat), a*math.sin(lat)])
        tang = tang/np.linalg.norm(tang)
        tang = rotY(lng,tang)
        r_vec = np.array([0.0, b*math.sin(lat), a*math.cos(lat)])
        r_vec = rotY(lng,r_vec)
        self.rmag = np.linalg.norm(r_vec)
        GM = np.interp(r,planet.layerProperty[planet.config.LP['R']],planet.layerProperty[planet.config.LP['GM']])
        self.g_static = GM/self.rmag**2
        try:
            self.gamma = nsl*math.acos(np.dot(r_vec,norm)/self.rmag)  # don't need to worry about direction here, since norm etc defined
        except ValueError:
            arg = np.dot(r_vec,norm)/self.rmag
            #print norm,np.linalg.norm(norm),r_vec,self.rmag
            print 'gamma error:  [%f,%f,%f,%f,%f]' % (arg,a,b,pclat,delta_lng),
            print '...but proceeding anyway by setting gamma=0.0'
            self.gamma = 0.0
        self.pclat = lat
        self.pglat = lat + self.gamma
        self.delta_lng = lng
        self.r = r_vec
        self.n = norm
        self.t = tang
        self.g = [self.g_static,self.g_static]
        self.gmag = self.g_static

        if plot3d:
            _x = []; _y = []; _z = []; _p = []; _r = []; _g = []
            lat = abs(pclat)
            latstep = 0.1
            for vlat in np.arange(-lat,lat+latstep,latstep):
                nsl=1.0
                if vlat<0.0:
                    nsl=-1.0
                r_vec = np.array([0.0, b*math.sin(vlat*np.pi/180.0), a*math.cos(vlat*np.pi/180.0)])
                r_vec = rotY(lng,r_vec)
                norm = np.array([0.0,  a*math.sin(vlat*np.pi/180.0), b*math.cos(vlat*np.pi/180.0)])
                norm = norm/np.linalg.norm(norm)
                norm = rotY(lng,norm)
                _x.append(r_vec[_X])
                _y.append(r_vec[_Y])
                _z.append(r_vec[_Z])
                #-debug-#_r.append(np.linalg.norm(r_vec))
                #-debug-#_p.append(vlat)
                #-debug-#try:
                #-debug-#    _g.append(nsl*math.acos(np.dot(r_vec,norm)/np.linalg.norm(r_vec)))
                #-debug-#except ValueError:
                #-debug-#    _g.append( 0.0 )
            colorm = clr.colorConverter.to_rgb(color)
            if USE_MAYAVI:
                mlab.plot3d(_x,_y,_z,color=colorm,opacity=0.5,tube_radius=250)
            #-debug-#s = 'test'+string.capitalize(gtype)+'.out'
            #-debug-#fp = open(s,'w')
            #-debug-#for i in range(len(_x)):
            #-debug-#    s = '%f\t%f\t%f\t%f\t%f\t%f\n' % (_x[i],_y[i],_z[i],_p[i],_r[i],_g[i])
            #-debug-#    fp.write(s)
            #-debug-#fp.close()
            del _x, _y, _z, _p, _r, _g
        return self.rmag
        
def plotMethods(planet,r,latstep=0.1,delta_lng=0.0, gtypes=['ellipse','sphere','reference','gravity'],color=None):
    if color == None:
        colors = ['k','r','g','b','y','m','c']
    else:
        colors = []
        for i in range(len(gtypes)):
            colors.append(color)
    if type(latstep) == str:
        latstep = geoid.default_latstep
    for i,gtype in enumerate(gtypes):
        geoid = Shape(gtype)
        r = geoid.calcShape(planet, r, pclat=90.0, delta_lng=delta_lng, plot3d=True, latstep=latstep, color=colors[i])
    del geoid

def rotX(x,V):
    Rx = np.array([[1.0,       0.0,        0.0],
                   [0.0, np.cos(x), -np.sin(x)],
                   [0.0, np.sin(x),  np.cos(x)]])
    return np.dot(Rx,V)
def rotY(y,V):
    Ry = np.array([[ np.cos(y), 0.0, np.sin(y)],
                   [ 0.0,       1.0,       0.0],
                   [-np.sin(y), 0.0, np.cos(y)]])
    return np.dot(Ry,V)
def rotZ(z,V):
    Rz = np.array([[np.cos(z), -np.sin(z), 0.0],
                   [np.sin(z),  np.cos(z), 0.0],
                   [      0.0,        0.0, 1.0]])
    return np.dot(Rz,V)
