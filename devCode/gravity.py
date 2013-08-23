import numpy as np
import matplotlib.pyplot as plt
import scipy.special as scisp
import math
import planetGlobal as pg

H = {'psi':0, 'lat':1, 'r':2, 'n':3, 't':4}

def doall(b=1.0):
    J = pg.Jupiter()
    S = pg.Saturn()
    U = pg.Uranus()
    N = pg.Neptune()
    P = [J,S,U,N]
    for planet in P:
        calcGeoid(planet,b)

def calcGeoid(planet, b=1.0):
    latstep = 0.1
    pclat = np.arange(0.0,90.0+latstep,latstep)
    grav = []
    psi = []
    psi_ellipse = []
    eGlobal = []
    eLocal = []
    df = []
    r_geoid = []
    r_globalEllipse = []
    if pclat[0]<-89.0:
        r = planet.Rpol*b
    else:
        r = planet.Req*b
    for ii,v in enumerate(pclat):
        radv = v*np.pi/180.0
        vw = np.interp(v,planet.vwlat,planet.vwdat)/1000.0
        omega = planet.omega_m + vw/(r*np.cos(radv))
        g, geoid, ellipse = gravity(v, r, planet.RJ, planet.Req, planet.Rpol, planet.GM, omega, planet.Jn)
        grav.append(g)
        
        psi.append((180.0/np.pi)*geoid[H['psi']])
        psiel = math.atan( (planet.Req/planet.Rpol)*math.tan(radv) )*180.0/np.pi - v
        psi_ellipse.append(psiel)
        df.append( (180.0/np.pi)*geoid[H['psi']] - psiel)
        if not ii%20:
            print '%5.2f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f' % \
                    (v, r, g[0]*1000.0, g[1]*1000.0, g[2]*1000.0, geoid[H['psi']]*180.0/np.pi,psiel,ellipse[0],ellipse[1])
        eGlobal.append([planet.Req*b*math.cos(radv),planet.Rpol*b*math.sin(radv)])
        eLocal.append([r*math.cos(radv),r*math.sin(radv)])
        dlat = latstep*np.pi/180.0
        dr_vec = r*dlat*np.array([geoid[H['t']][0], geoid[H['t']][1]])
        r_vec = geoid[H['r']]+dr_vec
        r = math.sqrt(r_vec[0]**2 + r_vec[1]**2)
        ###r = math.sqrt((ellipse[0]*math.cos(radv + dlat))**2 + (ellipse[1]*math.sin(radv+dlat))**2) #Doesn't work!
        r_geoid.append(r)
        r_globalEllipse.append( math.sqrt( (planet.Req*b*math.cos(radv))**2 + (planet.Rpol*b*math.sin(radv))**2 ) )

    grav = np.array(grav)
    plt.figure(2)
    plt.plot(pclat,grav[:,0])  #gr
    plt.plot(pclat,grav[:,1])  #gphi
    plt.plot(pclat,grav[:,2])  #gtotal
    plt.figure(3)
    plt.plot(pclat,psi)
    plt.plot(pclat,psi_ellipse)
    plt.plot(pclat,df)
    plt.figure(5)
    eGlobal = np.array(eGlobal)
    eLocal = np.array(eLocal)
    plt.plot(eGlobal[:,0],eGlobal[:,1])
    plt.plot(eLocal[:,0],eLocal[:,1])
    r_geoid = np.array(r_geoid)
    r_globalEllipse = np.array(r_globalEllipse)
    plt.figure(6)
    plt.plot(pclat,r_geoid)
    plt.plot(pclat,r_globalEllipse)
    plt.plot(pclat,r_geoid-r_globalEllipse)

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
    psi = math.atan2(gp, gr)
    norm = [math.cos(lat+psi), math.sin(lat+psi)]
    tang = [-math.sin(lat + psi), math.cos(lat + psi)]
    r_vec = [r*math.cos(lat), r*math.sin(lat)]
    geoid = [psi, lat, r_vec, norm, tang]

    # This option doesn't work - too much curvature
    if psi < 1.0e-6:
        rba = Rpol/Req
    else:
        rba = math.tan(lat)/math.tan(lat+psi)
    #rba = Rpol/Req
    a = r/math.sqrt(math.cos(lat)**2 + (rba*math.sin(lat))**2)
    b = a*rba
    ellipse = [a, b]
    
    return g, geoid, ellipse


def test(lat=45.0, b = 0.8, a = 1.0):
    eps = b/a
    latitude = lat*np.pi/180.0
    ds = np.sqrt(np.cos(lat)**2 + (eps*np.sin(lat))**2)
    print 'ds(ellipse) = ',ds
    t = np.arange(0,90)*np.pi/180.0
    Ximpact = a*np.cos(latitude) + np.array([0.0, a/5.0])
    Yimpact = b*np.sin(latitude) + np.array([0.0,0.0])
    plt.plot(Ximpact,Yimpact,'g')
    plt.axis('image')
    plt.grid(1)
    Xplancen = [0.0, a*np.cos(latitude)]
    Yplancen = [0.0, b*np.sin(latitude)]
    plt.plot(Xplancen,Yplancen,'g--')

    # ellipse
    unit = np.sqrt(a**2 + b**2)
    xe = a*np.cos(t)
    ye = b*np.sin(t)
    plt.plot(xe,ye,'g')
    xne = a*np.cos(latitude) + np.array([0.0,b*np.cos(latitude)/unit])
    yne = b*np.sin(latitude) + np.array([0.0,a*np.sin(latitude)/unit])
    plt.plot(xne,yne,'g--')
    fp = open('drawEllipse.scr','a')
    for i in range(len(xe)):
        if not i:
            s = 'PLINE %f,%f\n' % (xe[i],ye[i])
        else:
            s = '%f,%f\n' % (xe[i],ye[i])
        fp.write(s)
    fp.write('\n')
    fp.close()
    # circle
    unit = np.sqrt(a**2 + a**2)
    xc = a*np.cos(t)
    yc = a*np.sin(t)
    plt.plot(xc,yc,'b')
    plt.axis('image')
    xnc = a*np.cos(latitude) + np.array([0.0,a*np.cos(latitude)/unit])
    ync = a*np.sin(latitude) + np.array([0.0,a*np.sin(latitude)/unit])
    plt.plot(xnc,ync,'b--')

                                        

def checkLegendre(i=2):
    spvec = np.arange(-1.0,1.0,0.1)
    for sp in spvec:
        if i==6:
            lin = P6(sp)
        elif i==4:
            lin = P4(sp)
        else:
            lin = P2(sp)
        pyt = scisp.legendre(i)(sp)
        print lin, pyt

def P2(sp):
    return 0.5*(3.0*sp**2 - 1.0)

def P4(sp):
    return (1.0/8.0)*(35.0*pow(sp,4.0) - 30.0*sp**2 + 3.0)

def P6(sp):
    return (1.0/16.0)*(231.0*pow(sp,6.0) - 315.0*pow(sp,4.0) + 105.0*sp**2 - 5.0)
        
