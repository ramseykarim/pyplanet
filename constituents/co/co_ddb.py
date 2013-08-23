import math
import os.path

# Some constants
coef = 7.244E+21     # coefficient from GEISA_PP.TEX eq. 14
T0 = 296.0           # reference temperature in K
hck = 1.438396       # hc/k  [K cm]
GHz = 29.9792458     # conversion from cm^-1 to GHz

#Set data arrays
f0 = []
I0 = []
E = []
def readInputFiles(path,verbose=False):
    """If needed this reads in the data files for co"""
    useLinesUpTo = 26   # index number
    global nlin
    nlin = 0
    if verbose:
        print "Reading co lines"
    filename = os.path.join(path,'co.lin')
    ifp = open(filename,'r')
    for line in ifp:
        if nlin >= useLinesUpTo:
            break
        nlin+=1
        data = line.split()
        if len(data) == 3:
            f0.append(float(data[0]))
            I0.append(float(data[1]))
            E.append(float(data[2]))
        else:
            break
    ifp.close()
    if verbose:
        print '   '+str(nlin)+' lines'
    return nlin

def alpha(freq,T,P,X,P_dict,otherPar,units='dBperkm',path='./',verbose=False):

    ###Voigt coefficients from Janssen p67
    PLimits = [0.001,0.1]
    avoigt = [122.60793178,214.38238869,181.92853309, 93.15558046, 30.18014220, 5.91262621, 0.56418958,0.0]
    bvoigt = [122.60793178,352.73062511,457.33447878,348.70391772,170.35400182,53.99290691,10.47985711,1.0]
    
    # Read in data if needed
    if len(f0)==0:
        readInputFiles(path,verbose)

    P_h2 = P*X[P_dict['H2']]
    P_he = P*X[P_dict['HE']]
    P_co = P*X[P_dict['CO']]

    GH2 = 1.960
    GHe = 1.200
    GCO = 6.000
    n_dvl = 0.7
    n_int = 3.0/2.0
    gamma = pow((T0/T),n_dvl)*(GH2*P_h2 + GHe*P_he + GCO*P_co)
    g2 = gamma**2
    zeta = 0.0
    z2 = zeta**2
    delta = 0.0

    alpha_co = []
    for f in freq:
        f2 = f**2
        alpha = 0.0
        for i in range(nlin):
            ITG = I0[i]*math.exp(-((1.0/T)-(1.0/T0))*E[i]*hck)
            w = (P-PLimits[0])/(PLimits[1]-PLimits[0])
            if w<0.0:
                w=0.0
            elif w>1.0:
                w=1.0
            shape_Voigt = 0.0
            if P<=PLimits[1] or otherPar=='voigt' or otherPar=='diff':
                ###Doppler broadening Janssen p59
                betaD = 4.3e-7*math.sqrt(T/28.0)*f
                ###Voigt Janssen p67
                num = 0.0 + 0.0j
                den = 0.0 + 0.0j
                xi = gamma/betaD + (1.0j)*(f-f0[i])/betaD
                for jjj in range(len(avoigt)):
                    num+=avoigt[jjj]*(xi**jjj)
                    den+=bvoigt[jjj]*(xi**jjj)
                val = num/den
                shape_Voigt = GHz*( 1.0/(math.sqrt(math.pi)*betaD) )*val.real
            shape_VVW = 0.0
            if P>=PLimits[0] or otherPar=='vvw' or otherPar=='diff':
                num = (gamma-zeta)*f2 + (gamma+zeta)*( pow(f0[i]+delta,2.0) + g2 - z2)
                den = pow((f2 - pow(f0[i]+delta,2.0) - g2 + z2),2.0) + 4.0*f2*g2
                shape_VVW = GHz*2.0*pow(f/f0[i],2.0)*num/(math.pi*den)
            #print 'co_ddb.py: L94 hardcode w=0.0'
            #w = 0.0
            shape = w*shape_VVW + (1.0-w)*shape_Voigt

            if otherPar == 'voigt':
                alpha+=shape_Voigt
            elif otherPar == 'vvw':
                alpha+=shape_VVW
            elif otherPar == 'diff':
                alpha+=(shape_Voigt-shape_VVW)
            else:
                alpha+=shape*ITG
        a = coef*(P_co/T0)*pow((T0/T),n_int+2)*alpha
        if otherPar=='vvw' or otherPar=='voigt' or otherPar=='diff':
            print 'testtesttestP = ',P
            a = alpha
        if units=='dBperkm':
            a*=434294.5
        alpha_co.append(a)
        
    return alpha_co

