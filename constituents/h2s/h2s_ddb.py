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
GH2S = []
def readInputFiles(path,verbose=False):
    """If needed this reads in the data files for h2s"""
    useLinesUpTo = 200   # index number
    global nlin
    nlin = 0
    if verbose:
        print "Reading h2s lines"
    filename = os.path.join(path,'h2s.lin')
    ifp = open(filename,'r')
    for line in ifp:
        if nlin >= useLinesUpTo:
            break
        nlin+=1
        data = line.split()
        if len(data) == 4:
            f0.append(float(data[0]))
            I0.append(float(data[1]))
            E.append(float(data[2]))
            GH2S.append(float(data[3]))
        else:
            break
    ifp.close()
    if verbose:
        print '   '+str(nlin)+' lines'
    return nlin

def alpha(freq,T,P,X,P_dict,otherPar,units='dBperkm',path='./',verbose=False):

    # Read in data if needed
    if len(f0)==0:
        readInputFiles(path,verbose)

    P_h2 = P*X[P_dict['H2']]
    P_he = P*X[P_dict['HE']]
    P_h2s= P*X[P_dict['H2S']]
    GH2 = 1.960
    GHe = 1.200
    ZH2 = 0.000
    ZHe = 0.000
    ZH2S= 0.000
    D = 1.28
    C = 1.0
    n_dvl = 0.7
    n_int = 3.0/2.0
    delta = D*P_h2s

    alpha_h2s = []
    for f in freq:
        f2 = f**2
        alpha = 0.0
        for i in range(nlin):
            gamma = pow((T0/T),n_dvl)*(GH2*P_h2 + GHe*P_he + GH2S[i]*P_h2s)
            g2 = gamma**2
            zeta = gamma
            z2 = zeta**2
            ITG = I0[i]*math.exp(-((1.0/T)-(1.0/T0))*E[i]*hck)
            num = (gamma-zeta)*f2 + (gamma+zeta)*( pow(f0[i]+delta,2.0) + g2 - z2)
            den = pow((f2 - pow(f0[i]+delta,2.0) - g2 + z2),2.0) + 4.0*f2*g2
            shape = GHz*2.0*pow(f/f0[i],2.0)*num/(math.pi*den)
            alpha += shape*ITG

        a = coef*(P_h2s/T0)*pow((T0/T),n_int+2)*alpha
        if units=='dBperkm':
            a*=434294.5
        alpha_h2s.append(a)
        
    return alpha_h2s

