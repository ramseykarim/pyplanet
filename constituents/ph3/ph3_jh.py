import math
import os.path

# Some constants
coef = 7.244E+21        # coefficient from GEISA_PP.TEX eq. 14
T0 = 300.0              # reference temperature in K
hck = 1.438396          # hc/k  [K cm]
GHz = 29.9792458        # conversion from cm^-1 to GHz
PI = math.pi

# Set data arrays and values
f0 = []
I0 = []
E = []
WgtI0 = []
WgtFGB = []
WgtSB = []
def readInputFiles(path,verbose=False):
    """If needed this reads in the data files for ph3"""

    # Read line info
    useLinesUpTo = 320
    filename = os.path.join(path,'ph3jh.lin')
    if verbose:
        print "Reading ph3 data "+filename
    ifp = open(filename,'r')
    global nlin
    nlin = 0
    for line in ifp:
        if nlin >= useLinesUpTo:
            break
        nlin+=1
        data = line.split()
        if len(data) == 3:
            f0.append(float(data[0]))
            I0.append(float(data[1])/(GHz*1.0E17))
            E.append(float(data[2]))
        else:
            break
    ifp.close()
    if verbose:
        print '   '+str(nlin)+' lines'

    # Read weight info
    filename = os.path.join(path,'PH3WGT.dat')
    ifp = open(filename,'r')
    numWgt = 40
    global nwgt
    nwgt = 0
    for line in ifp:
        if nwgt >= numWgt:
            break
        nwgt+=1
        data = line.split()
        if len(data) == 4:
            WgtI0.append(float(data[1]))
            WgtFGB.append(float(data[2]))
            WgtSB.append(float(data[3]))
        else:
            break
    ifp.close()
    if verbose:
        print '    '+str(nwgt)+' weights'
    return nlin

def alpha(freq,T,P,X,P_dict,otherPar,units='dBperkm',path='./',verbose=False):
    """Computes the absorption due to ph3"""

    # Read in data if needed
    if len(f0)==0:
        readInputFiles(path,verbose=verbose)

    P_h2 = P*X[P_dict['H2']]
    P_he = P*X[P_dict['HE']]
    P_ph3= P*X[P_dict['PH3']]
    
    # Line parameters
    GH2   = 3.2930
    GHe   = 1.6803
    GPH3  = 4.2157
    n_dvl = 2.0/3.0
    n_int = 3.0/2.0
    zeta  = 0.0
    z2 = zeta**2
    delta = 0.0

    alpha_ph3 = []
    for f in freq:
        f2 = f**2
        alpha = 0.0
        for i in range(nlin):
            if i<nwgt:
                gamma = pow((T0/T),n_dvl)*(GH2*P_h2 + GHe*P_he)*WgtFGB[i] + pow((T0/T),1.0)*GPH3*P_ph3*WgtSB[i]
                g2 = gamma**2
                ITG = I0[i]*WgtI0[i]*math.exp(-((1.0/T)-(1.0/T0))*E[i]*hck)
            else:
                gamma = pow((T0/T),n_dvl)*(GH2*P_h2 + GHe*P_he) + pow((T0/T),1.0)*GPH3*P_ph3
                g2 = gamma**2
                ITG = I0[i]*math.exp(-((1.0/T)-(1.0/T0))*E[i]*hck);
            num = (gamma-zeta)*f2 + (gamma+zeta)*( (f0[i]+delta)**2 + g2 - z2)
            den = (f2 - (f0[i]+delta)**2 - g2 + z2)**2 + 4.0*f2*g2
            shape = GHz*2.0*((f/f0[i])**2)*num/(PI*den)
            alpha += shape*ITG
        a= coef*(P_ph3/T0)*pow((T0/T),n_int+2)*alpha
        if units=='dBperkm':
            a*=434294.5
        alpha_ph3.append(a)

    return alpha_ph3
