import planet
import numpy as np
import scipy.signal
import matplotlib.pyplot as plt
import os
import math

#from https://science.nrao.edu/facilities/vla/docs/manuals/oss/performance/resolution
vla = {} #band, GHz, arcsec
vla['4'] = [0.074, 24.00, 80.00, 260.00, 850.0]
vla['P'] = [0.350,  5.60, 18.50,  60.00, 200.0]
vla['L'] = [1.500,  1.30,  4.30,  14.00,  46.0]
vla['S'] = [3.000,  0.65,  2.10,   7.00,  23,0]
vla['C'] = [6.000,  0.33,  1.00,   3.50,  12.0]
vla['X'] = [10.00,  0.20,  0.60,   2.01,   7.2]
vla['Ku']= [15.00,  0.13,  0.42,   1.40,   4.6]
vla['K'] = [22.00,  0.089, 0.28,   0.95,   3.1]
vla['Ka']= [33.00,  0.059, 0.19,   0.63,   2.1]
vla['Q'] = [45.00,  0.043, 0.14,   0.47,   1.5]
f = 0
A = 1
B = 2
C = 3
D = 4

band = 'C'
arrayConfig = A

angle = 90.0
delta_b = 0.01

freq = vla[band][f]
beamwidth = vla[band][arrayConfig]
print 'Beamwidth:  ',beamwidth
beamwidth = beamwidth/3600.0*(math.pi/180.0)
radius = 71492.
distance = 777908927.
res = math.atan(delta_b*radius/distance)*(180.0/math.pi)*3600.0   #just FYI
b05 = math.tan(beamwidth/2.0)*(distance/radius)
print 'Resolution = ',res
print 'b05 = ',b05
kerexp = math.log(2.0)/b05**2


genFiles = False
if genFiles:
    j = planet.planet('jupiter')
    bv = [angle]
    for v in np.arange(0.0,1.2,delta_b):
        bv.append(v)
    j.run(freqs=freq,b=bv)

else:
    shapefiles = ['jup0shape_C.dat']#,'jup45shape_C.dat','jup90shape_C.dat']
    secfiles = ['jup0sec_C.dat','jup45sec_C.dat','jup90sec_C.dat']
    cols = ['b','r','k','g','y','m']

    jdat = []
    for i,fil in enumerate(shapefiles):

        # Do shape files
        label = fil.split('.')[0]
        tmpData = np.loadtxt(fil)
        b = tmpData[:,2]
        b = np.concatenate([-np.flipud(b),b])
        np.delete(b,len(b)/2)

        Ttmp = tmpData[:,3]
        centerTemp = Ttmp[0]
        Ttmp = np.concatenate([np.flipud(Ttmp),Ttmp])
        np.delete(Ttmp,len(Ttmp)/2)
        jdat.append(Ttmp)
        ker = np.exp(-kerexp*b**2)
        ker = ker/np.sum(ker)
        cnvshape = scipy.signal.convolve(Ttmp,ker,'same')

        plt.figure(1)
        plt.plot(b,cnvshape,color=cols[i],linestyle='-',label=label)

        # Do power law version
        g = np.where(abs(b) < 1.0)
        costh = np.zeros(np.shape(b))
        costh[g] = np.sqrt(1.0 - b[g]**2)**0.16
        scTmp = centerTemp*costh
        cnvshape = scipy.signal.convolve(scTmp,ker,'same')
#        plt.plot(b,cnvshape,color='g',linewidth=4)
        plt.plot(b,scTmp,color='g',linewidth=4)

        # Do sec files
        label = secfiles[i].split('.')[0]
        tmpData = np.loadtxt(secfiles[i])
        b = tmpData[:,2]
        b = np.concatenate([-np.flipud(b),b])
        np.delete(b,len(b)/2)

        sTtmp = tmpData[:,3]
        sTtmp = np.concatenate([np.flipud(sTtmp),sTtmp])
        np.delete(sTtmp,len(sTtmp)/2)
        jdat.append(sTtmp)
        ker = np.exp(-kerexp*b**2)
        ker = ker/np.sum(ker)
        cnvsec = scipy.signal.convolve(sTtmp,ker,'same')

        plt.plot(b,cnvsec,color=cols[i],linestyle='--',label=label)
	plt.plot(b,scTmp,color='g')

        plt.figure(2)
        plt.plot(b,ker)

        # Plot unconvolved versions
        plt.figure(3)
        plt.plot(b,Ttmp)
        plt.plot(b,scTmp)
        plt.plot(b,sTtmp)
