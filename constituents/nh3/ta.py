import nh3_kd
import nh3_bg
import matplotlib.pyplot as plt
import math
import numpy as np

jupiterTestP = [1.,10.,20.,100.,200.,400.,1000.,2000.,4000.,7500.,18000.]
jupiterTestT = [150.,300.,400.,700.,800.,1000.,1200.,1600.,2000.,2400.,3200.,]

f = []
fmin = 1.0
fmax = 500.0
fstep = 1.0
n = int(math.ceil((fmax-fmin)/fstep))+1
for i in range(n):
    f.append(fmin + i*fstep)
otherPar = []

##print 'Fig. 5.2'
##P = 1.009
##T = 216.4
##X_nh3 = 0.0095
##X_he = 0.1347
##X_h2 = 1.0 - (X_nh3+X_he)
##X_partial=[X_h2,X_he,X_nh3]
##P_dict = {'H2':0,'HE':1,'NH3':2}
##a = nh3_kd.alpha(f,T,P,X_partial,P_dict,otherPar)
##plt.semilogy(f,a)
##plt.axis([1,25,.01,1000])
##
##
##plt.figure('vsP')
##f = [0.5]
##T = 200.0
##Pvals = np.arange(1.0,100.0,1.0)
##alpha = []
##for P in Pvals:
##    alpha.append(nh3_kd.alpha(f,T,P,X_partial,P_dict,otherPar))
##plt.plot(Pvals,alpha)
##
##plt.figure('vsT')
##f = [0.5]
##P = 200.0
##Tvals = np.arange(100.0,1000.0,10.0)
##alpha = []
##for T in Tvals:
##    alpha.append(nh3_kd.alpha(f,T,P,X_partial,P_dict,otherPar))
##plt.plot(Tvals,alpha)


f = [4.0, 8.0,20.0]
fclr = ['b','r','g']


jupiterTestP = [1.,2.,4.,10.,20.,100.,200.,400.,1000.,2000.,4000.,7500.,18000.]
jupiterTestT = [150.,200.,250.,300.,400.,700.,800.,1000.,1200.,1600.,2000.,2400.,3200.,]
X_nh3 = 0.005
X_he = 0.10
X_h2 = 1.0 - (X_nh3+X_he)
X_partial=[X_h2,X_he,X_nh3]
P_dict = {'H2':0,'HE':1,'NH3':2}

a_kd = []
a_bg = []
for i in range(len(jupiterTestP)):
    T = jupiterTestT[i]
    P = jupiterTestP[i]
    a_kd.append(nh3_kd.alpha(f,T,P,X_partial,P_dict,otherPar))
    a_bg.append(nh3_bg.alpha(f,T,P,X_partial,P_dict,otherPar))
    #print '%5.0f\t%4.0f\t%.3f\t%.3f' % (P,T,a_kd[i],a_bg[i])
a_kd = np.array(a_kd)
a_bg = np.array(a_bg)
plt.figure('samples')
for i in range(len(f)):
    s = '%.0f GHz' % (f[i])
    print s
    plt.plot(jupiterTestP,a_kd[:,i],'bo')
    plt.plot(jupiterTestP,a_kd[:,i],'b',label=s)
    plt.plot(jupiterTestP,a_bg[:,i],'ro')
    plt.plot(jupiterTestP,a_bg[:,i],'r',label=s)
plt.xscale('log')
plt.legend()
plt.xlabel('Pressure [bars]')
plt.ylabel(r'$\alpha$ [dB/km]')
s = r'%.1f GHz:  %.3f H$_2$, %.3f He, %.3f NH$_3$' % (f[0],X_h2,X_he,X_nh3)
plt.title(s)
print X_partial

if True:
    fp = open('juptest.dat','w')
    f = [4.0, 8.0,20.0]
    fclr = ['b','r','g']
    jupdat = np.loadtxt('jupiter.paulSolar')
    jupiterTestP = jupdat[:,2]
    jupiterTestT = jupdat[:,1]
    X_nh3 = jupdat[:,6]
    X_he = jupdat[:,4]
    X_h2 = jupdat[:,3]
    a_kd = []
    a_bg = []
    for i in range(len(jupiterTestP)):
        T = jupiterTestT[i]
        P = jupiterTestP[i]
        X_partial = [X_h2[i],X_he[i],X_nh3[i]]
        a_kd.append(nh3_kd.alpha(f,T,P,X_partial,P_dict,otherPar))
        a_bg.append(nh3_bg.alpha(f,T,P,X_partial,P_dict,otherPar))
        fp.write('%f\t%f\t%f\t%f\t%f\t%s\t%s\n' % (P,T,X_partial[0],X_partial[1],X_partial[2],str(a_kd[i]),str(a_bg[i])))
    a_kd = np.array(a_kd)
    a_bg = np.array(a_bg)
    plt.figure('ssampless')
    for i in range(len(f)):
        plt.plot(jupiterTestP,a_kd[:,i],fclr[i])
        plt.plot(jupiterTestP,a_bg[:,i],fclr[i]+'--')
    plt.xscale('log')
    plt.xlabel('Pressure [bars]')
    plt.ylabel(r'$\alpha$ [dB/km]')
    #s = r'%.1f GHz:  %.3f H$_2$, %.3f He, %.3f NH$_3$' % (f[0],X_h2,X_he,X_nh3)
    #plt.title(s)
    print X_partial
    fp.close()
