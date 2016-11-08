from ta_NH3setup import *

if False:
    print 'Fig. 5.2'
    f = []
    fmin = 1.0
    fmax = 500.0
    fstep = 1.0
    n = int(math.ceil((fmax-fmin)/fstep))+1
    for i in range(n):
        f.append(fmin + i*fstep)
    P = 1.009
    T = 216.4
    X_nh3 = 0.0095
    X_he = 0.1347
    X_h2 = 1.0 - (X_nh3+X_he)
    X_partial=[X_h2,X_he,X_nh3]
    a = nh3_kd.alpha(f,T,P,X_partial,P_dict,otherPar)
    plt.semilogy(f,a)
    plt.axis([1,25,.01,1000])

    plt.figure('vsP')
    f = [0.5]
    T = 200.0
    Pvals = np.arange(1.0,100.0,1.0)
    alpha = []
    for P in Pvals:
        alpha.append(nh3_kd.alpha(f,T,P,X_partial,P_dict,otherPar))
    plt.plot(Pvals,alpha)

    plt.figure('vsT')
    f = [0.5]
    P = 200.0
    Tvals = np.arange(100.0,1000.0,10.0)
    alpha = []
    for T in Tvals:
        alpha.append(nh3_kd.alpha(f,T,P,X_partial,P_dict,otherPar))
    plt.plot(Tvals,alpha)


f = [2.0, 8.0,20.0,40]
fclr = ['b','r','g','m']
ltyp = ['--','-','-','--']

a_kd = []
a_bg = []
a_sjs = []
a_sjsd = []
a_dbs = []
a_dbs_sjs = []
usenh3 = {'kd':True,'bg':True,'sjs':True,'sjsd':True,'dbs':True,'dbs_sjs':True}
for i in range(len(P_Jup)):
    T = T_Jup[i]
    P = P_Jup[i]
    X_partial = X_Jup[i]
    if P<0.01:
        continue
    if usenh3['kd']:
        a_kd.append(nh3_kd.alpha(f,T,P,X_partial,P_dict,otherPar))
    if usenh3['bg']:
        a_bg.append(nh3_bg.alpha(f,T,P,X_partial,P_dict,otherPar))
    if usenh3['sjs']:
        a_sjs.append(nh3_sjs.alpha(f,T,P,X_partial,P_dict,otherPar))
    if usenh3['sjsd']:
        a_sjsd.append(nh3_sjsd.alpha(f,T,P,X_partial,P_dict,otherPar))
    if usenh3['dbs']:
        a_dbs.append(nh3_dbs.alpha(f,T,P,X_partial,P_dict,otherPar))
    if usenh3['dbs_sjs']:
        a_dbs_sjs.append(nh3_dbs_sjs.alpha(f,T,P,X_partial,P_dict,otherPar))
    #print '%5.0f\t%4.0f\t%.3f\t%.3f' % (T,P,a_kd[i][0],a_bg[i][0])
a_kd = np.array(a_kd)
a_bg = np.array(a_bg)
a_sjs = np.array(a_sjs)
a_sjsd = np.array(a_sjsd)
a_dbs = np.array(a_dbs)
a_dbs_sjs = np.array(a_dbs_sjs)
plt.figure('samples')
for i in range(len(f)):
    ########------------KD------------########
    if usenh3['kd']:
        if i==1:
            s = 'kd'
        else:
            s = None
        ll = 'b'+ltyp[i]
        plt.plot(P_Jup,a_kd[:,i],ll,linewidth=2,label=s)
        Y = a_kd[10,i]
    ########------------BG------------########
    if usenh3['bg']:
        if i==1:
            s = 'b-g'
        else:
            s = None
        ll = 'r'+ltyp[i]
        plt.plot(P_Jup,a_bg[:,i],ll,linewidth=2,label=s)
        Y = a_bg[10,i]
    ########------------SJS------------########
    if usenh3['sjs']:
        if i==1:
            s = 'ts+jj/ps'
        else:
            s = None
        ll = 'k'+ltyp[i]
        plt.plot(P_Jup,a_sjs[:,i],ll,linewidth=2,label=s)
        Y = a_sjs[10,i]
    ########------------SJSD------------########
    if usenh3['sjsd']:
        if i==1:
            s = 'kd+(ts+jj/ps)'
        else:
            s = None
        ll = 'g'+ltyp[i]
        plt.plot(P_Jup,a_sjsd[:,i],ll,linewidth=2,label=s)
        Y = a_sjsd[10,i]
    ########------------DBS------------########
    if usenh3['dbs']:
        if i==1:
            s = 'ab/ps'
        else:
            s = None
        ll = 'm'+ltyp[i]
        plt.plot(P_Jup,a_dbs[:,i],ll,linewidth=2,label=s)
        Y = a_dbs[10,i]
    ########-----------DBSSJS-----------########
    if usenh3['dbs_sjs']:
        if i==1:
            s = 'ab/ps+(ts+jj/ps)'
        else:
            s = None
        ll = 'c'+ltyp[i]
        plt.plot(P_Jup,a_dbs_sjs[:,i],ll,linewidth=2,label=s)
        Y = a_dbs_sjs[10,i]
    s = '%.0f GHz' % (f[i])
    plt.text(P_Jup[10],Y,s)
    print s
plt.xscale('log')
plt.yscale('log')
#plt.legend()
plt.xlabel('Pressure [bars]')
plt.ylabel(r'$\alpha$ [dB/km]')
s = 'jupiter.paulSolar'
plt.title(s)

if False:
    f = [1.4,4.0, 8.0,20.0,40.0]
    fclr = ['c','b','r','g']
    a_kd = []
    a_bg = []
    a_sjs = []
    a_sjsd = []
    for i in range(len(P_Jup)):
        T = T_Jup[i]
        P = P_Jup[i]
        X_partial = X_Jup
        a_kd.append(nh3_kd.alpha(f,T,P,X_partial,P_dict,otherPar))
        a_bg.append(nh3_bg.alpha(f,T,P,X_partial,P_dict,otherPar))
        a_sjs.append(nh3_sjs.alpha(f,T,P,X_partial,P_dict,otherPar))
        a_sjsd.append(nh3_sjsd.alpha(f,T,P,X_partial,P_dict,otherPar))
        fp.write('%f\t%f\t%f\t%f\t%f\t%s\t%s\n' % (P,T,X_partial[0],X_partial[1],X_partial[2],str(a_kd[i]),str(a_bg[i])))
    a_kd = np.array(a_kd)
    a_bg = np.array(a_bg)
    a_sjs = np.array(a_sjs)
    a_sjsd = np.array(a_sjsd)
    plt.figure('ssampless')
    for i in range(len(f)):
        plt.plot(P_Jup,a_kd[:,i],fclr[i])
        plt.plot(P_Jup,a_bg[:,i],fclr[i]+'--')
        plt.plot(P_Jup,a_sjs[:,i],fclr[i]+':')
        plt.plot(P_Jup,a_sjs[:,i],'k')
    plt.xscale('log')
    plt.xlabel('Pressure [bars]')
    plt.ylabel(r'$\alpha$ [dB/km]')
    fp.close()
