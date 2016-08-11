### This is the file to calculate the radiometric properties of the planets
import math
import matplotlib.pyplot as plt
import numpy as np
import scipy.special as ss
import utils
import sys
import raypath as ray
import os.path

class brightness():
    def __init__(self,log=None,verbose=False,plot=False):
        """This calculates the brightness temperature of the planets.
           It must be used with atmosphere and alpha"""
        self.verbose = verbose
        self.plot = plot
        self.log = utils.setupLogFile(log)
        print '\n---Brightness---\n'

    def layerAbsorption(self,freqs,atm,alpha,verbose=None,plot=None):
        self.layerAlpha = self.__layerAbsorp__(freqs,atm,alpha,verbose=verbose,plot=plot)
        if plot:
            P = atm.gas[atm.config.C['P']]
            plt.figure('alpha')
            for i,f in enumerate(freqs):
                label = '%.1f GHz' % (f)
                plt.loglog(self.layerAlpha[i],P,label=label)
            v = list(plt.axis())
            v[2] = 100.0*math.ceil(atm.gas[atm.config.C['P']][-1]/100.0)
            v[3] = 1.0E-7*math.ceil(atm.gas[atm.config.C['P']][0]/1E-7)
            plt.axis(v)
            plt.xlabel(utils.alphaUnits)
            plt.ylabel('P [bars]')
            plt.legend()
            #lgd=plt.legend(loc='upper left',bbox_to_anchor=(1,1))
            #lgd.set_visible(True)  # This is just to remind me...
        
    def __layerAbsorp__(self,freqs,atm,alpha,verbose=None,plot=None):
        alphaUnits = 'invcm'
        if verbose == None:
            verbose = self.verbose
        if plot == None:
            plot = self.plot
        self.freqs = freqs
        numLayers = len(atm.gas[0])
        layerAlp=[]
        P = atm.gas[atm.config.C['P']]
        T = atm.gas[atm.config.C['T']]
        utils.log(self.log,'%d layers' % (numLayers),True)
        #print '\t Computing absorption in layers...'
        for layer in range(numLayers):
            print '\r\tAbsorption in layer %d   ' % (layer+1),
            sys.stdout.flush()
            layerAlp.append(alpha.getAlpha(freqs,T[layer],P[layer],atm.gas[:,layer],atm.config.C,atm.cloud[:,layer],atm.config.Cl,units=alphaUnits,verbose=verbose))
        layerAlp = np.array(layerAlp).transpose()
        print ' '
        return layerAlp 

    def single(self, freqs, atm, b, alpha, orientation=None, taulimit=20.0, plot= None, verbose = None, discAverage = False, normW4plot=True):
        """This computes the brightness temperature along one ray path"""

        if verbose is None:
            verbose = self.verbose
        if plot is None:
            plot = self.plot
        # get path lengths (ds_layer) vs layer number (num_layer) - currently frequency independent refractivity
        self.path = ray.compute_ds(atm,b,orientation,gtype=None,verbose=verbose,plot=plot)
        if self.path.ds == None:
            print 'Off planet'
            self.Tb = []
            for j in range(len(freqs)):
                self.Tb.append(utils.T_cmb)
            return self.Tb
        
        # these profiles are saved
        self.tau = []
        self.W = []
        self.Tb_lyr = []
        
        # temporary arrays
        taus = []
        Tbs = []
        Ws = []

        # initialize
        for j in range(len(freqs)):
            taus.append(0.0)
            Tbs.append(0.0)
            Ws.append(0.0)
        self.tau.append(taus)
        self.W.append(Ws)
        self.Tb_lyr.append(Tbs)

        if alpha.config.Doppler:
            P = atm.gas[atm.config.C['P']]
            T = atm.gas[atm.config.C['T']]
            alphaUnits = 'invcm'
            #--debug--#self.debugDoppler = []
            print ''
        
        for i in range( len(self.path.ds)-1 ):
            ds = self.path.ds[i]*utils.Units[utils.processingAtmLayerUnit]/utils.Units['cm']
            taus = []
            Ws = []
            Tbs = []
            ii = self.path.layer4ds[i]
            ii1= self.path.layer4ds[i+1]
            
            for j,f in enumerate(freqs):
                if not alpha.config.Doppler:
                    a1 = self.layerAlpha[j][ii1]
                    a0 = self.layerAlpha[j][ii]
                else:
                    fshifted=[[f/self.path.doppler[i]],[f/self.path.doppler[i+1]]]
                    #--debug--#self.debugDoppler.append(fshifted[0][0])
                    print '\rdoppler corrected frequency at layer',i,
                    a1 = alpha.getAlpha(fshifted[0],T[ii1],P[ii1],atm.gas[:,ii1],atm.config.C,atm.cloud[:,ii1],atm.config.Cl,units=alphaUnits,verbose=False)
                    a0 = alpha.getAlpha(fshifted[1],T[ii],P[ii],atm.gas[:,ii],atm.config.C,atm.cloud[:,ii],atm.config.Cl,units=alphaUnits,verbose=False)
                dtau = (a0 + a1)*ds/2.0
                taus.append(self.tau[i][j] + dtau)         # this is tau_(i+1)
                T1 = atm.gas[atm.config.C['T']][ii1]
                T0 = atm.gas[atm.config.C['T']][ii]
                
                if discAverage==True:
                    Ws.append( 2.0*a1*ss.expn(2,taus[j]) )     # this is W_(i+1) for disc average
                #    dTb = ( T1*ss.expn(2,taus[j])/scriptR(T1,freqs[j]) + T0*ss.expn(2,self.tau[ii][j])/scriptR(T0,freqs[j]) )*dtau
                #    Tbs.append( self.Tb_lyr[i][j] + dTb )
                else:
                    Ws.append( a1*math.exp(-taus[j]) )                  # this is W_(i+1) for non disc average
                dTb = ( T1*Ws[j]/scriptR(T1,freqs[j]) + T0*self.W[i][j]/scriptR(T0,freqs[j]) )*ds/2.0
                Tbs.append( self.Tb_lyr[i][j] + dTb)
            self.tau.append(taus)
            self.W.append(Ws)
            self.Tb_lyr.append(Tbs)
        print ''

        # final spectrum
        self.Tb = []
        for j in range(len(freqs)):
            top_Tb_lyr = self.Tb_lyr[-1][j]
            if top_Tb_lyr < utils.T_cmb:
                top_Tb_lyr = utils.T_cmb
            self.Tb.append(top_Tb_lyr)
        self.tau = np.array(self.tau).transpose()
        self.W = np.array(self.W).transpose()
        self.Tb_lyr = np.array(self.Tb_lyr).transpose()

        try:
            if plot:
                # save a local copy of
                self.P = atm.gas[atm.config.C['P']][0:len(self.W[0])]
                self.z = atm.gas[atm.config.C['Z']][0:len(self.W[0])]
                    
                #####-----Weigthing functions
                plt.figure('radtran')
                plt.subplot(121)
                for i,f in enumerate(freqs):
                    #label=r'$\tau$: %.1f GHz' % (f)
                    #plt.semilogy(self.tau[i],self.P,label=label)
                    if normW4plot:
                        wplot = self.W[i]/np.max(self.W[i])
                    else:
                        wplot = self.W[i]
                    label=r'$W$: %.1f GHz' % (f)
                    label=r'%.1f cm' % (30.0/f)
                    #label=r'%.0f$^o$' % ((180.0/math.pi)*math.asin(b[0]))
                    plt.semilogy(wplot,self.P,label=label,linewidth=3)
                    #label=r'Tlyr$_b$: %.1f GHz' % (f)
                    #plt.semilogy(self.Tb_lyr[i],self.P,label=label)
                plt.legend()
                plt.axis(ymin=100.0*math.ceil(self.P[-1]/100.0), ymax=1.0E-7*math.ceil(self.P[0]/1E-7))
                #plt.xlabel('units')
                plt.ylabel('P [bars]')
                #####-----Alpha
                plt.figure('alpha')
                for i,f in enumerate(freqs):
                    label=r'$\alpha$: %.1f GHz' % (f)
                    label=r'%.1f cm' % (30.0/f)
                    pl = list(self.layerAlpha[i])
                    del pl[0]
                    #delete because alpha is at the layer boundaries, so there are n+1 of them
                    plt.loglog(pl,self.P,label=label)
                plt.legend()
                v = list(plt.axis())
                v[2] = 100.0*math.ceil(self.P[-1]/100.0)
                v[3] = 1.0E-7*math.ceil(self.P[0]/1E-7)
                plt.axis(v)
                #plt.legend()
                #plt.xlabel('units')
                plt.ylabel('P [bars]')
                #####-----Brightness temperature
                plt.figure('brightness')
                lt = '-'
                if (len(self.Tb)==1):
                    lt = 'o'
                plt.plot(freqs,self.Tb,lt)
                plt.xlabel('Frequency [GHz]')
                plt.ylabel('Brightness temperature [K]')
        except:
            print 'Plotting broke'
            
        del taus, Tbs, Ws

        return self.Tb

    def savertm(self,tag=None,path='Output'):
        if tag==None:
            filename = None
        else:
            filename = 'alpha_'+tag+'.out'
        self.saveAlpha(filename,path)
        if tag==None:
            filename = None
        else:
            filename = 'wgt_'+tag+'.out'
        self.saveWeight(filename,path)
        if tag==None:
            filename = None
        else:
            filename = 'tau_'+tag+'.out'
        self.saveTau(filename,path)
        if tag==None:
            filename = None
        else:
            filename = 'tblayer_'+tag+'.out'
        self.saveTblayer(filename,path)

    def saveAlpha(self,filename=None,path='.'):
        if filename == None:
            filename = 'alpha.out'
        os.path.join(path,filename)
        fp = open(filename,'w')
        s = '#P  \tz  \t'
        for f in self.freqs:
            s+='%.2f\t' % (f)
        s+='GHz\n'
        fp.write(s)
        for j in range(len(self.P)):
            s = '%s\t%.2f\t' % (repr(self.P[j]),self.z[j])
            for i in range(len(self.freqs)):
                s+='%s\t' % (repr(self.layerAlpha[i][j]))
            s+='\n'
            fp.write(s)
        s = '%s (%d x %d)' % (filename,i+1,j+1)
        return s

    def saveWeight(self,norm=False,filename=None,path='.'):
        if filename == None:
            filename = 'wgt.out'
        os.path.join(path,filename)
        fp = open(filename,'w')
        s = '#P  \tz  \t'
        for f in self.freqs:
            s+='%.2f\t' % (f)
        s=s.strip()+'GHz\n'
        fp.write(s)
        scale = []
        for i in range(len(self.freqs)):
            if norm:
                scale.append(np.max(self.W[i]))
            else:
                scale.append(1.0)
        for j in range(len(self.P)):
            s = '%s\t%.2f\t' % (repr(self.P[j]),self.z[j])
            for i in range(len(self.freqs)):
                s+='%s\t' % (repr(self.W[i][j]/scale[i]))
            s=s.strip()+'\n'
            fp.write(s)
        s = '%s (%d x %d)' % (filename,i+1,j+1)
        return s

    def saveTau(self,filename=None,path='.'):
        if filename == None:
            filename = 'tau.out'
        os.path.join(path,filename)
        fp = open(filename,'w')
        s = '#P  \tz  \t'
        for f in self.freqs:
            s+='%.2f\t' % (f)
        s+='GHz\n'
        fp.write(s)
        for j in range(len(self.P)):
            s = '%s\t%.2f\t' % (repr(self.P[j]),self.z[j])
            for i in range(len(self.freqs)):
                s+='%s\t' % (repr(self.tau[j]))
            s+='\n'
            fp.write(s)
        s = '%s (%d x %d)' % (filename,i+1,j+1)
        return s

    def saveTblayer(self,filename=None,path='.'):
        if filename == None:
            filename = 'tblayer.out'
        os.path.join(path,filename)
        fp = open(filename,'w')
        s = '#P  \tz  \t'
        for f in self.freqs:
            s+='%.2f\t' % (f)
        s+='GHz\n'
        fp.write(s)
        for j in range(len(self.P)):
            s = '%s\t%.2f\t' % (repr(self.P[j]),self.z[j])
            for i in range(len(self.freqs)):
                s+='%s\t' % (repr(self.Tb_lyr[j]))
            s+='\n'
            fp.write(s)
        s = '%s (%d x %d)' % (filename,i+1,j+1)
        return s
                       
def scriptR(T,freq):
    """See Janssen pg 7"""
    a = (utils.hP*freq*utils.Units[utils.processingFreqUnit]) / (utils.kB*T)
    R = (math.exp(a) - 1.0)/a
    return R
