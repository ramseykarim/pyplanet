import math
import string
import matplotlib.pyplot as plt
import numpy as np
import utils
import sys
import os
import os.path
import config as pcfg

###local imports
import raypath as ray
try:
    pyPlanetPath = os.getenv('PYPLANETPATH')
except:
    print 'No PYPLANETPATH environment variable'
    pyPlanetPath = './'
sys.path.append(os.path.join(pyPlanetPath,'constituents/'))
import properties
import regrid

planetDictionary = {'Jupiter':0,'Saturn':1,'Uranus':2,'Neptune':3}

class Atmosphere:
    def __init__(self,planet,config='config.par',path=None,log=None,verbosity=False,plot=True):
        """reads/computes atmospheres.  This should return:
               self.gas
               self.cloud
               self.layerProperty
            on the appropriate grid
            Note that the research is in the input files and modifying the tweak modules
            All of the default config parameters are hard-coded here:  see __init__, setConfig, showConfig."""

        planet = string.capitalize(planet)
        self.planet = planet
        self.verbosity=verbosity
        self.plot=plot
        self.logFile = utils.setupLogFile(log)
        self.batch = False
        
        print '\n---Atmosphere of %s---' % (planet)
        if type(config) == str:
            config = pcfg.planetConfig(self.planet,configFile=config,log=log,verbosity=verbosity)
        self.config = config

        ###Create function dictionaries
        self.gasGen = {}
        self.gasGen['read'] = self.readGas
        self.gasGen['compute'] = self.computeGas
        self.cloudGen = {}
        self.cloudGen['read'] = self.readCloud
        self.cloudGen['compute'] = self.computeCloud
        self.propGen = {}
        self.propGen['read'] = self.readProp
        self.propGen['compute'] = self.computeProp

        self.tweak_function = lambda gas, cloud, C, Cl: "There was no tweak defined!!"

        print 'Planet '+self.planet
        if self.config.gasType == 'read':  # this assumes that cloudType is then also 'read'
            utils.log(self.logFile,'\tReading from: '+self.config.path,True)
            utils.log(self.logFile,'\tAtmosphere file:  '+self.config.gasFile,True)
            utils.log(self.logFile,'\tCloud file:  '+self.config.cloudFile,True)
        if verbosity:
            print self.config.show()
        # print '\tIf in interactive mode, change any parameters in setpar() before run()'

    def run(self,Pmin=None,Pmax=None,regridType=None,
            gasType=None,cloudType=None,otherType=None,
            plot=None,verbosity=None,
            tweak=True, tweak_fn=None):
        """This is the standard pipeline"""
        ###Set run defaults
        if Pmin == None:
            Pmin = self.config.pmin
        if Pmax == None:
            Pmax = self.config.pmax
        if regridType == None:
            regridType = self.config.regridType
        if verbosity==None:
            verbosity = self.verbosity
        if plot==None:
            plot = self.plot
        if gasType==None:
            gasType = self.config.gasType
        if cloudType==None:
            cloudType = self.config.cloudType
        if otherType==None:
            otherType = self.config.otherType
        self.nAtm = 0

        ### Generate gas profile (gasType is 'read' or 'compute')
        if not self.gasGen.has_key(gasType):
            print 'Error:  No such gasType: '+gasType
            return 0
        else:
            self.gasGen[gasType](verbosity=verbosity)

        if not self.batch:
            ### Generate cloud profile (cloudType is 'read' or 'compute')
            if not self.cloudGen.has_key(cloudType):
                print 'Error:  No such cloudType: '+cloudType
                return 0
            else:
                self.cloudGen[cloudType](verbosity=verbosity)

            if tweak:  # This loads and calls the module 'tweakFile'
                if tweak_fn is None:
                    self.tweak_function = self.obtain_tweak_module()
                else:
                    self.tweak_function = tweak_fn
                self.tweakAtm()

            ### Compute other parameters that are needed
            if not self.propGen.has_key(self.config.otherType):
                print 'Error:  no such otherTpe: '+otherType
                return 0
            else:
                self.propGen[otherType](verbosity=verbosity)

            ### Put onto common grid
            regridded = regrid.regrid(self,regridType=regridType,Pmin=Pmin,Pmax=Pmax)
            self.nAtm = len(self.gas[0])

            angularDiameter = 2.0*math.atan(self.layerProperty[self.config.LP['R']][0]/self.config.distance)
            if verbosity:
                print 'angular radius = %f arcsec' % ((180.0/np.pi)*3600.0*angularDiameter/2.0)
        
            ### Plot data
            if plot:
                self.plotTP()
                self.plotGas()
                self.plotCloud()
                self.plotProp()

        return self.nAtm

    def getval(self,val=None,vtype='gas'):
        """Returns one of the constituent or cloud profiles"""
        if val == None:
            print "Usage:  getVal('v',['gas'/'cloud'/'other'])"
            print "    'gas' is default"
            print 'These are the gas values:'
            print self.config.C
            print 'These are the cloud values:'
            print self.config.Cl
            print 'These are the layerProperty (other) values:'
            print self.config.LP
            return None
        v = string.upper(val)
        vt = string.lower(vtype)
        rv = 0
        found = False
        
        if vt=='gas':
            if self.config.C.has_key(v):
                rv = self.gas[self.config.C[v]]
                print 'Found '+val+' in gas'
                found = True
        elif vt=='cloud':
            if self.config.Cl.has_key(v):
                rv = self.cloud[self.config.Cl[v]]
                print 'Found '+val+' in cloud'
                found = True
        elif vt=='other':
            if self.config.LP.has_key(v):
                rv = self.layerProperty[self.config.LP[v]]
                print 'Found '+val+' in layerProperty'
                found = True
        print vt+'  '+val,
        if found:
            print ':  found'
        else:
            print ':  not found'
            
        return rv

    def plotTP(self,plot='auto'):
        """Plot the T-P profile"""
        if plot=='auto':
            plt.figure(planetDictionary[self.planet])
        plt.title(self.planet+':  T-P profile')
        plt.loglog(self.gas[self.config.C['T']],self.gas[self.config.C['P']])
        v = list(plt.axis())
        v[2] = 100.0*math.ceil(self.gas[self.config.C['P']][-1]/100.0)
        v[3] = 1.0E-7*math.ceil(self.gas[self.config.C['P']][0]/1E-7)
        plt.axis(v)
        plt.ylabel('P [bars]')
        plt.xlabel('T [K]')

    def plotCloud(self,dontPlot=['Z','P','T','DZ'],plot='auto'):
        """Plots the clouds"""
        if plot=='auto':
            plt.figure(planetDictionary[self.planet]+11)
        plt.title(self.planet+': clouds')
        for cloud in self.config.Cl:
            present, cl = self.__isPresent__(self.cloud[self.config.Cl[cloud]])
            if cloud in dontPlot or not present:
                continue
            plt.loglog(cl,self.cloud[self.config.Cl['P']], label=cloud)
        v = list(plt.axis())
        if v[0] < 1E-10:
            v[0] = 1E-10
        v[2] = 100.0*math.ceil(self.gas[self.config.C['P']][-1]/100.0)
        v[3] = 1.0E-7*math.ceil(self.gas[self.config.C['P']][0]/1E-7)
        plt.axis(v)
        plt.ylabel('P [bars]')
        plt.xlabel(r'Density [g/cm$^3$]')
        plt.legend()

    def plotGas(self,dontPlot=['Z','P','T','DZ'],plot='auto'):
        """Plots the constituents"""
        if plot=='auto':
            plt.figure(planetDictionary[self.planet]+10)
        plt.title(self.planet+': gas')
        for gas in self.config.C:
            present, g = self.__isPresent__(self.gas[self.config.C[gas]])
            if gas in dontPlot or not present:
                continue
            plt.loglog(g,self.gas[self.config.C['P']], label=gas)
        v = list(plt.axis())
        if v[0] < 1E-10:
            v[0] = 1E-10
        v[2] = 100.0*math.ceil(self.gas[self.config.C['P']][-1]/100.0)
        v[3] = 1.0E-7*math.ceil(self.gas[self.config.C['P']][0]/1E-7)
        plt.axis(v)
        plt.ylabel('P [bars]')
        plt.xlabel('Fractional Abundance')
        plt.legend()

    def plotProp(self,dontPlot=['Z','P','T'],plot='auto'):
        if plot=='auto':
            plt.figure(planetDictionary[self.planet]+12)
        plt.title(self.planet+': other')
        for other in self.config.LP:
            present, g = self.__isPresent__(self.layerProperty[self.config.LP[other]])
            if other in dontPlot or not present:
                continue
            plt.loglog(g,self.gas[self.config.C['P']], label=other)
        v = list(plt.axis())
        if v[0] < 1E-10:
            v[0] = 1E-10
        v[2] = 100.0*math.ceil(self.gas[self.config.C['P']][-1]/100.0)
        v[3] = 1.0E-7*math.ceil(self.gas[self.config.C['P']][0]/1E-7)
        plt.axis(v)
        plt.ylabel('P [bars]')
        plt.xlabel('Property')
        plt.legend()

    def readGas(self,gasFile=None,numHeaderLines=None,Cdict=None,verbosity=False):
        """Reads gas profile file as self.gas"""

        if gasFile is None:
            gasFile=self.config.gasFile
        if Cdict is None:
            Cdict = self.config.C
        if gasFile == 'batch':
            self.batch = True
            return 0
        else:
            self.batch = False
        gasFile = os.path.join(self.config.path,gasFile)
        if numHeaderLines is None:
            numHeaderLines = self.config.gasFileHdr

        print 'Reading '+gasFile+'  (Header:  '+str(numHeaderLines)+')'
        self.gas = []
        print '\tUsing atmsopheric component:  ',
        for k in Cdict:
            print k+'  ',
            self.gas.append([])
        self.nConstituent = len(Cdict)
        try:
            fp = open(gasFile,"r")
        except IOError:
            print gasFile+' was not found - returning no gas profile\n\n'
            raise IOError
        i=0
        lineno = 0
        pastHeader = False
        print ' '
        for line in fp:
            lineno+=1
            data = line.split()
            skip_row = line[0]=='!' or len(data) < 4 or lineno<=numHeaderLines
            ######put in something to eliminate the need for numHeaderLines..........
            if skip_row:
                if verbosity:
                    print '\tHEADER: '+line,
                continue
            pastHeader = True
            if len(line) < 3:
                if pastHeader:
                    break
                else:
                    continue
            if i==0:
                print '\tCompiling '+str(self.nConstituent)+' components (data-file: '+str(len(data))+')'
            for n in range(self.nConstituent): ## Initialize all of the constituents to 0.0
                self.gas[n].append(0.0)
            columns = Cdict.values()
            columns.sort()
            for n,v in enumerate(columns): ## ...now read in data
                if v < len(data):
                    self.gas[n][i] = float(data[v])
                else:
                    self.gas[n][i] = 0.0
            if (self.nConstituent-n)!=1:
                print 'line '+str(i+1)+' has incorrect number of components'
                print '['+line+']'
            i+=1
        self.nGas = i
        ###Redo the constituent dictionary for the self.gas index positions
        nid,sk = utils.invertDictionary(Cdict)
        for i,k in enumerate(sk):
            Cdict[nid[k]] = i
        self.config.C = Cdict
        
        print '\tRead '+str(self.nGas)+' lines'
        fp.close()
        self.gas = np.array(self.gas)
        ### Check that P is monotonically increasing
        monotonic = np.all(np.diff(self.gas[self.config.C['P']])>0.0)
        if not monotonic:
            self.gas = np.fliplr(self.gas)
            monotonic = np.all(np.diff(self.gas[self.config.C['P']])>0.0)
        if not monotonic:
            print "Error in "+gasFile+".  Pressure not monotonically increasing"

        ### Renormalize so that deepest z is 0
        zDeep = self.gas[self.config.C['Z']][-1]
        if abs(zDeep) > 1E-6:
            for i in range(len(self.gas[self.config.C['Z']])):
                self.gas[self.config.C['Z']][i]-=zDeep
            print 'Deepest levels set from %.1f to 0' % (zDeep)

        ### Set dz
        dz = np.abs(np.diff(self.gas[self.config.C['Z']]))*1.0E5  #convert from km to cm (so no unit checking!!!)
        self.gas[self.config.C['DZ']] = np.append(np.array([0.0]),dz)
        return self.nGas

    def writeGas(self,outputFile='gas.dat'):
        fp = open(outputFile,'w')
        gdat = np.transpose(self.gas)
        fp.write(str(len(gdat))+'\n')
        for data in gdat:
            s = ''
            for i in range(11):
                s+= (str(data[i])+'\t')
            s+='\n'
            fp.write(s)
        fp.close()
                

    def readCloud(self,cloudFile=None,numHeaderLines=None,Cldict=None,verbosity=False):
        """Reads in cloud data if we have it..."""

        if cloudFile == None:
            cloudFile=self.config.cloudFile
        if Cldict == None:
            Cldict = self.config.Cl
        if cloudFile == 'batch':
            self.batch = True
            return 0
        else:
            self.batch = False
        cloudFile = os.path.join(self.config.path,cloudFile)
        if numHeaderLines is None:
            numHeaderLines = self.config.cloudFileHdr

        print 'Reading '+cloudFile+'  (Header:  '+str(numHeaderLines)+')'
        self.cloud = []
        print '\tUsing cloud component:  ',
        for k in Cldict:
            print k+'  ',
            self.cloud.append([])
        self.nParticulate = len(Cldict)
        try:
            fp = open(cloudFile,"r")
        except IOError:
            print cloudFile+' was not found - returning no clouds\n\n'
            raise IOError
        i=0
        lineno = 0
        pastHeader = False
        print ' '
        for line in fp:
            lineno+=1
            if line[0]=='!' or lineno<=numHeaderLines:
                if verbosity:
                    print '\tHEADER: '+line,
                continue
            pastHeader = True
            if len(line) < 3:
                if pastHeader:
                    break
                else:
                    continue
            data = line.split()
            if i==0:
                print '\tReading '+str(self.nParticulate)+' cloud particles (data-file: '+str(len(data))+')'
            for n in range(self.nParticulate): ## Initialize all of the particulates to 0.0
                self.cloud[n].append(0.0)
            columns = Cldict.values()
            columns.sort()
            for n,v in enumerate(columns): ## ...now read in data
                if v < len(data):
                    self.cloud[n][i] = float(data[v])
                else:
                    self.cloud[n][i] = 0.0
            if (self.nParticulate-n)!=1:
                print 'line '+str(i+1)+' has incorrect number of particulates'
                print '['+line+']'
            i+=1
        self.nCloud = i
        ###Redo the particulate dictionary for the self.cloud index positions
        nid,sk = utils.invertDictionary(Cldict)
        for i,k in enumerate(sk):
            Cldict[nid[k]] = i
        self.config.Cl = Cldict
        
        print '\tRead '+str(self.nCloud)+' lines'
        fp.close()
        self.cloud = np.array(self.cloud)
        ### Check that P is monotonically increasing
        monotonic = np.all(np.diff(self.cloud[self.config.Cl['P']])>0.0)
        if not monotonic:
            self.cloud = np.fliplr(self.cloud)
            monotonic = np.all(np.diff(self.cloud[self.config.Cl['P']])>0.0)
        if not monotonic:
            print "Error in "+cloudFile+".  Pressure not monotonically increasing"
        self.cloud[self.config.Cl['DZ']] = np.abs(np.append(np.diff(self.cloud[self.config.Cl['Z']]),0.0))
        return self.nCloud

    def readProp(self,otherFile=None,verbosity=False):
        """Reads in other property data if we have it..."""
        if otherFile == None:
            otherFile = self.config.otherFile
        otherFile = os.path.join(self.config.path,otherFile)
        print "readProp currently doesn't read in any other file..."
        print "it will probably not be needed, since it will get computed in computeProp"
        return False

    def tweakAtm(self):
        """Tweaks the atmosphere data..."""
        if self.batch:
            print 'Defer tweaking'
            return 0
        nAtm = len(self.gas[self.config.C['P']])
        if nAtm - self.nGas != 0:
            print 'Error in number of layers - check it out ('+str(nAtm)+'/'+str(self.nGas)+')'
            print 'Returned from tweakAtm'
            return 0

        # Run module then log
        self.tweakComment, self.gas, self.cloud = self.tweak_function(self.gas, self.cloud,
                                                                      self.config.C, self.config.Cl)
        print '---tweakComment'
        print self.tweakComment
        print '---'
        utils.log(self.logFile,self.tweakComment,False)
        _tf = os.path.join(self.config.path,self.config.tweakFile+'.py')
        _tp = open(_tf,'r')
        dt = _tp.read()
        utils.log(self.logFile,'======================'+_tf+'=====================',False)
        utils.log(self.logFile,dt,False)
        utils.log(self.logFile,'====================================================================',False)
        _tp.close()

        return nAtm

    def obtain_tweak_module(self):
        """
        Author: Ramsey Karim
        Returns the modifying function
        This function is meant to be overwritten if the tweak function
        is to be recreated multiple times in a single instance.
        -------

        """
        # Import tweakFile
        sys.path.append(self.config.path)
        try:
            __import__(self.config.tweakFile)
            tweakModule = sys.modules[self.config.tweakFile]
            return tweakModule.modify
        except SyntaxError:
            utils.log(self.logFile,"Syntax Error:  check "+self.config.tweakFile,True)
            raise SyntaxError
        except:
            utils.log(self.logFile,"non-syntax tweakAtm error",True)
            raise RuntimeError

    def computeProp(self,verbosity=False):
        """This module computes derived atmospheric properties (makes self.layerProperty)"""
        if self.batch:
            print 'Defer computing properties'
            return 0
        nAtm = len(self.gas[self.config.C['P']])
        self.layerProperty = []
        for op in self.config.LP:
            self.layerProperty.append([])
        zOffset = 0.0
        iOffset = 0
        psep = 1.0E6
        for i,zv in enumerate(self.gas[self.config.C['Z']]):     # find the nearest z value at p_ref
            P = self.gas[self.config.C['P']][i]
            if abs(P-self.config.p_ref) < psep:
                psep = abs(P-self.config.p_ref)
                iOffset = i
        zOffset = self.gas[self.config.C['Z']][iOffset]
        z_at_p_ref = self.config.Req
        if verbosity:
            print "z,P offset:  ",zOffset,self.gas[self.config.C['P']][iOffset]
        
        for i,zv in enumerate(self.gas[self.config.C['Z']]):
            T = self.gas[self.config.C['T']][i]
            P = self.gas[self.config.C['P']][i]
            self.layerProperty[self.config.LP['P']].append(P)
            self.layerProperty[self.config.LP['Z']].append(zv)
            rr = z_at_p_ref + zv - zOffset
            self.layerProperty[self.config.LP['R']].append(rr)   # note that this is the "actual" z along equator  referenced to planet center (aka radius)
            ###set mean amu
            amulyr = 0.0
            for key in properties.amu:
                if key in self.config.C:
                    amulyr+=properties.amu[key]*self.gas[self.config.C[key]][i]
            self.layerProperty[self.config.LP['AMU']].append(amulyr)
            ###set GM pre-calc (normalized further down) and get lapse rate
            if not i:
                self.layerProperty[self.config.LP['GM']].append(0.0)
                self.layerProperty[self.config.LP['LAPSE']].append(0.0)
            else:
                rho = (amulyr*P)/(properties.R*T)
                dr = abs(zv - self.gas[self.config.C['Z']][i-1])
                dV = 4.0*math.pi*(rr**2)*dr
                dM = 1.0e11*rho*dV
                GdM = self.layerProperty[self.config.LP['GM']][i-1] + properties.GravConst*dM    # in km3/s2
                self.layerProperty[self.config.LP['GM']].append( GdM )  # mass added as you make way into atmosphere by radius r (times G)
                dT = abs(T - self.gas[self.config.C['T']][i-1])
                self.layerProperty[self.config.LP['LAPSE']].append(dT/dr)
            ###set refractivity and index of refraction
            refrlyr = 0.0
            for key in properties.refractivity:
                if key in self.config.C:
                    refrlyr+=properties.refractivity[key]*self.gas[self.config.C[key]][i]
            refrlyr = refrlyr*P*(293.0/T)
            self.layerProperty[self.config.LP['REFR']].append(refrlyr)
            nlyr = refrlyr/1.0E6 + 1.0
            self.layerProperty[self.config.LP['N']].append(nlyr)
 
        ###Now need to normalize GM to planet and calculate scale height (H)
        GMnorm = self.layerProperty[self.config.LP['GM']][iOffset]  # G*(Mass added by p_ref)
        for i, mv in enumerate(self.layerProperty[self.config.LP['GM']]):
            gm = self.config.GM_ref - (mv-GMnorm)
            self.layerProperty[self.config.LP['GM']][i] = gm
            little_g = gm/self.layerProperty[self.config.LP['R']][i]**2
            m_bar = self.layerProperty[self.config.LP['AMU']][i]
            T = self.gas[self.config.C['T']][i]
            self.layerProperty[self.config.LP['H']].append( (properties.R*T)/(little_g*m_bar)/1000.0 )
            self.layerProperty[self.config.LP['g']].append( little_g )
        self.layerProperty = np.array(self.layerProperty)
            
    def computeCloud(self,verbosity=False):
        """This computes cloud stuff"""
        print 'computeCloud does nothing yet.  This probably wont do anything since gas/cloud/other will get computed in the same tcm'
    def computeGas(self,verbosity=False):
        """Computes an atmosphere given stuff"""
        print 'computeGas does nothing yet.  This probably wont do anything since gas/cloud/other will get computed in the same tcm'
        print 'This will probably just call a external tcm module'

    def __isPresent__(self,c,tiny=1.0E-30):
        """This checks to see if a constituent is there and sets 0.0 or negative values to tiny.  This is generally for log plotting"""
        vsum = 0.0
        present = False
        vnew = []
        for v in c:
            vsum+=v
            if v==0.0:
                vnew.append(tiny)
            elif v<0.0:
                vnew.append(tiny)
            else:
                vnew.append(v)
                present = True
                
        return present, vnew
