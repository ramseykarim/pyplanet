import math
import string
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import numpy as np
import utils
import sys
###local imports
import raypath as ray
sys.path.append('constituents')
import properties

planetDictionary = {'Jupiter':0,'Saturn':1,'Uranus':2,'Neptune':3}
commentChars = ['!','#','$','%','&','*']

##########Use all upper case to help getVal##########
###This 'C'onstituent dictionary has atm layer gas data.   Needs to correspond to datafile if not computing.
CJupiter = {'Z':0,'T':1,'P':2,'H2':3,'HE':4,'CH4':5,'NH3':6,'H2O':7,'H2S':8,'SOLN':9,'OTHER':10, 'DZ':11}
CSaturn = {}
CUranus = {}
CNeptune = {'Z':0,'T':1,'P':2,'H2':3,'HE':4,'CH4':5,'NH3':6,'H2O':7,'H2S':8,'SOLN':9,'OTHER':10,\
           'PH3':11,'CO':12,'CO13':13,'HCN':14,'DZ':15}

###This 'Cl'oud dictionary has the cloud parameters.  Needs to correspond to datafile if not computing.
ClJupiter = {}
ClSaturn = {}
ClUranus = {}
ClNeptune = {'Z':0,'T':1,'P':2,'SOLN':3,'H2O':4,'NH4SH':5,'NH3':6,'H2S':7,'CH4':8,'AR':9,'PH3':10, 'DZ':11}

###This 'lyr' dictionary has other atmospheric layer properties (but not alpha or tau)
lyrProp = {'Z':0,'DS':1,'G':2,'MU':3,'REFR':4,'N':5}

class atmosphere:
    def __init__(self,planet,config=None,log=None,verbose=False,plot=False):
        """reads/computes atmospheres.  This should return:
               self.gas
               self.cloud
               self.layerProperty
            on the appropriate grid
            Note that the research is in the input files and modifying the tweak modules"""

        planet = string.capitalize(planet)
        self.planet = planet
        self.verbose=verbose
        self.plot=plot
        self.logFile = utils.setupLogFile(log)
        
        self.LP = lyrProp
        ###Set default input files and Constituent(C)/Particulate(Cl) dictionaries
        if planet=='Neptune':
            self.path = 'Neptune/'
            self.gasFile = 'neptune.paulCO_cloud21_fletcher_best_dry'
            self.cloudFile = 'nepcloud_CO.cloud21_fletcher_best_dry'
            self.otherFile = None
            self.C = CNeptune
            self.Cl= ClNeptune
            self.gasType = 'read'     ### read vs compute
            self.cloudType = 'read'   ###     "
            self.otherType = 'compute'   ###     "
            self.g_ref = 1130.0         #[cm/s^2]
            self.r_ref = 2.45e9         #[cm]
            self.p_ref = 1.0            #[bars]
        elif planet=='Jupiter':
            self.path = 'Jupiter/'
            self.gasFile = 'jupiter.paulSolar'
            self.cloudFile = 'jupiter.paulclSolar'
            self.otherFile = None
            self.C = CJupiter
            self.Cl= ClJupiter
            self.gasType = 'read'       ### read vs compute
            self.cloudType = 'read'     ###     "
            self.otherType = 'compute'  ###     "
            self.g_ref = 2417.0         #[cm/s^2]
            self.r_ref = 1.0e12         #[cm]
            self.p_ref = 1.0            #[bars]
        else:
            self.gasFile = None
            self.cloudFile = None
            self.otherFile = None
            self.path = None
            self.C = None
            self.Cl= None
            print 'No planet values set for '+planet
        self.tweakType = self.planet               ### should be a planet name (default); if not, tweaking gets skipped
        self.regridType = 'z simple lin 1.0 km'    ### see the regrid module to default
        configFile = 'none'
        if config!=None:
            configFile = self.path+config
            self.readConfig(configFile)
        utils.log(self.logFile,planet,False)
        pars = self.dispPar()
        utils.log(self.logFile,configFile,False)
        utils.log(self.logFile,pars,False)

        ###Create function dictionaries
        self.gasGen = {}
        self.gasGen['read'] = self.readGas
        self.gasGen['compute'] = self.computeGas
        self.cloudGen = {}
        self.cloudGen['read'] = self.readCloud
        self.cloudGen['compute'] = self.computeCloud
        self.otherGen = {}
        self.otherGen['read'] = self.readOther
        self.otherGen['compute'] = self.computeOther
        self.tweakAtm = {}
        self.tweakAtm['Neptune'] = self.tweakAtmNeptune
        self.tweakAtm['Jupiter'] = self.tweakAtmJupiter

        print 'Planet '+self.planet
        if self.gasType == 'read':  # this assumes that cloudType is then also 'read'
            utils.log(self.logFile,'\tReading from: '+self.path,True)
            utils.log(self.logFile,'\tAtmosphere file:  '+self.gasFile,True)
            utils.log(self.logFile,'\tCloud file:  '+self.cloudFile,True)
        if verbose:
            print self.dispPar()
        print '\tIf in interactive mode, change any parameters in setpar() before run()'

    def run(self,regridType=None,gasType=None,cloudType=None,otherType=None,tweakType=None,plot=None,verbose=None):
        """This is the standard pipeline"""
        if regridType == None:
            regridType = self.regridType
        if verbose==None:
            verbose = self.verbose
        if plot==None:
            plot = self.plot
        if gasType==None:
            gasType = self.gasType
        if cloudType==None:
            cloudType = self.cloudType
        if otherType==None:
            otherType = self.otherType
        if tweakType==None:
            tweakType = self.tweakType

        ### Generate gas profile (gasType is 'read' or 'compute')
        if not self.gasGen.has_key(gasType):
            print 'Error:  No such gasType: '+gasType
            return 0
        else:
            self.gasGen[gasType](verbose=verbose)

        ### Generate cloud profile (cloudType is 'read' or 'compute')
        if not self.cloudGen.has_key(cloudType):
            print 'Error:  No such cloudType: '+cloudType
            return 0
        else:
            self.cloudGen[cloudType](verbose=verbose)

        ### Compute other parameters that are needed
        if not self.otherGen.has_key(self.otherType):
            print 'Error:  no such otherTYpe: '+otherType
            return 0
        else:
            self.otherGen[otherType](verbose=verbose)

        ### Tweak atmosphere (tweakType is PlanetName or 'none')
        if self.tweakAtm.has_key(self.tweakType):
            self.tweakAtm[self.planet]()
        
        self.regrid(regridType=regridType)
        self.nAtm = len(self.gas[0])

        ### Plot data
        if plot:
            self.plotTP()
            self.plotGas()
            self.plotCloud()
            self.plotOther()
        return self.nAtm

    def readConfig(self,configFile):
        nTokMod = 0
        try:
            fp = open(configFile,'r')
        except IOError:
            print configFile+' not found.  Using defaults.'
            return None
        for line in fp:
            data = line.split()
            tok = data[0].lower()
            if tok[0] in commentChars or len(line)<8:
                continue
            nTokMod += 1 
            if tok == 'gasfile':
                self.gasFile = data[1]
            elif tok == 'cloudfile':
                self.cloudFile = data[1]
            elif tok == 'otherfile':
                self.otherFile = data[1]
            elif tok == 'constituents':
                for i,w in enumerate(data):
                    if i==0:
                        continue
                    w = w.upper()
                    self.C[w] = i-1
            elif tok == 'clouds':
                for i,w in enumerate(data):
                    if i==0:
                        continue
                    w = w.upper()
                    self.Cl[w] = i-1
            elif tok == 'gastype':
                self.gasType = data[1]
            elif tok == 'cloudtype':
                self.cloudType = data[1]
            elif tok == 'othertype':
                self.otherType = data[1]
            elif tok == 'tweaktype':
                self.tweakType = data[1]
            elif tok == 'regridtype':
                self.regridType = data[1]
            elif tok == 'g_ref':  # in cm/s^2
                self.g_ref = float(data[1])
            elif tok == 'r_ref':  # in cm
                self.r_ref = float(data[1])
            elif tok == 'p_ref':  # in bars
                self.p_ref = float(data[1])
            else:
                print 'Unknown token:  ',tok
                nTokMod-=1
        fp.close()
        return nTokMod

    def getval(self,val=None,vtype='gas'):
        """Returns one of the constituent or cloud profiles"""
        if val == None:
            print "Usage:  getVal('v',['gas'/'cloud'])"
            print "    'gas' is default"
            print 'These are the gas values:'
            print self.C
            print 'These are the cloud values:'
            print self.Cl
            return None
        v = string.upper(val)
        vt = string.lower(vtype)
        rv = 0
        found = False
        
        if vt=='gas':
            if self.C.has_key(v):
                rv = self.gas[self.C[v]]
                found = True
        elif vtype=='cloud':
            if self.C.has_key(v):
                rv = self.cloud[self.Cl[v]]
                found = True
        print vt+'  '+val,
        if found:
            print ':  found'
        else:
            print ':  not found'
            
        return rv

    def sp(self,gasType=None,cloudType=None,gasFile=None,cloudFile=None,tweakType=None,regridType=None):
        """Shortcut mimic for setPar"""
        self.setPar(gasType=gasType,cloudType=cloudType,gasFile=gasFile,cloudFile=cloudFile,tweakType=tweakType,regridType=regridType)
    def setPar(self,gasType=None,cloudType=None,gasFile=None,cloudFile=None,tweakType=None,regridType=None):
        """Sets atmosphere run parameters"""
        if gasType!=None:    ### Should be read or compute
            self.atmType = atmType
        if cloudType!=None:  ### Should be read or compute
            self.cloudType = cloudType
        if gasFile!=None:
            self.gasFile=gasFile
        if cloudFile!=None:
            self.cloudFile=cloudFile
        if tweakType!=None:  ### Should be a planet name else no tweaking
            if tweakType == 'auto':
                self.tweakType = self.planet
            else:
                self.tweakType=tweakType
        if regridType!=None:
            self.regridType=regridType

    def dp(self):
        """Shortcut mimic for dispPar but only prints, doesn't return"""
        d=self.dispPar()
        print d
    def dispPar(self):
        """Displays parameters.  Only returns doesn't print."""
        
        d = 'Run parameters for '+self.planet
        d+= '\n\tgasType:  '+self.gasType
        d+= '\n\tgasFile:  '+self.gasFile
        d+= '\n\tcloudType:  '+self.cloudType
        d+= '\n\tcloudFile:  '+self.cloudFile
        d+= '\n\totherType:  '+self.otherType
        d+= '\n\ttweakType:  '+self.tweakType
        d+= '\n\tregridType: '+self.regridType
        d+= '\n'
        return d

    def plotTP(self):
        """Plot the T-P profile"""
        plt.figure(planetDictionary[self.planet])
        plt.title(self.planet+':  T-P profile')
        plt.loglog(self.gas[self.C['T']],self.gas[self.C['P']])
        plt.axis([10,1000,1500,1E-7])
        plt.ylabel('P [bars]')
        plt.xlabel('T [K]')

    def plotCloud(self,dontPlot=['Z','P','T','DZ']):
        """Plots the clouds"""
        plt.figure(planetDictionary[self.planet]+11)
        plt.title(self.planet+': clouds')
        for cloud in self.Cl:
            present, cl = self.__isPresent__(self.cloud[self.Cl[cloud]])
            if cloud in dontPlot or not present:
                continue
            plt.loglog(cl,self.cloud[self.Cl['P']], label=cloud)
        plt.axis([1.0E-10,0.1,1500.0,1E-7])
        plt.ylabel('P [bars]')
        plt.xlabel('Density [g/cm^3]')
        plt.legend()

    def plotGas(self,dontPlot=['Z','P','T','DZ']):
        """Plots the constituents"""
        plt.figure(planetDictionary[self.planet]+10)
        plt.title(self.planet+': gas')
        for gas in self.C:
            present, g = self.__isPresent__(self.gas[self.C[gas]])
            if gas in dontPlot or not present:
                continue
            plt.loglog(g,self.gas[self.C['P']], label=gas)
        plt.axis([1.0E-10,1.0,1500.0,1E-7])
        plt.ylabel('P [bars]')
        plt.xlabel('Fractional Abundance')
        plt.legend()

    def plotOther(self,dontPlot=['Z']):
        plt.figure(planetDictionary[self.planet]+12)
        plt.title(self.planet+': other')
        for other in self.LP:
            present, g = self.__isPresent__(self.layerProperty[self.LP[other]])
            if other in dontPlot or not present:
                continue
            plt.loglog(g,self.gas[self.C['P']], label=other)
        plt.axis([1.0E-5,1.0E5,1500.0,1E-7])
        plt.ylabel('P [bars]')
        plt.xlabel('Property')
        plt.legend()

    def readGas(self,gasFile=None,verbose=False):
        """Reads gas profile file as self.gas"""

        if gasFile == None:
            gasFile=self.gasFile
        gasFile = self.path+gasFile

        print 'Reading '+self.gasFile
        self.gas = []
        print '\tUsing atmsopheric constituent:  ',
        for k in self.C:
            print k+'  ',
            self.gas.append([])

        try:
            fp = open(gasFile,"r")
        except:
            print gasFile+' was not found - returning no gas profile'
            return 0
        i=0
        print ' '
        for line in fp:
            if line[0]=='!':
                if verbose:
                    print '\tHEADER: '+line,
                continue
            elif len(line) < 3:
                continue
            data = line.split()
            if i==0:
                self.nConstituent = len(data)
                print '\tReading '+str(self.nConstituent)+' components'
            for n in range(len(self.C)): ## Initialize all of the constituents to 0.0
                self.gas[n].append(0.0)
            for n,v in enumerate(data): ## ...now read in data
                self.gas[n][i] = float(v)
            if (self.nConstituent-n)!=1:
                print 'line '+str(i+1)+' has incorrect number of components'
                print '['+line+']'
            i+=1
        self.nGas = i
        print '\tRead '+str(self.nGas)+' lines'
        fp.close()
        self.gas = np.array(self.gas)
        return self.nGas

    def readCloud(self,cloudFile=None,verbose=False):
        """Reads in cloud data if we have it..."""
        if cloudFile == None:
            cloudFile=self.cloudFile
        cloudFile = self.path+cloudFile

        print 'Reading '+cloudFile
        self.cloud = []
        print '\tUsing cloud component:  ',
        for k in self.Cl:
            print k+'  ',
            self.cloud.append([])

        try:
            fp = open(cloudFile,"r")
        except:
            print cloudFile+' was not found - returning no clouds'
        i=0
        pastHeader = False
        print ' '
        for line in fp:
            if line[0]=='!':
                if verbose:
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
                self.nParticulate = len(data)
                print '\tReading '+str(self.nParticulate)+' cloud particles'
            for n in range(len(self.Cl)): ## Initialize all of the particulates to 0.0
                self.cloud[n].append(0.0)
            for n,v in enumerate(data): ## ...now read in data
                self.cloud[n][i] = float(v)
            if (self.nParticulate-n)!=1:
                print 'line '+str(i+1)+' has incorrect number of particulates'
                print '['+line+']'
            i+=1
        self.nCloud = i
        print '\tRead '+str(self.nCloud)+' lines'
        fp.close()
        self.cloud = np.array(self.cloud)
        return self.nCloud

    def readOther(self,otherFile=None,verbose=False):
        """Reads in other property data if we have it..."""
        if otherFile == None:
            otherFile = self.otherFile
        otherFile = self.path+otherFile
        print "readOther currently doesn't read in any other file..."
        print "it will probably not be needed, since it will get computed in computeOther"

    def tweakAtmNeptune(self,gas=None,cloud=None):
        """Tweak the Neptune atmosphere data..."""
        self.tweakComment = 'You should really add a comment when you tweak the atmosphere of %s' % (self.tweakType)
        print 'TWEAK:  Make this an external module.'
        if self.tweakType:
            print '-------'
            print self.tweakComment
            print '-------'
            utils.log(self.logFile,self.tweakComment,False)
        globalGas=False
        if gas==None:
            globalGas=True
            gas=self.gas
        globalCloud = False
        if cloud==None:
            globalCloud=True
            cloud=self.cloud
            
        sol = 0.0
        zSol = 0.0
        nAtm = len(gas[self.C['P']])
        if nAtm - self.nGas != 0:
            print 'Error in number of layers - check it out ('+str(nAtm)+'/'+str(self.nGas)+')'
        for i in range(nAtm):
            ### Process CO & HCN
            if gas[self.C['P']][i] > 0.1585:
                gas[self.C['CO']][i] = 0.0
            else:
                gas[self.C['CO']][i] = 1.0E-6
            gas[self.C['CO13']][i] = 1.0E-2*gas[self.C['CO']][i]
            gas[self.C['HCN']][i] = 0.0
            gas[self.C['SOLN']][i] = 0.0
            gas[self.C['PH3']][i] = 0.0
            ### Process NH3
            if gas[self.C['T']][i] > 400.0:
                gas[self.C['NH3']][i] = 1.94E-4
            ### compute dz and integrated solution cloud and solution cloud height
            if i==0:
                gas[self.C['DZ']][i] = 0.0
            else:
                gas[self.C['DZ']][i] = abs(gas[self.C['Z']][i]-gas[self.C['Z']][i-1])*1.0E5
                sol+=gas[self.C['SOLN']][i]*gas[self.C['DZ']][i]
                if gas[self.C['SOLN']][i] > 0.0:
                    zSol+=gas[self.C['DZ']][i]
        ### if using self.gas, need to write atmosphere back to self
        if globalGas:
            self.sol = sol
            self.zSol = zSol
            self.gas=gas
        if globalCloud:
            self.cloud = cloud
        return 1

    def tweakAtmJupiter(self,gas=None):
        """Tweak the atmosphere data..."""

        print 'TWEAK:  Make this an external module.'
        self.tweakComment = 'You should really add a comment when you tweak the atmosphere!'
        print tweakComment
        print "This doesn't do anything yet"
        globalGas=False
        if gas==None:
            globalGas=True
            gas=self.gas

        sol = 0.0
        zSol = 0.0
        nGas = len(gas[self.C['P']])
        for i in range(nGas):
            sol = 0.0

        ### if using self.gas, need to write atmosphere back to self
        if globalGas:
            self.sol = sol
            self.zSol = zSol
            self.gas=gas
        return gas

    def computeOther(self,verbose=False):
        """This module computes derived atmospheric properties:
           ds = length ray travels
           mu = molecular weight [AMU]
           g = local gravity [cm/s^2]
           refr = refractivity """
        self.layerProperty = []
        for op in lyrProp:
            self.layerProperty.append([])
        zOffset = 0.0
        iOffset = 0
        psep = 1.0E6
        for i,zv in enumerate(self.gas[self.C['Z']]):     # find the nearest z value at p_ref
            P = self.gas[self.C['P']][i]
            if abs(P-self.p_ref) < psep:
                psep = abs(P-self.p_ref)
                iOffset = i
        zOffset = self.gas[self.C['Z']][iOffset]
        z_at_p_ref = self.r_ref/1.0E5
        if verbose:
            print "z,P offset:  ",zOffset,self.gas[self.C['P']][iOffset]
        
        for i,zv in enumerate(self.gas[self.C['Z']]):
            T = self.gas[self.C['T']][i]
            P = self.gas[self.C['P']][i]
            self.layerProperty[self.LP['Z']].append(z_at_p_ref + zv - zOffset)   # note that this is the "actual" z referenced to planet center
            ###set mu
            mutmp = properties.AMU_H2*self.gas[self.C['H2']][i] + properties.AMU_He*self.gas[self.C['HE']][i] + \
                    properties.AMU_H2S*self.gas[self.C['H2S']][i] + properties.AMU_NH3*self.gas[self.C['NH3']][i] + \
                    properties.AMU_H2O*self.gas[self.C['H2O']][i] + properties.AMU_CH4*self.gas[self.C['CH4']][i] + \
                    properties.AMU_PH3*self.gas[self.C['PH3']][i]
            self.layerProperty[self.LP['MU']].append(mutmp)
            ###set g
            gtmp = 2.0*math.log(P/self.p_ref)*(properties.R*T)/(self.r_ref*mutmp)
            btmp = self.g_ref+gtmp
            gtmp = 0.5*( btmp + math.sqrt(btmp*btmp - gtmp*gtmp))
            self.layerProperty[self.LP['G']].append(gtmp)
            ###set R
            refrtmp = properties.REFR_H2*self.gas[self.C['H2']][i]   + properties.REFR_He*self.gas[self.C['HE']][i] + \
                      properties.REFR_H2S*self.gas[self.C['H2S']][i] + properties.REFR_NH3*self.gas[self.C['NH3']][i] + \
                      properties.REFR_H2O(T)*self.gas[self.C['H2O']][i] + properties.REFR_CH4*self.gas[self.C['CH4']][i] + \
                      properties.REFR_PH3*self.gas[self.C['PH3']][i]
            refrtmp*=P*(293.0/T)
            self.layerProperty[self.LP['REFR']].append(refrtmp)
            refrtmp = refrtmp/1.0E6 + 1.0
            self.layerProperty[self.LP['N']].append(refrtmp)
        self.layerProperty[self.LP['DS']] = ray.compute_dS(0.0,self.layerProperty[self.LP['Z']],self.layerProperty[self.LP['N']],plot='limits')
        print 'Setting convenience arrays:  self.z, self.n, self.ds'
        self.z = self.layerProperty[self.LP['Z']]
        self.n = self.layerProperty[self.LP['N']]
        self.ds = self.layerProperty[self.LP['DS']]
            
    def computeCloud(self,verbose=False):
        """This computes cloud stuff"""
        print 'computeCloud oes nothing yet.  This probably wont do anything since gas/cloud/other will get computed in the same tcm'
    def computeGas(self,verbose=False):
        """Computes an atmosphere given stuff"""
        print 'computeGas does nothing yet.  This probably wont do anything since gas/cloud/other will get computed in the same tcm'
        print 'This will probably just call a external tcm module'

    def regrid(self,regridType=None,Pmin=None,Pmax=None,kpol='linear'):
        """This puts atm and cloud on the same grid used later for calculations
            regridType is a 2, 4 or 5 element string with the following regridding options:
                 1:  variable to regrid on ('P' or 'z')
                 2:  type ('simple' or 'spline')  OR filename containing the desired levels
                 3:  log or linear regridding ('log' or 'lin')
                 4:  value (either a float or an integer for number of levels)
                 5:  ['unit' - if present assumes that the grid is on values not number of steps]
                 e.g. P spline log 100            ==> 100 layers evenly spaced in log
                      z simple lin 1 km (default) ==> layers as needed, evenly spaced by 1 km
            Pmin/Pmax are optional - defaults are min/max in both atm and cloud.
            Don't mess with kpol."""

        regridType = 'none'
        print "!!!Hardcoding regridType within 'regrid' until it's working"

        if regridType==None:
            regridType = self.regridType
        if string.lower(regridType) == 'none' or regridType == None:
            print 'No regridding.  Note that there is a risk that not everything is on the same grid...'
            return 0
        if string.lower(regridType) == 'auto' or string.lower(regridType) == 'default':
            regridType = 'z simple lin 1 km'

        regrid = regridType.split()
        xvar = string.upper(regrid[0])
        loglin = string.lower(regrid[1])
        
        if xvar != 'Z':
            print 'only support Z right now'
            return 0

        if Pmin==None:
            m = [min(self.gas[self.C['P']]), min(self.cloud[self.Cl['P']])]
            Pmin = min(m)
        if Pmax==None:
            m = [max(self.gas[self.C['P']]), max(self.cloud[self.Cl['P']])]
            Pmax = max(m)
        m = [min(self.gas[self.C['Z']]), min(self.cloud[self.Cl['Z']])]
        zmin = min(m)
        m = [max(self.gas[self.C['Z']]), max(self.cloud[self.Cl['Z']])]
        zmax = max(m)

        ###Need to make the numpy arrays back into lists  (Do I REALLY need to?)
        self.gas = list(self.gas)
        self.cloud = list(self.cloud)

        ###Prep pressure to use in setting/checking limits
        Ptmp = []
        for pv in self.gas[self.C['P']]:
            if loglin == 'log':
                Ptmp.append(math.log(pv))
            else:
                Ptmp.append(pv)
        if loglin == 'log':
            Pmax = math.log(Pmax)
            Pmin = math.log(Pmin)
        elif loglin == 'lin':
            pass
        else:
            print loglin+' is not valid'
            return 0

        ###xarr is the interpolating abscissa
        ###note that the interpolation needs xarr increasing, so multiply by -1 for that purpose (given definition of z)
        ###    this is confusing - esp watch for zmin/zmax issues
        xarr = []
        for v in self.gas[self.C[xvar]]:
            xarr.append(-1.0*v)
        schk = interp1d(xarr,Ptmp,kind=kpol,bounds_error=True)  #want to check pressures
        ###xs is the desired abscissa (note it is 'backwards' like xarr)
        xs = []
        if len(regrid) == 4:
            zstep = float(regrid[2])
            zunit = regrid[3]   ### should check units but don't yet
            print 'interpolating on '+xvar+' at '+str(zstep)+' '+zunit+' ('+loglin+')'
            ztry = math.floor(xarr[0])
            while ztry < zmin:
                ztry+=zstep
                approxP = schk(ztry)
                if approxP < Pmin:
                    continue
                elif approxP > Pmax:
                    break
                xs.append(ztry)
        else:
            print 'nothing yet'
            return 0

        if regrid[2].lower() == 'simple':
            ###This does the 'simple' regridding
            print 'simple regridding just linearly interpolates the values onto the common grid between values'
        else:
            ###This does the spline regridding
            ###interpolate gas onto grid
            xarr = []  ###This duplicates above, but just in case that changes...
            for v in self.gas[self.C[xvar]]:
                xarr.append(-1.0*v)
            for yvar in self.C:
                present, logVer = self.__isPresent__(self.gas[self.C[yvar]])
                if yvar == xvar:
                    continue
                yarr = []
                i = 0
                for val in self.gas[self.C[yvar]]:
                    if loglin == 'log' and present:
                        yarr.append(math.log(logVer[i]))
                    else:
                        yarr.append(val)
                    i+=1
                s = interp1d(xarr,yarr,kind=kpol,bounds_error=True)
                gas = []
                if loglin == 'log' and present:
                    for g in s(xs):
                        gas.append(math.exp(g))
                else:
                    gas = s(xs)
                self.gas[self.C[yvar]] = np.array(gas)
            ###interpolate cloud onto grid:  note that ranges are likely contained within gas, so P,T,DZ are set to gas and
            ###    outside range for clouds is set to 0.0  (Z is handled separately)
            print 'DDB: (atmosphere.py:L519):  NOTE THAT CLOUD REGRIDDING DOES NOT WORK AT THE MOMENT'
            copyOver = ['P','T','DZ']
            for co in copyOver:
                self.cloud[self.Cl[co]] = list(self.gas[self.C[co]])
            xarr = []
            for v in self.cloud[self.Cl[xvar]]:
                xarr.append(-1.0*v)
            for yvar in self.Cl:
                present, logVer = self.__isPresent__(self.cloud[self.Cl[yvar]])
                if yvar == xvar or yvar in copyOver:
                    continue
                yarr = []
                i = 0
                print 'Cloud '+yvar+'  ',present
                for val in self.cloud[self.Cl[yvar]]:
                    if loglin == 'log' and present:
                        yarr.append(math.log(logVer[i]))
                    else:
                        yarr.append(val)
                    i+=1
                s = interp1d(xarr,yarr,kind=kpol,bounds_error=False,fill_value=0.0)
                cloud = []
                if loglin == 'log' and present:
                    for cl in s(xs):
                        cloud.append(math.exp(cl))
                else:
                    cloud = s(xs)
                self.cloud[self.Cl[yvar]] = list(cloud)
            ### change xs back to decreasing
            for i in range(len(xs)):
                xs[i]*=-1.
            self.gas[self.C[xvar]] = list(xs)
            self.cloud[self.Cl[xvar]] = list(xs)

        ### Make back to numpy arrays
        self.gas = np.array(self.gas)
        self.cloud = np.array(self.cloud)
        
        return 1


    def __isPresent__(self,c,tiny=1.0E-30):
        """This checks to see if a constituent is there and sets 0.0 or negative values to tiny.  This is generally for log plotting"""
        vsum = 0.0
        vnew = []
        for v in c:
            vsum+=v
            if v==0.0:
                vnew.append(tiny)
            elif v<0.0:
                vnew.append(tiny)
            else:
                vnew.append(v)
        if vsum<tiny:
            present = False
        else:
            present = True
        return present, vnew
