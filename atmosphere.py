import math
import string
import matplotlib.pyplot as plt
import numpy as np
import utils
import sys
import os
import os.path

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


##########Use all upper case to help getVal##########
###This 'C'onstituent dictionary has atm layer gas data.   Needs to correspond to datafile if not computing.
CJupiter = {'Z':0,'T':1,'P':2,'H2':3,'HE':4,'CH4':5,'NH3':6,'H2O':7,'H2S':8,'SOLN':9,'OTHER':10, 'DZ':11}
CSaturn = {}
CUranus = {}
CNeptune = {'Z':0,'T':1,'P':2,'H2':3,'HE':4,'CH4':5,'NH3':6,'H2O':7,'H2S':8,'SOLN':9,'OTHER':10,
           'PH3':11,'CO':12,'CO13':13,'HCN':14,'DZ':15}

###This 'Cl'oud dictionary has the cloud parameters.  Needs to correspond to datafile if not computing.
ClJupiter = {'Z':0,'T':1,'P':2,'SOLN':3,'H2O':4,'NH4SH':5,'NH3':6,'H2S':7,'CH4':8,'AR':9,'PH3':10, 'DZ':11}
ClSaturn = {}
ClUranus = {}
ClNeptune = {'Z':0,'T':1,'P':2,'SOLN':3,'H2O':4,'NH4SH':5,'NH3':6,'H2S':7,'CH4':8,'AR':9,'PH3':10, 'DZ':11}

###This 'lyr' dictionary has other atmospheric layer properties
lyrProp = {'Z':0,'R':1,'P':2,'GM':3,'AMU':4,'REFR':5,'N':6}

class atmosphere:
    def __init__(self,planet,path=None,config=None,log=None,verbose=False,plot=True):
        """reads/computes atmospheres.  This should return:
               self.gas
               self.cloud
               self.layerProperty
            on the appropriate grid
            Note that the research is in the input files and modifying the tweak modules
            All of the default config parameters are hard-coded here:  see __init__, setConfig, showConfig."""

        planet = string.capitalize(planet)
        self.planet = planet
        self.verbose=verbose
        self.plot=plot
        self.logFile = utils.setupLogFile(log)
        print '\n---Atmosphere---\n'
        
        self.LP = lyrProp
        ###Set default input files and Constituent(C)/Particulate(Cl) dictionaries
        if path:
            self.path = os.path.join(path,planet)
        else:
            self.path = planet
        if planet=='Neptune':
            self.gasFile = 'neptune.paulCO_cloud21_fletcher_best_dry'
            self.cloudFile = 'nepcloud_CO.cloud21_fletcher_best_dry'
            self.tweakFile = 'NeptuneTweak'
            self.otherFile = None
            self.C = CNeptune
            self.Cl= ClNeptune
            self.gasType = 'read'        ### read vs compute
            self.cloudType = 'read'      ###     "
            self.otherType = 'compute'   ###     "
            self.regridType = 'z grid.dat'
            self.h2state = 'e'
            self.pmin = None
            self.pmax = None
            self.distance = 29.704*utils.Units['AU']/utils.Units[utils.processingAtmLayerUnit]
            self.p_ref = 1.0           # bars
            self.Req = 24766.0         # in km
            self.Rpol = 24342.0
            self.gtype = 'ellipse'
            self.orientation = [348.8274, -28.97]   # position angle and sub-earth latitude (planetographic)
            self.RJ = self.Req
            self.GM_ref = 0.6835096e7
            self.Jn = [0.0, 0.0, 0.3539e-2, 0.0, -0.28e-4, 0.0, 0.0]
            self.omega_m = 1.01237195e-4
            self.zonal = 'zonalNeptune.dat'
            self.Doppler = False
        elif planet=='Jupiter':
            self.gasFile = 'jupiter.paulSolar'
            self.cloudFile = 'jupiter.paulclSolar'
            self.tweakFile = None
            self.otherFile = None
            self.C = CJupiter
            self.Cl= ClJupiter
            self.gasType = 'read'       ### read vs compute
            self.cloudType = 'read'     ###     "
            self.otherType = 'compute'  ###     "
            self.regridType = 'z lin 1.0 km'
            self.h2state = 'e'
            self.pmin = None
            self.pmax = None
            self.p_ref = 1.0           # bars
            self.Req = 71492.0         # in km
            self.Rpol = 66854.0
            self.gtype = 'ellipse'
            self.orientation = [3.12, 0.0]
            self.RJ = self.Req
            self.GM_ref = 12.6686538e7
            self.Jn = [0.0, 0.0, 1.4697e-2, 0.0, -5.84e-4, 0.0, 0.31e-4]
            self.omega_m = 1.7585e-4
            self.zonal = 'Jupiter/zonalJupiter.dat'
            self.Doppler = False
        else:
            self.gasFile = None
            self.cloudFile = None
            self.tweakFile = None
            self.otherFile = None
            self.regridType = 'z lin 1.0 km'
            self.path = None
            self.C = None
            self.Cl= None
            self.pmin = None
            self.pmax = None
            print 'No planet values set for '+planet
        self.setConfig(config)
        pars = self.showConfig()
        utils.log(self.logFile,planet,False)
        utils.log(self.logFile,config,False)
        utils.log(self.logFile,pars,True)
        if self.gtype == 'geoid' or self.Doppler:
            try:
                filename = os.path.join(self.path,self.zonal)
                fp = open(filename,'r')
                self.vwlat = []
                self.vwdat = []
                for line in fp:
                    data = line.split()
                    self.vwlat.append(float(data[0]))
                    self.vwdat.append(float(data[1]))
                fp.close()
            except IOError:
                self.vwlat = [0.0,90.0]
                self.vw.dat = [0.0,0.0]
        else:
            self.vwlat = [0.0,90.0]
            self.vwdat = [0.0,0.0]

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

        print 'Planet '+self.planet
        if self.gasType == 'read':  # this assumes that cloudType is then also 'read'
            utils.log(self.logFile,'\tReading from: '+self.path,True)
            utils.log(self.logFile,'\tAtmosphere file:  '+self.gasFile,True)
            utils.log(self.logFile,'\tCloud file:  '+self.cloudFile,True)
        if verbose:
            print self.dispPar()
        # print '\tIf in interactive mode, change any parameters in setpar() before run()'

    def run(self,Pmin=None,Pmax=None,regridType=None,gasType=None,cloudType=None,otherType=None,tweak=True,plot=None,verbose=None):
        """This is the standard pipeline"""
        if Pmin == None:
            Pmin = self.pmin
        if Pmax == None:
            Pmax = self.pmax
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

        if tweak:  # This loads and calls the module 'tweakFile'
            self.tweakAtm()

        ### Compute other parameters that are needed
        if not self.propGen.has_key(self.otherType):
            print 'Error:  no such otherTpe: '+otherType
            return 0
        else:
            self.propGen[otherType](verbose=verbose)
        
        regridded = regrid.regrid(self,regridType=regridType,Pmin=Pmin,Pmax=Pmax)
        self.nAtm = len(self.gas[0])

        angularDiameter = 2.0*math.atan(self.layerProperty[self.LP['R']][0]/self.distance)
        print 'angular radius = %f arcsec' % ((180.0/np.pi)*3600.0*angularDiameter/2.0)
        
        ### Plot data
        if plot:
            self.plotTP()
            self.plotGas()
            self.plotCloud()
            self.plotProp()
        return self.nAtm

    def setConfig(self,configFile):
        """Reads in config files and updates after default set in __init__.  These are all shown in showConfig"""
        if configFile == None:
            print 'Atmosphere:  using default config'
            return 0
        nTokMod = 0
        try:
            fp = open(configFile,'r')
        except IOError:
            print configFile+' not found.  Using defaults.'
            return 0
        print 'Reading '+configFile
        for line in fp:
            validTok = True
            if line[0] in utils.commentChars or len(line)<4:
                continue
            data = line.split()
            tok = data[0].lower()
            del(data[0])
            if data[0].lower() == 'none':
                data[0] = None
            nTokMod += 1 
            if tok == 'gasfile':
                self.gasFile = data[0]
            elif tok == 'cloudfile':
                self.cloudFile = data[0]
            elif tok == 'otherfile':
                self.otherFile = data[0]
            elif tok == 'constituents':
                for i,w in enumerate(data):
                    if w[0] == '#':
                        break
                    w = w.upper()
                    self.C[w] = i
            elif tok == 'clouds':
                for i,w in enumerate(data):
                    if w[0] == '#':
                        break
                    w = w.upper()
                    self.Cl[w] = i
            elif tok == 'gastype':
                self.gasType = data[0]
            elif tok == 'cloudtype':
                self.cloudType = data[0]
            elif tok == 'othertype':
                self.otherType = data[0]
            elif tok == 'tweakfile':
                self.tweakFile = data[0]
            elif tok == 'regridtype':
                s=''
                for w in data:
                    if w[0]=='#':
                        break
                    else:
                        s+= w+' '
                self.regridType = s
            elif tok == 'pmin':
                try:
                    self.pmin = float(data[0])
                    try:
                        inpUnit = data[1]
                        if inpUnit not in utils.Units:
                            inpUnit = utils.processingPressureUnit
                    except:
                        inpUnit = utils.processingPressureUnit
                    self.pmin*= utils.Units[inpUnit]/utils.Units[utils.processingPressureUnit]
                except:
                    self.pmin = data[0]
            elif tok == 'pmax':
                try:
                    self.pmax = float(data[0])
                    try:
                        inpUnit = data[1]
                        if inpUnit not in utils.Units:
                            inpUnit = utils.processingPressureUnit
                    except:
                        inpUnit = utils.processingPressureUnit
                    self.pmax*=utils.Units[inpUnit]/utils.Units[utils.processingPressureUnit]
                except:
                    self.pmin = data[0]
            elif tok == 'omega':
                try:
                    self.omega = float(data[0])
                except:
                    self.omega = data[0]
            elif tok == 'jn':
                try:
                    self.Jn = []
                    for d in data:
                        self.Jn.append( float(d) )
                except:
                    self.Jn = data[0]
            elif tok == 'gm':  # #########################################################
                try:
                    self.GM_ref = float(data[0])
                    try:
                        inpUnit = data[1]
                        if inpUnit not in utils.Units:
                            print inpUnit+' not found - assuming '+utils.processingAccelUnit
                            inpUnit = utils.processingAccelUnit
                    except:
                        inpUnit = utils.processingAccelUnit
                    self.GM_ref*= utils.Units[inpUnit]/utils.Units[utils.processingAccelUnit]
                except:
                    self.GM_ref = data[0]
            elif tok == 'req': #########################################################
                try:
                    self.Req = float(data[0])
                    try:
                        inpUnit = data[1]
                        if inpUnit not in utils.Units:
                            print inpUnit+' no found - assuming '+utils.processingAtmLayerUnit
                    except:
                        inpUnit = utils.processingAtmLayerUnit
                    self.Req*= utils.Units[inpUnit]/utils.Units[utils.processingAtmLayerUnit]
                except:
                    self.Req = data[0]
            elif tok == 'rpol': #########################################################
                try:
                    self.Rpol = float(data[0])
                    try:
                        inpUnit = data[1]
                        if inpUnit not in utils.Units:
                            print inpUnit+' no found - assuming '+utils.processingAtmLayerUnit
                    except:
                        inpUnit = utils.processingAtmLayerUnit
                    self.Rpol*= utils.Units[inpUnit]/utils.Units[utils.processingAtmLayerUnit]
                except:
                    self.Rpol = data[0]
            elif tok == 'distance': #########################################################
                try:
                    self.distance = float(data[0])
                    try:
                        inpUnit = data[1]
                        if inpUnit not in utils.Units:
                            print inpUnit+' no found - assuming '+utils.processingAtmLayerUnit
                    except:
                        inpUnit = utils.processingAtmLayerUnit
                    self.distance*= utils.Units[inpUnit]/utils.Units[utils.processingAtmLayerUnit]
                except:
                    self.distance = data[0]
            elif tok == 'rj': #########################################################
                try:
                    self.RJ = float(data[0])
                    try:
                        inpUnit = data[1]
                        if inpUnit not in utils.Units:
                            print inpUnit+' no found - assuming '+utils.processingAtmLayerUnit
                    except:
                        inpUnit = utils.processingAtmLayerUnit
                    self.RJ*= utils.Units[inpUnit]/utils.Units[utils.processingAtmLayerUnit]
                except:
                    self.RJ = data[0]
            elif tok == 'p_ref': #########################################################
                try:
                    self.p_ref = float(data[0])
                    try:
                        inpUnit = data[1]
                        if inpUnit not in utils.Units:
                            print inpUnit+' no found - assuming '+utils.processingPressureUnit
                    except:
                        inpUnit = utils.processingPressureUnit
                    self.p_ref*= utils.Units[inpUnit]/utils.Units[utils.processingPressureUnit]
                except:
                    self.p_ref = data[0]
            elif tok == 'zonal':
                self.zonal = data[0]
            elif tok == 'gtype':
                self.gtype = data[0]
            elif tok == 'orientation':
                self.orientation = []
                try:
                    self.orientation = [float(data[0]),float(data[1])]
                except:
                    self.orientation = data[0]
            elif tok == 'h2state':
                self.h2state = data[0]
            elif tok=='doppler':
                self.Doppler = False
                if data[0] in utils.affirmative:
                    self.Doppler = True
                print 'Doppler setting:  ',self.Doppler
            else:
                print 'Unknown token:  ',tok
                nTokMod-=1
        fp.close()
        return nTokMod

    def showConfig(self):
        """Displays configuration and returns string.  See __init__ and setConfig."""
        s = 'Run parameters:\n'
        s+= '\tpath:  '+str(self.path)+'\n'
        s+= '\tgasFile:  '+str(self.gasFile)+'\n'
        s+= '\tcloudFile:  '+str(self.cloudFile)+'\n'
        s+= '\ttweakFile:  '+str(self.tweakFile)+'\n'
        s+= '\totherFile:  '+str(self.otherFile)+'\n'
        s+= '\tC:  '+str(self.C)+'\n'
        s+= '\tCl:  '+str(self.Cl)+'\n'
        s+= '\tgasType:  '+str(self.gasType)+'\n'
        s+= '\tcloudType:  '+str(self.cloudType)+'\n'
        s+= '\totherType:  '+str(self.otherType)+'\n'
        s+= '\tregridType:  '+str(self.regridType)+'\n'
        s+= '\th2state:  '+str(self.h2state)+'\n'
        s+= '\tpmin:  '+str(self.pmin)+'\n'
        s+= '\tpmax:  '+str(self.pmax)+'\n'
        s+= '\tdistance:  '+str(self.distance)+'\n'
        s+= '\tp_ref:  '+str(self.p_ref)+'\n'
        s+= '\tReq:  '+str(self.Req)+'\n'
        s+= '\tRpol:  '+str(self.Rpol)+'\n'
        s+= '\tgtype:  '+str(self.gtype)+'\n'
        s+= '\torientation:  '+str(self.orientation)+'\n'
        s+= '\tRJ:  '+str(self.RJ)+'\n'
        s+= '\tGM_ref:  '+str(self.GM_ref)+'\n'
        s+= '\tJn:  '+str(self.Jn)+'\n'
        s+= '\tomega_m:  '+str(self.omega_m)+'\n'
        s+= '\tzonal:  '+str(self.zonal)+'\n'
        return s

    def getval(self,val=None,vtype='all'):
        """Returns one of the constituent or cloud profiles"""
        if val == None:
            print "Usage:  getVal('v',['gas'/'cloud'/'other'/'all'])"
            print "    'gas' is default"
            print 'These are the gas values:'
            print self.C
            print 'These are the cloud values:'
            print self.Cl
            print 'These are the layerProperty (other) values:'
            print self.LP
            return None
        v = string.upper(val)
        vt = string.lower(vtype)
        rv = 0
        found = False
        
        if vt=='gas' or vt=='all':
            if self.C.has_key(v):
                rv = self.gas[self.C[v]]
                print 'Found '+val+' in gas'
                found = True
        elif vt=='cloud' or vt=='all':
            if self.Cl.has_key(v):
                rv = self.cloud[self.Cl[v]]
                print 'Found '+val+' in cloud'
                found = True
        elif vt=='other' or vt=='all':
            if self.layerProperty.has_key(v):
                rv = self.layerProperty[self.LP[v]]
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
        plt.loglog(self.gas[self.C['T']],self.gas[self.C['P']])
        v = list(plt.axis())
        v[2] = 100.0*math.ceil(self.gas[self.C['P']][-1]/100.0)
        v[3] = 1.0E-7*math.ceil(self.gas[self.C['P']][0]/1E-7)
        plt.axis(v)
        plt.ylabel('P [bars]')
        plt.xlabel('T [K]')

    def plotCloud(self,dontPlot=['Z','P','T','DZ'],plot='auto'):
        """Plots the clouds"""
        if plot=='auto':
            plt.figure(planetDictionary[self.planet]+11)
        plt.title(self.planet+': clouds')
        for cloud in self.Cl:
            present, cl = self.__isPresent__(self.cloud[self.Cl[cloud]])
            if cloud in dontPlot or not present:
                continue
            plt.loglog(cl,self.cloud[self.Cl['P']], label=cloud)
        v = list(plt.axis())
        if v[0] < 1E-10:
            v[0] = 1E-10
        v[2] = 100.0*math.ceil(self.gas[self.C['P']][-1]/100.0)
        v[3] = 1.0E-7*math.ceil(self.gas[self.C['P']][0]/1E-7)
        plt.axis(v)
        plt.ylabel('P [bars]')
        plt.xlabel(r'Density [g/cm$^3$]')
        plt.legend()

    def plotGas(self,dontPlot=['Z','P','T','DZ'],plot='auto'):
        """Plots the constituents"""
        if plot=='auto':
            plt.figure(planetDictionary[self.planet]+10)
        plt.title(self.planet+': gas')
        for gas in self.C:
            present, g = self.__isPresent__(self.gas[self.C[gas]])
            if gas in dontPlot or not present:
                continue
            plt.loglog(g,self.gas[self.C['P']], label=gas)
        v = list(plt.axis())
        if v[0] < 1E-10:
            v[0] = 1E-10
        v[2] = 100.0*math.ceil(self.gas[self.C['P']][-1]/100.0)
        v[3] = 1.0E-7*math.ceil(self.gas[self.C['P']][0]/1E-7)
        plt.axis(v)
        plt.ylabel('P [bars]')
        plt.xlabel('Fractional Abundance')
        plt.legend()

    def plotProp(self,dontPlot=['Z','P','T'],plot='auto'):
        if plot=='auto':
            plt.figure(planetDictionary[self.planet]+12)
        plt.title(self.planet+': other')
        for other in self.LP:
            present, g = self.__isPresent__(self.layerProperty[self.LP[other]])
            if other in dontPlot or not present:
                continue
            plt.loglog(g,self.gas[self.C['P']], label=other)
        v = list(plt.axis())
        if v[0] < 1E-10:
            v[0] = 1E-10
        v[2] = 100.0*math.ceil(self.gas[self.C['P']][-1]/100.0)
        v[3] = 1.0E-7*math.ceil(self.gas[self.C['P']][0]/1E-7)
        plt.axis(v)
        plt.ylabel('P [bars]')
        plt.xlabel('Property')
        plt.legend()

    def readGas(self,gasFile=None,verbose=False,numHeaderLines=1):
        """Reads gas profile file as self.gas"""

        if gasFile == None:
            gasFile=self.gasFile
        gasFile = os.path.join(self.path,gasFile)

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
        lineno = 0
        pastHeader = False
        print ' '
        for line in fp:
            lineno+=1
            data = line.split()
            if line[0]=='!' or len(data) < 4 or lineno<=numHeaderLines:
                if verbose:
                    print '\tHEADER: '+line,
                continue
            pastHeader = True
            if len(line) < 3:
                if pastHeader:
                    break
                else:
                    continue
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
        ### Check that P is monotonically increasing
        monotonic = np.all(np.diff(self.gas[self.C['P']])>0.0)
        if not monotonic:
            self.gas = np.fliplr(np.gas)
            monotonic = np.all(np.diff(self.gas[self.C['P']])>0.0)
        if not monotonic:
            print "Error in "+gasFile+".  Pressure not monotonically increasing"
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
                

    def readCloud(self,cloudFile=None,verbose=False,numHeaderLines=7):
        """Reads in cloud data if we have it..."""
        if cloudFile == None:
            cloudFile=self.cloudFile
        cloudFile = os.path.join(self.path,cloudFile)

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
        lineno = 0
        pastHeader = False
        print ' '
        for line in fp:
            lineno+=1
            if line[0]=='!' or lineno<=numHeaderLines:
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
        ### Check that P is monotonically increasing
        monotonic = np.all(np.diff(self.cloud[self.Cl['P']])>0.0)
        if not monotonic:
            self.cloud = np.fliplr(self.cloud)
            monotonic = np.all(np.diff(self.cloud[self.Cl['P']])>0.0)
        if not monotonic:
            print "Error in "+cloudFile+".  Pressure not monotonically increasing"
        return self.nCloud

    def readProp(self,otherFile=None,verbose=False):
        """Reads in other property data if we have it..."""
        if otherFile == None:
            otherFile = self.otherFile
        otherFile = os.path.join(self.path,otherFile)
        print "readProp currently doesn't read in any other file..."
        print "it will probably not be needed, since it will get computed in computeProp"
        return False

    def tweakAtm(self):
        """Tweaks the atmosphere data..."""
        nAtm = len(self.gas[self.C['P']])
        if nAtm - self.nGas != 0:
            print 'Error in number of layers - check it out ('+str(nAtm)+'/'+str(self.nGas)+')'
            print 'Returned from tweakAtm'
            return 0
        # Import tweakFile
        sys.path.append(self.path)
        try:
            __import__(self.tweakFile)
            tweakModule = sys.modules[self.tweakFile]
        except SyntaxError:
            utils.log(self.logFile,"Syntax Error:  check "+self.tweakFile,True)
            return 0
        except:
            utils.log(self.logFile,"non-syntax tweakAtm error",True)
            return 0

        # Run module then log
        self.tweakComment, self.gas, self.cloud = tweakModule.modify(self.gas,self.cloud,self.C,self.Cl)
        print '---tweakComment'
        print self.tweakComment
        print '---'
        utils.log(self.logFile,self.tweakComment,False)
        _tf = os.path.join(self.path,self.tweakFile+'.py')
        _tp = open(_tf,'r')
        dt = _tp.read()
        utils.log(self.logFile,'======================'+_tf+'=====================',False)
        utils.log(self.logFile,dt,False)
        utils.log(self.logFile,'====================================================================',False)
        _tp.close()

        return nAtm

    def computeProp(self,verbose=False):
        """This module computes derived atmospheric properties:
           amu = molecular weight [AMU]
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
        z_at_p_ref = self.Req
        if verbose:
            print "z,P offset:  ",zOffset,self.gas[self.C['P']][iOffset]
        
        for i,zv in enumerate(self.gas[self.C['Z']]):
            T = self.gas[self.C['T']][i]
            P = self.gas[self.C['P']][i]
            self.layerProperty[self.LP['P']].append(P)
            self.layerProperty[self.LP['Z']].append(zv)
            rr = z_at_p_ref + zv - zOffset
            self.layerProperty[self.LP['R']].append(rr)   # note that this is the "actual" z along equator  referenced to planet center (aka radius)
            ###set amu
            amutmp = properties.AMU_H2*self.gas[self.C['H2']][i] + properties.AMU_He*self.gas[self.C['HE']][i] + \
                    properties.AMU_H2S*self.gas[self.C['H2S']][i] + properties.AMU_NH3*self.gas[self.C['NH3']][i] + \
                    properties.AMU_H2O*self.gas[self.C['H2O']][i] + properties.AMU_CH4*self.gas[self.C['CH4']][i] + \
                    properties.AMU_PH3*self.gas[self.C['PH3']][i]
            self.layerProperty[self.LP['AMU']].append(amutmp)
            ###set GM
            if not i:
                self.layerProperty[self.LP['GM']].append(0.0)
            else:
                Gdm = self.layerProperty[self.LP['GM']][i-1] + (1.0e11)*properties.GravConst*(4.0*math.pi*amutmp*P*abs(zv-self.gas[self.C['Z']][i-1])*rr**2) / (properties.R*T)  # in km3/s2
                self.layerProperty[self.LP['GM']].append( Gdm )
            ###set R
            refrtmp = properties.REFR_H2*self.gas[self.C['H2']][i]   + properties.REFR_He*self.gas[self.C['HE']][i] + \
                      properties.REFR_H2S*self.gas[self.C['H2S']][i] + properties.REFR_NH3*self.gas[self.C['NH3']][i] + \
                      properties.REFR_H2O(T)*self.gas[self.C['H2O']][i] + properties.REFR_CH4*self.gas[self.C['CH4']][i] + \
                      properties.REFR_PH3*self.gas[self.C['PH3']][i]
            refrtmp*=P*(293.0/T)
            self.layerProperty[self.LP['REFR']].append(refrtmp)
            refrtmp = refrtmp/1.0E6 + 1.0
            self.layerProperty[self.LP['N']].append(refrtmp)
        GMnorm = self.layerProperty[self.LP['GM']][iOffset]
        for i, mv in enumerate(self.layerProperty[self.LP['GM']]):
            gm = self.GM_ref - (mv-GMnorm)
            self.layerProperty[self.LP['GM']][i] = gm
        self.layerProperty = np.array(self.layerProperty)
            
    def computeCloud(self,verbose=False):
        """This computes cloud stuff"""
        print 'computeCloud oes nothing yet.  This probably wont do anything since gas/cloud/other will get computed in the same tcm'
    def computeGas(self,verbose=False):
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
