import string
import utils
import os.path

class planetConfig:
    def __init__(self,planet,configFile='config.par',path=None,log=None,verbose=False,printHelp=False):
        """reads in config file"""
        self.toks = {'gasfile':['gasFile',str], 'cloudfile':['cloudFile',str], 'gasfilehdr':['gasFileHdr',int], 'cloudfilehdr':['cloudFileHdr',int],
                     'constituents':['C',str], 'clouds':['Cl',str], 'tweakfile':['tweakFile',str], 'regridtype':['regridType',str],
                     'pmin':['pmin',float], 'pmax':['pmax',float], 'omega':['omega_m',float], 'jn':['Jn',str], 'gm':['GM_ref',float],
                     'req':['Req',float], 'rpol':['Rpol',float], 'distance':['distance',float], 'rj':['RJ',float], 'p_ref':['p_ref',float], 'zonal':['zonal',str],
                     'gtype':['gtype',str], 'orientation':['orientation',str], 'h2state':['h2state',str], 'doppler':['Doppler',str], 'h2newset':['h2newset',str],
                     'water':['water_p',float],'ice':['ice_p',float],'nh4sh':['nh4sh_p',float],'nh3ice':['nh3ice_p',float],'h2sice':['h2sice_p',float],'ch4':['ch4_p',float],
                     'limb':['limb',str]}
        planet = string.capitalize(planet)
        self.planet = planet
        self.verbose=verbose
        self.logFile = utils.setupLogFile(log)
        if path:
            self.path = path
        else:
            self.path = planet
        configFile = os.path.join(self.path,configFile)

        print '\n---Setting config for %s---\n' % (planet)

        #Set "universal" defaults (i.e. things not currently set but still used internally)
        self.gasType = 'read'       # "read vs compute -- normally read (currently can't compute)"
        self.cloudType = 'read'     # "read vs compute -- normally read (currently can't compute)"
        self.otherType = 'compute'  # "read vs compute -- normally compute (currently can't read)"
        self.otherFile = None       # if otherType could 'read', this file would have the data...
        self.LP = {'Z':0,'R':1,'P':2,'GM':3,'AMU':4,'REFR':5,'N':6,'H':7,'LAPSE':8,'g':9}    #...and this would hold its dictionary properties

        #These are the current config values -- get seeded with help text
        #   note that the tok names largely differ from variable name, hence the tok dictionary (although generally just all lowercase)
        self.gasFile = 'atmospheric constituent file - column order set by C (string)'
        self.gasFileHdr = 'number of header lines (to ignore) in gasFile (int)'
        self.cloudFile = 'atmospheric cloud file - column order set by Cl (string)'
        self.cloudFileHdr = 'number of header lines (to ignore) in cloudFile (int)'
        self.tweakFile = 'module that tweaks the read atmosphere (string)'
        self.C  = "This 'C'onstituent dictionary has atm layer gas data.   Needs to correspond to datafile if not computing {string}"
        self.Cl = "This 'Cl'oud dictionary has the cloud parameters.  Needs to correspond to datafile if not computing {string}"
        self.regridType = 'instructions or file to regrid data (see regrid.py) (string/int)'
        self.h2state = 'hydrogen state [e or n] (char)'
        self.pmin = 'pmin that is used in the regrid (none uses file min) [bars] (float)'
        self.pmax = 'pmax that is used in the regrid (none uses file max) [bars] (float)'
        self.distance = 'distance to the planet [AU] (float)'
        self.p_ref = 'the pressure where the radii are specified (Req/Rpol) [bars] (float)'
        self.Req = 'equatorial radius at p_ref [km] (float)'
        self.Rpol = 'polar radius at p_ref [km] (float)'
        self.gtype = 'type for planet shape [ellipse/reference/gravity/sphere] (string)'
        self.orientation = 'position angle and sub-earth latitude (planetographic) [float,float]'
        self.RJ = 'radius for gravity terms (usually Req) [km] (float)'
        self.GM_ref = 'GM at radius RJ [km3/s2] (float)'
        self.Jn = 'gravity terms (float,float,float,float,float,float)'
        self.omega_m = 'rotation velocity [rad/s] (float)'
        self.zonal = 'file with zonal winds (string)'
        self.h2newset = 'related to h2_orton - can be deprecated? (bool)'
        self.water_p = 'water particle size [um?] (float)'
        self.ice_p = 'ice particle size [um?] (float)'
        self.nh4sh_p = ' nh4sh particle size [um?] (float)'
        self.nh3ice_p = 'ammonia-ice particle size [um?] (float)'
        self.h2sice_p = 'h2s-ice particle size [um?] (float)'
        self.ch4_p = 'methane particle size [um?] (float)'
        self.Doppler = 'use Doppler or not (bool)'
        self.limb = 'limb type - used in compute_ds to test limb darkening [shape/sec] (str)'
        if printHelp:
            print 'eventually will print out tok help...'

        self.planetaryDefaults()   # this sets to base defaults, which get overwritten if valid config file
        self.setConfig(configFile)
        pars = self.show()
        utils.log(self.logFile,planet,False)
        utils.log(self.logFile,configFile,False)
        utils.log(self.logFile,pars,True)


    def setConfig(self,configFile):
        """Reads in config files and updates after default set in __init__.  These are all shown in showConfig"""
        if configFile == None:
            print 'Atmosphere:  using default config'
            return 0
        try:
            fp = open(configFile,'r')
        except IOError:
            print configFile+' not found.  Using defaults.'
            return 0
        print 'Reading '+configFile

        nTokMod = 0
        for line in fp:
            if line[0] in utils.commentChars or len(line)<4:
                continue
            if '#' in line:
                line = line[:line.index('#')]
            data = line.split()
            tok = data[0].lower()
            if tok not in self.toks:
                print 'token %s not found' % (tok)
                continue
            nTokMod += 1
            
            del(data[0])
            line = ''
            for w in data:
                line+=(w+' ')

            # First pass
            try:
                self.__dict__[self.toks[tok][0]] = self.toks[tok][1](data[0])
            except:
                self.__dict__[self.toks[tok][0]] = str(data[0])
            if type(self.__dict__[self.toks[tok][0]]) == str:
                b = self.__dict__[self.toks[tok][0]].lower()
                if b == 'none':
                    self.__dict__[self.toks[tok][0]] = None
                elif b == 'true':
                    self.__dict__[self.toks[tok][0]] = True
                elif b == 'false':
                    self.__dict__[self.toks[tok][0]] = False

            # Second pass for some
            if tok == 'gasfile':
                try:
                    self.gasFileHdr = int(data[1])
                except:
                    self.gasFileHdr = 0
            elif tok == 'cloudfile':
                try:
                    self.cloudFileHdr = int(data[1])
                except:
                    self.cloudFileHdr = 0
            elif tok == 'constituents':
                self.C = {}
                for i,w in enumerate(data):
                    w = w.upper()
                    self.C[w] = i
            elif tok == 'clouds':
                self.Cl = {}
                for i,w in enumerate(data):
                    w = w.upper()
                    self.Cl[w] = i
            elif tok == 'regridtype':
                self.regridType = line
            elif tok == 'pmin':
                self.pmin = self.rescale(self.pmin,data,utils.processingPressureUnit)
            elif tok == 'pmax':
                self.pmax = self.rescale(self.pmax,data,utils.processingPressureUnit)
            elif tok == 'jn':
                try:
                    self.Jn = []
                    for d in data:
                        self.Jn.append( float(d) )
                except:
                    self.Jn = data[0]
            elif tok == 'gm':
                self.GM_ref = self.rescale(self.GM_ref,data,utils.processingAccelUnit)
            elif tok == 'req':
                self.Req = self.rescale(self.Req,data,utils.processingAtmLayerUnit)
            elif tok == 'rpol':
                self.Rpol = self.rescale(self.Rpol,data,utils.processingAtmLayerUnit)
            elif tok == 'distance':
                self.distance = self.rescale(self.distance,data,utils.processingAtmLayerUnit)
            elif tok == 'rj':
                self.RJ = self.rescale(self.RJ,data,utils.processingAtmLayerUnit)
            elif tok == 'p_ref':
                self.p_ref = self.rescale(self.p_ref,data,utils.processingPressureUnit)
            elif tok == 'orientation':
                self.orientation = []
                try:
                    self.orientation = [float(data[0]),float(data[1])]
                except:
                    self.orientation = data[0]
            elif tok=='doppler':
                if type(self.Doppler) is not bool:
                    if data[0] in utils.affirmative:
                        self.Doppler = True
                    else:
                        self.Doppler = False
            elif tok=='h2newset':
                if type(self.h2newset) is not bool:
                    if data[0] in utils.affirmative:
                        self.h2newset = True
                    else:
                        self.h2newset = False
        fp.close()

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
            self.vwdat = [0.0,0.0]

        
        return nTokMod

    def rescale(self,v,d,procUnit):
        inpUnit = procUnit
        if len(d) > 1:
            if d[1] in utils.Units:
                inpUnit = d[1]
        if type(v) == float:
            v*= utils.Units[inpUnit]/utils.Units[procUnit]
        return v       

    def show(self):
        """Displays configuration and returns string.  See __init__ and setConfig."""
        s = 'Run parameters:\n'
        #---Below gives all variables defined in class---
        #variables =  vars(self)
        #for key in variables:
        #    s+='\t%s:  %s\n' %(key,str(variables[key]))
        keys = self.toks.keys()
        keys.sort()
        for key in keys:
            s+='\t%15s:  %15s \t%s\n' % (key, str(type(self.__dict__[self.toks[key][0]])), str(self.__dict__[self.toks[key][0]]))
        return s

    def planetaryDefaults(self):
        """This sets some defaults (mainly to have them in a different place than just the directory config.par file, which gets changed)"""
        # Set defaults
        if self.planet=='Neptune':
            self.gasFile = 'neptune.paulCO_cloud21_fletcher_best_dry'
            self.gasFileHdr = 0
            self.cloudFile = 'nepcloud_CO.cloud21_fletcher_best_dry'
            self.cloudFileHdr = 0
            self.tweakFile = 'NeptuneTweak'
            self.otherFile = None
            self.C = {'Z':0,'T':1,'P':2,'H2':3,'HE':4,'CH4':5,'NH3':6,'H2O':7,'H2S':8,'SOLN':9,'OTHER':10,
                      'PH3':11,'CO':12,'CO13':13,'HCN':14,'DZ':15}
            self.Cl= {'Z':0,'T':1,'P':2,'SOLN':3,'H2O':4,'NH4SH':5,'NH3':6,'H2S':7,'CH4':8,'AR':9,'PH3':10, 'DZ':11}
            self.regridType = 'grid.dat'
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
            self.h2newset = True
            self.limb = 'shape'
        elif self.planet=='Jupiter':
            self.gasFile = 'jupiter.paulSolar'
            self.gasFileHdr = 0
            self.cloudFile = 'jupiter.paulclSolar'
            self.cloudFileHdr = 0
            self.tweakFile = None
            self.otherFile = None
            self.C = {'Z':0,'T':1,'P':2,'H2':3,'HE':4,'CH4':5,'NH3':6,'H2O':7,'H2S':8,'SOLN':9,'OTHER':10,
                      'PH3':11,'CO':12, 'CO13':13,'HCN':14,'DZ':15}
            self.Cl= {'Z':0,'T':1,'P':2,'SOLN':3,'H2O':4,'NH4SH':5,'NH3':6,'H2S':7,'CH4':8,'AR':9,'PH3':10, 'DZ':11}
            self.regridType = 'grid.dat'
            self.h2state = 'e'
            self.pmin = None
            self.pmax = None
            self.distance = 5.2*utils.Units['AU']/utils.Units[utils.processingAtmLayerUnit]
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
            self.h2newset = True
            self.limb = 'shape'

