import string
import utils
import os.path

class planetConfig:
    def __init__(self,planet,configFile=None,path=None,log=None,verbose=False,printHelp=False):
        """reads in config file"""

        deprecated_toks = {'gastype':'gasType','cloudtype':'cloudType','othertype':'otherType','otherfile':'otherFile'}
        toks = {'gasfile':'gasFile','cloudfile':'cloudFile','constituents':'C','clouds':'Cl','tweakfile':'tweakFile','regridtype':'regridType',
                'pmin':'pmin','pmax':'pmax','omega':'omega_m','jn':'Jn','gm':'GM_ref','req':'Req','rpol':'Rpol','distance':'distance','rj':'RJ',
                'p_ref':'p_ref','zonal':'zonal','gtype':'gtype','orientation':'orientation','h2state':'h2state','doppler':'Doppler',
                'h2newset':'h2newset','water':'water_p','ice':'ice_p','nh4sh':'nh4sh_p','nh3ice':'nh3ice_p','h2sice':'h2sice_p','ch4':'ch4_p'}
        toksFloat = {'pmin':'pmin','pmax':'pmax','omega':'omega_m','jn':'Jn','gm':'GM_ref','req':'Req','rpol':'Rpol','distance':'distance','rj':'RJ',
                'p_ref':'p_ref','zonal':'zonal','gtype':'gtype','orientation':'orientation','h2state':'h2state','doppler':'Doppler',
                'h2newset':'h2newset','water':'water_p','ice':'ice_p','nh4sh':'nh4sh_p','nh3ice':'nh3ice_p','h2sice':'h2sice_p','ch4':'ch4_p'}
        toks = {'gasfile':'gasFile','cloudfile':'cloudFile','constituents':'C','clouds':'Cl','tweakfile':'tweakFile','regridtype':'regridType',
                'pmin':'pmin','pmax':'pmax','omega':'omega_m','jn':'Jn','gm':'GM_ref','req':'Req','rpol':'Rpol','distance':'distance','rj':'RJ',
                'p_ref':'p_ref','zonal':'zonal','gtype':'gtype','orientation':'orientation','h2state':'h2state','doppler':'Doppler',
                'h2newset':'h2newset','water':'water_p','ice':'ice_p','nh4sh':'nh4sh_p','nh3ice':'nh3ice_p','h2sice':'h2sice_p','ch4':'ch4_p'}
        toks = {'gasfile':'gasFile','cloudfile':'cloudFile','constituents':'C','clouds':'Cl','tweakfile':'tweakFile','regridtype':'regridType',
                'pmin':'pmin','pmax':'pmax','omega':'omega_m','jn':'Jn','gm':'GM_ref','req':'Req','rpol':'Rpol','distance':'distance','rj':'RJ',
                'p_ref':'p_ref','zonal':'zonal','gtype':'gtype','orientation':'orientation','h2state':'h2state','doppler':'Doppler',
                'h2newset':'h2newset','water':'water_p','ice':'ice_p','nh4sh':'nh4sh_p','nh3ice':'nh3ice_p','h2sice':'h2sice_p','ch4':'ch4_p'}
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
        #   note that the tok values can differ, hence the tok dictionary (although generally just all lowercase)
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
        self.h2newset = 'forget... (bool)'
        self.water_p = 'water particle size [um?] (float)'
        self.ice_p = 'ice particle size [um?] (float)'
        self.nh4sh_p = ' nh4sh particle size [um?] (float)'
        self.nh3ice_p = 'ammonia-ice particle size [um?] (float)'
        self.h2sice_p = 'h2s-ice particle size [um?] (float)'
        self.ch4_p = 'methane particle size [um?] (float)'
        self.Doppler = 'use Doppler or not (bool)'
        if printHelp:
            print 'eventually will print out tok help...'

        self.planetaryDefaults()
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
                try:
                    self.cloudFileHdr = float(data[1])
                except:
                    self.cloudFileHdr = 0
            elif tok == 'cloudfile':
                self.cloudFile = data[0]
                try:
                    self.cloudFileHdr = float(data[1])
                except:
                    self.cloudFileHdr = 0
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
                    self.omega_m = float(data[0])
                except:
                    self.omega_m = data[0]
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
            elif tok == 'h2newset':
                if data[0][0].lower() == 'f':
                    self.h2newset = False
                else:
                    self.h2newset = True
            elif tok=='water':
                self.water_p = float(data[0])
                #  should also check for units!!!
            elif tok=='ice':
                self.ice_p = float(data[0])
                #  should also check for units!!!
            elif tok=='nh4sh':
                self.nh4sh_p = float(data[0])
                #  should also check for units!!!
            elif tok=='nh3ice':
                self.nh3ice_p = float(data[0])
                #  should also check for units!!!
            elif tok=='h2sice':
                self.h2sice_p = float(data[0])
                #  should also check for units!!!
            elif tok=='ch4':
                self.ch4_p = float(data[0])
            else:
                print 'Unknown token:  ',tok
                nTokMod-=1
        fp.close()

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
        
        return nTokMod

    def show(self):
        """Displays configuration and returns string.  See __init__ and setConfig."""
        s = 'Run parameters:\n'
        variables =  vars(self)
        for key in variables:
            s+='\t%s:  %s\n' %(key,str(variables[key]))
        print 'CONFIG.PY:360 JUST TO REMIND HOW TO ACCESS VARIABLES ==>',self.__dict__['Doppler']
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

