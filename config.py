import string
import utils

class planetConfig:
    def __init__(self,planet,configFile=None,path=None,log=None,verbose=False):
        """reads in config file"""
        planet = string.capitalize(planet)
        self.planet = planet
        self.verbose=verbose
        self.logFile = utils.setupLogFile(log)
        print '\n---Reading config for %s---\n' % (planet)

        """These are the config parameters (see self.configTokLists)"""
        if path:
            self.path = os.path.join(path,planet)
        else:
            self.path = planet
        self.gasFile = 'atmospheric constituent file - column order set by C'
        self.cloudFile = 'atmospheric cloud file - column order set by Cl'
        self.tweakFile = 'module that tweaks the read atmosphere'
        self.otherFile = 'this would have the other properties, but normally is None'
        self.C  = {}  #This 'C'onstituent dictionary has atm layer gas data.   Needs to correspond to datafile if not computing.
        self.Cl = {}  #This 'Cl'oud dictionary has the cloud parameters.  Needs to correspond to datafile if not computing.
        self.LP = {}  #This 'lyr' dictionary has other atmospheric layer properties
        self.gasType = "read vs compute -- normally read (currently can't compute)"
        self.cloudType = "read vs compute -- normally read (currently can't compute)"
        self.otherType = "read vs compute -- normally compute (currently can't read"
        self.regridType = 'instructions or file to regrid data'
        self.h2state = 'hydrogen state e or n'
        self.pmin = None
        self.pmax = None
        self.distance = 1.0*utils.Units['AU']/utils.Units[utils.processingAtmLayerUnit] # distance to planet
        self.p_ref = 1.0     # bars
        self.Req = 1.0       # in km
        self.Rpol = 1.0
        self.gtype = 'ellipse'
        self.orientation = [0.0, 0.0]   # position angle and sub-earth latitude (planetographic)
        self.RJ = self.Req
        self.GM_ref = 1.0e7
        self.Jn = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        self.omega_m = 0.0
        self.zonal = 'file with zonal winds'
        self.h2newset = True
        self.water_p = 0.0
        self.ice_p = 0.0
        self.nh4sh_p = 0.0
        self.nh3ice_p = 0.0
        self.h2sice_p = 0.0
        self.ch4_p = 0.0
        self.Doppler = False

        self.LP = {'Z':0,'R':1,'P':2,'GM':3,'AMU':4,'REFR':5,'N':6,'H':7,'LAPSE':8,'g':9}
        # Set defaults
        if planet=='Neptune':
            self.gasFile = 'neptune.paulCO_cloud21_fletcher_best_dry'
            self.cloudFile = 'nepcloud_CO.cloud21_fletcher_best_dry'
            self.tweakFile = 'NeptuneTweak'
            self.otherFile = None
            self.C = {'Z':0,'T':1,'P':2,'H2':3,'HE':4,'CH4':5,'NH3':6,'H2O':7,'H2S':8,'SOLN':9,'OTHER':10,
                      'PH3':11,'CO':12,'CO13':13,'HCN':14,'DZ':15}
            self.Cl= {'Z':0,'T':1,'P':2,'SOLN':3,'H2O':4,'NH4SH':5,'NH3':6,'H2S':7,'CH4':8,'AR':9,'PH3':10, 'DZ':11}
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
            self.C = {'Z':0,'T':1,'P':2,'H2':3,'HE':4,'CH4':5,'NH3':6,'H2O':7,'H2S':8,'SOLN':9,'OTHER':10,
                      'PH3':11,'CO':12, 'CO13':13,'HCN':14,'DZ':15}
            self.Cl= {'Z':0,'T':1,'P':2,'SOLN':3,'H2O':4,'NH4SH':5,'NH3':6,'H2S':7,'CH4':8,'AR':9,'PH3':10, 'DZ':11}
            self.gasType = 'read'       ### read vs compute
            self.cloudType = 'read'     ###     "
            self.otherType = 'compute'  ###     "
            self.regridType = 'z lin 1.0 km'
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

        self.setConfig(configFile)
        pars = self.showConfig()
        utils.log(self.logFile,planet,False)
        utils.log(self.logFile,configFile,False)
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

    def configTokLists(self):
        self.ctl = ['path','gasfile','cloudfile','otherfile','constituents','clouds','gastype','cloudtype','othertype','tweakfile','regridtype','pmin','pmax',
                    'omega','jn','gm','req','rpol','distance','rj','p_ref','zonal','gtype','orientation','h2state','doppler']
        print 'These get read in at atmosphere.py:'
        print self.ctl
        self.atl = ['h2state','h2newset','water','ice','nh4sh','nh3ice','h2sice','doppler']
        print 'These get read in at alpha.py'
        print self.atl

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
        s+= '\tL:  '+str(self.LP)+'\n'
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
        s+= '\th2state:  '+str(self.h2state)+'\n'
        s+= '\tdoppler:  '+str(self.Doppler)+'\n'
        s+= '\th2newset:  '+str(self.h2newset)+'\n'
        s+= '\twater_p:  '+str(self.water_p)+'\n'
        s+= '\tice_p:  '+str(self.ice_p)+'\n'
        s+= '\tnh4sh_p:  '+str(self.nh4sh_p)+'\n'
        s+= '\tnh3ice_p:  '+str(self.nh3ice_p)+'\n'
        s+= '\th2sice_p:  '+str(self.h2sice_p)+'\n'
        s+= '\tch4_p:  '+str(self.ch4_p)+'\n'
        return s

