      ###This is the 'executive' class for planets
import math
import string
import matplotlib.pyplot as plt
import numpy as np
import atmosphere as atm
import config as  pcfg
import alpha
import brightness as bright
import data_handling
import utils
import datetime
import os
import copy

version = '1.0'

class Planet:
    def __init__(self, name, freqs=None, b=None, freqUnit='GHz', config='config.par', log='auto', verbosity=0, plot=True):
        """This is the 'executive function class to compute overall planetary emission
           Arguments here set defaults, however often get set specifically in run. See pyPlanet.pdf for documentation.
           Inputs:
                name:  'Jupiter', 'Saturn', 'Uranus', 'Neptune' [or 'functions' if you just want to load without running]
                freqs: options are:
                    - int/float:  does that one frequency
                    - list of length 3:  assumes it is [start,stop,step]
                    - list not of length 3:   does those frequencies
                b: 'impact parameter' b=1 is the radius of the maximum projected disc.
                   Determines outType from 'spectrum','profile','image' (along with freqs to some extent)
                    - doublet list is one position, [0,0] is the center
                    - float will generate a grid at that spacing, may need to set blocks during run
                    - list of length > 2, assumes a line of those locations at angle of first entry (deg)
                         if the length is four it assumes [angle,start,stop,step]
                    - 'disc' for disc-averaged
                    - 'stamp' for postage stamp (queries values)
                    - list of doublet lists, evaluate at those locations
               freqUnit: unit for above
               config:  config file name [config.par], 'manual' [equivalent none]
               log:  log data from run, either a file name, a 'no', or 'auto' (for auto filename)
               verbosity:  True/False
               plot:  True/False"""

        if name.lower()[0:4] == 'func':
            return

        planetList = ['Jupiter','Saturn','Neptune','Uranus']
        self.planet = string.capitalize(name)

        runStart = datetime.datetime.now()
        self.header = {}

        print 'Planetary modeling  (ver '+version+')\n'
        print "PLANET.PY_L51:  In alpha, clouds_idp need otherPar['refr'] - still?"

        if self.planet in planetList:
            ### Set up log file
            if string.lower(log)=='auto':
                self.logFile = '%s_%d%02d%02d_%02d%02d.log' % (self.planet,runStart.year,runStart.month,runStart.day,runStart.hour,runStart.minute)
            elif string.lower(log)=='no':
                self.logFile=None
            else:
                self.logFile = log
            self.log=utils.setupLogFile(self.logFile,path='Logs/')
            utils.log(self.log,self.planet+' start '+str(runStart),True)
            self.plot = plot
            self.verbosity = verbosity

            ### Some convenience values for the specific Neptune observations
            self.fvla_old = [4.86,8.46,14.94,22.46,43.34]
            self.fvla_new = [1.5,3.0,6.0,10.,15.,22.,33.,45.]
            self.fvla = [3.0, 6.0, 10.0, 15.0, 33.0]
            anglecap = 13.24
            bvalcap = [0.5,0.6,0.7,0.8,0.9,0.925,0.95]
            self.bvla = []
            for bval in bvalcap:
                self.bvla.append([-bval*math.sin(math.pi*anglecap/180.0),-bval*math.cos(math.pi*anglecap/180.0)])
            
            ### Get frequencies
            if freqs != None:
                freqs = self.__freqRequest__(freqs, freqUnit)
            else:
                self.freqUnit = freqUnit
            
            ### Get viewing
            self.imRow = False
            if b!= None:
                b = self.__bRequest__(b,[1,1])

            ### Get config
            if config == 'manual' or config=='none':
                config = None
            self.config = pcfg.planetConfig(self.planet,configFile=config,log=self.log,verbosity=verbosity)

            ### Create atmosphere:  attributes are self.atm.gas, self.atm.cloud and self.atm.layerProperty
            self.atm = atm.Atmosphere(self.planet,config=self.config,log=self.log,verbosity=verbosity,plot=plot)
            self.atm.run()
            self.log.flush()

            ### Read in absorption modules:  to change absorption, edit files under /constituents'
            self.alpha = alpha.Alpha(config=self.config,log=self.log,verbosity=verbosity,plot=plot)
            #self.alpha.test(f=0,verbosity=True,plot=True)
            self.log.flush()

            ### Next compute radiometric properties - initialize bright and return data class
            self.bright = bright.Brightness(log=self.log,verbosity=verbosity,plot=plot)
            self.data_return = data_handling.DataReturn()

    def run(self, freqs=[1.0,10.0,1.0],b=[0.0,0.0], freqUnit='GHz', orientation=None, block=[1,1], outputType='frequency',
            tweak_function=None):
        """Runs the model to produce the brightness temperature, weighting functions etc etc
           b = 0.04175 is a good value for Neptune images (don't remember why at the moment...)
           outputType sets whether frequency, wavelength, or both are written
           if outputType is 'batch', it forces the outType to be 'Spectrum' and outputType to be 'frequency'"""

        if not tweak_function is None:
            self.atm.run(tweak_fn=tweak_function)
            self.log.flush()


        self.imSize = None
        self.outType = None
        self.bType = None

        ###Set freqs
        reuse = False
        if freqs == None and freqUnit == None:
            freqs = self.freqs
            freqUnit = self.freqUnit
        elif freqs == 'reuse':
            freqs = self.freqs
            reuse = True
        else:
            if freqUnit == None:
                freqUnit = 'GHz'
                self.freqUnit = freqUnit
            if freqs == None:
                freqs = self.freqs
            else:
                freqs = self.__freqRequest__(freqs, freqUnit)
            freqUnit = self.freqUnit
        self.freqs = freqs
        self.freqUnit = freqUnit
        self.data_return.f = self.freqs
        
        ###Set b
        if b == None:
            b = self.b
        else:
            b = self.__bRequest__(b,block)
        if self.outType is None or self.outType == 'Spectrum' or self.outType == 'Profile':
            if len(b) > 5*len(freqs):
                self.outType = 'Profile'
            else:
                self.outType = 'Spectrum'
        if outputType == 'batch':
            self.outType = 'Spectrum'
            outputType = 'frequency'
        self.b = b
        self.data_return.b = self.b
        if self.verbosity>2:
            print 'outType = '+self.outType
        if self.outType == 'Image' and len(freqs) > 1:
            print 'Warning:  Image must be at only one frequency'
            print 'Using %f %s' % (freqs[0],self.freqUnit)
            self.freqs = list(freqs[0])
            freqs = self.freqs

        ###Start
        runStart = datetime.datetime.now()
        self.Tb=[]
        hit_b=[]
        btmp = ''
        self.rNorm = None; self.tip = None; self.rotate = None
        if self.outType == 'Image':  ##We now treat it as an image at one frequency
            if self.verbosity > 1:
                print 'imgSize = %d x %d' % (self.imSize[0], self.imSize[1])
            self.Tb_img = []
            imtmp = []
            if abs(block[1])>1:
                btmp = '_%02dof%02d'%(block[0],abs(block[1]))
            else:
                btmp = ''

        if not reuse:
            if not self.config.Doppler:
                self.bright.layerAbsorption(freqs,self.atm,self.alpha)
            
        for i,bv in enumerate(b):
            if self.verbosity>1:
                print '%d of %d (view [%.4f, %.4f])  ' % (i+1,len(b),bv[0],bv[1]),
            Tbt = self.bright.single(freqs,self.atm,bv,self.alpha,orientation,plot=not(self.outType=='Image'),discAverage=(self.bType=='disc'))
            if self.bright.path != None and self.rNorm == None:
                self.rNorm = self.bright.path.rNorm
            if self.bright.path != None and self.tip == None:
                self.tip = self.bright.path.tip
            if self.bright.path != None and self.rotate == None:
                self.rotate = self.bright.path.rotate
            if Tbt == None:  #I've now done away with returning None by returning T_cmb in brightness.py (at least I thought so...)
                Tbt = []
                for i in range(len(freqs)):
                    Tbt.append(utils.T_cmb)
            else:            #   ... so should always go to 'else'
                hit_b.append(bv)
            self.Tb.append(Tbt)
            self.data_return.Tb.append(Tbt)
            if self.outType == 'Image':
                imtmp.append(Tbt[0])
                if not (i+1)%self.imSize[0]:
                    self.Tb_img.append(imtmp)
                    imtmp = []
        self.log.flush()
        self.data_return.Tb = self.Tb

        ###Write output files (this needs to be compatible with TBfile  -- eventually should incorporate it in there###
        datFile = 'Output/%s_%s%s_%d%02d%02d_%02d%02d.dat' % (self.planet,self.outType,btmp,runStart.year,runStart.month,runStart.day,runStart.hour,runStart.minute)
        if self.verbosity>2:
            print '\nWriting '+self.outType+' data to ',datFile
        df = open(datFile,'w')
        self.__setHeader__(self.rNorm)
        self.__writeHeader__(df)
        if self.outType == 'Image':
            for data0 in self.Tb_img:
                s = ''
                for data1 in data0:
                    s+= '%7.2f\t' % (data1)
                s+='\n'
                df.write(s)
        elif self.outType == 'Spectrum':
            fp_lineoutput = open('specoutputline.dat','w')
            if outputType == 'frequency':
                s = '# %s  K@b  \t' % (self.freqUnit)
            elif outputType == 'wavelength':
                s = '# cm   K@b  \t'
            else:
                s = '# %s    cm    K@b  \t' % (self.freqUnit)
            for i,bv in enumerate(hit_b):
                s+='(%5.3f,%5.3f)\t' % (bv[0],bv[1])
            s=s.strip('\t')
            s+='\n'
            df.write(s)
            for i,f in enumerate(freqs):
                wlcm = 100.0*utils.speedOfLight / (f*utils.Units[self.freqUnit])
                if outputType == 'frequency':
                    s = '%.2f\t  ' % (f)
                elif outputType == 'wavelength':
                    s = '%.4f\t  ' % (wlcm)
                else:
                    s = '%.2f     %.4f \t ' % (f,wlcm)
                for j in range(len(hit_b)):
                    s+='  %7.2f  \t' % (self.Tb[j][i])
                s=s.strip()
                s+='\n'
                fp_lineoutput.write(s)
                df.write(s)
            fp_lineoutput.close()
        elif self.outType == 'Profile':
            if outputType == 'frequency':
                s = '# b  K@%s \t'  % (self.freqUnit)
            elif outputType == 'wavelength':
                s = '# b  K@cm  \t'
            else:
                s = '# b  K@GHz,cm  \t'
            for i,fv in enumerate(freqs):
                wlcm = 100.0*utils.speedOfLight / (fv*utils.Units[self.freqUnit])
                if outputType == 'frequency':
                    s+='  %9.4f   \t' % (fv)
                elif outputType == 'wavelength':
                    s+='  %.4f   \t' % (wlcm)
                else:
                    s+=' %.2f,%.4f\t' % (fv,wlcm)
            s.strip('\t')
            s+='\n'
            df.write(s)
            bs = []
            for i,bv in enumerate(hit_b):
                s = '%5.3f %5.3f\t' % (bv[0],bv[1])
                bs.append(math.sqrt(bv[0]**2 + bv[1]**2))
                for j in range(len(freqs)):
                    s+=' %7.2f\t ' % (self.Tb[i][j])
                s+='\n'
                df.write(s)
            if plot:
                plt.figure("Profile")
                Tbtr = np.transpose(self.Tb)
                for j in range(len(freqs)):
                    frqs = ('%.2f %s' % (self.freqs[j],self.freqUnit))
                    plt.plot(bs,Tbtr[j],label=frqs)
                plt.legend()
                plt.xlabel('b')
                plt.ylabel('$T_B$ [K]')
        else:
            print 'Invalid outType:  '+self.outType
        df.close()
        return self.data_return

    """
    BEGIN Ramsey's Edit (the main one)
    """

    def multi_run(self, parameter_list, default_P_cutoff=2,
                  freqs=[1.0, 10.0, 1.0], b=[0.0, 0.0], freqUnit='GHz', orientation=None, block=[1, 1],
                  outputType='frequency'):
        """
        Runs over a grid of pressure profiles in a 'brute force' attempt at fitting to a model.
        Does this by adjusting the Tweak function between every run.

        Make sure parameter_list is a list of lists!

        Runs the model to produce the brightness temperature, weighting functions etc etc
        b = 0.04175 is a good value for Neptune images (don't remember why at the moment...)
        outputType sets whether frequency, wavelength, or both are written
        if outputType is 'batch', it forces the outType to be 'Spectrum' and outputType to be 'frequency'
        """
        runStart = datetime.datetime.now()
        folder_path = 'Output/%s_%s_%d%02d%02d_%02d%02d' % (
            self.planet, outputType, runStart.year,
            runStart.month, runStart.day, runStart.hour,
            runStart.minute
        )
        os.mkdir(folder_path)

        original_b = copy.deepcopy(b)

        print "BEGINNING ITERATIVE RUN!"
        runtime = len(parameter_list)

        for num_run, parameter in enumerate(parameter_list):
            print "===================================="
            print "ITERATION %d / %d\n" % (num_run+1, runtime)
            if isinstance(parameter, list) or isinstance(parameter, tuple):
                saturation_parameter = parameter[0]
                cutoff_parameter = parameter[1]
            else:
                saturation_parameter = parameter
                cutoff_parameter = default_P_cutoff

            def tweak_function(gas, cloud, C, Cl):
                comment = "Jupiter tweaking"
                nAtm = len(gas[C['P']])
                for i in range(nAtm):
                    Plyr = gas[C['P']][i]
                    Tlyr = gas[C['T']][i]
                    #  Process H2S
                    #  Process NH3
                    if Plyr <= cutoff_parameter:
                        gas[C['NH3']][i] *= saturation_parameter
                    #  Process CO
                    gas[C['CO']][i] = 0.0
                    #  Process CO13
                    gas[C['CO13']][i] = 1.0E-2 * gas[C['CO']][i]
                    #  Process HCN, PH3
                    gas[C['HCN']][i] = 0.0
                    gas[C['PH3']][i] = 0.0
                return comment, gas, cloud

            self.atm.run(tweak_fn=tweak_function,
                         plot=False, verbosity=False)
            self.log.flush()

            self.imSize = None
            self.outType = None
            self.bType = None

            ###Set freqs
            reuse = False
            if freqs == None and freqUnit == None:
                freqs = self.freqs
                freqUnit = self.freqUnit
            elif freqs == 'reuse':
                freqs = self.freqs
                reuse = True
            else:
                if freqUnit == None:
                    freqUnit = 'GHz'
                    self.freqUnit = freqUnit
                if freqs == None:
                    freqs = self.freqs
                else:
                    freqs = self.__freqRequest__(freqs, freqUnit)
                freqUnit = self.freqUnit
            self.freqs = freqs
            self.freqUnit = freqUnit
            self.data_return.f = self.freqs

            ###Set b
            b = copy.deepcopy(original_b)
            if b == None:
                b = self.b
            else:
                b = self.__bRequest__(b, block)
            if self.outType is None or self.outType == 'Spectrum' or self.outType == 'Profile':
                if len(b) > 5 * len(freqs):
                    self.outType = 'Profile'
                else:
                    self.outType = 'Spectrum'
            if outputType == 'batch':
                self.outType = 'Spectrum'
                outputType = 'frequency'
            self.b = b
            self.data_return.b = self.b
            if self.verbosity > 2:
                print 'outType = ' + self.outType
            if self.outType == 'Image' and len(freqs) > 1:
                print 'Warning:  Image must be at only one frequency'
                print 'Using %f %s' % (freqs[0], self.freqUnit)
                self.freqs = list(freqs[0])
                freqs = self.freqs

            ###Start
            runStart = datetime.datetime.now()
            self.Tb = []
            hit_b = []
            btmp = ''
            self.rNorm = None;
            self.tip = None;
            self.rotate = None
            if self.outType == 'Image':  ##We now treat it as an image at one frequency
                if self.verbosity > 1:
                    print 'imgSize = %d x %d' % (self.imSize[0], self.imSize[1])
                self.Tb_img = []
                imtmp = []
                if abs(block[1]) > 1:
                    btmp = '_%02dof%02d' % (block[0], abs(block[1]))
                else:
                    btmp = ''

            if not reuse:
                if not self.config.Doppler:
                    self.bright.layerAbsorption(freqs, self.atm, self.alpha)

            for i, bv in enumerate(b):
                if self.verbosity > 1:
                    print '%d of %d (view [%.4f, %.4f])  ' % (i + 1, len(b), bv[0], bv[1]),
                Tbt = self.bright.single(freqs, self.atm, bv, self.alpha, orientation, plot=not (self.outType == 'Image'),
                                         discAverage=(self.bType == 'disc'))
                if self.bright.path != None and self.rNorm == None:
                    self.rNorm = self.bright.path.rNorm
                if self.bright.path != None and self.tip == None:
                    self.tip = self.bright.path.tip
                if self.bright.path != None and self.rotate == None:
                    self.rotate = self.bright.path.rotate
                if Tbt == None:  # I've now done away with returning None by returning T_cmb in brightness.py (at least I thought so...)
                    Tbt = []
                    for i in range(len(freqs)):
                        Tbt.append(utils.T_cmb)
                else:  # ... so should always go to 'else'
                    hit_b.append(bv)
                self.Tb.append(Tbt)
                self.data_return.Tb.append(Tbt)
                if self.outType == 'Image':
                    imtmp.append(Tbt[0])
                    if not (i + 1) % self.imSize[0]:
                        self.Tb_img.append(imtmp)
                        imtmp = []
            self.log.flush()
            self.data_return.Tb = self.Tb

            ###Write output files (this needs to be compatible with TBfile  -- eventually should incorporate it in there###
            datFile = folder_path + '/%s%05d.dat' % (
                self.planet, num_run
            )
            if self.verbosity > 2:
                print '\nWriting ' + self.outType + ' data to ', datFile
            df = open(datFile, 'w')
            df.write("# %.2f %.4f\n" % (saturation_parameter, cutoff_parameter))
            df.write("# ITERATION %02d PARAMETERS:\n" % num_run)
            df.write("# NH3 Saturation: %.2f pct\n" % (saturation_parameter * 100))
            df.write("# NH3 P Cutoff (lower): %.4f bar\n" % cutoff_parameter)
            self.__setHeader__(self.rNorm)
            self.__writeHeader__(df)
            if self.outType == 'Image':
                for data0 in self.Tb_img:
                    s = ''
                    for data1 in data0:
                        s += '%7.2f\t' % (data1)
                    s += '\n'
                    df.write(s)
            elif self.outType == 'Spectrum':
                fp_lineoutput = open('specoutputline.dat', 'w')
                if outputType == 'frequency':
                    s = '# %s  K@b  \t' % (self.freqUnit)
                elif outputType == 'wavelength':
                    s = '# cm   K@b  \t'
                else:
                    s = '# %s    cm    K@b  \t' % (self.freqUnit)
                for i, bv in enumerate(hit_b):
                    s += '(%5.3f,%5.3f)\t' % (bv[0], bv[1])
                s = s.strip('\t')
                s += '\n'
                df.write(s)
                for i, f in enumerate(freqs):
                    wlcm = 100.0 * utils.speedOfLight / (f * utils.Units[self.freqUnit])
                    if outputType == 'frequency':
                        s = '%.2f\t  ' % (f)
                    elif outputType == 'wavelength':
                        s = '%.4f\t  ' % (wlcm)
                    else:
                        s = '%.2f     %.4f \t ' % (f, wlcm)
                    for j in range(len(hit_b)):
                        s += '  %7.2f  \t' % (self.Tb[j][i])
                    s = s.strip()
                    s += '\n'
                    fp_lineoutput.write(s)
                    df.write(s)
                fp_lineoutput.close()
            elif self.outType == 'Profile':
                if outputType == 'frequency':
                    s = '# b  K@%s \t' % (self.freqUnit)
                elif outputType == 'wavelength':
                    s = '# b  K@cm  \t'
                else:
                    s = '# b  K@GHz,cm  \t'
                for i, fv in enumerate(freqs):
                    wlcm = 100.0 * utils.speedOfLight / (fv * utils.Units[self.freqUnit])
                    if outputType == 'frequency':
                        s += '  %9.4f   \t' % (fv)
                    elif outputType == 'wavelength':
                        s += '  %.4f   \t' % (wlcm)
                    else:
                        s += ' %.2f,%.4f\t' % (fv, wlcm)
                s.strip('\t')
                s += '\n'
                df.write(s)
                bs = []
                for i, bv in enumerate(hit_b):
                    s = '%5.3f %5.3f\t' % (bv[0], bv[1])
                    bs.append(math.sqrt(bv[0] ** 2 + bv[1] ** 2))
                    for j in range(len(freqs)):
                        s += ' %7.2f\t ' % (self.Tb[i][j])
                    s += '\n'
                    df.write(s)
                if plot:
                    plt.figure("Profile")
                    Tbtr = np.transpose(self.Tb)
                    for j in range(len(freqs)):
                        frqs = ('%.2f %s' % (self.freqs[j], self.freqUnit))
                        plt.plot(bs, Tbtr[j], label=frqs)
                    plt.legend()
                    plt.xlabel('b')
                    plt.ylabel('$T_B$ [K]')
            else:
                print 'Invalid outType:  ' + self.outType
            df.close()
            print "===================================="
        return 0

    """
    END Ramsey's Edit
    """

    def __setHeader__(self,intercept):
        if not intercept:  # didn't intercept the planet
            self.header['res'] = '# res not set\n'
            self.header['orientation'] = '# orientation not set\n'
            self.header['aspect'] = '# aspect tip, rotate not set\n'
            self.header['rNorm'] = '# rNorm not set\n'
        else:
            self.header['orientation'] = '# orientation:   %s\n' % (repr(self.config.orientation))
            self.header['aspect'] = '# aspect tip, rotate:  %.4f  %.4f\n' % ((180.0/math.pi)*self.tip, (180.0/math.pi)*self.rotate)
            self.header['rNorm'] = '# rNorm: %f\n' % self.rNorm
            if self.bType:
                self.header['bType'] = '# bType:  %s\n' % self.bType
            if self.outType:
                self.header['outType'] = '# outType:  %s\n' % (self.outType)
                if self.outType == 'Image':
                    self.header['imgSize'] = '# imgSize: %.0f, %.0f\n' % (self.imRow,self.imCol)
                    resolution = (180.0/math.pi)*3600.0*math.atan(abs(self.b[1][0]-self.b[0][0])*self.rNorm/self.config.distance)
                    print 'resolution = ',resolution
                    self.header['res'] = '# res:  %f arcsec\n' % (resolution)
        self.header['gtype'] = '# gtype: %s\n' % (self.config.gtype)
        self.header['radii'] = '# radii:  %.1f  %.1f  km\n' % (self.config.Req,self.config.Rpol)
        self.header['distance'] = '# distance:  %f km\n' % (self.config.distance)
    def __writeHeader__(self,fp):
        for hdr in self.header:
            fp.write(self.header[hdr])

    def __bRequest__(self,b,block):
        """b has a number of options for different bType:
               'points':  discrete number of points
               'line':  radial lines
               'image':  full image
               'stamp':  small image of region
               'disc':  disc-averaged
           b, bType and outType get set"""
        self.bType = None
        self.header['b'] ='# b request:  '+str(b)+'  '+str(block)+'\n'
        
        self.imSize = None
        if type(b) == float:  ## this generates a grid at that spacing and blocking
            self.bType = 'image'
            self.outType = 'Image'
            pb = []
            grid = -1.0*np.flipud(np.arange(b,1.5+b,b))
            grid = np.concatenate( (grid,np.arange(0.0,1.5+b,b)) )
            # get blocks
            bsplit = len(grid)/abs(block[1])
            lastRow = block[0]/abs(block[1])
            if abs(block[1])==1:
                lastRow = 0
            for i in range(bsplit+lastRow):
                ii = i+(block[0]-1)*bsplit
                vrow = grid[ii]
                for vcol in grid:
                    pb.append([vcol,vrow])
            b = pb
            self.imSize = [len(grid),len(b)/len(grid)]
        elif type(b) == str:
            if b.lower() == 'disc':
                b = [[0.0,0.0]]
                self.bType = 'disc'
                self.outType = 'Spectrum'
                print 'Setting to disc-average'
            elif b.lower() == 'stamp':
                self.bType = 'stamp'
                self.outType = 'Image'
                print 'Setting to postage stamp'
                try:
                    bres = float(raw_input('...Input postage stamp resolution in b-units:  '))
                except ValueError:
                    bres = None
                bxmin,bxmax = raw_input('...Input bx_min, bx_max:  ').split(',')
                try:
                    bxmin = float(bxmin)
                    bxmax = float(bxmax)
                except ValueError:
                    bres = None
                bymin,bymax = raw_input('...Input by_min, by_max:  ').split(',')
                try:
                    bymin = float(bymin)
                    bymax = float(bymax)
                except ValueError:
                    bres = None
                if bres:
                    pb = []
                    for x in np.arange(bxmin,bxmax+bres/2.0,bres):
                        for y in np.arange(bymin,bymax+bres/2.0,bres):
                            pb.append([y,x])
                    b = pb
                    xbr = len(np.arange(bymin,bymax+bres/2.0,bres))
                    self.imSize = [xbr,len(b)/xbr]
            else:
                self.bType = b
                self.outType = None
                b = None
        elif len(np.shape(b)) == 1:     ## this makes:
            pb = []
            if len(b) == 2:
                self.bType = 'points'
                self.outType = 'Spectrum'
                pb.append(b)            ##...data at one location
            else:
                self.bType = 'line'
                self.outType = 'Profile'
                angle = b[0]*math.pi/180.0
                del b[0]
                if len(b) == 3:
                    b = np.arange(b[0],b[1]+b[2]/2.0,b[2]) ##...this generates the points as start,stop,step
                for v in b:
                    pb.append([v*math.cos(angle),v*math.sin(angle)])  ##...a line at given angle (angle is first entry)
            b = pb
        else:
            self.bType = b
            self.outType = None
            b = None
        if not b:
            print 'Invalid bType ('+self.bType+'), but soldiering on'
            self.bType = None
        self.b = b
        return b
    
    def __freqRequest__(self,freqs, freqUnit):
        """ Internal processing of frequency list.
               if freqs is a string, it reads frequencies from that file
               if freqs is a scalar, it is made to a list of length 1
               if freqs is a list with first entry a string:
                    'r'ange:  sets range with 3 following values in list as start,stop,step
                        if the step is negative, it is assumed as a log step
                    that is currently the only option
               otherwise the list is used"""
        self.header['freqs'] = '# freqs request: '+str(freqs)+' '+freqUnit+'\n'
        ### Process frequency range "request"
        if type(freqs)==str:
            print 'Loading frequencies from file:  ',freqs
            freqs = np.loadtxt(freqs)
            freqs = list(freqs)  #convert numpy array back to list (for no good reason really)
        elif type(freqs) != list:
            freqs = [freqs]
        elif type(freqs[0])==str:
            if freqs[0][0].lower()=='r':
                #this used to be how to make a range...#elif len(freqs)==3:
                print 'Computing frequency range (len(freqs)==3)'
                fstart=freqs[1]
                fstop=freqs[2]
                fstep=freqs[3]
                freqs=[]
                f=fstart
                if (fstep<0.0):  # do log steps
                    fstep = abs(fstep)
                    while f<fstop:
                        freqs.append(f)
                        f*=fstep
                else:
                    while f<=fstop:
                        freqs.append(f)
                        f+=fstep
            else:
                print 'No other string allowed here yet...'
                freqs=None
        for i in range(len(freqs)):
            freqs[i]*=utils.Units[freqUnit]/utils.Units[utils.processingFreqUnit]
        if len(freqs) > 1:
            utils.log(self.log,self.planet+' in '+str(len(freqs))+' frequency steps ('+str(freqs[0])+' - '+str(freqs[-1])+' '+utils.processingFreqUnit+')',True)
        else:
            utils.log(self.log,self.planet+' at '+str(freqs[0])+' '+utils.processingFreqUnit,True)
        self.freqs = freqs
        self.freqUnit = freqUnit
        return freqs
