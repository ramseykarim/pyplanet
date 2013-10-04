###This is the 'executive' class for planets
import math
import string
import matplotlib.pyplot as plt
import numpy as np
from scipy import special
import atmosphere as atm
import alpha
import brightness as bright
import utils
import datetime
import os.path

version = 'pre-release'
header = {'freqs':'', 'b':'', 'res':'', 'gtype':'', 'orientation':'', 'aspect':'', 'radii':'', 'rNorm':'', 'distance':''}

class planet:
    def __init__(self, name, freqs=None, b=None, freqUnit='GHz', config='config.par', log='auto', verbose=False, plot=True):
        """This is the 'executive function class to compute overall planetary emission
           Inputs:
               planet:  'Jupiter', 'Saturn', 'Uranus', 'Neptune'
               freqs: options are:
                  a value:  does that one frequency
                  a triple:  assumes it is [start,stop,step]
                  a list: does those frequencies - don't use a length 3 list (see above)
               freqUnit: unit for above
               config:  config file name [config.par], 'manual' [equivalent none]
               log:  either a file name, a 'no', or 'auto' (for auto filename)
               verbose:  True/False
               plot:  True/False"""

        if name == 'functions':
            return

        self.planet = string.capitalize(name)
        runStart = datetime.datetime.now()

        print 'Planetary modeling  (ver '+version+')\n'
        print "In alpha, clouds_idp need otherPar['refr']"

        ### Set up log file
        if string.lower(log)=='auto':
            self.logFile = '%s_%d%02d%02d_%02d%02d.log' % (self.planet,runStart.year,runStart.month,runStart.day,runStart.hour,runStart.minute)
        elif string.lower(log)=='no':
            self.logFile=None
        else:
            self.logFile = log
        self.log=utils.setupLogFile(self.logFile,path='Logs/')
        utils.log(self.log,self.planet+' start '+str(runStart),True)

        ### Some convenience values
        self.fvla = [4.86,8.46,14.94,22.46,43.34]
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

        if config == 'manual' or config=='none':
            config = None
        else:
            config = os.path.join(self.planet,config)
        self.config = config
        self.plot = plot
        self.verbose = verbose

        ### Create atmosphere:  outputs are self.atm.gas, self.atm.cloud and self.atm.layerProperty
        self.atm = atm.atmosphere(self.planet,config=config,log=self.log,verbose=verbose,plot=plot)
        self.atm.run()
        self.log.flush()

    def run(self, freqs=[1.0,10.0,1.0], b=0.04175, freqUnit='GHz', orientation=None, block=[1,1], verbose=None, plot=None):
        if freqs == None and freqUnit == None:
            freqs = self.freqs
            freqUnit = self.freqUnit
        else:
            if freqUnit == None:
                freqUnit = 'GHz'
                self.freqUnit = freqUnit
            if freqs == None:
                freqs = self.freqs
            else:
                freqs = self.__freqRequest__(freqs, freqUnit)
            freqUnit = self.freqUnit
        if b == None:
            b = self.b
        else:
            b = self.__bRequest__(b,block)
        if verbose == None:
            verbose = self.verbose
        if plot == None:
            plot = self.plot
        runStart = datetime.datetime.now()

        ### Read in absorption modules:  to change absorption, edit files under /constituents'
        self.alpha = alpha.alpha(config=self.config,log=self.log,verbose=verbose,plot=plot)
        #self.alpha.test(f=0,verbose=True,plot=True)
        self.log.flush()

        ### Next compute radiometric properties
        self.bright = bright.brightness(log=self.log,verbose=verbose,plot=plot)
        if not self.alpha.Doppler:
            self.bright.layerAbsorption(freqs,self.atm,self.alpha)
        self.Tb=[]
        hit_b=[]
        isImg = False
        self.rNorm = None; self.tip = None; self.rotate = None
        if len(freqs)==1 and self.imRow:
            isImg = True
            imCol = len(b)/self.imRow
            print 'imgSize = %d x %d' % (self.imRow, imCol)
            self.Tb_img = []
            imtmp = []
        for i,bv in enumerate(b):
            print '%d of %d:  ' % (i,len(b)),
            Tbt = self.bright.single(freqs,self.atm,bv,self.alpha,orientation,isImage=isImg,discAverage=self.discAverage)
            if self.bright.path != None and self.rNorm == None:
                self.rNorm = self.bright.path.rNorm
            if self.bright.path != None and self.tip == None:
                self.tip = self.bright.path.tip
            if self.bright.path != None and self.rotate == None:
                self.rotate = self.bright.path.rotate
            if Tbt == None:
                Tbt = [0.0]
            else:
                hit_b.append(bv)
                self.Tb.append(Tbt)
            if isImg:
                #if Tbt[0]>250.0:
                #    Tbt[0]=0.0
                imtmp.append(Tbt[0])
                if not (i+1)%self.imRow:
                    self.Tb_img.append(imtmp)
                    imtmp = []
        self.log.flush()

        global header
        if isImg:
            if abs(block[1])>1:
                btmp = '%02dof%02d'%(block[0],abs(block[1]))
            else:
                btmp = ''
            datFile = 'Output/%s_Image%s_%d%02d%02d_%02d%02d.dat' % (self.planet,btmp,runStart.year,runStart.month,runStart.day,runStart.hour,runStart.minute)
            print '\nWriting image data to ',datFile
            self.__setHeader__(self.rNorm)
            df = open(datFile,'w')
            for hdr in header:
                df.write(header[hdr])
            for data0 in self.Tb_img:
                s = ''
                for data1 in data0:
                    s+= '%.4f\t' % (data1)
                s+='\n'
                df.write(s)
            df.close()
        else:
            datFile = 'Output/%s_%d%02d%02d_%02d%02d.dat' % (self.planet,runStart.year,runStart.month,runStart.day,runStart.hour,runStart.minute)
            print '\nWriting data to ',datFile
            matchFile = 'match.dat'
            print '\nWriting data to ',matchFile
            matchfp = open(matchFile,'w')
            self.__setHeader__(None)
            df = open(datFile,'w')
            for hdr in header:
                df.write(header[hdr])
            s = 'GHz \tK@km'
            sm = 'b:  '
            for i,bv in enumerate(hit_b):
                s+='(%.0f,%.0f)\t' % (self.rNorm*bv[0],self.rNorm*bv[1])
                cma = ', '
                if i==0:
                    cma = ''
                sm+='%s%.3f'% (cma,math.sqrt(bv[0]**2. + bv[1]**2.))
            s+='\n'
            sm+='\n'
            df.write(s)
            matchfp.write(sm)
            for i,f in enumerate(freqs):
                s = '%.9f\t' % (f)
                sm = '%.5f:  ' % (f)
                for j in range(len(b)):
                    s+='  %.2f\t  ' % (self.Tb[j][i])
                    cma = ', '
                    if j==0:
                        cma=''
                    sm+='%s%.2f' % (cma,self.Tb[j][i])
                print s
                s+='\n'
                sm+='\n'
                df.write(s)
                matchfp.write(sm)
            df.close()
            matchfp.close()
    def __setHeader__(self,intercept):
        global header
        if not intercept:  # didn't intercept the planet
            header['res'] = '# res not set\n'
            header['orientation'] = '# orientation not set\n'
            header['aspect'] = '# aspect tip, rotate not set\n'
            header['rNorm'] = '# rNorm not set\n'
        else:
            resolution = (180.0/math.pi)*3600.0*math.atan(abs(self.b[1][0]-self.b[0][0])*self.rNorm/self.atm.distance)
            print 'resolution = ',resolution
            header['res'] = '# res:  %f arcsec\n' % (resolution)
            header['orientation'] = '# orientation:   %s\n' % (repr(self.atm.orientation))
            header['aspect'] = '# aspect tip, rotate:  %.4f  %.4f\n' % ((180.0/math.pi)*self.tip, (180.0/math.pi)*self.rotate)
            header['rNorm'] = '# rNorm: %f\n' % self.rNorm
        header['gtype'] = '# gtype: %s\n' % (self.atm.gtype)
        header['radii'] = '# radii:  %.1f  %.1f  km\n' % (self.atm.Req,self.atm.Rpol)
        header['distance'] = '# distance:  %f km\n' % (self.atm.distance)
        
    def bRequest(self,b,block):
        b=self.__bRequest__(b,block)
        return b
    def __bRequest__(self,b,block):
        self.discAverage = False
        global header
        header['b'] ='# b request:  '+str(b)+'  '+str(block)+'\n'
        self.imRow = None
        if type(b) == float:  ## this generates a grid at that spacing
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
            self.imRow = len(grid)
        elif type(b) == str:     ## this assumes a disc averaged value is desired
            b = [[0.0,0.0]]
            self.discAverage = True
            print 'Setting to disc-average'
        elif len(np.shape(b)) == 1:   ## this makes a line along the equator
            pb = []
            if len(b) == 2:
                pb.append(b)
            else:
                for v in b:
                    pb.append([v,0.0])
            b = pb
        else:
            if len(b) > 8:  ## I don't remember...
                self.imRow = int(math.sqrt(len(b)))
        self.b = b
        if block[1]<0:
            print '---------------------------------------------------------'
            print b
            print '---------------------------------------------------------'
        return b
    def __freqRequest__(self,freqs, freqUnit):
        """ Internal processing of frequency list.
               if there is a scalar, it is made to a list of length 1
               if there is a list of three it is assumed to be start, stop, step
                    if the step is negative, it is assumed as a log step
               if it is a list!=3, the list is used"""
        global header
        header['freqs'] = '# freqs request: '+str(freqs)+' '+freqUnit+'\n'
        ### Process frequency range "request"
        if type(freqs) != list:
            freqs = [freqs]
        elif len(freqs)==3:
            fstart=freqs[0]
            fstop=freqs[1]
            fstep=freqs[2]
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
        for i in range(len(freqs)):
            freqs[i]*=utils.Units[freqUnit]/utils.Units[utils.processingFreqUnit]
        if len(freqs) > 1:
            utils.log(self.log,self.planet+' in '+str(len(freqs))+' frequency steps ('+str(freqs[0])+' - '+str(freqs[-1])+' '+utils.processingFreqUnit+')',True)
        else:
            utils.log(self.log,self.planet+' at '+str(freqs[0])+' '+utils.processingFreqUnit,True)
        self.freqs = freqs
        self.freqUnit = freqUnit
        return freqs
