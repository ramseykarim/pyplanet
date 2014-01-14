import os.path
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import os
import numpy as np

commentChars = ['!','#','$','%','&','*']
affirmative = [1,'1','y','Y','t','T']

Units = {'Frequency':1,'Hz':1.0,'kHz':1.0E3,'MHz':1.0E6,'GHz':1.0E9,
         'Length':1, 'm':1.0, 'km':1.0E3,'cm':1.0E-2,'AU':149597870691.0,
         'Pressure':1, 'bars':1.0, 'atm':1.01325,
         'Time':1, 'sec':1.0, 'min':60.0, 'hr':3600.0, 'day':86400.0, 'year':31536000.0,
         'Acceleration':1, 'mpersec2':1.0, 'cmpersec2':0.01}
processingFreqUnit = 'GHz'  	#alpha assumes this unit
processingAtmLayerUnit = 'km'   #alpha assumes this unit for layers
processingPressureUnit = 'bars'
processingAccelUnit = 'mpersec2'
speedOfLight = 2.9979E8     	# m/s
kB = 1.3806503E-23          	# m2 kg s-2 K-1  Boltzmann's constant
hP = 6.626068E-34       	# m2 kg / s	 Planck's constant

rfBands = {'HF':[0.003,0.03],'VHF':[0.03,0.3],'UHF':[0.3,1.0],'L':[1.0,2.0],'S':[2.0,4.0],'C':[4.0,8.0],'X':[8.0,12.0],'Ku':[12.0,18.0],
           'K':[18.0,26.5],'Ka':[26.5,40.0],'Q':[40.0,50.0],'V':[50.0,75.0],'W':[75.0,110.0]}
def getRFband(freq,unit='GHz'):
    freq = freq*Units[unit]/Units['GHz']
    for bnd in rfBands:
        if rfBands[bnd][0] <= freq < rfBands[bnd][1]:
            return bnd
    return None

def setupLogFile(log,path='Logs/'):
    if type(log) == file:
        logfp = log
    elif type(log) == str:
        lf = os.path.join(path,log)
        try:
            logfp = open(lf,'a')
        except IOError:
            print lf+' not found.  No logging (and save the whales).'
            return None
        print 'Logging to '+lf
    else:
        log = None
        logfp = None
    return logfp

def log(logfp,msg,printOut=True):
    if logfp:
        logfp.write(msg+'\n')
        if printOut:
            print msg

def close(logfp):
    logfp.close()

def ls(directory='Output', show=True, returnList=False):
    """Generates file list for plotTB and writeWavel"""
    filelist = os.listdir(directory)
    files = []
    i = 0
    for fff in filelist:
        if fff[0] != '.':
            files.append(os.path.join(directory,fff))
            if show:
                print '%d:  %s' % (i,files[i])
            i+=1
    if returnList:
        return files

def bang():
    files = ls(show=False,returnList=True)
    for f in files:
        plotTB(f,xaxis='wavel',xlog=True,justFreq=True)

def plotTB(fn=None,xaxis='Frequency',xlog=False, justFreq=False, directory='Output',distance=4377233696.68):
    """plots brightness temperature against frequency and disc location:
           fn = filename to read (but then ignores directory) | '?', '.' or None | integer [None]
           xaxis = 'f[requency]' | 'w[avelength' ['freq']
           xlog = True | False [False]
           justFreq = True | False [False]
           directory = subdirectory for data (not used if filename given) ['Output']
           distance = distance for angular size plot in km [4377233696 km for Neptune]"""

    filename,Tb,f,wavel,b,xlabel,ylabels = readTB(fn=fn,directory=directory)
    title = filename.split('/')

    ## Frequency plot
    plt.figure('TB')
    for i in range(len(b)):
        if xaxis[0].lower() == 'f':
            plotx = f
        else:
            plotx = wavel
            xlabel='Wavelength [cm]'
        if xlog:
            plt.semilogx(plotx,Tb[i],label=ylabels[i])
        else:
            plt.plot(plotx,Tb[i],label=ylabels[i])
        #plt.plot(f,Tb[i],'x')
    #plt.legend()
    plt.xlabel(xlabel)
    plt.ylabel('Brightness Temperature [K]')
    plt.title(title[-1])

    if justFreq:
        return len(f)

    ## b plot
    plt.figure('b')
    for i in range(len(f)):
        plt.plot(b[:,0],Tb[:,i],label=str(f[i]))
        plt.plot(b[:,0],Tb[:,i],'o')
    plt.legend()
    plt.title(title[-1])
    plt.xlabel('km')
    plt.ylabel('Brightness Temperature [K]')

    ## b plot vs angle
    angle = []
    for r in b[:,0]:
        angle.append( (r/distance)*(180.0/np.pi)*3600.0 )
    plt.figure('b_vs_angle')
    for i in range(len(f)):
        plt.plot(angle,Tb[:,i],label=str(f[i]))
        plt.plot(angle,Tb[:,i],'o')
    plt.legend()
    plt.title(title[-1])
    plt.xlabel('arcsec')
    plt.ylabel('Brightness Temperature [K]')

    return len(b)

def writeWavel(fn=None,outputFile=None,directory='Output'):

    filename,Tb,f,wavel,b,xlabel,ylabels = readTB(fn=fn,directory=directory)
    title = filename.split('/')[-1].split('.')[0]


    ## Write file
    title+='_wavel.dat'
    if outputFile == None:
        outputFile = title
    print 'Writing to '+outputFile
    fp=open(outputFile,'w')
    for i in range(len(wavel)):
        s = '%f\t' % (wavel[i])
        fp.write(s)
        for j in range(len(b)):
            s = '%f\t' % (Tb[j][i])
            fp.write(s)
        fp.write('\n')
    fp.close()
            
    return n

def readTB(fn=None,directory='Output'):
    """reads brightness temperature file:
           fn = filename to read (but then ignores directory) | '?', '.' or None | integer [None]
           directory = subdirectory for data (not used if filename given) ['Output']
           """
    print '==>I should make TB files a class...'
    filename = fn
    if fn=='?' or fn=='.' or fn==None:
        files = ls(directory=directory, show=True, returnList=True)
        i = input('File number:  ')
        filename = files[i]
    elif type(fn) == int:
        files = ls(directory=directory, show=False, returnList=True)
        filename = files[fn]
    print '\nOpening '+filename
    
    try:
        fp = open(filename,'r')
    except IOError:
        print filename+' not found'
        return 0

    ## Get past any header and get first line
    headerText = []
    line ='# file:  '+filename+'\n'
    while line[0]=='#':
        headerText.append(line)
        line = fp.readline()
    print '=============Header (note update to use img.parseHeader)============'
    for hdr in headerText:
        print hdr,
    print '\n\n'
        
    bvals = []
    labels = line.split()
    xlabel = labels[0]
    del(labels[0])
    print 'b = ',
    for b in labels:
        print ' '+b,
        bb = b.split('(')[1].strip(')').split(',')
        bb = [float(bb[0]),float(bb[1])]
        bvals.append(bb)
    b = np.array(bvals)
    print ''
    ylabels = labels

    ## Rest of data
    print '>>>>>>>>>>>>>>>>wavelength assumes cm/GHz for now<<<<<<<<<<<<<<<<'
    f = []
    Tb = []
    n = 0
    wavel = []
    global funit
    for line in fp:
        data = line.split()
        f.append(float(data[0]))
        wavel.append((speedOfLight/1E7)/float(data[0]))
        del(data[0])
        t = []
        for d in data:
            t.append(float(d))
        Tb.append(t)
        n+=1
    Tb = np.array(Tb).transpose()
    print 'Freq:  %.3f - %.3f %s  (%d points)' % (f[0],f[-1],xlabel,n)

    return filename,Tb,f,wavel,b,xlabel,ylabels
