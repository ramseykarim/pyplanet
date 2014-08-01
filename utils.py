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
T_cmb = 2.725

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

def concatdat(files,directory='Output'):
    """Given a list of utils.ls indices, returns TB, freq, wavel"""
    aTB = []
    bf = []
    cwavel = []
    for fil in files:
        data = readTB(fn=fil)
        a = list(data[1][0])
        b = list(data[2])
        c = list(data[3])
        aTB+=a
        bf+=b
        cwavel+=c
    sa = np.argsort(np.array(bf))
    TB = []
    f = []
    wavel = []
    for i in sa:
        TB.append(aTB[i])
        f.append(bf[i])
        wavel.append(cwavel[i])
    
    return TB,f,wavel

def plotObs(fn,cols=[0,1,2], color='b',marker='o',delimiter=None,comline='!'):
    try:
        fp=open(fn,'r')
    except IOError:
        print fn+' not found'
        return 0
    data = []
    for line in fp:
        if comline in line[0:5]:
            continue
        dline = line.split(delimiter)
        if len(dline)<max(cols):
            continue
        drow = []
        for c in cols:
            drow.append(float(dline[c]))
        data.append(drow)
    data = np.array(data)
    plt.semilogx(data[:,0],data[:,1],color=color,marker=marker)
    plt.errorbar(data[:,0],data[:,1],yerr=data[:,2],color=color,marker=marker,ls='none')
    
        
