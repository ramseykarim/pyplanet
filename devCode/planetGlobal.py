import matplotlib.pyplot as plt

class Jupiter:
    def __init__(self):
        self.Req = 71492.0         # in km
        self.Rpol = 66854.0
        self.RJ = self.Req
        self.GM = 12.6686538e7
        oblate = 0.06487
        self.Jn = [0.0, 0.0, 1.4697e-2, 0.0, -5.84e-4, 0.0, 0.31e-4]
        self.omega_m = 1.7585e-4
        fp = open('zonalJupiter.dat','r')
        self.vwlat = []
        self.vwdat = []
        for line in fp:
            data = line.split()
            self.vwlat.append(float(data[0]))
            self.vwdat.append(float(data[1]))
        plt.plot(self.vwlat,self.vwdat,label='Jupiter')
        plt.title('Zonal Winds')
        plt.xlabel('Lat [deg]')
        plt.ylabel('Speed [m/s]')

class Saturn:
    def __init__(self):
        self.Req = 60268.0         # in km
        self.Rpol = 54364.0
        self.RJ = self.Req
        self.GM = 3.7931272e7
        oblate = 0.09796
        self.Jn = [0.0, 0.0, 1.6332e-2, 0.0, -9.19e-4, 0.0, 1.04e-4]
        self.omega_m = 1.6378499e-4
        fp = open('zonalSaturn.dat','r')
        self.vwlat = []
        self.vwdat = []
        for line in fp:
            data = line.split()
            self.vwlat.append(float(data[0]))
            self.vwdat.append(float(data[1]))
        plt.plot(self.vwlat,self.vwdat,label='Saturn')
        plt.title('Zonal Winds')
        plt.xlabel('Lat [deg]')
        plt.ylabel('Speed [m/s]')

class Uranus:
    def __init__(self):
        self.Req = 25559.0         # in km
        self.Rpol = 24973.0
        self.RJ = self.Req
        self.GM = 0.5793947e7
        oblate = 0.02293
        self.Jn = [0.0, 0.0, 0.3516e-2, 0.0, -0.354e-4, 0.0, 0.0]
        self.omega_m = 1.0833825e-4
        fp = open('zonalUranus.dat','r')
        self.vwlat = []
        self.vwdat = []
        for line in fp:
            data = line.split()
            self.vwlat.append(float(data[0]))
            self.vwdat.append(float(data[1]))
        plt.plot(self.vwlat,self.vwdat,label='Uranus')
        plt.title('Zonal Winds')
        plt.xlabel('Lat [deg]')
        plt.ylabel('Speed [m/s]')

class Neptune:
    def __init__(self):
        self.Req = 24766.0         # in km
        self.Rpol = 24342.0
        self.RJ = self.Req
        self.GM = 0.6835096e7
        oblate = 0.0171
        self.Jn = [0.0, 0.0, 0.3539e-2, 0.0, -0.28e-4, 0.0, 0.0]
        self.omega_m = 1.01237195e-4
        fp = open('zonalNeptune.dat','r')
        self.vwlat = []
        self.vwdat = []
        for line in fp:
            data = line.split()
            self.vwlat.append(float(data[0]))
            self.vwdat.append(float(data[1]))
        plt.plot(self.vwlat,self.vwdat,label='Neptune')
        plt.title('Zonal Winds')
        plt.xlabel('Lat [deg]')
        plt.ylabel('Speed [m/s]')

