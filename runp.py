#! /usr/bin/env python
#script planet processing
import planet
import sys

oldFreqs = [8.46,14.94,22.46,43.34]
newFreqs = [2.44,3.44,4.42,5.45,6.45,7.45,5.90,3.04,1.46,8.58,9.61,10.38,11.41,13.18,14.21,15.18,16.21,17.38,22.45,23.45,24.45,25.45]
newFreqs.sort()
freqs = oldFreqs

print 'Reading input file ',sys.argv[1]

j = planet.planet('jupiter')
j.atm.readGas(sys.argv[1])
j.atm.readCloud(sys.argv[1])

j.atm.tweakAtm()
j.atm.computeProp()
j.atm.regrid.regrid(j.atm,regridType=j.atm.config.regridType,Pmin=j.atm.config.Pmin,Pmax=j.atm.config.Pmax)
j.atm.nAtm=len(j.atm.gas[0])

j.run(freqs)
