#! /usr/bin/env python
#script planet processing
import planet
import sys

freqs = [1.0]
j = planet.planet('jupiter')
print '==================='
j.atm.readGas()
print '++++++++++++++++++++'
j.atm.readCloud()
print '11111111111111111111'
j.atm.tweakAtm()
j.atm.computeProp()
###need to make sure regrid called correctly...
j.atm.regrid.regrid()


#j.run(freqs)
