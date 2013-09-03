import planet
n = planet.planet('neptune')
frq = 43.34 
bstep = 0.01
sizeofblock = 15
#thisrun = [11]
#thisrun = [1,2,3,4]
thisrun = [12,13,14,15]

for i in thisrun:
    n.run(freqs=frq,b=bstep,block=[i,sizeofblock])

