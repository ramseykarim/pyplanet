###This tests the refractive bending math.  It assumes a perfectly spherical planet, which sets the local normal.
import math
import matplotlib.pyplot as plt

def refract(layers):
    """Refurns refractive index at layers"""
    n = []
    for r in layers:
        v = 1.0 + (50000.0/r)**2
        v = 1
        n.append(v)
    return n

def layers(rmin=10000.0,rmax=20000.0,nlyr=100):
    """Returns layer mid-points and edges.
       Note that:
           layer 0 is free space
           layer n+1 is interior-most layer continued
        so, oddly, there is one more layer than edge"""
    dr = (rmax-rmin)/nlyr
    edges = []
    for i in range(nlyr):
        r = rmax - i*dr
        edges.append(r)
    mid = [rmax+dr]
    for i in range(nlyr-1):
        r = (edges[i]+edges[i+1])/2.0
        mid.append(r)
    mid.append(r)
    return mid, edges

#######################################################
#make layers
mid, edge = layers(rmin=10000.0,rmax=20000.0,nlyr=100)
n = refract(mid)
#######################################################

#######################
# set impact parameter
b = 0.75*edge[0]
#######################

# initialize x,y position of ray outside atmosphere
path = [ [1.001*edge[0]], [b] ]                 # this is path[][0]
#  ... now at the atmosphere edge
path[0].append(math.sqrt(edge[0]**2 - b**2))    # this is path[0][1]
path[1].append(b)                               # this is path[1][1]

# initialize angles
sinThetaInc = b/edge[0]
cosThetaInc = math.sqrt(1.0 - sinThetaInc**2)
sinThetaTran=sinThetaInc    # do this to get the correct first inc angle in loop
dr = 0.0                    # for the first one we want the value at the edge
incAngle = [math.asin(sinThetaInc)]
# initialize s-vector and norm-vector
s = [ [-1.0], [0.0] ]                            # this is s[][0]
na = math.atan2(path[1][1],path[0][1])
norm = [ [math.cos(na)], [math.sin(na)] ]        # this is n[][0]

for i,r in enumerate(edge):
    # get refractive index ratio
    nratio = n[i]/n[i+1]
    print 'x = ',path[0][i+1]
    # get angles
    #sinThetaInc = sinThetaTran*(r+dr)/(r)         # Cassini's model for spherical atmosphere
    sinThetaInc = math.sin(math.atan2(path[1][i+1],path[0][i+1]))
    cosThetaInc = math.sqrt(1.0 - sinThetaInc**2)
    incAngle.append(math.asin(sinThetaInc))
    print 'incident angle ', incAngle[i+1]*180.0/math.pi
    sinThetaTran = nratio*sinThetaInc
    cosThetaTran = math.sqrt(1.0 - sinThetaTran**2)
    # update norm-vector
    na = math.atan2(path[1][i+1],path[0][i+1])  
    norm[0].append(math.cos(na))                  # this is norm[0][i+1]
    norm[1].append(math.sin(na))                  # this is norm[1][i+1]
    # update s-vector
    A = nratio*cosThetaInc - cosThetaTran
    s[0].append(nratio*s[0][i] + A*norm[0][i])    # this is s[0][i+1]
    s[1].append(nratio*s[1][i] + A*norm[1][i])    # this is s[1][i+1]


    ##############################################################
    # update the path
    # need to check for when the ray crosses the tangent point and do something appropriate
    try:
        dr = r - edge[i+1]
    except:
        dr = 0.0
    if dr > 0:
        bq = path[0][i+1]*s[0][i+1] + path[1][i+1]*s[1][i+1]
        try:
            ds = -bq - math.sqrt(bq**2 - 2.0*r*dr + dr**2)
        except:
            ds = 0.0
    else:
        ds = 0.0
        dr = 0.0
    path[0].append(path[0][i+1] + ds*s[0][i+1])    # this is path[0][i+2]
    path[1].append(path[1][i+1] + ds*s[1][i+1])    # this is path[1][i+2]
    ##############################################################

plt.plot(path[0],path[1],'-o')
v = plt.axis()


circle = []
for i,r in enumerate(edge):
    circle.append([[],[]])
    for th in range(180):
        angle = (th-30.0)*math.pi/180.0
        circle[i][0].append(r*math.cos(angle))
        circle[i][1].append(r*math.sin(angle))
    plt.plot(circle[i][0],circle[i][1],color='black')


#plt.axis([v[0],v[1],v[2]*0.99,v[3]*1.01])
        
