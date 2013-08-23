import nh3_kd
import matplotlib.pyplot as plt
import math

f = []
fmin = 1.0
fmax = 500.0
fstep = 1.0
n = int(math.ceil((fmax-fmin)/fstep))+1
for i in range(n):
    f.append(fmin + i*fstep)

print 'Fig. 5.2'
P = 1.009
T = 216.4
X_nh3 = 0.0095
X_he = 0.1347
X_h2 = 1.0 - (X_nh3+X_he)
X_partial=[X_h2,X_he,X_nh3]
P_dict = {'H2':0,'HE':1,'NH3':2}
a = nh3_kd.alpha(f,T,P,X,P_dict)
plt.semilogy(f,a)
axis([1,25,.01,1000])
