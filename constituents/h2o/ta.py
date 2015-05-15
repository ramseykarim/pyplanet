import h2o_ddb
import matplotlib.pyplot as plt

f = []
for i in range(1000):
    f.append(float(i)/2.0)

P = 6.0
T = 300.0
X_h2o = 3.3E-5
X_he = 0.1
X_h2 = 1.0 - (X_h2o+X_he)
X = [X_h2, X_he, X_h2o]
P_dict = {'H2':0,'HE':1,'H2O':2}
otherpar = []
a = h2o_ddb.alpha(f,T,P,X,P_dict,otherpar)
plt.plot(f,a)
