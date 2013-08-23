import h2s_ddb
import matplotlib.pyplot as plt

f = []
for i in range(1000):
    f.append(float(i)/2.0)

P = 6.0
T = 300.0
X_h2s = 3.3E-5
X_he = 0.1
X_h2 = 1.0 - (X_h2s+X_he)
X = [X_h2,X_he,X_h2s]
P_dict = {'H2':0,'HE':1,'H2S':2}
a = h2s_ddb.alpha(f,T,P,X,P_dict)
plt.plot(f,a)
