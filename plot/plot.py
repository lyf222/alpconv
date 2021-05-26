import matplotlib.pyplot as plt
import numpy as np

ee,pp = np.loadtxt('C:\\Users\\Lenovo\\Desktop\\codes\\alp\\rs_par1.txt').T
plt.clf()
plt.plot(ee,pp)
plt.xscale('log')
#plt.ylim(0.75,1.1)
plt.show()

r1,br,r2,bth,r3,bph = np.loadtxt('Bicm_reg.txt').T
plt.clf()
plt.plot(r1,br,label='br')
plt.plot(r2,bth,label='btheta')
plt.plot(r3,bph,label='bphi')
plt.legend()
plt.show()

d = np.loadtxt('matrixoutput.txt').T
xx,bu,bv,bb = d
plt.clf()
plt.plot(xx,bb)
plt.show()


dd = np.loadtxt('rm.txt')
bins = np.linspace(-7500,7500,15)
yy,xx = np.histogram(dd,bins=bins)
x0 = (xx[:-1]+xx[1:])*0.5
plt.clf()
plt.plot(x0,yy)
plt.show()


plt.clf()
d = np.loadtxt('dd.txt').T
plt.plot(d[0],d[1])
plt.plot(d[0],d[2])
plt.plot(d[0],d[3])
plt.show()
