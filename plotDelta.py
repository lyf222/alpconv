
import numpy as np

kl=2*np.pi/35
kh=2*np.pi/0.7
q=-2.8
rmax=500.
sigamb=10.0
yita=0.5

nk =int(20*(np.log10(kh/kl)+3))  # nk is the total number of spacings in k; where k is the wave number
stepk=pow(10,(np.log10(kh/kl)+3)/nk)
stepk2 = (np.log10(kh/kl)+3)/nk

kmin = 1.e-3*kl
lmin = kmin

print(kmin,kl,kh,nk,stepk)

for i in range(nk):
    k=lmin*pow(stepk,i)
    k2=10**(np.log10(lmin)+stepk2*i)
    deltk=k*(stepk-1)
    print(i,k,k2,deltk,np.log10(k))

# ma(neV),e(KeV),ne(cm^-3) gag(10^-11 GeV^-1),b(10^-6 G) deta (Kpc)-1
ga = 1 #22.62
ma = 1 #6.1

ne = 0.1
b = 1.
e = np.logspace(np.log10(0.1),np.log10(1.e4),100) # 0.1 GeV - 10 TeV
e = e*1.e6  # GeV to keV
dpl = -1.1e-4/(e/100)*ne/1.0e-7*1.0e-3
dqed = 4.1e-16*(e/100)*(b/1.0e-3)*(b/1.0e-3)*1.0e-3
dgg = 8.0e-2*(e/1.0e9)*1.0e-3

dag=1.52e-2*(ga)*(b/1.0e-3)*1.0e-3
da=-7.8e-3*ma*ma/1.0e-8/(e/100)*1.0e-3

#print dpl,dqed,dgg

import matplotlib.pyplot as plt
plt.plot(e,-dpl,label='-dpl')
plt.plot(e,dqed,label='dqed')
plt.plot(e,dgg,label='dgg')
plt.plot(e,dag*np.ones_like(e),label='dag')
plt.plot(e,-da,label='da')
plt.xscale('log')
plt.yscale('log')
plt.ylim(1.e-10,10)
plt.legend()
plt.savefig('test.png')