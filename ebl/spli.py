import numpy as np
import argparse
parser = argparse.ArgumentParser(description='''get opacity for z=0-1 
  for example python spli.py -z=0.1 -filename=z=0.01.dat
''',formatter_class=argparse.RawTextHelpFormatter) 
parser.add_argument('-z', type=float, help='red-shift')
parser.add_argument('-filename', type=str, help='outfilename')
args = parser.parse_args()

ifile=open("Gamma-gamma-opacity-z=0-1.dat")
z=args.z
filename=args.filename
xy=np.zeros((50,2))
xy1=np.zeros((50,2))
zt=0
zt1=0
for i in range(1000):
    xy=np.copy(xy1)
    zt=zt1
    line=ifile.readline().split()
    zt1=float(line[3])
    ifile.readline()
    for j in range(50):
        line=ifile.readline().split()
        xy1[j,0]=line[0]
        xy1[j,1]=line[3]
    if(abs(zt1-z)<0.00001):
        print("ppp")
        np.savetxt(filename,xy1,fmt='%.8e',header='z= %f' %(z))
        exit()
    if(zt1>z):
        res=xy
        y=xy[:,1]
        y1=xy1[:,1]
        res[:,1]=y+(y1-y)*((z-zt)/(zt1-zt))
        np.savetxt(filename,res,fmt='%.8e',header='z= %f' %(z))
        print(zt1,z,abs(zt1-z))
        exit()
