# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 10:15:11 2020

@author: farad
"""
import numpy as np
import matplotlib.pyplot as plt
from skimage import io
from os import listdir
from os.path import isfile, join
from scipy.ndimage.morphology import binary_dilation
from mpl_toolkits.mplot3d import Axes3D
from itertools import product

def myint(s,b):
	if s==b:return 0
	else: return 1


########################## MAIN 

mypath = input("What is the path to your masks? ")
cratio = float(input("What is the slice to pixel ratio? "))
d = int(input("What is the size of the region of interest? "))
onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f)) and f[-3::]=='tif' ]
ncells=len(onlyfiles)


C=[] 
Cb=[] 
k = np.ones((3,3),dtype=int)
for f in onlyfiles:
	c=[] 
	b=[] 
	im = io.imread(mypath+'/'+f)
	for z in range(im.shape[0]): 
		I=im[z][:][:]
		I=[[myint(x,y[0]) for x in y] for y in I]	
		I=np.array(I)
		b.append(I)
		Il=binary_dilation(I==0, k) & I
		c.append(Il)
	C.append(c)
	Cb.append(b)


Allx,Ally,Allz,AllN=[],[],[],[]
for n in range(ncells):
    Ly,Lx=len(C[n][0]),len(C[n][0][0])
    X,Y,Z,N=[],[],[],[]
    for z in range(im.shape[0]):
        for i in range(Ly):
            for j in range(Lx):
                if C[n][z][i][j]==1:
                    X.append(j)
                    Y.append(i)
                    Z.append(z)
                    N.append(n)
    Allx=Allx+X
    Ally=Ally+Y
    Allz=Allz+Z
    AllN=AllN+N

T=[]
maxcc=0
for k in range(len(Allx)):
    xx,yy,zz,nn=Allx[k],Ally[k],Allz[k],AllN[k]
    cc=0
    for n in range(ncells):
    	for z,i,j in product(range(max(0,zz-1),min(im.shape[0],zz+1)),range(max(yy-d,0),min(Ly,yy+d)),range(max(0,xx-d),min(Lx,xx+d))):
            if C[n][z][i][j]==1:
                cc=cc+1
                break
    if cc>2:
    	T.append([j,i,z,cc])
    	if cc>maxcc:maxcc=cc


for n in range(3,maxcc+1):
    fi=open(mypath+'/Junctions_'+str(d)+'_'+str(n)+'.txt', 'w')
    for k in range(len(T)):
          if T[k][3]==n: fi.write(str(T[k][2]*cratio)+';'+str(T[k][1])+';'+str(T[k][0])+';'+str(T[k][3])+'\n')
    fi.close()

            
