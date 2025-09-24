import os, math, re, sys
import struct
import numpy as np
from matplotlib import pyplot as plt
from scipy import optimize
import matplotlib

datadir = "data2/"

HeaderFile = open(datadir+"Header.mat", 'rb')

LX=struct.unpack('=i', HeaderFile.read(4))[0]

LY=struct.unpack('=i', HeaderFile.read(4))[0]

LZ=struct.unpack('=i', HeaderFile.read(4))[0]

ndim=struct.unpack('=i', HeaderFile.read(4))[0]

t_zero = 0
tstart = 0

tend = 350000
struct.unpack('=i', HeaderFile.read(4))[0]
tinc = struct.unpack('=i', HeaderFile.read(4))[0]

numcomp = 3

def calc_R(xc, yc,x,y):
    """ calculate the distance of each 2D points from the center (xc, yc) """
    return np.sqrt((x-xc)**2 + (y-yc)**2)

def f_2(c,x,y):
    """ calculate the algebraic distance between the data points and the mean circle centered at c=(xc, yc) """
    Ri = calc_R(*c,x,y)
    #print(Ri)
    return Ri - Ri.mean()

def circleFit(c,LX,LY,h1,h2,w1,w2,offset):
    x=np.where(np.logical_and(c[w1:w2,h1:h2]>=0.45, c[w1:w2,h1:h2]<=0.55))[0]
    y=np.where(np.logical_and(c[w1:w2,h1:h2]>=0.45, c[w1:w2,h1:h2]<=0.55))[1]

    x_m=np.mean(x)
    y_m=np.mean(y)
    #print(x_m,y_m)
    center_estimate = x_m, y_m
    center_2, ier = optimize.leastsq(f_2, center_estimate,args=(x,y))

    xc_2, yc_2 = center_2
    Ri_2       = calc_R(*center_2,x,y)
    R_2        = Ri_2.mean()
    residu_2   = sum((Ri_2 - R_2)**2)
    #print(xc_2,yc_2)
    theta=90-np.arcsin(abs(yc_2-offset)/R_2)*180/np.pi
    print(theta)
    return theta,R_2,x,y,xc_2,yc_2

outDirName = "figures"
os.system("mkdir -p %s"%outDirName)
v = np.zeros((LX,LY,LZ,ndim))
for t in range(tstart,tend+1,tinc):

    t_file =t+t_zero

    file_name = datadir+"Velocity_t%li.mat"%t_file
    File2 = open(file_name, 'rb')
    v = np.zeros((LX,LY,LZ,ndim))
    dat=File2.read()
    v = np.ndarray((LX,LY,ndim),'=d',dat,0,(ndim*8*LY,ndim*8,8))

    C = np.zeros((LX,LY,2))
    for iC in range(2):
        file_name = f"{datadir}OrderParameter{iC}_t{t_file}.mat"
        File = open(file_name, 'rb')
        dat=File.read()
        File.close()
        C[:,:,iC] = np.frombuffer(dat, '=d').reshape(LX, LY)

    liquid = C[:,:,0] + 0.5*C[:,:,1]
    print("v",np.amax(v[:,:,0]))
    print("c",np.sum(C[:,:,0]))
    fig,ax=plt.subplots(1,1,figsize=(6,6))
    output = "%s/component_plot_%012d.png"%(outDirName,t)

    rgbv = np.zeros((LY,LX))
    rgbv[:,:] = np.flip(liquid).transpose()

    im=ax.imshow(rgbv,interpolation='nearest',origin='upper')
    fig.colorbar(im)
    plt.savefig(output, dpi=400, format='png')
    plt.close(fig)
