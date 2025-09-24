import os, math, re, sys
import struct
import numpy as np
from matplotlib import pyplot as plt
from scipy import optimize
import matplotlib
import glob

datadir = "data/diffusivity_0.36/pDiff_0.0001/lx_100/ly_32/lz_1/Csat_0.7/"

HeaderFile = open(datadir+"Header.mat", 'rb')

LX=struct.unpack('=i', HeaderFile.read(4))[0]

LY=struct.unpack('=i', HeaderFile.read(4))[0]

LZ=struct.unpack('=i', HeaderFile.read(4))[0]

ndim=struct.unpack('=i', HeaderFile.read(4))[0]

t_zero = 0
tstart = 0

tend = struct.unpack('=i', HeaderFile.read(4))[0]
tinc = struct.unpack('=i', HeaderFile.read(4))[0]

outDirName = "figures"
os.system("mkdir -p %s"%outDirName)
v = np.zeros((LX,LY,LZ,ndim))

plot=True
plotnameheight=False

for t in range(tstart,tend+1,tinc):

    print("t=%s"%t)
    t_file =t+t_zero

    file_name = datadir+"Density_t%li.mat"%t_file

    File = open(file_name, 'rb')

    file_name = datadir+"Concentration_t%li.mat"%t_file

    File0 = open(file_name, 'rb')

    file_name = datadir+"BoundaryLabels_t%li.mat"%t_file

    FileSolid = open(file_name, 'rb')
    
    file_name = datadir+"Velocity_t%li.mat"%t_file

    File2 = open(file_name, 'rb')

    print(file_name)

    def coord_k(k, LY, LZ):
        """From a k value, determines its xk, yk, and zk."""    

        xk = math.floor(k/(LY*LZ))
        yk = math.floor((k - xk*LZ*LY)/LZ)
        zk = k - xk*LZ*LY - yk*LZ
        return xk, yk, zk

    NLatt=LX*LY*LZ

    rho = np.zeros((LX,LY,LZ))
    rho2 = np.zeros((LX,LY,LZ))
    v = np.zeros((LX,LY,LZ,ndim))
    solid = np.zeros((LX,LY,LZ))

    dat=File.read()
    rho = np.ndarray((LX,LY,LZ),'=d',dat,0,(8*LY*LZ,8*LZ,8))

    dat=File0.read()
    concentration = np.ndarray((LX,LY,LZ),'=d',dat,0,(8*LY*LZ,8*LZ,8))

    dat=FileSolid.read()
    solid = np.ndarray((LX,LY,LZ),'=i',dat,0,(4*LY*LZ,4*LZ,4))
    slicezz = (1*LZ)//4
    liquid = np.array(rho[:,:,slicezz])
    liquid[np.where(np.logical_or(solid[:,:,slicezz]==1,solid[:,:,slicezz]==-1))[0],np.where(np.logical_or(solid[:,:,slicezz]==1,solid[:,:,slicezz]==-1))[1]] = 0

    concentration = np.array(concentration[:,:,:])
    concentration[np.where(solid!=0)[0],np.where(solid!=0)[1],np.where(solid!=0)[2]] = 0.1

    dat=File2.read()
    v = np.ndarray((LX,LY,LZ,ndim),'=d',dat,0,(ndim*8*LY*LZ,ndim*8*LZ,ndim*8,8))

    File.close()
    File0.close()
    File2.close()

    output = "%s/component_plot_%012d.png"%(outDirName,t)

    if plot==True:
        
        rgbv = np.zeros((LY,LX))
        rgbv[:,:] = np.flip(concentration[:,:LY,slicezz]).transpose()

        if not(os.path.exists("figures/") and os.path.isdir("figures/")):
            os.mkdir("figures/")
        
        fig,ax=plt.subplots(1,1,figsize=(6,int(6*LY//LX)+1))
        im=ax.imshow(rgbv,interpolation='nearest',origin='upper')
        
        stepx=1
        stepy=1
        print(LX,LY)
        
        X,Y=np.meshgrid(np.linspace(0,LX-1,int((LX)/stepx)),np.linspace(0,LY-1,int((LY)/stepy)))
        ax.quiver(X.T,Y.T,np.flip(-v[0:LX:stepx,0:LY:stepy,0,0]),np.flip(v[0:LX:stepx,0:LY:stepy,0,1]),width=0.00008,headwidth=7.5,headlength=7.5)
        fig.colorbar(im)
        
        plt.savefig(outDirName + "/concentration_imshow_%012d.png"%(t), dpi=1200, format='png')
        plt.close(fig)

