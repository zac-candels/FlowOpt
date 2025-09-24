import os, math, re, sys
import struct
import numpy as np
from matplotlib import pyplot as plt
import matplotlib

#sys.path.append("../../../pyutils/")
#import NewLBMCode.TernaryLBM.examples.ternary.lens.fit_interface as fit_interface # Script to fit Neumann angle.
import fit_interface

datadir = "data/" # Directory with simulation output data.

HeaderFile = open(datadir+"Header.mat", 'rb') # Open file which contains simulation parameters.

LX=struct.unpack('=i', HeaderFile.read(4))[0] # Size in X direction.
LY=struct.unpack('=i', HeaderFile.read(4))[0] # Size in Y direction.
LZ=struct.unpack('=i', HeaderFile.read(4))[0] # Size in Z direction.

ndim=struct.unpack('=i', HeaderFile.read(4))[0] # Number of directions (1D, 2D or 3D).

tstart = 0 # Time to start reading from.
tend = struct.unpack('=i', HeaderFile.read(4))[0] # Time to stop reading from.
tinc = struct.unpack('=i', HeaderFile.read(4))[0] # Saving interval.
print(LX,LY,LZ,ndim,tend,tinc)

plot = True # Plot concentrations?

numcomp = 3 # 3 immiscible fluid components.

outDirName = "figures" # Directory to save figures.
os.system("mkdir -p %s"%outDirName)

for t in range(tstart,tend+1,tinc):

    C = np.zeros((LX,LY,2)) # Array to store concentration fluid components 0 and 1. Component 3 is calculated as 1 - conceltration0 - concentration1.
    for iC in range(2):
        file_name = str(datadir)+"OrderParameter"+str(iC)+"_t"+str(t)+".mat"
        File = open(file_name, 'rb') # Open file given by file name above.
        dat=File.read()
        File.close()
        C[:,:,iC] = np.frombuffer(dat, '=d').reshape(LX, LY) # Fill array with fluid concentration of component "iC".

    CNew=C[:,LY//2:,0]
    C06=np.where(CNew>0.5)
    point1y=np.amax(C06[1])
    point1x=LX//2
    point2x=np.amin(C06[0])
    point2y=0
    H=point1y
    W=point1x-point2x
    print("Angle12= ",2*2*np.arctan(H/W)*180/np.pi)

    if plot==True:

        fluid_visualisation = C[:,:,0] + 0.5*C[:,:,1] # Save weighted sum of fluid concentrations to visualise
        fig,ax=plt.subplots(1,1,figsize=(6,6)) # Make figure
        output = "%s/component_plot_%012d.png"%(outDirName,t)

        rgbv = np.zeros((LY,LX))
        rgbv[:,:] = np.flip(fluid_visualisation).transpose() # Fill array to plot

        ax.imshow(rgbv,interpolation='nearest') # Plot concentration
        plt.savefig(output, dpi=400, format='png')  # Save figure
        plt.close(fig)