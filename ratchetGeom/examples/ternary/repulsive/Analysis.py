import os, math, re, sys
import struct
import numpy as np
from matplotlib import pyplot as plt
import matplotlib


datadir = "data/" # Directory with simulation output data.

HeaderFile = open(datadir+"Header.mat", 'rb') # Open file which contains simulation parameters.

LX=struct.unpack('=i', HeaderFile.read(4))[0] # Size in X direction.
LY=struct.unpack('=i', HeaderFile.read(4))[0] # Size in Y direction.
LZ=struct.unpack('=i', HeaderFile.read(4))[0] # Size in Z direction.

ndim=struct.unpack('=i', HeaderFile.read(4))[0] # Number of directions (1D, 2D or 3D).

tstart = 0 # Time to start reading from.
tend = struct.unpack('=i', HeaderFile.read(4))[0] # Time to stop reading from.
tinc = struct.unpack('=i', HeaderFile.read(4))[0] # Saving interval.

plot = True # Plot concentrations?

surfacetension01=struct.unpack('=d', HeaderFile.read(8))[0] # Surface tension between droplet and lower fluid.
surfacetension02=struct.unpack('=d', HeaderFile.read(8))[0] # Surface tension between droplet and upper fluid.
surfacetension12=struct.unpack('=d', HeaderFile.read(8))[0] # Surface tension between upper and lower fluid.

numcomp = 3 # 3 immiscible fluid components.

outDirName = "figures" # Directory to save figures.
os.system("mkdir -p %s"%outDirName)

vmag = np.array([])

for t in range(tstart,tend+1,tinc):

    C = np.zeros((LX,LY,2)) # Array to store concentration fluid components 0 and 1. Component 3 is calculated as 1 - conceltration0 - concentration1.
    v = np.zeros((LX,LY,2))
    for iC in range(2):
        file_name = str(datadir)+"OrderParameter"+str(iC)+"_t"+str(t)+".mat"
        File = open(file_name, 'rb') # Open file given by file name above.
        dat=File.read()
        File.close()
        C[:,:,iC] = np.frombuffer(dat, '=d').reshape(LX, LY) # Fill array with fluid concentration of component "iC".

    file_name = str(datadir)+"Velocity_t"+str(t)+".mat"
    File = open(file_name, 'rb') # Open file given by file name above.
    dat=File.read()
    File.close()
    v[:,:,:] = np.frombuffer(dat, '=d').reshape(LX, LY,2)

    print("MAX VELOCITY", np.max(np.sqrt(v[:,:,0]**2 + v[:,:,1]**2))) # Print maximum velocity.
    vmag=np.append(vmag,np.max(np.sqrt(v[:,:,0]**2 + v[:,:,1]**2))) # Save maximum velocity.
    i2 = np.argmax(np.logical_and((C[0,5:LY,1] < 0.9),(C[0,5:LY,1] >= 0.5))) # Find point just above interface of components 1 and 2
    i1 = i2 - 1  # Find point just below interface of components 1 and 2
    c1_1 = C[0,5+i1,1] # Concentration of fluid 1 at i1
    c1_2 = C[0,5+i2,1] # Concentration of fluid 1 at i2
    h = 5+i1 + (c1_1 - 0.5) / (c1_1 - c1_2) # Interpolate to find interface

    if plot==True:

        fluid_visualisation = C[:,:,0] + 0.5*C[:,:,1] # Save weighted sum of fluid concentrations to visualise
        fig,ax=plt.subplots(1,1,figsize=(6,6)) # Make figure
        output = "%s/component_plot_%012d.png"%(outDirName,t)

        rgbv = np.zeros((LY,LX))
        rgbv[:,:] = np.flip(fluid_visualisation).transpose() # Fill array to plot

        ax.imshow(rgbv,interpolation='nearest',origin='lower') # Plot concentration
        stepx=1
        stepy=1
        print(LX,LY)
        X,Y=np.meshgrid(np.linspace(0,LX-1,int((LX)/stepx)),np.linspace(0,LY-1,int((LY)/stepy)))
        ax.quiver(X.T,Y.T,np.flip(-v[0:LX:stepx,0:LY:stepy,0]),np.flip(v[0:LX:stepx,0:LY:stepy,1]),width=0.0008,headwidth=7.5,headlength=7.5)
        
        plt.savefig(output, dpi=400, format='png')  # Save figure
        plt.close(fig)

plt.figure()
plt.plot(vmag)
plt.savefig("vmag.png")

fit_interface.measure_angle2d(C[:,LY//2:,0],None,[0,[0,1]],True,[LX//2,0,20]) # Plot neumann angle fit