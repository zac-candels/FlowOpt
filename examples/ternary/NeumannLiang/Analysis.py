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

plot = True # Plot concentrations?

surfacetension01=struct.unpack('=d', HeaderFile.read(8))[0] # Surface tension between droplet and lower fluid.
surfacetension02=struct.unpack('=d', HeaderFile.read(8))[0] # Surface tension between droplet and upper fluid.
surfacetension12=struct.unpack('=d', HeaderFile.read(8))[0] # Surface tension between upper and lower fluid.

AnalyticalNeumannTop = 180*np.arccos((surfacetension01**2+surfacetension12**2-surfacetension02**2) # Analytical angle above the halfway point.
                                     /(2*surfacetension01*surfacetension12))/np.pi
AnalyticalNeumannBottom = 180*np.arccos((surfacetension02**2+surfacetension12**2-surfacetension01**2) # Analytical angle below the halfway point.
                                        /(2*surfacetension02*surfacetension12))/np.pi

AnalyticalNeumann12 = AnalyticalNeumannTop+AnalyticalNeumannBottom # Analytical Neumann angle between phase 1 and 2.
AnalyticalNeumann01 = 180-AnalyticalNeumannTop # Analytical Neumann angle between phase 0 and 2.
AnalyticalNeumann02 = 180-AnalyticalNeumannBottom # Analytical Neumann angle between phase 0 and 2.

print("Analytical Neumann Angle between phase 0 and 1: ",str(round(AnalyticalNeumann01,2))+"°")
print("Analytical Neumann Angle between phase 0 and 2: ",str(round(AnalyticalNeumann02,2))+"°")
print("Analytical Neumann Angle between phase 1 and 2: ",str(round(AnalyticalNeumann12,2))+"°")

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

    i2 = np.argmax(np.logical_and((C[0,5:LY,1] < 0.9),(C[0,5:LY,1] >= 0.5))) # Find point just above interface of components 1 and 2
    i1 = i2 - 1  # Find point just below interface of components 1 and 2
    c1_1 = C[0,5+i1,1] # Concentration of fluid 1 at i1
    c1_2 = C[0,5+i2,1] # Concentration of fluid 1 at i2
    h = 5+i1 + (c1_1 - 0.5) / (c1_1 - c1_2) # Interpolate to find interface

    TopAngle=fit_interface.measure_angle2d(C[:,LY//2:,0],None,[h-LY//2,[0,1]],False,[LX//2,0,20]) # Measured angle above the halfway point.
    BottomAngle=fit_interface.measure_angle2d(C[:,:LY//2,0],None,[h,[0,1]],False,[LX//2,LY//2,20]) # Measured angle below the halfway point.

    MeasuredNeumann01=TopAngle[0] # Measured angle between phases 0 and 1
    MeasuredNeumann02=BottomAngle[0] # Measured angle between phases 0 and 2
    MeasuredNeumann12=360-BottomAngle[0]-TopAngle[0] # Measured angle between phases 1 and 2
    
    print("t="+str(t))
    print("Measured angle between phase 0 and 1: ",str(round(MeasuredNeumann01,2))+"°")
    print("Difference to analytical solution: ",str(round(MeasuredNeumann01-AnalyticalNeumann01,2))+"°")
    print("Measured angle between phase 0 and 2: ",str(round(MeasuredNeumann02,2))+"°")
    print("Difference to analytical solution: ",str(round(MeasuredNeumann02-AnalyticalNeumann02,2))+"°")
    print("Measured angle between phase 1 and 2: ",str(round(MeasuredNeumann12,2))+"°")
    print("Difference to analytical solution: ",str(round(MeasuredNeumann12-AnalyticalNeumann12,2))+"°")

    if plot==True:

        fluid_visualisation = C[:,:,0] + 0.5*C[:,:,1] # Save weighted sum of fluid concentrations to visualise
        fig,ax=plt.subplots(1,1,figsize=(6,6)) # Make figure
        output = "%s/component_plot_%012d.png"%(outDirName,t)

        rgbv = np.zeros((LY,LX))
        rgbv[:,:] = np.flip(fluid_visualisation).transpose() # Fill array to plot

        ax.imshow(rgbv,interpolation='nearest',origin='lower') # Plot concentration
        plt.savefig(output, dpi=400, format='png')  # Save figure
        plt.close(fig)

fit_interface.measure_angle2d(C[:,LY//2:,0],None,[0,[0,1]],True,[LX//2,0,20]) # Plot neumann angle fit