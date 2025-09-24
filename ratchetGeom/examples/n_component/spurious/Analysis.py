import os, math, re, sys
import struct
import numpy as np
from matplotlib import pyplot as plt
from scipy import optimize
import matplotlib

datadir = "data/s12_0.02464/s02_0.015/s01_0.015/tau3_0.505/tau1_0.505/evaporationrate_0.0001/timesteps_450000/lz_101/ly_101/lx_300/"
#datadir = "data/s12_0.0004928/s02_0.0003/s01_0.0003/tau3_1.0/tau1_1.0/evaporationrate_0.0002/timesteps_250000/lz_101/ly_101/lx_300/"
#datadir = "data/volmore/s12_0.02464/s02_0.015/s01_0.015/tau3_0.505/tau1_0.505/evaporationrate_0.0001/timesteps_450000/lz_101/ly_101/lx_300/"
datadir = "data/volmore/s12_0.0004928/s02_0.0003/s01_0.0003/tau3_1.0/tau1_1.0/evaporationrate_0.0001/timesteps_250000/lz_101/ly_101/lx_300/"
datadir = "data/volmore/s12_0.008213/s02_0.005/s01_0.005/tau3_0.55/tau1_0.55/evaporationrate_0.001/timesteps_100000/lz_101/ly_101/lx_300/"
datadir = "data/volmore/s12_0.0005598/s02_0.0003/s01_0.0003/tau3_1.0/tau1_1.0/evaporationrate_1e-05/timesteps_400000/lz_101/ly_101/lx_300/"
datadir = "data/test15/s12_0.0055/"
#datadir = "data/volmore/s12_0.02464/s02_0.015/s01_0.015/tau3_0.505/tau1_0.505/evaporationrate_0.0001/timesteps_450000/lz_101/ly_101/lx_300/"
#datadir = "data/volmore/s12_0.02464/s02_0.015/s01_0.015/tau3_0.505/tau1_0.505/evaporationrate_0.0001/timesteps_450000/lz_101/ly_101/lx_300/"
#datadir = "data/volmore/s12_0.008213/s02_0.005/s01_0.005/tau3_0.55/tau1_0.55/evaporationrate_0.0001/timesteps_450000/lz_101/ly_101/lx_300/"
#datadir = "data/lx_100/ly_100/radius_20/"
#datadir = "data/s12_0.02464/s02_0.015/s01_0.015/tau3_0.505/tau1_0.505/evaporationrate_0.0001/timesteps_450000/lz_101/ly_101/lx_300/"
HeaderFile = open(datadir+"Header.mat", 'rb')

LX=struct.unpack('=i', HeaderFile.read(4))[0]

LY=struct.unpack('=i', HeaderFile.read(4))[0]

LZ=struct.unpack('=i', HeaderFile.read(4))[0]

ndim=struct.unpack('=i', HeaderFile.read(4))[0]

t_zero = 0
tstart = 000

tend = 7450000
tend2 = 0
struct.unpack('=i', HeaderFile.read(4))[0]
tinc = 25000
struct.unpack('=i', HeaderFile.read(4))[0]

plot = False
slicepos=0

sliceaxis=1
if LY==1:
    sliceaxis=1
elif LX==1:
    sliceaxis=0

height = np.array([])
vol=np.array([])
vol2=np.array([[]])
fluxsum=np.array([])
fluxsum2=np.array([])
densityG = 1
vapourFrac = 0.5
diffusion = 0.008
densityL = 1
H = 98

def mlfit(x):
    return (1-(x/np.amax(x))**2)**(1/6-0.5)

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

cangs=np.empty(int((tend-tstart)/tinc)+1)
vols = np.array([])
plot3d = False

print(tend,flush=True)
outDirName = "figures"
os.system("mkdir -p %s"%outDirName)
v = np.zeros((LX,LY,LZ,ndim))
flux=[0.001]
va=np.array([])
vE=np.array([])
vE2=np.array([])
for i in flux:
    #datadir = "data/s12_0.02464/s02_0.015/s01_0.015/tau3_0.505/tau1_0.505/evaporationrate_"+str(i)+"/timesteps_450000/lz_101/ly_101/lx_300/"
    #datadir = "data/s12_0.0004928/s02_0.0003/s01_0.0003/tau3_1.0/tau1_1.0/evaporationrate_"+str(i)+"/timesteps_250000/lz_101/ly_101/lx_300/"
    #datadir = "data/volmore/s12_0.02464/s02_0.015/s01_0.015/tau3_0.505/tau1_0.505/evaporationrate_"+str(i)+"/timesteps_450000/lz_101/ly_101/lx_300/"
    #datadir = "data/volmore/s12_0.0004928/s02_0.0003/s01_0.0003/tau3_1.0/tau1_1.0/evaporationrate_"+str(i)+"/timesteps_250000/lz_101/ly_101/lx_300/"
    datadir = "data/volmore/s12_0.008213/s02_0.005/s01_0.005/tau3_0.55/tau1_0.55/evaporationrate_0.001/timesteps_100000/lz_101/ly_101/lx_300/"
    datadir = "data/volmore/s12_0.0005598/s02_0.0003/s01_0.0003/tau3_1.0/tau1_1.0/evaporationrate_1e-05/timesteps_400000/lz_101/ly_101/lx_300/"
    datadir = "data/test15/s12_0.0055/"
    #datadir = "data/volmore/s12_0.02464/s02_0.015/s01_0.015/tau3_0.505/tau1_0.505/evaporationrate_0.0001/timesteps_450000/lz_101/ly_101/lx_300/"
    #datadir = "data/volmore/s12_0.02464/s02_0.015/s01_0.015/tau3_0.505/tau1_0.505/evaporationrate_0.0001/timesteps_450000/lz_101/ly_101/lx_300/"
    #datadir = "data/s12_0.02464/s02_0.015/s01_0.015/tau3_0.505/tau1_0.505/evaporationrate_0.0001/timesteps_450000/lz_101/ly_101/lx_300/"
    #datadir = "data/volmore/s12_0.008213/s02_0.005/s01_0.005/tau3_0.55/tau1_0.55/evaporationrate_0.0001/timesteps_450000/lz_101/ly_101/lx_300/"
    for t in range(tstart,tend+1,tinc):
        print("t=%s"%t,flush=True)
        t_file =t+t_zero

        #file_name = "data/"+"OrderParameter_t%li.mat"%t_file
        #file_name = datadir+"Pressure_t%li.mat"%t_file
        #file_name = datadir+"Density_t%li.mat"%t_file
        #file_name = "data/"+"Humidity_t%li.mat"%t_file
        #file_name = datadir+"OrderParameter_t%li.mat"%t_file
        #file_name = datadir+"ChemicalPotential_t%li.mat"%t_file
        file_name = datadir+"OrderParameter0_t%li.mat"%t_file
        #file_name = "data/"+"BoundaryLabels_t%li.mat"%t_file

        File = open(file_name, 'rb')

        file_name = datadir+"Density_t%li.mat"%t_file
        #file_name = "data/"+"BoundaryLabels_t%li.mat"%t_file

        FileDensity = open(file_name, 'rb')

        #file_name = "data/"+"OrderParameter_t%li.mat"%t_file
        #file_name = datadir+"Pressure_t%li.mat"%t_file
        #file_name = datadir+"Density_t%li.mat"%t_file
        #file_name = "data/"+"Humidity_t%li.mat"%t_file
        #file_name = datadir+"OrderParameter_t%li.mat"%t_file
        #file_name = datadir+"ChemicalPotential_t%li.mat"%t_file
        file_name = datadir+"OrderParameter1_t%li.mat"%t_file
        #file_name = "data/"+"BoundaryLabels_t%li.mat"%t_file

        FileOP1 = open(file_name, 'rb')

        #file_name = "data/"+"OrderParameter_t%li.mat"%t_file
        #file_name = datadir+"Pressure_t%li.mat"%t_file
        #file_name = datadir+"Density_t%li.mat"%t_file
        #file_name = "data/"+"Humidity_t%li.mat"%t_file
        #file_name = datadir+"OrderParameter_t%li.mat"%t_file
        #file_name = datadir+"ChemicalPotential_t%li.mat"%t_file
        file_name = datadir+"OrderParameter2_t%li.mat"%t_file
        #file_name = "data/"+"BoundaryLabels_t%li.mat"%t_file

        FileOP2 = open(file_name, 'rb')

        file_name = datadir+"BoundaryLabels_t%li.mat"%t_file
        #file_name = "data/"+"Pressure_t%li.mat"%t_file
        #file_name = "data/"+"Density_t%li.mat"%t_file
        #file_name = "data/"+"Humidity_t%li.mat"%t_file
        #file_name = "data/"+"LaplacianOrderParameter_t%li.mat"%t_file
        #file_name = "data/"+"ChemicalPotential_t%li.mat"%t_file
        #file_name = "data/"+"BoundaryLabels_t%li.mat"%t_file

        FileSolid = open(file_name, 'rb')
    
        
        file_name = datadir+"Velocity_t%li.mat"%t_file

        File2 = open(file_name, 'rb')

        #file_name = datadir+"MixedGradientPressure_t%li.mat"%t_file

        #FileGP = open(file_name, 'rb')

        #file_name = "data/"+"OrderParameter_t%li.mat"%t_file
        #file_name = datadir+"Pressure_t%li.mat"%t_file
        #file_name = datadir+"Density_t%li.mat"%t_file
        #file_name = datadir+"ChemicalPotential_t%li.mat"%t_file
        #file_name = datadir+"OrderParameter_t%li.mat"%t_file
        file_name = datadir+"ChemicalPotential0_t%li.mat"%t_file
        #file_name = datadir+"LaplacianOrderParameter_t%li.mat"%t_file
        #file_name = "data/"+"BoundaryLabels_t%li.mat"%t_file

        File3 = open(file_name, 'rb')

        file_name = datadir+"ChemicalPotential1_t%li.mat"%t_file
        #file_name = datadir+"LaplacianOrderParameter_t%li.mat"%t_file
        #file_name = "data/"+"BoundaryLabels_t%li.mat"%t_file

        File4 = open(file_name, 'rb')

        file_name = datadir+"ChemicalPotential2_t%li.mat"%t_file
        #file_name = datadir+"LaplacianOrderParameter_t%li.mat"%t_file
        #file_name = "data/"+"BoundaryLabels_t%li.mat"%t_file

        File5 = open(file_name, 'rb')

        file_name = datadir+"ChemicalPotential3_t%li.mat"%t_file
        #file_name = datadir+"LaplacianOrderParameter_t%li.mat"%t_file
        #file_name = "data/"+"BoundaryLabels_t%li.mat"%t_file

        File6 = open(file_name, 'rb')

        #file_name = "data/"+"Density_t%li.mat"%t_file

        #File3 = open(file_name, 'rb')
        print(file_name)

        def coord_k(k, LY, LZ):
            """From a k value, determines its xk, yk, and zk."""    

            xk = math.floor(k/(LY*LZ))
            yk = math.floor((k - xk*LZ*LY)/LZ)
            zk = k - xk*LZ*LY - yk*LZ
            return xk, yk, zk

        NLatt=LX*LY*LZ

        rho = np.zeros((LX,LY,LZ))
        v = np.zeros((LX,LY,LZ,ndim))
        gradp = np.zeros((LX,LY,LZ,ndim))
        solid = np.zeros((LX,LY,LZ))

        dat=File.read()
        rho = np.ndarray((LX,LY,LZ),'=d',dat,0,(8*LY*LZ,8*LZ,8))
        dat=FileOP1.read()
        rho2 = np.ndarray((LX,LY,LZ),'=d',dat,0,(8*LY*LZ,8*LZ,8))
        dat=FileSolid.read()
        solid = np.ndarray((LX,LY,LZ),'=i',dat,0,(4*LY*LZ,4*LZ,4))
        dat=FileOP2.read()
        rho3 = np.ndarray((LX,LY,LZ),'=d',dat,0,(8*LY*LZ,8*LZ,8))

        dat=File3.read()
        cp = np.array(np.ndarray((LX,LY,LZ),'=d',dat,0,(8*LY*LZ,8*LZ,8)))
        cp[np.where(solid!=0)]=np.average(cp[np.where(solid==0)])

        dat=File4.read()
        cp2 = np.array(np.ndarray((LX,LY,LZ),'=d',dat,0,(8*LY*LZ,8*LZ,8)))
        cp2[np.where(solid!=0)]=np.average(cp2[np.where(solid==0)])

        dat=File5.read()
        cp3 = np.array(np.ndarray((LX,LY,LZ),'=d',dat,0,(8*LY*LZ,8*LZ,8)))
        cp3[np.where(solid!=0)]=np.average(cp3[np.where(solid==0)])
        
        dat=File6.read()
        cp4 = np.array(np.ndarray((LX,LY,LZ),'=d',dat,0,(8*LY*LZ,8*LZ,8)))
        cp4[np.where(solid!=0)]=np.average(cp4[np.where(solid==0)])
        

        dat=FileDensity.read()
        density = np.ndarray((LX,LY,LZ),'=d',dat,0,(8*LY*LZ,8*LZ,8))

        #liquid = np.array(-lap[:,3,:,1]-lap[:,3,:,0])
        #liquid = rho[:,5,:,0]
        #####liquid = rho[:,5,:,0]
        #liquid = humid[:,:,0]
        #liquid = cp[:,5,:,0]+3*(0.003-0.002464)*lap[:,5,:,0]#1-rho[:,3,:,1]-rho[:,3,:,0]
        #liquid = np.array([np.gradient(liquid,axis=0),np.gradient(liquid,axis=1)])
        #liquid = np.array([np.gradient(liquid[0],axis=0),np.gradient(liquid[1],axis=1)])
        #liquid = liquid[0]+liquid[1]

        dat=File2.read()
        v = np.ndarray((LX,LY,LZ,ndim),'=d',dat,0,(ndim*8*LY*LZ,ndim*8*LZ,ndim*8,8))

        #dat=FileGP.read()
        #gradp = np.ndarray((LX,LY,LZ,ndim),'=d',dat,0,(ndim*8*LY*LZ,ndim*8*LZ,ndim*8,8))
        if (plot3d):
            if(t==tstart):
                liquid = np.array(rho[:,:,:,1])
                #liquid = np.concatenate((np.flip(np.array(liquid[:,:LY//2-1,:]),axis=1),np.array(liquid),np.flip(np.array(liquid[:,LY//2-1:,:]),axis=1)),axis=1)
                #liquid = np.concatenate((np.flip(liquid[:,:,0:LZ//2-1],axis=2),liquid,np.flip(liquid[:,:,LZ//2-1:],axis=2)),axis=2)
                liquid[:1,np.where(liquid[:1,:,:]>0.5)[1],np.where(liquid[:1,:,:]>0.5)[2]]=0

                liquid[np.where(liquid[:,:1,:]>0.5)[0],:1,np.where(liquid[:,:1,:]>0.5)[2]]=0
                liquid[np.where(liquid[:,:,:1]>0.5)[0],np.where(liquid[:,:,:1]>0.5)[1],:1]=0
                liquid[LX-2:,np.where(liquid[LX-2:,:,:]>0.5)[1],np.where(liquid[LX-2:,:,:]>0.5)[2]]=0

                liquid[np.where(liquid[:,LY-2:,:]>0.5)[0],LY-2:,np.where(liquid[:,LY-2:,:]>0.5)[2]]=0
                liquid[np.where(liquid[:,:,LZ-2:]>0.5)[0],np.where(liquid[:,:,LZ-2:]>0.5)[1],LZ-2:]=0

            liquid2 = np.array(rho[:,:,:,0])
            liquid2[:1,np.where(liquid2[:1,:,:]>0.5)[1],np.where(liquid2[:1,:,:]>0.5)[2]]=0

            liquid2[np.where(liquid2[:,:1,:]>0.5)[0],:1,np.where(liquid2[:,:1,:]>0.5)[2]]=0
            liquid2[np.where(liquid2[:,:,:1]>0.5)[0],np.where(liquid2[:,:,:1]>0.5)[1],:1]=0
            liquid2[LX-2:,np.where(liquid2[LX-2:,:,:]>0.5)[1],np.where(liquid2[LX-2:,:,:]>0.5)[2]]=0

            liquid2[np.where(liquid2[:,LY-2:,:]>0.5)[0],LY-2:,np.where(liquid2[:,LY-2:,:]>0.5)[2]]=0
            liquid2[np.where(liquid2[:,:,LZ-2:]>0.5)[0],np.where(liquid2[:,:,LZ-2:]>0.5)[1],LZ-2:]=0

            #liquid2[np.where(liquid2[:,:,:2]>0.333)[0],np.where(liquid2[:,:,:2]>0.333)[1],:2]=0
            print(np.sum(liquid2))
            vols = np.append(vols, np.sum(liquid2))
        else:
            #liquid = np.array(152000./907./50./(1./3.)/50.*cp2[:,:,LZ//2]-66000./907./50./(1./3.)/50.*cp[:,:,LZ//2]-8000./907./50./(1./3.)/50.*cp3[:,:,LZ//2]-78000./907./50./(1./3.)/50.*cp4[:,:,LZ//2])#-cp4[:,:,LZ//2])
            #liquid = np.array(148000./907./50./(1./3.)/50.*cp[:,:,LZ//2]-66000./907./50./(1./3.)/50.*cp2[:,:,LZ//2]-92000./907./50./(1./3.)/50.*cp3[:,:,LZ//2]+10000./907./50./(1./3.)/50.*cp4[:,:,LZ//2])
            #liquid = np.array(167500./907./50./(1./3.)/50.*cp3[:,:,LZ//2]-92000./907./50./(1./3.)/50.*cp[:,:,LZ//2]-8000./907./50./(1./3.)/50.*cp2[:,:,LZ//2]-67500./907./50./(1./3.)/50.*cp4[:,:,LZ//2])
            #print(np.amax(liquid))
            #liquid = np.array(cp[:,:,LZ//2])
            #liquid=np.array(1-rho3-rho2-rho)
            #liquid = np.array(rho3[:,:,LZ//2])
            #liquid=np.array(rho[:,:,LZ//2])*np.array(rho2[:,:,LZ//2])*(1-np.array(rho[:,:,LZ//2])-np.array(rho2[:,:,LZ//2]))
            liquid = np.array(rho[:,:,LZ//2])+ 0.666*np.array(rho2[:,:,LZ//2])+0.333*np.array(rho3[:,:,LZ//2])#*np.array(gh[:,LY//2,:,0])**2+np.array(gh[:,LY//2,:,1])**2+np.array(gh[:,LY//2,:,2])**2)
            #
            #liquid = np.sqrt(np.array(gh[:,LY//2,:,0])**2+np.array(gh[:,LY//2,:,1])**2+np.array(gh[:,LY//2,:,2])**2)
            #liquid = np.array(solid[:,LY//2,:])
            #liquid = (np.sign(nTernaryLeeExtraPotentialp.array(go[:,LY//2,:,0,0])*(np.array(go[:,LY//2,:,0,0])+np.array(go[:,LY//2,:,1,0]))+np.array(go[:,LY//2,:,0,1])*(np.array(go[:,LY//2,:,0,1])+np.array(go[:,LY//2,:,1,1]))+np.array(go[:,LY//2,:,0,2])*(np.array(go[:,LY//2,:,0,2])+np.array(go[:,LY//2,:,1,2]))))*np.sqrt(np.abs(np.array(go[:,LY//2,:,0,0])*(np.array(go[:,LY//2,:,0,0])+np.array(go[:,LY//2,:,1,0]))+np.array(go[:,LY//2,:,0,1])*(np.array(go[:,LY//2,:,0,1])+np.array(go[:,LY//2,:,1,1]))+np.array(go[:,LY//2,:,0,2])*(np.array(go[:,LY//2,:,0,2])+np.array(go[:,LY//2,:,1,2]))))
            #liquid = np.sign(np.array(go[:,LY//2,:,0,2]))*np.sqrt(np.abs(go[:,LY//2,:,0,2]*(np.array(go[:,LY//2,:,0,2])+np.array(go[:,LY//2,:,1,2]))))*gh[:,LY//2,:,2]
        thetafit=0#circleFit(liquid,LX,LY,12,LY-1,0,LX-1,0)
        print(thetafit)
        #dat=File3.read()
        #gh = np.ndarray((LX,LY,ndim),'=d',dat,0,(ndim*8*LY,ndim*8,8))
        print("LIQUID VOLUME",np.sum(rho[:,:,LZ//2]))
        File.close()
        FileSolid.close()
        File2.close()
        File3.close()
        print(np.array(rho2[LX//2,LY//2,LZ//2]))
        #####fig,ax=plt.subplots(1,1,figsize=(6,6))
        #fig = plt.figure()
        #ax = plt.axes(projection='3d')
        output = "%s/component_plot_%012d.png"%(outDirName,t)
        #rho3=2*0.01*(rho2-0.2)*(rho2-1)*(2*rho2-0.2-1)-0.0128*rho4
        if (not plot3d):
            rgbv = np.zeros((LY,LX))
            rgbv[:,:] = np.flip(liquid).transpose()
        #rgbv[:,:] = np.flip(liquid[:,:]).transpose()
        #print(np.sum(liquid2))
        #rgbv[:,:] = np.flip(rho.take(indices=slicepos,axis=sliceaxis)).transpose()
        #rgbv[:,:,1] = np.flip(rho.take(indices=slicepos,axis=sliceaxis)).transpose()
        #rgbv[:,:,2] = np.flip(rho.take(indices=slicepos,axis=sliceaxis)).transpose()
        
        #im=ax.imshow(np.flip(rho0.take(indices=slicepos,axis=sliceaxis)).transpose(),interpolation='nearest',origin='upper')

        #print(v[228-5,195,0])
        #pts = np.where(1-rho[:,:,:,0]-rho[:,:,:,1]>0.5)
        #print(np.shape(pts))
        #print(pts)
        #X,Y,Z=np.meshgrid(np.linspace(0,101,102), np.linspace(0,300,301), np.linspace(0,101,102))
        #ax.voxels(X,Y,Z,rho[:,:,:,0])
        
        if plot == True:
            fig,ax=plt.subplots(1,1,figsize=(6,6))
            im=ax.imshow(rgbv,interpolation='nearest',origin='upper')
            fig.colorbar(im)
            stepx=1
            stepz=1
            ax.set_aspect('equal')
            X,Z=np.meshgrid(np.linspace(0,LX-1,int((LX)/stepx)),np.linspace(0,LY-1,int((LY)/stepz)))
            ax.quiver(X.T,Z.T,np.flip(np.flip(-v[0:LX:stepx,0:LY:stepz,0,0],0),1),np.flip(v[0:LX:stepx,0:LY:stepz,0,1]),width=0.0008,headwidth=7.5,headlength=7.5)
            #plt.scatter(thetafit[2],80-(thetafit[3]+10))
            plt.savefig(output, dpi=400, format='png')
            plt.close(fig)
        #plt.scatter(61,92,s=0.5,color="r")
        #print("here ",rho[161,5,92,0],rho[161,5,92,1],1-rho[161,5,92,0]-rho[161,5,92,1])
        #print("here2 ",liquid[161,92])
        #print("here3 ",3*(0.003-0.002464)*lap[161,5,92,0])
        #ax.contour(np.flip((1-rho[:,2,:,0]-rho[:,2,:,1])).T, levels=[0.5], colors="k", zorder=1, linewidths=0.75)
        #im=ax.imshow(np.sqrt((gh.take(indices=0,axis=2).take(indices=slicepos,axis=sliceaxis))**2+(gh.take(indices=1,axis=2).take(indices=slicepos,axis=sliceaxis))**2),interpolation='nearest',origin='upper')
        #im=ax.imshow(np.sqrt((go[:,:,1,:].take(indices=0,axis=2))**2+(go[:,:,1,:].take(indices=1,axis=2))**2),interpolation='nearest',origin='upper')
        #im=ax.imshow(np.sqrt((gh.take(indices=0,axis=2).take(indices=slicepos,axis=sliceaxis))**2),interpolation='nearest',origin='upper')
        #im=ax.imshow(np.flip(np.flip(np.sqrt((v.take(indices=0,axis=2))**2+(v.take(indices=1,axis=2))**2).T,1),0),interpolation='nearest',origin='upper')
        #print(np.flip(rho.take(indices=slicepos,axis=sliceaxis)).transpose()[70,70])
        #ax.scatter(70,70)
        stepx=1
        stepz=1
        #####X,Z=np.meshgrid(np.linspace(0,LX-1,int((LX)/stepx)),np.linspace(0,LY-1,int((LY)/stepz)))
        #print(np.sum(np.sqrt((gh.take(indices=0,axis=3).take(indices=slicepos,axis=sliceaxis))**2+(gh.take(indices=1,axis=3).take(indices=slicepos,axis=sliceaxis))**2)))
        #ax.quiver(X.T,Y.T,np.flip(v[:,:,0].take(indices=slicepos,axis=sliceaxis)),np.flip(-v[:,:,1].take(indices=slicepos,axis=sliceaxis)),width=0.001,headwidth=2.5,headlength=1.5)
        #ax.quiver(X.T,Y.T,np.flip(-v[0:-1:stepx,0:-1:stepy,0]),np.flip(v[0:-1:stepx,0:-1:stepy,1]),width=0.0002,headwidth=7.5,headlength=7.5)
        #####ax.quiver(X.T,Z.T,np.flip(np.flip(-v[0:LX:stepx,5,0:LY:stepz,0],0),1),np.flip(v[0:LX:stepx,5,0:LY:stepz,1]),width=0.0008,headwidth=7.5,headlength=7.5)
        #print(np.sum(rho[int(LX/2),:])/(0.002/(100-h)*np.log(1/(1-0.1))))
        print("V ",np.average(np.sqrt(v[:,:,:,0]**2+v[:,:,:,1]**2)))#np.sqrt(v[:,:,0]**2+v[:,:,1]**2)))
        #print("HERE ",np.sum(rho))
        #print("HERE ",np.sum(rho>0.5))
        #print("HERE ",rho[LX//8,LY//2])
        #print("HERE ",rho[LX//2,LY//6])
        #####fig.colorbar(im)
        #ax.scatter(49,49)
        #####plt.savefig(output, dpi=400, format='png')
        #####plt.close(fig)
        #plt.figure()
        #plt.plot(rho[int(LX/2),:])
        #plt.plot(rho[:,4,LZ//2])
        #plt.plot(rho[:,4,LZ//2])
        #plt.savefig("test_%012d.png"%(t), dpi=200, format='png')
        #plt.close()
        if t==tstart:
            vol0=np.sum(rho[1:LX-1,1:LY-1,1:LZ-1])
            vol=np.append(vol,0)
        else:
            vol=np.append(vol,(np.sum(rho[1:LX-1,1:LY-1,1:LZ-1]))-vol0)
        #print(gh[int(LX/2),51,1])
        va = np.append(va,np.amax(np.sqrt(v[:,:,:,0]**2+v[:,:,:,1]**2)))
        vE = np.append(vE,np.sum(density*(v[:,:,:,0]**2+v[:,:,:,1]**2)))
    #np.append(va,np.amax(np.amax(np.sqrt(v[:,:,:,0]**2+v[:,:,:,1]**2))-np.average(v[:,:,:,0])))#np.amax(v[:,:,:,0]))

fig, ax = plt.subplots(1, 1, figsize = (6, 5))
plt.plot(rho[:,LY//2,LZ//2])
plt.plot(rho2[:,LY//2,LZ//2])
plt.plot(rho3[:,LY//2,LZ//2])
plt.plot(1-rho[:,LY//2,LZ//2]-rho2[:,LY//2,LZ//2]-rho3[:,LY//2,LZ//2])
#plt.plot(np.linspace(tstart,tend2,(tend2-tstart)//tinc+1),vE2)
plt.savefig("C.png", dpi=400, format='png')
plt.close(fig)

fig, ax = plt.subplots(1, 1, figsize = (6, 5))
ax.set_yscale('log')
ax.set_xscale('log')
plt.ylim([10**-30,10**-0])
plt.xlim([10**4,3*10**7])
plt.ylabel("Kinetic Energy")
plt.xlabel("Simulation Iterations")
plt.plot(np.linspace(tstart,tend,(tend-tstart)//tinc+1),vE)
#plt.plot(np.linspace(tstart,tend2,(tend2-tstart)//tinc+1),vE2)
plt.savefig("vmax.png", dpi=400, format='png')
plt.close(fig)

print(vol2)
plt.figure()
#plt.plot(np.sqrt(np.diagonal(v[:,:,0])**2+np.diagonal(v[:,:,1])**2))
#plt.plot(v[LX//2,:,0]/3)

plt.plot((vol2[0,:])-(vol2[0,-1]),fluxsum2[0,:]*2,label="0.0001")
#plt.plot(-(vol2[1,:]),fluxsum2[1,:],label="0.0001")
#plt.plot(-(vol2[2,:]),fluxsum2[2,:],label="0.00005")
#plt.plot(-(vol2[3,:]),fluxsum2[3,:],label="0.00002")
plt.xlabel("Volume $μm^3$")
plt.ylabel("Area $μm^2$")
plt.xlim([440000,-20000])
#plt.plot(fluxsum2[0,:],label="0.0002")
#plt.plot(fluxsum2[1,:],label="0.0001")
#plt.plot(fluxsum2[2,:],label="0.00005")
#plt.plot(fluxsum2[3,:],label="0.00002")
#plt.legend()
plt.savefig("test6.png", dpi=200, format='png')

plt.figure()
#plt.plot(np.sqrt(np.diagonal(v[:,:,0])**2+np.diagonal(v[:,:,1])**2))
#plt.plot(v[LX//2,:,0]/3)

plt.plot((vol2[0,:])-(vol2[0,-1]),label="0.0001")
plt.ylabel("Volume $μm^3$")
plt.xlabel("time")
plt.savefig("test7.png", dpi=200, format='png')


plt.figure()
plt.xlabel("Time (lattice Units)")
plt.ylabel("Height (lattice Units)")
# plt.ylim([0,1.1])

t = tinc*np.linspace(0,height.size-1,height.size)

plt.scatter(t,height[0:],label="Data")

def height_func(t, h0, k):

    a = np.array(2*k*t + (H-h0)**2)
    a[a<0] = 0
    return H - np.sqrt(a)

def height_func2(t, h0, k):
    dt = np.diff(t, prepend=[0])
    intGradVapour = np.cumsum(gradVapour * dt)
    return h0 - k * intGradVapour

#popt, _ = curve_fit(height_func, t, height, p0=[70,-1e-3])
#print(popt)
# print(height_func(t, popt[0], popt[1]))
#plt.plot(t, height_func(t, popt[0], popt[1]), 'r-', label="Fit")

k = (densityG*1/(1-vapourFrac)) *vapourFrac * diffusion/densityL

h0 = height[0]
print(h0, k)
#print(popt[1]/k)
height_an = height_func(t, h0, k)
#height_an = height_func2(t, h0, 0.5*k/vapourFrac)
#height_an2 = height_func2(t, h0, k/vapourFrac)
plt.plot(t, height_an, 'k-', label="Theory")
#plt.plot(humidity[0,:])
#plt.plot(masssink[0,:]/np.amax(masssink[0,:]))
#plt.plot(0.5-0.5*np.tanh(0.5*(np.linspace(0,humidity[0,:].size-1,humidity[0,:].size)-h)))
#plt.plot(liquid[0,:])
#plt.plot(np.linspace(height[-1],399,height.size),np.amax(humidity[0,:])*(1-np.linspace(0,1,height.size)))
#plt.scatter(height[-1],0)
#plt.plot(t, height_an2, 'r-', label="Matching mass loss")
#plt.plot([0.1483, 0.16367, 0.16935, 0.16904, 0.16935])
print((height[0]-height[-1])/(height_an[0]-height_an[-1]))
plt.legend()
plt.savefig("massloss.png", dpi=500, format='png')

def parabola(x,a,b,c):
    return a*x**2 + b*x + c

def analysePoiseuille():

    parabfit=np.empty([3])

    fit_params, pcov = optimize.curve_fit(parabola, np.linspace(7,LY-8,LY-14), v[math.trunc(LX/2),7:LY-7,0]/3)
    parabfit[:]=fit_params
    print(parabfit)
    plt.figure()
    print(np.amax(v[:,:,0]),0.0000001/(2*(1/3*0.5))*((LY-12.5)/2)**2*3)
    plt.plot(np.linspace(0,LY-1,LY),parabola(np.linspace(0,LY-1,LY),*parabfit[:]))
    plt.plot(v[LX//2,:,0]/3)
    plt.xlabel("Length Across Z Direction (lattice units)")
    plt.ylabel("Parabola fit")
    plt.tight_layout()
    plt.savefig("poiseuille.png", dpi=500, format='png')
    plt.close(fig)

analysePoiseuille()
