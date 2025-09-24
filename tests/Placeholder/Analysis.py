import os, math, re, sys
import struct
import numpy as np
from matplotlib import pyplot as plt
from scipy import optimize
import matplotlib

#datadir = "data/inflow3/diffusivity_0.0008/inflowmomentum_-0.0002/lx_228/ly_200/postwidth_162/offsety_-17/theta_150/"
#datadir = "data/inflowsimple7/surfacetension_0.0001/theta_90/offsety_-36.0/postwidth_486.0/ly_600.0/lx_684.0/inflowmomentum_-0.0006/diffusivity_0.024/"
#datadir = "data/inflowsimple2/inflowmomentum_-0.0002/lx_228/ly_200/postwidth_162/offsety_-17/theta_120/"
#datadir = "data/inflow3/inflowmomentum_-0.002/lx_228/ly_300/postwidth_162/offsety_-17/theta_30/"
#datadir = "data/inflowsimple7/surfacetension_0.0001/theta_150/offsety_-36/postwidth_270/ly_399/lx_468/inflowmomentum_-0.0005/diffusivity_0.012/"
datadir = "data/inflowsimple7/surfacetension_0.0001/theta_30/offsety_-12/postwidth_90/ly_133/lx_156/inflowmomentum_-0.0006/diffusivity_0.004/"

HeaderFile = open(datadir+"Header.mat", 'rb')

LX=struct.unpack('=i', HeaderFile.read(4))[0]

LY=struct.unpack('=i', HeaderFile.read(4))[0]

LZ=struct.unpack('=i', HeaderFile.read(4))[0]

ndim=struct.unpack('=i', HeaderFile.read(4))[0]

t_zero = 0
tstart = 0

tend = struct.unpack('=i', HeaderFile.read(4))[0]
tinc = struct.unpack('=i', HeaderFile.read(4))[0]

slicepos=0

sliceaxis=1
if LY==1:
    sliceaxis=1
elif LX==1:
    sliceaxis=0

height = np.array([])

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
def fitfunc(x, a, b, c):
    return a+b*x**c
def derivfunc(x, b, c):
    return c*b*x**(c-1)
print(tend)
outDirName = "figures"
os.system("mkdir -p %s"%outDirName)
v = np.zeros((LX,LY,LZ,ndim))
#thetas = [150,120,90,60,30]
thetas = [30]

dml = {}
dheightavg = {}
dvol = {}
#a=["a","aa","aaa","b","bb","bbb","c","cc","ccc","d","dd","ddd","e","ee","eee","f","ff","fff","g","gg","ggg","h","hh","hhh","i","ii","iii","j","jj","jjj","k","kk","kkk","l","ll","lll","m","mm","mmm","n","nn","nnn","o","oo","ooo","p","pp","ppp","q","qq","qqq","r","rr","rrr","s","ss","sss","t","tt","ttt","u","uu","uuu","v","vv","vvv","w","ww","www","x","xx","xxx","y","yy","yyy","z","zz""zzz"]
for th in thetas:
    heightavg = np.array([])
    ml = np.array([])
    vol = np.array([])
    #datadir = "data/inflowsimple/inflowmomentum_-0.0002/lx_228/ly_200/postwidth_162/offsety_-17/theta_"+str(th)+"/"
    #datadir = "data/inflow3/diffusivity_0.0008/inflowmomentum_-0.0002/lx_228/ly_200/postwidth_162/offsety_-17/theta_"+str(th)+"/"
    #datadir = "data/inflowsimple7/surfacetension_0.0001/theta_"+str(th)+"/offsety_-48/postwidth_360/ly_532/lx_624/inflowmomentum_-0.0006/diffusivity_0.016/"
    datadir = "data/inflowsimple7/surfacetension_0.0001/theta_"+str(th)+"/offsety_-12/postwidth_90/ly_133/lx_156/inflowmomentum_-0.0005/diffusivity_0.004/"
    #datadir = "data/inflowsimple7/surfacetension_0.0001/theta_"+str(th)+"/offsety_-36/postwidth_270/ly_399/lx_468/inflowmomentum_-0.0005/diffusivity_0.012/"
    #datadir = "data/inflowsimple2/inflowmomentum_-0.0002/lx_228/ly_200/postwidth_162/offsety_-17/theta_"+str(th)+"/"data/inflowsimple7/surfacetension_0.0001/theta_30/offsety_-12/postwidth_90/ly_133/lx_156/inflowmomentum_-0.0005/diffusivity_0.004/
    #datadir = "data/inflow3/inflowmomentum_-0.002/lx_228/ly_300/postwidth_162/offsety_-17/theta_"+str(th)+"/"
    for t in range(tstart,tend+1,tinc):
        print("t=%s"%t)
        t_file =t+t_zero

        #file_name = "data/"+"OrderParameter_t%li.mat"%t_file
        #file_name = datadir+"Pressure_t%li.mat"%t_file
        #file_name = datadir+"Density_t%li.mat"%t_file
        #file_name = "data/"+"Humidity_t%li.mat"%t_file
        file_name = datadir+"OrderParameter_t%li.mat"%t_file
        #file_name = datadir+"ChemicalPotential_t%li.mat"%t_file
        #file_name = "data/"+"BoundaryLabels_t%li.mat"%t_file

        File = open(file_name, 'rb')

        #file_name = datadir+"MassSink_t%li.mat"%t_file
        #file_name = datadir+"Pressure_t%li.mat"%t_file
        #file_name = datadir+"Density_t%li.mat"%t_file
        file_name = datadir+"Humidity_t%li.mat"%t_file
        #file_name = "data/"+"LaplacianOrderParameter_t%li.mat"%t_file
        #file_name = "data/"+"ChemicalPotential_t%li.mat"%t_file
        #file_name = "data/"+"BoundaryLabels_t%li.mat"%t_file

        File0 = open(file_name, 'rb')

        file_name = datadir+"MassSink_t%li.mat"%t_file
        #file_name = "data/"+"Pressure_t%li.mat"%t_file
        #file_name = "data/"+"Density_t%li.mat"%t_file
        #file_name = "data/"+"Humidity_t%li.mat"%t_file
        #file_name = "data/"+"LaplacianOrderParameter_t%li.mat"%t_file
        #file_name = "data/"+"ChemicalPotential_t%li.mat"%t_file
        #file_name = "data/"+"BoundaryLabels_t%li.mat"%t_file

        File00 = open(file_name, 'rb')

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

        #file_name = datadir+"GradientHumidity_t%li.mat"%t_file

        #File3 = open(file_name, 'rb')

        file_name = datadir+"OrderParameter_t%li.mat"%t_file

        File4 = open(file_name, 'rb')

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
        rho2 = np.zeros((LX,LY,LZ))
        v = np.zeros((LX,LY,LZ,ndim))
        solid = np.zeros((LX,LY,LZ))

        dat=File.read()
        rho = np.ndarray((LX,LY),'=d',dat,0,(8*LY,8))

        dat=File0.read()
        rho2 = np.ndarray((LX,LY),'=d',dat,0,(8*LY,8))
        #rho = np.ndarray((LX,LY),'=i',dat,0,(4*LY,4))

        dat=File00.read()
        humidity = np.ndarray((LX,LY),'=d',dat,0,(8*LY,8))

        dat=FileSolid.read()
        solid = np.ndarray((LX,LY),'=i',dat,0,(4*LY,4))

        liquid = np.array(rho)
        liquid[np.where(np.logical_or(solid==1,solid==-1))[0],np.where(np.logical_or(solid==1,solid==-1))[1]] = 0.5

        mlnosolid = np.array(humidity)
        mlnosolid[np.where(np.logical_or(solid==1,solid==-1))[0],np.where(np.logical_or(solid==1,solid==-1))[1]] = 0

        dat=File2.read()
        v = np.ndarray((LX,LY,ndim),'=d',dat,0,(ndim*8*LY,ndim*8,8))

        #dat=File3.read()
        #gh = np.ndarray((LX,LY,ndim),'=d',dat,0,(ndim*8*LY,ndim*8,8))

        dat=File4.read()
        c = np.ndarray((LX,LY),'=d',dat,0,(8*LY,8))

        File.close()
        File2.close()
        #File3.close()
        File4.close()
        
        

        #output = "%s/"%(outDirName)+str(a[t//25000])+"component_plot_%012d.png"%(t)
        output = "%s/component_plot_%012d.png"%(outDirName,t)
        #rho3=2*0.01*(rho2-0.2)*(rho2-1)*(2*rho2-0.2-1)-0.0128*rho4
        rgbv = np.zeros((LY,LX))
        rgbv[:,:] = np.flip(mlnosolid).transpose()
        #rgbv[:,:] = np.flip(rho.take(indices=slicepos,axis=sliceaxis)).transpose()
        #rgbv[:,:,1] = np.flip(rho.take(indices=slicepos,axis=sliceaxis)).transpose()
        #rgbv[:,:,2] = np.flip(rho.take(indices=slicepos,axis=sliceaxis)).transpose()
        
        #im=ax.imshow(np.flip(rho0.take(indices=slicepos,axis=sliceaxis)).transpose(),interpolation='nearest',origin='upper')
        i2 = np.argmax(np.logical_and((c[0,:] <= 0.5),(c[0,:] > 0.1)))

        #print("sum",np.sum([0.5*(liquid[0,i+1]-liquid[0,i-1]) for i in range(i2,H-2)]))
        i1 = i2 - 1
        c1_1 = c[0,i1]
        c1_2 = c[0,i2]
        h = i1 + (c1_1 - 0.5) / (c1_1 - c1_2) + 1
        print(h)
        #print(v[228-5,195,0])

        
        #im=ax.imshow(np.sqrt((gh.take(indices=0,axis=2).take(indices=slicepos,axis=sliceaxis))**2+(gh.take(indices=1,axis=2).take(indices=slicepos,axis=sliceaxis))**2),interpolation='nearest',origin='upper')
        #im=ax.imshow(np.sqrt((gh.take(indices=0,axis=2))**2+(gh.take(indices=1,axis=2))**2),interpolation='nearest',origin='upper')
        #im=ax.imshow(np.sqrt((gh.take(indices=0,axis=2).take(indices=slicepos,axis=sliceaxis))**2),interpolation='nearest',origin='upper')
        #im=ax.imshow(np.sqrt((v.take(indices=0,axis=3).take(indices=slicepos,axis=sliceaxis))**2+(v.take(indices=1,axis=3).take(indices=slicepos,axis=sliceaxis))**2),interpolation='nearest',origin='upper')
        #print(np.flip(rho.take(indices=slicepos,axis=sliceaxis)).transpose()[70,70])
        #ax.scatter(70,70)
        
        fig,ax=plt.subplots(1,1,figsize=(6,6))
        im=ax.imshow(rgbv,interpolation='nearest',origin='upper')
        ax.contour(np.flip(liquid).T, levels=[0.5], colors="k", zorder=1, linewidths=0.75)
        stepx=1
        stepy=1
        X,Y=np.meshgrid(np.linspace(0,LX-1,int((LX)/stepx)),np.linspace(0,LY-1,int((LY)/stepy)))
        ax.quiver(X.T,Y.T,np.flip(-v[0:LX:stepx,0:LY:stepy,0]),np.flip(v[0:LX:stepx,0:LY:stepy,1]),width=0.0008,headwidth=7.5,headlength=7.5)
        fig.colorbar(im)
        plt.savefig(output, dpi=400, format='png')
        plt.close(fig)
        
        #print(np.sum(np.sqrt((gh.take(indices=0,axis=3).take(indices=slicepos,axis=sliceaxis))**2+(gh.take(indices=1,axis=3).take(indices=slicepos,axis=sliceaxis))**2)))
        #ax.quiver(X.T,Y.T,np.flip(v[:,:,0].take(indices=slicepos,axis=sliceaxis)),np.flip(-v[:,:,1].take(indices=slicepos,axis=sliceaxis)),width=0.001,headwidth=2.5,headlength=1.5)
        #ax.quiver(X.T,Y.T,np.flip(-v[0:-1:stepx,0:-1:stepy,0]),np.flip(v[0:-1:stepx,0:-1:stepy,1]),width=0.0002,headwidth=7.5,headlength=7.5)
        
        #print(np.sum(rho[int(LX/2),:])/(0.002/(100-h)*np.log(1/(1-0.1))))
        #print("V ",np.amax(v))
        #print("HERE ",np.sum(rho))
        #print("HERE ",np.sum(rho>0.5))
        #print("HERE ",rho[LX//8,LY//2])
        #print("HERE ",rho[LX//2,LY//6])
        
        #ax.scatter(49,49)

        print("V",np.amax(v))
        print("V",np.average(v[LX//4,:]))
        print("V",np.average(v[LX//2,:]))
        print("V",np.average(v[3*LX//4,:]))
        
        pw=80
        ym=35
        ymin=15
        bulkfreenergy = liquid[pw:(228-pw),ymin:ym]**2*(1-liquid[pw:(228-pw),ymin:ym])**2
        gradx = np.gradient(liquid[pw:(228-pw),ymin:ym], axis = 0)
        grady = np.gradient(liquid[pw:(228-pw),ymin:ym], axis = 1)
        gradx[np.where(np.logical_or(solid[pw:(228-pw),ymin:ym]==-1,solid[pw:(228-pw),ymin:ym]>=1))[0],np.where(np.logical_or(solid[pw:(228-pw),ymin:ym]==-1,solid[pw:(228-pw),ymin:ym]>=1))[1]] = 0
        grady[np.where(np.logical_or(solid[pw:(228-pw),ymin:ym]==-1,solid[pw:(228-pw),ymin:ym]>=1))[0],np.where(np.logical_or(solid[pw:(228-pw),ymin:ym]==-1,solid[pw:(228-pw),ymin:ym]>=1))[1]] = 0
        surfacearea = np.sum(4*bulkfreenergy + 0.75 * 3 * (gradx**2 + grady**2))
        #plt.figure()
        #plt.plot(rho[int(LX/2),:])
        #plt.plot(rho[:,int(LY/2),0])
        #plt.savefig("test_%012d.png"%(t), dpi=200, format='png')
        height = np.append(height, h)
        #print(gh[int(LX/2),51,1])
        vol=np.append(vol,np.sum(liquid))#>=0.5))
        heightavg = np.append(heightavg,np.average(np.where(np.logical_and((liquid[:,:] <= 0.8),(liquid[:,:] > 0.2)))[1]))

        #ml=np.append(ml,surfacearea)
        if (th<90):
            ml=np.append(ml,np.sum(mlnosolid[:,:])/surfacearea)
            #ml=np.append(ml,surfacearea)
        else:
            ml=np.append(ml,np.sum(mlnosolid[:,:])/surfacearea)
            #ml=np.append(ml,surfacearea)
        print(np.average(v[4,:,0]))
        print(np.average(v[LX-8,:,0]))
        print(np.sum(liquid))
    print(ml)
    fit=optimize.curve_fit(fitfunc, np.linspace(tstart,tend,int((tend-tstart)/tinc)+1), vol, p0=[10000,-10,0.5],maxfev=1000000)[0]
    #fit=optimize.curve_fit(derivfunc, np.linspace(tstart,tend,int((tend-tstart)/tinc)+1), ml, p0=[0.01949235,0.79008422],maxfev=500000)[0]
    print(fit)
    print(vol)
    dheightavg[th] = heightavg
    dml[th] = -derivfunc(np.linspace(tstart,tend,int((tend-tstart)/tinc)+1),fit[1],fit[2])#/ml#derivfunc(np.linspace(tstart,tend,int((tend-tstart)/tinc)+1),fit[0],fit[1])#ml#-derivfunc(np.linspace(tstart,tend,int((tend-tstart)/tinc)+1),fit[0],fit[1],fit[2])#ml
    dvol[th] = vol

plt.figure()
for th in thetas:
    plt.plot(dheightavg[th],dml[th],label="Contact Angle: "+str(th))
plt.xlabel("Height")
plt.ylabel("dV/dt")# per unit area")
plt.title("Height Measured as Average Height")
plt.gca().invert_xaxis()
plt.legend()
plt.tight_layout()
#plt.savefig("ContactLine_fluxvsheightNOTNORMALISED.PNG",dpi=200,format='png')
plt.savefig("0.0002.PNG",dpi=200,format='png')

plt.figure()
plt.plot(rho2[5*int(LX/6),LY//3+5:LY-10])
plt.savefig("test.png", dpi=200, format='png')


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
