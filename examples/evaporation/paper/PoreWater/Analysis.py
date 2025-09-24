import os, math, re, sys
import struct
import numpy as np
from matplotlib import pyplot as plt
from scipy import optimize
import scipy.ndimage.filters as filters
import matplotlib
import matplotlib.patches as mpatches
import matplotlib.path as mpath
from scipy import interpolate


height = np.array([])

def mlfit(x):
    return (1-(x/np.amax(x))**2)**(1/6-0.5)

def fitfunc(x, a, b, c):
    return a+b*x**c
def derivfunc(x, b, c,d):
    return c*b*x**(c-1) + d
def derivfunc2(aa, x,y):
    #print("test",aa[1]*aa[0]*x**(aa[1]-1)+aa[2] - y)
    return aa[1]*aa[0]*x**(aa[1]-1)+aa[2] - y

outDirName = "figures"
os.system("mkdir -p %s"%outDirName)

#thetas = [150,120,90,60,30]
thetas = [150]
factors = [1]
plot=True
dml = {}
dheightavg = {}
dvol = {}
fml = {}
fheightavg = {}
fvol = {}

"""fig,ax=plt.subplots(1,1,figsize=(6,6))
flxs = (np.array([[1.09060120017596,1.084191929488,1.06924895305013,1.0411049227347],[1.09347989372926,1.08374013792761,1.06882588385361,1.04614637858097],[1.089214646939,1.08181165322233,1.06320238108245,1.04147010784058]])-1)*100
im=ax.imshow(flxs,interpolation='nearest',origin='upper')
ax.set_xticks([0,1,2,3])
ax.set_xticklabels(["80%","65%","45%","30%"])
ax.set_xlabel("Pore width/unit cell width")
ax.set_yticks([2,1,0])
ax.set_yticklabels(["15%","35%","55%"])
ax.set_ylabel("Outflow distance/pore height")
cbar = fig.colorbar(im)
cbar.set_label(r"% increase in flux for $30\degree$ vs $90\degree$")
plt.savefig("flxs.png", dpi=400, format='png')
plt.close(fig)
sys.exit()"""
fig,ax=plt.subplots(1,1,figsize=(6,6))
flxs = (np.array([[0.00125745108378601,0.00125745108378601,0.00125378867768036
],[0.000820247656824959,0.000848380836081792,0.000850900302529716
],[0.000555180717467066,0.000619647921112395,0.000651638217006107
]]))
im=ax.imshow(flxs,interpolation='nearest',origin='upper')
ax.set_xticks([0,1,2])
ax.set_xticklabels(["18%","65%","125%"])
ax.set_xlabel("Post width/pore width")
ax.set_yticks([0,1,2])
ax.set_yticklabels(["120%","165%","230%"])
ax.set_ylabel("Outflow distance/pore height")
cbar = fig.colorbar(im)
cbar.set_label("Flux")
plt.savefig("flxs.png", dpi=400, format='png')
plt.close(fig)
#sys.exit()
#a=["a","aa","aaa","b","bb","bbb","c","cc","ccc","d","dd","ddd","e","ee","eee","f","ff","fff","g","gg","ggg","h","hh","hhh","i","ii","iii","j","jj","jjj","k","kk","kkk","l","ll","lll","m","mm","mmm","n","nn","nnn","o","oo","ooo","p","pp","ppp","q","qq","qqq","r","rr","rrr","s","ss","sss","t","tt","ttt","u","uu","uuu","v","vv","vvv","w","ww","www","x","xx","xxx","y","yy","yyy","z","zz""zzz"]
wid = 12
hei = 195
for fact in factors:
    datadir = "data/diffusive7/postheight_46/surfacetension_0.001/theta_30/offsety_-17/postwidth_42/ly_100/lx_108/diffusivity_0.02141/factor_"+str(fact)+"/"
    datadir = "data/diffusive8/postheight_46/surfacetension_0.001/theta_150/offsety_-17/postwidth_82/ly_150/lx_148/diffusivity_0.02141/factor_2/"
    datadir = "data/diffusive7/postheight_46/surfacetension_0.001/theta_30/offsety_-17/postwidth_42/ly_100/lx_108/diffusivity_0.02141/factor_"+str(fact)+"/"
    datadir = "data/diffusive8/postheight_46/surfacetension_0.001/theta_90/offsety_-17/postwidth_82/ly_120/lx_148/diffusivity_0.02141/factor_"+str(fact)+"/"
    datadir = "data/diffusive9/postheight_46/surfacetension_0.001/theta_150/offsety_-17/postwidth_82/ly_106/lx_148/diffusivity_0.02141/factor_1/"
    datadir = "data/diffusive9/postheight_46/surfacetension_0.001/theta_30/offsety_-17/postwidth_"+str(wid)+"/ly_"+str(hei)+"/lx_"+str(66+wid)+"/diffusivity_0.02141/factor_1/"
    datadir = "data/diffusive9/postheight_46/surfacetension_0.001/theta_150/offsety_-17/postwidth_"+str(wid)+"/ly_"+str(hei)+"/lx_"+str(66+wid)+"/diffusivity_0.02141/factor_"+str(fact)+"/"
    datadir = "data/paper/wb/factor_1/diffusivity_0.02141/lx_"+str(66+wid)+"/ly_"+str(hei)+"/postwidth_"+str(wid)+"/offsety_-15/theta_30/surfacetension_0.001/postheight_76/outflowoffset_0/"
    datadir= "data/paper/phasediagram/factor_2/diffusivity_0.02141/lx_138/ly_195/postwidth_72/offsety_0/theta_150/surfacetension_0.001/postheight_76/outflowoffset_180/"
    datadir = "data/outflowoffset_190/postheight_132/surfacetension_0.00111/theta_30/offsety_-15/postwidth_24/ly_370/lx_156/diffusivity_0.16666667/factor_1/"
    #datadir ="data/paper/cap/factor_1/diffusivity_0.02141/lx_138/ly_195/postwidth_72/offsety_-35/theta_150/surfacetension_0.001/postheight_76/outflowoffset_0/capwidth_56/capthickness_10/"
    #datadir = "data/paper/pinnedkinetic/equilibriumtimesteps_10000/postheight_96/surfacetension_0.01/theta_30/offsety_0/postwidth_42/ly_200/lx_108/diffusivity_0.02141/factor_1/"
    #datadir = "data/diffusive8/postheight_46/surfacetension_0.001/theta_90/offsety_-17/postwidth_82/ly_120/lx_148/diffusivity_0.02141/factor_"+str(fact)+"/"
    #datadir = "data/diffusive6/postheight_96/surfacetension_0.001/theta_30/offsety_-17/postwidth_162/ly_200/lx_228/diffusivity_0.02141/factor_"+str(fact)+"/"
    #datadir = "data/diffusive6/postheight_46/surfacetension_0.001/theta_30/offsety_-17/postwidth_82/ly_100/lx_148/diffusivity_0.02141/factor_1/"
    #postheight_96/surfacetension_0.001/theta_150/offsety_-17/postwidth_42/ly_110/lx_228/diffusivity_0.08141/factor_1/
    #datadir = "data/inflowsimple7/surfacetension_0.0001/theta_90/offsety_-12/postwidth_90/ly_133/lx_156/pDiff_1e-05/diffusivity_0.006/"
    #datadir = "data/inflowsimple7/surfacetension_0.0001/theta_90/offsety_-36/postwidth_270/ly_399/lx_468/pDiff_3.3333333333333337e-06/diffusivity_0.018000000000000002/"
    HeaderFile = open(datadir+"Header.mat", 'rb')

    LX=struct.unpack('=i', HeaderFile.read(4))[0]

    LY=struct.unpack('=i', HeaderFile.read(4))[0]

    LZ=struct.unpack('=i', HeaderFile.read(4))[0]

    ndim=struct.unpack('=i', HeaderFile.read(4))[0]

    t_zero = 0
    tstart =150000

    tend = [150000,15000,15000,700000]
    struct.unpack('=i', HeaderFile.read(4))[0]
    tinc = 10000
    struct.unpack('=i', HeaderFile.read(4))[0]
    v = np.zeros((LX,LY,LZ,ndim))
    for th in thetas:
        heightavg = np.array([])
        ml = np.array([])
        vol = np.array([])
        #datadir = "data/inflowsimple/inflowmomentum_-0.0002/lx_228/ly_200/postwidth_162/offsety_-17/theta_"+str(th)+"/"
        #datadir = "data/inflow3/diffusivity_0.0008/inflowmomentum_-0.0002/lx_228/ly_200/postwidth_162/offsety_-17/theta_"+str(th)+"/"
        #datadir = "data/inflowsimple7/surfacetension_0.0001/theta_"+str(th)+"/offsety_-48/postwidth_360/ly_532/lx_624/inflowmomentum_-0.0006/diffusivity_0.016/"
        #datadir = "data/inflowsimple7/surfacetension_0.0001/theta_"+str(th)+"/offsety_-12/postwidth_90/ly_133/lx_156/inflowmomentum_-0.0006/diffusivity_0.004/"
        
        #datadir = "data/inflow8/surfacetension_0.0001/theta_"+str(th)+"/offsety_-24/postwidth_180/ly_266/lx_312/pDiff_5e-06/diffusivity_0.012/"
        datadir = "data/diffusive/surfacetension_0.001/theta_"+str(th)+"/offsety_-17/postwidth_162/ly_200/lx_228/diffusivity_0.08141/"
        datadir = "data/diffusive3/postheight_96/surfacetension_0.001/theta_"+str(th)+"/offsety_-17/postwidth_162/ly_150/lx_228/diffusivity_0.08141/factor_"+str(fact)+"/"
        datadir = "data/diffusive7/postheight_46/surfacetension_0.001/theta_"+str(th)+"/offsety_-17/postwidth_42/ly_100/lx_108/diffusivity_0.02141/factor_"+str(fact)+"/"
        datadir = "data/diffusive8/postheight_46/surfacetension_0.001/theta_"+str(th)+"/offsety_-17/postwidth_82/ly_120/lx_148/diffusivity_0.02141/factor_"+str(fact)+"/"
        datadir = "data/diffusive9/postheight_46/surfacetension_0.001/theta_150/offsety_-17/postwidth_82/ly_106/lx_148/diffusivity_0.02141/factor_1/"
        datadir = "data/diffusive9/postheight_46/surfacetension_0.001/theta_"+str(th)+"/offsety_-17/postwidth_"+str(wid)+"/ly_"+str(hei)+"/lx_"+str(66+wid)+"/diffusivity_0.02141/factor_1/"
        datadir = "data/diffusive9/postheight_46/surfacetension_0.001/theta_"+str(th)+"/offsety_-17/postwidth_"+str(wid)+"/ly_"+str(hei)+"/lx_"+str(66+wid)+"/diffusivity_0.02141/factor_"+str(fact)+"/"
        datadir = "data/paper/wb/factor_1/diffusivity_0.02141/lx_"+str(66+wid)+"/ly_"+str(hei)+"/postwidth_"+str(wid)+"/offsety_-15/theta_"+str(th)+"/surfacetension_0.001/postheight_76/outflowoffset_0/"
        datadir = "data/paper/phasediagram/factor_2/diffusivity_0.02141/lx_138/ly_195/postwidth_72/offsety_0/theta_"+str(th)+"/surfacetension_0.001/postheight_76/outflowoffset_180/"
        datadir = "data/outflowoffset_190/postheight_132/surfacetension_0.00111/theta_30/offsety_-15/postwidth_24/ly_370/lx_156/diffusivity_0.16666667/factor_1/"
        #datadir = "data/paper/pinnedkinetic/equilibriumtimesteps_10000/postheight_96/surfacetension_0.01/theta_"+str(th)+"/offsety_0/postwidth_42/ly_200/lx_108/diffusivity_0.02141/factor_1/"
        #datadir = "data/diffusive8/postheight_46/surfacetension_0.001/theta_"+str(th)+"/offsety_-17/postwidth_82/ly_120/lx_148/diffusivity_0.02141/factor_"+str(fact)+"/"
        #datadir = "data/diffusive8/postheight_46/surfacetension_0.001/theta_150/offsety_-17/postwidth_82/ly_150/lx_148/diffusivity_0.02141/factor_2/"
        #datadir = "data/diffusive6/postheight_96/surfacetension_0.001/theta_"+str(th)+"/offsety_-17/postwidth_162/ly_200/lx_228/diffusivity_0.02141/factor_"+str(fact)+"/"
        #datadir = "data/diffusive6/postheight_46/surfacetension_0.001/theta_30/offsety_-17/postwidth_82/ly_100/lx_148/diffusivity_0.02141/factor_1/"
        #datadir = "data/inflowsimple7/surfacetension_0.0001/theta_"+str(th)+"/offsety_-12/postwidth_90/ly_133/lx_156/pDiff_1e-05/diffusivity_0.006/"
        #datadir = "data/inflowsimple7/surfacetension_0.0001/theta_"+str(th)+"/offsety_-36/postwidth_270/ly_399/lx_468/pDiff_3.3333333333333337e-06/diffusivity_0.018000000000000002/"

        #datadir = "data/inflowsimple2/inflowmomentum_-0.0002/lx_228/ly_200/postwidth_162/offsety_-17/theta_"+str(th)+"/"
        #datadir = "data/inflow3/inflowmomentum_-0.002/lx_228/ly_300/postwidth_162/offsety_-17/theta_"+str(th)+"/"
        #datadir = "data/paper/cap/factor_1/diffusivity_0.02141/lx_138/ly_195/postwidth_72/offsety_-15/theta_"+str(th)+"/surfacetension_0.001/postheight_76/outflowoffset_0/capwidth_56/capthickness_10/"
        for t in range(tstart,tend[fact-1]+1,tinc):
            #if th==150 and t>=1125000:
            #    break
            print("t=%s"%t)
            t_file =t+t_zero

            file_name = datadir+"OrderParameter_t%li.mat"%t_file
            #file_name = datadir+"Pressure_t%li.mat"%t_file
            #file_name = datadir+"Density_t%li.mat"%t_file
            #file_name = "data/"+"Humidity_t%li.mat"%t_file
            #file_name = datadir+"Pressure_t%li.mat"%t_file
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

            file_name = datadir+"GradientHumidity_t%li.mat"%t_file

            File3 = open(file_name, 'rb')

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
            
            solid = np.zeros((LX,LY,LZ))

            dat=File.read()
            rho = np.ndarray((LX,LY),'=d',dat,0,(8*LY,8))

            dat=FileSolid.read()
            solid = np.ndarray((LX,LY),'=i',dat,0,(4*LY,4))

            dat=File0.read()
            rho2 = np.array(np.ndarray((LX,LY),'=d',dat,0,(8*LY,8)))
            rho2[np.where(np.logical_and(solid!=0,solid!=5))[0],np.where(np.logical_and(solid!=0,solid!=5))[1]] = 0
            
            #rho = np.ndarray((LX,LY),'=i',dat,0,(4*LY,4))

            dat=File00.read()
            humidity = np.ndarray((LX,LY),'=d',dat,0,(8*LY,8))

            

            liquid = np.array(rho)
            liquid[np.where(np.logical_or(solid==1,solid==-1))[0],np.where(np.logical_or(solid==1,solid==-1))[1]] = 0

            mlnosolid = np.array(humidity)
            mlnosolid[np.where(np.logical_or(solid==1,solid==-1))[0],np.where(np.logical_or(solid==1,solid==-1))[1]] = 0
            
            dists = -0.5 * 3.0 * np.arctanh((liquid - 0.5) / 0.5)#

            #idcs=np.where(np.logical_and(-0.5 * 3.0 * np.arctanh((liquid - 0.5) / 0.5)>7,-0.5 * 3.0 * np.arctanh((liquid - 0.5) / 0.5)<7.5))
            """print(idcs)
            plt.figure()
            plt.plot(mlnosolid[idcs[0],idcs[1]])
            plt.savefig("testdist.png")
            plt.close()"""
            v = np.zeros((LX,LY,LZ,ndim))
            dat=File2.read()
            v = np.ndarray((LX,LY,ndim),'=d',dat,0,(ndim*8*LY,ndim*8,8))

            dat=File3.read()
            gh = np.array(np.ndarray((LX,LY,ndim),'=d',dat,0,(ndim*8*LY,ndim*8,8)))
            gh[np.where(solid!=0)[0],np.where(solid!=0)[1],:] = 0

            dat=File4.read()
            c = np.ndarray((LX,LY),'=d',dat,0,(8*LY,8))

            File.close()
            File2.close()
            File3.close()
            File4.close()
            
            

            #output = "%s/"%(outDirName)+str(a[t//25000])+"component_plot_%012d.png"%(t)
            output = "%s/component_plot_%012d.png"%(outDirName,t)
            #rho3=2*0.01*(rho2-0.2)*(rho2-1)*(2*rho2-0.2-1)-0.0128*rho4
            if plot==True:
                #rgbv = np.zeros((2*LY//3-21,LX))
                rgbv = np.zeros((LY,LX))
                #lap = filters.laplace(mlnosolid[:,:LY//3*3])[1:LX-1,LY//3+2:LY-5]
                #rgbv[:,:] = np.flip(lap).transpose()
                #rgbv2 = np.zeros((LY//3-3,111-46))
                #lap = filters.laplace(mlnosolid[:,:LY//3*3])[46:111,5:LY//3+2]
                #rgbv2[:,:] = np.flip(lap).transpose()
                #rgbv[:,:] = np.flip(np.sqrt(gh[:,:,0]**2+gh[:,:,1]**2)).transpose()
                #rgbv[:,:] = 0.03*np.flip(mlnosolid[:,LY//3+10:LY//3*3-10]).transpose()#+np.flip(gh[:,LY//3+10:LY//3*3-10,0]*v[:,LY//3+10:LY//3*3-10,0]+gh[:,LY//3+10:LY//3*3-10,1]*v[:,LY//3+10:LY//3*3-10,1]).transpose()
                #rgbv[:,:] = np.flip(rho2[:,:]).transpose()
                rgbv[:,:] = np.flip(liquid[:,:]*(liquid[:,:]<0)).transpose()
            #rgbv[:,:] = np.flip(rho.take(indices=slicepos,axis=sliceaxis)).transpose()0.004081596000898668,0.003980521608340414,0.003899587801806342,0.0039410080368059265,0.0038240387164630843
            #im=ax.imshow(np.flip(rho0.take(indices=slicepos,axis=sliceaxis)).transpose(),interpolation='nearest',origin='upper')
            i2 = np.argmax(np.logical_and((c[0,:] <= 0.5),(c[0,:] > 0.1)))
            print("SUM",np.sum(liquid))
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
            if plot==True:
                poly_verts = [
                    (0,47),
                    (0, 70),
                    (LX-1, 70),
                    (LX-1, 47),
                    (87,47),
                    (87,2),
                    (21,2),
                    (21,2),
                    (21,47),
                    (0,47)
                ]
                poly_codes = [mpath.Path.MOVETO] + (len(poly_verts) - 2) * [mpath.Path.LINETO] + [mpath.Path.CLOSEPOLY]
                path = mpath.Path(poly_verts, poly_codes)
                patch = mpatches.PathPatch(path, facecolor='none', edgecolor='k')
                fig,ax=plt.subplots(1,1,figsize=(6,6))
                im=ax.imshow(rgbv,interpolation='nearest',origin='upper')
                plt.colorbar(im, orientation='vertical', pad=0.01, aspect=50, label="Humidity")
                #im = ax.contour(np.flip(rgbv),colors="k")
                #ax.add_patch(patch)  ## TRY COMMENTING THIS OUT
                #for col in im.collections:
                #    col.set_clip_path(patch)
                
                """print(im.levels)
                p = im.collections[1].get_paths()[0]
                print(p)
                vv = p.vertices#.vertices
                x = vv[:,0]
                y = vv[:,1]"""

                #print(x[65],y[65])
                #print(np.interp(x*LY+y,np.linspace(0,LX*LY-1,LX*LY),mlnosolid.flatten()))

                """plt.figure()
                plt.plot(np.round(x),np.interp(np.round(x)*LY+np.round(y),np.linspace(0,LX*LY-1,LX*LY),mlnosolid.flatten()))
                print(th,''.join([str(i)+"," for i in np.round(x)]),"\n",''.join([str(i)+"," for i in np.interp(np.round(x)*LY+np.round(y),np.linspace(0,LX*LY-1,LX*LY),mlnosolid.flatten())]))
                #print(interpolate.CubicSpline(np.linspace(0,LX*LY-1,LX*LY),mlnosolid.flatten()))
                #plt.plot(interpolate.CubicSpline(np.linspace(0,LX*LY-1,LX*LY),mlnosolid.flatten())(x*LY+y))
                plt.savefig("testdist.png")
                plt.close()"""
                #plt.clabel(im, fontsize=6)
                #plt.axis("scaled")
                #ax.scatter(np.where(np.logical_and((liquid[:,:] <= 0.6),(liquid[:,:] > 0.4)))[0],np.where(np.logical_and((liquid[:,:] <= 0.6),(liquid[:,:] > 0.4)))[1],s=1)
                #im=ax.contourf(rgbv,interpolation='nearest',origin='upper')
                #im=ax.imshow(rgbv2,interpolation='nearest',origin='upper',extent=[46,111, 0, LY//3+3])
                im=ax.imshow(rgbv,interpolation='nearest',origin='upper')#,extent=[0, LX, LY//3+3, LY])
                #ax.set_ylim([0,LY])
                #print(np.average(np.where(np.logical_and((liquid[21,:] <= 0.7),(liquid[21,:] > 0.3)))))
                #liq=ax.contour(np.flip(np.flip(liquid).T), levels=[0.5], colors="r",linestyles = 'dashed')#,extent=[1, LX-1, int(np.average(np.where(np.logical_and((liquid[21,:] <= 0.7),(liquid[21,:] > 0.3))))), LY//2])#, zorder=1, linewidths=0.75,extent=[1, LX-1, 0, LY//2])
                liq=ax.contour(np.flip(solid).T, levels=[0.5], colors="k",linestyles = 'dashed')
                #for col in liq.collections:
                #    col.set_clip_path(patch)
                stepx=1
                stepy=1
                #print(LX,LY)
                X, Z = np.meshgrid(np.linspace(0, LX - 1, int((LX) / stepx)), np.linspace(0, LY - 1, int((LY) / stepy)))
                # Plot velocity arrows
                #ax.quiver(X.T, Z.T, np.flip(np.flip(-v[0:LX:stepx, 0:LY:stepy,0], 0), 1),np.flip(v[0:LX:stepx, 0:LY:stepy, 1]))#, width=0.0008, headwidth=7.5, headlength=7.5)

                plt.savefig(output, dpi=400, format='png')
                plt.close(fig)
            print("HEIGHT ",np.sum(liquid)/(66*fact))#np.average(np.where(np.logical_and((liquid[:,:] <= 0.8),(liquid[:,:] > 0.2)))[1]))
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
            #bulkfreenergy = liquid[pw:(228-pw),ymin:ym]**2*(1-liquid[pw:(228-pw),ymin:ym])**2
            #gradx = np.gradient(liquid[pw:(228-pw),ymin:ym], axis = 0)
            #grady = np.gradient(liquid[pw:(228-pw),ymin:ym], axis = 1)
            #gradx[np.where(np.logical_or(solid[pw:(228-pw),ymin:ym]==-1,solid[pw:(228-pw),ymin:ym]>=1))[0],np.where(np.logical_or(solid[pw:(228-pw),ymin:ym]==-1,solid[pw:(228-pw),ymin:ym]>=1))[1]] = 0
            #grady[np.where(np.logical_or(solid[pw:(228-pw),ymin:ym]==-1,solid[pw:(228-pw),ymin:ym]>=1))[0],np.where(np.logical_or(solid[pw:(228-pw),ymin:ym]==-1,solid[pw:(228-pw),ymin:ym]>=1))[1]] = 0
            #surfacearea = np.sum(4*bulkfreenergy + 0.75 * 3 * (gradx**2 + grady**2))
            #plt.figure()
            #plt.plot(rho[int(LX/2),:])
            #plt.plot(rho[:,int(LY/2),0])
            #plt.savefig("test_%012d.png"%(t), dpi=200, format='png')
            #if (np.isnan(mlnosolid[:,:]).any()):
            #    break
            treached = t
            if (np.sum(liquid)/(66*fact)<20):
                break
            height = np.append(height, h)
            #print(gh[int(LX/2),51,1])
            
            vol=np.append(vol,np.sum(liquid))
            heightavg = np.append(heightavg,np.sum(liquid)/(66*fact))#np.mean(np.where(np.logical_and((liquid[:,:] <= 0.8),(liquid[:,:] > 0.2)))[1]))
            
            #ml=np.append(ml,surfacearea)
            if (th<90):
                ml=np.append(ml,np.sum(mlnosolid[:,:]*(mlnosolid[:,:]>0)))#/surfacearea)
                #ml=np.append(ml,surfacearea)
            else:
                ml=np.append(ml,np.sum(mlnosolid[:,:]*(mlnosolid[:,:]>0)))#/surfacearea)
                #ml=np.append(ml,surfacearea)
            print(np.average(v[4,:,0]))
            print(np.average(v[LX-8,:,0]))
            print(np.sum(liquid))
        
        plt.figure()
        #plt.plot(mlnosolid[44,:])
        #plt.plot(mlnosolid[45,2:])
        #plt.plot(np.diff(mlnosolid[45,2:]))
        #plt.plot(mlnosolid[46,:])
        #plt.plot(gh[45,2:,1])
        plt.plot(vol)
        plt.savefig("test.png", dpi=200, format='png')
        print("MASSLOSS",ml)
        print(np.amin(liquid)*(77.87604605263631-0.08396117887190388)+0.08396117887190388)
        print(np.amin(liquid))
        #fit=optimize.curve_fit(fitfunc, np.linspace(tstart,treached,int((treached-tstart)/tinc)+1), vol, p0=[10000,-10,0.5],maxfev=1000000)[0]
        #fit=optimize.curve_fit(derivfunc, np.linspace(tstart,tend,int((tend-tstart)/tinc)+1), ml, p0=[0.01949235,0.79008422,0],maxfev=500000)[0]
        #fit=optimize.least_squarerho2s(derivfunc2, [0.01,0.01,0], loss='cauchy',max_nfev=500000, args = (np.linspace(tstart,tend,int((tend-tstart)/tinc)+1),ml)).x
        #print(fit)
        print(vol)
        dheightavg[th] = heightavg
        dml[th] = ml#derivfunc(np.linspace(tstart,tend,int((tend-tstart)/tinc)+1),fit[0],fit[1],fit[2])#ml#derivfunc(np.linspace(tstart,tend,int((tend-tstart)/tinc)+1),fit[0],fit[1],fit[2])#-derivfunc(np.linspace(tstart,treached,int((treached-tstart)/tinc)+1),fit[1],fit[2])#/ml#derivfunc(np.linspace(tstart,tend,int((tend-tstart)/tinc)+1),fit[0],fit[1])#ml#-derivfunc(np.linspace(tstart,tend,int((tend-tstart)/tinc)+1),fit[0],fit[1],fit[2])#ml
        dvol[th] = vol
        #plt.figure()
        #plt.plot(rho2[int(LX/2),:])#,0])
        #plt.savefig("test"+str(th)+".png", dpi=200, format='png')
    #dheightavg2 = dheightavg.copy()
    fheightavg[fact] = dheightavg.copy()
    fml[fact] = dml.copy()#derivfunc(np.linspace(tstart,tend,int((tend-tstart)/tinc)+1),fit[0],fit[1],fit[2])#ml#derivfunc(np.linspace(tstart,tend,int((tend-tstart)/tinc)+1),fit[0],fit[1],fit[2])#-derivfunc(np.linspace(tstart,treached,int((treached-tstart)/tinc)+1),fit[1],fit[2])#/ml#derivfunc(np.linspace(tstart,tend,int((tend-tstart)/tinc)+1),fit[0],fit[1])#ml#-derivfunc(np.linspace(tstart,tend,int((tend-tstart)/tinc)+1),fit[0],fit[1],fit[2])#ml
    fvol[fact] = dvol.copy()
plt.figure()
print(np.amin(liquid)*(1000-1)+1)
plt.plot(liquid[13,:])#,0])
plt.savefig("phi.png", dpi=200, format='png')
#print("HEREE",'\n'.join([''.join([str(i)+"," for i in dheightavg[30]]),''.join([str(i)+"," for i in dml[30]]),''.join([str(i)+"," for i in dheightavg[90]]),''.join([str(i)+"," for i in dml[90]]),''.join([str(i)+"," for i in dheightavg[150]]),''.join([str(i)+"," for i in dml[150]])]))
#print("AVERAGERATE",np.average(fml[1][90]))
#sys.exit()
def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

h30=[21.0,22.0,22.0,23.0,23.0,24.0,24.0,25.0,25.0,26.0,26.0,27.0,27.0,28.0,28.0,29.0,29.0,30.0,31.0,31.0,32.0,32.0,33.0,33.0,34.0,35.0,35.0,36.0,37.0,37.0,38.0,39.0,39.0,40.0,41.0,42.0,42.0,43.0,44.0,45.0,45.0,46.0,47.0,48.0,49.0,50.0,51.0,52.0,53.0,54.0,55.0,56.0,57.0,58.0,59.0,60.0,61.0,62.0,63.0,63.0,64.0,65.0,66.0,66.0,67.0,68.0,69.0,69.0,70.0,71.0,71.0,72.0,73.0,73.0,74.0,75.0,75.0,76.0,76.0,77.0,77.0,78.0,79.0,79.0,80.0,80.0,81.0,81.0,82.0,82.0,83.0,83.0,84.0,84.0,85.0,85.0,86.0,86.0,87.0]
f30=[0.12561889383399058,0.1421398812174986,0.1421398812174986,0.1536602051562387,0.1536602051562387,0.16245759236426144,0.16245759236426144,0.16834394944139547,0.16834394944139547,0.17269597689805313,0.17269597689805313,0.17592690121949273,0.17592690121949273,0.17864574809218217,0.17864574809218217,0.18095535613878524,0.18095535613878524,0.18308461626349845,0.18084925577825087,0.18084925577825087,0.18300067022041455,0.18300067022041455,0.18505020314478943,0.18505020314478943,0.18332098801347185,0.18543399802509986,0.18543399802509986,0.18749818405804808,0.1861372101154809,0.1861372101154809,0.18827845021450865,0.18712681337393597,0.18712681337393597,0.1860771293978384,0.18837400875651972,0.18749902265133922,0.18749902265133922,0.18670633885234283,0.18913879935078898,0.1885015931179556,0.1885015931179556,0.18793831478381876,0.18744579690619764,0.1900758110149338,0.18972060166134647,0.1894307983339843,0.1892061799056976,0.18904727661814522,0.18895213492009802,0.18892053540638915,0.1889523256295686,0.1890476602870914,0.18920676022603436,0.18943158404985486,0.18972159705648628,0.19007701411814631,0.18744719387941272,0.18793993081629273,0.18850342925864527,0.18850342925864527,0.18914084762099445,0.18670855262757427,0.18750143151905796,0.18750143151905796,0.1883765814873688,0.18607971468644421,0.18712948151160733,0.18712948151160733,0.1882812189048036,0.18613973681356874,0.18613973681356874,0.18750081431783178,0.18543638649361985,0.18543638649361985,0.1833232972823184,0.18505288251503793,0.18505288251503793,0.18300367550968039,0.18300367550968039,0.1808528882345651,0.1808528882345651,0.18308932470676864,0.1809610971253888,0.1809610971253888,0.17865293865508314,0.17865293865508314,0.17593563441470456,0.17593563441470456,0.17270658053573687,0.17270658053573687,0.16835714329580184,0.16835714329580184,0.16247189829907838,0.16247189829907838,0.15367384152234093,0.15367384152234093,0.1421539898910962,0.1421539898910962,0.1256308556992925]
h90=[19.0,20.0,21.0,22.0,23.0,24.0,25.0,26.0,27.0,28.0,29.0,30.0,31.0,32.0,33.0,34.0,35.0,36.0,37.0,38.0,39.0,40.0,41.0,42.0,43.0,44.0,45.0,46.0,47.0,48.0,49.0,50.0,51.0,52.0,53.0,54.0,55.0,56.0,57.0,58.0,59.0,60.0,61.0,62.0,63.0,64.0,65.0,66.0,67.0,68.0,69.0,70.0,71.0,72.0,73.0,74.0,75.0,76.0,77.0,78.0,79.0,80.0,81.0,82.0,83.0,84.0,85.0,86.0,87.0,88.0,89.0]
f90=[0.17665464365555558,0.17665464365555558,0.17665464365555558,0.17669486064069345,0.17677417954387495,0.17688229583384055,0.1770248338862673,0.17716242682716446,0.17732732797402337,0.17750654339192518,0.17768971653351223,0.17788101245057986,0.17807195225008066,0.17826357321267788,0.1784513733962093,0.1786361544011204,0.1788140255861049,0.17898644483217438,0.179150143067022,0.179306169113822,0.17945311519341306,0.1795908909701086,0.1797194078604124,0.17983819529412623,0.17994749615655942,0.18004700363404164,0.18013690748812752,0.1802171145315265,0.18028773240527332,0.18034878120383743,0.18040032689887486,0.18044241996029295,0.18047510714188117,0.18049842920165815,0.18051241550020575,0.18051708373072026,0.18051244113161968,0.18049848030302276,0.18047518339185845,0.18044252088043597,0.1804004518584691,0.180348929425324,0.180287902971686,0.18021730639476327,0.18013711947910174,0.18004723447265697,0.17994774446454806,0.17983845960682815,0.17971968664165153,0.1795911826259099,0.17945341808768336,0.179306481586782,0.17915046344687616,0.17898677146701908,0.17881435683899397,0.17863648871373028,0.1784517092429574,0.1782639092462809,0.17807228709019027,0.1778813452113625,0.1776900458147875,0.17750686933287713,0.17732764905408513,0.17716274377123314,0.1770251491968274,0.17688260422825405,0.17677448505217996,0.1766951638082561,0.17665494575132853,0.17665494575132853,0.17665494575132853]
h150=[21.0,21.0,22.0,22.0,23.0,23.0,23.0,24.0,24.0,25.0,25.0,26.0,27.0,27.0,28.0,28.0,29.0,29.0,30.0,31.0,31.0,32.0,32.0,33.0,34.0,34.0,35.0,36.0,36.0,37.0,38.0,39.0,39.0,40.0,41.0,42.0,42.0,43.0,44.0,45.0,46.0,46.0,47.0,48.0,49.0,50.0,51.0,52.0,53.0,54.0,55.0,56.0,57.0,58.0,59.0,60.0,61.0,62.0,62.0,63.0,64.0,65.0,66.0,66.0,67.0,68.0,69.0,69.0,70.0,71.0,72.0,72.0,73.0,74.0,74.0,75.0,76.0,76.0,77.0,77.0,78.0,79.0,79.0,80.0,80.0,81.0,81.0,82.0,83.0,83.0,84.0,84.0,85.0,85.0,85.0,86.0,86.0,87.0,87.0]
f150=[0.1942207030746758,0.1942207030746758,0.19273236959038031,0.19273236959038031,0.19121960733854718,0.189057125720646,0.189057125720646,0.18726531649071296,0.18726531649071296,0.18540284067854412,0.18540284067854412,0.18344907925334633,0.18447647038782256,0.18447647038782256,0.18250704467361836,0.18250704467361836,0.18039176941899218,0.18039176941899218,0.17812126749269971,0.1794703538426652,0.1794703538426652,0.17705148565854031,0.17705148565854031,0.1784894506497919,0.17587036290460245,0.17587036290460245,0.17733071085648341,0.17447553930814882,0.17447553930814882,0.1758953892056348,0.17277664821104954,0.17410882490185786,0.17410882490185786,0.17071250574742994,0.17191438858132937,0.1730711117618249,0.1730711117618249,0.16928550050094296,0.1702670445486644,0.1711795717541084,0.17201660790425158,0.17201660790425158,0.17277093801818155,0.1683466483535743,0.1688865205857854,0.1693334128903819,0.16968418626947265,0.16993642690841174,0.17008846462533062,0.17013968421688672,0.1700884830419105,0.16993646358536701,0.169684241074196,0.16933348530887796,0.168886610142405,0.16834675439942492,0.17277107015692897,0.17201675600772237,0.17201675600772237,0.17117973451460355,0.17026722094966162,0.16928568881955303,0.17307132691809274,0.17307132691809274,0.17191461361139188,0.17071273862991712,0.17410908241341547,0.17410908241341547,0.17277691049127475,0.17589567426888572,0.1744758251442179,0.1744758251442179,0.17733101547941313,0.1758706661745573,0.1758706661745573,0.17848976872060093,0.17705179909235197,0.17705179909235197,0.17947067797414812,0.17947067797414812,0.1781215823911423,0.18039208912034865,0.18039208912034865,0.1825073653965842,0.1825073653965842,0.18447678686258204,0.18447678686258204,0.18344938370157965,0.1854031382369231,0.1854031382369231,0.1872656045665972,0.1872656045665972,0.18905740217144965,0.18905740217144965,0.19121987672491808,0.19273262202175498,0.19273262202175498,0.19422093839252985,0.19422093839252985]

plt.figure()
plt.plot(h30[4:len(h30)-4],smooth(f30,9)[4:len(h30)-4],label="30 degrees")
plt.plot(h90[4:len(h90)-4],smooth(f90,9)[4:len(h90)-4],label="90 degrees")
plt.plot(h150[4:len(h150)-4],smooth(f150,9)[4:len(h150)-4],label="150 degrees")
plt.legend()
plt.ylabel("Vapor concentration fraction")
plt.savefig("testdist.png")
plt.close()
def fitfunc3(x, a, b,c):
    return a-b*(-x+c)**-1
fit90=optimize.curve_fit(fitfunc3, fheightavg[1][90], fml[1][90], p0=[0,-0.2,100],maxfev=500000)[0]
curve90 = fitfunc3(20,fit90[0],fit90[1],fit90[2])
plt.figure(figsize=(5,4))
for fact in factors:
    for th in thetas:
        #plt.plot(fheightavg[fact][th]/fact,fml[fact][th]/fact/fact,label="Resolution: "+str(fact)+"x, Contact Angle: "+str(th))#,0])
        fit = optimize.curve_fit(fitfunc3, fheightavg[1][th], fml[1][th], p0=[0,-0.2,100],maxfev=500000,)[0]
        #plt.plot(dheightavg[90],fitfunc3(dheightavg[90],fit[0],fit[1],fit[2]),label="Contact Angle: "+str(th))
        plt.plot(fheightavg[fact][th]/fact,fml[fact][th]/fact/fact/fml[1][90],label="Contact Angle: "+str(th),ls="dashed")#,0])
        print(np.amax((((fitfunc3(20,fit[0],fit[1],fit[2])-curve90)/curve90))))
        print(fit[0],fit[1],fit[2],th)
"""
fit30=optimize.curve_fit(fitfunc3, fheightavg[1][30], fml[1][30], p0=[0.05,0.0017,80],maxfev=500000,bounds = ([0,0,80],[0.05,0.0017,100]))[0]
plt.plot(dheightavg[90],fitfunc3(dheightavg[90],fit30[0],fit30[1],fit30[2]),label="Contact Angle: "+str(30))
fit150=optimize.curve_fit(fitfunc3, fheightavg[1][150], fml[1][150], p0=[0.05,0.0017,80],maxfev=500000,bounds = ([0,0,80],[0.05,0.0017,100]))[0]
plt.plot(dheightavg[90],fitfunc3(dheightavg[90],fit150[0],fit150[1],fit150[2]),label="Contact Angle: "+str(150))
"""
#plt.xlim([1.05*max(i for v in fheightavg[1].values() for i in v),0.8*min(i for v in fheightavg[1].values() for i in v)])

"""plt.plot(np.array([80.10897435897436,79.49324324324324,78.4407894736842,77.45,77.25324675324676,76.27564102564102,75.48993288590604,75.19620253164557,74.34868421052632,73.54861111111111,73.2012987012987,72.34868421052632,71.67333333333333,71.17307692307692,70.1554054054054,70.07236842105263,69.18243243243244,68.55333333333333,68.10759493670886,67.38513513513513,67.01973684210526,66.1,65.55333333333333,64.99358974358974,64.18,63.955128205128204,63.11643835616438,62.74342105263158,62.0328947368421,61.42763157894737,60.967532467532465,60.479729729729726,59.955128205128204,59.08904109589041,58.99342105263158,58.060810810810814,57.74342105263158,57.0328947368421,56.577922077922075,55.993506493506494,55.55102040816327,54.967532467532465,54.38513513513514,53.98026315789474,53.358108108108105,52.94,52.43333333333333,51.88815789473684,51.439189189189186,50.824675324675326,50.363013698630134,49.824675324675326,49.3445945945946,48.908496732026144,48.358108108108105,47.98026315789474,47.358108108108105,46.9671052631579,46.50657894736842,45.967532467532465,45.69032258064516,44.953947368421055,44.81875,44.1,43.92948717948718]),np.array([0.0020101434804481637,0.001969507815555941,0.001963572669547593,0.00197445783271317,0.0019056179152960055,0.0018905698747088667,0.0019075208699599491,0.0018497736746207738,0.0018397720633752802,0.0018603946803995692,0.0017889650201258424,0.001795730750502422,0.0017916827829727528,0.0017574333241059865,0.0017591656200548688,0.001718066598691627,0.0017078562387163022,0.0017269310166759212,0.0016904797166410945,0.0016814566221235202,0.0016503195568116516,0.0016409032200900048,0.0016624300162091145,0.0016280375381803054,0.001622951356263875,0.0015861072011476513,0.0015872814579638206,0.0015571055507254713,0.0015583657134636106,0.0015773967436554783,0.0015417096515510723,0.0015499150492269296,0.0015186575787184933,0.0015221362958163587,0.0014895275255623255,0.0014933606158162875,0.001470275269278798,0.0014766916283139463,0.0014468809465434777,0.0014560850375324227,0.001420521775798775,0.001444409777838967,0.001444591853979257,0.0014242819118844272,0.0014402397888297986,0.001406538782383586,0.0014220394346268937,0.0013893380362525194,0.0014059777784483432,0.0013733155263087925,0.0013912061200396643,0.0013596386426444913,0.0013703959714737061,0.0013486247762694594,0.0013169820804185453,0.0013354235477109682,0.0013039723935119535,0.0013193772134565463,0.0012941518526441253,0.0013098595926332873,0.001286808600774511,0.0012997242852265847,0.001273235024967886,0.0012890818757213224,0.0012662783894562574]),label="1x Contact Angle: "+str(30))
plt.plot(np.array([78.5,77.5,76.86466165413533,76.31060606060606,75.5,74.78358208955224,74.20769230769231,73.5,72.79545454545455,72.3,71.5,70.85606060606061,70.5,69.5,68.97727272727273,68.5,67.73484848484848,67.3923076923077,66.5,65.96268656716418,65.5,64.78358208955224,64.5,63.5,63.09848484848485,62.5,61.93283582089552,61.5,60.82575757575758,60.5,59.62307692307692,59.5,58.5,58.21969696969697,57.5,57.06818181818182,56.5,56.00757575757576,55.5,54.946969696969695,54.5,53.93283582089552,53.5,52.916666666666664,52.5,51.916666666666664,51.5,50.93283582089552,50.5,49.946969696969695,49.5,48.97727272727273,48.5,48.03787878787879,47.5,47.13740458015267,46.5,46.404761904761905,45.5,45.5,44.67424242424242,44.5,43.82575757575758,43.5,42.916666666666664]),np.array([0.001916565394888632,0.0018988932387627466,0.0018747436489064627,0.0018599101336823796,0.0018338520811376042,0.00181624009012411,0.001800945573194718,0.0017811485969783214,0.001766211325402805,0.0017533062322236212,0.001734717752535893,0.001719559068867583,0.0016916455579696185,0.0016917543800958668,0.0016772062480378278,0.0016598779278227175,0.0016526157437795528,0.0016368029816795989,0.0016255839556347366,0.0016131712362248782,0.0015974543918912019,0.0015885314502546892,0.0015638305528848854,0.0015668432090640834,0.001555379966412947,0.0015445104437418829,0.001534344170384496,0.001520718645616329,0.0015135810932957257,0.001495360841414578,0.0014963030782728832,0.0014688456879871099,0.0014754240220934104,0.001467628063052854,0.0014572616321589182,0.0014495793643002177,0.0014394061625165745,0.001430778848173073,0.0014220657504124582,0.0014139615575741887,0.001405252128961308,0.0013981272630527075,0.0013891194479521828,0.0013808356221943923,0.001373741277910475,0.0013657989097442163,0.0013591105701809092,0.0013529288415938003,0.0013451473039200888,0.0013372679586055273,0.0013317126791545174,0.0013242526246654402,0.0013186298106614598,0.0013161440474112334,0.0013057617563781848,0.0012985130242527355,0.0012931321958971594,0.0012844820225835735,0.0012805527662862067,0.0012638884060403683,0.0012701533562142495,0.0012570297192154318,0.001256713922339274,0.0012490827045338372,0.0012451891181226257]),label="1x Contact Angle: "+str(90))
plt.plot(np.array([76.34810126582279,75.53125,74.83766233766234,74.31333333333333,73.21153846153847,72.59375,72.16447368421052,71.37820512820512,70.53125,70.00649350649351,69.37820512820512,68.59375,68.26973684210526,67.21875,66.66013071895425,66.31333333333333,65.35714285714286,64.83766233766234,64.27564102564102,63.72784810126582,63.23026315789474,62.375,61.83766233766234,61.24342105263158,60.75625,60.31333333333333,59.436708860759495,59.098684210526315,58.34415584415584,57.83766233766234,57.3051948051948,56.71974522292994,56.234177215189874,55.72784810126582,55.31333333333333,54.63461538461539,54.203947368421055,53.357142857142854,53.166666666666664,52.357142857142854,52.098684210526315,51.375,50.94078947368421,50.40384615384615,49.76973684210526,49.3525641025641,48.76973684210526,48.3525641025641,47.846666666666664,47.435064935064936,47.03333333333333,46.375,46.166666666666664,45.357142857142854,45.166666666666664,44.40909090909091,44.203947368421055,43.63461538461539,43.263513513513516,42.72784810126582,42.25,41.76875,41.27564102564103,40.69871794871795,40.26543209876543]),np.array([0.0019566800351382774,0.001931301984388222,0.0018993816596799974,0.0018816450916216629,0.0018560008697800943,0.0018432010358552242,0.0018195589946556393,0.00180463978251043,0.0017916356777578093,0.0017673858977675465,0.001753638238913207,0.0017454933048439815,0.0017227167341861956,0.0017056259864042466,0.001691898576987754,0.0016811712927168173,0.0016683545991136046,0.0016530503007765203,0.0016375341586637943,0.001631090710115604,0.0016086737899696987,0.0016036339154670677,0.0015901860843992453,0.0015747949829867135,0.0015684830554405144,0.0015591937339310216,0.0015452079655697683,0.0015321122332154246,0.0015211995956737328,0.0015121031179437167,0.0015013088280361668,0.0014946228834057202,0.001483417683229397,0.0014721650989234316,0.0014686331474438369,0.0014574485849733316,0.0014426264611962073,0.001436793212762439,0.001428691036843696,0.0014242974091816312,0.0014111499862649588,0.001403686730497485,0.00139374054506856,0.001388041316847095,0.0013792037545022622,0.0013729916068107066,0.0013635499940740745,0.0013581718693849692,0.0013496513100906327,0.0013437814543323127,0.0013343736711623777,0.001329372988577138,0.001321378361400167,0.0013184043617951431,0.0013076573462515,0.0013011422718415546,0.001294521510809609,0.0012922289127372298,0.0012859519974933351,0.0012747659245292801,0.001274860208897291,0.0012663115156836856,0.0012588352299381118,0.0012532211381190124,0.0012489257266965139]),label="1x Contact Angle: "+str(150))
plt.plot(np.array([159.66993464052288,158.21241830065358,157.28819444444446,156.75850340136054,155.85069444444446,155.09931506849315,154.6541095890411,153.7533783783784,153.09121621621622,152.69256756756758,151.65666666666667,150.98993288590603,150.55,149.61986301369862,148.85664335664336,148.35,147.4966216216216,146.87152777777777,146.263698630137,145.75503355704697,145.1858108108108,144.32993197278913,143.80872483221478,143.19256756756758,142.38013698630138,141.6768707482993,141.14802631578948,140.53040540540542,140.1510067114094,139.14930555555554,138.58965517241379,138.11744966442953,137.5,137.10855263157896,136.38486842105263,135.42123287671234,135.08503401360545,134.47278911564626,133.95,133.31,132.55369127516778,132.03716216216216,131.40753424657535,130.9046052631579,130.37828947368422,129.75986842105263,129.19127516778522,128.5758620689655,128.00335570469798,127.1858108108108,126.8841059602649,126.34121621621621,125.67880794701986,125.3476821192053,124.41156462585035,124.03960396039604,123.1734693877551,122.8758389261745,122.27364864864865,121.6217105263158,121.15771812080537,120.69463087248322,120.2687074829932,119.63175675675676,119.1241610738255,118.27181208053692,117.9496644295302,117.23469387755102,116.89869281045752,116.02413793103449,115.7312925170068,115.17465753424658,114.69463087248322,114.16780821917808,113.66554054054055,113.09121621621621,112.50337837837837,112.01700680272108,111.5728476821192,111.0641891891892,110.53973509933775,110.03666666666666,109.51324503311258,108.89666666666666,108.28523489932886,107.8125,107.24149659863946,106.89735099337749,106.19932432432432,105.85099337748345,105.0641891891892,104.85666666666667])/2,np.array([0.004081596000898668,0.003980521608340414,0.003899587801806342,0.0039410080368059265,0.0038240387164630843,0.0038104722935392225,0.003814835281115909,0.003790414769934954,0.003828825428131254,0.003740502084305284,0.0037399568133378563,0.003756443234973223,0.0036837210342987423,0.003684107823408906,0.0036440303136445616,0.0036492070084495776,0.003638973179631296,0.003598744180229738,0.0035993628185323904,0.0036106654631561466,0.003550469568096312,0.003553496665280626,0.003499921095382902,0.0035241465386719865,0.003511626187486861,0.003460281779400573,0.0034759783711107643,0.0034868310501829794,0.0034347917125400742,0.0034293115139993465,0.00339321051573943,0.003396754876788498,0.0032174497330362596,0.0033540925810726886,0.0033660851226044097,0.0033242396322643378,0.003331106894395256,0.003283200114160657,0.0032964984012089573,0.003310714466664944,0.0032575872338202858,0.0032707556172231,0.003229350251222035,0.00323604600716428,0.0031898495982211235,0.0032072463036132053,0.0032166859920868263,0.0031753814175406363,0.0031829361786521945,0.0031434288367072354,0.003154257066075739,0.003108709057861909,0.0031300374290563126,0.0031423949456562726,0.0031030956058483,0.003102088285997284,0.003068412403567757,0.003083506204414881,0.003048796471769996,0.003061892375932575,0.003010932211378355,0.003035589880403361,0.0029806251527113165,0.003009184028615464,0.0030271561018520273,0.002987756588094399,0.002995096799071043,0.0029662050227986323,0.0029790067996201238,0.0029366222114952318,0.002963166288443649,0.002912271420025697,0.0029383524643048714,0.002891931694198943,0.002918228547678161,0.0028691787281341745,0.00289691920155456,0.0028503566450343095,0.0028747900872305022,0.0028317448502222513,0.0028520310299561752,0.002812689785391318,0.0028370194585958107,0.002797050178374528,0.0028196981395309326,0.0027823255223034674,0.0028019268262340566,0.0027664594456107,0.0027885814601438354,0.0027520890530639544,0.0027724791056833215,0.0027353855053134637])/2,label="2x Contact Angle: "+str(30))
plt.plot(np.array([157.5,156.6060606060606,155.90943396226416,155.1780303030303,154.45283018867926,153.81132075471697,153.08712121212122,152.43984962406014,151.75378787878788,151.04166666666666,150.45283018867926,149.73584905660377,149.0340909090909,148.5,147.75849056603772,147.06439393939394,146.5,145.81439393939394,145.14772727272728,144.5,143.8901515151515,143.33969465648855,142.61278195488723,141.98872180451127,141.5,140.7781954887218,140.12830188679246,139.5,138.92045454545453,138.5,137.71804511278197,137.0719696969697,136.5,135.90530303030303,135.5,134.72348484848484,134.10227272727272,133.5,132.9436090225564,132.5,131.79323308270676,131.24809160305344,130.58712121212122,130.02651515151516,129.5,128.90530303030303,128.5,127.7781954887218,127.25378787878788,126.61742424242425,126.04887218045113,125.5,124.95075757575758,124.5,123.85984848484848,123.5,122.76315789473684,122.27862595419847,121.64772727272727,121.08712121212122,120.5,120.01136363636364,119.5,118.95075757575758,118.5,117.89015151515152,117.5,116.82954545454545,116.5,115.7781954887218,115.4423076923077,114.72348484848484,114.26335877862596,113.6780303030303,113.16287878787878,112.61509433962264,112.10227272727273,111.5,111.0719696969697,110.5,110.04887218045113,109.5,109.03383458646617,108.5,108.02651515151516,107.5,107.02651515151516,106.5,106.02651515151516,105.5,105.02651515151516,104.5])/2,np.array([0.003792485069236947,0.00376557160470397,0.003729432359711212,0.0037054897582287374,0.003683421993899841,0.003655658660022809,0.0036403513259391385,0.0036189766266599716,0.003595581129182025,0.0035826168954983598,0.0035617386583615476,0.0035429397283582914,0.003535806626593197,0.0035046327642737833,0.0034941469645788257,0.0034853864356004963,0.0034617467678994835,0.0034493816973045995,0.003437398601460159,0.003420069358217467,0.003402571748436032,0.003389321837094789,0.003374917656289239,0.0033631959498317374,0.0033385706360666177,0.003331457650045414,0.003322010489977044,0.003306961102680436,0.0032940484857307975,0.003262337437556004,0.0032641962684034015,0.003255032194066484,0.003240911403238954,0.0032296369479611897,0.0031984990373120237,0.003200697936737089,0.0031940100997188362,0.0031799982098012396,0.003167855111555646,0.003145935646663685,0.0031424667315961856,0.0031303899753484937,0.003123390131059356,0.0031122687422299247,0.003097460256805702,0.0030901116924226604,0.003064230165611052,0.0030643441530450216,0.0030519980815785505,0.00304195652685143,0.0030356980742407094,0.0030240363725470595,0.0030124067773704805,0.0029981506982694793,0.0029951270521186677,0.0029673900412741897,0.0029717609886668354,0.002959724246463165,0.0029512424777403466,0.002944171218369186,0.0029350655664590647,0.002924846181336054,0.0029157209572475643,0.002908330056563513,0.002893917924556031,0.0028895819606376306,0.0028707822126471274,0.0028682364340077756,0.002846514592123643,0.0028520323511985663,0.0028402280360443296,0.0028327402251334824,0.00282560761040236,0.0028174865886131054,0.0028156863304789007,0.0028022446307885644,0.002793377089608353,0.0027917291957759156,0.0027793534444520373,0.0027701232733771883,0.0027633196865226705,0.0027547949743791343,0.002747884552615735,0.002739478096871691,0.002730707453252987,0.0027242368369078066,0.002717665992764213,0.002709231099141791,0.002703005704865253,0.002694477891729036,0.002686303898675603,0.002679891069605528])/2,label="2x Contact Angle: "+str(90))
plt.plot(np.array([155.19281045751634,154.39144736842104,153.40255591054313,152.36217948717947,151.85576923076923,151.4950495049505,150.3603896103896,149.35534591194968,148.8202614379085,148.31935483870967,147.34516129032258,146.82903225806453,146.29738562091504,145.4056603773585,144.6487341772152,144.3476821192053,143.38387096774193,142.5891719745223,142.1979865771812,141.2753164556962,140.5891719745223,140.14144736842104,139.4108280254777,138.5891719745223,138.26333333333332,137.1867088607595,136.66772151898735,136.31045751633988,135.47151898734177,134.81612903225806,134.33653846153845,133.21974522292993,133.37171052631578,132.40822784810126,131.48417721518987,131.12091503267973,130.39171974522293,129.859477124183,129.27243589743588,128.59935897435898,128.10333333333332,127.34713375796179,126.86858974358974,126.31290322580645,125.60576923076923,125.05629139072848,124.3407643312102,124.18211920529801,123.11858974358974,122.48417721518987,122.21612903225807,121.21974522292993,121.19155844155844,120.33647798742139,119.84935897435898,119.2275641025641,118.5828025477707,118.12091503267973,117.16346153846153,117.13398692810458,116.2515923566879,116.27777777777777,115.25483870967741,114.55448717948718,114.31089743589743,113.43464052287581,113.1078431372549,112.16346153846153,112.02960526315789,111.34713375796179,111.19666666666667,110.19811320754717,110.31935483870967,109.3301282051282,108.63548387096775,108.24193548387096,107.390625,107.17628205128206,106.50974025974025,106.16558441558442,105.48064516129033,105.04276315789474,104.33647798742139,104.16558441558442,103.21698113207547,103.1590909090909,102.34713375796179,102.10855263157895,101.28525641025641,101.10264900662251,100.34713375796179,100.1470588235294])/2,np.array([0.00394449279189614,0.00391303252693742,0.0038765263140046772,0.0038395186272325217,0.003806736669939202,0.0037891268180896332,0.0037735561689572093,0.003744661384247152,0.0037151632382720116,0.003711075097023609,0.0036844642608353574,0.0036590519597087825,0.003660357954993183,0.003633926195507629,0.0036116879672898625,0.0035989183793233865,0.0035855627896492987,0.0035629639314019605,0.0035447207613817206,0.003533687254303863,0.003512123814771436,0.00349628890590263,0.0034869503561421916,0.0034685266198199952,0.0034516321055747005,0.0034363143631264874,0.0034210910965300927,0.003413035371250648,0.003388582601205859,0.0033683329722124476,0.0033629741761479703,0.003341232820584866,0.0033259307834035924,0.0033220247115257834,0.003309291755343104,0.0032929465221542475,0.0032818601484967615,0.003260440140561486,0.0032554147544155226,0.003238232825558074,0.003227416672555723,0.003217327975610092,0.0031927287985306003,0.0031934838483197695,0.003174453424326185,0.0031667549661751225,0.0031587209151676125,0.0031365482441175056,0.003135616614779962,0.00312292429483805,0.0031108210677700753,0.003094820330430037,0.00308418456553108,0.0030824388122026956,0.0030620744304450617,0.003056307519711477,0.0030438720579655373,0.0030361062695222893,0.003018544373593964,0.0030075663187828015,0.003004333561520949,0.0029844885094264537,0.002984261626330994,0.002977642628152268,0.002963576345656501,0.0029466558719589463,0.0029440673903149504,0.002931289226716878,0.002925161794481404,0.0029164248255189617,0.0029027139404853002,0.0029021832339709576,0.0028848640915452623,0.002881884371680647,0.002874198822742552,0.00286168692070396,0.0028504164563546808,0.002841977153123897,0.0028301915353920472,0.0028262472215335004,0.002813485333461976,0.002807615685079275,0.0027992365946054344,0.002789971842843815,0.002790321419836843,0.0027749280008577534,0.0027704151407626433,0.0027575808824552043,0.0027532141668574804,0.0027411347734510386,0.0027391661958169992,0.002726173241870135])/2,label="2x Contact Angle: "+str(150))
plt.plot(np.array([238.94155844155844,237.33630289532294,236.2873831775701,235.9579646017699,235.07432432432432,234.25681818181818,233.81627906976743,232.73076923076923,231.91822429906543,231.29908675799086,230.69954128440367,229.8732718894009,229.20094562647753,228.61036036036037,227.8490990990991,227.13720930232557,226.63004484304932,225.7186046511628,224.97477064220183,224.62217194570135,223.59490740740742,222.718009478673,222.35550458715596,221.77546296296296,220.88443396226415,220.60416666666666,219.69495412844037,219.33715596330276,218.38235294117646,217.47453703703704,217.21428571428572,216.59318181818182,215.7476851851852,215.21330275229357,214.69655172413792,213.94819819819818,213.14350797266516,212.45560747663552,212.02045454545456,211.44170403587444,211.12045454545455,210.02522935779817,209.2523148148148,208.97085201793723,208.33710407239818,207.86805555555554,207.1027397260274,206.2206896551724,205.85227272727272,205.204128440367,204.71461187214612,204.0829596412556,203.24654377880185,202.9609195402299,202.162100456621,201.72706422018348,201.12276785714286,200.4004524886878,199.9284064665127,198.9322429906542,198.70871559633028,198.042600896861,197.57623318385652,196.99321266968326,196.07110091743118])/3,np.array([0.006052388237388816,0.005943818258493665,0.005862495097506974,0.005829086721670471,0.0057284670940725875,0.005792270494664331,0.0057400130402240615,0.005747701565804319,0.005638071240837787,0.0056367054452610196,0.0056310514767126145,0.005614714408396606,0.0055498219151223225,0.0055681453420894185,0.00555654759061901,0.005587336057117117,0.00551318710287097,0.005501781796540694,0.005445295595699915,0.005465924264445117,0.0054542447706031096,0.005397990694523861,0.005424417486067322,0.00542766160905962,0.005360693946396297,0.005368187948752266,0.005402125299661485,0.0053219271482226,0.0053220426187357605,0.005266507701054051,0.0052786588041100115,0.005289895994831181,0.005229270955033763,0.005243399450779177,0.005164726506347905,0.005195878751640195,0.005195543521977119,0.005140276839837374,0.005155677415849227,0.005181055203908317,0.005119295802348821,0.0051193870210789385,0.005070380994100202,0.005084271077265608,0.005030486632606129,0.005041179472282564,0.0050429694568350746,0.005009318464341138,0.005021670989872029,0.004962880722887304,0.004979007428225751,0.004996460915120617,0.004940883014720851,0.004953957268177423,0.004900835760693877,0.004912534373810455,0.0049368738472077185,0.004882172523710941,0.004889913759142435,0.004840132475283747,0.004862660903760452,0.004801991029340581,0.004820113085035787,0.004845642488290494,0.004802281586187425])/3,label="3x Contact Angle: "+str(30))
plt.plot(np.array([236.5,235.60353535353536,234.92131979695432,234.10353535353536,233.43401015228426,232.7550505050505,232.08333333333334,231.32997481108313,230.63888888888889,230.02272727272728,229.2348484848485,228.58333333333334,227.98232323232324,227.18434343434345,226.55527638190955,225.96733668341707,225.19444444444446,224.55583756345177,223.97222222222223,223.20854271356785,222.58838383838383,222.00252525252526,221.2286432160804,220.65736040609136,220.05778894472363,219.31565656565655,218.7550505050505,218.1281407035176,217.4318181818182,216.89698492462313,216.21464646464648,215.5580808080808,214.99748743718592,214.28535353535352,213.7348484848485,213.11363636363637,212.4238578680203,211.91161616161617,211.24874371859298,210.63888888888889,210.05778894472363,209.3391959798995,208.87121212121212,208.21859296482413,207.61675126903555,207.03768844221105,206.35929648241205,205.86683417085428,205.2286432160804,204.64720812182742,204.06313131313132,203.39646464646464,202.90151515151516,202.28535353535352,201.74683544303798,201.1180904522613,200.47848101265822,199.97738693467338,199.38790931989925,198.83668341708542,198.2247474747475,197.69597989949747,197.08333333333334,196.4570707070707,195.95969773299748])/3,np.array([0.0056452426446632785,0.005611911104864525,0.005577610755841604,0.005537343849508287,0.005512482349958429,0.005488515623676816,0.005456027636852745,0.005438468394039372,0.005418673922444886,0.005391316582214499,0.005378097260780927,0.005359793338803735,0.005332847151176855,0.005324426463643504,0.005305776830627906,0.005277196873193162,0.0052736920352773,0.005253814693612503,0.005227283672435727,0.005223440918851192,0.0052046742636401785,0.005180243976642528,0.005174948744555776,0.005159429755847772,0.005135289996625707,0.0051313797696757105,0.005112376074608813,0.005092059667167751,0.005079984410132267,0.005055526578874339,0.005050167541952284,0.0050369423899202116,0.005015493773463768,0.0050105050993870665,0.004997527173152812,0.004976276976745577,0.004969361712767362,0.004948424252453146,0.004938849600835546,0.004926765744202741,0.004907225577602557,0.004902858644317571,0.004876622111417049,0.0048713635244670254,0.0048611300790915996,0.0048422600973608825,0.004838378716139686,0.004812738648663195,0.004807460184304525,0.004797622746044148,0.004780615785064928,0.004776221322100148,0.004752122883634088,0.004745285045460682,0.004726027275646202,0.004721290315560685,0.004717401112104482,0.00469494151353857,0.004687378897795689,0.00466935387386089,0.004663042861739914,0.00464462024130061,0.004643301426964674,0.004635843879939399,0.004615855377741816])/3,label="3x Contact Angle: "+str(90))
plt.plot(np.array([233.1351931330472,232.13771186440678,231.20464135021098,230.98922413793105,230.13829787234042,229.11181434599155,228.34110169491527,228.13695652173914,227.1645299145299,226.26150627615064,225.5255319148936,224.9625550660793,224.06932773109244,223.5214592274678,223.23390557939913,222.0363247863248,221.5168776371308,220.49789915966386,220.02813852813853,219.0699152542373,218.64439655172413,218.03333333333333,217.0609243697479,216.34978540772534,216.16960352422907,215.1276150627615,214.4785407725322,214.13377192982455,213.0483193277311,212.57296137339057,212.13377192982455,211.33958333333334,210.56437768240343,209.9069264069264,209.16596638655463,208.4636752136752,207.97173913043477,207.23126338329766,206.46008403361344,206.07974137931035,205.5,204.9636752136752,204.13713080168776,203.60991379310346,203.12715517241378,202.08958333333334,201.32553191489362,200.86909871244634,200.34689507494647,199.9636752136752,199.27330508474577,198.59442060085837,198.0194805194805,197.37076271186442,196.59533898305085,195.98504273504273,195.72435897435898,195.0926724137931,194.3886554621849,193.41489361702128,192.741452991453,192.69574468085105,192.10897435897436,191.3886554621849,190.52350427350427])/3,np.array([0.00590188656658659,0.005841897160537602,0.005807560192412338,0.005769491208211345,0.005739206252216643,0.005722572464487128,0.005686412778035406,0.0056496605333558915,0.005636362809450443,0.005625752526791934,0.005600191732910781,0.005576052865550846,0.005555620934137162,0.005532039702197623,0.005507519785400079,0.005486499539007502,0.005485981039084816,0.005468771106274706,0.005442264757388487,0.005431885031870871,0.005404818640841514,0.005385833316235176,0.005368741069577952,0.0053478927156334595,0.00532726596221569,0.005314383741436591,0.0052952117016579125,0.0052758780295351915,0.005267291550363631,0.005247721258915395,0.005226730130662149,0.005217973331509328,0.005200854092059287,0.005183128495268153,0.005167606884115635,0.005155077170875487,0.005145981164635049,0.005143921428327916,0.005117450426774838,0.0051025384410151664,0.005089207163642449,0.005069721571281341,0.005055507571035174,0.005040276892677522,0.005024480474984507,0.005009671463454395,0.004994357166464681,0.004985925222760681,0.004976132544343849,0.004959117493390527,0.004942369228720871,0.00492473579560349,0.004910844996970725,0.004899786351826471,0.004896882013079541,0.004882325703159155,0.0048681475049870415,0.004849411020832213,0.004837063520531266,0.004827970436506929,0.004815509787787792,0.004805065710670457,0.004791011753913391,0.004775307820464191,0.004769175677709095])/3,label="3x Contact Angle: "+str(150))
plt.xlim([85,35])"""
plt.xlabel("Liquid Height/Post Height (%)")
plt.ylabel("dV/dt")
plt.tight_layout()
#plt.ylim([0.00135,0.00135*4.5])
#plt.xticks(ticks = [60,50,40,30,20], labels = [round(60/74,2)*100,round(50/74,2)*100,round(40/74,2)*100,round(30/74,2)*100,round(20/74,2)*100])
plt.legend()
plt.savefig("cangs_w"+str(wid)+"_h"+str(hei)+".png", dpi=200, format='png')

def fitfunc2(x, a, b,c):
    return a-b*(-x+c)**0.5


plt.figure()
for th in thetas:
    plt.plot(fvol[1][th],label=th)
plt.xlabel("t")
plt.ylabel("V")# per unit area")
#plt.xlim([0,62])
#plt.title("Height Measured as Average Height")
#plt.gca().invert_xaxis()
plt.legend()
plt.tight_layout()
plt.savefig("volvst.PNG",dpi=200,format='png')

fit30=optimize.curve_fit(fitfunc2, fheightavg[1][30], fml[1][30], p0=[0.05,0.0017,80],maxfev=500000,bounds = ([0,0,80],[0.05,0.0017,100]))[0]
#print(fit30, "HEEEEEEEEEEEEEEEEEEEEEEEEERE")
fit90=optimize.curve_fit(fitfunc2, fheightavg[1][90], fml[1][90], p0=[0.05,0.0017,80],maxfev=500000,bounds = ([0,0,80],[0.05,0.0017,100]))[0]
#fit150=optimize.curve_fit(fitfunc2, dheightavg[150], dml[150], p0=[0,0.000079008422,75],maxfev=500000)[0]
#fit30=optimize.curve_fit(fitfunc2, dheightavg[30], dml[30], p0=[0,0.000079008422],maxfev=500000)[0]
#fit150=optimize.curve_fit(fitfunc2, dheightavg[150], dml[150], p0=[0,0.000079008422],maxfev=500000)[0]
plt.figure()
#for th in thetas:
    #print("Here",dheightavg[th])
    #plt.plot(dheightavg[th],(dml[th]-dml[90])/(dml[90]),label="Contact Angle: "+str(th))
    #plt.plot(dheightavg[th],dml[th],label="Contact Angle: "+str(th))
plt.plot(dheightavg[90],fitfunc2(dheightavg[90],fit30[0],fit30[1],fit30[2]),label="Contact Angle: "+str(30))
plt.plot(dheightavg[90],fitfunc2(dheightavg[90],fit90[0],fit90[1],fit90[2]),label="Contact Angle: "+str(30))
#plt.plot(dheightavg[90],(dml[90]/2-dml[90]/2)/(dml[90]/2),label="Contact Angle: "+str(90))
#plt.plot(fheightavg[1][90],(fitfunc2(fheightavg[1][90],fit30[0],fit30[1],fit30[2])),label="Contact Angle: "+str(30))
#plt.plot(fheightavg[1][th]/1,fml[1][th]/1,label="Resolution: "+str(1)+"x, Contact Angle: "+str(th))#,0])
#plt.plot(fheightavg[1][90],(fitfunc2(fheightavg[1][90],fit90[0],fit90[1],fit90[2])),label="Contact Angle: "+str(90))
#plt.plot(dheightavg[90],(fitfunc2(dheightavg[90],fit150[0],fit150[1],fit150[2],fit150[3])/2-dml[90]/2)/(dml[90]/2),label="Contact Angle: "+str(90))
plt.xlabel("Height")
plt.ylabel("dV/dt")# per unit area")
plt.title("Height Measured as Average Height")
plt.gca().invert_xaxis()
plt.legend()
#plt.ylim([0.06,0.12])
plt.xlim([80,0])
plt.tight_layout()
#plt.savefig("ContactLine_fluxvsheightNOTNORMALISED.PNG",dpi=200,format='png')
plt.savefig("0.0002.PNG",dpi=200,format='png')


#print("!!!!!!!",fit)
plt.figure()
#for th in thetas:
#    print("Here",dheightavg[th])
#    #plt.plot(dheightavg[th],(dml[th]/2-dml[90]/2)/(dml[90]/2),label="Contact Angle: "+str(th))
#    plt.plot(dheightavg[th],dml[th]/2,label="Contact Angle: "+str(th))
#plt.plot(dheightavg[30],dml[30]/2,label="Contact Angle: "+str(30))
heightz = np.append(dheightavg[90],[31,28,25,22,19,16,13,10,7,4,1])
plt.plot(heightz,fitfunc2(heightz,fit30[0],fit30[1],fit30[2])/2,label="TEST: "+str(30))
plt.plot(heightz,fitfunc2(heightz,fit90[0],fit90[1],fit90[2])/2,label="TEST: "+str(90))
plt.plot(heightz,fitfunc2(heightz,fit150[0],fit150[1],fit150[2])/2,label="TEST: "+str(150))
#plt.plot(dheightavg[90],dml[90]/2,label="Contact Angle: "+str(90))
#plt.plot(dheightavg[150],dml[150]/2,label="Contact Angle: "+str(150))
plt.xlabel("Height")
plt.ylabel("dV/dt")# per unit area")
plt.xlim([0,62])
plt.title("Height Measured as Average Height")
plt.gca().invert_xaxis()
plt.legend()
plt.tight_layout()
#plt.savefig("ContactLine_fluxvsheightNOTNORMALISED.PNG",dpi=200,format='png')
plt.savefig("mlvsh.PNG",dpi=200,format='png')

plt.figure()
#for th in thetas:
#    print("Here",dheightavg[th])
#    plt.plot(dheightavg[th],(dml[th]/2-dml[90]/2)/(dml[90]/2),label="Contact Angle: "+str(th))
#    #plt.plot(dheightavg[th],dml[th]/2,label="Contact Angle: "+str(th))
#heightz = np.append(dheightavg[90],[31,28,25,22,19,16,13,10,7,4,1])
plt.plot(heightz,(fitfunc2(heightz,fit30[0],fit30[1],fit30[2])/2-fitfunc2(heightz,fit90[0],fit90[1],fit90[2])/2)/(fitfunc2(heightz,fit90[0],fit90[1],fit90[2])/2),label="Contact Angle: "+str(30))
plt.plot(heightz,(fitfunc2(heightz,fit90[0],fit90[1],fit90[2])/2-fitfunc2(heightz,fit90[0],fit90[1],fit90[2])/2)/(fitfunc2(heightz,fit90[0],fit90[1],fit90[2])/2),label="Contact Angle: "+str(90))
plt.plot(heightz,(fitfunc2(heightz,fit150[0],fit150[1],fit150[2])/2-fitfunc2(heightz,fit90[0],fit90[1],fit90[2])/2)/(fitfunc2(heightz,fit90[0],fit90[1],fit90[2])/2),label="Contact Angle: "+str(150))

plt.xlabel("Height")
plt.ylabel("dV/dt")# per unit area")
plt.title("Height Measured as Average Height")
plt.gca().invert_xaxis()
plt.legend()
plt.tight_layout()
#plt.savefig("ContactLine_fluxvsheightNOTNORMALISED.PNG",dpi=200,format='png')
plt.savefig("diff.PNG",dpi=200,format='png')


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
#print(h0, k)
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
#print((height[0]-height[-1])/(height_an[0]-height_an[-1]))
plt.legend()
plt.savefig("massloss.png", dpi=500, format='png')

