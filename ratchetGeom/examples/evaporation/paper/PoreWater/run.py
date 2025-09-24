import os
import numpy as np
import random
import shutil
import sys
import fileinput
import math
import time

owd = os.getcwd() # Directory containing the python file
datadir="data" # All data directories start with the "data" folder
os.system("mkdir "+datadir) # Make the "data" folder (if it doesn't exist already)
os.chdir(datadir) # Enter it
factor = 1

params={"factor":np.array([factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor]),"diffusivity":np.array([0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141]),"lx":np.array([78,108,138,168,198,78,108,138,168,198,78,108,138,168,198,78,108,138,168,198,78,108,138,168,198,78,108,138,168,198,78,108,138,168,198,78,108,138,168,198,78,108,138,168,198,78,108,138,168,198,78,108,138,168,198,78,108,138,168,198,78,108,138,168,198,78,108,138,168,198,78,108,138,168,198]),"ly":np.array([195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195]),"postwidth":np.array([12,42,72,102,132,12,42,72,102,132,12,42,72,102,132,12,42,72,102,132,12,42,72,102,132,12,42,72,102,132,12,42,72,102,132,12,42,72,102,132,12,42,72,102,132,12,42,72,102,132,12,42,72,102,132,12,42,72,102,132,12,42,72,102,132,12,42,72,102,132,12,42,72,102,132]),"offsety":np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]),"theta":np.array([30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,90,90,90,90,90,90,90,90,90,90,90,90,90,90,90,90,90,90,90,90,90,90,90,90,90,150,150,150,150,150,150,150,150,150,150,150,150,150,150,150,150,150,150,150,150,150,150,150,150,150]),"surfacetension":np.array([0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01]),"postheight":np.array([76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76]),"outflowoffset":np.array([105,105,105,105,105,90,90,90,90,90,60,60,60,60,60,30,30,30,30,30,0,0,0,0,0,105,105,105,105,105,90,90,90,90,90,60,60,60,60,60,30,30,30,30,30,0,0,0,0,0,105,105,105,105,105,90,90,90,90,90,60,60,60,60,60,30,30,30,30,30,0,0,0,0,0])*factor}
params={"factor":np.array([factor,factor,factor]),"diffusivity":np.array([0.02141,0.02141,0.02141]),"lx":np.array([78,78,78]),"ly":np.array([195,195,195]),"postwidth":np.array([12,12,12]),"offsety":np.array([-30,-30,-30]),"theta":np.array([30,90,150]),"surfacetension":np.array([0.01,0.01,0.01]),"postheight":np.array([76,76,76]),"outflowoffset":np.array([120,120,120])*factor}
params={"factor":np.array([factor]),"diffusivity":np.array([0.16666667]),"lx":np.array([78])*2,"ly":np.array([185])*2,"postwidth":np.array([12])*2,"offsety":np.array([-15]),"theta":np.array([30]),"surfacetension":np.array([0.00111]),"postheight":np.array([66])*2,"outflowoffset":np.array([95])*2}
#params={"factor":np.array([factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor]),"diffusivity":np.array([0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141]),"lx":np.array([78,108,138,168,198,78,108,138,168,198,78,108,138,168,198]),"ly":np.array([195,195,195,195,195,195,195,195,195,195,195,195,195,195,195]),"postwidth":np.array([12,42,72,102,132,12,42,72,102,132,12,42,72,102,132]),"offsety":np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]),"theta":np.array([30,30,30,30,30,30,90,90,90,90,90,150,150,150,150,150]),"surfacetension":np.array([0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01]),"postheight":np.array([76,76,76,76,76,76,76,76,76,76,76,76,76,76,76]),"outflowoffset":np.array([180,180,180,180,180,180,180,180,180,180,180,180,180,180,180])}
#wids=np.array([108,78,108,138,78,198,78,168,198,198,78,168,198,168,78,108,138,168,198,78,108,138,168,198,78,168,78,108,138,168,198,138,168,198,108,168,78,138,78,78])
#spac=wids-66
#params={"factor":np.array([factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor]),"diffusivity":np.array([0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141]),"lx":wids,"ly":np.array([195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195]),"postwidth":spac,"offsety":np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]),"theta":np.array([30,30,30,30,30,30,30,30,30,30,30,30,30,30,90,90,90,90,90,90,90,90,90,90,90,90,90,90,90,90,90,90,90,90,90,150,150,150,150,150]),"surfacetension":np.array([0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01]),"postheight":np.array([76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76]),"outflowoffset":np.array([0,30,30,30,60,60,90,90,90,110,180,180,180,220,0,0,0,0,0,30,30,30,30,30,60,60,90,90,90,90,90,110,110,110,220,30,90,90,180,220])}

"""
params={"factor":np.array([factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor,factor]),"diffusivity":np.array([0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141,0.02141]),"lx":np.array([78,108,138,168,198,138,168,198,78,108,138,168,198,78,108,138,168,198,78,108,138,168,198,78,108,138,168,198,78,108,138,168,198,78,108,138,168,198,78,108,138,168,198,78,108,138,168,198,78,108,138,168,198,78,108,138,168,198]),"ly":np.array([95,95,95,95,95,195,195,195,95,95,95,95,95,105,105,105,105,105,135,135,135,135,135,165,165,165,165,165,195,195,195,195,195,95,95,95,95,95,105,105,105,105,105,135,135,135,135,135,165,165,165,165,165,195,195,195,195,195]),"postwidth":np.array([12,42,72,102,132,72,102,132,12,42,72,102,132,12,42,72,102,132,12,42,72,102,132,12,42,72,102,132,12,42,72,102,132,12,42,72,102,132,12,42,72,102,132,12,42,72,102,132,12,42,72,102,132,12,42,72,102,132]),"offsety":np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]),"theta":np.array([30,30,30,30,30,30,30,30,90,90,90,90,90,90,90,90,90,90,90,90,90,90,90,90,90,90,90,90,90,90,90,90,90,150,150,150,150,150,150,150,150,150,150,150,150,150,150,150,150,150,150,150,150,150,150,150,150,150]),"surfacetension":np.array([0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01]),"postheight":np.array([76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76,76])}
"""
#params={"factor":np.array([factor]),"diffusivity":np.array([0.02141]),"lx":np.array([78]),"ly":np.array([195]),"postwidth":np.array([12]),"offsety":np.array([-5]),"theta":np.array([30]),"surfacetension":np.array([0.01]),"postheight":np.array([76]),"outflowoffset":np.array([0])}
#params={"factor":np.array([factor]),"diffusivity":np.array([0.02141,0.0002141]),"lx":np.array([108]),"ly":np.array([200]),"postwidth":np.array([42]),"offsety":np.array([0]),"theta":np.array([5,30,60,90]),"surfacetension":np.array([0.01]),"postheight":np.array([96]),"equilibriumtimesteps":np.array([10000])}
#params={"factor":np.array([factor]),"diffusivity":np.array([0.02141]),"lx":np.array([148]),"ly":np.array([106]),"postwidth":np.array([82]),"offsety":np.array([-17]),"theta":np.array([150]),"surfacetension":np.array([0.01]),"postheight":np.array([46])}
#params={"factor":np.array([factor]),"diffusivity":np.array([0.02141]),"lx":np.array([148]),"ly":np.array([120]),"postwidth":np.array([82]),"offsety":np.array([-17]),"theta":np.array([30,90,150]),"surfacetension":np.array([0.01]),"postheight":np.array([46])}

##################################################################################
# Recursive loop function. Pass it a dictionary d of strings and arrays, with keys
# y and length n. It will run a simulation for EVERY combination of parameters.
##################################################################################
def loop_rec(d,y,n,li:list,dirr):
    if n >= 1: # n will reduce with each loop until it gets to zero after we have reached the last parameter
        new=li
        for x in d[y[n-1]]: # Iterate through values in the array stored with the key y[n-1]
            os.system("mkdir "+y[n-1]+"_"+str(x)) # Make a new directory formatted as "{parameter name}_{parameter_value}"
            #print(y[n-1]+"_"+str(x))
            os.chdir(y[n-1]+"_"+str(x)) # Enter it
            #if (y[n-1]=="ly"):
            #    print(x)
            new.append(x) # Add parameter x from this loop
            #print(new,n)
            #print(new)
            loop_rec(d,y, n - 1,new,dirr) # Run a new loop with n now reduced by 1 and a list "new" containing the parameter "x" from this loop
            #print(new)
            new.pop() # Reset new for the next loop
            #print(new)
            os.chdir("..") # Exit the directory we just created so we can repeat this process
    else:
        
        curdir= os.getcwd() # Save the directory we came from so we can go back to it, as we are within a for loop
        os.chdir(dirr) # Go to directory containing this file
        datapath="" # Empty string
        for i in range(len(li)):
            #print(y[len(li)-1-i]+"_"+str(li[i])+"/")
            datapath+=y[len(li)-1-i]+"_"+str(li[i])+"/" # Construct the data path in the same format as we created the directories (this could be done way better)
        seed=random.randint(100000, 999999) # Random seed to make simulation specific input file
        shutil.copy("input.txt","input/input"+str(seed)+".txt") # Create simulation specific input file

        for line in fileinput.input(["input/input"+str(seed)+".txt"], inplace=True): # Find the line corresponding to the data directory and change it
            if line.strip().startswith("datadir"):
                    line = "datadir=\""+datadir+"/"+datapath+"\"\n"
            sys.stdout.write(line)

        for i in range(len(li)):
            
            print(li)
            for line in fileinput.input(["input/input"+str(seed)+".txt"], inplace=True): # Change each parameter to their desired value
                if line.strip().startswith(y[len(li)-1-i]+"="):
                    line = y[len(li)-1-i]+"="+str(li[i])+"\n"
                sys.stdout.write(line)

        for line in fileinput.input(["submitscP"], inplace=True): # Set the seed in the slurm file so that it can be passed to the simulation (so it knows what input file to read)
            if line.strip().startswith("seed"):
                    line = "seed=\""+str(seed)+"\"\n"
            sys.stdout.write(line)
        os.system("sbatch submitscP") # Run submit script
        os.chdir(curdir) # Go back to original directory

##################################################################################
# Normal loop function. Pass it a dictionary d of strings and arrays, with keys y 
# and length n. It will run a simulation for for parameters in the order they are
# placed in the arrays. For example, if I want to change the radius for the values
# [a,b,c] and the contact angle with [d,e,f] it will run simulations with
# {radius:a, contact angle:d} then {radius:b, contact angle:e} then
# {radius:c, contact angle:f}.
##################################################################################
def loop(d,y,n,dirr):
    for i in range(n):
        os.chdir(dirr+"/"+datadir) # Start in data directory
        li=[] # Reset list
        for key in y: # Iterate through parameters
            temp=key
            os.system("mkdir "+temp+"_"+str(d[key][i])) # Make directory corresponding to the parameter and its value at the index i
            os.chdir(temp+"_"+str(d[key][i])) # Enter it
            li.append(d[key][i]) # Add to list
        
        os.chdir(dirr) # Pretty much the same as above
        datapath=""
        for j in range(len(li)):
            datapath+=y[j]+"_"+str(li[j])+"/"
        print(datapath) # Print the path to the data directory
        seed=random.randint(100000, 999999)
        shutil.copy("input.txt","input/input"+str(seed)+".txt")

        for line in fileinput.input(["input/input"+str(seed)+".txt"], inplace=True):
            if line.strip().startswith("datadir"):
                    line = "datadir=\""+datadir+"/"+datapath+"\"\n"
            sys.stdout.write(line)

        for i in range(len(li)):
            print(li)
            for line in fileinput.input(["input/input"+str(seed)+".txt"], inplace=True):
                if line.strip().startswith(y[i]+"="):
                    if isinstance(li[i], str):
                        line = y[i]+"=\""+str(li[i])+"\"\n"
                    else:
                        line = y[i]+"="+str(li[i])+"\n"
                sys.stdout.write(line)

        for line in fileinput.input(["submitscP"], inplace=True):
            if line.strip().startswith("seed"):
                    line = "seed=\""+str(seed)+"\"\n"
            sys.stdout.write(line)
        os.system("sbatch submitscP")
        new=[]
        #time.sleep(0.5)
        



loop_rec(params,list(params),len(params),[],owd) # Run the loop
#loop(params,list(params),len(params[list(params)[0]]),owd)