import os
import numpy as np
import random
import shutil
import sys
import fileinput

owd = os.getcwd() # Directory containing the python file
datadir="data/inflowsimple7" # All data directories start with the "data" folder
os.system("mkdir "+datadir) # Make the "data" folder (if it doesn't exist already)
os.chdir(datadir) # Enter it

#params3={"inflowmomentum":np.array([-0.0002,-0.0002,-0.0002,-0.0002,-0.0002,-0.00004,-0.00004,-0.00004,-0.00004,-0.00004]),"lx":np.array([228,228,228,228,228,228,228,228,228,228]),"ly":np.array([200,200,200,200,200,200,200,200,200,200]),"postwidth":np.array([162,162,162,162,162,162,162,162,162,162]),"offsety":np.array([-17,-17,-17,-17,-17,-17,-17,-17,-17,-17]),"theta":np.array([30,60,90,120,150,30,60,90,120,150])}
#params3={"inflowmomentum":np.array([-0.0002,-0.0002,-0.0002,-0.0002,-0.0002]),"lx":np.array([228,228,228,228,228]),"ly":np.array([200,200,200,200,200]),"postwidth":np.array([162,162,162,162,162]),"offsety":np.array([-17,-17,-17,-17,-17]),"theta":np.array([30,60,90,120,150])}


#params3={"diffusivity":np.array([0.0008,0.0008,0.0008,0.0008,0.0008,0.0008,0.0008,0.0008,0.0008,0.0008]),"inflowmomentum":np.array([-0.0002,-0.00002,-0.0002,-0.00002,-0.0002,-0.00002,-0.0002,-0.00002,-0.0002,-0.00002]),"lx":np.array([228,228,228,228,228,228,228,228,228,228]),"ly":np.array([200,200,200,200,200,200,200,200,200,200]),"postwidth":np.array([162,162,162,162,162,162,162,162,162,162]),"offsety":np.array([-17,-17,-17,-17,-17,-17,-17,-17,-17,-17]),"theta":np.array([30,30,60,60,90,90,120,120,150,150])}
#params3={"diffusivity":np.array([0.0008,0.0008,0.0008,0.0008,0.0008]),"inflowmomentum":np.array([-0.0002,-0.0002,-0.0002,-0.0002,-0.0002,]),"lx":np.array([228,228,228,228,228]),"ly":np.array([200,200,200,200,200]),"postwidth":np.array([162,162,162,162,162]),"offsety":np.array([-12,-12,-12,-12,-12]),"theta":np.array([30,90,150,60,120])}
#params3={"diffusivity":np.array([0.0008]),"inflowmomentum":np.array([-0.0002]),"lx":np.array([228,228,228,228,228,228,228,228,228,228]),"ly":np.array([200]),"postwidth":np.array([162]),"offsety":np.array([-17]),"theta":np.array([30])}

params3={"diffusivity":np.array([0.0008,0.008,0.08]),"inflowmomentum":np.array([-0.0002]),"lx":np.array([228]),"ly":np.array([200]),"postwidth":np.array([162]),"offsety":np.array([-12]),"theta":np.array([30,90,150,60,120]),"surfacetension":np.array([0.00125,0.0005,0.0001,0.00005])}
params3={"diffusivity":np.array([0.0008,0.008,0.08]),"inflowmomentum":np.array([-0.0002]),"lx":np.array([228]),"ly":np.array([200]),"postwidth":np.array([162]),"offsety":np.array([-12]),"theta":np.array([30,90,150,60,120]),"surfacetension":np.array([0.00001])}
params3={"diffusivity":np.array([0.024]),"inflowmomentum":np.array([-0.0006]),"lx":np.array([228*3-72*3]),"ly":np.array([200*2]),"postwidth":np.array([162*3]),"offsety":np.array([-12*3]),"theta":np.array([30,90,150]),"surfacetension":np.array([0.0001])}
#params3={"diffusivity":np.array([0.012]),"inflowmomentum":np.array([-0.0006,-0.0006,-0.0006,-0.0003,-0.0003,-0.0003,-0.0002,-0.0002,-0.0002,-0.00015,-0.00015,-0.00015]),"lx":np.array([228*1-72*1,228*1-72*1,228*1-72*1,228*2-72*2,228*2-72*2,228*2-72*2,228*3-72*3,228*3-72*3,228*3-72*3,228*4-72*4,228*4-72*4,228*4-72*4]),"ly":np.array([133*1,133*2,133*3,133*4]),"postwidth":np.array([162*1,162*2,162*3,162*4]),"offsety":np.array([-12*1,12*2,12*3,12*4]),"theta":np.array([30,90,150,30,90,150,30,90,150,30,90,150]),"surfacetension":np.array([0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001])}
#params3={"diffusivity":np.array([0.004,0.004,0.004,0.008,0.008,0.008,0.012,0.012,0.012,0.016,0.016,0.016]),"inflowmomentum":np.array([-0.0006,-0.0006,-0.0006,-0.0006,-0.0006,-0.0006,-0.0006,-0.0006,-0.0006,-0.0006,-0.0006,-0.0006]),"lx":np.array([228*1-72*1,228*1-72*1,228*1-72*1,228*2-72*2,228*2-72*2,228*2-72*2,228*3-72*3,228*3-72*3,228*3-72*3,228*4-72*4,228*4-72*4,228*4-72*4]),"ly":np.array([133*1,133*1,133*1,133*2,133*2,133*2,133*3,133*3,133*3,133*4,133*4,133*4]),"postwidth":np.array([162*1,162*1,162*1,162*2,162*2,162*2,162*3,162*3,162*3,162*4,162*4,162*4]),"offsety":np.array([-12*1,-12*1,-12*1,-12*2,-12*2,-12*2,-12*3,-12*3,-12*3,-12*4,-12*4,-12*4]),"theta":np.array([30,90,150,30,90,150,30,90,150,30,90,150]),"surfacetension":np.array([0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001])}
factor=1
params3={"diffusivity":np.array([0.004*factor]),"inflowmomentum":np.array([-0.0005]),"lx":np.array([228*factor-72*factor]),"ly":np.array([133*factor]),"postwidth":np.array([162*factor-72*factor]),"offsety":np.array([-12*factor]),"theta":np.array([30]),"surfacetension":np.array([0.0001])}
#factor=1
#params3={"diffusivity":np.array([0.004*factor]),"inflowmomentum":np.array([-0.0005]),"lx":np.array([228*factor-72*factor]),"ly":np.array([133*factor]),"postwidth":np.array([162*factor-72*factor]),"offsety":np.array([-12*factor]),"theta":np.array([150]),"surfacetension":np.array([0.0001])}


#params3={"diffusivity":np.array([0.004*factor]),"inflowmomentum":np.array([-0.001]),"lx":np.array([228*factor-72*factor]),"ly":np.array([133*factor]),"postwidth":np.array([162*factor-72*factor]),"offsety":np.array([-12*factor]),"theta":np.array([90]),"surfacetension":np.array([0.0001])}
#params3={"inflowmomentum":np.array([0,0,0,0,0]),"lx":np.array([228,228,228,228,228]),"ly":np.array([200,200,200,200,200]),"postwidth":np.array([162,162,162,162,162]),"offsety":np.array([-17,-17,-17,-17,-17]),"theta":np.array([30,60,90,120,150])}
#params3={"inflowmomentum":np.array([-0.0002]),"lx":np.array([228]),"ly":np.array([200]),"postwidth":np.array([162]),"offsety":np.array([-17]),"theta":np.array([150])}
#params3={"lx":np.array([228,228]),"ly":np.array([200,200]),"postwidth":np.array([162,162]),"offsety":np.array([-17,-17]),"theta":np.array([30,90])}


##################################################################################
# Recursive loop function. Pass it a dictionary d of strings and arrays, with keys
# y and length n. It will run a simulation for EVERY combination of parameters.
##################################################################################
def loop_rec(d,y,n,li:list,dirr):
    if n >= 1: # n will reduce with each loop until it gets to zero after we have reached the last parameter
        new=li
        for x in d[y[n-1]]: # Iterate through values in the array stored with the key y[n-1]
            os.system("mkdir "+y[n-1]+"_"+str(x)) # Make a new directory formatted as "{parameter name}_{parameter_value}"
            os.chdir(y[n-1]+"_"+str(x)) # Enter it
            new.append(x) # Add parameter x from this loop
            print(new)
            loop_rec(d,y, n - 1,new,dirr) # Run a new loop with n now reduced by 1 and a list "new" containing the parameter "x" from this loop
            new.remove(x) # Reset new for the next loop
            os.chdir("..") # Exit the directory we just created so we can repeat this process
    else:
        
        curdir= os.getcwd() # Save the directory we came from so we can go back to it, as we are within a for loop
        os.chdir(dirr) # Go to directory containing this file
        datapath="" # Empty string
        for i in range(len(li)):
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
        



#loop_rec(params1,list(params1),len(params1),[])
loop_rec(params3,list(params3),len(params3),[],owd) # Run the loop
#loop(params3,list(params3),len(params3[list(params3)[0]]),owd)