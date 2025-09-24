import os
import numpy as np
import random
import shutil
import sys
import fileinput

owd = os.getcwd() # Directory containing the python file
datadir="data" # All data directories start with the "data" folder
os.system("mkdir "+datadir) # Make the "data" folder (if it doesn't exist already)
os.chdir(datadir) # Enter it

#params3={"lx":np.array([300]),"ly":np.array([101]),"lz":np.array([101]),"timesteps":np.array([450000]),"evaporationrate":np.array([0.00003]),"tau1":np.array([0.502]),"tau3":np.array([0.502]),"s01":np.array([0.015]),"s02":np.array([0.015]),"s12":np.array([0.02464])}



#params3={"timesteps":np.array([2000000]),"s12":np.array([0.0015,0.00176047,0.00201303,0.00225,0.00246418,0.00264907,0.00279904,0.00290954,0.00297721]),"channeltheta":np.array([30])}#,0.00225])}
params3={"lx":np.array([200,245,283,316,346]),
         "ly":np.array([200,245,283,316,346]),
         "ncomp":np.array([2,3,4,5,6]),
         "prop1":np.array([0.5,0.3333334,0.25,0.2,0.16666667]),
         "prop2":np.array([0.5,0.3333333,0.25,0.2,0.16666667]),
         "prop3":np.array([0.0,0.3333333,0.25,0.2,0.16666667]),
         "prop4":np.array([0.0,0.0,0.25,0.2,0.16666667]),
         "prop5":np.array([0.0,0.0,0.0,0.2,0.16666666]),
         "prop6":np.array([0.0,0.0,0.0,0.0,0.16666666]),
         "s12":np.array([0.005,0.005,0.005,0.005,0.005])}#,0.00225])}

params3={"lx":np.array([245,245,245,245,245,283,283,283,283,283,346,346,346,346,346]),
         "ly":np.array([245,245,245,245,245,283,283,283,283,283,346,346,346,346,346]),
         "ncomp":np.array([3,3,3,3,3,4,4,4,4,4,6,6,6,6,6]),
         "prop1":np.array([0.1,0.2,0.5,0.7,0.55,0.1,0.1,0.1,0.1,0.2,0.1,0.05,0.5,0.05,0.05]),
         "prop2":np.array([0.45,0.4,0.25,0.15,0.3,0.1,0.3,0.1,0.2,0.2,0.18,0.19,0.1,0.1,0.05]),
         "prop3":np.array([0.45,0.4,0.25,0.15,0.15,0.4,0.3,0.1,0.3,0.2,0.18,0.19,0.1,0.15,0.05]),
         "prop4":np.array([0.0,0.0,0.0,0.0,0.0,0.4,0.3,0.7,0.5,0.4,0.18,0.19,0.1,0.2,0.85/3.0]),
         "prop5":np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.18,0.19,0.1,0.25,0.85/3.0]),
         "prop6":np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.18,0.19,0.1,0.25,0.85/3.0]),
         "s12":np.array([0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005])}

params3={"lx":np.array([245*2,283*2,316*2,346*2]),#,374*2]),
         "ly":np.array([245*2,283*2,316*2,346*2]),#,374*2]),
         "ncomp":np.array([3,4,5,6]),#,7]),
         "prop1":np.array([0.3333334,0.25,0.2,0.16666667]),#,0.14285714285]),
         "prop2":np.array([0.3333333,0.25,0.2,0.16666667]),#,0.14285714285]),
         "prop3":np.array([0.3333333,0.25,0.2,0.16666667]),#,0.14285714285]),
         "prop4":np.array([0.0,0.25,0.2,0.16666667]),#,0.14285714285]),
         "prop5":np.array([0.0,0.0,0.2,0.16666666]),#,0.14285714285]),
         "prop6":np.array([0.0,0.0,0.0,0.16666666]),#,0.14285714285]),
         "prop7":np.array([0.0,0.0,0.0,0.0]),#,0.14285714290]),
         #"s12":np.array([0.01,0.01,0.01,0.01,0.01,0.01])}#,0.00225])}
         "s12":np.array([0.005,0.005,0.005,0.005,0.005])}#,0.0025])}#,0.00225])}

params3={"lx":np.array([100]),#,374*2]),
         "ly":np.array([100]),#,374*2]),
         "ncomp":np.array([3]),#,7]),
         "prop1":np.array([0.3333334]),#,0.14285714285]),
         "prop2":np.array([0.3333333]),#,0.14285714285]),
         "prop3":np.array([0.3333333]),#,0.14285714285]),
         "prop4":np.array([0.0]),#,0.14285714285]),
         "prop5":np.array([0.0]),#,0.14285714285]),
         "prop6":np.array([0.0]),#,0.14285714285]),
         "prop7":np.array([0.0]),#,0.14285714290]),
         #"s12":np.array([0.01,0.01,0.01,0.01,0.01,0.01])}#,0.00225])}
         "s12":np.array([0.005])}

params3={"lx":np.array([139,147,155,162]),
         "ly":np.array([139,147,155,162]),
         "lz":np.array([139,147,155,162]),
         "ncomp":np.array([5,6,7,8]),
         "prop1":np.array([0.2,0.16666667,0.14285714285,0.0125]),
         "prop2":np.array([0.2,0.16666667,0.14285714285,0.0125]),
         "prop3":np.array([0.2,0.16666667,0.14285714285,0.0125]),
         "prop4":np.array([0.2,0.16666667,0.14285714285,0.0125]),
         "prop5":np.array([0.2,0.1666666,0.14285714285,0.0125]),
         "prop6":np.array([0.0,0.1666666,0.14285714285,0.0125]),
         "prop7":np.array([0.0,0.0,0.14285714290,0.0125]),
         "prop8":np.array([0.0,0.0,0,0.0125]),
         #"s12":np.array([0.01,0.01,0.01,0.01,0.01,0.01])}#,0.00225])}
         "surfacetension":np.array([0.005,0.005,0.005,0.005])}#,0.00225])}

params3={"lx":np.array([155,162]),
         "ly":np.array([155,162]),
         "lz":np.array([155,162]),
         "ncomp":np.array([7,8]),
         "prop1":np.array([0.14285714285,0.0125]),
         "prop2":np.array([0.14285714285,0.0125]),
         "prop3":np.array([0.14285714285,0.0125]),
         "prop4":np.array([0.14285714285,0.0125]),
         "prop5":np.array([0.14285714285,0.0125]),
         "prop6":np.array([0.14285714285,0.0125]),
         "prop7":np.array([0.14285714290,0.0125]),
         "prop8":np.array([0,0.0125]),
         #"s12":np.array([0.01,0.01,0.01,0.01,0.01,0.01])}#,0.00225])}
         "surfacetension":np.array([0.005,0.005])}

params3={"lx":np.array([186]),
         "ly":np.array([186]),
         "lz":np.array([186]),
         "ncomp":np.array([12]),
         "prop1":np.array([0.08333333333 ]),
         "prop2":np.array([0.08333333333 ]),
         "prop3":np.array([0.08333333333 ]),
         "prop4":np.array([0.08333333333 ]),
         "prop5":np.array([0.08333333333 ]),
         "prop6":np.array([0.08333333333 ]),
         "prop7":np.array([0.08333333333 ]),
         "prop8":np.array([0.08333333333 ]),
         "prop9":np.array([0.08333333334 ]),
         "prop10":np.array([0.08333333334 ]),
         "prop11":np.array([0.08333333334 ]),
         "prop12":np.array([0.08333333334 ]),
         #"s12":np.array([0.01,0.01,0.01,0.01,0.01,0.01])}#,0.00225])}
         "surfacetension":np.array([0.005])}

params3={"lx":np.array([int(117*3),int(135*3),int(151*3),int(165*3)]),
         "ly":np.array([int(117*3),int(135*3),int(151*3),int(165*3)]),
         "lz":np.array([25,25,25,25]),
         "ncomp":np.array([3,4,5,6]),
         "prop1":np.array([0.33333333,0.25,0.2,0.16666667]),
         "prop2":np.array([0.33333333,0.25,0.2,0.16666667]),
         "prop3":np.array([0.33333334,0.25,0.2,0.16666667]),
         "prop4":np.array([0.0,0.25,0.2,0.16666667]),
         "prop5":np.array([0.0,0.0,0.2,0.16666667]),
         "prop6":np.array([0.0,0.0,0.0,0.16666667]),
         "prop7":np.array([0.0,0.0,0.0,0.0]),
         "prop8":np.array([0.0,0.0,0.0,0.0]),
         #"s12":np.array([0.01,0.01,0.01,0.01,0.01,0.01])}#,0.00225])}
         "surfacetension":np.array([0.005,0.005,0.005,0.005])}#,0.00225])}


"""params3={"lx":np.array([int(117*1.5),int(135*1.5),int(151*1.5),int(165*1.5)]),
         "ly":np.array([int(117*1.5),int(135*1.5),int(151*1.5),int(165*1.5)]),
         "lz":np.array([75,75,75,75]),
         "ncomp":np.array([3,4,5,6]),
         "prop1":np.array([0.33333333,0.25,0.2,0.16666667]),
         "prop2":np.array([0.33333333,0.25,0.2,0.16666667]),
         "prop3":np.array([0.33333334,0.25,0.2,0.16666667]),
         "prop4":np.array([0.0,0.25,0.2,0.16666667]),
         "prop5":np.array([0.0,0.0,0.2,0.16666667]),
         "prop6":np.array([0.0,0.0,0.0,0.16666667]),
         "prop7":np.array([0.0,0.0,0.0,0.0]),
         "prop8":np.array([0.0,0.0,0.0,0.0]),
         #"s12":np.array([0.01,0.01,0.01,0.01,0.01,0.01])}#,0.00225])}
         "surfacetension":np.array([0.005,0.005,0.005,0.005])}"""
"""
params3={"lx":np.array([int(117*1.5),int(135*1.5),int(151*1.5),int(165*1.5)]),
         "ly":np.array([int(117*1.5),int(135*1.5),int(151*1.5),int(165*1.5)]),
         "lz":np.array([5,5,5,5]),
         "ncomp":np.array([3,4,5,6]),
         "prop1":np.array([0.33333333,0.25,0.2,0.16666667]),
         "prop2":np.array([0.33333333,0.25,0.2,0.16666667]),
         "prop3":np.array([0.33333334,0.25,0.2,0.16666667]),
         "prop4":np.array([0.0,0.25,0.2,0.16666667]),
         "prop5":np.array([0.0,0.0,0.2,0.16666667]),
         "prop6":np.array([0.0,0.0,0.0,0.16666667]),
         "prop7":np.array([0.0,0.0,0.0,0.0]),
         "prop8":np.array([0.0,0.0,0.0,0.0]),
         #"s12":np.array([0.01,0.01,0.01,0.01,0.01,0.01])}#,0.00225])}
         "surfacetension":np.array([0.005,0.005,0.005,0.005])}"""

"""params3={"lx":np.array([102,117,129]),
         "ly":np.array([102,117,129]),
         "lz":np.array([102,117,129]),
         "ncomp":np.array([2,3,4]),
         "prop1":np.array([0.5,0.33333334,0.25]),
         "prop2":np.array([0.5,0.33333333,0.25]),
         "prop3":np.array([0.0,0.33333333,0.25]),
         "prop4":np.array([0.0,0.0,0.25]),
         "prop5":np.array([0.0,0.0,0.0]),
         "prop6":np.array([0.0,0.0,0.0]),
         "prop7":np.array([0.0,0.0,0.0]),
         "prop8":np.array([0.0,0.0,0.0]),
         #"s12":np.array([0.01,0.01,0.01,0.01,0.01,0.01])}#,0.00225])}
         "surfacetension":np.array([0.005,0.005,0.005])}"""#,0.00225])}

"""params3={"lx":np.array([155]),
         "ly":np.array([155]),
         "lz":np.array([155]),
         "ncomp":np.array([7]),
         "prop1":np.array([0.14285714285]),
         "prop2":np.array([0.14285714285]),
         "prop3":np.array([0.14285714285]),
         "prop4":np.array([0.14285714285]),
         "prop5":np.array([0.14285714285]),
         "prop6":np.array([0.14285714285]),
         "prop7":np.array([0.14285714290]),
         "prop8":np.array([0.0]),
         #"s12":np.array([0.01,0.01,0.01,0.01,0.01,0.01])}#,0.00225])}
         "surfacetension":np.array([0.005])}#,0.00225])}"""
"""
params3={"lx":np.array([130]),
         "ly":np.array([130]),
         "lz":np.array([130]),
         "ncomp":np.array([7]),
         "prop1":np.array([0.14285714285]),
         "prop2":np.array([0.14285714285]),
         "prop3":np.array([0.14285714285]),
         "prop4":np.array([0.14285714285]),
         "prop5":np.array([0.14285714285]),
         "prop6":np.array([0.14285714285]),
         "prop7":np.array([0.14285714290]),
         #"s12":np.array([0.01,0.01,0.01,0.01,0.01,0.01])}#,0.00225])}
         "surfacetension":np.array([0.01])}"""

"""params3={"lx":np.array([245,283,316,346]),#,374*2]),
         "ly":np.array([245,283,316,346]),#,374*2]),
         "ncomp":np.array([3,4,5,6]),#,7]),
         "prop1":np.array([0.3333334,0.25,0.2,0.16666667]),#,0.14285714285]),
         "prop2":np.array([0.3333333,0.25,0.2,0.16666667]),#,0.14285714285]),
         "prop3":np.array([0.3333333,0.25,0.2,0.16666667]),#,0.14285714285]),
         "prop4":np.array([0.0,0.25,0.2,0.16666667]),#,0.14285714285]),
         "prop5":np.array([0.0,0.0,0.2,0.16666666]),#,0.14285714285]),
         "prop6":np.array([0.0,0.0,0.0,0.16666666]),#,0.14285714285]),
         "prop7":np.array([0.0,0.0,0.0,0.0]),#,0.14285714290]),
         #"s12":np.array([0.01,0.01,0.01,0.01,0.01,0.01])}#,0.00225])}
         "s12":np.array([0.005,0.005,0.005,0.005,0.005])}"""

"""params3={"lx":np.array([500]),
         "ly":np.array([500]),
         "ncomp":np.array([4]),
         "prop1":np.array([0.25]),
         "prop2":np.array([0.25]),
         "prop3":np.array([0.25]),
         "prop4":np.array([0.25]),
         "prop5":np.array([0.0]),
         "prop6":np.array([0.0]),
         "prop7":np.array([0.0]),
         "s12":np.array([0.005])}"""

"""params3={"lx":np.array([374]),
         "ly":np.array([374]),
         "ncomp":np.array([7]),
         "prop1":np.array([0.14285714285]),
         "prop2":np.array([0.14285714285]),
         "prop3":np.array([0.14285714285]),
         "prop4":np.array([0.14285714285]),
         "prop5":np.array([0.14285714285]),
         "prop6":np.array([0.14285714285]),
         "prop7":np.array([0.14285714285]),
         "s12":np.array([0.005])}"""

#params3={"lx":np.array([228]),"ly":np.array([200]),"postwidth":np.array([162]),"offsety":np.array([-17]),"theta":np.array([30])}
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
#loop_rec(params3,list(params3),len(params3),[],owd) # Run the loop
loop(params3,list(params3),len(params3[list(params3)[0]]),owd)
