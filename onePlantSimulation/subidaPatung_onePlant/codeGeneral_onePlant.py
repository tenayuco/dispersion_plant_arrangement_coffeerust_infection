#!/usr/bin/env python3
# -*- coding: utf-8
"""
This codes instantiates the functions defined 
in __init_py (kept in the .nwg file) 
using the different values defined in the submitfile
"""


import numpy as np
import argparse
from os import path
import pickle
import pandas as pd

""" 
We first import the values of the initial conditions in the submitfile
"""


parser = argparse.ArgumentParser(description = "arguments that have to be read by the code")
parser.add_argument('--code', help = 'initial matrix', required= True)
parser.add_argument('--subdir', help= 'subdirectory', required= True)
args = parser.parse_args()




with open(path.join(args.subdir, args.code), "rb") as out:
        condiciones = pickle.load(out)


from nwg import *  # here we import the __init__.py script


"""
Fixed parameters
"""

time= 300 #days
sizeStep = 0.01 #para la integracion
steps= int(time/sizeStep) # esto es para que calcule los pasos necesario para llegar al tiempo Adim
mode= "mono"


"""
Here all the parameters of the ODE equation can be modifyed. For Figure 3 we only modifyed the g and io
"""

r= condiciones[0]
b1= condiciones[1]
b2= condiciones[2]
a= condiciones[3]
mu= condiciones[4]
g= condiciones[5]
I0= condiciones[6]
condIni = [1,0,I0]
parameters =[r, b1,b2, g, a, mu]
dic_ODE = {"cond_Ini": condIni, "para_":parameters, "steps_": steps ,"size_step": sizeStep, "mode_": mode}

"""
General dynamics
"""

matrizSim = creadorMatriz_sim(dic_ODE)

liga =  "/srv/home/emilio/sim_onePlant"+ "/salida/matriz_onePlant"

with open(path.join(liga, "matrizLocal_%s" %(args.code)), "wb") as output:
        pickle.dump(matrizSim, output)



ro= ((b2*mu)+(b1*a))/(mu*g) 


"""
Here we run some procedure to create the data frame
"""
    

ano2_Ini= int(300/sizeStep-1)

Imax1 = max(matrizSim[0:ano2_Ini, 1])
Ieq= max(matrizSim[ano2_Ini:ano2_Ini+1, 1])

Idif= Imax1-Ieq
#aqui se podria poner e valor tambien, pero da igual.   
pasosImax1= np.where(matrizSim[0:ano2_Ini,1]==Imax1)[0][0]
tiempoAnio1 = sizeStep*pasosImax1

dicLocal= {"r":[r],"b1":[b1], "b2":[b2], "a":[a],"mu":[mu], "g":[g], "R0":[ro], "I0":[I0], "Imax1":[Imax1],"Imax1-Ieq":[Idif],  "tiempoImax1":[tiempoAnio1]}


df_local_resumen = pd.DataFrame(data=dicLocal)

#df_temporal_resumen.append(dict(zip(df_temporal_resumen.columns,[r,b1, b2, a,mu, g, ro, I0, Imax1,Idif,  tiempoAnio1])), ignore_index=True)


liga2 =  "/srv/home/emilio/sim_onePlant"+ "/salida/dfOnePlant"  #patung directory

#uso este y no un pickle porque la salida ya la uqiero poder usar en r (en pycle te deja todo en python)
df_local_resumen.to_csv(str(liga2)+ "dfLocalresumen_%s.csv"%(args.code), index = False)            
            
   