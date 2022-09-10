#!/usr/bin/env python3
# -*- coding: utf-8 

"""
Here we have all the necessary functions without instantiation. This script should be saved in .nwg for the condor 9.0.11 structure.


0. Import the used libraries
"""

from operator import add
import numpy as np
import seaborn as sns; sns.set()

"""
I. GENERAL DYNAMICS dinaMig

This function runs the general dynamics. It uses other functions enlisted below
 
It receives: 
    
a) necessary conditions for ODE to run (dic_ODE)


It exits: 
Generates a 2-dimensional matrix. One of the entries is a vector of 1*3 entries: (steps)* (s, i, x) values

It calls:
euler_red() to call the ODE in each of the cell of the lattice 
"""


def creadorMatriz_sim(dic_ODE): 

    pasos = dic_ODE["steps_"]
    matrix_zvalues = np.zeros((pasos,len(dic_ODE["cond_Ini"])))
    matrix_zvalues[0,:]= dic_ODE["cond_Ini"]    
    for n in range(pasos-1):
        #después llama a la función Euler
        z1= euler_red(matrix_zvalues[n], dic_ODE["para_"], dic_ODE["size_step"], dic_ODE["mode_"])
        matrix_zvalues[n+1]= z1
    return(matrix_zvalues)
        
"""
II. euler_red

This function creates the euler integration using the differential (  Z_t+h = f(s,i,p)*h + Z_t)

They receive:
    a) (s,i,x) of the n-1 step for one specific coordinate
    b) the parameters for the ode
    c) the size step
    d) the mode of the growth function (mono or log)
    
They exit: a list of (s,i, x) for the n step

They call:
    a) the SIX_eu_normal to calculate the differential
"""

def euler_red(z_pre, para_, size_step, mode_): #z es para referir el vector sip que le metes
    h= size_step
    suma= [i * h for i in SIX_eu_normal(z_pre, para_, mode_)]#aqui hago el calcylo
    znu= list(map(add,z_pre, suma))     # esto es la z más la suma de la diferencial pero por si quiero ver varios pasos en algun momento
    return(znu) #renglon


"""
III. six_eu_normal
This functions creates the differential

It receives:
    a) (s,i,x) of the n-1 step for one specific coordinate
    b) the parameters for the ode
    c) the mode of the growth function (mono or log)
    
It exits: a list of (s,i, x) differential 

It does not call any other function
"""
def SIX_eu_normal(z_pre, para_, mode_):
   #parametros del modelo
    
    para= para_
    z= z_pre
    
    r= para[0]
    b1=para[1]
    b2 = para[2]
    g= para[3]
    a= para[4]
    mu= para[5]

    s= z[0]
    i= z[1]
    p= z[2]

    
    #ecuaciones
    if str(mode_)== "log":
        dsdt = r*s*(1-(s))-b1*s*p-b2*s*i  
    
    elif str(mode_) == "mono":
        dsdt = r*(1-(s))-b1*s*p-b2*s*i
    
    else: 
        print ("ERROR DE MODO")
    
    didt = b1*s*p+b2*s*i-g*i
    dpdt= a*i -mu*p 
    dzdt = [dsdt, didt, dpdt]    
    return dzdt

