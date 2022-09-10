#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
Here we have all the necessary functions to run specific examples (for Fig.3. One plant simulation examples)



0. Import the used libraries
"""

from operator import add
import numpy as np
import seaborn as sns; sns.set()
import pandas as pd

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




"""
AQUI PARA CORRER TODO
"""
###script para correr la funcion una vez



"""
Primero para poder graficar las soluciones en el tiempo
"""

 #el chiste es controlar los parametros como si fuera un año
##########CARACRTETRISTICAS GENERALES


tiempoReal= 300 #dias
sizeStep = 0.01 #para la integracion
steps= int(tiempoReal/sizeStep) # esto es para que calcule los pasos necesario para llegar al tiempo Adim
mode= "mono"

b=0.011  
b1=0.035  
b2=0.035
a= 0.65
g= 0.056   
l= 100
mu= 0.2  
    
condIni = [1,0,0.01]

parametros =[b, b1,b2, g, a, mu]

dic_ODE = {"cond_Ini": condIni, "para_":parametros, "steps_": steps ,"size_step": sizeStep, "mode_": mode}

ro= ((b2*mu)+(b1*a))/(mu*g) 
print ("ro=",ro)
guardar= "FALSE"


"""
This codes generates some specific examples for the simulation as save it as examplesIndividualPlant
"""

dfMatrizSimTotal= pd.DataFrame(columns = ['S','I','X', "I0", "g"])
for Ires in [0.001, 0.01, 0.1]:
    for g in [0.015, 0.056, 0.12]:  #este cambio lo hice post_co. el 16 de agosto del 2022
         parametros =[b, b1,b2, g, a, mu]
         condIni = [1,0,Ires]
         dic_ODE = {"cond_Ini": condIni, "para_":parametros, "steps_": steps ,"size_step": sizeStep, "mode_": mode}
         matrizSim = creadorMatriz_sim(dic_ODE)
         dfMatrizSim = pd.DataFrame(matrizSim, columns = ['S','I','X'])
         dfMatrizSim["I0"]=Ires
         dfMatrizSim["g"]=g
         
         dfMatrizSimTotal= pd.concat([dfMatrizSimTotal, dfMatrizSim])


dfMatrizSimTotal.to_csv(" ./codes_MoraVanCauwelaert_2022/FiguresCodes/onePlant_data/examplesIndividualPlant.csv")      #this paths depends on where you put your general folder    
    
    

#
