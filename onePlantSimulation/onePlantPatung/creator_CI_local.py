#!/usr/bin/env python3
# -*- coding: utf-8
"""
Created on Mon Nov 23 12:54:23 2020

@author: emilio

This scripts will make a full bunch of data composed of one line ofr eeach combination of 
parameter

With this script I will do 2 tests, the sensibility test with al the parameters
And then, the specific test fixing the parameters that are not important to me




"""

import numpy as np
import pickle
import os
from os import path




"""
MEGA LOOP
"""

# I put the conditions 

#valoresOptimo = np.arange(0, 200, 20)

#aquiLosDatosParaAnalsisi de sensibilidad. 


valoresR = np.linspace(0.009, 0.013,1)
valoresB1 = np.linspace(0.03,0.04,1)
valoresB2 =np.linspace(0.03,0.04,1)
valoresA = np.linspace(0.65, 0.65,1)
valoresG = np.linspace(0.0152, 0.2,10) # para ver más variación (nno porque sea rea)
valoresMU = np.linspace(0.2, 0.2,1)
valoresI0 = np.linspace(0.001, 0.1, 10) #variable


linea = []


n=0

for r in valoresR:
    for b1 in valoresB1:
        for b2 in valoresB2:
            for a in valoresA:
                for mu in valoresMU:
                    for I0 in valoresI0:
                        for g in valoresG:
                            
                            linea.append(r)
                            linea.append(b1)
                            linea.append(b2)
                            linea.append(a)
                            linea.append(mu)
                            linea.append(g)
                            linea.append(I0)
    
                            liga = '/srv/home/emilio/sim_local/condicionesIniciales'
                            os.makedirs(liga, exist_ok= True)
    
                            with open(path.join(liga, "CI_%s_%s_%s_%s_%s_%s_%s.txt" %(r, b1, b2, a, mu, g, I0)), "wb") as output:
                                pickle.dump(linea, output)
                                         #output.write(lista)
                            linea= []



