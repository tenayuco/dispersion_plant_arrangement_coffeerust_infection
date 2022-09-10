#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This code creates a file (initialConditions) with all the textfiles that correspond to each of the initial conditions. In Condor language, this file will
be used to create the submit_file


@author: emilio
"""

import numpy as np
import pandas as pd
import pickle
import os
from os import path

#!/usr/bin/env python3
# -*- coding: utf-8 -*-


""" 
1. Initial conditions
Here we define the different values of the parameters

1.1. For the general figures we will fix  r,b1, b2, a and mu values. 
"""

valoresR = [0.011]
valoresB1 = [0.035]
valoresB2 = [0.035]
valoresA = [0.65]
valoresMU = [0.2]
sizeEuler = [0.01]

"""
Then we define the different values for m, g, io and the patterns
"""

valoresM= [0.001, 0.005, 0.01,0.05,0.1]
valoresG = [0.015, 0.056]
valoresIni = [0.001, 0.1]
randomList = list(np.repeat("random",30))
valoresPatrones = ["spaced", "aggregated", "rows"]
valoresPatrones = valoresPatrones + randomList

"""
1.2 If we want to do the sensitivy analysis we ran our initial procedure, but changing the value
of the other parameters at once.
"""

#valoresR = [0.009, 0.011, 0.013]
#valoresB1 = [0.03, 0.035, 0.04]
#valoresB2 = [0.03, 0.035, 0.04]
#valoresA = [0.1, 0.65, 1.2]
#valoresMU = [0.1, 0.2, 0.3]


"""
1.3 If we want to do the euler analysis we  we ran our initial procedure, but changing the value
of the step size.
"""

#sizeEuler = [0.01, 0.02]

"""
Here we run the general loop
"""


linea = []

patungDirectory = '/srv/home/emilio/sim_fullLattice'  #this changes depending the computer


sim = 0

for g in valoresG:
    for m in valoresM:
        for ini in valoresIni:
            for pat in valoresPatrones:
                
                if pat == "random":
                    sim = sim +1
                    
                else:
                    sim = 1
                
                for r in valoresR:
                    for b1 in valoresB1:
                        for b2 in valoresB2:
                            for a in valoresA:
                                for mu in valoresMU:
                                    for size in sizeEuler:
                                        linea.append(g)
                                        linea.append(m)
                                        linea.append(ini)
                                        linea.append(pat)
                                        linea.append(sim)
                                        linea.append(r)
                                        linea.append(b1)
                                        linea.append(b2)
                                        linea.append(a)
                                        linea.append(mu)
                                        linea.append(size)
                
                        
                                        liga = patungDirectory+ "/initialConditions"
                                        os.makedirs(liga, exist_ok= True)
                        
                                        with open(path.join(liga, "CI_%s_%s_%s_%s_%s_%s_%s_%s_%s_%s_%s.txt" %(g,m, ini, pat, sim, r, b1, b2, a, mu, size)), "wb") as output:
                                            pickle.dump(linea, output)
                                 #output.write(lista)
                                        linea= []

