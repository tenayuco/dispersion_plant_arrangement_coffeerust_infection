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
Here we fix the different values of the parameters
"""


actMigra = 1
valoresM= [0.001, 0.005, 0.01,0.05,0.1]
valoresG = [0.015, 0.056]
valoresIni = [0.001, 0.01, 0.1]
randomList = list(np.repeat("random",30))
valoresPatrones = ["spaced", "aggregated", "rows"]
valoresPatrones = valoresPatrones + randomList

"""
Here we run the loop
"""


linea = []


sim = 0

for g in valoresG:
    for m in valoresM:
        for ini in valoresIni:
            for pat in valoresPatrones:
                
                if pat == "random":
                    sim = sim +1
                    
                else:
                    sim = 1
                        
                linea.append(g)
                linea.append(m)
                linea.append(ini)
                linea.append(pat)
                linea.append(sim)
                
                        
                liga = path.join("patung_superComputer_directory", "/initialConditions")
                os.makedirs(liga, exist_ok= True)
                        
                with open(path.join(liga, "CI_%s_%s_%s_%s_%s.txt" %(g,m, ini, pat, sim)), "wb") as output:
                    pickle.dump(linea, output)
                                 #output.write(lista)
                linea= []

