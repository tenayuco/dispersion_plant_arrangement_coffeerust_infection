#!/usr/bin/env python3
# -*- coding: utf-8
"""
This script opens the ListCI.txt (that list the names of the files in initialConditions file and creates the submit file for condor language)
In particular, each line of the submit file is paralelized in the supercomputer Patung (written in Condor 9.0.11)

It defines where to keep track of the log, output and errors

It also defines which code to read and execute (here the code is

codigoMa.py that instantiates the funcitons in __init__py)
with the initial conditions defined in the submit file

"""

with open('listCI_local.txt', "r") as f:
    condiciones = f.readlines()

subdir = "initialConditions"


head = """
getenv = true
executable = codigoMa.py

output = bitacora/local_$(Process).out
error = bitacora/local_$(Process).err
log = bitacora/local_$(Process).log

request_cpus = 3

"""

with open("submit_file_local.sub", "w") as out:
    out.write(head)
    for line in condiciones:
        out.write("arguments = --subdir %s --code %s\n" % (subdir,line.strip()))
        out.write("queue\n\n")
