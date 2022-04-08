"""
Here we have all the necessary functions without instantiation. This script should be saved in .nwg for the condor 9.0.11 structure.





0. Import the used libraries
"""

from operator import add #esto para crear nuevos directorios para guardar tu trabajo
import numpy as np  ##numpy paquete que acelera todo el manejo de matrices y demás
import random 
import seaborn as sns; sns.set() #parauete para graficar y hacer los mapeos y las rejillas
import pandas as pd 
from matplotlib import cm  #esto para definir paletas de colores manualmente
viridis = cm.get_cmap('viridis', 12) #esto pa usar la paleta de colores viridis




"""
PART A
These functions are called during the dynamics of the simulation

"""


"""
I. GENERAL DYNAMICS dinaMig

This function runs the general dynamics. It uses other functions enlisted below
 
It receives: 
    
a) necessary conditions for ODE to run (dic_ODE)
b) Necessary conditions to define the space parameters (dic_Lattice)
c) Necessary values for migration/diffusion (dic_Migratio)
d) Necessary values for an external perurbation (dic_Perturbation) NON USED


It exits: 
Generates a 4-dimensional matrix. One of the entries is a vector of 1*3 entries: (steps)*(xcoor)*(ycoor)* (s, i, x) values

It calls:
marcadoresPosibles() to define the initial planting arrangement
genera_pob_ini() to generate the initial lattice 
migracionRand() o migracion() according to the call to define a migration matrix
euler_red() to call the ODE in each of the cell of the lattice
 
"""



def dinaMig(dic_ODE, dic_Lattice, dic_Migration, dic_Perturbation):

    marcadores = marcadoresPosibles(dic_Lattice)  
    pasos= dic_ODE["steps"]    
    latticeX= dic_Lattice["lattice_size"][0]
    latticeY= dic_Lattice["lattice_size"][1]
    
    
    #we define the general matrix
    matrizGeneral = np.zeros((pasos, latticeX,latticeY,3))
    #no we have a dynamics matrix that will generate each step of the simulation
    matrizDin = genera_pob_ini(latticeX,latticeY, dic_ODE["condIni"], dic_Lattice["Ires"],marcadores)    
    matrizGeneral[0,:,:,:]= matrizDin[:,:,:] 
    
    tasa_mig= dic_Migration["tasa_Mig"]  
    actuMigra= dic_Migration["actualizacion"]
    
    """
    General dynamic

    Each step takes as an entry an initial lattice (the state matrix) with one combination of initial conditions.
    It then creates a migration matrix where each entry of the matrix represents 
    the migration term to the corresponding entry of the state matrix. 
    Both matrices are integrated with the Euler method. 
    This creates a new state matrix and the cycle is repeated for 30 000 steps (representing 300 days).
    """    
    
    
    for n in range(pasos-1):  
        
        #we generate the migration matriz, using the last line of the dynamic matrix
        
        if n in np.arange(0, pasos-1, actuMigra):
            Pmigra = migracion(matrizDin, tasa_mig)  
           
        else: 
            Pmigra = np.zeros((len(matrizDin), len(matrizDin)))
            
    
        ## Now for each coordinate of the lattice, we define the borders as 0
        ## We run the ode integration for the cells with trees (defined by the marcadores) #we overwrite the matrizDin
        ## we run the ode integration external, for the cell without trees
        
        for i in range(len(matrizDin)):
            for j in range(len(matrizDin[i])):  
                
                
                if i== 0 or i==len(matrizDin)-1 or j==0 or j==len(matrizDin[i])-1: #a las fronteras no les llega nada
                    matrizDin[i,j,:]=0
                
                elif [i,j] in marcadores:
                    #print([i,j])
                    X_eu= euler_red(matrizDin[i,j,:], Pmigra[i,j], dic_ODE["para"], dic_ODE["sizeStep"], dic_ODE["mode"])
                    matrizDin[i,j,:] = X_eu
                else:
                    X_eu= euler_red_ext(matrizDin[i,j,:], Pmigra[i,j], dic_ODE["para"], dic_ODE["sizeStep"], dic_ODE["mode"])
                    matrizDin[i,j,:] = X_eu
                    
                    

        matrizGeneral[n+1,:,:,:]= matrizDin[:,:,:]       #we save the new dynamic matrix in the general matrix
    
    return ([matrizGeneral, marcadores])



##cree esta nueva funcion en 22 maro 2021

"""
II. MarcadoresPosibles

This function defines the planting arrangement
 
It receives: 
a) Necessary conditions to define the space parameters (dic_Lattice)

It exits: 
A list with N entries, where each entry is a [x,y]

It does not call any other function
"""


def marcadoresPosibles(dic_Lattice):
    indPosibles = []
    
    x_arb= dic_Lattice["lattice_size"][0]
    y_arb = dic_Lattice["lattice_size"][1]
    den = dic_Lattice["density"]  # we deinfe density 0.5 for 50 trees
    patron= dic_Lattice["patron"]
    
    numOfTrees = ((x_arb-2)*(y_arb-2))*den
    if patron == "rows":
        for j in np.arange(1, y_arb-1,1):
            for k in np.arange(1, x_arb-1,2):
                indPosibles.append([k,j])
        
        
    elif patron == "spaced":
        for j in np.arange(1, y_arb-1,1):
            ini = j%2
            for k in np.arange(ini+1, x_arb-1,2):
               indPosibles.append([k,j])
        
                
    elif patron == "random":#suponemos que el numero de plantas
        indPosibles = [[random.randint(1, x_arb-2), random.randint(1, y_arb-2)] for i in range(int(numOfTrees))] 
        
        indIni = [[5,6]]
        
        ##for the random pattern we hve the verify that the cell [5,6] contains a tree (the initially infected tree)
     
        if indIni[0] not in indPosibles:
            print("reemplazo de centro")
            eli = random.choice(indPosibles)
            indPosibles.remove(eli)
            indPosibles.append(indIni[0])
        
        
        
    elif patron == "aggregated":  #this code is larger, because it is defined algoritmically and not geometrically
        
        indNucleo = [[int((x_arb/2)-1),int((y_arb/2)-1)], [int((x_arb/2)-1),int((y_arb/2))], [int((x_arb/2)),int((y_arb/2)-1)], [int((x_arb/2)),int((y_arb/2))]] #esto es para crear un nucleo de donde va a crecer de manera uniforme
        
        n = len(indNucleo) #esto es el nucleo
        
        indicesVecinos = [[0,-1],[0,1],[-1,0],[1,0]]
        
        while n< numOfTrees:             
            
            for ind in indNucleo: 
                for i in indicesVecinos: 
                    propuestaInd = [ind[0]+i[0], ind[1]+i[1]]
                                   
                    if n< numOfTrees:                   
                        if propuestaInd not in indNucleo:
                            indNucleo.append(propuestaInd)
                            n = n+1
            
                    else: 
                        break
                                          
                else: 
                    pass                         
        indPosibles= indNucleo      
        
    else: 
        print("non exisisting pattern")  #when we have a tipo calling the function
        
        
    return(indPosibles)
    
    
    
"""
III. genera_pob_ini

This function generates the initial lattice and initial conditions    

It receives:
a) lattice characteristics (x, y) size
b) ODE intial 
c) li_Ires (initial porportion of infected trees)
d) marcadores (vector of trees indices, i.e. the planting arrangement)
    
It exits:
The initial 3dim- matrix [50 x coordinate, 50y coordinates, (s,i,x) initial values]    

Does not call any other function
"""


#def genera_pob_ini(x_arb, y_arb, s_ini, i_ini, x_ini, aleatoria, dic_Lattice):
def genera_pob_ini(x_arb, y_arb, condIni, li_Ires, marcadores):


    matrizParIni = np.zeros((x_arb, y_arb,3))
    matrizParIni[:,:,0]= condIni[0] #s_ini #this way I call all the cuadros
    matrizParIni[:,:,1]= condIni[1] #in:ii
    matrizParIni[:,:,2]= condIni[2] #x ini"
        
  
##we set the frontiers to 0
    matrizParIni[0,:,:]= 0
    matrizParIni[x_arb-1,:,:]= 0
    matrizParIni[:,0,:]= 0
    matrizParIni[:,y_arb-1,:]= 0
    
    indPosibles = marcadores
    
    # we set S leaves =1 (trees presence)
    for i in indPosibles:
        matrizParIni[i[0],i[1],0] = 1
    
    indIni = [[5,6]]
  
    
    for i in indIni:
        matrizParIni[i[0],i[1],1] = li_Ires[0]  #aqui todas inician con lo que le mande a la planta
    
    return (matrizParIni)

"""
IV. Migration/Diffusion matrix

It receives:
a) The dynamical matrix (last step of the previous integration)
b) the diffusion rate
    

It creates
a) A migration matrix that will integrate with the next ode integration

It does not all any other function
"""
def migracion(matriz_Din, tasa_mig):
    m= tasa_mig
    Pmig= np.zeros((len(matriz_Din), len(matriz_Din)))
    

    #the amount of matter that each cell has after the exchange with its neighbours depends on whether is is localted on the borders or in the center
    
    for i in range(len(matriz_Din)):
            for j in range(len(matriz_Din[i])):                
                
                if i== 0 or i==len(matriz_Din)-1 or j==0 or j==len(matriz_Din[i])-1: #a las fronteras no les llega nada
                    pass
                elif (i==1 or i==len(matriz_Din)-2) and (j==1 or j==len(matriz_Din[i])-2):
                    Pmig[i,j]= m*(matriz_Din[i-1,j,2]+matriz_Din[i+1,j,2]+matriz_Din[i,j-1,2]+matriz_Din[i,j+1,2]- 2*matriz_Din[i,j,2])
                
                elif i==1 or i==len(matriz_Din)-2 or j==1 or j==len(matriz_Din[i])-2:
                     Pmig[i,j]= m*(matriz_Din[i-1,j,2]+matriz_Din[i+1,j,2]+matriz_Din[i,j-1,2]+matriz_Din[i,j+1,2]- 3*matriz_Din[i,j,2])
                
                else:
                    Pmig[i,j]= m*(matriz_Din[i-1,j,2]+matriz_Din[i+1,j,2]+matriz_Din[i,j-1,2]+matriz_Din[i,j+1,2]- 4*matriz_Din[i,j,2])
   
    return (Pmig)

      

"""
V. euler_red, euler_red ext 

These functions create the euler integration step for cell with trees and without trees respectively (using the differential  ( Z_t+h = f(s,i,p)*h + Z_t))

They receive:
    a) (s,i,x) of the n-1 step for one specific coordinate
    b) the value of the migration matrix for one specific coordinate
    c) the parameters for the ode
    d) the size step
    e) the mode of the growth function (mono or log)
    
They exit: a list of (s,i, x) for the n step

They call:
    a) the SIX_eu_mig_Ext() and SIX_eu_mig()) to calculate the differential
"""


def euler_red(z_pre, migra_, para_, size_step, mode_): #z es para referir el vector sip que le metes
    h= size_step
    suma= [i * h for i in SIX_eu_mig(z_pre, para_, migra_, mode_)]#aqui hago el calcylo
    #print("int", suma)
    znu= list(map(add,z_pre, suma))     # esto es la z más la suma de la diferencial pero por si quiero ver varios pasos en algun momento
    return(znu) #renglon


def euler_red_ext(z_pre, migra_, para_, size_step, mode_): #z es para referir el vector sip que le metes
    h= size_step
    suma= [i * h for i in SIX_eu_mig_ext(z_pre, para_, migra_, mode_)]#aqui hago el calcylo
    #print("ext", suma)
    znu= list(map(add,z_pre, suma))     # esto es la z más la suma de la diferencial pero por si quiero ver varios pasos en algun momento
    return(znu) #renglon
 

"""
VI. six_eu_mig and six_eu_mig_ext

These functions create the differential

They receive:
    a) (s,i,x) of the n-1 step for one specific coordinate
    b) the value of the migration matrix for one specific coordinate
    c) the parameters for the ode
    d) the mode of the growth function (mono or log)
    
They exit: a list of (s,i, x) differential 

They do not call any other function
"""


def SIX_eu_mig(z_pre, para_, migra_, mode_):
    para= para_
    migra= migra_
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
        dsdt = r*(1-(s))-b1*s*p-b2*s*i   #aqui está    
    else: 
        print ("ERROR DE MODO")
    
    didt = b1*s*p+b2*s*i-g*i
    dpdt= a*i-mu*p+ migra
    dzdt = [dsdt, didt, dpdt]    
    return dzdt



def SIX_eu_mig_ext(z_pre, para_, migra_, mode_):
   #parametros del modelo
    para= para_
    migra= migra_
    z= z_pre
    
    mu= para[5]
    
    p= z[2]
    
    #ecuaciones
    if str(mode_)== "log":
        dsdt = 0
    
    elif str(mode_) == "mono":
        dsdt = 0  #aqui está    
    else: 
        print ("ERROR DE MODO")
    
    didt = 0
    dpdt= -mu*p+ migra #this is the main difference between the integration inside a cell with tree and a cell withut tree
    dzdt = [dsdt, didt, dpdt]    
    return dzdt


"""
PART B
These functions are called after the simulations (they use the general matrix generated during the simulations)
"""

"""
I. CalculadorIndice

This function calculates the H value for the used simulations. (in order to calculate the D value in another script)

It receives
a) The general indices of the trees
b) the dic_lattice 
c) the simulation number
    
It exits:
    The H value for each plant
    

"""

def calculadorIndice(marcadores, dic_Lattice, sim_):
    maxIndice = dic_Lattice["lattice_size"][0] #suponiendo una matriz cuadrada
    minIndice = 0
    
    
    df_Indices = pd.DataFrame({"coordenada":0, "1/D_c":0, "patron":dic_Lattice["patron"], "tamaño":dic_Lattice["lattice_size"], "sim":sim_})  

    
    for coor in marcadores:
        
        D_l = 0 
        D_r = 0
        D_d = 0
        D_u = 0
        n = coor[0]
        m= coor[1]
        
        l= 1
        r=1
        u= 1
        d= 1
        
        while n-l>minIndice:
            #print("A", n-l)
            if [n-l,m] in marcadores:
                D_l = 1/l
                break
            else:
                l= l+1
                
        while n+r<maxIndice:
            #print("B", n+r)
            if [n+r,m] in marcadores:
                D_r = 1/r
                break
            else:
                r= r+1
                
        while m-d>minIndice:
            #print("C", m-d)
            if [n,m-d] in marcadores:
                D_d = 1/d
                break
            else:
                d= d+1
                
        while m+u<maxIndice:
            #print("D", m+u)
            if [n,m+u] in marcadores:
                D_u = 1/u
                break
            else:
                u= u+1
           
        InvD = (D_l + D_d + D_u + D_r)/4    
        df_Indices.loc[len(df_Indices.index)] = [coor, InvD,dic_Lattice["patron"], dic_Lattice["lattice_size"], sim_]
    
    df_Indices =  df_Indices.drop([0,1])
    
    return(df_Indices)




"""
II. creador DF_maPromedio 

This function create a general dataframe (dataFrame_promedios) for R using the general values of the simulation and the average values 

They receive:
    a) The general matrix (the 4-dimensional matrix)
    b) the general values of the simlations (dic_DE, dic_Lattice, di_maigration, sim, Dc, marcadores)
    c) the list of trees (marcadores)
    
It exit: 
   a) A general dataframe with the averaged values and the simulation specifications

It call the promediador_espacio to generate the average matrix
"""


def creadorDF_MaPromedio(matrizGen, dic_ODE, dic_Lattice, dic_Migration, sim, marcadores, Dc_):
    
    print("creadorDFPromedio")
    
    
    Dc = Dc_
    Ires = dic_Lattice["Ires"]
    Iini= Ires[0]
    #Ipro = Ires[1]
    Ipro = 0.01    
    patron =dic_Lattice["patron"]
    parameters = dic_ODE["para"]
    mig = dic_Migration["tasa_Mig"]
    actmigra = dic_Migration["actualizacion"]

    df_matrizPromedio = pd.DataFrame(promediador_espacio(matrizGen, marcadores), columns=["S_prom", "I_prom", "X_prom", "S_sd", "I_sd", "X_sd"])  
    
    
    
    df_matrizPromedio["g"] = parameters[3]  
    df_matrizPromedio["m"] = mig 
    df_matrizPromedio["Iini"] = Iini
    df_matrizPromedio["IresPro"] = Ipro
    df_matrizPromedio["actMigra"] =actmigra
    df_matrizPromedio["sim"]=sim
    df_matrizPromedio["patron"] = patron
    df_matrizPromedio["Dc"] = Dc
    return(df_matrizPromedio)




"""
III. promediador_espacio

This function averages the general values for s,i, x for each step, for all the trees 

it receives:
    a) The general matrix (the 4-dimensional matrix)
    b) the list of trees (marcadores)
    
It exits: 
   a) An average matrix of n steps and 6 values for each step (an average and stand dev for s, i and x)

It does not call any other function
"""


def promediador_espacio(matrizGen, marcadores):
    indPosibles= marcadores


    matriz1 = matrizGen
    #creador NAs
    for i in range(len(matriz1[0])):  #dimensiones x
        for j in range(len(matriz1[0][0])): #dimensiones y
            if [i, j] not in indPosibles:
                matriz1[:,i,j,:] =np.nan

    N= len(matriz1)
    
    MatPromedio= np.zeros((N,6))
    
    
    for n in range(N):
    
        MatPromedio[n,0] = np.nanmean(matriz1[n,:,:,0])
        MatPromedio[n,1] = np.nanmean(matriz1[n,:,:,1])
        MatPromedio[n,2] = np.nanmean(matriz1[n,:,:,2])
        MatPromedio[n,3] = np.nanstd(matriz1[n,:,:,0])
        MatPromedio[n,4] = np.nanstd(matriz1[n,:,:,1])
        MatPromedio[n,5] = np.nanstd(matriz1[n,:,:,2])

    return(MatPromedio)


"""
IV. df_calculadorPicosPlant

This function calculates the maximum value of infection for each plant and the time to get to that value and creates a general data
frame with these values


It receives:
    a) The general matrix (the 4-dimensional matrix)
    b) the list of trees (marcadores)
    c) the General conditions of the simulation (all the dictionaries)
    d) the D value
    
It exit: 
   a) a dataframe with the maxmum infection and time to max infection for each plant

It does not call any other function
"""





def df_calculadorPicosPlanta(matrizGen, dic_ODE, dic_Lattice, dic_Migration, sim, marcadores, Dc_):
    print("creadorPicos")
    
    dataFrame_planta = pd.DataFrame(columns=["Sim", "Plant", "Imax", "TiempoImax", "Iini", "Ipro", "Patron", "Plot", "tasaMig", "g"])

    Ires = dic_Lattice["Ires"]
    Iini= Ires[0]
    Ipro = 1
   # Ipro = Ires[1]
    patron =dic_Lattice["patron"]
    parameters = dic_ODE["para"]
    g= parameters[3]
    mig = dic_Migration["tasa_Mig"]
    Dc = Dc_
    
    
    indPosibles= marcadores

    for a in indPosibles:

        x=a[0]
        y=a[1]
       # print(x,y)
        
        Imax =  np.max(matrizGen[:,x,y,1], axis= None)
        #print(Imax)
        pasosImax= np.where(matrizGen[:,x,y,1]==Imax)[0][0]
        stepSize= 0.01
        tiempoImax = stepSize*pasosImax
        dataFrame_planta = dataFrame_planta.append({"Sim":sim, "Dc": Dc, "Plant":[x,y], "Imax": Imax ,"TiempoImax":tiempoImax, "Ini":Iini , "Ipro":Ipro, "Patron":patron, "tasaMig":mig,"g":g},ignore_index = True)

    return(dataFrame_planta)




