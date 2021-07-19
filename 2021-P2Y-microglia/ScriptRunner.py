import sys
import pickle as pk
import numpy as np
import math
import analyzeGotran as ao
import subprocess as sb
from subprocess import PIPE
import shlex
import time

def gotranMicroglia(sim_time        = 1000,
                    ATP             = 0, # in uM
                    period          = 3,
                    data_name1      = 'Cai',
                    data_name2      = None,
                    data_name3      = None,
                    data_name4      = None,
                    data_name5      = None,
                    data_name6      = None,
                    data_name7      = None,
                    data_name8      = None,
                    data_name9      = None,
                    data_name10     = None,
                    output_name     = 'p2xp2y',
                    output_switch   = 0,
                    ode_file_name   = 'p2xp2y',
                    removePickle    = 0,
                    timePrint       = 1,
                    **kwargs):
    
    start = time.time()

    ######### Name of files #############################
    nameOfODEfile= ode_file_name+'.ode'
    simTime = str(sim_time) # in milliseconds
    outputName = output_name
    ######################################################

    ######### Do not touch ############
    #mainCommand = 'singularity exec /home/bending456/singularity-img/gtr_env_new.img python2 dcBen.py -odeName'
    mainCommand = 'python2 dcBen.py -odeName'
    basicSetup = '-dt 0.001 -dSr 1 -jit -iters 1 -T'
    extra1 = ['-name']
    #####################################

    ######### Assign new parameters if necessary #########
    variables = ['-var','stim_amplitude', str(ATP), # ATP concentration in [nM]
                 '-var','stim_period',    str(period)] 
    ######################################################
    
    ######### Creating command for variable adjustment #######
    addedArg = []

    for key, value in kwargs.items():
        command = '-var ' + key + ' ' + str(value)
        addedArg = addedArg + shlex.split(command)
    ##########################################################

    ##### Executing calculation ####################################################################################################
    inputArg = shlex.split(mainCommand) + [nameOfODEfile] + shlex.split(basicSetup) + [simTime] + variables\
               + addedArg\
               + extra1 + [outputName]
    out = sb.Popen(inputArg,stdout=PIPE).communicate()[0]
    ################################################################################################################################
    
    ############ Printing Output #########################
    if output_switch == 1:
        print(out)
    ######################################################

    ############ Storing Data ######################
    with open(outputName+"_cat.pickle", 'rb') as f:
        py2data = pk.load(f, encoding='latin1')    
    temp1 = ao.GetData(py2data,data_name1)
    y = np.vstack([temp1.t,temp1.valsIdx])
    if data_name2 != None:
        temp2 = ao.GetData(py2data,data_name2)
        y = np.vstack([y,temp2.valsIdx])
    if data_name3 != None:
        temp3 = ao.GetData(py2data,data_name3)
        y = np.vstack([y,temp3.valsIdx])
    if data_name4 != None:
        temp4 = ao.GetData(py2data,data_name4)
        y = np.vstack([y,temp4.valsIdx])
    if data_name5 != None:
        temp5 = ao.GetData(py2data,data_name5)
        y = np.vstack([y,temp5.valsIdx])
    if data_name6 != None:
        temp6 = ao.GetData(py2data,data_name6)
        y = np.vstack([y,temp6.valsIdx])
    if data_name7 != None:
        temp7 = ao.GetData(py2data,data_name7)
        y = np.vstack([y,temp7.valsIdx])
    if data_name8 != None:
        temp8 = ao.GetData(py2data,data_name8)
        y = np.vstack([y,temp8.valsIdx])
    if data_name9 != None:
        temp9 = ao.GetData(py2data,data_name9)
        y = np.vstack([y,temp9.valsIdx])
    if data_name10 != None:
        temp10 = ao.GetData(py2data,data_name10)
        y = np.vstack([y,temp10.valsIdx])
    ################################################
    
    ### Cleaning pickle files from the parent directory ### 
    if removePickle == 1:
        sb.call(['mv',outputName+'_cat.pickle','./Data'])
        sb.call(['rm','-f',outputName+'_1.pickle'])
    
    #######################################################
    if timePrint == 1:
        print(" -------------- %s seconds --------------" % (time.time() - start))
        print(" ------------ End of Simulation -----------")
    
    return y