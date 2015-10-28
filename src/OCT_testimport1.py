# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 19:26:37 2015

@author: OHNS
"""

from ctypes import *
import numpy.ctypeslib as np_ctypes
from enum import Enum
import numpy as np
import matplotlib.pyplot as plt


class OCTSetup(Enum):
    OCT_50kHz = 0
    OCT_200kHz = 1
    
class DLLInterface:
    # load in DLL
    # see header file LV_OCT_Raw_Data.h for C function definitions
    pass

def InitFPGA(setup):
    err = 1
    return err
    
def ConfigAcq(numTriggers, samplesPerTrigger):
    err = 2
    return err
    
def GrabData(numTriggers, samplesPerTrigger):    
    err = 3    

    t=np.linspace(0,samplesPerTrigger/50000,samplesPerTrigger)
    ch0_data=np.random.random((numTriggers, samplesPerTrigger))       
    ch1_data=np.random.random((numTriggers, samplesPerTrigger))
    samplesRemaining=0
    timeElapsed=0
    return (err, ch0_data, ch1_data, samplesRemaining, timeElapsed)

def CloseFPGA():
    err = 4
    return err

if __name__ == "__main__":
    plotColors = ['-r', '-b', '-g', '-c', '-m', '-y', '-k']
    numTriggers = 10
    samplesPerTrigger = 2048
    err = InitFPGA(0)
    print("Init FPGA: err = %d" % err)
    
       
    try:
        err = ConfigAcq(numTriggers, samplesPerTrigger)
        print("Config Acq: err = %d" % err)
        (err, ch0_data, ch1_data, samplesRemaining, timeElapsed) = GrabData(numTriggers, samplesPerTrigger)
        print("Grab Data: err = %d" % err)
        print("Grab Data: ch0_data.shape = %s ch1_data.shape= %s" % (repr(ch0_data.shape), repr(ch1_data.shape)))
        print("Grab Data: Samples Remaing= %d TimeElapsed= %d" % (samplesRemaining, timeElapsed))
        plt.figure(1)
        plt.clf()
        for n in range(0, numTriggers):
            clr = plotColors[n % len(plotColors)]
            plt.subplot(2, 1, 1)
            plt.plot(ch0_data[n, :], clr)
            plt.subplot(2, 1, 2)
            plt.plot(ch1_data[n, :], clr)
        plt.subplot(2, 1, 1)
        plt.title("Channel 0")
        plt.subplot(2, 1, 2)
        plt.title("Channel 1")
    except Exception as ex:
        raise ex
    finally:
        err = CloseFPGA()
        print("Close FPGA: err = %d" % err)