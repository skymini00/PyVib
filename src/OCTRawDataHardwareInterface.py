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
    lvlib = CDLL('LV_OCT_Raw_Data')
    
    # setup prototypes to functions in library
    _init_fpga = lvlib.InitFPGARaw
    _init_fpga.argtypes = [c_uint16]
    _init_fpga.resttype = c_int32   
    
    _config_acq = lvlib.ConfigureFPGAAcquisitionRaw
    _config_acq.argtypes = [c_uint32, c_uint16]
    _config_acq.resttype = c_int32   

     # AcquireFPGADataRaw(uint32_t NumTriggers, int16_t ch0Data[], 
	#int16_t ch1Data[], int32_t data_len_in, int32_t *data_len_out, 
	#uint32_t *SamplesRemaining, uint32_t *TimeElapsedMs);    

    data_t = np_ctypes.ndpointer(c_int16, flags="C_CONTIGUOUS")
    len_ptr_t = POINTER(c_int32)
    uint32_ptr_t = POINTER(c_uint32)
    _acquire_data = lvlib.AcquireFPGADataRaw
    _acquire_data.argtypes = [c_uint32, data_t, data_t, c_int32, len_ptr_t, uint32_ptr_t, uint32_ptr_t]
    _acquire_data.resttype = c_int32   
    
    _close_fpga = lvlib.CloseFPGARaw
    _close_fpga.argtypes = []
    _close_fpga.resttype = c_int32   

def InitFPGA(setup):
    err = DLLInterface._init_fpga(c_uint16(setup))
    return err
    
def ConfigAcq(numTriggers, samplesPerTrigger):
    err = DLLInterface._config_acq(c_uint32(numTriggers), c_uint16(samplesPerTrigger))
    return err
    
def GrabData(numTriggers, samplesPerTrigger):
    numSamples = c_int32(numTriggers*samplesPerTrigger)
    numSamplesOut = c_int32(0)
    timeElapsed = c_uint32(0)
    samplesRemaining = c_uint32(0)
    ch0_data = np.zeros((numTriggers, samplesPerTrigger), dtype=np.int16)
    ch1_data = np.zeros((numTriggers, samplesPerTrigger), dtype=np.int16)
    err = DLLInterface._acquire_data(c_uint32(numTriggers), ch0_data, ch1_data, numSamples, byref(numSamplesOut), byref(samplesRemaining), byref(timeElapsed))
    print()
    
    return (err, ch0_data, ch1_data, samplesRemaining.value, timeElapsed.value)

def CloseFPGA():
    err = DLLInterface._close_fpga()
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