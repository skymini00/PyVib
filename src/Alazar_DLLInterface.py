# -*- coding: utf-8 -*-
"""
Created on Thu 02/04/2016

@author: Sangmin Kim
"""

import numpy as np
from ctypes import *
import matplotlib.pyplot as plt

class Alazar_DLLInterface:
    
    def InitInterface(self):
        # load DLLs
        alazar_led_on = CDLL('..\\dll\\AlazarBD_DLL_LED_ON') # for checking AlazarCard
        alazar_led_off = CDLL('..\\dll\\AlazarBD_DLL_LED_OFF') # for checking AlazarCard
        alazar_initconfig = CDLL('..\\dll\\AlazarBD_DLL_InitConfigBoard')
        alazar_acq_raw = CDLL('..\\dll\\AlazarBD_DLL_Acquisition')
        
        # change DLL names
        self.led_on_alazar = alazar_led_on.board_led_on
        self.led_off_alazar = alazar_led_off.board_led_off
        self.initconfig_alazar = alazar_initconfig.InitConfigBoard
        self.acquire_alazar = alazar_acq_raw.Acquisition
        
    def LedOn(self):
        err = self.led_on_alazar()
        if err > 0:
            err = -1
            return err
        
        return err
    
    def LedOff(self):
        err = self.led_off_alazar()
        if err > 0:
            err = -1
            return err
        
        return err
    
    def InitConfigAlazar(self):
        err = self.initconfig_alazar(c_uint(0), c_double(0))
        if err == 0:
            err = -1
            return err
        else:
            err = 0
            return err
    
    def AcquireOCTDataRaw(self, numTriggers, samplesPerTrig):
        # declare variables
        total_length = numTriggers * samplesPerTrig
        raw_data = (c_uint16 * total_length)(0)
        
        err = self.acquire_alazar(c_uint32(samplesPerTrig), c_uint32(numTriggers), c_uint32(1), raw_data, 0) 
        
        ch0_data = np.zeros(total_length, np.float32)
        ch0_data[:] = raw_data[:]
        ch0_data = ch0_data.reshape(numTriggers, samplesPerTrig)
                  
        if err == 0:
            err = -1
            return err, None, None
        else:
            err = 0
            return err, ch0_data, None
                            
        
        
    
         
        

