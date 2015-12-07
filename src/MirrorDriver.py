# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 12:44:20 2015

@author: OHNS
"""

from enum import Enum
import re   # regular expressions

class MirrorType(Enum):
    DISSECTING_MICROSCOPE = 0
    MEMS_MICROSCOPE = 1
    MEMS_ENDOSCOPE = 2

class MirrorDriver:
    def __init__(self):
        self.mirrorType = MirrorType.DISSECTING_MICROSCOPE  
        self.voltsPerMillimeter = 2.55      # how many volts to apply to to displace the beam  by one milliemeter
        self.voltRange = (-10, 10)          # range 
        self.settleTime = 1e-3              # time if takes for mirror to settle to small command changes
        self.flybackTime = 10e-3            # time it takes for mirro to flyback across its length
        self.reverseTime = 2e-3             # time for mirror to reverse position
        self.X_daqChan = "Dev1/ao2"
        self.Y_daqChan = "Dev1/ao3"
        self.trig_daqChan = "/Dev1/PXI_Trig2"
        self.DAQoutputRate = 500e3
        self.DAQdevice = 'Dev1'
        self.fastScanMaxFreq = 80e3
        self.newField = ''
        self.skew = 1
        self.phaseAdjust = 0
        self.angularScanFreq = 0
        self.volScanFreq = 0
        self.LPFcutoff = 10
        self.voltsPerMillimeterResonant = 1
        
    # return the  move the mirro to a single (x,y) point, where x and y in millimeters
    def makeMirrorCommand(self, x, y):
        x_cmd = self.voltsPerMillimeter*x
        y_cmd = self.voltsPerMillimeter*y
        return (x_cmd, y_cmd)
        
    def encodeToString(self):
        s = ""
        s = s + "\nType= %s" % repr(self.mirrorType)
        s = s + "\nVolts per millimeter= %f" % (self.voltsPerMillimeter)
        s = s + "\nVolt Range= %f %f" % (self.voltRange[0], self.voltRange[1]) 
        s = s + "\nSettle Time= %f" % (self.settleTime)
        s = s + "\nFlyback Time= %f" % (self.flybackTime)
        s = s + "\nReverse Time= %f" % (self.reverseTime)
        s = s + "\nX DAQ Chan= %s" % (self.X_daqChan)
        s = s + "\nY DAQ Chan= %s" % (self.Y_daqChan)
        s = s + "\nOCT Trigger DAQ Chan= %s"% (self.trig_daqChan)
        s = s + "\nOutput rate= %f" % (self.DAQoutputRate)
        s = s + "\nDAQ Device= %s" % self.DAQdevice
        s = s + "\nFast Scan Max Freqe= %f" % self.fastScanMaxFreq
        
        return s
        
    def decodeFromString(self, mirror_str):
        lines = re.split('\n', mirror_str)    # break up lines into array
        for s in lines:
            x = re.split('=', s)
            if(len(x) < 2):
                continue
            fld = x[0].rstrip()
            val = x[1]
            if(fld == "Volts per millimeter"):
                self.voltsPerMillimeter = float(val)
            elif(fld == "Volt Range"):
                val2 = re.split(' ', val)
                self.voltRange = (float(val2[1]), float(val2[2]))
                print('voltrange', self.voltRange)
            elif(fld == "Settle Time"):
                self.settleTime = float(val)
            elif(fld == "Flyback Time"): 
                self.flybackTime = float(val)
            elif(fld == "Reverse Time"):
                self.reverseTime = float(val)
            elif(fld == "X DAQ Chan"): 
                self.X_daqChan = val
            elif(fld == "Y DAQ Chan"):
                self.Y_daqChan = val
            elif(fld == "OCT Trigger DAQ Chan"):
                self.trig_daqChan = val
            elif(fld == "Output rate"):
                self.DAQoutputRate = float(val)
            elif(fld == "DAQ Device"):
                self.DAQdevice = val
            elif(fld == "Fast Scan Max Freq"):
                self.fastScanMaxFreq = float(val)
            elif(fld == "Type"):
                self.mirrorType = MirrorType(int(val))
            elif(fld == "skew"):
                self.skew = float(val)
            elif(fld == "phaseAdjust"):
                self.phaseAdjust = float(val)
            elif(fld == "angularScanFreq"):
                self.angularScanFreq = float(val)
            elif(fld == "volScanFreq"):
                self.volScanFreq = float(val)
            elif(fld == "LPF cutoff"):
                self.LPFcutoff = float(val)
            elif(fld == "voltsPerMillimeterResonant"):
                self.voltsPerMillimeterResonant = float(val)
    def __repr__(self):        
        return self.encodeToString()
        
def readMirrorDriverConfig(filepath):    
    g = MirrorDriver()
    
    f = open(filepath, "r")
    txt = f.read()
    f.close()
    
    g.decodeFromString(txt) 
    return g