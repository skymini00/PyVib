# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 16:27:38 2015

@author: OHNS
"""

import numpy as np

from DebugLog import DebugLog
from ctypes import byref
import time
import PyDAQmx as daqmx

class DAQHardware:
    def __init__(self):
        highSamples = 1000
        numpts = 3 * highSamples
        #dn = 100000
        doData = np.zeros((numpts,), dtype=np.uint32)
        doData[highSamples:2*highSamples] =  2**32 - 1
        self.doTrigData = doData
        
        
    # write single value to maultiple channels 
    def writeValues(self, chanNames, data):
        DebugLog.log("DAQhardware.writeValue(): chanNames= %s val= %s" % (repr(chanNames), repr(data)))
        
        self.analog_output = None   # ensure the output task is closed
        samplesWritten = daqmx.int32()
        analog_output = daqmx.Task()
        data = np.vstack((data, data))
        data = data.transpose()
        data = np.require(data, np.double, ['C', 'W'])
        numSamples = 2
        outputRate = 1000
        
        for chanName in chanNames:
            analog_output.CreateAOVoltageChan(chanName,"",-10.0,10.0, daqmx.DAQmx_Val_Volts, None)
        analog_output.CfgSampClkTiming("",outputRate, daqmx.DAQmx_Val_Rising, daqmx.DAQmx_Val_FiniteSamps, numSamples)
        analog_output.WriteAnalogF64(numSampsPerChan=numSamples, autoStart=True,timeout=1.0, dataLayout=daqmx.DAQmx_Val_GroupByChannel, writeArray=data, reserved=None, sampsPerChanWritten=byref(samplesWritten))
        DebugLog.log("DAQhardware.setupAnalogOutput(): Wrote %d samples" % samplesWritten.value)
        
        # wait until write is completeled
        isDone = False
        isDoneP = daqmx.c_ulong()
        while not isDone:
            err = analog_output.IsTaskDone(byref(isDoneP))
            isDone = isDoneP.value != 0
                
        analog_output = None
        
        
    def setupAnalogOutput(self, chanNames, digTrigChan, outputRate, data):
        numSamples = data.shape[0]

        DebugLog.log("DAQhardware.setupAnalogOutput(): chanNames= %s digTrigChan= %s outputRate= %f numSamples= %f" % (chanNames, digTrigChan, outputRate, numSamples))
        
        self.clearAnalogOutput()  # ensure the output task is closed
        samplesWritten = daqmx.int32()
        analog_output = daqmx.Task()
        data = np.require(data, np.float, ['C', 'W'])
        for chanName in chanNames:
            analog_output.CreateAOVoltageChan(chanName,"",-10.0,10.0, daqmx.DAQmx_Val_Volts, None)
        analog_output.CfgSampClkTiming("",outputRate, daqmx.DAQmx_Val_Rising, daqmx.DAQmx_Val_FiniteSamps, numSamples)
        analog_output.CfgDigEdgeStartTrig(digTrigChan, daqmx.DAQmx_Val_Rising) 
        #analog_output.WriteAnalogF64(numSampsPerChan=numSamples, autoStart=False,timeout=3.0, dataLayout=daqmx.DAQmx_Val_GroupByChannel, writeArray=data, reserved=None, sampsPerChanWritten=byref(samplesWritten))
        analog_output.WriteAnalogF64(numSampsPerChan=numSamples, autoStart=False,timeout=3.0, dataLayout=daqmx.DAQmx_Val_GroupByScanNumber, writeArray=data, reserved=None, sampsPerChanWritten=byref(samplesWritten))
        DebugLog.log("DAQhardware.setupAnalogOutput(): Wrote %d samples" % samplesWritten.value)
        
        self.analog_output = analog_output

    def startAnalogOutput(self):
        self.analog_output.StartTask()
        
    def stopAnalogOutput(self):
        self.analog_output.StopTask()

    def setupAnalogInput(self, chanNames, digTrigChan, inputRate, numInputSamples): 
        # ensure old task has been cosed
        self.clearAnalogInput()
        
        DebugLog.log("DAQhardware.setupAnalogOutput(): chanNames= %s digTrigChan= %s inputRate= %f numInputSamples= %f" % (chanNames, digTrigChan, inputRate, numInputSamples))
        ## DAQmx Configure Code
        # analog_input.CreateAIVoltageChan("Dev1/ai0","",DAQmx_Val_Cfg_Default,-10.0,10.0,DAQmx_Val_Volts,None)
        analog_input = daqmx.Task()
        for chanName in chanNames:
            analog_input.CreateAIVoltageChan(chanName,"", daqmx.DAQmx_Val_Cfg_Default,-10.0,10.0, daqmx.DAQmx_Val_Volts,None)
            
        analog_input.CfgSampClkTiming("", inputRate, daqmx.DAQmx_Val_Rising, daqmx.DAQmx_Val_FiniteSamps, numInputSamples)
        analog_input.CfgDigEdgeStartTrig(digTrigChan, daqmx.DAQmx_Val_Rising) 
        
        numCh = len(chanNames)
        self.dataIn = np.zeros((numInputSamples*numCh,))
        #
        ## DAQmx Start Code
        self.analog_input = analog_input

    def startAnalogInput(self):
        self.analog_input.StartTask()
        
    def stopAnalogInput(self):
        self.analog_input.StopTask()
        
    def sendDigTrig(self, trigOutLine):
        # setup the digital trigger
        DebugLog.log("DAQhardware.sendDigTrig(): trigOutLine= %s " % trigOutLine)
        dig_out = daqmx.Task()
        dig_out.CreateDOChan(trigOutLine, "", daqmx.DAQmx_Val_ChanForAllLines)

        doSamplesWritten = daqmx.int32()
        doData = self.doTrigData
        numpts = len(doData)
        dig_out.WriteDigitalU32(numSampsPerChan=numpts, autoStart=True, timeout=1.0, dataLayout=daqmx.DAQmx_Val_GroupByChannel, writeArray=doData, reserved=None, sampsPerChanWritten=byref(doSamplesWritten))
        DebugLog.log("DAQhardware.sendDigTrig():  Wrote %d samples" % doSamplesWritten.value)
        dig_out.ClearTask()
        
    def readAnalogInput(self, timeout=3.0): 
        ## DAQmx Read Code
        read = daqmx.int32()
        numSamplesIn = len(self.dataIn)
        self.analog_input.ReadAnalogF64(numSamplesIn, timeout, daqmx.DAQmx_Val_GroupByChannel, self.dataIn, numSamplesIn, byref(read), None)
        DebugLog.log("DAQhardware.sendDigTrig(): Read %s samples" % repr(read))
        data = self.dataIn
        
        return data
        
    def sendDigOutCmd(self, outLines, outCmd, timeBetweenPts=0.0001):
        dig_out = daqmx.Task()
        dig_out.CreateDOChan(outLines, "", daqmx.DAQmx_Val_ChanForAllLines)
        # print(dig_out)
        numpts = len(outCmd)
        doSamplesWritten = daqmx.c_int32()
        # dig_out.WriteDigitalU32(numSampsPerChan=numpts, autoStart=True, timeout=1.0, dataLayout=daqmx.DAQmx_Val_GroupByChannel, writeArray=outCmd, reserved=None, sampsPerChanWritten=byref(doSamplesWritten))
        # DebugLog.log("Digital output: Wrote %d samples" % doSamplesWritten.value)
        for n in range(0, numpts):
            outCmd_n = np.array(outCmd[n])
            dig_out.WriteDigitalU32(numSampsPerChan=1, autoStart=True, timeout=0.01, dataLayout=daqmx.DAQmx_Val_GroupByChannel, writeArray=outCmd_n, reserved=None, sampsPerChanWritten=byref(doSamplesWritten))            
            time.sleep(timeBetweenPts)

    def waitDoneTask(task, timeout):
        err = 0
        isDone = False
        isDoneP = daqmx.c_ulong()
        tElapsed = 0
        lastTime = time.time()
        while not isDone:
            task.IsTaskDone(byref(isDoneP))
            DebugLog.log("waitDoneTask(): isDoenP= %s" % repr(isDoneP))
            isDone = isDoneP.value != 0
            timeNow = time.time()
            tElapsed += timeNow - lastTime
            lastTime = timeNow
            if timeout >= 0 and tElapsed > timeout:
                if not isDone:
                    err = -1
                isDone = True
            
        return err
        
    def waitDoneOutput(self, timeout=-1):
        err = DAQHardware.waitDoneTask(self.analog_output, timeout)        
        return err
        
    def waitDoneInput(self, timeout=-1):
        err = DAQHardware.waitDoneTask(self.analog_input, timeout)        
        return err
        
    def clearAnalogOutput(self):
        self.analog_output = None
        
    def clearAnalogInput(self):
        self.analog_input = None
        

        
        
import matplotlib.pyplot as plt
import copy

def runAOTest(daqHW,):
    chanNames = ['Dev1/ao0', 'Dev1/ao3']    
    chanNames = ['PXI1Slot2/ao0', 'PXI1Slot2/ao1']    

    data = np.zeros(2)
    daqHW.writeValues(chanNames, data)
    
    outputRate = 200000
    amp=1
    amp2 = amp/2
    outT = 100e-3
    npts = int(outputRate*outT)
    t = np.linspace(0, outT, npts)
    freq = 2000
    freq2 = freq*2
    cmd_x = amp*np.sin(2*np.pi*freq*t)
    cmd_y = amp2*np.sin(2*np.pi*freq2*t)
    
#    plt.figure(1);
#    plt.clf()
#    plt.plot(cmd_x, '-r')
#    plt.plot(cmd_y, '-b')
#    plt.show()
#    plt.draw()
    
    outData = np.vstack((cmd_x, cmd_y))
    #outData = cmd_x
    trigChan = '/Dev1/PFI0'
    trigOutLine = 'Dev1/port0/line0'
    chanName = chanNames[1]
    
#    numSamples = npts
#    samplesWritten = daqmx.int32()
#    analog_output = daqmx.Task()
#    outData = np.require(outData, np.float, ['C', 'W'])

#    analog_output.CreateAOVoltageChan(chanName,"",-10.0,10.0, daqmx.DAQmx_Val_Volts, None)
    #for chanName in chanNames:
    #    analog_output.CreateAOVoltageChan(chanName,"",-10.0,10.0, daqmx.DAQmx_Val_Volts, None)
#    analog_output.CfgSampClkTiming("",outputRate, daqmx.DAQmx_Val_Rising, daqmx.DAQmx_Val_FiniteSamps, numSamples)
#    analog_output.CfgDigEdgeStartTrig(trigChan, daqmx.DAQmx_Val_Rising) 
#    analog_output.WriteAnalogF64(numSampsPerChan=numSamples, autoStart=False,timeout=3.0, dataLayout=daqmx.DAQmx_Val_GroupByChannel, writeArray=outData, reserved=None, sampsPerChanWritten=byref(samplesWritten))

    dataOut = copy.copy(outData)
    daqHW.setupAnalogOutput(chanNames, trigChan, outputRate, dataOut.transpose())
    for n in range(0, 20):
        daqHW.startAnalogOutput()
        print("n= ", n)
        #nalog_output.StartTask()
        daqHW.sendDigTrig(trigOutLine)
        
        # wait unti output is finsihed and clean up tasks
        err = daqHW.waitDoneOutput(3)
        if err < 0:
            DebugLog.log("waitDoneOutput() err = %s" % repr(err))
        
        daqHW.stopAnalogOutput()
#        isDone = False
#        isDoneP = daqmx.c_ulong()
#        while not isDone:
#            err = analog_output.IsTaskDone(byref(isDoneP))
#            isDone = isDoneP.value != 0
#            
#        analog_output.StopTask()
        time.sleep(0.1)
        
#    analog_output.ClearTask()
    daqHW.clearAnalogOutput()
    
if __name__ == "__main__":

    daqHW = DAQHardware()
    digOutLines = 'PXI1Slot2/port0/line1:2'    
    numpts = 100
    digSig = np.zeros(numpts, dtype=np.uint32)
    digSig[0::2] = 1
    digSig[1::2] = 2
    
    daqHW.sendDigOutCmd(digOutLines, digSig)

    # runAOTest(daqHW)
