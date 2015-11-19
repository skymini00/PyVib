# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 12:45:03 2015

@author: OHNS
"""
import numpy as np
import re   # regular expressions
from DebugLog import DebugLog

class AudioHardware:
    def __init__(self):
        self.micVoltsPerPascal = 1
        self.speakerOutputRng = (10e-3, 2.82)
        self.speakerL_daqChan = "Dev1/ao0"
        self.speakerR_daqChan = "Dev1/ao1"
        self.daqTrigChanOut = "Dev1/port0/line0"
        self.daqTrigChanIn = "/Dev1/PFI1"
        self.mic_daqChan = "Dev1/ai0"
        self.usingAttenuator = True  
        self.attenL_daqChan = "Dev1/line1:3"
        self.attenR_daqChan = "Dev1/line4:6"
        self.DAQOutputRate = 500e3       # output rate of the DAQ in Hz
        self.DAQInputRate = 250e3       # output rate of the DAQ in Hz
        self.maxAtten = 60                # maximum attenuation level possible
        self.speakerCalVolts = 0.1        # output voltage of the speaker calibration
        
         # speaker calibration frequencies, in Hz, 2D array, first index is left speaker, second index is right speaker
        self.speakerCalFreq = np.array([[20, 100e3], [20, 100e3]]) 
        
         # the speaker calibration, 2D array of dB responses, first index is left speaker, second idex is right speaker 
        self.speakerCal = np.array([[60, 60], [60, 60]]) 
        
         # the speaker calibration, 2D array of dB responses, first index is left speaker, second idex is right speaker 
        self.speakerCalPhase = np.array([[0, 0], [0, 0]]) 
        
        self.DAQdevice = 'Dev1'

    # load in the speaker callibration from the processsed data - check OCTProcessing.SpeakerCalData for format
    def loadSpeakerCalFromProcData(self, spCal):
        self.speakerCalVolts = spCal.voltsOut
        self.speakerCalFreq = spCal.freq
        self.speakerCal = spCal.magResp
        try:
            self.speakerCalPhase = spCal.phaseResp
        except:
            self.speakerCalPhase = None
        
    def getCalibratedOutputVoltageAndAttenLevel(self, freq, ampdB, speakerNum):
        freqArray = self.speakerCalFreq[speakerNum]
        calArray = self.speakerCal[speakerNum]
        DebugLog.log("AudioHardware.getCalibratedOutputVoltageAndAttenLevel freq= %f freqArray= %s calArray= %s" % (freq, repr(freqArray), repr(calArray)))
        caldBarr = np.interp([freq], freqArray, calArray)
        caldB = caldBarr[0]
        dBdiff = ampdB - caldB
        DebugLog.log("AudioHardware.getCalibratedOutputVoltageAndAttenLevel ampdB= %f caldB= %f dBdiff= %f" % (ampdB, caldB, dBdiff))
        outV = self.speakerCalVolts*(10 ** (dBdiff/20))
        minV = self.speakerOutputRng[0]
        maxV = self.speakerOutputRng[1]
        attenLevel = 0
        if outV > maxV:
            outV = 0
        elif outV < minV:
            attenLevel = np.ceil(20 * np.log10(minV/outV))
            attenLevel = int(attenLevel)
            outV = minV
        if attenLevel > self.maxAtten:
            outV = 0
   
        return (outV, attenLevel)
       
    def getCalibratedOutputPhase(self, freq, speakerNum):
        if self.speakerCalPhase is None:
            return 0
            
        freqArray = self.speakerCalFreq[speakerNum]
        phaseArray = self.speakerCal[speakerNum]
        phaseArr = np.interp([freq], freqArray, phaseArray)
        return phaseArr[0]
        
    def encodeToString(self, prefix=''):
        s = ""
        s = s + "\n" + prefix + "Speaker Output Range= %f %f" % (self.speakerOutputRng[0], self.speakerOutputRng[1]) 
        s = s + "\n" + prefix + "Speaker Left DAQ Channel= %s" % (self.speakerL_daqChan)
        s = s + "\n" + prefix + "Speaker Right DAQ Channel= %s" % (self.speakerR_daqChan)
        s = s + "\n" + prefix + "DAQ TrigChanOut= %s" % (self.daqTrigChanOut)
        s = s + "\n" + prefix + "DAQ TrigChanIn= %s" % (self.daqTrigChanIn)
        s = s + "\n" + prefix + "Mic volts per Pascal= %f" % (self.micVoltsPerPascal)
        s = s + "\n" + prefix + "Mic DAQ Chan= %s" % (self.mic_daqChan)
        s = s + "\n" + prefix + "Output Rate= %f" % (self.DAQOutputRate)
        s = s + "\n" + prefix + "Input Rate= %f" % (self.DAQInputRate)
        s = s + "\n" + prefix + "Using Attenuator= " + repr(self.usingAttenuator)
        s = s + "\n" + prefix + "Left Atten DAQ Chans= %s" % (self.attenL_daqChan) 
        s = s + "\n" + prefix + "Right Atten DAQ Chans= %s" % (self.attenR_daqChan) 
        s = s + "\n" + prefix + "Max Attenuation= %d" % (self.maxAtten)
        s = s + "\n" + prefix + "DAQ Device= %s" % (self.DAQdevice)
        return s
        
    def decodeFromString(self, audio_str):
        lines = re.split('\n', audio_str)  # break up lines into array
        for s in lines:
            x = re.split('=', s)
            if(len(x) < 2):
                continue
            fld = x[0]
            val = x[1]
            if(fld == 'Speaker Output Range'):
                val2 = re.split(' ', val)
                self.speakerOutputRng = (float(val2[1]), float(val2[2]))
            elif(fld == 'Speaker Left DAQ Channel'):
                self.speakerL_daqChan = val
            elif(fld == 'Speaker Right DAQ Channel'):
                self.speakerR_daqChan = val
            elif(fld == 'DAQ TrigChanOut'):
                self.daqTrigChanOut = val
            elif(fld == 'DAQ TrigChanIn'):
                self.daqTrigChanIn = val
            elif(fld == 'Mic volts per Pascal'):
                self.micVoltsPerPascal = float(val)
            elif(fld == 'Mic DAQ Chan'):
                self.mic_daqChan = val
            elif(fld == 'Using Attenuator'):
                self.usingAttenuator = (val == ' True')
            elif(fld == 'Left Atten DAQ Chans'):
                self.attenL_daqChan = val
            elif(fld == 'Right Atten DAQ Chans'):
                self.attenR_daqChan = val
            elif(fld == 'Output Rate'):
                self.DAQOutputRate = float(val)
            elif(fld == 'Input Rate'):
                self.DAQInputRate = float(val)
            elif(fld == 'Max Attenuation'):
                self.maxAtten = int(val)
            elif(fld == "DAQ Device"):
                self.DAQdevice = val
    
    def saveSpeakerCal(self, fileName):
        numPts = self.speakerCalFreq.shape[1]
        left_freq = self.speakerCalFreq[0, :]
        right_freq = self.speakerCalFreq[1, :]
        left_resp = self.speakerCal[0, :]
        right_resp = self.speakerCal[1, :]
        
        # build the frequency and response strings
        left_freq_str = ''
        right_freq_str = ''
        left_resp_str = ''
        right_resp_str = ''
        for n in range(0, numPts):
            left_freq_str = left_freq_str + " " + repr(left_freq[n])
            right_freq_str = right_freq_str + " " + repr(right_freq[n])
            left_resp_str = left_resp_str + " " + repr(left_resp[n])
            right_resp_str = right_resp_str + " " + repr(right_resp[n])
            
        left_ph_str = ''
        right_ph_str = ''
        if self.speakerCalPhase is not None:
            left_phase = self.speakerCalPhase[0, :]
            right_phase = self.speakerCalPhase[1, :]
            for n in range(0, numPts):
                left_ph_str = left_ph_str + " " + repr(left_phase[n])
                right_ph_str = right_ph_str + " " + repr(right_phase[n])
                
        f = open(filename, 'w')            
        f.write('\nCalibration Volts= %f' % (self.speakerCalVolts))
        f.write('\nNum points= %d' % (numPts))
        f.write("\nLeft freq= %s" % (left_freq_str))
        f.write("\nRight freq= %s" % (right_freq_str))
        f.write("\nLeft resp dB= %s" % (left_freq_str))
        f.write("\nRight resp dB= %s" % (right_freq_str))
        if left_ph_str != '':
            f.write("\nLeft phase (rad)= %s" % (left_ph_str))
            f.write("\nRight phase (rad)= %s" % (right_ph_str))
        f.close()
        

    def loadSpeakerCal(self, fileName):
        f = open(filename, 'r')
        s = f.read()   # read the entire file
        f.close()

        lines = re.split('\n', s)        
        for s in lines:
            x = re.split('=', s)
            fld = x[0]
            val = x[1]
            if(fld == 'Calibration Volts'):
                self.speakerCalVolts = float(val)
            elif(fld == 'Num points'):
                numPts = int(val)
                self.speakerCalFreq = np.zeros(2, numPts)
                self.speakerCaal = np.zeros(2, numPts)
                
        for s in lines:
            x = re.split('=', s)
            fld = x[0]
            val = x[1]
            if(fld == 'Left freq'):
                v = rep.split(" ", val)
                v = v[1:]  # get rid of leading whitespace character
                for n in range(0, numPts):
                    self.speakerCalFreq[0, n] = v[n]
            elif(fld == 'Right freq'):
                v = rep.split(" ", val)
                v = v[1:]  # get rid of leading whitespace character
                for n in range(0, numPts):
                    self.speakerCalFreq[1, n] = v[n]
            elif(fld == 'Left resp dB'):
                v = rep.split(" ", val)
                v = v[1:]  # get rid of leading whitespace character
                for n in range(0, numPts):
                    self.speakerCal[0, n] = v[n]
            elif(fld == 'Right resp dB'):
                v = rep.split(" ", val)
                v = v[1:]  # get rid of leading whitespace character
                for n in range(0, numPts):
                    self.speakerCal[1, n] = v[n]
            elif(fld == 'Left phase'):
                v = rep.split(" ", val)
                v = v[1:]  # get rid of leading whitespace character
                for n in range(0, numPts):
                    self.speakerCalPhase[0, n] = v[n]
            elif(fld == 'Right phase'):
                v = rep.split(" ", val)
                v = v[1:]  # get rid of leading whitespace character
                for n in range(0, numPts):
                    self.speakerCalPhase[1, n] = v[n]                    
                    
                    
def readAudioHWConfig(filepath):    
    hw = AudioHardware()
    
    f = open(filepath, "r")
    txt = f.read()
    f.close()
    hw.decodeFromString(txt)
    
    return hw

# make the digital signal to set the attenuator level for given dB level
def makeLM1972AttenSig(attenLvl):
    numpts = 16*3 + 3
    loadSig = np.zeros(numpts, dtype=np.uint8)
    
    # load signal 
    loadSig[0] = 1 
    loadSig[numpts-1] = 1  
    clkSig = np.zeros(numpts, dtype=np.uint8)
    dataSig =  np.zeros(numpts, dtype=np.uint8)
    
    for n in range(0, 16):
        clkSig[n*3 + 2] = 1
        
    # generate number in binary by repeated divinsion by 2
    for n in range(7, -1, -1):
        idx = n*3 + 2
        r = attenLvl % 2
        dataSig[idx-1:idx+2] = r
        attenLvl = attenLvl // 2
        
    # generate signal output as 32-bit number because thats what NIDAQ uses for writing digital 
    sig = np.zeros(numpts, dtype=np.uint32)
    loadMask = 1 << 2
    dataMask = 1 << 1
    clkMask = 1 << 0
    for n in range(0, numpts):
        sig[n] = (loadMask*loadSig[n]) | (dataMask*dataSig[n]) | (clkMask*clkSig[n])
        
    return sig

    
if __name__ == "__main__":
    from DAQHardware import DAQHardware
    
    sig = makeLM1972AttenSig(30)
    for n in range(0, len(sig)):
        print("%10x" % sig[n])    
        
    daqHW = DAQHardware()
    audioHW = AudioHardware()
    
    # outLines = audioHW.attenL_daqChan
    # lf.attenL_daqChan = "Dev1/line1:3"
    outLines = "PXI1Slot2/port0/line1:3"
    daqHW.sendDigOutCmd(outLines, sig)
    