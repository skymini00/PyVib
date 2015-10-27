# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 12:22:44 2015

@author: OHNS
"""

import OCTCommon
from DebugLog import DebugLog

from OCTProtocolParams import *
import numpy as np
from PyQt4 import QtCore, QtGui, uic
import AudioHardware
import pickle
import os
import sys
import traceback

class SpeakerCalData:
    def __init__(self, audioParams):
        self.voltsOut = 0.1
        self.freq = audioParams.freq
        numFreq = audioParams.getNumFrequencies()
        self.magResp = np.zeros((2, numFreq))
        self.magResp[:, :] = np.NaN
        self.phaseResp = np.zeros((2, numFreq))
        self.phaseResp[:, :] = np.NaN


class MicData:
    def __init__(self):
        self.raw = None   # raw mic response
        self.t = None   # time
        self.fft_mag = None
        self.fft_phase = None
        self.fft_freq = None
        self.stim_freq_mag = None
        self.stim_freq_phase = None
        
def makeSpeakerCalibrationOutput(freq, audioHW, audioParams):
    outV = 100e-3
    outputRate = audioHW.DAQOutputRate
    trialDur = 1e-3*audioParams.getTrialDuration(80)
    trialPts = np.ceil(trialDur * outputRate)
    stimEnv = 1e-3*audioParams.stimEnvelope
    envPts = np.ceil(stimEnv * outputRate)
    t = np.linspace(0, trialDur, trialPts)
    sig = outV*np.sin(2*np.pi*1000*freq*t)
    envFcn = np.ones((trialPts))
    envFcn[0:envPts] = np.linspace(0, 1, envPts)
    envFcn[trialPts-envPts:] = np.linspace(1, 0, envPts)
    sig = sig*envFcn
    
    return sig
   
def processSpkCalData(mic_data, freq, freq_idx, audioParams, inputRate, speakerCalIn, spkNum):
    # print("SpeakerCalProtocol: processData: mic_data=" + repr(mic_data))
    # ensure data is 1D
    if len(mic_data.shape) > 1:
        mic_data = mic_data[:, 0]
        
    numpts = len(mic_data)
    DebugLog.log("SpeakerCalProtocol: processData: numpts= %d" % (numpts))

    t = np.linspace(0, numpts/inputRate, numpts)
    numfftpts = numpts*2
    winfcn = np.hanning(numpts)
    mic_fft = np.fft.fft(winfcn*mic_data, numfftpts)
    endIdx = np.ceil(numfftpts/2)
    mic_fft = mic_fft[0:endIdx]
    mic_fft_mag = 2*np.abs(mic_fft)
    
    # convert to dB, correctting for RMS and FFT length
    fftrms_corr = 2/(numpts*np.sqrt(2))
    mic_fft_mag = 20*np.log10(fftrms_corr*mic_fft_mag/20e-6)   # 20e-6 pa
    
    mic_fft_phase = np.angle(mic_fft)
    mic_freq = np.linspace(0, inputRate/2, endIdx)
    fIdx = int(np.floor(freq*numfftpts/inputRate))
    DebugLog.log("SpeakerCalibration: processData: freq= %f fIdx= %d" % (freq, fIdx))

    stim_freq_mag = np.NAN
    stim_freq_phase = np.NAN

    try:            
        mag_rgn = mic_fft_mag[fIdx-1:fIdx+1]
        phase_rgn = mic_fft_phase[fIdx-1:fIdx+1]
        maxIdx = np.argmax(mag_rgn)
        stim_freq_mag = mag_rgn[maxIdx]
        stim_freq_phase = phase_rgn[maxIdx]
    except Exception as ex:
        DebugLog.log(ex)
    
    DebugLog.log("SpeakerCalibration: processData: stim_freq_mag= %f stim_freq_phase= %f" % (stim_freq_mag, stim_freq_phase))
    micData = MicData()
    micData.raw = mic_data
    micData.t = t
    micData.fft_mag = mic_fft_mag
    micData.fft_phase = mic_fft_phase
    micData.fft_freq = mic_freq
    micData.stim_freq_mag = stim_freq_mag
    micData.stim_freq_phase = stim_freq_phase

    speakerCalIn.magResp[spkNum, freq_idx] = stim_freq_mag
    speakerCalIn.phaseResp[spkNum, freq_idx] = stim_freq_phase
        
    return micData, speakerCalIn
    
    # save the processed data of this protocol
def saveSpeakerCal(spCalData, saveDir):
    filepath = os.path.join(saveDir, 'speaker_cal_last.pickle')
    f = open(filepath, 'wb')
    pickle.dump(spCalData, f)
    f.close()
   
def loadSpeakerCal(filepath):
    f = open(filepath, 'rb')
    spCal = pickle.load(f)
    
    f.close()
    
    return spCal
    
def runSpeakerCal(appObj, testMode=False):
    DebugLog.log("runSpeakerCal")
    appObj.tabWidget.setCurrentIndex(1)
    appObj.doneFlag = False
    appObj.isCollecting = True
    # trigRate = octfpga.GetTriggerRate()
    audioHW = appObj.audioHW
    outputRate = audioHW.DAQOutputRate
    inputRate = audioHW.DAQInputRate
    
    if testMode:
        testDataDir = os.path.join(appObj.basePath, 'exampledata', 'Speaker Calibration')
        filePath = os.path.join(testDataDir, 'AudioParams.pickle')
        f = open(filePath, 'rb')
        audioParams = pickle.load(f)
        f.close()
    else:
        audioParams = appObj.getAudioParams()
    numSpk = audioParams.getNumSpeakers()
    
    if not testMode:
        from DAQHardware import DAQHardware
        daq = DAQHardware()

    chanNamesIn= [ audioHW.mic_daqChan]
    micVoltsPerPascal = audioHW.micVoltsPerPascal

    spCal = SpeakerCalData(audioParams)
    try:
        frameNum = 0
        isSaveDirInit = False
        saveOpts = appObj.getSaveOpts()
        for spkNum in range(0, numSpk):
            chanNameOut = audioHW.speakerL_daqChan 
            attenLines = audioHW.attenL_daqChan
            spkIdx = 0
                
            if (numSpk == 1 and audioParams.speakerSel == Speaker.RIGHT) or spkNum == 2:
                chanNameOut = audioHW.speakerR_daqChan
                attenLines = audioHW.attenR_daqChan
                spkIdx = 1
    
            freq_array = audioParams.freq[spkIdx, :]
            DebugLog.log("freq_array=" + repr(freq_array))
            freq_idx = 0
            for freq in freq_array:
                spkOut = makeSpeakerCalibrationOutput(freq, audioHW, audioParams)    
                npts = len(spkOut)
                t = np.linspace(0, npts/outputRate, npts)
                
                pl = appObj.plot_spkOut
                pl.clear()
                endIdx = int(5e-3 * outputRate)        # only plot first 5 ms
                pl.plot(t[0:endIdx], spkOut[0:endIdx], pen='b')
                        
                attenSig = AudioHardware.makeLM1972AttenSig(0)
                numInputSamples = int(inputRate*len(spkOut)/outputRate) 
                
                if testMode:
                    mic_data = OCTCommon.loadRawData(testDataDir, frameNum, dataType=3)                    
                else:
                    daq.sendDigOutCmd(attenLines, attenSig)
                    # setup the output task
                    daq.setupAnalogOutput([chanNameOut], audioHW.daqTrigChanIn, int(outputRate), spkOut)
                    daq.startAnalogOutput()
                    
                    # setup the input task
                    daq.setupAnalogInput(chanNamesIn, audioHW.daqTrigChanIn, int(inputRate), numInputSamples) 
                    daq.startAnalogInput()
                
                    # trigger the acquiisiton by sending ditital pulse
                    daq.sendDigTrig(audioHW.daqTrigChanOut)
                    
                    mic_data = daq.readAnalogInput()
                    mic_data = mic_data/micVoltsPerPascal

                
                if not testMode:
                    daq.stopAnalogInput()
                    daq.stopAnalogOutput()
                    daq.clearAnalogInput()
                    daq.clearAnalogOutput()
                
                npts = len(mic_data)
                t = np.linspace(0, npts/inputRate, npts)
                pl = appObj.plot_micRaw
                pl.clear()
                pl.plot(t, mic_data, pen='b')
                
                labelStyle = appObj.xLblStyle
                pl.setLabel('bottom', 'Time', 's', **labelStyle)
                labelStyle = appObj.yLblStyle
                pl.setLabel('left', 'Response', 'Pa', **labelStyle)
                
                micData, spCal = processSpkCalData(mic_data, freq*1000, freq_idx, audioParams, inputRate, spCal, spkIdx)
                
                pl = appObj.plot_micFFT
                pl.clear()
                df = micData.fft_freq[1] - micData.fft_freq[0]
                nf = len(micData.fft_freq)
                i1 = int(1000*freq_array[0]*0.9/df)
                i2 = int(1000*freq_array[-1]*1.1/df)
                DebugLog.log("SpeakerCalibration: df= %0.3f i1= %d i2= %d nf= %d" % (df, i1, i2, nf))
                pl.plot(micData.fft_freq[i1:i2], micData.fft_mag[i1:i2], pen='b')
                labelStyle = appObj.xLblStyle
                pl.setLabel('bottom', 'Frequency', 'Hz', **labelStyle)
                labelStyle = appObj.yLblStyle
                pl.setLabel('left', 'Magnitude', 'db SPL', **labelStyle)
                
                pl = appObj.plot_micMagResp
                pl.clear()
                pl.plot(1000*spCal.freq[spkIdx, :], spCal.magResp[spkIdx, :], pen="b", symbol='o')
                labelStyle = appObj.xLblStyle
                pl.setLabel('bottom', 'Frequency', 'Hz', **labelStyle)
                labelStyle = appObj.yLblStyle
                pl.setLabel('left', 'Magnitude', 'db SPL', **labelStyle)
                
                freq_idx += 1
                
                if appObj.getSaveState():
                    if not isSaveDirInit:
                        saveDir = OCTCommon.initSaveDir(saveOpts, 'Speaker Calibration', audioParams=audioParams)
                        isSaveDirInit = True
    
                    if saveOpts.saveRaw:
                        OCTCommon.saveRawData(mic_data, saveDir, frameNum, dataType=3)
                    
                QtGui.QApplication.processEvents() # check for GUI events, such as button presses
                
                # if done flag, break out of loop
                if appObj.doneFlag:
                    break
                
                frameNum += 1

                
            # if done flag, break out of loop
            if appObj.doneFlag:
                break
                
        if not appObj.doneFlag:
            saveDir = appObj.settingsPath
            saveSpeakerCal(spCal, saveDir)
            appObj.audioHW.loadSpeakerCalFromProcData(spCal)
            appObj.spCal = spCal            
            
    except Exception as ex:
        traceback.print_exc(file=sys.stdout)
        QtGui.QMessageBox.critical (appObj, "Error", "Error during scan. Check command line output for details")           
        
    8# update the audio hardware speaker calibration                     
    appObj.isCollecting = False
    QtGui.QApplication.processEvents() # check for GUI events, such as button presses
    appObj.finishCollection()


    