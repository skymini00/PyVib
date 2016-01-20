# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 12:22:44 2015

@author: OHNS
"""

import OCTCommon
from DebugLog import DebugLog
import VolumeScan
import AudioHardware
import BScan
import JSOraw

# import OCTFPGAProcessingInterface as octfpga

from PyQt4 import QtCore, QtGui, uic
from ROIImageGraphicsView import *
from OCTProtocolParams import *
from DebugLog import DebugLog
import PIL
import copy
import traceback
import sys
import os
import tifffile
import struct
import time
import pickle
import multiprocessing as mproc
import queue 
from OCTProtocolParams import *
import scipy.signal


# Mscan processing options
class MscanProcOpts:
    def __init__(self):
        self.FFTZeroPadFactor = 5
        self.singlePt_zROI_indices = [0]        # for single point or points list mscans the Z point and spread at which to 
        self.singlePt_zROI_spread = 0           #    calculate TD and FFT (looking for max SNR in this region )
                                                  #    singlePt_zROIidx could be array in cae of points list scan types
        self.noiseBandwidth = 100            # bandwidth for noise calculations in Hz
        self.spikeRemoval = False            # remove spikes from time domain signal  
        self.spikeThreshold_nm = 50            # remove spikes from time domain signal  
        # self.FFTindex = 0
        
        self.zRes = 0                    # z resolution in microns (this 0 value will be changed to what is imported from oct.txt)
        self.correctAspectRatio = True      # whether to correct the aspect ratio in the image
        
        self.logAline = True
        self.bscanNormLow = 20               # range to normalize bscan after taking 20*log10 of magnitude
        self.bscanNormHigh = 150
        self.refractiveIndex = 0            # changed based on oct.txt
        self.centerWavelength = 0          # center wavelength in nm; changed based on oct.txt
        # [1.0590e-06 -1.0900e-10 4.4570e-15 -1.9210e-18]
      
      
class MScanData:
    def __init__(self):
        self.avgAline = None                # average of aline data
        self.posLenStep = None              # scan pos step
        self.posWidthStep = None            # scan width step
        self.stimFreqIdx = None             # sound stimulation frequency index in AudioOutputParams frequency arrray
        self.stimAmpIdx = None              # sound stimulation amplitude index in AudioOutputParams frequency arrray
        self.trialTime = None               # time of trial
        self.mscanTD = None                 # the mscan in time domain, at desired point (1D array)
        self.mscanFFT = None                # FFT of time doamin, at desired point (1D array)
        self.maxFFTFreq = None
        self.mscanSampleRate = None         # effective sample rate of mscan 
        self.stimRespMag = None             # magnitude resposne at stimulus frequency
        self.stimRespPhase = None           # phase resposne at stimulus frequency    
        
        self.TDnoiseSD = None                 # time domain (TD) noise 
        self.FDnoiseAboveStimFreqMean = None  # noise calculations for FFT noise (FD = Fourier domain) above and below stim frequency
        self.FDnoiseAboveStimFreqSD = None    #   Mean= mean noise, SD = standard deviation
        self.FDnoiseBelowStimFreqMean = None  #   num pts is deermined by noise bandiwth in mscan processing opts
        self.FDnoiseBelowStimFreqSD = None  
        self.frameNum = -1

# aggregate tuning curve at single point
class MscanTuningCurve:
    def __init__(self, audioParams):
        self.freq = audioParams.freq[0, :]
        self.amp = audioParams.amp             # stim amp array (db SPL)  
        numFreq = audioParams.getNumFrequencies()
        numAmp = len(self.amp)
        self.magResp = np.zeros((numFreq, numAmp))
        self.magResp[:, :] = np.NaN
        self.phaseResp = np.zeros((numFreq, numAmp))
        self.phaseResp[:, :] = np.NaN
        self.phaseRespUnwrapped = np.zeros((numFreq, numAmp))
        self.phaseRespUnwrapped[:, :] = np.NaN
        self.TDnoise = np.zeros((numFreq, numAmp))
        self.TDnoise[:, :] = np.NaN
        self.FDnoiseAboveStimFreqMean = np.zeros((numFreq, numAmp))  
        self.FDnoiseAboveStimFreqMean[:, :] = np.NaN
        self.FDnoiseAboveStimFreqSD = np.zeros((numFreq, numAmp))
        self.FDnoiseAboveStimFreqSD[:, :] = np.NaN
        self.FDnoiseBelowStimFreqMean = np.zeros((numFreq, numAmp))
        self.FDnoiseBelowStimFreqMean[:, :] = np.NaN
        self.FDnoiseBelowStimFreqSD = np.zeros((numFreq, numAmp))
        self.FDnoiseBelowStimFreqSD[:, :] = np.NaN
        
# aggregate data for mscan over a fairly region (B Mscan, Volume Mscan with possible ROI mask)
class MscanRegionData:
    def __init__(self, audioParams, scanParams, procOpts, numZPts):
        hsl = scanParams.length/2
        hsw = scanParams.width/2
        numFreq = audioParams.getNumFrequencies()
        numAmp = len(audioParams.amp)
        len_steps = scanParams.lengthSteps
        width_steps = scanParams.widthSteps
        
        self.posLen = np.linspace(-hsl, hsl, len_steps)     # 1D array of length positions
        self.posWidth = np.linspace(-hsw, hsw, width_steps)   # 1D array of width positions
        # numZPts = (procOpts.zROIend - procOpts.zROIstart + 1)        
        maxZ =procOpts.zRes*numZPts
        self.posDepth = np.linspace(0, maxZ, numZPts)              # 1D array of depth positions
        self.audioOutputParams = audioParams       # audio output data for this  Mscan
        self.scanParams = scanParams
        
        # 5D array of mscan magnitude data (pos z, pos x, pos y, freq, amplitude)        
        self.magResp = np.zeros((numZPts, width_steps, len_steps, numFreq, numAmp))  
        self.magResp[:, :, :, :, :] = np.NaN
        
         # 5D array of mscan phase data
        self.phaseResp = np.zeros((numZPts, width_steps, len_steps, numFreq, numAmp))  
        self.phaseResp[:, :, :, :, :] = np.NaN
        
        # noise calculations for FFT noise (FD = Fourier domain) above and below stim frequency
        #   Mean= mean noise, SD = standard deviation
        #   num pts is deermined by noise bandiwth in mscan processing opts
        self.FDnoiseAboveStimFreqMean = np.zeros((numZPts, width_steps, len_steps, numFreq, numAmp))    
        self.FDnoiseAboveStimFreqMean[:, :, :, :, :] = np.NaN
        
        self.FDnoiseAboveStimFreqSD = np.zeros((numZPts, width_steps, len_steps, numFreq, numAmp))  
        self.FDnoiseAboveStimFreqSD[:, :, :, :, :] = np.NaN
        
        self.FDnoiseBelowStimFreqMean = np.zeros((numZPts, width_steps, len_steps, numFreq, numAmp))  
        self.FDnoiseBelowStimFreqMean[:, :, :, :, :] = np.NaN
        
        self.FDnoiseBelowStimFreqSD = np.zeros((numZPts, width_steps, len_steps, numFreq, numAmp))  
        self.FDnoiseBelowStimFreqSD[:, :, :, :, :] = np.NaN
        
        # last valid indeixes 
        self.freqIdx = -1 
        self.ampIdx = -1 
        self.widthStep = -1
        self.lengthStep = -1
        
        self.xRes = np.NaN
        self.yRes = np.NaN
        self.zRes = np.NaN
        
        
# class reprenesenting the psoition, frequency and intensity of the mscan
class MscanPosAndStim:
    def __init__(self):
        self.posLenStep = 0    # the length position, in step #
        self.posWidthStep = 0  # the width position, in step #
        self.freqIdx = 0      # index of stimulus frequency in audio paras frequency array
        self.ampIdx = 0       # index of stimlus amplitude in audio paras amltide array
        self.stimFreq = 1e3   # stimulus frequency
        self.numFreq = 0  
        self.numAmp = 0 

# spke removal for  mscan
def spikeRemoval(mscanTD, spikeThresh):
    td_diff = np.diff(mscanTD)
    
    spike_idx = np.where(np.abs(td_diff) > spikeThresh)
    spike_idx2 = spike_idx[0][:] + 1
    spike_corr = np.zeros(len(mscanTD))
    spike_corr[spike_idx2] = -td_diff[spike_idx]
    spike_corr_cum = np.cumsum(spike_corr)
    mscanTD = mscanTD + spike_corr_cum
    
    return mscanTD
    
# process mscna data give
#    oct_data: 3D array of complex FFT indiexed by (trigger, zpos, trial)
def processMscanData(oct_data, mscanPosAndStim, scanParams, audioParams, procOpts, OCTtrigRate, mscanTuningCurveList, mscanRegionData, volData):
    scanP = scanParams

    mag = None
    phase = None
    rawdata = oct_data
    shp = rawdata.shape
    
    # rawdata = trigsPerPt = np.ceil(galv.settleTime * OCTtriggerRate)
    numTrigs = shp[0]
    numZPts = shp[1]
    numTrials = shp[2]
    mscan = MScanData()
    if DebugLog.isLogging:        
        DebugLog.log("Mscan.processMscanData(): numTrigs %d numTrials= %d OCTtrigRate= %g" % (numTrigs, numTrials, OCTtrigRate))

    t1 = time.time()
    if mag is None:
        mag = np.abs(rawdata)
        
    mag = np.clip(mag, 1, np.inf)
    mag = 20*np.log10(mag)
        
    if DebugLog.isLogging:        
        DebugLog.log("Mscan.processMscanData(): (before averaging) mag.shape= " + repr(mag.shape))          
    
    # average all trials
    mag = np.mean(mag, 2)
    # average all triggers
    mag = np.mean(mag, 0)
    avgAlineTime = time.time() - t1
    if DebugLog.isLogging:        
        DebugLog.log("Mscan.processMscanData(): (after averaging) avg aline time= %0.4f mag.shape= %s" % (avgAlineTime, repr(mag.shape))) 
        
    mscan.avgAline = mag
    mscan.posLenStep = mscanPosAndStim.posLenStep
    mscan.posWidthStep = mscanPosAndStim.posWidthStep
    
    freqIdx = mscanPosAndStim.freqIdx 
    ampIdx = mscanPosAndStim.ampIdx
    # if single pt mscan, then select pt of maximum SNR to perform FFT
    DebugLog.log("Mscan.processMscanData(): posLenStep= %d posWidthStep= %d freqIdx= %d ampIdxs= %d" % (mscan.posLenStep, mscan.posWidthStep, freqIdx, ampIdx))

    mscan.stimFreqIdx = freqIdx 
    mscan.stimAmpIdx = ampIdx
    trigRate = OCTtrigRate
    mscan.trialTime = numTrigs / trigRate
    mscan.mscanSampleRate = trigRate
    mscan.maxFFTFreq = trigRate/2

    zeroPad = procOpts.FFTZeroPadFactor
    numfftpts = int(2 ** np.round(np.log2(numTrigs*zeroPad)))
    numfftpts_2 = numfftpts // 2
    winfcn = np.hanning(numTrigs)
    
    stimFreq = mscanPosAndStim.stimFreq
    mscan.stimFreq = stimFreq
    
    rad_to_nm = procOpts.centerWavelength / (4 * np.pi * procOpts.refractiveIndex)
    win_magcorr = 2
    
    scanP = scanParams
    # initialize aggregate data if necesary
    if (scanP.lengthSteps == 1) and (scanP.widthSteps == 1) or scanP.pattern == ScanPattern.ptsList:
        if mscanTuningCurveList is None:   # initialize list of tuning curves
            # points list mscan
            audioP = audioParams
    
            tcurves = []
            numPts = len(scanP.ptsList)
            numPts = max((numPts, 1))
            for n in range(0, numPts):
                tcurves.append(MscanTuningCurve(audioP))
                
            mscanTuningCurveList = tcurves
    # volume or B mscan
    else:  
        if volData is None:
            mscanRegionData = MscanRegionData(audioParams, scanParams, procOpts, numZPts)
            volData = VolumeScan.VolumeData()
            volData.scanParams = scanParams
            volData.volumeImg = np.zeros((scanP.widthSteps, scanP.lengthSteps, numZPts))
            volData.xRes = 1e3*scanP.length/scanP.lengthSteps
            volData.yRes = 1e3*scanP.width/scanP.widthSteps
            volData.zRes = procOpts.zRes

    if ((scanP.lengthSteps == 1) and (scanP.widthSteps == 1)) or (scanP.pattern == ScanPattern.ptsList):
        zroi_indices = procOpts.singlePt_zROI_indices
        sprd = procOpts.singlePt_zROI_spread
        if len(zroi_indices) > 1:
            idx = zroi_indices[frameNum]
        else:
            idx = zroi_indices[0]

        # find the z point where the magnitude response is maximum (max SNR)
        idx1 = max(idx - sprd, 0)
        idx2 = min(idx + sprd, len(mag) - 1)
        
        if DebugLog.isLogging:        
            DebugLog.log("Mscan.processMscanData(): idx1= %d idx2= %d" % (idx1, idx2))
        
        maxSNRidx = np.argmax(mag[idx1:idx2])
        
        mscan.zROIidx = (idx1, idx2)
        mscan.maxSNR_Zidx = idx1 + maxSNRidx
        
        # calculate time domain phase
        if phase is None:
            dataRgn = rawdata[:, idx1 + maxSNRidx, :]
            ph = np.angle(dataRgn)
            if DebugLog.isLogging:        
                DebugLog.log("Mscan.processMscanData(): dataRgn.shape= " + repr(dataRgn.shape) + " ph.shape= " + repr(ph.shape) + " maxSNRidx= " + repr(maxSNRidx))
        else:
            ph = phase[:, idx1 + maxSNRidx, :]
            ph = ph / (2*np.pi)
            if DebugLog.isLogging:        
                DebugLog.log("Mscan.processMscanData(): ph.shape= " + repr(ph.shape) + " maxSNRidx= " + repr(maxSNRidx))
        
        
        # unwrap  
        #plt.plot(np.real(dataRgn[:, 0]), '-b', np.imag(dataRgn[:, 0]), '-r')
        #plt.show()
        
        unwrapThreshold = np.pi
        ph = np.unwrap(ph, discont=unwrapThreshold, axis=0)
        
        if DebugLog.isLogging:        
            DebugLog.log("Mscan.processMscanData(): ph.shape= " + repr(ph.shape))
        #plt.plot(ph)
        #plt.show()
        
        
        # remove DC/LF components by subtracting from mean 
        ph_mean = np.mean(ph, 0) 
        shp = ph.shape
        ph_mean = np.tile(ph_mean, (shp[0], 1))                
#        ph_mean = np.transpose(ph_mean)
        if DebugLog.isLogging:        
            DebugLog.log("Mscan.processMscanData(): ph.shape= %s ph_mean.shape= %s " % (repr(ph.shape), repr(ph_mean.shape)))

        ph = ph - ph_mean
        
        # average all trials
        ph = np.mean(ph, 1)
        mscan.TDnoiseSD = np.std(ph)
        
        # convert to nanometers
        ph = ph * rad_to_nm
        if DebugLog.isLogging:        
            DebugLog.log("Mscan.processMscanData(): ph.shape= " + repr(ph.shape))
        
        if procOpts.spikeRemoval:
            ph = spikeRemoval(ph, procOpts.spikeThreshold_nm)
        
        # take the FFT
        mscan_fft = np.fft.fft(ph * winfcn, numfftpts)
        mscan_fft = mscan_fft[0:numfftpts_2] * win_magcorr / numTrigs
        if DebugLog.isLogging:        
            DebugLog.log("Mscan.processMscanData(): numfftpts= %d mscan_fft.shape= %s " % (numfftpts, repr(mscan_fft.shape)))
        
        fftmag = np.abs(mscan_fft)
        fftphase = np.angle(mscan_fft)
        
        # find the magnitude and phase response at the stimulus frequency
        idx = np.floor(numfftpts*stimFreq/trigRate)
        i1 = idx - 2*zeroPad
        i1 = max(i1, 0)
        i2 = idx + 2*zeroPad
        i2 = min(i2, len(fftmag))
        
        if DebugLog.isLogging:        
            DebugLog.log("Mscan.processMscanData(): i1= %d i2= %d maxIdx= " % (i1, i2))            
            
        maxIdx = np.argmax(fftmag[i1:i2])
        
        if DebugLog.isLogging:        
            DebugLog.log("Mscan.processMscanData(): i1= %d i2= %d maxIdx= %d" % (i1, i2, maxIdx))
            
        mscan.stimRespMag = fftmag[i1 + maxIdx]
        mscan.stimRespPhase = fftphase[i1 + maxIdx]
        
        # calculate noise
        noisePts = np.ceil(numTrigs * procOpts.noiseBandwidth / trigRate)
        noisePts = noisePts*zeroPad
        noiseRgn = fftmag[i2:(i2 + noisePts)]
        mscan.FDnoiseAboveStimFreqMean = np.mean(noiseRgn)
        mscan.FDnoiseAboveStimFreqSD = np.std(noiseRgn)
        noiseRgn = fftmag[(i1 - noisePts):i1]
        mscan.FDnoiseBelowStimFreqMean = np.mean(noiseRgn)
        mscan.FDnoiseBelowStimFreqSD = np.std(noiseRgn)
        
        # fill in aggregate data
        scanPtNum = mscanPosAndStim.posLenStep
        tcurve = mscanTuningCurveList[scanPtNum]
        
        tcurve.magResp[freqIdx, ampIdx] = mscan.stimRespMag
        tcurve.phaseResp[freqIdx, ampIdx] = mscan.stimRespPhase
        tcurve.phaseRespUnwrapped[freqIdx, ampIdx] = mscan.stimRespPhase
        tcurve.phaseRespUnwrapped[:, ampIdx] = np.unwrap(tcurve.phaseRespUnwrapped[:, ampIdx], axis=0)
        
        tcurve.FDnoiseAboveStimFreqMean[freqIdx, ampIdx] = mscan.FDnoiseAboveStimFreqMean
        tcurve.FDnoiseAboveStimFreqSD[freqIdx, ampIdx] = mscan.FDnoiseAboveStimFreqSD
        tcurve.FDnoiseBelowStimFreqMean[freqIdx, ampIdx] = mscan.FDnoiseBelowStimFreqMean
        tcurve.FDnoiseBelowStimFreqSD[freqIdx, ampIdx] = mscan.FDnoiseBelowStimFreqSD
    else:  # B mscan or Volume M scan
        # calculate time domain phase
        t1 = time.time()
        if phase is None:
            ph = np.angle(rawdata)
            if DebugLog.isLogging:        
                DebugLog.log("Mscan.processMscanData(): rawdata.shape= " + repr(rawdata.shape) + " ph.shape= " + repr(ph.shape))
        else:
            ph = phase / (2 * np.pi)
            if DebugLog.isLogging:        
                DebugLog.log("Mscan.processMscanData(): ph.shape= " + repr(ph.shape))
        
        phaseCalcTime= time.time() - t1
        # unwrap  
        #plt.plot(np.real(dataRgn[:, 0]), '-b', np.imag(dataRgn[:, 0]), '-r')
        #plt.show()
        
        unwrapThreshold = np.pi
        t1 = time.time()
        ph = np.unwrap(ph, discont=unwrapThreshold, axis=0)
        unwrapTime = time.time() - t1
        if DebugLog.isLogging:        
            DebugLog.log("Mscan.processMscanData(): unwrap time=%0.4f ph.shape= %s" % (unwrapTime, repr(ph.shape)))
        #plt.plot(ph)
        #plt.show()
        
        # remove DC/LF components by subtracting from mean 
        t1 = time.time()
        ph_mean = np.mean(ph, 0) 
        if DebugLog.isLogging:        
            DebugLog.log("Mscan.processMscanData(): ph_mean.shape= " + repr(ph_mean.shape))
        ph_mean = np.tile(ph_mean, (numTrigs, 1, 1))
        if DebugLog.isLogging:        
            DebugLog.log("Mscan.processMscanData(): ph_mean.shape= " + repr(ph_mean.shape))
        ph = ph - ph_mean
        if DebugLog.isLogging:        
            DebugLog.log("Mscan.processMscanData(): ph.shape= " + repr(ph.shape))
        
        # average all trials
        ph = np.mean(ph, 2)
        phaseMeanSubTime = time.time() - t1
        
        mscan.TDnoiseSD = np.std(ph)
        
        # convert to nanometers
        ph = ph * rad_to_nm
        
        #if procOpts.mscan.spikeRemoval:
        #    ph = spikeRemoval(ph, procOpts.mscan.spikeThreshold_nm)
        
        # take the FFT
        t1 = time.time()
        winfcn = np.tile(winfcn, (numZPts, 1))
        winfcn = np.transpose(winfcn)
        win_ph = ph * winfcn 
        winFcnTime = time.time() - t1
        
        if DebugLog.isLogging:        
            DebugLog.log("Mscan.processMscanData(): ph.shape= %s winfcn.shape= %s win_ph.dtype= %s" % (repr(ph.shape), repr(winfcn.shape), repr(win_ph.dtype)))
        
        
        #mscan_fft = np.fft.fft(win_ph, numfftpts, 0)
        #mscan_fft = mscan_fft[0:numfftpts_2, :] * win_magcorr / numTrigs

        # preallocate result array, this makes a HUGE speed diffrence when using multiprocesing
        fft_output_pts = numfftpts // 2
        if numfftpts % 2 == 0:
            fft_output_pts += 1
        mscan_fft = np.zeros((fft_output_pts, win_ph.shape[1]), dtype=np.complex)
        t1 = time.time()
        mscan_fft[:, :] = np.fft.rfft(win_ph, n=numfftpts, axis=0)
        
        fftTime = time.time() - t1            
        
        
        if DebugLog.isLogging:        
            DebugLog.log("Mscan.processMscanData(): fft time= %0.4f numfftpts= %d mscan_fft.shape= %s " % (fftTime, numfftpts, repr(mscan_fft.shape)))
            
        mscan_fft = mscan_fft * win_magcorr / numTrigs
        
        t1 = time.time()
        #fftmag = np.abs(mscan_fft)
        #fftphase = np.angle(mscan_fft)
        magPhaseCalcTime = time.time() - t1            
        
        # find the magnitude and phase response at the stimulus frequency
        t1 = time.time()
        idx = np.floor(numfftpts*stimFreq/trigRate)
        i1 = idx - 2*zeroPad
        i1 = max(i1, 0)
        i2 = idx + 2*zeroPad
        #i2 = min(i2, fftmag.shape[0])
        i2 = min(i2, mscan_fft.shape[0])
        
        if DebugLog.isLogging:        
            DebugLog.log("Mscan.processMscanData(): i1= %d i2= %d " % (i1, i2))            
        numZ = mscan_fft.shape[1]
        mscan.stimRespMag = np.zeros(numZ)
        mscan.stimRespPhase = np.zeros(numZ)
        
        for n in range(0, numZ):
            fftmag = np.abs(mscan_fft[i1:i2, n])
            #maxIdx = np.argmax(fftmag[i1:i2, n], 0)
            maxIdx = np.argmax(fftmag, 0)
            # DebugLog.log("MscanProtocol.processData(): i1= %d i2= %d maxIdx= %d" % (i1, i2, maxIdx))
            mscan.stimRespMag[n] = fftmag[maxIdx]
            mscan.stimRespPhase[n] = np.angle(mscan_fft[i1+maxIdx, n])
        
        maxMagCalcTime = time.time() - t1        
        if DebugLog.isLogging:        
            DebugLog.log("Mscan.processMscanData(): mscan.stimRespMag.shape= %s time= %0.4f" % (repr(mscan.stimRespMag.shape), maxMagCalcTime))
        
        # calculate noise
        t1 = time.time()
        noisePts = np.ceil(numTrigs * procOpts.noiseBandwidth / trigRate)
        noisePts = noisePts*zeroPad
        #noiseRgn = fftmag[i2:(i2 + noisePts), :]
        noiseRgn = np.abs(mscan_fft[i2:(i2 + noisePts), :])
        mscan.FDnoiseAboveStimFreqMean = np.mean(noiseRgn, 0)
        mscan.FDnoiseAboveStimFreqSD = np.std(noiseRgn, 0)
        #noiseRgn = fftmag[(i1 - noisePts):i1, :]
        noiseRgn = np.abs(mscan_fft[(i1 - noisePts):i1, :])
        mscan.FDnoiseBelowStimFreqMean = np.mean(noiseRgn, 0)
        mscan.FDnoiseBelowStimFreqSD = np.std(noiseRgn, 0)
        noiseCalcTime = time.time() - t1 
        if DebugLog.isLogging:        
            DebugLog.log("Mscan.processMscanData(): mscan.stimRespMag.FDnoiseAboveStimFreqMean= %s" % (repr(mscan.FDnoiseAboveStimFreqMean.shape)))
        
        t1 = time.time()
        rgnData = mscanRegionData
        # rgnData = self.aggregateData.mscanRegionData
        posWidthStep = mscan.posWidthStep
        posLenStep = mscan.posLenStep
        freqIdx = mscan.stimFreqIdx
        ampIdx = mscan.stimAmpIdx
        rgnData.magResp[:, posWidthStep, posLenStep, freqIdx, ampIdx] = mscan.stimRespMag
        rgnData.phaseResp[:, posWidthStep, posLenStep, freqIdx, ampIdx] = mscan.stimRespPhase
        
        rgnData.FDnoiseAboveStimFreqMean[:, posWidthStep, posLenStep, freqIdx, ampIdx] = mscan.FDnoiseAboveStimFreqMean
        rgnData.FDnoiseAboveStimFreqSD[:, posWidthStep, posLenStep, freqIdx, ampIdx] = mscan.FDnoiseAboveStimFreqSD
        rgnData.FDnoiseBelowStimFreqMean[:, posWidthStep, posLenStep, freqIdx, ampIdx] = mscan.FDnoiseBelowStimFreqMean
        rgnData.FDnoiseBelowStimFreqSD[:, posWidthStep, posLenStep, freqIdx, ampIdx] = mscan.FDnoiseBelowStimFreqSD
        
        rgnData.freqIdx = freqIdx
        rgnData.ampIdx = ampIdx
        rgnData.widthStep = posWidthStep
        rgnData.lengthStep = posLenStep
        
        rgnData.xRes = 1e3*scanP.length/scanP.lengthSteps
        rgnData.yRes = 1e3*scanP.width/scanP.widthSteps
        rgnData.zRes = procOpts.zRes
        
        
        nL = procOpts.bscanNormLow
        nH = procOpts.bscanNormHigh
        nRng = nH - nL
        dataAssignTime = time.time() - t1
        if DebugLog.isLogging:        
            DebugLog.log("Mscan.processMscanData(): nL= %g nH=%g nRng= %g" % (nL, nH, nRng))
        
        t1 = time.time()
        # remap range to 0...1
        aline = (mscan.avgAline - nL)/nRng  
    
        # remap range to 0 ... to 2^16 - 1
        aline16b = aline*65535
        aline16b = np.clip(aline16b, 0, 65535)
        aline16b = np.require(aline16b, 'uint16')
        volData.volumeImg[posWidthStep, posLenStep, :] = aline16b
        alineCalcTime = time.time() - t1

        if DebugLog.isLogging:        
            DebugLog.log("Mscan.processMscanData(): \n\tavgAlineTime= %0.4f \n\tphaseCalcTime= %0.4f \n\tphaseMeanSubTime= %0.4f" % (avgAlineTime, phaseCalcTime, phaseMeanSubTime))
            DebugLog.log("\tunwrapTime= %0.4f winFcnTime= %0.4f \n\tFFTtime= %0.4f  \n\tmagPhaseCalcTime= %0.4f" % (unwrapTime, winFcnTime, fftTime, magPhaseCalcTime))
            DebugLog.log("\tmaxMagCalcTime= %0.4f noiseCalcTime= %0.4f \n\tdataAssignTime= %0.4f \n\talineCalcTime= %0.4f" % (maxMagCalcTime, noiseCalcTime, dataAssignTime, alineCalcTime))
        
    mscan.mscanTD = ph
    mscan.mscanFFT = mscan_fft
    
    return mscan, mscanTuningCurveList, mscanRegionData, volData


def makeVolMScanBoxROIScanParams(appObj):
    xyScanP = appObj.volDataLast.scanParams
    DebugLog.log("makeVolMScanBoxROIScanParams: xyScanP= %s" % (repr(xyScanP)))
    
    ul = appObj.vol_plane_proj_gv.ROIBox_pt1  # upper left point
    lr = appObj.vol_plane_proj_gv.ROIBox_pt2  # lower right point
    # exchange box points if necessary
    if ul[1] > lr[1]:
        tmp = ul
        ul = lr
        lr = tmp
        
    DebugLog.log("makeVolMScanBoxROIScanParams: ul= %s lr= %s" % (repr(ul), repr(lr)))
    
    if ul is not None and lr is not None:
        roi_dw = np.abs(lr[0] - ul[0])   # roi delta width
        roi_dh = np.abs(lr[1] - ul[1])   # roi delta height
        roi_ow = ul[0] + (roi_dw // 2)
        roi_oh = ul[1] + (roi_dh // 2)
    else: 
        return None

    (img_w, img_h) = appObj.vol_plane_proj_gv.getImageWidthHeight()
    DebugLog.log("makeVolMScanBoxROIScanParams: img_w= %d img_h= %d" % (img_w, img_h))
    scanP = ScanParams()
    scanP.length = roi_dw * xyScanP.length / img_w 
    scanP.width = roi_dh * xyScanP.width / img_h
    scanP.lengthOffset = xyScanP.lengthOffset + roi_ow * xyScanP.length / img_w
    scanP.widthOffset = xyScanP.widthOffset + roi_oh * xyScanP.width / img_h
    scanP.rotation_Z = xyScanP.rotation_Z
    
    xyRes = appObj.volMscan_sampling_dblSpinBox.value()*1e-3
    scanP.lengthSteps = int(round(scanP.length / xyRes))
    scanP.widthSteps = int(round(scanP.width / xyRes))
    DebugLog.log("makeVolMScanBoxROIScanParams: lengthSteps= %d widthSteps= %d" % (scanP.lengthSteps, scanP.widthSteps))
    xyScanP.continuousScan = False
    
    return scanP

def makeVolMScanPolyROIScanParams(appObj, poly = None):
    scanP = ScanParams()
    
    if poly is None:
        poly = appObj.vol_plane_proj_gv.ROI_poly
    
    rect = poly.boundingRect()  # get the bounding rectangle of the polygon

    ul = rect.topLeft()
    lr = rect.bottomRight()
    # exchange box points if necessary
    if ul.y() > lr.y():
        tmp = ul
        ul = lr
        lr = tmp
        
    if ul is not None and lr is not None:
        roi_dw = np.abs(lr.x() - ul.x())   # roi delta width
        roi_dh = np.abs(lr.y() - ul.y())   # roi delta height
        roi_ow = ul.x() + (roi_dw // 2)
        roi_oh = ul.y() + (roi_dh // 2)
    else:
        return None

    (img_w, img_h) = appObj.vol_plane_proj_gv.getImageWidthHeight()
    scanP = ScanParams()
    scanP.length = roi_dw * xyScanP.length / img_w 
    scanP.width = roi_dh * xyScanP.width / img_h 
    scanP.lengthOffset = xyScanP.lengthOffset + roi_ow * xyScanP.length / img_w
    scanP.widthOffset = xyScanP.widthOffset + roi_oh * xyScanP.width / img_h
    scanP.rotation_Z = xyScanP.rotation_Z
    
    xyRes = appObj.volMscan_sampling_dblSpinBox.value()*1e-3
    len_steps = int(round(scanP.length / xyRes))
    width_steps = int(round(scanP.width / xyRes))
    scanP.lengthSteps = len_steps 
    scanP.widthSteps = width_steps
    DebugLog.log("makeVolMScanBoxROIScanParams: lengthSteps= %d widthSteps= %d" % (scanP.lengthSteps, scanP.widthSteps))
    xyScanP.continuousScan = False
    
    # create ROI mask 
    shp = (len_steps, width_steps)
    scanP.boxROIMaskXY = np.zeros(shp, np.bool)
    for l in range(0, len_steps):
        for w in range(0, width_steps):
            x = roi_dw * l / len_steps + ul.x()
            y = roi_dh * w / len_steps + ul.y()
            ptf = QtCore.QPointF(x, y)
            scanP.boxROIMaskXY[l, w] = poly.contains(ptf)
    
    return scanP    

def makeVolMScanFreeROIScanParams(appObj):
    # since the free draw ROI also consists of a polygon, we can use same method
    # as the polygon ROI
    poly = appObj.vol_plane_proj_gv.freedraw_poly
    scanP = makeVolMScanPolyROIScanParams(appObj, poly)
    
    return scanP
    
def makeMscanScanParamsAndZROI(appObj):
    scanP = ScanParams()
    
    zROIIndices = []
    zROIspread = appObj.mscan_zROIlspread_spinBox.value()
    
    roiBegin = -1
    roiEnd = -1

    if appObj.mscan_single_pt_button.isChecked():
        scanP.length = 0
        scanP.width = 0 
        scanP.lengthSteps = 1
        scanP.widthSteps = 1 
        
        imgScanP = appObj.imgDataScanParams
        DebugLog.log ("makeMscanScanParamsAndZROI imgScanP= %s" % (repr(imgScanP)))

        # pt = self.bscan_img_gv.ptsList[0]
        pt = appObj.bscan_img_gv.singlePt
        
        img_w = appObj.imgdata_8b.shape[1]
        img_h = appObj.imgdata_8b.shape[0]
        DebugLog.log ("makeMscanScanParamsAndZROI pt= %s  img_w=%d img_h= %d " % (repr(pt), img_w, img_h))
        scanP.lengthOffset = (pt[0] - (img_w / 2) ) * imgScanP.length / img_w  + imgScanP.lengthOffset
        scanP.widthOffset = imgScanP.widthOffset
        scanP.rotation_Z = imgScanP.rotation_Z
        
        roiBegin = appObj.imgdata_zROI[0]
        roiEnd = appObj.imgdata_zROI[1]
        
        zROIIndices.append(pt[1])
    elif appObj.bmscan_box_region_button.isChecked():
        imgScanP = appObj.imgDataScanParams
        DebugLog.log ("makeMscanScanParamsAndZROI imgScanP= %s" % (repr(imgScanP)))
        #pt1 = self.bscan_img_gv.ROIBox_pt1
        #pt2 = self.bscan_img_gv.ROIBox_pt2
        (pt1, pt2) = appObj.bmscan_box_rgn
        
        # calculate (x1, y1) as upper left corner and (x2, y2) as lower right corner (with (0,0) being top left corner)
        x1 = min(pt1[0], pt2[0])
        x2 = max(pt1[0], pt2[0])
        y1 = min(pt1[1], pt2[1])
        y2 = max(pt1[1], pt2[1])
        
        img_w = appObj.imgdata_8b.shape[1]            
        img_h = appObj.imgdata_8b.shape[0]
        
        DebugLog.log ("makeMscanScanParamsAndZROI pt1= %s pt2= %s img_w=%d img_h= %d " % (repr(pt1), repr(pt2), img_w, img_h))
        scanP.width = 0 
        # scanP.lengthSteps = int(x2 - x1)
        
        scanP.lengthSteps = appObj.BMscan_numSteps_spinBox.value()
        scanP.widthSteps = 1 

        dw = np.abs(x2 - x1)
        scanP.length = imgScanP.length*dw/img_w
        x_mid = (x1 + x2)/2
        scanP.lengthOffset = (x_mid - (img_w / 2) ) * imgScanP.length / img_w  + imgScanP.lengthOffset
        scanP.widthOffset = imgScanP.widthOffset
        scanP.rotation_Z = imgScanP.rotation_Z 
        
        roiBegin = int(y1 + appObj.imgdata_zROI[0])
        roiEnd = int(y2 + appObj.imgdata_zROI[0])
        roiSize = roiEnd - roiBegin 
        # ensure ROI is aligned to multiple of 4
        # this is required if data is 16-bit 
        if roiSize % 4 != 0:
            roiEnd += 4 - roiSize % 4
    else:
        volBox = appObj.vol_boxROI_pushButton.isChecked()
        volPoly = appObj.vol_polyROI_pushButton.isChecked()
        volFree = appObj.vol_freeROI_pushButton.isChecked()
        
        if volBox:
            scanP = makeVolMScanBoxROIScanParams(appObj)
        elif volPoly:
            scanP = makeVolMScanPolyROIScanParams(appObj)
        elif volFree:
            scanP = makeVolMScanFreeROIScanParams(appObj)
        else:
            scanP = None  # this case indicates that no valid region has been selected
        
        roiBegin = appObj.imgdata_zROI[0]
        roiEnd = appObj.imgdata_zROI[1]
        roiSize = roiEnd - roiBegin 
        # ensure ROI is aligned to multiple of 4
        # this is required if data is 16-bit 
        if roiSize % 4 != 0:
            roiEnd += 4 - roiSize % 4
            
    DebugLog.log ("makeMscanScanParamsAndZROI scanP= %s" % (repr(scanP)))
    DebugLog.log ("makeMscanScanParamsAndZROI roiBegin= %d roiEnd= %d zROIindices= %s" % (roiBegin, roiEnd, repr(zROIIndices)))
    return (scanP, roiBegin, roiEnd, zROIIndices, zROIspread)
    
def getXYPos(lenStep, widthStep, scanParams):
    scanP = scanParams
    if scanP.pattern == ScanPattern.ptsList:
        pt = scanP.ptsList[lenStep]
        xPos = pt[0]
        yPos = pt[1]
    else:
        rotRad = np.pi*scanP.rotation_Z/180
        # rotation effect on the length offset
        cos_rot = np.cos(rotRad)
        sin_rot = np.sin(rotRad)
        
        # rotation effect on the width offset
        wo_rot = rotRad + np.pi/2
        wo_cos_rot = np.cos(wo_rot)
        wo_sin_rot = np.sin(wo_rot)
        
        hsl = scanP.length / 2       # half scan length
        hsw = scanP.width / 2        # half scan width
       
        # where scan is along the length before accounting for offset and rotation
        len_offset = -hsl + lenStep * scanP.length / scanP.lengthSteps 
        # x and y offsets contributing by the length step
        len_xOffset =  len_offset*cos_rot + scanP.lengthOffset*cos_rot
        len_yOffset =  len_offset*sin_rot + scanP.lengthOffset*sin_rot
        
        # where scan is along the widthh before accounting for offset and rotation
        width_offset = -hsw + widthStep * scanP.width / scanP.widthSteps
        # x and y offsets contributing by the width step
        width_xOffset = width_offset*wo_cos_rot
        width_yOffset = width_offset*wo_sin_rot
        
        xPos = len_xOffset + width_xOffset 
        yPos = len_yOffset + width_yOffset 
        
    return (xPos, yPos)


def makeAudioOutput(audioParams, audioHW, spkNum, f, a):
    outputRate = audioHW.DAQOutputRate
    trialDur = 1e-3*audioParams.getTrialDuration(a)
    trialPts = np.ceil(trialDur * outputRate)
    
    spkOut = []
    DebugLog.log("makeAudioOutput: spkNum=%d trialDur=%f trialPts=%d stimType= %s" % (spkNum, trialDur, trialPts, audioParams.stimType))
    if spkNum == 1 and (audioParams.stimType == AudioStimType.TONE_LASER):
        spkOut = np.zeros((trialPts))    
        laserStimDur = audioParams.stimDuration*1e-3
        stimPts = int(np.ceil(laserStimDur * outputRate ))
        DebugLog.log("makeAudioOutput: laserStimDurr=%f stimPts=%d" % (laserStimDur, stimPts))
        i1 = int(np.ceil(outputRate*1e-3*audioParams.stimOffset))
        i2 = i1 + stimPts
        spkOut[i1:i2] = 5
        attenLvl = 0
    else:
        (outV, attenLvl) = audioHW.getCalibratedOutputVoltageAndAttenLevel(f, a, spkNum)
        
        if(outV > 0):
            stimDur = 1e-3*audioParams.stimDuration
            stimOffset = 1e-3*audioParams.stimOffset
            
            if audioParams.stimType == AudioStimType.TONE_LASER:
                stimDur = trialDur
                stimOffset = 0
                
            stimEnv = 1e-3*audioParams.stimEnvelope
            offsetPts = np.ceil(stimOffset * outputRate)
            stimPts = np.ceil(stimDur * outputRate)
            
             # in the case that stim + offset will excdeed trial duration, we must trim the stim 
            if (stimPts + offsetPts) > trialPts:  
                stimPts = trialPts - offsetPts
    
            envPts = np.ceil(stimEnv * outputRate)
            spkOut = np.zeros((trialPts))
            t = np.linspace(0, stimDur, stimPts)
            sig = outV*np.sin(2*np.pi*1000*f*t)
            envFcn = np.ones((stimPts))
            envFcn[0:envPts] = np.linspace(0, 1, envPts)
            envFcn[stimPts-envPts:] = np.linspace(1, 0, envPts)
            sig = sig*envFcn
            spkOut[offsetPts:offsetPts+stimPts] = sig

            
            
    return (spkOut, attenLvl)
        
        
def makeCSVTableString(data, dataFmt, colHdr, colHdrFmt, rowHdr, rowHdrFmt):
    s = ''
    fmt = ', ' + colHdrFmt
    # write teh column headers
    for n in range(0, len(colHdr)):
        s = s + fmt % (colHdr[n])
        
    s = s + '\n'
    for i in range(0, len(rowHdr)):
        s = s + rowHdrFmt % (rowHdr[i])
        for k in range(0, len(colHdr)):
            fmt = ', ' + dataFmt
            s = s + fmt % (data[i, k])
            
        s = s + '\n'
        
    return s
    
# save the time domain mscan data
def saveMscanDataTimeDomain(mscanData, freq, amp, frameNum, saveDir):
    fileName = 'Mscan TD Phase %0.2f kHz %0.1f dB frame %d' % (freq, amp, frameNum)
    filePath = os.path.join(saveDir, fileName)
    f = open(filePath, 'wb')
    f.write(struct.pack('d', mscanData.trialTime))
    f.write(struct.pack('I', mscanData.mscanTD.shape[0]))
    f.write(struct.pack('%sd' % (len(mscanData.mscanTD)), *mscanData.mscanTD))
    
    # pickle.dump(spCal, f)
    f.close()

def writeExcelFreqAmpHeader(ws, freq, amp, row=0, col=1):
    # write amplitude header in first row
    r = row
    c = col
    for a in amp:
        c += 1
        ws.write(r, c, a)


    #write frequency header down col
    c = col
    for f in freq:
        r += 1
        ws.write(r, c, f)
        
def writeExcel2DData(ws, data, row=1, col=1):
    for r in range(0, data.shape[0]):
        for c in range(0, data.shape[1]):
            ws.write(r+row, c+col, data[r, c])

"""
    saveMscanTuningExcel
        save mscan tuning curve data as Excel workbook (XLSX format)
"""
        
def saveMscanTuningCurveExcel(mscanTuningCurve, audioParams, ptNum, saveDir):
    try: 
        import xlsxwriter
    except:
        DebugLog.log("could not import xlsxwriter, not saving in Excel")
        return
        
    fileName = 'Mscan Tuning PtNum %d.xlsx' % (ptNum)
    filepath = os.path.join(saveDir, fileName)
    tcurve = mscanTuningCurve
    
    amp = audioParams.amp
    freq = audioParams.freq[0, :]
    numFreq = audioParams.getNumFrequencies()
    numAmp = len(amp)

    row = 0
    col = 0
    workbook = xlsxwriter.Workbook(filepath)
    ws = workbook.add_worksheet('Mag')
    ws.write(row, col, 'Stim Mag Resp')
    writeExcelFreqAmpHeader(ws, freq, amp, row+1, col)
    writeExcel2DData(ws, tcurve.magResp, row+2, col+1)
    
    chart = workbook.add_chart({'type' :'scatter', 'subtype' : 'straight_with_markers'})
    catStr = '=Mag!$A$3:$A$' + str(numFreq+2)
    # Configure the chart. In simplest case we add one or more data series.
    for a_i in range(0, numAmp):
        colStr = chr(ord('B') + a_i)
        rngStr = colStr + '3:' + colStr + str(2+numFreq)
        valStr = '=Mag!' + rngStr
        print("saveMscanTuningCurveExcel(): valStr= ", valStr)
        chart.add_series({'values': valStr, 'name' : str(amp[a_i]), 'categories': catStr})
    
    chart.set_title({'name' : 'Magnitude Response at Stim Frequency'})
    chart.set_x_axis({'name' : 'Stim Frequency (kHz)'})
    chart.set_y_axis({'name' : 'Response (nm)'})
    
    # create the location string where chart will be located
    locStr = chr(ord('B') + numAmp + 1) + '1'
    # Insert the chart into the worksheet.
    ws.insert_chart(locStr, chart)

    ws = workbook.add_worksheet('Phase')
    ws.write(row, col, 'Phase Resp Unwrapped')
    writeExcelFreqAmpHeader(ws, freq, amp, row+1, col)
    writeExcel2DData(ws, tcurve.phaseRespUnwrapped, row+2, col+1)

    row = row + numFreq + 3
    ws.write(row, col, 'Phase Resp Raw')
    writeExcelFreqAmpHeader(ws, freq, amp, row+1, col)
    writeExcel2DData(ws, tcurve.phaseResp, row+2, col+1)
    
    ws = workbook.add_worksheet('Mag Noise')
    row = 0
    ws.write(row, col, 'Mean Above Stim Freq')
    writeExcelFreqAmpHeader(ws, freq, amp, row+1, col)
    writeExcel2DData(ws, tcurve.phaseRespUnwrapped, row+2, col+1)
    
    row = row + numFreq + 3
    ws.write(row, col, 'St. dev. Above Stim Freq')
    writeExcelFreqAmpHeader(ws, freq, amp, row+1, col)
    writeExcel2DData(ws, tcurve.phaseRespUnwrapped, row+2, col+1)
 
    col = numAmp + 2
    row = 0
    ws.write(row, col, 'Mean Below Stim Freq')
    writeExcelFreqAmpHeader(ws, freq, amp, row+1, col)
    writeExcel2DData(ws, tcurve.phaseRespUnwrapped, row+2, col+1)
    
    row = row + numFreq + 3
    ws.write(row, col, 'St. dev. Below Stim Freq')
    writeExcelFreqAmpHeader(ws, freq, amp, row+1, col)
    writeExcel2DData(ws, tcurve.phaseRespUnwrapped, row+2, col+1)
    
    workbook.close()

"""
    saveMscanTuningCurve
        save mscan tuning curve data as CSV and Excel if import is successful
"""
def saveMscanTuningCurve(mscanTuningCurve, audioParams, ptNum, saveDir):
    tcurve = mscanTuningCurve
    fileName = 'Mscan Tuning PtNum %d.csv' % (ptNum)
    filePath = os.path.join(saveDir, fileName)
    f = open(filePath, 'w') 
    amp = audioParams.amp
    freq = audioParams.freq[0, :]

    datafmt = '%0.4f'
    rowfmt = '%0.2f'
    colfmt = '%0.1f'
    f.write('Stim Mag Resp\n')
    f.write(makeCSVTableString(tcurve.magResp, datafmt, amp, colfmt, freq, rowfmt))

    f.write('\nStim Phase Resp Unwrapped\n')    
    f.write(makeCSVTableString(tcurve.phaseRespUnwrapped, datafmt, amp, colfmt, freq, rowfmt))
        
    f.write('\nStim Phase Resp Raw\n')    
    f.write(makeCSVTableString(tcurve.phaseResp, datafmt, amp, colfmt, freq, rowfmt))
        
    f.write('\n')     # write new line
    f.close()
    
    fileName = 'Mscan Noise PtNum %d.csv' % (ptNum)
    filePath = os.path.join(saveDir, fileName)
    f = open(filePath, 'w')
    
    f.write('FFT Noise Above Stim Freq Mean\n')
    f.write(makeCSVTableString(tcurve.FDnoiseAboveStimFreqMean, datafmt, amp, colfmt, freq, rowfmt))
    
    f.write('\nFFT Noise Above Stim Freq SD\n')
    f.write(makeCSVTableString(tcurve.FDnoiseAboveStimFreqSD, datafmt, amp, colfmt, freq, rowfmt))
    
    f.write('\nFFT Noise Below Stim Freq Mean\n')
    f.write(makeCSVTableString(tcurve.FDnoiseBelowStimFreqMean, datafmt, amp, colfmt, freq, rowfmt))
    
    f.write('\nFFT Noise Below Stim Freq SD\n')
    f.write(makeCSVTableString(tcurve.FDnoiseBelowStimFreqSD, datafmt, amp, colfmt, freq, rowfmt))
    
    f.close()
    
    try: 
        saveMscanTuningCurveExcel(mscanTuningCurve, audioParams, ptNum, saveDir)
    except:
        DebugLog.log("exception attempting to save as Excel file")
        traceback.print_exc()
    
        
# save vibratorery response images
def saveMscanRegionData(mscanRegionData, volData, saveDir):
    volImg = volData.volumeImg
    
    shp = volImg.shape
    fileName = 'Avg Aline'
    filePath = os.path.join(saveDir, fileName)
    numpts = np.prod(shp)
    DebugLog.log("Mscan.saveProcessedData() volume numpts= %s shp=%s" % (repr(numpts), repr(shp)))
    f = open(filePath, 'wb')
    f.write(struct.pack('III', shp[0], shp[1], shp[2]))
    fmt_str = '%dH' % numpts
    DebugLog.log("Mscan.saveProcessedData() volume numpts= %s shp=%s fmt_str= %s" % (repr(numpts), repr(shp), repr(fmt_str)))
    volImg = volImg.reshape(numpts)
    volImg = np.require(volImg, np.uint16)
    b = struct.pack(fmt_str, *volImg)
    f.write(b)
    f.close()
            
    # save vibrotory response data
    rgnData = mscanRegionData
    fileName = 'Mscan Stim Freq Resp Mag'
    filePath = os.path.join(saveDir, fileName)
    f = open(filePath, 'wb')
    shp = rgnData.magResp.shape
    numpts = np.prod(shp)
    DebugLog.log("MscanProtocol.saveProcessedData() mscanRegionData magResp.shp=%s numpts= %s" % (repr(shp), repr(numpts)))
    f.write(struct.pack('IIIII', shp[0], shp[1], shp[2], shp[3], shp[4]))
    f.write(struct.pack('%df' % numpts, *rgnData.magResp.reshape(numpts)))
    f.close()
    
    fileName = 'Mscan Stim Freq Resp Phase'
    filePath = os.path.join(saveDir, fileName)
    f = open(filePath, 'wb')
    f.write(struct.pack('IIIII', shp[0], shp[1], shp[2], shp[3], shp[4]))
    f.write(struct.pack('%df' % numpts, *rgnData.phaseResp.reshape(numpts)))
    f.close()
    
    fileName = 'Mscan FD Noise Above Stim Mean'
    filePath = os.path.join(saveDir, fileName)
    f = open(filePath, 'wb')
    f.write(struct.pack('IIIII', shp[0], shp[1], shp[2], shp[3], shp[4]))
    f.write(struct.pack('%df' % numpts, *rgnData.FDnoiseAboveStimFreqMean.reshape(numpts)))
    f.close()
    
    fileName = 'Mscan FD Noise Above Stim SD'
    filePath = os.path.join(saveDir, fileName)
    f = open(filePath, 'wb')
    f.write(struct.pack('IIIII', shp[0], shp[1], shp[2], shp[3], shp[4]))
    f.write(struct.pack('%df' % numpts, *rgnData.FDnoiseAboveStimFreqSD.reshape(numpts)))
    f.close()
    
    fileName = 'Mscan FD Noise Below Stim Mean'
    filePath = os.path.join(saveDir, fileName)
    f = open(filePath, 'wb')
    f.write(struct.pack('IIIII', shp[0], shp[1], shp[2], shp[3], shp[4]))
    f.write(struct.pack('%df' % numpts, *rgnData.FDnoiseBelowStimFreqMean.reshape(numpts)))
    f.close()
    
    fileName = 'Mscan FD Noise Below Stim SD'
    filePath = os.path.join(saveDir, fileName)
    f = open(filePath, 'wb')
    f.write(struct.pack('IIIII', shp[0], shp[1], shp[2], shp[3], shp[4]))
    f.write(struct.pack('%df' % numpts, *rgnData.FDnoiseBelowStimFreqSD.reshape(numpts)))
    f.close()
    

class VibDataType(Enum):
    STIM = 0
    F1 = 1
    DP = 2
    
class VibNoiseRegion(Enum):
    ABOVE = 0
    BELOW = 1
    
# make an indiexed vibratory response (magnitude and phase) images overlayed on top of the reflectivity
# only the pixels that above given threshold will 
# this is accomplished by splitting the colormap into two halves
# index values 0...127 represent the image (relection intenisty), while values 128 to 255 represent vibration
# an appropriate color map that has      

# returns a tuple (magImg, phaseImg, ,minMag, maxMag)
def makeVibRespImg(mscanRgnData, volData, magAutoNorm=True, magNorms=(0, 100), intThresh=0, vibDataType=VibDataType.STIM, vibNoiseRgn=VibNoiseRegion.ABOVE, noiseNumSD=3):
    widthStep = mscanRgnData.widthStep
    freqIdx= mscanRgnData.freqIdx
    ampIdx = mscanRgnData.ampIdx
    
    # find the indices where the intensity is below threshold
    intImg = volData.volumeImg[widthStep, :, :]
    intImg = intImg.transpose()
    if DebugLog.isLogging:
        DebugLog.log("makeVibRespImg intImg.shape= %s intImg max= %f min= %f" % (repr(intImg.shape), np.max(intImg), np.min(intImg)))
    intThreshIdx = np.where(intImg < intThresh)
    
    # make the intensity image
    intImg = VolumeScan.makeImgSliceFromVolume(volData, widthStep, correctAspectRatio=False)
    intImg = intImg // 2  # since intensity image maps from 0.255

    if(vibDataType == VibDataType.STIM):
        magData = mscanRgnData.magResp[:, widthStep, :, freqIdx, ampIdx]      
        phaseData = mscanRgnData.phaseResp[:, widthStep, :, freqIdx, ampIdx]      
    else:
        magData = mscanRgnData.magResp[:, widthStep, :, freqIdx, ampIdx]      
        phaseData = mscanRgnData.phaseResp[:, widthStep, :, freqIdx, ampIdx]      
    # TODO add more data types when they become avaiable
    
    
    # since we are going to be modifying magData, we need to make acopy
    magData = magData.copy()  
    phaseData = phaseData.copy()
    
    if(vibNoiseRgn == VibNoiseRegion.ABOVE):
        noiseMean = mscanRgnData.FDnoiseAboveStimFreqMean[:, widthStep, :, freqIdx, ampIdx]
        noiseSD = mscanRgnData.FDnoiseAboveStimFreqSD[:, widthStep, :, freqIdx, ampIdx]
    else:
        noiseMean = mscanRgnData.FDnoiseBelowStimFreqMean[:, widthStep, :, freqIdx, ampIdx]
        noiseSD = mscanRgnData.FDnoiseBelowStimFreqSD[:, widthStep, :, freqIdx, ampIdx]
        
    magNoiseThresh = noiseMean + noiseNumSD*noiseSD
    vibThreshIdx = np.where(magData < magNoiseThresh)
    
    # set points below threshold to 0 
    magData[intThreshIdx] = 0
    magData[vibThreshIdx] = 0
    
    # normalize the maggnitude and phase data
    if(magAutoNorm):
        nL = np.min(magData)
        nH = np.max(magData)
    else:
        nL = norms(0)
        nH = norms(1)

    minMag = nL
    maxMag = nH

    magImg = 127*(magData - nL)/(nH-nL)
    magImg = 128 + np.clip(magImg, 0, 127)
    
    nL = -np.pi
    nH = np.pi
    phaseImg = 127*(phaseData - nL)/(nH-nL)
    phaseImg = 128 + np.clip(phaseImg, 0, 127)

    # replace secintos which where below threshold with image data
    #magImg = magImg.transpose()
    #phaseImg = phaseImg.transpose()
    if DebugLog.isLogging:
        DebugLog.log("bRespImg magImg.shape= %s minMag= %f maxMag= %f" % (repr(magImg.shape), minMag, maxMag))
    
    magImg[intThreshIdx] = intImg[intThreshIdx]
    magImg[vibThreshIdx] = intImg[vibThreshIdx]
    phaseImg[intThreshIdx] = intImg[intThreshIdx]
    phaseImg[vibThreshIdx] = intImg[vibThreshIdx]
    
    xRes = mscanRgnData.xRes
    zRes = mscanRgnData.zRes
    
    magImg = BScan.correctImageAspectRatio(magImg, xRes, zRes, PIL.Image.NEAREST)
    phaseImg = BScan.correctImageAspectRatio(phaseImg, xRes, zRes, PIL.Image.NEAREST)

    return (magImg, phaseImg, minMag, maxMag)          
                    
def displayMscanDataSinglePt(appObj, mscanData, tuningCurve):
    if mscanData is not None:
        pl = appObj.plot_mscan_avgAline
        pl.clear()
        pl.plot(mscanData.avgAline, pen="b")
        al_max = np.max(mscanData.avgAline)
        al_min = np.min(mscanData.avgAline)
        if hasattr(mscanData, 'zROIidx'):
            idx1 = mscanData.zROIidx[0]
            idx2 = mscanData.zROIidx[1]
            pl.plot([idx1, idx1], [al_min, al_max], pen="m")
            pl.plot([idx2, idx2], [al_min, al_max], pen="m")
            # zROIpl1.sigClicked.connect(self.zROIplot1Clicked())
        if hasattr(mscanData, 'maxSNR_Zidx'):
            idx = mscanData.maxSNR_Zidx
            DebugLog.log("displayMscanDataSinglePt: mscanData.maxSNR_Zidx = %d" % (idx))
            pl.plot([idx, idx], [al_min, al_max], pen="r")
            
        if len(mscanData.stimRespMag.shape) == 1:
            appObj.mscan_stim_mag_resp.setText("%0.3g nm" % (mscanData.stimRespMag))
    
        if len(mscanData.mscanTD.shape) == 1:
            pl = appObj.plot_mscan_TD
            pl.clear()
            t = np.linspace(0, mscanData.trialTime, len(mscanData.mscanTD))
            
            pl.plot(t, 1e-9*mscanData.mscanTD, pen="b")
            labelStyle = appObj.xLblStyle
            pl.setLabel('bottom', 'Time', 's', **labelStyle)
            labelStyle = appObj.yLblStyle
            pl.setLabel('left', 'Displacement', 'm', **labelStyle)
        
        DebugLog.log("displayMscanDataSinglePt: mscanData.maxFFTFreq = %g" % (mscanData.maxFFTFreq))
        if len(mscanData.mscanFFT.shape) == 1:
            freq = np.linspace(0, mscanData.maxFFTFreq, len(mscanData.mscanFFT))
            
            # don't plot below 80 Hz because of large LF component that throws off scaling
            idx = int(np.floor(80 / freq[1]))
            pl = appObj.plot_mscan_FFT
            pl.clear()
            pl.plot(freq[idx:], 1e-9*np.abs(mscanData.mscanFFT[idx:]), pen="b")
            labelStyle = appObj.xLblStyle
            pl.setLabel('bottom', 'Frequency', 'Hz', **labelStyle)
            labelStyle = appObj.yLblStyle
            pl.setLabel('left', 'Displacement', 'm', **labelStyle)
            
            noiseNumSD = 3
            noiseAbove = mscanData.FDnoiseAboveStimFreqMean + noiseNumSD * mscanData.FDnoiseAboveStimFreqSD
            noiseBelow = mscanData.FDnoiseBelowStimFreqMean + noiseNumSD * mscanData.FDnoiseBelowStimFreqSD
            
            appObj.mscan_noise_below_lineEdit.setText("%0.3g nm" % (noiseAbove))
            appObj.mscan_noise_above_lineEdit.setText("%0.3g nm" % (noiseBelow))
            
            clr = QtGui.QColor(128, 0, 0)     # dark red
            qpen = QtGui.QPen(clr)
            qbrush = QtGui.QBrush(clr)
        
            stimFreq = [ mscanData.stimFreq ]
            stimResp = [ mscanData.stimRespMag*1e-9 ]
            pl.plot(stimFreq, stimResp, symbol='o', pen=qpen, symbolBrush=qbrush)
    
    if tuningCurve is not None:
        DebugLog.log("displayMscanDataSinglePt() plotting mscan tuning curves")
        mag_plt = appObj.plot_mscan_mag_tuning
        phase_plt = appObj.plot_mscan_phase_tuning
        mag_plt.clear()
        phase_plt.clear()
        
        #yAxis = mag_plt.getAxis('left')
        #yAxis.setLogMode(True)
        mag_plt.setLogMode(y=True)
    
        # penColors = ['b', 'r', 'g', 'y']
        tcurve = tuningCurve
        numAmp = len(tcurve.amp)
        for aIdx in range(0, numAmp):
            # mag = tcurve.magResp[:, aIdx] * 1e-9
            mag = tcurve.magResp[:, aIdx] 
            pn = appObj.penArray[aIdx % len(appObj.penArray)]
            br = appObj.brushArray[aIdx % len(appObj.brushArray)]
            mag_plt.plot(tcurve.freq*1e3, mag, symbol='o', pen=pn, symbolBrush=br)
            
            ph = tcurve.phaseRespUnwrapped[:, aIdx]
            phase_plt.plot(tcurve.freq*1e3, ph, symbol='o', pen=pn, symbolBrush=br)
    
        labelStyle = appObj.xLblStyle
        mag_plt.setLabel('bottom', 'Frequency', 'Hz', **labelStyle)
        phase_plt.setLabel('bottom', 'Frequency', 'Hz', **labelStyle)
        labelStyle = appObj.yLblStyle
        mag_plt.setLabel('left', 'Displacement (nm)', '', **labelStyle)
        phase_plt.setLabel('left', 'Phase', 'rad', **labelStyle)


    
def displayMscanRegionData(mscanRegionData, volumeData, appObj, useLastFreqAmpIdx=True):
    intThreshVal = 65535*appObj.mscan_intensityThresh_slider.value()/100
    numSD = appObj.mscan_noise_numSD_dblSpinBox.value()
    
    if not useLastFreqAmpIdx:
        # TODO set frequency, amp and width sstep indices based off UI elements (probably fromcombo boxes )
        mscanRegionData.freqIdx = appObj.mscan_vol_freq_comboBox.currentIndex()
        mscanRegionData.ampIdx = appObj.mscan_vol_amp_comboBox.currentIndex()
    
    DebugLog.log("MScan.displayMscanRegionData() intThreshVal= %d numSD= %f" % (intThreshVal, numSD))
    # def makeVibRespImg(mscanRgnData, volData, magAutoNorm=True, magNorms=(0, 100), intThresh=0, vibDataType=VibDataType.STIM, vibNoiseRgn=VibNoiseRegion.ABOVE, noiseNumSD=3):
    (magImg, phaseImg, minMag, maxMag) = makeVibRespImg(mscanRegionData, volumeData, intThresh=intThreshVal, noiseNumSD=numSD)
    #(magImg, phaseImg) = makeVibRespImg(aggData.mscanRegionData, self.volumeData)
    
    magClrMap = ROIImageGraphicsView.COLORMAP_JET_BW
    phaseClrMap = ROIImageGraphicsView.COLORMAP_HSV_BW
    # rset = (appObj.frameNumLast == 0)  
    rset = True
    appObj.mscan_img_mag_roi_gv.setImage(magImg, magClrMap, resetTransform=rset)
    appObj.mscan_img_phase_roi_gv.setImage(phaseImg, phaseClrMap, resetTransform=rset)
    appObj.mscan_img_mag_min_lbl.setText("%0.1f" % (minMag))
    appObj.mscan_img_mag_max_lbl.setText("%0.1f" % (maxMag))
    appObj.mscan_img_mag_qtr_max_lbl.setText("%0.1f" % (maxMag/4))
    appObj.mscan_img_mag_half_max_lbl.setText("%0.1f" % (maxMag/2))
    appObj.mscan_img_mag_3qtr_max_lbl.setText("%0.1f" % (3*maxMag/4))
    magClrMap = ROIImageGraphicsView.COLORMAP_JET
    phaseClrMap = ROIImageGraphicsView.COLORMAP_HSV
    sz = appObj.mscan_img_mag_clrbar_lbl.size()
    w = sz.width()
    h = sz.height()
    imgData = np.round(np.linspace(255, 0, h))
    imgData = np.tile(imgData, (w, 1))
    imgData = imgData.transpose()
    DebugLog.log("MScan.displayMscanRegionData() imgData.shape= %s minMag= %f maxMag= %f" % (repr(imgData.shape), minMag, maxMag))
    imgData = np.require(imgData, np.uint8, 'C')
    qImg = QtGui.QImage(imgData.data, w, h, w, QtGui.QImage.Format_Indexed8)
    qImg.setColorTable(magClrMap)
    qPixMap = QtGui.QPixmap.fromImage(qImg)
    appObj.mscan_img_mag_clrbar_lbl.setPixmap(qPixMap)
    
    sz = appObj.mscan_img_phase_clrbar_lbl.size()
    w = sz.width()
    h = sz.height()
    imgData = np.round(np.linspace(255, 0, h))
    imgData = np.tile(imgData, (w, 1))
    imgData = imgData.transpose()
    imgData = np.require(imgData, np.uint8, 'C')
    qImg = QtGui.QImage(imgData.data, w, h, w, QtGui.QImage.Format_Indexed8)
    qImg.setColorTable(phaseClrMap)
    qPixMap = QtGui.QPixmap.fromImage(qImg)
    appObj.mscan_img_phase_clrbar_lbl.setPixmap(qPixMap)
    
    rgnData = mscanRegionData
    widthStep = rgnData.widthStep
    imgData = VolumeScan.makeImgSliceFromVolume(volumeData, widthStep)
    DebugLog.log("MScan.displayMscanRegionData()() imgData.shape= %s" % (repr(imgData.shape)))
    
    rset = (rgnData.widthStep, rgnData.lengthStep, rgnData.freqIdx, rgnData.ampIdx) == (0, 0, 0, 0)
    appObj.mscan_img_vol_roi_gv.setImage(imgData, ROIImageGraphicsView.COLORMAP_HOT, rset)
    if rset:
        gvs = (appObj.mscan_img_vol_roi_gv, appObj.mscan_img_mag_roi_gv, appObj.mscan_img_phase_roi_gv)
        for gv in gvs:
            hscroll = gv.horizontalScrollBar()
            hscroll.setSliderPosition(-500)
            vscroll = gv.verticalScrollBar()
            vscroll.setSliderPosition(-500)
                    
def GetTestData(frameNum):
    return oct_data

# this class is used to define raw data for a given frame
class MscanRawData:
    def __init__(self):
        self.frameNum = None
        self.oct_data = None
        self.mic_data = None
        self.mscanPosAndStim = None
        self.lastFrame = False
        
def MscanGetStepFromFrameNum(frameNum, scanParams, audioParams):
    numAmpSteps = len(audioParams.amp)
    numFreqSteps = audioParams.getNumFrequencies()
    numLenSteps = scanParams.lengthSteps
    numWidthSteps = scanParams.widthSteps
    numSoundSteps = numAmpSteps*numFreqSteps

    posWidthStep = frameNum // (numSoundSteps * numLenSteps)
    posLenStep = (frameNum // numSoundSteps) % numLenSteps
    freqStep = (frameNum // numAmpSteps) % numFreqSteps 
    ampStep = frameNum % numAmpSteps
    
    if posWidthStep >= numWidthSteps:
        return -1, -1, -1, -1
        
    return posLenStep, posWidthStep, freqStep, ampStep
    
# collect mscan data for a given frame number
# oct_hw is a LV_DLL_Interface
def MscanCollectFcn(oct_hw, frameNum, trigRate, extraArgs):
    t1 = time.time()
    scanParams = extraArgs[0]
    audioParams = extraArgs[1]
    mirrorDriver = extraArgs[2]
    audioHW = extraArgs[3]
    zROI = extraArgs[4]
    testDataDir =  extraArgs[5]
    
    numAmpSteps = len(audioParams.amp)
    numFreqSteps = audioParams.getNumFrequencies()
    numLenSteps = scanParams.lengthSteps
    numWidthSteps = scanParams.widthSteps
    posLenStep, posWidthStep, freqStep, ampStep = MscanGetStepFromFrameNum(frameNum, scanParams, audioParams)
    # if beyond last step, then return nothing
    if posLenStep < 0:
        return None, None
    
    mirrChanNames = [mirrorDriver.X_daqChan, mirrorDriver.Y_daqChan]
    mirrOutData = np.zeros(2)
    
    outputRate = audioHW.DAQOutputRate
    inputRate = audioHW.DAQInputRate
    trigChan = mirrorDriver.trig_daqChan
    numSpk = audioParams.getNumSpeakers()
    if not oct_hw.IsOCTTestingMode():
        from DAQHardware import DAQHardware
        daq = DAQHardware()
    chanNamesIn= [ audioHW.mic_daqChan]
    
    chanNamesOut = [audioHW.speakerL_daqChan]
    attenLines = audioHW.attenL_daqChan
    spkNum = 0
    if audioParams.speakerSel == Speaker.RIGHT:
        chanNamesOut = [audioHW.speakerR_daqChan]
        spkNum = 1
        attenLines = audioHW.attenR_daqChan
    elif audioParams.speakerSel == Speaker.BOTH:
        chanNamesOut = [audioHW.speakerl_daqChan, audioHW.speakerR_daqChan]
    
    mscanPosAndStim = MscanPosAndStim()
    
    # set mirror position
    (xPos, yPos) = getXYPos(posLenStep, posWidthStep, scanParams)
    (x_cmd, y_cmd) = mirrorDriver.makeMirrorCommand(xPos, yPos)
    mirrOutData[0] = x_cmd
    mirrOutData[1] = y_cmd
    if not oct_hw.IsDAQTestingMode():
        daq.writeValues(mirrChanNames, mirrOutData)
    
    # set up audio output
    freq = audioParams.freq[spkNum, freqStep]
    amp = audioParams.amp[ampStep]
    endIdx = int(5e-3 * outputRate)        # only plot first 5 ms
    
    if audioParams.speakerSel == Speaker.BOTH:
        audioOutputL, attenLvlL = makeAudioOutput(audioParams, audioHW, 0, freq, amp)
        audioOutputR, attenLvlR = makeAudioOutput(audioParams, audioHW, 1, freq, amp)
        numOutputSamples = len(audioOutputL)
        # t = np.linspace(0, numOutputSamples/outputRate, numOutputSamples)
        # pl.plot(t[0:endIdx], audioOutputL[0:endIdx], pen='b')
        # pl.plot(t[0:endIdx], audioOutputR[0:endIdx], pen='r')
        
        audioOutput = np.vstack((audioOutputL, audioOutputR))
        DebugLog.log("Mscan CollectionProcess(): attenLvlL= %s attenLvlR= %s" % (repr(attenLvlL), repr(attenLvlR)))
        # set attenuator level
        if not oct_hw.IsDAQTestingMode():
            attenSig = AudioHardware.makeLM1972AttenSig(attenLvlL)
            appObj.oct_hw.SetAttenLevel(attenLvlL, attenLines)
            #daq.sendDigOutCmd(audioHW.attenL_daqChan, attenSig)
            attenSig = AudioHardware.makeLM1972AttenSig(attenLvlR)
            appObj.oct_hw.SetAttenLevel(attenLvlR, attenLines)
            #daq.sendDigOutCmd(audioHW.attenL_daqChan, attenSig)
    else:
        audioOutput, attenLvl = makeAudioOutput(audioParams, audioHW, spkNum, freq, amp)
        DebugLog.log("Mscan CollectionProcess(): attenLvL= %s " %  repr(attenLvl))

        numOutputSamples = len(audioOutput)
#        t = np.linspace(0, numOutputSamples/outputRate, numOutputSamples)
#        pl.plot(t[0:endIdx], audioOutput[0:endIdx], pen='b')
        if not oct_hw.IsDAQTestingMode():
            attenSig = AudioHardware.makeLM1972AttenSig(attenLvl)
            # daq.sendDigOutCmd(attenLines, attenSig)
            oct_hw.SetAttenLevel(attenLvl, attenLines)
    
    numInputSamples = int(inputRate*numOutputSamples/outputRate) 
    if not oct_hw.IsDAQTestingMode():
        daq.setupAnalogOutput(chanNamesOut, trigChan, outputRate, audioOutput.transpose())
    
        # setup the input task
        daq.setupAnalogInput(chanNamesIn, trigChan, int(inputRate), numInputSamples) 

    numTrials = audioParams.numTrials
    oct_data = None
    for n in range(0, numTrials):
        if not oct_hw.IsDAQTestingMode():
            daq.startAnalogOutput()
            daq.startAnalogInput()
        
        # setup and grab the OCT data
        startTrigOffset = 0
        numTrigs = int(np.floor(trigRate*numOutputSamples/outputRate))
        if oct_hw.IsOCTTestingMode():
            oct_data_tmp = OCTCommon.loadRawData(testDataDir, frameNum, dataType=0)
            oct_data_tmp = oct_data_tmp[:, :, 0]
        else:
            err, oct_data_tmp = oct_hw.AcquireOCTDataFFT(numTrigs, zROI, startTrigOffset)
        
        if oct_data is None:
            shp = oct_data_tmp.shape
            oct_data = np.zeros((shp[0], shp[1], numTrials), np.complex)
            
        oct_data[:, :, n] = oct_data_tmp
        
        if not oct_hw.IsDAQTestingMode():
            mic_data = daq.readAnalogInput()
            mic_data = mic_data/audioHW.micVoltsPerPascal
            
            daq.waitDoneOutput()
            daq.stopAnalogOutput()
            daq.waitDoneInput()
            daq.stopAnalogInput()
        else:
            mic_data = OCTCommon.loadRawData(testDataDir, frameNum, dataType=3)

    # clear DAQ OUTPUT
    if not oct_hw.IsDAQTestingMode():
        daq.clearAnalogOutput()
        daq.clearAnalogInput()
        
    # copy values to object
    mscanPosAndStim.ampIdx = ampStep
    mscanPosAndStim.freqIdx = freqStep
    mscanPosAndStim.posLenStep = posLenStep
    mscanPosAndStim.posWidthStep = posWidthStep
    mscanPosAndStim.stimFreq = freq*1000
    mscanPosAndStim.numAmp = numAmpSteps
    mscanPosAndStim.numFreq = numFreqSteps    
    rawData = MscanRawData()
    
    rawData.frameNum = frameNum
    rawData.oct_data = oct_data
    rawData.mic_data = mic_data
    rawData.mscanPosAndStim = mscanPosAndStim
    rawData.collectTime = time.time() - t1

    return rawData, audioOutput
        
class MScanRegionVolData():
    def __init__(self, mscanRegionData, volData):
        self.mscanRegionData = mscanRegionData
        self.volData = volData
        
def MscanProcessingProcess(audioParams, scanParams, zROI, regionMscan, procOpts, trigRate, framesPerScan, rawDataQ, procDataQ, procRawDataQ, msgQ, statusQ):
    shutdown = False
    mscanTuningCurveList = None
    numZPts = zROI[1] - zROI[0] + 1
    mscanRegionData = None
    volData = None
    mscanData = None
    
    frameNum = 0
    putTimeout = False  # indicates there was a timeout attempting to send data to 
    while not shutdown and frameNum < framesPerScan:
        try:
            if not putTimeout:
                if not rawDataQ.empty():
                    mscanData = None
                    rawData = rawDataQ.get(timeout=0.25)
                    
                    if rawData is not None and isinstance(rawData, MscanRawData) and not putTimeout:
                        # convet to cocorret type
                        #if isinstance(data, MScanData):
                        #    mscanData = data
                        #rawData.frameNum
                        #rawData.mic_data
                        DebugLog.log("MscanProcessingProcess(): got raw data")
                        # process the data
                        t1 = time.time()
                        frameNum = rawData.frameNum
                        posLenStep, posWidthStep, freqStep, ampStep = MscanGetStepFromFrameNum(frameNum, scanParams, audioParams)
                        mscanData, mscanTuningCurveList, mscanRegionData, volData = processMscanData(rawData.oct_data, rawData.mscanPosAndStim, scanParams, audioParams, procOpts, trigRate, mscanTuningCurveList, mscanRegionData, volData)
                        mscanRgnVolData = MScanRegionVolData(mscanRegionData, volData)
                        mscanData.frameNum = frameNum
                        mscanData.collectTime = rawData.collectTime
                        mscanData.processTime = time.time() - t1
                        
                        
            # send processsed data to main program
            if mscanData is not None:
                try:
                    if procRawDataQ is not None:
                        procRawDataQ.put(rawData, timeout=0.25)
                    if regionMscan:
                        procDataQ.put(mscanRgnVolData, timeout=0.25)
                    else:
                        ptNum = posLenStep
                        tCurve = mscanTuningCurveList[ptNum]
                        procDataQ.put(tCurve, timeout=0.25)
                        #(mscanPosAndStim.ampIdx == (numAmpSteps - 1)) and (mscanPosAndStim.freqIdx == (numFreqSteps - 1))
                        
                    procDataQ.put(mscanData, timeout=0.25)   # send this data last because it contains frame number which client uss to detect whether acquisition is complete                                
                    putTimeout = False
                except queue.Full:
                    putTimeout = True
        except Exception as ex:
            traceback.print_exc(file=sys.stdout)
            statusMsg = OCTCommon.StatusMsg(OCTCommon.StatusMsgSource.PROCESSING, OCTCommon.StatusMsgType.ERROR)
            statusMsg.param = ex
            try:
                statusQ.put(statusMsg, False)
            except queue.Full as ex:
                pass
            shutdown = True
        
        # chedk for shutdown messages 
        if not msgQ.empty():
            msg = msgQ.get()
            if msg == 'shutdown':
                shutdown = True

def handleStatusMessage(statusMsg):
    err = False
    # if statusMsg.msgSrc == OCTCommon.StatusMsgSource.COLLECTION:
    if statusMsg.msgType == OCTCommon.StatusMsgType.ERROR:
        err = True
    elif statusMsg.msgType == OCTCommon.StatusMsgType.DAQ_OUTPUT:
        pass
        
    return err
    
def runMscanMultiProcess(appObj, scanParams, zROI, procOpts, trigRate, testDataDir, regionMscan):
    mirrorDriver = appObj.mirrorDriver
    saveOpts = appObj.getSaveOpts()
    audioParams = appObj.getAudioParams()

    rset = True
    isSaveDirInit = False
    audioHW = appObj.audioHW
    oct_hw = appObj.oct_hw
    
    #procDataQ = mproc.Queue(4)        
#    procMsgQ = mproc.Queue(4)        
 #   dataQ = oct_hw.
    # start up the processing process
  #  procProc = mproc.Process(target=MscanProcessingProcess, args=(audioParams, scanParams, zROI, regionMscan, procOpts, dataQ, procDataQ, procMsgQ), daemon=True)
   # procProc.start()
    startTime = time.time()    
    DebugLog.log("runBScanMultiProcess: new acquisiiton")
    oct_hw.NewAcquisition()
    oct_hw.SetSendExtraInfo(False)   # do not send audio output
    DebugLog.log("runBScanMultiProcess: setting acquire function")
    oct_hw.SetAcqFunction(MscanCollectFcn)
    extraArgs = [scanParams, audioParams, mirrorDriver, audioHW, zROI, testDataDir]
    DebugLog.log("runBScanMultiProcess: setting acquire functiona args")
    oct_hw.SetAcqFunctionArgs(extraArgs)
    DebugLog.log("runBScanMultiProcess: SetAcqFunctionArgs() exit")
    
    mscanTuningCurveList = None
    mscanRegionData = None
    volData = None
    startAcq = True
    
    numAmpSteps = len(audioParams.amp)
    numFreqSteps = audioParams.getNumFrequencies()
    numLenSteps = scanParams.lengthSteps
    numWidthSteps = scanParams.widthSteps
    
    framesPerScan = numAmpSteps*numFreqSteps*numLenSteps*numWidthSteps
        
    DebugLog.log("runBScanMultiProcess: cleaning status message log")
    statusMsg = oct_hw.GetStatus()
    while statusMsg is not None:
        DebugLog.log("runBScanMultiProcess: got status message type=" + repr(statusMsg.msgType))
        err = handleStatusMessage(statusMsg)
        statusMsg = oct_hw.GetStatus()        
    
    procDataQ = mproc.Queue(10)
    procRawDataQ = None
    if saveOpts.saveRaw:   # if saving raw data, create a raw data queue so that mscan processing process will resned raw data to this function
        procRawDataQ = mproc.Queue(10)
    msgQ = mproc.Queue(10)
    rawDataQ = oct_hw.rawDataQ
    statusQ = oct_hw.statusQ
    procProcess = mproc.Process(target=MscanProcessingProcess, args=[audioParams, scanParams, zROI, regionMscan, procOpts, trigRate, framesPerScan, rawDataQ, procDataQ, procRawDataQ, msgQ, statusQ], daemon=True)
    DebugLog.log("runBScanMultiProcess(): starting processing process")
    procProcess.start()
    frameNum = -1
        
    while not appObj.doneFlag and frameNum < framesPerScan-1:
        # update parameters in background process
        # start the acquisitio on first loop iteration
        # we don't just do this outside the loop because we haven't loaded function args yet
        if startAcq:  
            DebugLog.log("runBScanMultiProcess: starting acquisition")
            oct_hw.StartAcquisition() 
            startAcq = False
        
        
        if not procDataQ.empty():
            DebugLog.log("runBScanMultiProcess: grabbing data")

            data = procDataQ.get()      
            mscanRegionData = None
            mscanData = None
            tuningCurve = None
            if isinstance(data, MScanData):
                DebugLog.log("runBScanMultiProcess: received mscan data")
                frameNum = data.frameNum
                mscanData = data
                appObj.acquisition_progressBar.setValue(round(100*(frameNum+1)/framesPerScan))                
                appObj.mscanCollectionTime_label.setText("%0.1f ms" % (mscanData.collectTime * 1000))
                appObj.mscanProcessTime_label.setText("%0.1f ms" % (mscanData.processTime * 1000))

                if not regionMscan:
                    displayMscanDataSinglePt(appObj, mscanData, tuningCurve)
                    
            elif isinstance(data, MscanTuningCurve):
                DebugLog.log("runBScanMultiProcess: received tuning curve data")
                tuningCurve = data
                displayMscanDataSinglePt(appObj, mscanData, tuningCurve)
            elif isinstance(data, MScanRegionVolData):
                DebugLog.log("runBScanMultiProcess: received mscan region data")
                mscanRegionData = data.mscanRegionData
                volData = data.volData
                displayMscanRegionData(mscanRegionData, volData, appObj, useLastFreqAmpIdx=True)          

            # save the mscan tuning curve
            if appObj.getSaveState():
                if not isSaveDirInit:
                    saveDir = OCTCommon.initSaveDir(saveOpts, 'MScan', scanParams=scanParams, audioParams=audioParams)
                    isSaveDirInit = True
                if regionMscan:
                    if mscanRegionData is not None:
                        saveMscanRegionData(mscanRegionData, volData, saveDir)
                else:
                    if mscanData is not None:
                        saveMscanDataTimeDomain(mscanData, freq, amp, frameNum, saveDir)                   
                    if tuningCurve is not None:
                        saveMscanTuningCurve(tuningCurve, audioParams, posLenStep, saveDir)
                        
                if saveOpts.saveRaw:
                    if not procRawDataQ.empty():
                        rawData = procRawDataQ.get()
                        OCTCommon.saveRawData(rawData.oct_data, saveDir, frameNum, dataType=0)
                        OCTCommon.saveRawData(rawData.mic_data, saveDir, frameNum, dataType=3)
                    
            statusMsg = oct_hw.GetStatus()
            while statusMsg is not None:
                DebugLog.log("runBScanMultiProcess: got status message type=" + repr(statusMsg.msgType))
                err = handleStatusMessage(statusMsg)
                if err:
                    appObj.doneFlag = True  # if error occured, stop pcollecting
                statusMsg = oct_hw.GetStatus()
            

        tElapsed = time.time() - startTime
        tMins = int(np.floor(tElapsed / 60))
        tSecs = int(tElapsed - 60*tMins)
        appObj.timeElapsed_label.setText("%d mins %d secs" % (tMins, tSecs))

        # check for GUI events, particularly the "done" flag
        QtGui.QApplication.processEvents() 
        time.sleep(0.005)
    
    msgQ.put('shutdown')  # tell processing process to stop
    DebugLog.log("runBScanMultiProcess: finishd acquiring data")        
    oct_hw.PauseAcquisition()        
    appObj.isCollecting = False
    QtGui.QApplication.processEvents() # check for GUI events
    appObj.finishCollection()    
    

def runMScan(appObj, multiProcess=False):
    DebugLog.log("runMscan")
    oct_hw = appObj.oct_hw
    
    try: 
        appObj.doneFlag = False
        appObj.isCollecting = True    
        startTime = time.time()
        
        if not appObj.oct_hw.IsOCTTestingMode():
            from DAQHardware import DAQHardware
            daq = DAQHardware()
            
        scanParams, roiBegin, roiEnd, zROIIndices, zROIspread = makeMscanScanParamsAndZROI(appObj)
        if scanParams is None:
            QtGui.QMessageBox.critical (appObj, "Error", "No point or region has been set for the Mscan.  Please select a point or region.")
            appObj.tabWidget.setCurrentIndex(0)   # B-Mscan, Volume Mscan
            appObj.Mscan_pushButton.setChecked(False)
            appObj.isCollecting = False
            return

        trigRate = appObj.octSetupInfo.getTriggerRate()
        mirrorDriver = appObj.mirrorDriver
        saveOpts = appObj.getSaveOpts()
        isSaveDirInit = False
        audioParams = appObj.getAudioParams()
        locPixMap = QtGui.QPixmap.grabWidget(appObj.bscan_img_gv)

        if (scanParams.lengthSteps == 1 and scanParams.widthSteps == 1):
            regionMscan = False
            appObj.tabWidget.setCurrentIndex(3)
        else:   # region mscan
            regionMscan = True
            appObj.tabWidget.setCurrentIndex(4)
        
        DebugLog.log("Mscan runMscan(): regionMscan= %s" % repr(regionMscan))
        
        # if in testing mode, load proper paramaeters instead of getting them from GUI
        if oct_hw.IsOCTTestingMode():
          
            processMode = OCTCommon.ProcessMode(appObj.processMode_comboBox.currentIndex())
            if processMode == OCTCommon.ProcessMode.SOFTWARE:
                appObj.savedDataBuffer.loadData(appObj)
            else:
                if regionMscan:
                    testDataDir = os.path.join(appObj.basePath, 'exampledata', 'MScan B-Mscan')
                else:
                    testDataDir = os.path.join(appObj.basePath, 'exampledata', 'MScan single pt')
                filePath = os.path.join(testDataDir, 'ScanParams.pickle')
                f = open(filePath, 'rb')
                scanParams = pickle.load(f)
                f.close()
                filePath = os.path.join(testDataDir, 'AudioParams.pickle')
                zROIIndices = [85]
                f = open(filePath, 'rb')
                audioParams = pickle.load(f)
                f.close()
                trigRate = 49.9598e3
                        
        procOpts = MscanProcOpts()
        procOpts.bscanNormLow = appObj.normLow_spinBox.value()
        procOpts.bscanNormHigh = appObj.normHigh_spinBox.value()
        procOpts.zRes = appObj.octSetupInfo.zRes
        procOpts.singlePt_zROI_indices = zROIIndices
        procOpts.singlePt_zROI_spread = zROIspread
        zROI = [roiBegin, roiEnd]
        procOpts.refractiveIndex= appObj.octSetupInfo.refractiveIndex
        procOpts.centerWavelength= appObj.octSetupInfo.centerWavelength
        
        if multiProcess:
            runMscanMultiProcess(appObj, scanParams, zROI, procOpts, trigRate, testDataDir, regionMscan)
            return
        
        rset = True
        
        frameNum = 0
        posLenStep = 0
        posWidthStep = 0
        freqStep = 0
        ampStep = 0
        
        numAmpSteps = len(audioParams.amp)
        numFreqSteps = audioParams.getNumFrequencies()
        numLenSteps = scanParams.lengthSteps
        numWidthSteps = scanParams.widthSteps
        frameNum = 0
        
        mirrChanNames = [mirrorDriver.X_daqChan, mirrorDriver.Y_daqChan]
        mirrOutData = np.zeros(2)
        
        
        audioHW = appObj.audioHW
        outputRate = audioHW.DAQOutputRate
        inputRate = audioHW.DAQInputRate
        trigChan = mirrorDriver.trig_daqChan
        numSpk = audioParams.getNumSpeakers()
        chanNamesIn= [ audioHW.mic_daqChan]
        
        chanNamesOut = [audioHW.speakerL_daqChan]
        attenLines = audioHW.attenL_daqChan
        spkNum = 0
        if audioParams.speakerSel == Speaker.RIGHT:
            chanNamesOut = [audioHW.speakerR_daqChan]
            spkNum = 1
            attenLines = audioHW.attenR_daqChan
        elif audioParams.speakerSel == Speaker.BOTH:
            chanNamesOut = [audioHW.speakerL_daqChan, audioHW.speakerR_daqChan]
        
        mscanPosAndStim = MscanPosAndStim()
        
        mscanTuningCurveList = None
        mscanRegionData = None
        volData = None
        
        while not appObj.doneFlag and posWidthStep < numWidthSteps:
            # set mirror position
            (xPos, yPos) = getXYPos(posLenStep, posWidthStep, scanParams)
            (x_cmd, y_cmd) = mirrorDriver.makeMirrorCommand(xPos, yPos)
            mirrOutData[0] = x_cmd
            mirrOutData[1] = y_cmd
            if not oct_hw.IsDAQTestingMode():
                daq.writeValues(mirrChanNames, mirrOutData)
        
            processMode = OCTCommon.ProcessMode(appObj.processMode_comboBox.currentIndex())
        
            # set up audio output
            freq = audioParams.freq[spkNum, freqStep]
            amp = audioParams.amp[ampStep]
            pl = appObj.plot_spkOut
            pl.clear()
            endIdx = int(5e-3 * outputRate)        # only plot first 5 ms
            
            if audioParams.speakerSel == Speaker.BOTH:
                audioOutputL, attenLvlL = makeAudioOutput(audioParams, audioHW, 0, freq, amp)
                audioOutputR, attenLvlR = makeAudioOutput(audioParams, audioHW, 1, freq, amp)
                DebugLog.log("Mscan runMscan(): len(audioOutputL)= %d len(audioOutputR) = %d" % (len(audioOutputL), len(audioOutputR)))
                numOutputSamples = len(audioOutputL)
                t = np.linspace(0, numOutputSamples/outputRate, numOutputSamples)
                pl.plot(t[0:endIdx], audioOutputL[0:endIdx], pen='b')
                pl.plot(t[0:endIdx], audioOutputR[0:endIdx], pen='r')
                audioOutput = np.vstack((audioOutputL, audioOutputR))
                DebugLog.log("Mscan runMscan(): attenLvlL= %s attenLvlR= %s" % (repr(attenLvlL), repr(attenLvlR)))
                # set attenuator level
                if not oct_hw.IsDAQTestingMode():
                    attenSig = AudioHardware.makeLM1972AttenSig(attenLvlL)
                    #daq.sendDigOutCmd(audioHW.attenL_daqChan, attenSig)
                    appObj.oct_hw.SetAttenLevel(attenLvlL, attenLines)
                    attenSig = AudioHardware.makeLM1972AttenSig(attenLvlR)
                    #daq.sendDigOutCmd(audioHW.attenL_daqChan, attenSig)
                    appObj.oct_hw.SetAttenLevel(attenLvlR, attenLines)
            else:
                audioOutput, attenLvl = makeAudioOutput(audioParams, audioHW, spkNum, freq, amp)
                DebugLog.log("Mscan runMscan(): attenLvL= %s " %  repr(attenLvl))

                numOutputSamples = len(audioOutput)
                t = np.linspace(0, numOutputSamples/outputRate, numOutputSamples)
                pl.plot(t[0:endIdx], audioOutput[0:endIdx], pen='b')
                if not oct_hw.IsDAQTestingMode():
                    attenSig = AudioHardware.makeLM1972AttenSig(attenLvl)
                    #daq.sendDigOutCmd(attenLines, attenSig)
                    appObj.oct_hw.SetAttenLevel(attenLvl, attenLines)
                                
            numInputSamples = int(inputRate*numOutputSamples/outputRate) 
            if not oct_hw.IsDAQTestingMode():
                daq.setupAnalogOutput(chanNamesOut, trigChan, outputRate, audioOutput.transpose())
                daq.setupAnalogInput(chanNamesIn, trigChan, int(inputRate), numInputSamples) 

            numTrials = audioParams.numTrials
            oct_data = None
            for n in range(0, numTrials):
                if not oct_hw.IsDAQTestingMode():
                    daq.startAnalogOutput()
                    daq.startAnalogInput()
                
                # setup and grab the OCT data
                startTrigOffset = 0
                numTrigs = int(np.floor(trigRate*numOutputSamples/outputRate))
                oct_data_tmp = None
                if processMode == OCTCommon.ProcessMode.FPGA:
                    if appObj.oct_hw.IsOCTTestingMode():
                        oct_data = OCTCommon.loadRawData(testDataDir, frameNum, dataType=0)
                        oct_data_tmp = oct_data[:, :, 0]
                    else:
                        err, oct_data_tmp = appObj.oct_hw.AcquireOCTDataFFT(numTrigs, zROI, startTrigOffset)
                        
                    dataIsRaw = False
                elif processMode == OCTCommon.ProcessMode.SOFTWARE:
                    if appObj.oct_hw.IsOCTTestingMode():
                        ch0_data,ch1_data=JSOraw.getSavedRawData(numTrigs,appObj.dispData.requestedSamplesPerTrig,appObj.savedDataBuffer)
                    else:
                        # def AcquireOCTDataRaw(self, numTriggers, samplesPerTrig=-1, Ch0Shift=-1, startTrigOffset=0):
                        samplesPerTrig = appObj.oct_hw.fpgaOpts.SamplesPerTrig*2
                        err, ch0_data,ch1_data = appObj.oct_hw.AcquireOCTDataRaw(numTrigs, samplesPerTrig, startTrigOffset=startTrigOffset)
                    dataIsRaw = True
                    oct_data_tmp, klin = JSOraw.softwareProcessing(ch0_data,ch1_data,zROI,appObj)
                else:
                    QtGui.QMessageBox.critical (appObj, "Error", "Unsuppoted processing mode for current hardware")
                    break
                                
                if oct_data is None:
                    shp = oct_data_tmp.shape
                    oct_data = np.zeros((shp[0], shp[1], numTrials), np.complex)
                    
                oct_data[:, :, n] = oct_data_tmp                
                if not oct_hw.IsDAQTestingMode():
                    mic_data = daq.readAnalogInput()
                    mic_data = mic_data/audioHW.micVoltsPerPascal
                    
                    daq.waitDoneOutput()
                    daq.stopAnalogOutput()
                    daq.waitDoneInput()
                    daq.stopAnalogInput()
                else:
                    if processMode == OCTCommon.ProcessMode.SOFTWARE:
                        mic_data=scipy.signal.resample(audioOutput,numInputSamples)                                                
                    else:
                        mic_data = OCTCommon.loadRawData(testDataDir, frameNum, dataType=3)
                
                # check for done flag
                QtGui.QApplication.processEvents() 
                if appObj.doneFlag:
                    break

            if appObj.doneFlag or oct_data is None:
                break
            
            mscanPosAndStim.ampIdx = ampStep
            mscanPosAndStim.freqIdx = freqStep
            mscanPosAndStim.posLenStep = posLenStep
            mscanPosAndStim.posWidthStep = posWidthStep
            mscanPosAndStim.stimFreq = freq*1000
            mscanPosAndStim.numAmp = numAmpSteps
            mscanPosAndStim.numFreq = numFreqSteps
            
            # process the data
            mscanData, mscanTuningCurveList, mscanRegionData, volData = processMscanData(oct_data, mscanPosAndStim, scanParams, audioParams, procOpts, trigRate, mscanTuningCurveList, mscanRegionData, volData)
                
            if regionMscan:
                displayMscanRegionData(mscanRegionData, volData, appObj, useLastFreqAmpIdx=True)          
            else:
                tuningCurve = mscanTuningCurveList[posLenStep]
                displayMscanDataSinglePt(appObj, mscanData, tuningCurve)
                
            npts = len(mic_data)
            t = np.linspace(0, npts/inputRate, npts)
            pl = appObj.plot_micRaw
            pl.clear()
            pl.plot(t, mic_data, pen='b')
            
            labelStyle = appObj.xLblStyle
            pl.setLabel('bottom', 'Time', 's', **labelStyle)
            labelStyle = appObj.yLblStyle
            pl.setLabel('left', 'Response', 'Pa', **labelStyle)            
            
            if not oct_hw.IsDAQTestingMode():
                daq.clearAnalogOutput()
                daq.clearAnalogInput()

            # increment the frameNum and step            
            frameNum += 1
            ampStep += 1
            saveTC = False   # whether or not to save the tuning curve 
            if ampStep == numAmpSteps:
                ampStep = 0
                freqStep += 1
                if freqStep == numFreqSteps:
                    freqStep = 0
                    saveTC = True    # save tuning curve after all freq/amp positions reached
                    posLenStep += 1
                    if posLenStep == numLenSteps:
                        posLenStep = 0
                        posWidthStep += 1

            # save data
            if appObj.getSaveState():
                if not isSaveDirInit:
                    saveDir = OCTCommon.initSaveDir(saveOpts, 'MScan', scanParams=scanParams, audioParams=audioParams)
                    fName = 'location.png'
                    fPath = os.path.join(saveDir, fName)
                    locPixMap.save(fPath, 'PNG')
                    isSaveDirInit = True
                    # TODO add save code here
                if regionMscan:
                    saveMscanRegionData(mscanRegionData, volData, saveDir)
                else:
                    saveMscanDataTimeDomain(mscanData, freq, amp, frameNum, saveDir)                   
                    if saveTC:
                        tuningCurve = mscanTuningCurveList[posLenStep]
                        saveMscanTuningCurve(tuningCurve, audioParams, posLenStep, saveDir)                        
                        
                if saveOpts.saveRaw:
                    if dataIsRaw:
                        outfile = os.path.join(saveDir, 'RawData %d.npz' % (frameNum-1))
                        np.savez_compressed(outfile, ch0_data=ch0_data, ch1_data=ch1_data)
                    else:
                        OCTCommon.saveRawData(oct_data, saveDir, frameNum-1, dataType=0)
                    
                    OCTCommon.saveRawData(mic_data, saveDir, frameNum-1, dataType=3)

            framesPerScan = numAmpSteps*numFreqSteps*numLenSteps*numWidthSteps
            appObj.acquisition_progressBar.setValue(round(100*frameNum/framesPerScan))
            tElapsed = time.time() - startTime
            tMins = int(np.floor(tElapsed / 60))
            tSecs = int(tElapsed - 60*tMins)
            appObj.timeElapsed_label.setText("%d mins %d secs" % (tMins, tSecs))
            
            # check for GUI events, particularly the "done" flag
            QtGui.QApplication.processEvents() 
    except Exception as ex:
        # raise ex
        traceback.print_exc(file=sys.stdout)
        QtGui.QMessageBox.critical (appObj, "Error", "Error during scan. Check command line output for details")
    finally:
        appObj.isCollecting = False
        QtGui.QApplication.processEvents() # check for GUI events
        appObj.finishCollection()    
