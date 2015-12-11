# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 12:22:44 2015

@author: OHNS
"""

import OCTCommon
from DebugLog import DebugLog
import scipy.signal
# import OCTFPGAProcessingInterface as octfpga

from PyQt4 import QtCore, QtGui, uic
from DebugLog import DebugLog
from ROIImageGraphicsView import *
from OCTProtocolParams import *
import os
import pickle
import traceback
import sys
import JSOraw
#class DispersionData:
#    def __init__(self):
#        self.magWin = []
#        self.phaseCorr = []
#        self.phDiode = []
#        self.uncorrAline = []
#        self.corrAline = []
#        self.phDiode_background = None

class DispersionData:  # this class holds all the dispersion compensation data
    def __init__(self, fpgaOpts=None):
        self.magWin = None
        self.phaseCorr = None
        self.phDiode_background = None
        self.requestedSamplesPerTrig = -1
        self.startSample = -1
        self.endSample = -1
        self.numKlinPts = -1
        self.numShiftPts = 0
        self.sampleOffset = 0
        self.filterWidth = 0
        self.PDfilterCutoffs=[] 
        self.mziFilter=[]
        self.magWin_LPfilterCutoff=[]
        self.k0Reference = []
        self.uncorrAline = []
        self.corrAline = []
        self.Klin=[]
        self.dispCode=-1
        self.dispMode=-1
            
        
        if fpgaOpts is not None:
            self.requestedSamplesPerTrig = fpgaOpts.SamplesPerTrig*2
            self.startSample = fpgaOpts.klinRoiBegin
            self.endSample = fpgaOpts.klinRoiEnd
            self.numKlinPts = fpgaOpts.numKlinPts
            self.numShiftPts = fpgaOpts.Ch0Shift*2
            self.sampleOffset = fpgaOpts.SampleOffset*2
        
def loadDispData(appObj, fileName='dispersionLast.pickcle'):
    infile=os.path.join(appObj.configPath, 'Dispersion', fileName)  
    file2 = open(infile,'rb')
    dispData = pickle.load(file2)
    file2.close()
    
    return dispData
    
# process dispersion data
# interpPD = 2D numpy array, trigger# is first indiex, ptNum is second
# numklinpts = num klinpts
# PDfilterCutoffs = array like , 2 elements, with low/high filter cuoffs for the photodiode dieata, should be between 0 and 0.5
# magWin_LPfilterCutoff = low pass filter cutoffs for the magnitude window, hould be between 0 and 0.5
def processData(interpPD, dispData, numklinpts, PDfilterCutoffs, magWin_LPfilterCutoff, pd_background=None, bg_collect=False):
    # scanP = self.scanParams

    numTrigs = interpPD.shape[0]
    numPts = interpPD.shape[1]
    DebugLog.log("Dispersion.processData numTrigs= %d numPts= %d" % (numTrigs, numPts))
    
    pd = np.mean(interpPD, 0)
    dispData.phDiode = pd
    
    winFcn = np.hanning(numklinpts)
    k = np.linspace(0, 1000, numklinpts)

    if bg_collect:
        pd = np.real(interpPD[:, 0:numklinpts])
        pd_background = np.mean(pd, 0)
    
    magWin = np.zeros((numTrigs, numklinpts))
    phaseCorr = np.zeros((numTrigs, numklinpts))
    
    # HP filter to get rid of LF components
    # filterCutoff = procOpts.dispersion_PD_HPfilterCutoff

    # (b, a) = scipy.signal.butter(2, filterCutoff, 'highpass')
    lpfc = PDfilterCutoffs[0]
    hpfc = PDfilterCutoffs[1]
    Wn = [hpfc, lpfc]
    (b, a) = scipy.signal.butter(2, Wn=Wn, btype='bandpass')
    (b2, a2) = scipy.signal.butter(2, 0.4, 'lowpass')
    # subtract background
    sigdata = np.real(interpPD[:, 0:numklinpts])
    if pd_background is not None:
        bg = np.tile(pd_background, (numTrigs, 1))
        sigdata = sigdata - bg
        
    for n in range(0, numTrigs):
        sig = sigdata[n, :]
        #sig = scipy.signal.lfilter(b, a, sig)
        
        # filter signal perfectly by using FFT-IFFT method
        fftsig = np.fft.rfft(sig)
        ln = len(sig)
        # calculate indices in array by multipling by filter cuoffs
        idx0 = np.round(hpfc*ln)
        idx1 = np.round(lpfc*ln)
        if n == 0:
            DebugLog.log("DispersonProtocol.processData() len(sig)= %d len(fftsig)= %d idx0= %d idx1= %d" % (len(sig), len(fftsig), idx0, idx1))
            
        fftsig[0:idx0] = 0
        fftsig[idx1:] = 0
        sig = np.fft.irfft(fftsig, len(sig))
        
        # HP filter to get rid of LF components
        sig_ht = scipy.signal.hilbert(sig)
        
        #mag = np.complex(sig, np.imag(sig_ht))
        mag = np.abs(sig_ht)
        mag0 = mag[0]
        mag = mag - mag0
        mag = scipy.signal.lfilter(b2, a2, mag)
        mag = mag + mag0
#            mag = mag 
        # print("mag.shape= %s winFcn.shape= %s" % (repr(mag.shape), repr(winFcn.shape)))
        magWin[n, :] = winFcn / mag
        # magWin[n, :] = mag 
        
        ph = np.angle(sig_ht)
        ph_unwr = np.unwrap(ph)
        pcof = np.polyfit(k, ph_unwr, 1)
        fity = np.polyval(pcof, k)
        phaseCorr[n, :] = ph_unwr - fity

    magWin = np.mean(magWin, 0)
    
    # low pass filter to get rid of ripple
    filterCutoff = magWin_LPfilterCutoff
    (b, a) = scipy.signal.butter(2, filterCutoff, 'lowpass')
    magWin = scipy.signal.lfilter(b, a, magWin)
    
    # renomalze to 0...1
    minWin = np.min(magWin)
    maxWin = np.max(magWin)
    magWin = (magWin - minWin)/(maxWin - minWin)
    
    #magWin = magWin[0, :]
    phaseCorr = np.mean(phaseCorr, 0)
    dispData.magWin = magWin
    dispData.phaseCorr = phaseCorr
    
    winFcnTile = np.tile(winFcn, numTrigs).reshape(numTrigs, numklinpts)
    
    # make uncorrected aline
    fft_sig_u = np.fft.fft(winFcnTile * sigdata, 2048, 1)
    dispData.uncorrAline = 20*np.log10(np.mean(np.abs(fft_sig_u) + 1, 0))
    
    # make orrected aline
    magWinTile = np.tile(magWin, numTrigs).reshape(numTrigs, numklinpts)
    phaseCorrTile = np.tile(phaseCorr, numTrigs).reshape(numTrigs, numklinpts)
    sigDataComplex = magWinTile * sigdata * (np.cos(-phaseCorrTile) + 1j * np.sin(-phaseCorrTile))
    fft_sig_c = np.fft.fft(sigDataComplex, 2048, 1)
    dispData.corrAline = 20*np.log10(np.mean(np.abs(fft_sig_c) + 1, 0))
    
    dispData.numTrigs = numTrigs
    dispData.numklinpts = numklinpts
    dispData.magWin_LPfilterCutoff = magWin_LPfilterCutoff
    dispData.PDfilterCutoffs = PDfilterCutoffs
    dispData.phDiode_background = pd_background
    
    return dispData

def processUniqueDispersion(pd_data, dispData, pd_background=None, bg_collect=False):
    numTrigs = pd_data.shape[0]
    numklinpts = pd_data.shape[1]
    DebugLog.log("Dispersion.processUniqueDispersion: numklinpts= %d" % numklinpts)
    winFcn = np.hanning(numklinpts)

    pd = np.mean(pd_data, 0)
    dispData.phDiode = pd
    
    if bg_collect:
        pd = np.real(pd_data[:, 0:numklinpts])
        pd_background = np.mean(pd, 0)
    
    
    magWin = np.zeros((numTrigs, numklinpts))
    phaseCorr = np.zeros((numTrigs, numklinpts))
    
    # set up filter coefficients
    idx0 = np.round(dispData.PDfilterCutoffs[1]*numklinpts)
    idx1 = np.round(dispData.PDfilterCutoffs[0]*numklinpts)            
    (b, a) = scipy.signal.butter(2, dispData.magWin_LPfilterCutoff, 'lowpass')
    
    # subtract background
    sigdata = np.real(pd_data[:, 0:numklinpts])
    if pd_background is not None:
        bg = np.tile(pd_background, (numTrigs, 1))
        sigdata = sigdata - bg
        
    for n in range(0, numTrigs):
        sig = sigdata[n, :]
        
        # first filter the interferograph using an ideal filter with the FFT-IFFT method
        fftsig = np.fft.rfft(sig)
        fftsig[0:idx0] = 0
        fftsig[idx1:] = 0
        sig = np.fft.irfft(fftsig, len(sig))
        
        """
        Hilbert Transform
        A perfect reflector like a mirror should make the interferogram like a
        sine wave. The phase of the Hilbert Transform of a sine wave should be
        a linear ramp. Any deviation from this must therefore be from dispersion.
        """
        sig_ht = scipy.signal.hilbert(sig)
        mag = np.abs(sig_ht)
        magWin[n, :] = winFcn / mag       
        ph = np.angle(sig_ht)
        ph_unwr = np.unwrap(ph)
        phaseCorr[n, :] = scipy.signal.detrend(ph_unwr)

    
    magWin = np.mean(magWin, 0)
    magWin = scipy.signal.lfilter(b, a, magWin)
    minWin = np.min(magWin)
    maxWin = np.max(magWin)
    dispData.magWin = (magWin - minWin)/(maxWin - minWin)
    phaseCorr = np.mean(phaseCorr, 0)
    dispData.phaseCorr = phaseCorr
    dispData.phDiode_background = pd_background
    
    winFcnTile = np.tile(winFcn, numTrigs).reshape(numTrigs, numklinpts)
    
    # make uncorrected aline
    fft_sig_u = np.fft.fft(winFcnTile * sigdata, 2048, 1)
    dispData.uncorrAline = 20*np.log10(np.mean(np.abs(fft_sig_u) + 1, 0))
    
    DebugLog.log("Dispersion.processUniqueDispersion: len(ph_unwr)= %d phaseCorr.shape= %s" % (len(ph_unwr), repr(phaseCorr.shape)))
    DebugLog.log("Dispersion.processUniqueDispersion: magWin.shape= %s" % repr(magWin.shape))


    # make orrected aline
    magWinTile = np.tile(magWin, numTrigs).reshape(numTrigs, numklinpts)
    phaseCorrTile = np.tile(phaseCorr, numTrigs).reshape(numTrigs, numklinpts)
    sigDataComplex = magWinTile * sigdata * (np.cos(-phaseCorrTile) + 1j * np.sin(-phaseCorrTile))
    fft_sig_c = np.fft.fft(sigDataComplex, 2048, 1)
    dispData.corrAline = 20*np.log10(np.mean(np.abs(fft_sig_c) + 1, 0))
    
    return dispData
    
def saveDispersionData(dispData, basePath):
    fileName = 'dispersionLast'
    saveDir = os.path.join(basePath, 'Dispersion')
        
    filepath = os.path.join(saveDir, '%s.pickle' % fileName)
    f = open(filepath, 'wb')
    pickle.dump(dispData, f)
    f.close()
    
    filepath = os.path.join(saveDir, '%s.txt' % fileName)
    f = open(filepath, 'w')
    s = "%0.5f" % dispData.magWin[0]
    magWin = dispData.magWin[1:]
    for xi in magWin:
        s = s + ", %0.5f" % xi
    f.write(s + '\n')           
    
    s = "%0.5f" % dispData.phaseCorr[0]
    phaseCorr = dispData.phaseCorr[1:]
    for xi in phaseCorr:
        s = s + ", %0.5f" % xi
        
    f.write(s)           
    f.close()
    
def plotDispData(appObj, dispData, PDfiltCutoffs):
    pl = appObj.plot_disp_phdiode
    pl.clear()
    pl.plot(dispData.phDiode, pen='b')
    
    pl = appObj.plot_disp_magwin
    pl.clear()
    pl.plot(dispData.magWin, pen='b')

    pl = appObj.plot_disp_phasecorr            
    pl.clear()
    pl.plot(dispData.phaseCorr, pen='b')
    
    pl = appObj.plot_disp_uncorrected_aline
    pl.clear()
    pl.plot(dispData.uncorrAline, pen='b')
    x1 = 2*1024*PDfiltCutoffs[0]
    x1 = np.array([x1, x1])
    x2 = 2*1024*PDfiltCutoffs[1]
    x2 = np.array([x2, x2])
    y = np.array([np.min(dispData.uncorrAline), np.max(dispData.uncorrAline)])
    pl.plot(x1, y, pen='r')
    pl.plot(x2, y, pen='r')
    
    pl = appObj.plot_disp_corrected_aline
    pl.clear()
    pl.plot(dispData.corrAline, pen='b')    
    
# collect dispersion data
# appObj is an OCTWindowClass as defined in OCTnew.py
def runDispersion(appObj):
    DebugLog.log("runDispersion")
    appObj.tabWidget.setCurrentIndex(5)
    appObj.doneFlag = False
    appObj.isCollecting = True
    # trigRate = octfpga.GetTriggerRate()
    mirrorDriver = appObj.mirrorDriver
    
    # set the mirror position to (0,0)
    chanNames = [mirrorDriver.X_daqChan, mirrorDriver.Y_daqChan]
    data = np.zeros(2)
    if not appObj.oct_hw.IsDAQTestingMode():
        from DAQHardware import DAQHardware
        daq = DAQHardware()
        daq.writeValues(chanNames, data)
    
    pd_background = None
    
    fpgaOpts = appObj.oct_hw.fpgaOpts
    numklinpts = fpgaOpts.numKlinPts
    if fpgaOpts.InterpDownsample > 0:
        numklinpts =  numklinpts // 2
    
    # keep looping until we are signlaed to stop by GUI (flag set in appObj)
    try:
        frameNum = 0
        saveDirInit = False
        testDataDir = os.path.join(appObj.basePath, 'exampledata', 'Dispersion')
        dispData = DispersionData(fpgaOpts)
        klin = None # initialze klin to None so it will be computed first iteration

        while not appObj.doneFlag: 
            # setup and grab the OCT data - this will also fire the mirror output
            numTrigs = appObj.disp_numTrigs_spinBox.value()
            processMode = OCTCommon.ProcessMode(appObj.processMode_comboBox.currentIndex())
                
            # get proessing optiosn from GUI        
            PD_LP_fc = appObj.disp_pd_lpfilter_cutoff_dblSpinBox.value()
            PD_HP_fc = appObj.disp_pd_hpfilter_cutoff_dblSpinBox.value()
            PDfiltCutoffs = [PD_LP_fc, PD_HP_fc]
            magWin_LPfilterCutoff = appObj.disp_magwin_lpfilter_cutoff_dblSpinBox.value()
            dispData.mziFilter = appObj.mziFilter.value()
            dispData.magWin_LPfilterCutoff = magWin_LPfilterCutoff
            dispData.PDfilterCutoffs = PDfiltCutoffs
            
            collectBG = appObj.disp_collectBG_pushButton.isChecked()
            
            pd_background = dispData.phDiode_background                
    
            if processMode == OCTCommon.ProcessMode.FPGA:
                if appObj.oct_hw.IsOCTTestingMode():
                    pd_data = OCTCommon.loadRawData(testDataDir, frameNum % 19, dataType=1)
                    numklinpts = 1400
                else:
                    err, pd_data = appObj.oct_hw.AcquireOCTDataInterpPD(numTrigs)
                    DebugLog.log("runBScan(): AcquireOCTDataInterpPD() err = %d" % err)
                    # process the data
                    dispData = processData(pd_data, dispData, numklinpts, PDfiltCutoffs, magWin_LPfilterCutoff, pd_background, collectBG)

            elif processMode == OCTCommon.ProcessMode.SOFTWARE:
                if appObj.oct_hw.IsOCTTestingMode():
                    ch0_data,ch1_data=JSOraw.getSavedRawData(numTrigs,appObj.dispData.requestedSamplesPerTrig,appObj.savedDataBuffer)
                else:
                    # def AcquireOCTDataRaw(self, numTriggers, samplesPerTrig=-1, Ch0Shift=-1, startTrigOffset=0):
                    samplesPerTrig = fpgaOpts.SamplesPerTrig*2 + fpgaOpts.Ch0Shift*2
                    err, ch0_data,ch1_data = appObj.oct_hw.AcquireOCTDataRaw(numTrigs, samplesPerTrig)

                pdData,mziData,actualSamplesPerTrig = JSOraw.channelShift(ch0_data,ch1_data,dispData)    # shift the two channels to account for delays in the sample data compared to the MZI data 
                mzi_hilbert, mzi_mag, mzi_ph, k0 = JSOraw.processMZI(mziData, dispData)                # calculate k0 from the phase of the MZI data
                k0Cleaned = JSOraw.cleank0(k0, dispData) # Adjust the k0 curves so that the unwrapping all starts at the same phase    
                pd_data, klin = JSOraw.processPD(pdData, k0Cleaned, dispData, klin)  # Interpolate the PD data based upon the MZI data
                dispData.Klin = klin
                dispData = processUniqueDispersion(pd_data, dispData, pd_background, collectBG)
            else:
                QtGui.QMessageBox.critical (appObj, "Error", "Unsuppoted processing mode for current hardware")
                
            # plot the data
            plotDispData(appObj, dispData, PDfiltCutoffs)
            
            if appObj.getSaveState():
                saveDispersionData(dispData, appObj.settingsPath)
            
                saveOpts = appObj.getSaveOpts()
                if saveOpts.saveRaw:
                    if not saveDirInit:
                        saveDir = OCTCommon.initSaveDir(saveOpts, 'Dispersion')
                        saveDirInit = True
                    
                    OCTCommon.saveRawData(pd_data, saveDir, frameNum, dataType=1)
                appObj.dispData = dispData
                
            frameNum += 1
            QtGui.QApplication.processEvents() # check for GUI events, particularly the "done" flag
    except Exception as ex:
        # raise ex
        traceback.print_exc(file=sys.stdout)
        QtGui.QMessageBox.critical (appObj, "Error", "Error during scan. Check command line output for details")
    finally:
        appObj.isCollecting = False
        QtGui.QApplication.processEvents() # check for GUI events
        appObj.finishCollection()        
    