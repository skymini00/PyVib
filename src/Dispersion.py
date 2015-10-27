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

class DispersionData:
    def __init__(self):
        self.magWin = []
        self.phaseCorr = []
        self.phDiode = []
        self.uncorrAline = []
        self.corrAline = []
        self.phDiode_background = None

# process dispersion data
# interpPD = 2D numpy array, trigger# is first indiex, ptNum is second
# numklinpts = num klinpts
# PDfilterCutoffs = array like , 2 elements, with low/high filter cuoffs for the photodiode dieata, should be between 0 and 0.5
# magWin_LPfilterCutoff = low pass filter cutoffs for the magnitude window, hould be between 0 and 0.5
def processData(interpPD, numklinpts, PDfilterCutoffs, magWin_LPfilterCutoff, pd_background=None, bg_collect=False):
    # scanP = self.scanParams

    numTrigs = interpPD.shape[0]
    numPts = interpPD.shape[1]
    DebugLog.log("Dispersion.processData numTrigs= %d numPts= %d" % (numTrigs, numPts))
    
    dispData = DispersionData()
    
    dispData.magWin = []
    dispData.phaseCorr = []
    pd = np.mean(interpPD, 0)
    dispData.phDiode = pd
    dispData.uncorrAline = []
    dispData.corrAline = []
    
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
        while not appObj.doneFlag: 
            # setup and grab the OCT data - this will also fire the mirror output
            numTrigs = appObj.disp_numTrigs_spinBox.value()
            
            if appObj.oct_hw.IsOCTTestingMode():
                pd_data = OCTCommon.loadRawData(testDataDir, frameNum % 19, dataType=1)
                numklinpts =  1400
            else:
                err, pd_data = appObj.oct_hw.AcquireOCTDataInterpPD(numTrigs)
                DebugLog.log("runBScan(): AcquireOCTDataInterpPD() err = %d" % err)
    
            # get proessing optiosn from GUI        
            PD_LP_fc = appObj.disp_pd_lpfilter_cutoff_dblSpinBox.value()
            PD_HP_fc = appObj.disp_pd_hpfilter_cutoff_dblSpinBox.value()
            PDfiltCutoffs = [PD_LP_fc, PD_HP_fc]
            magWin_LPfilterCutoff = appObj.disp_magwin_lpfilter_cutoff_dblSpinBox.value()
            
            collectBG = appObj.disp_collectBG_pushButton.isChecked()
            
            # process the data
            dispData = processData(pd_data, numklinpts, PDfiltCutoffs, magWin_LPfilterCutoff, pd_background, collectBG)
            pd_background = dispData.phDiode_background
                
            # plot the data
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
    
            pl = appObj.plot_disp_corrected_aline
            pl.clear()
            pl.plot(dispData.corrAline, pen='b')
            
            if appObj.getSaveState():
                saveDispersionData(dispData, appObj.settingsPath)
            
                saveOpts = appObj.getSaveOpts()
                if saveOpts.saveRaw:
                    if not saveDirInit:
                        saveDir = OCTCommon.initSaveDir(saveOpts, 'Dispersion')
                        saveDirInit = True
                    
                    OCTCommon.saveRawData(pd_data, saveDir, frameNum, dataType=1)
            
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
    