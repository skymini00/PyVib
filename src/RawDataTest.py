# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 12:25:50 2015

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
import copy

def runRawDataTest(appObj):
    try:
        appObj.doneFlag = False
        while not appObj.doneFlag:
            numTrigs = appObj.rd_numTrigs_spinBox.value()
            samplesPerTrig = appObj.rd_samplesPerTrig_spinBox.value()
            err, ch0_data, ch1_data = appObj.oct_hw.AcquireOCTDataRaw(numTrigs, samplesPerTrig)

            pd_data = ch0_data
            mzi_data = ch1_data
            
                
            sampleOffset = appObj.sampleOffset_spinBox.value()*2
            ch0shift = appObj.ch0shift_spinBox.value()*2
            pd_data = pd_data[:, sampleOffset+ch0shift:]
            mzi_data = mzi_data[:, sampleOffset:-ch0shift]
            mzi_data = copy.copy(mzi_data)
            
            (b, a) = scipy.signal.butter(2, 0.005, 'highpass')
            mzi_data=scipy.signal.lfilter(b, a, mzi_data,axis=-1)    

            appObj.pd_plot.clear()
            appObj.mzi_plot.clear()
            appObj.mzi_phase_plot.clear()
            appObj.k0_plot.clear()
            
            nTrigs = min(numTrigs, 10)
            for i in range(0, nTrigs):
                pn = (i, nTrigs)
                
                appObj.pd_plot.plot(pd_data[i, :], pen=pn)
                appObj.mzi_plot.plot(mzi_data[i, :], pen=pn)
            
            mzi_hilbert = scipy.signal.hilbert(mzi_data, axis=-1)
            mzi_mag = np.abs(mzi_hilbert)
            mzi_ph = np.angle(mzi_hilbert)
            k0 = np.unwrap(mzi_ph,axis=-1)    

            mzi_ph = mzi_ph[:, 15:]
            nTrigs = min(numTrigs, 50)                
            for i in range(0, nTrigs):
                pn = (i, nTrigs)
                
                appObj.mzi_phase_plot.plot(mzi_ph[i, 0:100], pen=pn)
                appObj.k0_plot.plot(k0[i, 0:100], pen=pn)
            
            QtGui.QApplication.processEvents() # check for GUI events, particularly the "done" flag
    except Exception as ex:
        traceback.print_exc(file=sys.stdout)
        QtGui.QMessageBox.critical (appObj, "Error", "Error during scan. Check command line output for details")
    finally:
        appObj.isCollecting = False
        QtGui.QApplication.processEvents() # check for GUI events
        appObj.finishCollection()        

        