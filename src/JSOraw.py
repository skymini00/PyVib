# -*- coding: utf-8 -*-
"""
Created on Fri Sep  4 10:15:58 2015

@author: OHNS
"""
import sys
import os
from PyQt4 import QtCore, QtGui, uic
import pyqtgraph as pg
import scipy.signal 
import numpy as np
import datetime
import matplotlib.pyplot as plt
from DebugLog import DebugLog

import OCT_testimport1 as OCTRaw                   # to import test data
#import OCTRawDataHardwareInterface as OCTRaw       # to import real data


class blankClass:
    def __init__(self):
        pass

def processMZI(mzi_data):
    filter=0    
    if filter==1:       # filtering seems to reduce sidebands caused by interpolation
        (b, a) = scipy.signal.butter(2, 0.005, 'highpass')
        mzi_data=scipy.signal.lfilter(b, a, mzi_data,axis=-1)    
    
    mzi_hilbert = scipy.signal.hilbert(mzi_data, axis=-1)
    mzi_mag = np.abs(mzi_hilbert)
    mzi_ph = np.angle(mzi_hilbert)
    mzi_hilbert = np.imag(mzi_hilbert)
    k0 = np.unwrap(mzi_ph,axis=-1)    
    return mzi_hilbert, mzi_mag, mzi_ph, k0
    
def processPD(pd_data, k0, klin_idx=[15, 2030], numklinpts=2000, klin=None):
    if klin is None:
        klin = np.linspace(k0[0,klin_idx[0]], k0[0,klin_idx[1]], numklinpts)
    pd_interp=np.zeros([pd_data.shape[0],numklinpts])
    for i in range(pd_data.shape[0]):
        pd_interp[i,:] = np.interp(klin, k0[i,:], pd_data[i,:])    
    return pd_interp, klin  

def dispersionCorrection(pd,dispMode='None'):
    numKlinPts=pd.shape[1]
    if dispMode=='Hanning':
        windowFunctionMag = 2*np.hanning(numKlinPts)
        windowFunctionPh = np.zeros(numKlinPts)
    elif dispMode=='Hamming':
        windowFunctionMag = 2*np.hamming(numKlinPts)
        windowFunctionPh = np.zeros(numKlinPts)
    elif dispMode=='Unique':
        windowFunctionMag,windowFunctionPh = processUniqueDispersion(pd)
    else:
        windowFunctionMag = np.ones(numKlinPts)
        windowFunctionPh = np.zeros(numKlinPts)
    pd_interpDisp = np.abs(windowFunctionMag * pd * (np.cos(-1*windowFunctionPh) + 1j * np.sin(-1*windowFunctionPh)))
#    pd_interpDisp = windowFunctionMag * pd
    return pd_interpDisp, windowFunctionMag, windowFunctionPh

def calculateAline(pd):
    numPts=pd.shape[1] 
    pd_fft = np.fft.fft(pd, n=2048,axis=-1)/numPts
    alineMag = 20*np.log10(np.abs(pd_fft) + 1)
    alinePhase=np.unwrap(np.angle(pd_fft),axis=-1)
    return pd_fft, alineMag, alinePhase

def channelShift(ch0_data,ch1_data,numShiftPts):
    actualSamplesPerTrig=ch0_data.shape[1]-numShiftPts
    pdData=ch0_data[:,numShiftPts:]
    mziData=ch1_data[:,0:actualSamplesPerTrig]
    return pdData,mziData,actualSamplesPerTrig

def processUniqueDispersion(pd_data):
    numTrigs = pd_data.shape[0]
    numklinpts = pd_data.shape[1]
    
    (b2, a2) = scipy.signal.butter(2, 0.6, 'lowpass')
    (b, a) = scipy.signal.butter(2, 0.01, 'highpass')
    pd_dataFilt=scipy.signal.lfilter(b, a, pd_data,axis=-1) #this smooths out the hilbert phase plot
    pd_hilbert=scipy.signal.hilbert(pd_dataFilt, axis=-1)
    pd_hilbertMag = scipy.signal.lfilter(b2, a2, np.abs(pd_hilbert), axis=-1)
    pd_hilbertPh=np.unwrap(np.angle(pd_hilbert), axis=-1)    
    phaseDetrend=scipy.signal.detrend(pd_hilbertPh,axis=-1)
    
    winFcn = np.hanning(numklinpts)  
    magWin=winFcn/np.mean(pd_hilbertMag,0)
    minWin = np.amin(magWin)    # renomalze to 0...1
    maxWin = np.amax(magWin)
    magWinOut = (magWin - minWin)/(maxWin - minWin)    
    phaseCorrOut=np.mean(phaseDetrend,0)
    
    plt.figure(1)
    plt.clf()
    plotColors = ['-r', '-b', '-g', '-c', '-m', '-y', '-k']
    for n in range(0, 1):
        clr = plotColors[n % len(plotColors)]
        plt.subplot(3, 1, 1)
        plt.plot(pd_hilbertMag[n, :], clr)
        plt.subplot(3, 1, 2)
        plt.plot(pd_hilbertPh[n, :], clr)
        plt.subplot(3, 1, 3)
        plt.plot(phaseDetrend[n, :], clr)
    plt.subplot(3, 1, 1)
    plt.title("hilbert magnitude after filter")
    plt.subplot(3, 1, 2)
    plt.title("hilbert phase")
    plt.subplot(3, 1, 3)
    plt.title("hilbert phase after detrend")

    return magWinOut, phaseCorrOut


def processUniqueDispersion_old(pd_data, PD_HPfilterCutoff=0.05, magWin_LPfilterCutoff = 0.05):
    numTrigs = pd_data.shape[0]
    numklinpts = pd_data.shape[1]
    
    winFcn = np.hanning(numklinpts)
    k = np.linspace(0, numklinpts, numklinpts)    
    magWin = np.zeros((numTrigs, numklinpts))
    phaseCorr = np.zeros((numTrigs, numklinpts))
    
    (b, a) = scipy.signal.butter(2, PD_HPfilterCutoff, 'highpass')
    (b2, a2) = scipy.signal.butter(2, 0.4, 'lowpass')
    for n in range(0, numTrigs):
        sig = pd_data[n, :]
        sig = scipy.signal.lfilter(b, a, sig)
        sig_ht = scipy.signal.hilbert(sig)        
        mag = np.abs(sig_ht)
        mag0 = mag[0]
        mag = mag - mag0
        mag = scipy.signal.lfilter(b2, a2, mag)
        mag = mag + mag0
        magWin[n, :] = winFcn / mag        
        ph = np.angle(sig_ht)
        ph_unwr = np.unwrap(ph)
        pcof = np.polyfit(k, ph_unwr, 1)
        fity = np.polyval(pcof, k)
        phaseCorr[n, :] = ph_unwr - fity
    magWin = np.mean(magWin, 0)    
    (b, a) = scipy.signal.butter(2, magWin_LPfilterCutoff, 'lowpass')
    magWin = scipy.signal.lfilter(b, a, magWin)    
    minWin = np.min(magWin)    # renomalze to 0...1
    maxWin = np.max(magWin)
    magWin = (magWin - minWin)/(maxWin - minWin)    
    phaseCorr = np.mean(phaseCorr, 0)       
    return magWin, phaseCorr








       
  
def getSavedRawData(self,numTrigs,requestedSamplesPerTrig):
    # oct_data is a 2D arary of complex numbers
    # When used to carry raw data:
    # Photodiode data (interferogram from the sample) is encoded as imaginary part 
    # MZI data (interferogram from MZI) is encoded as imaginary part        print(self.count,self.oct_data_all.shape[0])
#        print('size/numTrigs',self.oct_data_file.shape[0],numTrigs)        
    ch0_data=np.zeros([numTrigs,self.ch1_data_file.shape[1]])        
    ch1_data=np.zeros([numTrigs,self.ch1_data_file.shape[1]])        
    self.requestedSamplesPerTrig.setValue(self.ch1_data_file.shape[1])
   
    for i in range(numTrigs):
        if self.count==self.ch0_data_file.shape[0]:
            self.count=0           
        ch0_data[i,:]=self.ch0_data_file[self.count,:] 
        ch1_data[i,:]=self.ch1_data_file[self.count,:] 
        self.count=self.count+1             
    return ch0_data, ch1_data
    
def getNewRawData(self,numTrigs,requestedSamplesPerTrig):                                 
    try:
        err = OCTRaw.ConfigAcq(numTrigs, requestedSamplesPerTrig)
        if err>0:
            print("Config Acq: err = %d" % err)
        (err, ch0_data, ch1_data, samplesRemaining, timeElapsed) = OCTRaw.GrabData(numTrigs, requestedSamplesPerTrig)
        if err>0:
            print("Grab Data: err = %d" % err)
            print("Grab Data: ch0_data.shape = %s ch1_data.shape= %s" % (repr(ch0_data.shape), repr(ch1_data.shape)))
            print("Grab Data: Samples Remaing= %d TimeElapsed= %d" % (samplesRemaining, timeElapsed))
    except Exception as ex:
        print('Error collecting data')
        raise ex
        err = OCTRaw.CloseFPGA()
        if err>0:
            print("Close FPGA: err = %d" % err)
    
    if self.saveData_checkBox.isChecked()==True:      # get new data and save to disk for later use
        outfile='testData.npz'
        np.savez_compressed(outfile, ch0_data=ch0_data, ch1_data=ch1_data)
        print('saved data in :',outfile)
        self.saveData_checkBox.setChecked(False)                     
    return ch0_data, ch1_data
    
def rawDataTest_pushButton_on(self):                        
    self.doneFlag = False
           
    peakXPos=np.array([0],dtype=int)       
    peakYPos=np.array([0],dtype=float)       
    print('peak-init',peakXPos,peakYPos)            
               
    while self.doneFlag == False:
        # read data analysis settings from the GUI
        numTrigs=self.numTrig.value()
        requestedSamplesPerTrig=self.requestedSamplesPerTrig.value()
        startSample=self.startSample.value()
        endSample=self.endSample.value()
        klin_idx = [startSample, endSample]
        numKlinPts=self.numKlinPts.value()
        numShiftPts=self.numShiftPts.value()
        numTrigs=self.numTrig.value()
        dispCode=self.dispersionCompAlgorithm_comboBox.currentIndex()            
        if dispCode==0:
            dispMode='None'
        elif dispCode==1:
            dispMode='Hanning'
        elif dispCode==2:
            dispMode='Hamming'
        elif dispCode==3:
            dispMode='Unique'
        else:
            dispMode='None'
                 
        # Get data using one of several methods
        if self.getDataMethod=='saved raw data':
            ch0_data,ch1_data=self.getSavedRawData(numTrigs,requestedSamplesPerTrig)
        elif self.getDataMethod=='new raw data':
            ch0_data,ch1_data=self.getNewRawData(numTrigs,requestedSamplesPerTrig)
        
        # delay the MZI to account for it having a shorter optical path than the sample/reference arm path
        pdData,mziData,actualSamplesPerTrig=channelShift(ch0_data,ch1_data,numShiftPts)    
        textString='Actual samples per trigger: {actualSamplesPerTrig}'.format(actualSamplesPerTrig=actualSamplesPerTrig)            
        self.actualSamplesPerTrig_label.setText(textString)         
        
        # Process the data        
        mzi_hilbert, mzi_mag, mzi_ph, k0 = processMZI(mziData) 
        pd_interpRaw, klin = processPD(pdData, k0, klin_idx, numKlinPts)   
        pd_interpDispComp,windowFunctionMag,windowFunctionPh = dispersionCorrection(pd_interpRaw,dispMode)
        self.windowFunctionMag=windowFunctionMag       
        self.windowFunctionPh=windowFunctionPh       
        pd_fftNoInterp, alineMagNoInterp, alinePhaseNoInterp = calculateAline(pdData[:,startSample:endSample])
        pd_fftRaw, alineMagRaw, alinePhaseRaw = calculateAline(pd_interpRaw)
        pd_fftDispComp, alineMagDispComp, alinePhaseDispComp = calculateAline(pd_interpDispComp)
       
        #scale k0 and the MZI to the same range to plot them so they overlap
        k0Ripple= scipy.signal.detrend(k0[0,500:700],axis=-1)
        k0RippleNorm=k0Ripple/k0Ripple.max()
        mziDataRipple= scipy.signal.detrend(mziData[0,500:700],axis=-1)
        mziDataNorm=mziDataRipple/mziDataRipple.max()
       
        # Find the peak of the A-line within a range and calculate the phase noise
        rangePeak=[100, 900]
        alineAve=np.average(alineMagDispComp,axis=0) 
#            alineAve[501]=1.1
        peakXPos[0]=np.argmax(alineAve[rangePeak[0]:rangePeak[1]])+rangePeak[0]
        peakYPos[0]=alineAve[peakXPos[0]]              
        phaseNoiseTD=np.unwrap(alinePhaseDispComp[:,peakXPos[0]])
        phaseNoiseTD=phaseNoiseTD-phaseNoiseTD[0]
        phaseNoiseFFT = np.fft.fft(phaseNoiseTD, n=2048)
        phaseNoiseFD = 20*np.log10(np.abs(phaseNoiseFFT) + 1)        
                    
        # Clear all of the plots
        self.mzi_plot.clear() 
        self.pd_plot.clear()
        self.mzi_mag_plot.clear()
        self.mzi_phase_plot.clear()
        self.k0_plot.clear()
        self.interp_pdRaw_plot.clear()
        self.interp_pdDispComp_plot.clear()
        self.alineNoInterp_plot.clear()
        self.alineRaw_plot.clear()
        self.alineDispComp_plot.clear()
        self.phaseNoiseTD_plot.clear()
        self.phaseNoiseFD_plot.clear()
        self.dispWnfcMag_plot.clear()
        self.dispWnfcPh_plot.clear()
        
        # Plot all the data
        if self.plotFirstOnly_checkBox.isChecked()==True:
            i=0
            self.pd_plot.plot(pdData[i,:], pen='r')            
            self.mzi_plot.plot(mziData[i,:], pen='r')            
            self.mzi_mag_plot.plot(mzi_mag[i,:], pen='r')            
            self.k0_plot.plot(k0[i,:], pen='r')
            sampleNum=np.linspace(startSample,endSample,numKlinPts)
            self.k0_plot.plot(sampleNum,klin, pen='b')                      
            self.interp_pdRaw_plot.plot(pd_interpRaw[i,:], pen='r')           
            self.interp_pdDispComp_plot.plot(pd_interpDispComp[i,:], pen='r')           
            self.alineNoInterp_plot.plot(alineMagNoInterp[i,:], pen='r')
            self.alineRaw_plot.plot(alineMagRaw[i,:], pen='r')
            self.alineDispComp_plot.plot(alineMagDispComp[i,:], pen='r')
        else:                   
            for i in range(numTrigs):
                self.pd_plot.plot(pdData[i,:], pen=(i,numTrigs))            
                self.mzi_plot.plot(mziData[i,:], pen=(i,numTrigs))            
                self.mzi_mag_plot.plot(mzi_mag[i,:], pen=(i,numTrigs))            
                self.mzi_phase_plot.plot(mzi_ph[i,:], pen=(i,numTrigs))            
                self.k0_plot.plot(k0[i,:], pen=(i,numTrigs))                      
                self.interp_pdRaw_plot.plot(pd_interpRaw[i,:], pen=(i,numTrigs))           
                self.interp_pdDispComp_plot.plot(pd_interpDispComp[i,:], pen=(i,numTrigs))           
                self.alineNoInterp_plot.plot(alineMagNoInterp[i,:], pen=(i,numTrigs))
                self.alineRaw_plot.plot(alineMagRaw[i,:], pen=(i,numTrigs))
                self.alineDispComp_plot.plot(alineMagDispComp[i,:], pen=(i,numTrigs))
            
        self.alineDispComp_plot.plot(peakXPos,peakYPos, pen=None, symbolBrush='k', symbolPen='b')
        self.phaseNoiseTD_plot.plot(phaseNoiseTD, pen='r')
        self.phaseNoiseFD_plot.plot(phaseNoiseFD, pen='r')
        self.mzi_phase_plot.plot(mziDataNorm, pen='b')            
        self.mzi_phase_plot.plot(k0RippleNorm, pen='r')            
        
        # if you want to align the pd and the Mzi data
#            plotPDPhase.plot(pdData[0,:], pen='r')
#            plotPDPhase.plot(mziData[0,:], pen='b')

        self.dispWnfcMag_plot.plot(windowFunctionMag, pen='b')
        self.dispWnfcPh_plot.plot(windowFunctionPh, pen='b')
        
        
#            # Now create a bscan image from the 1 aline, but sweep the shift value between the mzi and pd to see what works best
#            nShift=201        
#            bscan=np.zeros([nShift, alineMag.shape[0]])        
#            for i in range(nShift):
#                shift=i-(nShift-1)/2
#                if shift<0:
#                    mzi_data_temp=mzi_data[-1*shift:]
#                    pd_data_temp=pd_data[0:mzi_data_temp.shape[0]]
#    #                print(mzi_data_temp.shape,pd_data_temp.shape)
#                elif shift>0:
#                    pd_data_temp=pd_data[shift:]
#                    mzi_data_temp=mzi_data[0:pd_data_temp.shape[0]]
#    #                print(mzi_data_temp.shape,pd_data_temp.shape)
#                elif shift==0:
#                    pd_data_temp=pd_data
#                    mzi_data_temp=mzi_data
#                    
#                mzi_hilbert, mzi_mag, mzi_ph, k0 = processMZI(mzi_data_temp)
#                pd_interpRaw, pd_interpHanning, pd_fft, alineMag, alinePhase, klin = processPD(pd_data_temp, k0, klin_idx, numklinpts)
#                bscan[i,:]=alineMag
#    
#            pl = self.bscan_plot
#            pl.setImage(bscan)            
        
#            print('alineMagDispComp ',alineMagDispComp.shape)
        self.bscan_plot.setImage(alineMagDispComp)
        QtGui.QApplication.processEvents() # check for GUI events                       print('finished loop')
                   
def rawDataTest_pushButton_clicked(self):  # CtoF button event handler
    if self.rawDataTest_pushButton.isChecked():          
        self.rawDataTest_pushButton_on() # start collecting if button is turned on
    else:
        self.doneFlag = True  # stop collecting if button is turned off

def saveDispersion_pushButton_clicked(self):                
    outfile=datetime.datetime.now().strftime("dispComp-%Y_%m_%d-%H_%M_%S.npz")
    windowFunctionMag=self.windowFunctionMag
    windowFunctionPh=self.windowFunctionPh
    np.savez_compressed(outfile, windowFunctionMag=windowFunctionMag, windowFunctionPh=windowFunctionPh)
    print('saved data in :',outfile)       
    self.dispCompFilename_label.setText(outfile)
          
def runJSOraw(appObj):
    DebugLog.log("runJSOraw")
    appObj.tabWidget.setCurrentIndex(7)
    appObj.doneFlag = False
    appObj.isCollecting = True
    testDataDir = os.path.join(appObj.basePath, 'exampledata', 'JSOraw')
    print('test data dir',testDataDir)
    
    JSOrawInfo=blankClass()
#        self.getDataMethod='new raw data'
    JSOrawInfo.doneFlag=False
    JSOrawInfo.getDataMethod='saved raw data'
#        self.getDataMethod='new raw data'
    JSOrawInfo.singleProcess=True
    JSOrawInfo.laserSweepFreq=50000
    
    if JSOrawInfo.getDataMethod=='new raw data':      # get new data            
        err = OCTRaw.InitFPGA(0)
        if err>0:
            print("Init FPGA: err = %d" % err)

    elif JSOrawInfo.getDataMethod=='saved raw data':    #load in the saved data and use it instead
#            outfile='testData2.npz'
        outfile=os.path.join(testDataDir,'testData3.npz')
        print('outfile',outfile)
        x=np.load(outfile)             
        JSOrawInfo.ch0_data_file=x['ch0_data']
        JSOrawInfo.ch1_data_file=x['ch1_data']
        JSOrawInfo.count=0
        JSOrawInfo.saveRawData=0
        appObj.requestedSamplesPerTrig.setValue(JSOrawInfo.ch1_data_file.shape[1])
        print('loaded data from :',outfile)

 