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
import scipy.cluster
import scipy.ndimage
import numpy as np
import datetime
import matplotlib.pyplot as plt
import traceback
from DebugLog import DebugLog
import pickle
import Dispersion
import time

#class DispersionData:  # this class holds all the dispersion compensation data
#    def __init__(self):
#        self.magWin =[]
#        self.phaseCorr = []
#        self.phDiode_background = None
#        self.requestedSamplesPerTrig=[]
#        self.startSample=[]
#        self.endSample=[]
#        self.numKlinPts=[]
#        self.Klin=[]
#        self.numShiftPts=[]
#        self.filterWidth=[]
#        self.PDfilterCutoffs=[]
#        self.mziFilter=[]
#        self.magWin_LPfilterCutoff=[]
#        self.k0Reference=[]

                
class SavedDataBuffer: 
    def __init__(self,appObj):
        self.ch0_data_file=[]
        self.ch1_data_file=[]
        self.count=0
        self.saveRawData=0
        self.testDataDir = os.path.join(appObj.basePath, 'exampledata', 'JSOraw')
        
    def listRawDataFilenames(self,appObj):        
        appObj.rawDataFilename_comboBox.clear()
        count=0
        for entry in os.listdir(self.testDataDir):
            if entry.endswith('.npz'):
                appObj.rawDataFilename_comboBox.insertItem(count,entry)
                count=count+1
                
    def loadData(self,appObj):
        #testDataDir = os.path.join(appObj.basePath, 'exampledata', 'JSOraw')
        #outfile=os.path.join(testDataDir,'testData.npz')
        fileName=appObj.rawDataFilename_comboBox.currentText()
        outfile=os.path.join(self.testDataDir, fileName)
        x=np.load(outfile)                    
        self.ch0_data_file=x['ch0_data']
        self.ch1_data_file=x['ch1_data']
        self.count=0
        self.saveRawData=0
        appObj.requestedSamplesPerTrig.setValue(self.ch1_data_file.shape[1])
        print('loaded saved data into buffer from :',outfile)
    
    def saveData(self,appObj,dataToSave,fileName):
        dateName=datetime.datetime.now().strftime("-%Y_%m_%d-%H_%M_%S.npz")
        fName=fileName + dateName
        print('fName',fName)
        outfile=os.path.join(self.testDataDir,fName)
        ch0_data = dataToSave[0]
        ch1_data = dataToSave[1]
        np.savez_compressed(outfile, ch0_data=ch0_data, ch1_data=ch1_data)       
        print('data was saved to: ', fName)
        self.listRawDataFilenames(appObj)
        
def processMZI(mzi_data, dispData, FIRcoeff=[0+0j]):
    t1 = time.time()
    complex_FIR = 1
    if complex_FIR == 0:
        #normal way using fft/ifft approach to do the hilbert transform
    
        # filtering seems to reduce sidebands created during the interpolation process
        (b, a) = scipy.signal.butter(2, dispData.mziFilter, 'highpass')
        mzi_data=scipy.signal.lfilter(b, a, mzi_data,axis=-1)        
    
        # mean subtraction is much faster, and we only need to eliminate the zero-frequency componment
        # mzi_data = mzi_data - np.mean(mzi_data, keepdims=True)
    
        mzi_complex = scipy.signal.hilbert(mzi_data, axis=-1)
    elif complex_FIR==1:      
        # Approximate Hilbert Transform using FIR filter
        """
        # code to design the filter coefficients (This is not needed since they are pasted in below, but the code included for reference)
        filter_order = 35
        stop_freq=0.25
        pass_freq=0.15
        prt_lpf=scipy.signal.remez(filter_order, [0, pass_freq, stop_freq, 0.5], [1, 0])
        exp_len=np.arange(-filter_order/2,filter_order/2,1)
        FIRcoeff=prt_lpf*np.exp(np.complex(0,1)*2*np.pi*0.25*exp_len)
        
        for i in range(filter_order):
            print(FIRcoeff[i],',')
        """    
        # these filter coefficients come from the above code
        FIRcoeff=np.array([
        (-0.00030351655287-0.00030351655287j) ,
        (0.00113140195234-0.00113140195234j) ,
        (-9.35755583151e-05-9.35755583151e-05j) ,
        (0.00252252948249-0.00252252948249j) ,
        (0.00218108432482+0.00218108432482j) ,
        (0.0032263770108-0.0032263770108j) ,
        (0.00690123992871+0.00690123992871j) ,
        (0.0001322186872-0.0001322186872j) ,
        (0.0122523152659+0.0122523152659j) ,
        (-0.00969906300727+0.00969906300727j) ,
        (0.0128955480118+0.0128955480118j) ,
        (-0.0267326545724+0.0267326545724j) ,
        (0.000205842190372+0.000205842190372j) ,
        (-0.0472714248792+0.0472714248792j) ,
        (-0.0409251974438-0.0409251974438j) ,
        (-0.0643104412456+0.0643104412456j) ,
        (-0.212334405586-0.212334405586j) ,
        (0.282605062261-0.282605062261j) ,
        (0.212334405586+0.212334405586j) ,
        (-0.0643104412456+0.0643104412456j) ,
        (0.0409251974438+0.0409251974438j) ,
        (-0.0472714248792+0.0472714248792j) ,
        (-0.000205842190372-0.000205842190372j) ,
        (-0.0267326545724+0.0267326545724j) ,
        (-0.0128955480118-0.0128955480118j) ,
        (-0.00969906300727+0.00969906300727j) ,
        (-0.0122523152659-0.0122523152659j) ,
        (0.0001322186872-0.0001322186872j) ,
        (-0.00690123992871-0.00690123992871j) ,
        (0.0032263770108-0.0032263770108j) ,
        (-0.00218108432482-0.00218108432482j) ,
        (0.00252252948249-0.00252252948249j) ,
        (9.35755583151e-05+9.35755583151e-05j) ,
        (0.00113140195234-0.00113140195234j) ,
        (0.00030351655287+0.00030351655287j)
        ])
#        mzi_complex=np.zeros(mzi_data.shape, dtype=complex)
#        for i_th in range(0, mzi_data.shape[0]):
#            mzi_complex[i_th,:]=np.convolve(mzi_data[i_th,:],FIRcoeff, mode='same')  # this works but is slow

#        mzi_complex=scipy.ndimage.filters.convolve1d(mzi_data,FIRcoeff)  # This gives errors - doesn't work with complex numbers correctly

        mzi_complex=np.apply_along_axis(lambda m: np.convolve(m, FIRcoeff, mode='same'), axis=1, arr=mzi_data) # this works and is the fastest way to do it that I could figure out

        # mzi_complex=scipy.signal.lfilter(FIRcoeff,[1],mzi_data,axis=1)  # this works, but is also very slow
        
        
    DebugLog.log("JSOraw.processMZI() hilbert trasnform time= %0.4f" % (time.time() - t1))
    
    mzi_mag = np.abs(mzi_complex)
    t1 = time.time()    
    mzi_ph = np.angle(mzi_complex)
    DebugLog.log("JSOraw.processMZI() angle calc time= %0.4f" % (time.time() - t1))
        
    mzi_hilbert = np.imag(mzi_complex)
    
    t1 = time.time()    
    k0 = np.unwrap(mzi_ph,axis=-1)    
    DebugLog.log("JSOraw.processMZI() unwraptime= %0.4f" % (time.time() - t1))
    return mzi_hilbert, mzi_mag, mzi_ph, k0
 
def cleank0(k0,dispData):  # for Disp Compensation: look at all of the k0 tracings and shift them 2pi rads so they overlap
    k0Cleaned=np.copy(k0)    
    k0Init=k0Cleaned[:,dispData.startSample]
    # If there are phase jumps going on in the data, find the bad k0 curves and shift them appropriately by 2pi
    while (np.abs(np.max(k0Init)-np.min(k0Init)))>(1.5*np.pi): 
        k0InitUW=np.unwrap(k0Init)
        diff=np.abs(k0InitUW-k0Init)
        indexGrp1=np.argwhere(diff>(1.5*np.pi))
        indexGrp2=np.argwhere(diff<(1.5*np.pi))
        shift=np.sign(k0Init[indexGrp1[0]]-k0Init[indexGrp2[0]])*2*np.pi
        if len(indexGrp1)>=len(indexGrp2): 
            k0Cleaned[indexGrp2,:]=k0Cleaned[indexGrp2,:]+shift
        else:
            k0Cleaned[indexGrp1,:]=k0Cleaned[indexGrp1,:]-shift
        k0Init=k0Cleaned[:,dispData.startSample]
    dispData.k0Reference=np.mean(k0Cleaned,axis=0)              
    return k0Cleaned

def cleank0Run(k0,dispData):  # for Runtime: use the cleaned k0Reference data and adjust all new k0 data by 2pi rads to overlap it
    k0Cleaned=np.copy(k0)    
    k0Init=k0Cleaned[:,dispData.startSample]
    DebugLog.log('JSOraw.cleank0Run: len(k0reference)= %d startSample= %d' % (len(dispData.k0Reference), dispData.startSample))
    k0RefInit=np.tile(dispData.k0Reference[dispData.startSample],len(k0Init))    
    
    # If there are phase jumps going on in the data, find the bad k0 curves and shift them appropriately by 2pi
    diff=k0Init-k0RefInit
    while np.max(np.abs(diff))>(1.5*np.pi): 
        indexGrp1=np.argwhere(diff>(1.5*np.pi))
        indexGrp2=np.argwhere(diff<(-1.5*np.pi))        
        k0Cleaned[indexGrp1,:]=k0Cleaned[indexGrp1,:]-2*np.pi
        k0Cleaned[indexGrp2,:]=k0Cleaned[indexGrp2,:]+2*np.pi
        k0Init=k0Cleaned[:,dispData.startSample]
        diff=k0Init-k0RefInit
    return k0Cleaned
    
def processPD(pd_data, k0, dispData, klin=None):
    numklinpts = dispData.numKlinPts
    # this section is only run during the dispersion compensation algorithm. Then, klin is saved with the dispersion file and used from there on
    if klin is None or (len(klin) != numklinpts):   
        klin = np.linspace(k0[0,dispData.startSample], k0[0,dispData.endSample], numklinpts)

    # if num klinpts is > 2048, need to use downsampling to get interp points below 2048
    if numklinpts > 2048:
        dsf = numklinpts // 2048 + 1
        num_klin_ds = numklinpts // dsf
        DebugLog.log("JSOraw.processPD numklinpts= %d dsf= %d len(klin)= %d pd_data.shape= %s" % (numklinpts, dsf, len(klin), repr(pd_data.shape)))
        
        pd_interp = np.zeros((pd_data.shape[0], numklinpts // dsf))
        for i in range(pd_data.shape[0]):
            interp_pd = np.interp(klin, k0[i,:], pd_data[i,:])    
            interp_pd = np.reshape(interp_pd, (num_klin_ds, dsf))
            interp_pd = np.mean(interp_pd, 1)
            pd_interp[i,:] = interp_pd
    else:
        pd_interp=np.zeros((pd_data.shape[0], numklinpts))
        for i in range(pd_data.shape[0]):
            pd_interp[i,:] = np.interp(klin, k0[i,:], pd_data[i,:])    
        
    return pd_interp, klin  

def dispersionCorrection(pd,dispData):
    if dispData.dispMode=='Hanning':
        dispData.magWin = 2*np.hanning(dispData.numKlinPts)
        dispData.phaseCorr = np.zeros(dispData.numKlinPts)
    elif dispData.dispMode=='Hamming':
        dispData.magWin = 2*np.hamming(dispData.numKlinPts)
        dispData.phaseCorr = np.zeros(dispData.numKlinPts)
    elif dispData.dispMode=='Mirror': 
        # this routine is based upon the concept that a mirror should give a single sharp peak with the phase of the Hilbert transform being a linear ramp 
        processMirrorDispersion(pd,dispData)
    elif dispData.dispMode=='Brian':
        #Brian put a call in to your dispersion subroutine here
        pass
    else:
        dispData.magWin = np.ones(dispData.numKlinPts)
        dispData.phaseCorr = np.zeros(dispData.numKlinPts)
    
def calculateAline(pd, returnPhase=True, unwrapPhase=True):
    numPts=pd.shape[1] 

    t1 = time.time()    
    pd_fft = (np.fft.fft(pd, n=2048,axis=-1)/numPts)/100
    DebugLog.log("JSOraw.calculateAline() FFT time= %0.4f" % (time.time() - t1))
    
    t1 = time.time()    
    alineMag = np.abs(pd_fft)
    DebugLog.log("JSOraw.calculateAline() abs time= %0.4f" % (time.time() - t1))
    
    t1 = time.time()    
    alineMag = 20*np.log10(alineMag + 1)
    DebugLog.log("JSOraw.calculateAline() log10 time= %0.4f" % (time.time() - t1))
    
    alinePhase = None
    if returnPhase:
        alinePhase= np.angle(pd_fft)
        if unwrapPhase:
            alinePhase= np.unwrap(alinePhase,axis=-1)
        
    return pd_fft, alineMag, alinePhase

def channelShift(ch0_data,ch1_data,dispData):
    actualSamplesPerTrig=ch0_data.shape[1]-dispData.numShiftPts
    pdData=ch0_data[:,dispData.numShiftPts:]
    mziData=ch1_data[:,0:actualSamplesPerTrig]
    return pdData,mziData,actualSamplesPerTrig

def processMirrorDispersion(pd_data, dispData, pd_background=None, bg_collect=False):
    numTrigs = pd_data.shape[0]
    numklinpts = pd_data.shape[1]
    winFcn = np.hanning(numklinpts)
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
    dispData.phaseCorr = np.mean(phaseCorr, 0)
    
def getSavedRawData(numTrigs,requestedSamplesPerTrig,JSOrawSavedData):
    # oct_data is a 2D arary of complex numbers
    # When used to carry raw data:
    # Photodiode data (interferogram from the sample) is encoded as imaginary part 
    # MZI data (interferogram from MZI) is encoded as imaginary part        print(self.count,self.oct_data_all.shape[0])
#        print('size/numTrigs',self.oct_data_file.shape[0],numTrigs)        
    DebugLog.log("JSOraw.getSavedRawData(): data file count= %d " % (JSOrawSavedData.count))
    ch0_data=np.zeros([numTrigs,JSOrawSavedData.ch1_data_file.shape[1]])        
    ch1_data=np.zeros([numTrigs,JSOrawSavedData.ch1_data_file.shape[1]])        
   
    for i in range(numTrigs):
        if JSOrawSavedData.count==JSOrawSavedData.ch0_data_file.shape[0]:
            JSOrawSavedData.count=0           
        ch0_data[i,:]=JSOrawSavedData.ch0_data_file[JSOrawSavedData.count,:] 
        ch1_data[i,:]=JSOrawSavedData.ch1_data_file[JSOrawSavedData.count,:] 
        JSOrawSavedData.count=JSOrawSavedData.count+1             
    return ch0_data, ch1_data
    
def getNewRawData(numTrigs,requestedSamplesPerTrig,appObj):                                 
    try:
        print('requestedSamplesPerTrig= ', requestedSamplesPerTrig)
        err, ch0_data, ch1_data = appObj.oct_hw.AcquireOCTDataRaw(numTrigs, requestedSamplesPerTrig)
        print('ch0_data shape= ', repr(ch0_data.shape))
    except Exception as ex:
        print('Error collecting data')
        raise ex
        
    return ch0_data, ch1_data
      
def saveDispersion_pushButton_clicked(appObj):
    dispData=appObj.dispData
    dispData.sampleOffset = 2
    filename=datetime.datetime.now().strftime("dispComp-%Y_%m_%d-%H_%M_%S.pickle")
    outfile=os.path.join(appObj.configPath, 'Dispersion',filename)
    file1=open(outfile,'wb')
    pickle.dump(dispData,file1)
    file1.close()
    appObj.dispCompFilename_label.setText(outfile)
    # appObj.loadDispersionIntoFPGA(outfile, appObj.oct_hw.fpgaOpts)
    
def updateDispersionGUI(appObj, dispData):
        # Clear all of the plots
    appObj.mzi_plot_2.clear() 
    appObj.pd_plot_2.clear()
    appObj.mzi_mag_plot_2.clear()
    appObj.mzi_phase_plot_2.clear()
    appObj.k0_plot_2.clear()
    appObj.interp_pdRaw_plot.clear()
    appObj.interp_pdDispComp_plot.clear()
    appObj.alineNoInterp_plot.clear()
    appObj.alineRaw_plot.clear()
    appObj.alineDispComp_plot.clear()
    appObj.phaseNoiseTD_plot.clear()
    appObj.phaseNoiseFD_plot.clear()
    appObj.dispWnfcMag_plot.clear()
    appObj.dispWnfcPh_plot.clear()
    appObj.k0_plot_3.clear()
    appObj.k0_plot_4.clear()
    appObj.k0_plot_5.clear()
    
    # Plot the dispersion window functions and update the processing values
    appObj.dispWnfcMag_plot.plot(dispData.magWin, pen='b')
    appObj.dispWnfcPh_plot.plot(dispData.phaseCorr, pen='b')
    appObj.k0_plot_5.plot(dispData.k0Reference, pen='b') 
    
    appObj.requestedSamplesPerTrig.setValue(dispData.requestedSamplesPerTrig)
    appObj.startSample.setValue(dispData.startSample)
    appObj.endSample.setValue(dispData.endSample)
    appObj.numKlinPts.setValue(dispData.numKlinPts)
    appObj.numShiftPts.setValue(dispData.numShiftPts)
    appObj.filterWidth.setValue(dispData.filterWidth)
    appObj.mziFilter.setValue(dispData.mziFilter)
    appObj.dispMagWindowFilter.setValue(dispData.magWin_LPfilterCutoff)
    try:  # this field may not exist, so wrap around try/except
        appObj.dispersionCompAlgorithm_comboBox.setCurrentIndex(dispData.dispCode) 
    except:
        pass
    
def loadDispersion_pushButton_clicked(appObj):
    loadpath=os.path.join(appObj.configPath, 'Dispersion')  
    w = QtGui.QWidget()    
    infile = QtGui.QFileDialog.getOpenFileName(w, 'Open File', loadpath) 
    file2=open(infile,'rb')
    dispData=pickle.load(file2)
    file2.close()
    appObj.dispData=dispData
    appObj.dispCompFilename_label.setText(infile) 
    updateDispersionGUI(appObj, appObj.dispData)
    
#    appObj.loadDispersionIntoFPGA(infile, appObj.oct_hw.fpgaOpts)
    
def loadDispersion_onStartup(appObj):   
#    infile=os.path.join(appObj.configPath, 'Dispersion','dispComp-initial.pickle')  
    infile = appObj.octSetupInfo.dispFilename
    dispData = Dispersion.loadDispData(appObj, infile)
    appObj.dispData=dispData
    DebugLog.log('JSORaw.loadDispersion_onStartup: set appObj.dispData')
    appObj.dispCompFilename_label.setText(infile)
    updateDispersionGUI(appObj, appObj.dispData)

def softwareProcessing(ch0_data,ch1_data,zROI,appObj, returnPhase=False, unwrapPhase=False):
    # This routine can be called from any other routine to do software processing of the raw data. 
    # It needs appObj.dispData so you must have loaded a dispersion file already for this to work.

    dispData=appObj.dispData
    pdData,mziData,actualSamplesPerTrig=channelShift(ch0_data,ch1_data,dispData)    # shift the two channels to account for delays in the sample data compared to the MZI data 
        
    t1 = time.time()
    mzi_hilbert, mzi_mag, mzi_ph, k0 = processMZI(mziData, dispData)                # calculate k0 from the phase of the MZI data
    DebugLog.log("JSOraw.softwareProcessing(): process MZI time= %0.4f" % (time.time() - t1))
    DebugLog.log("JSOraw.softwareProcessing(): numTriggers collected= %d" % (k0.shape[0]))

    k0Cleaned=cleank0Run(k0,dispData) # Adjust the k0 curves so that the unwrapping all starts at the same phase    
    Klin=dispData.Klin
    
    t1 = time.time()
    pd_interpRaw, klin = processPD(pdData, k0Cleaned, dispData, Klin)  # Interpolate the PD data based upon the MZI data    
    DebugLog.log("JSOraw.softwareProcessing(): process PD time= %0.4f" % (time.time() - t1))
    
    t1 = time.time()    
    pd_interpDispComp = dispData.magWin * pd_interpRaw * (np.cos(-1*dispData.phaseCorr) + 1j * np.sin(-1*dispData.phaseCorr))  # perform dispersion compensation
    pd_fftDispComp, alineMagDispComp, alinePhaseDispComp = calculateAline(pd_interpDispComp, returnPhase, unwrapPhase) # calculate the a-line
    DebugLog.log("JSOraw.softwareProcessing(): calc Aline time= %0.4f" % (time.time() - t1))
    oct_data = pd_fftDispComp[:, zROI[0]:zROI[1]] 
    return oct_data, klin
    
def calibrateScanMirror(appObj):
    DebugLog.log("calibrateScanMirror")    
    appObj.tabWidget.setCurrentIndex(7)
    appObj.doneFlag = False
    appObj.isCollecting = True
    appObj.JSOsaveDispersion_pushButton.setEnabled(True)
    appObj.JSOloadDispersion_pushButton.setEnabled(False)

    if not appObj.oct_hw.IsOCTTestingMode():     # prepare to get new data            
        from DAQHardware import DAQHardware
        daq = DAQHardware()
    audioHW=appObj.audioHW
    mirrorDriver = appObj.mirrorDriver    
    chanNames = [mirrorDriver.X_daqChan, mirrorDriver.Y_daqChan]
    trigChan = audioHW.daqTrigChanIn   #use the audio trigger to start the scan
    outputRate = mirrorDriver.DAQoutputRate

    while appObj.doneFlag == False:        # keep running until the button is turned off 
        scanParams = appObj.getScanParams()
       #    create scan pattern to drive the mirrors
        mode=appObj.scanShape_comboBox.currentIndex()
        print('mode',mode)
        if mode==0:   # create a spiral scan using fast (resonant) scanning
            Vmaxx=mirrorDriver.voltRange[1] # maximum voltage for MEMS mirror for x-axis
            Vmaxy=mirrorDriver.voltRange[1] # maximum voltage for MEMS mirror for y-axis
            xAdjust = 1    
            yAdjust = scanParams.skewResonant
            phaseShift = scanParams.phaseAdjust
            fr = mirrorDriver.resonantFreq  # angular scan rate (frequency of one rotation - resonant frequency)
            fv = scanParams.volScanFreq     # plotParam scan frequency, which scans in and then out, which is actually two volumes
            DebugLog.log("freq of one rotation (fr)= %d; scan frequency (fv)= %d" % (fr, fv))
            diameter = scanParams.length
            voltsPerMM = mirrorDriver.voltsPerMillimeterResonant
            A1=(Vmaxx/2)/xAdjust
            A2=(Vmaxy/2)/yAdjust
            A3=voltsPerMM*diameter/2 
            A=np.min([A1,A2,A3])           
            fs=mirrorDriver.DAQoutputRate   # galvo output sampling rate
            t=np.arange(0,np.around(fs/fv))*1/fs  # t is the array of times for the DAQ output to the mirrors
            r=1/2*(1-np.cos(2*np.pi*fv*t))            
            x=xAdjust*A*r*np.cos(2*np.pi*fr*t) # x and y are the coordinates of the laser at each point in time
            y=yAdjust*A*r*np.sin(2*np.pi*fr*t+phaseShift*np.pi/180)
            mirrorOut= np.vstack((x,y))
            
        elif mode==1:   # create a square scan using slow parameters
            Vmaxx=mirrorDriver.voltRange[1] # maximum voltage for MEMS mirror for x-axis
            Vmaxy=mirrorDriver.voltRange[1] # maximum voltage for MEMS mirror for y-axis
            xAdjust = 1    
            yAdjust = scanParams.skewNonResonant
            diameter = scanParams.length
            voltsPerMMX = mirrorDriver.voltsPerMillimeter*xAdjust
            voltsPerMMY = mirrorDriver.voltsPerMillimeter*yAdjust
            if ((diameter/2)*voltsPerMMX)>Vmaxx:
                diameter=2*Vmaxx/voltsPerMMX
            if ((diameter/2)*voltsPerMMY)>Vmaxy:
                diameter=2*Vmaxy/voltsPerMMY
            freq = appObj.cal_freq_dblSpinBox.value()
            if freq>mirrorDriver.LPFcutoff:  # can't go faster than the maximum scan rate
                appObj.cal_freq_dblSpinBox.setValue(mirrorDriver.LPFcutoff)
            fs=mirrorDriver.DAQoutputRate   # galvo output sampling rate
            t1=np.arange(0,np.around(fs/freq))*1/fs  
            n=np.around(t1.shape[0]/4)-1   # number of points in each 4th of the cycle (reduce by 1 to make it easy to shorten t1)
            t=t1[0:4*n]  # t is the array of times for the DAQ output to the mirrors
            cornerX=(diameter/2)*voltsPerMMX     # voltage at each corner of the square            
            cornerY=(diameter/2)*voltsPerMMY     # voltage at each corner of the square            
            
            # x and y are the coordinates of the laser at each point in time
            x=np.zeros(t.shape)            
            y=np.zeros(t.shape)            
            x[0:n]=np.linspace(-cornerX,cornerX,n)
            y[0:n]=-cornerY
            x[n:2*n]=cornerX
            y[n:2*n]=np.linspace(-cornerY,cornerY,n)
            x[2*n:3*n]=np.linspace(cornerX,-cornerX,n)
            y[2*n:3*n]=cornerY
            x[3*n:4*n]=-cornerX
            y[3*n:4*n]=np.linspace(cornerY,-cornerY,n)
            mirrorOut1= np.vstack((x,y))
            if mirrorDriver.MEMS==True:
                mirrorOut=scipy.signal.filtfilt(mirrorDriver.b_filt,mirrorDriver.a_filt,mirrorOut1)           
            else:
                mirrorOut=mirrorOut1    

        # plot mirror commands to GUI 
        pl = appObj.JSOmisc_plot1
        npts = mirrorOut.shape[1]
        t = np.linspace(0, npts/outputRate, npts)
        pl.clear()
        pl.plot(t, mirrorOut[0, :], pen='b')  
        pl.plot(t, mirrorOut[1, :], pen='r')  
        labelStyle = appObj.xLblStyle
        pl.setLabel('bottom', 'Time', 's', **labelStyle)
        labelStyle = appObj.yLblStyle
        pl.setLabel('left', 'Output', 'V', **labelStyle)
    
        pl2=appObj.JSOmisc_plot2
        pl2.clear()
        pl2.plot(mirrorOut[0, :],mirrorOut[1, :], pen='b')
        labelStyle = appObj.xLblStyle
        pl2.setLabel('bottom', 'X galvo', 'V', **labelStyle)
        labelStyle = appObj.yLblStyle
        pl2.setLabel('left', 'Y galvo', 'V', **labelStyle)
     
        if not appObj.oct_hw.IsDAQTestingMode():
            # setup the analog output DAQ device
            daq.setupAnalogOutput(chanNames, trigChan, outputRate, mirrorOut.transpose())        
            daq.startAnalogOutput()
            
            #start trigger and wait for output to finish 
            daq.sendDigTrig(audioHW.daqTrigChanOut)
            daq.waitDoneOutput(timeout=3, stopAndClear=True)
            
            QtGui.QApplication.processEvents() # check for GUI events
        else:
            appObj.doneFlag = True      # just run one time through if in test mode
            appObj.CalibrateScanMirror_pushButton.setChecked(False)
                  
    # when testing is over, set the mirror position to (0,0)
    if not appObj.oct_hw.IsDAQTestingMode():
        chanNames = [mirrorDriver.X_daqChan, mirrorDriver.Y_daqChan]
        data = np.zeros(2)
        daq.writeValues(chanNames, data)
    
    appObj.JSOsaveDispersion_pushButton.setEnabled(False)    
    appObj.JSOloadDispersion_pushButton.setEnabled(True)        
    appObj.isCollecting = False
    appObj.finishCollection()
    
def runJSOraw(appObj):
    DebugLog.log("runJSOraw")
    try:
        appObj.tabWidget.setCurrentIndex(7)
        appObj.doneFlag = False
        appObj.isCollecting = True
        appObj.JSOsaveDispersion_pushButton.setEnabled(True)
        appObj.JSOloadDispersion_pushButton.setEnabled(False)
        dispData = appObj.dispData             # this class holds all the dispersion compensation data    
        if dispData is None:
            dispData = Dispersion.DispersionData()
            
        laserSweepFreq=appObj.octSetupInfo.getTriggerRate()
        mirrorDriver = appObj.mirrorDriver
        
        if not appObj.oct_hw.IsOCTTestingMode():     # prepare to get new data            
            # set the mirror position to (0,0)
            chanNames = [mirrorDriver.X_daqChan, mirrorDriver.Y_daqChan]
            data = np.zeros(2)
            from DAQHardware import DAQHardware
            daq = DAQHardware()
            daq.writeValues(chanNames, data)
        else:
            appObj.savedDataBuffer.loadData(appObj)
    
        peakXPos=np.array([0],dtype=int)       
        peakYPos=np.array([0],dtype=float)       
        peakXPos1=np.array([0],dtype=int)       
        peakYPos1=np.array([0],dtype=float)       
                   
        while appObj.doneFlag == False:
            # read data analysis settings from the GUI
            numTrigs=appObj.numTrig.value()
            dispData.requestedSamplesPerTrig=appObj.requestedSamplesPerTrig.value()
            dispData.startSample=appObj.startSample.value()
            dispData.endSample=appObj.endSample.value()
            dispData.numKlinPts=appObj.numKlinPts.value()
            dispData.Klin=np.zeros(dispData.numKlinPts)         
            dispData.numShiftPts=appObj.numShiftPts.value()
            dispData.filterWidth=appObj.filterWidth.value()
            dispData.mziFilter=appObj.mziFilter.value()
            dispData.magWin_LPfilterCutoff=appObj.dispMagWindowFilter.value()
            dispData.PDfilterCutoffs=[0,0]
            dispData.dispCode=appObj.dispersionCompAlgorithm_comboBox.currentIndex() 
            dispData.dispMode=appObj.dispersionCompAlgorithm_comboBox.currentText()
                     
            # Get data using one of several methods
            if appObj.oct_hw.IsOCTTestingMode():
                ch0_data,ch1_data=getSavedRawData(numTrigs,dispData.requestedSamplesPerTrig,appObj.savedDataBuffer)
            else:
                ch0_data,ch1_data=getNewRawData(numTrigs,dispData.requestedSamplesPerTrig,appObj)
            
            if appObj.saveData_checkBox.isChecked()==True:      # save data to disk for later use if desired
                fileName='Mirror_Raw'
                dataToSave = (ch0_data, ch1_data)              
                appObj.savedDataBuffer.saveData(appObj,dataToSave,fileName)                                
                appObj.saveData_checkBox.setChecked(False)                     
                
            # delay the MZI to account for it having a shorter optical path than the sample/reference arm path, then calculate k0 as the MZI phase
            pdData,mziData,actualSamplesPerTrig=channelShift(ch0_data,ch1_data,dispData)    
            textString='Actual samples per trigger: {actualSamplesPerTrig}'.format(actualSamplesPerTrig=actualSamplesPerTrig)            
            appObj.actualSamplesPerTrig_label.setText(textString)         
            
            import time
            t1 = time.time()
            mzi_hilbert, mzi_mag, mzi_ph, k0 = processMZI(mziData, dispData) 
            mzi_proc_time = time.time() - t1
            print("MZI processing time = %0.4f ms" % (mzi_proc_time*1000))
    
            # Adjust the k0 curves so that the unwrapping all starts at the same phase
            appObj.k0_plot_3.clear()
            appObj.k0_plot_4.clear()
            appObj.k0_plot_5.clear()
            t1 = time.time()
            k0Cleaned=cleank0(k0,dispData)              
            k0clean_time = time.time() - t1
            print("k0 cleaning time = %0.4f ms" % (k0clean_time*1000))
            
            for i in range(numTrigs):
                appObj.k0_plot_3.plot(k0[i,:2*dispData.startSample], pen=(i,numTrigs)) 
            startMZIdata1=k0[:,dispData.startSample]
            appObj.k0_plot_4.plot(startMZIdata1, pen='r') 
            startMZIdata2=k0Cleaned[:,dispData.startSample]
            appObj.k0_plot_4.plot(startMZIdata2, pen='b') 
            for i in range(numTrigs):
                appObj.k0_plot_5.plot(k0Cleaned[i,:2*dispData.startSample], pen=(i,numTrigs)) 
            k0=k0Cleaned
            
            # Interpolate the PD data based upon the MZI data and calculate the a-lines before dispersion compensation      
            t1 = time.time()
            pd_interpRaw, klin = processPD(pdData, k0, dispData)
            interpPD_time = time.time() - t1
            print("Interp PD time = %0.4f ms" % (interpPD_time*1000))
            dispData.Klin=klin
            pd_fftNoInterp, alineMagNoInterp, alinePhaseNoInterp = calculateAline(pdData[:,dispData.startSample:dispData.endSample])

            t1 = time.time()
            pd_fftRaw, alineMagRaw, alinePhaseRaw = calculateAline(pd_interpRaw)
            alineCalc_time = time.time() - t1
            print("Aline calc time = %0.4f ms" % (alineCalc_time*1000))
            
            # find the mirror in the a-line to determine the filter settings, and then perform the dispersion compensatsion 
            rangePeak1=[100, 900]
            alineAve1=np.average(alineMagRaw,axis=0) 
            peakXPos1[0]=np.argmax(alineAve1[rangePeak1[0]:rangePeak1[1]])+rangePeak1[0]
            peakYPos1[0]=alineAve1[peakXPos1[0]]     
            width=dispData.filterWidth*(rangePeak1[1]-rangePeak1[0])/2         
            dispData.PDfilterCutoffs[0]=(peakXPos1[0]+width)/2048
            dispData.PDfilterCutoffs[1]=(peakXPos1[0]-width)/2048
        
            dispersionCorrection(pd_interpRaw,dispData)
            appObj.dispData=dispData      #store the local variable in the overall class so that it can be saved when the save button is pressed
            
            # now correct the data using dispersion compensation and then process the a-lines
            pd_interpDispComp = dispData.magWin * pd_interpRaw * (np.cos(-1*dispData.phaseCorr) + 1j * np.sin(-1*dispData.phaseCorr))
            pd_fftDispComp, alineMagDispComp, alinePhaseDispComp = calculateAline(pd_interpDispComp)
                  
            #scale k0 and the MZI to the same range to plot them so they overlap
            k0Ripple= scipy.signal.detrend(k0[0,500:700],axis=-1)
            k0RippleNorm=k0Ripple/k0Ripple.max()
            mziDataRipple= scipy.signal.detrend(mziData[0,500:700],axis=-1)
            mziDataNorm=mziDataRipple/mziDataRipple.max()
           
            # Find the peak of the A-line within a range and calculate the phase noise
            rangePeak=[100, 900]
            alineAve=np.average(alineMagDispComp,axis=0) 
            peakXPos[0]=np.argmax(alineAve[rangePeak[0]:rangePeak[1]])+rangePeak[0]
            peakYPos[0]=alineAve[peakXPos[0]]              
            t=np.arange(numTrigs)/laserSweepFreq
            phaseNoiseTD=np.unwrap(alinePhaseDispComp[:,peakXPos[0]])
            phaseNoiseTD=phaseNoiseTD-np.mean(phaseNoiseTD)
            phaseNoiseTD=phaseNoiseTD*1310e-9/(4*np.pi*1.32)
            phaseNoiseFFT = np.abs(np.fft.rfft(phaseNoiseTD))/(numTrigs/2)
#            phaseNoiseFD = 20*np.log10(np.abs(phaseNoiseFFT))        
            freq = np.fft.rfftfreq(numTrigs)*laserSweepFreq
    
            # Clear all of the plots
            appObj.mzi_plot_2.clear() 
            appObj.pd_plot_2.clear()
            appObj.mzi_mag_plot_2.clear()
            appObj.mzi_phase_plot_2.clear()
            appObj.k0_plot_2.clear()
            appObj.interp_pdRaw_plot.clear()
            appObj.interp_pdDispComp_plot.clear()
            appObj.alineNoInterp_plot.clear()
            appObj.alineRaw_plot.clear()
            appObj.alineDispComp_plot.clear()
            appObj.phaseNoiseTD_plot.clear()
            appObj.phaseNoiseFD_plot.clear()
            appObj.dispWnfcMag_plot.clear()
            appObj.dispWnfcPh_plot.clear()
           
            # Plot all the data
            if appObj.plotFirstOnly_checkBox.isChecked()==True:
                i=0
                appObj.pd_plot_2.plot(pdData[i,:], pen='r')            
                appObj.mzi_plot_2.plot(mziData[i,:], pen='r')            
                appObj.mzi_mag_plot_2.plot(mzi_mag[i,:], pen='r')            
                appObj.k0_plot_2.plot(k0[i,:], pen='r')
                sampleNum=np.linspace(dispData.startSample,dispData.endSample,dispData.numKlinPts)
                appObj.k0_plot_2.plot(sampleNum,klin, pen='b')                      
                appObj.interp_pdRaw_plot.plot(pd_interpRaw[i,:], pen='r')           
                appObj.interp_pdDispComp_plot.plot(np.abs(pd_interpDispComp[i,:]), pen='r')           
                appObj.alineNoInterp_plot.plot(alineMagNoInterp[i,:], pen='r')
                appObj.alineRaw_plot.plot(alineMagRaw[i,:], pen='r')
                appObj.alineDispComp_plot.plot(alineMagDispComp[i,:], pen='r')
            else:
                # limit plotting to first 10 or so triggers, otherwise this will freeze up
                nTrigs = min((numTrigs, 10))
                
                for i in range(nTrigs):
                    pen=(i,nTrigs)
                    appObj.pd_plot_2.plot(pdData[i,:], pen=pen)            
                    appObj.mzi_plot_2.plot(mziData[i,:], pen=pen)            
                    appObj.mzi_mag_plot_2.plot(mzi_mag[i,:], pen=pen)            
                    appObj.mzi_phase_plot_2.plot(mzi_ph[i,:], pen=pen)            
                    appObj.k0_plot_2.plot(k0[i,:], pen=pen)            
                    appObj.interp_pdRaw_plot.plot(pd_interpRaw[i,:], pen=pen)            
                    appObj.interp_pdDispComp_plot.plot(np.abs(pd_interpDispComp[i,:]), pen=pen)            
                    appObj.alineNoInterp_plot.plot(alineMagNoInterp[i,:], pen=pen)            
                    appObj.alineRaw_plot.plot(alineMagRaw[i,:], pen=pen)            
                    appObj.alineDispComp_plot.plot(alineMagDispComp[i,:], pen=pen)            
                
            appObj.alineRaw_plot.plot(peakXPos1,peakYPos1, pen=None, symbolBrush='k', symbolPen='b')
            appObj.alineDispComp_plot.plot(peakXPos,peakYPos, pen=None, symbolBrush='k', symbolPen='b')
            appObj.phaseNoiseTD_plot.plot(t,phaseNoiseTD, pen='r')
            appObj.phaseNoiseFD_plot.plot(freq,phaseNoiseFFT, pen='r')
            appObj.mzi_phase_plot_2.plot(mziDataNorm, pen='b')            
            appObj.mzi_phase_plot_2.plot(k0RippleNorm, pen='r')            
            
           
            # if you want to align the pd and the Mzi data
    #            plotPDPhase.plot(pdData[0,:], pen='r')
    #            plotPDPhase.plot(mziData[0,:], pen='b')
    
            appObj.dispWnfcMag_plot.plot(dispData.magWin, pen='b')
            appObj.dispWnfcPh_plot.plot(dispData.phaseCorr, pen='b')
            
            # plot filter cutoff ranges on the raw Aline plot
            yy=[np.min(alineMagRaw[0,:]),np.max(alineMagRaw[0,:])]
            xx0=[alineMagRaw.shape[1]*dispData.PDfilterCutoffs[0],alineMagRaw.shape[1]*dispData.PDfilterCutoffs[0]]        
            xx1=[alineMagRaw.shape[1]*dispData.PDfilterCutoffs[1],alineMagRaw.shape[1]*dispData.PDfilterCutoffs[1]]    
            appObj.alineRaw_plot.plot(xx0,yy, pen='b')
            appObj.alineRaw_plot.plot(xx1,yy, pen='b')
            
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
            if ~np.all(np.isnan(alineMagDispComp)):  # only make the bscan plot if there is data to show (this prevents an error from occurring)
                appObj.bscan_plot.setImage(alineMagDispComp)
            QtGui.QApplication.processEvents() # check for GUI events  
    except:
        traceback.print_exc(file=sys.stdout)
        QtGui.QMessageBox.critical (appObj, "Error", "Error during scan. Check command line output for details")
        
    appObj.JSOsaveDispersion_pushButton.setEnabled(False)    
    appObj.JSOloadDispersion_pushButton.setEnabled(True)
    appObj.isCollecting = False
    QtGui.QApplication.processEvents() # check for GUI events
    appObj.finishCollection()
