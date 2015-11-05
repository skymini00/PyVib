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
import numpy as np
import datetime
import matplotlib.pyplot as plt
from DebugLog import DebugLog
import pickle

class blankClass:
    def __init__(self):
        pass
    
class DispersionData:  # this class holds all the dispersion compensation data
    def __init__(self):
        self.magWin = []
        self.phaseCorr = []
        self.phDiode_background = None
        self.startSample=[]
        self.endSample=[]
        self.numKlinPts=[]
        self.numShiftPts=[]
        self.filterWidth=[]
        self.PDfilterCutoffs=[]
        self.mziFilter=[]
        self.magWin_LPfilterCutoff=[]
            
def processMZI(mzi_data, dispData):
    # filtering seems to reduce sidebands created during the interpolation process
    (b, a) = scipy.signal.butter(2, dispData.mziFilter, 'highpass')
    mzi_data=scipy.signal.lfilter(b, a, mzi_data,axis=-1)    
    
    mzi_hilbert = scipy.signal.hilbert(mzi_data, axis=-1)
    mzi_mag = np.abs(mzi_hilbert)
    mzi_ph = np.angle(mzi_hilbert)
    mzi_hilbert = np.imag(mzi_hilbert)
    k0 = np.unwrap(mzi_ph,axis=-1)    
    return mzi_hilbert, mzi_mag, mzi_ph, k0
 
def cleank0(k0,dispData):
    k0Cleaned=k0    
    k0Init=k0Cleaned[:,dispData.startSample]
    while (np.max(k0Init)-np.min(k0Init))>(1.5*np.pi):
        # there are phase jumps going on in the data
        # need to put code in here to unwrap, find # points that are different, and shift to the most common point. Then re-run algorithm
        break
        
    return k0Cleaned

   
def processPD(pd_data, k0, dispData, klin=None):
    if klin is None:
        klin = np.linspace(k0[0,dispData.startSample], k0[0,dispData.endSample], dispData.numKlinPts)
    pd_interp=np.zeros([pd_data.shape[0],dispData.numKlinPts])
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
    elif dispData.dispMode=='Unique':
        processUniqueDispersion(pd,dispData)
    else:
        dispData.magWin = np.ones(dispData.numKlinPts)
        dispData.phaseCorr = np.zeros(dispData.numKlinPts)
    
def calculateAline(pd):
    numPts=pd.shape[1] 
    pd_fft = np.fft.fft(pd, n=2048,axis=-1)/numPts
    alineMag = 20*np.log10(np.abs(pd_fft) + 1)
    alinePhase=np.unwrap(np.angle(pd_fft),axis=-1)
    return pd_fft, alineMag, alinePhase

def channelShift(ch0_data,ch1_data,dispData):
    actualSamplesPerTrig=ch0_data.shape[1]-dispData.numShiftPts
    pdData=ch0_data[:,dispData.numShiftPts:]
    mziData=ch1_data[:,0:actualSamplesPerTrig]
    return pdData,mziData,actualSamplesPerTrig

def calcOCTDataFFT(pdData, mziData, MZI_PD_shift, klinROI_idx, numklinpts, klin, dispCorr_mag, dispCorr_ph, zROI):
    DebugLog.log('calcOCTDataFFT: pdData.shape= ' + repr(pdData.shape))
    pdData,mziData,actualSamplesPerTrig=channelShift(pdData, mziData, MZI_PD_shift)    
    DebugLog.log('calcOCTDataFFT: after channel Shift pdData.shape= ' + repr(pdData.shape))
    mzi_hilbert, mzi_mag, mzi_ph, k0 = processMZI(mziData) 
    DebugLog.log('calcOCTDataFFT: processed the MZI data pdData.shape= ' + repr(pdData.shape))
    pd_interpRaw, klin = processPD(pdData, k0, klinROI_idx, numklinpts)   
    if numklinpts > 2048:  # downsample if over 2048 pts
        idx = range(0, numklinpts, 2)
        pd_interpRaw = pd_interpRaw[:, idx]
    print('processed the PDdata')
    pd_interpDisp = np.abs(dispCorr_mag * pd_interpRaw * (np.cos(-1*dispCorr_ph) + 1j * np.sin(-1*dispCorr_ph)))

    oct_data = np.fft.fft(pd_interpDisp, 2048)
    oct_data = oct_data[:, zROI[0]:zROI[1]]
    oct_data = oct_data / numklinpts   
    return oct_data, klin
    
def processUniqueDispersion(pd_data, dispData, pd_background=None, bg_collect=False):
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
    print('getting saved raw data',JSOrawSavedData.count)
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
        err, ch0_data, ch1_data = appObj.oct_hw.AcquireOCTDataRaw(numTrigs, requestedSamplesPerTrig)
    except Exception as ex:
        print('Error collecting data')
        raise ex
    
    if appObj.saveData_checkBox.isChecked()==True:      # get new data and save to disk for later use
        outfile='testData.npz'
        np.savez_compressed(outfile, ch0_data=ch0_data, ch1_data=ch1_data)
        print('saved data in :',outfile)
        appObj.saveData_checkBox.setChecked(False)                     
    return ch0_data, ch1_data
      
def saveDispersion_pushButton_clicked(appObj):
    dispData=appObj.dispData
    filename=datetime.datetime.now().strftime("dispComp-%Y_%m_%d-%H_%M_%S.pickle")
    outfile=os.path.join(appObj.configPath, 'Dispersion',filename)
    print(appObj.dispData.numKlinPts,appObj.dispData.magWin.shape)
    file1=open(outfile,'wb')
    pickle.dump(dispData,file1)
    file1.close()
    appObj.dispCompFilename_label.setText(outfile)
    
def loadDispersion_pushButton_clicked(appObj):
    loadpath=os.path.join(appObj.configPath, 'Dispersion')  
    w = QtGui.QWidget()    
    infile = QtGui.QFileDialog.getOpenFileName(w, 'Open File', loadpath)
    
    file2=open(infile,'rb')
    dispData=pickle.load(file2)
    file2.close()

    appObj.dispData=dispData
    appObj.dispCompFilename_label.setText(infile) 

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
    appObj.startSample.setValue(dispData.startSample)
    appObj.endSample.setValue(dispData.endSample)
    appObj.numKlinPts.setValue(dispData.numKlinPts)
    appObj.numShiftPts.setValue(dispData.numShiftPts)

def runJSOraw(appObj):
    DebugLog.log("runJSOraw")
    appObj.tabWidget.setCurrentIndex(7)
    appObj.doneFlag = False
    appObj.isCollecting = True
    testDataDir = os.path.join(appObj.basePath, 'exampledata', 'JSOraw')
    
    dispData = DispersionData()             # this class holds all the dispersion compensation data    
    laserSweepFreq=appObj.oct_hw.GetTriggerRate()
    mirrorDriver = appObj.mirrorDriver
    
    if not appObj.oct_hw.IsOCTTestingMode():     # prepare to get new data            
        # set the mirror position to (0,0)
        chanNames = [mirrorDriver.X_daqChan, mirrorDriver.Y_daqChan]
        data = np.zeros(2)
        from DAQHardware import DAQHardware
        daq = DAQHardware()
        daq.writeValues(chanNames, data)
    else:    # load in the saved data and use it instead
        outfile=os.path.join(testDataDir,'testData.npz')
        x=np.load(outfile)             
        JSOrawSavedData=blankClass()        #This class stores the data from the disk file
        JSOrawSavedData.ch0_data_file=x['ch0_data']
        JSOrawSavedData.ch1_data_file=x['ch1_data']
        JSOrawSavedData.count=0
        JSOrawSavedData.saveRawData=0
        appObj.requestedSamplesPerTrig.setValue(JSOrawSavedData.ch1_data_file.shape[1])
        print('loaded data from :',outfile)
           
    peakXPos=np.array([0],dtype=int)       
    peakYPos=np.array([0],dtype=float)       
    peakXPos1=np.array([0],dtype=int)       
    peakYPos1=np.array([0],dtype=float)       
               
    while appObj.doneFlag == False:
        # read data analysis settings from the GUI
        numTrigs=appObj.numTrig.value()
        requestedSamplesPerTrig=appObj.requestedSamplesPerTrig.value()
        dispData.startSample=appObj.startSample.value()
        dispData.endSample=appObj.endSample.value()
        dispData.numKlinPts=appObj.numKlinPts.value()
        dispData.numShiftPts=appObj.numShiftPts.value()
        dispData.filterWidth=appObj.filterWidth.value()
        dispData.mziFilter=appObj.mziFilter.value()
        dispData.magWin_LPfilterCutoff=appObj.dispMagWindowFilter.value()
        dispData.PDfilterCutoffs=[0,0]
        dispCode=appObj.dispersionCompAlgorithm_comboBox.currentIndex()            
        if dispCode==0:
            dispData.dispMode='None'
        elif dispCode==1:
            dispData.dispMode='Hanning'
        elif dispCode==2:
            dispData.dispMode='Hamming'
        elif dispCode==3:
            dispData.dispMode='Unique'
        else:
            dispData.dispMode='None'
                 
        # Get data using one of several methods
        if appObj.oct_hw.IsOCTTestingMode():
            ch0_data,ch1_data=getSavedRawData(numTrigs,requestedSamplesPerTrig,JSOrawSavedData)
        else:
            ch0_data,ch1_data=getNewRawData(numTrigs,requestedSamplesPerTrig,appObj)
        
        # delay the MZI to account for it having a shorter optical path than the sample/reference arm path, then calculate k0 as the MZI phase
        pdData,mziData,actualSamplesPerTrig=channelShift(ch0_data,ch1_data,dispData)    
        textString='Actual samples per trigger: {actualSamplesPerTrig}'.format(actualSamplesPerTrig=actualSamplesPerTrig)            
        appObj.actualSamplesPerTrig_label.setText(textString)         
        mzi_hilbert, mzi_mag, mzi_ph, k0 = processMZI(mziData, dispData) 

        # k-means cluster the k0 curves to figure out how to start the unwrapping
        k0Cleaned=cleank0(k0,dispData)              

        # Interpolate the PD data based upon the MZI data and calculate the a-lines before dispersion compensation      
        pd_interpRaw, klin = processPD(pdData, k0, dispData)
        pd_fftNoInterp, alineMagNoInterp, alinePhaseNoInterp = calculateAline(pdData[:,dispData.startSample:dispData.endSample])
        pd_fftRaw, alineMagRaw, alinePhaseRaw = calculateAline(pd_interpRaw)
        
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
        phaseNoiseTD=np.unwrap(alinePhaseDispComp[:,peakXPos[0]])
        phaseNoiseTD=phaseNoiseTD-phaseNoiseTD[0]
        phaseNoiseFFT = np.fft.fft(phaseNoiseTD, n=2048)
        phaseNoiseFD = 20*np.log10(np.abs(phaseNoiseFFT) + 1)        
        
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
            for i in range(numTrigs):
                appObj.pd_plot_2.plot(pdData[i,:], pen=(i,numTrigs))            
                appObj.mzi_plot_2.plot(mziData[i,:], pen=(i,numTrigs))            
                appObj.mzi_mag_plot_2.plot(mzi_mag[i,:], pen=(i,numTrigs))            
                appObj.mzi_phase_plot_2.plot(mzi_ph[i,:], pen=(i,numTrigs))            
                appObj.k0_plot_2.plot(k0[i,:], pen=(i,numTrigs))                      
                appObj.interp_pdRaw_plot.plot(pd_interpRaw[i,:], pen=(i,numTrigs))           
                appObj.interp_pdDispComp_plot.plot(np.abs(pd_interpDispComp[i,:]), pen=(i,numTrigs))           
                appObj.alineNoInterp_plot.plot(alineMagNoInterp[i,:], pen=(i,numTrigs))
                appObj.alineRaw_plot.plot(alineMagRaw[i,:], pen=(i,numTrigs))
                appObj.alineDispComp_plot.plot(alineMagDispComp[i,:], pen=(i,numTrigs))
            
        appObj.alineRaw_plot.plot(peakXPos1,peakYPos1, pen=None, symbolBrush='k', symbolPen='b')
        appObj.alineDispComp_plot.plot(peakXPos,peakYPos, pen=None, symbolBrush='k', symbolPen='b')
        appObj.phaseNoiseTD_plot.plot(phaseNoiseTD, pen='r')
        appObj.phaseNoiseFD_plot.plot(phaseNoiseFD, pen='r')
        appObj.mzi_phase_plot_2.plot(mziDataNorm, pen='b')            
        appObj.mzi_phase_plot_2.plot(k0RippleNorm, pen='r')            
        
        for i in range(numTrigs):
            appObj.k0_plot_3.plot(k0[i,:25], pen=(i,numTrigs)) 
        startMZIdata1=k0[:,dispData.startSample]
        appObj.k0_plot_4.plot(startMZIdata1, pen='r') 
        startMZIdata2=k0Cleaned[:,dispData.startSample]
        appObj.k0_plot_4.plot(startMZIdata2, pen='b') 
        
        # if you want to align the pd and the Mzi data
#            plotPDPhase.plot(pdData[0,:], pen='r')
#            plotPDPhase.plot(mziData[0,:], pen='b')

        appObj.dispWnfcMag_plot.plot(dispData.magWin, pen='b')
        appObj.dispWnfcPh_plot.plot(dispData.phaseCorr, pen='b')
        
        # plot filter cutoff ranges on the raw Aline plot
        yy=[0,np.max(alineMagRaw[0,:])]
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
        appObj.bscan_plot.setImage(alineMagDispComp)
        QtGui.QApplication.processEvents() # check for GUI events  
        
    appObj.isCollecting = False
    QtGui.QApplication.processEvents() # check for GUI events
    appObj.finishCollection()
