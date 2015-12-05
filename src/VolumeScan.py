# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 12:22:44 2015

@author: OHNS
"""

import OCTCommon
import OCTFPGAProcessingInterface as octfpga
import BScan
import JSOraw

from DebugLog import DebugLog
import scipy

from OCTProtocolParams import *
from PyQt4 import QtCore, QtGui, uic
from ROIImageGraphicsView import *

import PIL
import copy
import numpy as np
import traceback
import sys
import os
import tifffile
import pickle
import time
import multiprocessing as mproc
import queue

class blankRecord():
    pass

# processing optiosn for the volume scan
class ProcOpts:
    def __init__(self):
        self.normLow = 20
        self.normHigh = 100
        self.spiralScanThresholdVal = 0.3
        self.biDirTrigVolAdj = 0
        self.zRes = 8.3e-6

# processed data of a single stim of single pt
class VolumeData:
    def __init__(self):
        self.scanParams = None
        self.volumeImg = None                 # should be a numpy.ndarray (of 3 dimensions)
        self.zPixSize = np.NaN     # z pixel size in um
        self.xPixSize = np.NaN     # x pixel size in um
        self.yPixSize = np.NaN     # y pixel size in um 
        self.scanNum = 0
        self.spiralScanData = None

class SpiralScanData():
    def __init__(self):
        self.bothPlots = None
        self.surfacePlot = None
        self.bscanPlot = None
        self.bscanPlot_16b = None

# raw data used for multi processing
class VolumeRawData:
    def __init__(self, frameNum, packedData=None, oct_data=None, isPackedMagOnly=True):
        self.packedData = packedData
        self.isPackedMagOnly = isPackedMagOnly
        self.oct_data = oct_data
        self.frameNum = frameNum

        
def setupScan(scanParams, mirrorDriver, zROI, OCTtrigRate, procOpts):
    """
    This function should be called to setup the output and reconstruction parameters prior to starting the protocol,
    so that the getMirrorOutput(), processData() functions will return correct values
    """
    DebugLog.log("VolumeScan.setupScan() zROI=%s OCTtrigRate=%d " % (repr(zROI), OCTtrigRate))
    # define plotting parameters and then setup plotting algorithms
    plotParam=blankRecord() # create an empty record to pass the plotting parameters in
    plotParam.zROI = zROI
    plotParam.zPixel=plotParam.zROI[1]-plotParam.zROI[0]+1              
    plotParam.zPixelsize=procOpts.zRes    
    plotParam.xPixel=scanParams.lengthSteps   # number of pixels in x dimension  
    plotParam.yPixel=scanParams.widthSteps   # number of pixels in y dimension        
    plotParam.xCenter=np.int(plotParam.xPixel/2)
    plotParam.yCenter=np.int(plotParam.yPixel/2)            
    plotParam.rangeCenter=((plotParam.xCenter+plotParam.yCenter)//2)//4   
    
    DebugLog.log("VolumeScan.setupScan() plotParam xPixel= %d yPixel= %d zPixel= %d xCenter= %d yCenter= %d rangeCenter= %d" % (plotParam.xPixel, plotParam.yPixel, plotParam.zPixel, plotParam.xCenter, plotParam.yCenter, plotParam.rangeCenter))
    #  setup galvo voltages and calculate which samples fit into which pixels                 
    scanDetails=blankRecord()  #create an empty record to get the spiral scan parameter
    calcuateAngledBScans(plotParam)
    if scanParams.pattern == ScanPattern.spiral:
        setupSpiralScan(scanParams, mirrorDriver, scanDetails, plotParam, OCTtrigRate)
    elif scanParams.pattern == ScanPattern.wagonWheel:
        setupWagonWheelScan(scanParams, mirrorDriver, scanDetails, plotParam, OCTtrigRate)
    elif scanParams.pattern == ScanPattern.zigZag:
        setupZigZagScan(scanParams, mirrorDriver, scanDetails, plotParam, OCTtrigRate)

    return plotParam, scanDetails
    
# TODO: need to implement this function
def setupWagonWheelScan(scanParams, mirrorDriver, scanDetails,plotParam, OCTtrigRate):
    pass

def setupZigZagScan(scanParams, mirrorDriver, scanDetails, plotParam, OCTtrigRate):
        """
        This function builds the zigzag scan voltages (scanDetails.mirrOut) 
        and the variables that are used when reformatting the collected A-lines
        into the proper 3D format(scanDetails.c,scanDetails.cTile, and scanDetails.c3)
        """
        Vmaxx=mirrorDriver.voltRange[1] #maximum voltage for MEMS mirror for x-axis
        Vmaxy=mirrorDriver.voltRange[1] #minimum voltage for MEMS mirror for y-axis
        MaxVinx = (Vmaxx / 7.5) - 9.333 #maximum input voltage to amplifier assuming
        MaxViny = (Vmaxy / 7.5) - 9.333 #maximum input voltage
        LPFcutoff=mirrorDriver.LPFcutoff  #recommended low pass cutoff frequency
        RFx=scanParams.angularScanFreq #resonant frequency used for x-axis 
        RFy=scanParams.volScanFreq # frequency to collect a two complete volumes along the y-axis 

        # System parameters
        plotParam.voltsPerMM=mirrorDriver.voltsPerMillimeter
        angle=0 #angle in degrees for the Bscan 
        fs= OCTtrigRate  #laser sweep rate
        fDAQ= self.mirrorDriver.DAQoutputRate  #DAQ sampling rate       
        FOV = 10 #field of view in mm
        plotParam.xPixelsize=(FOV/plotParam.xPixel)*1000  # size on one pixel in the x dimension in microns
        plotParam.yPixelsize=(FOV/plotParam.yPixel)*1000  # size on one pixel in the y dimension in microns         
        spatialsampling= 5  #spatial sampling in microns       
        Vpmmx=1*plotParam.voltsPerMM #volts per millimeter for x-axis
        Vpmmy=1*plotParam.voltsPerMM #volts per millimeter for y-axis
        if FOV>2*np.min([MaxVinx/Vpmmx,MaxViny/Vpmmy]):
            FOV=2*np.min([MaxVinx/Vpmmx,MaxViny/Vpmmy])
        n=round(FOV/(spatialsampling/1000))  #number of samples in both x and y
        # Create filtered sawtooth scan on the x-axis      
        faxis= fs/n #fast axis frequency
        cycles=51  #number of cycles to use in to generate filtered waveform, we will use the middle cycle to avoid artifacts at the edges
        ts=np.arange(0, 1/faxis*cycles, 1/fs) #generate time array for laser sweep clock
        tDAQ=np.arange(0, 1/faxis*cycles, 1/fDAQ) #generate time array for DAQ clock
        xraw=signal.sawtooth(2*np.pi*faxis*ts,0.5) #raw sawtooth signal for laser sweep
        xDAQraw=signal.sawtooth(2*np.pi*faxis*tDAQ,0.5) #raw sawtooth signal for DAQ
        
        b,a=signal.butter(3,LPFcutoff/fs/2,'low')  #note since using filtfilt a 3 pole is actually a 6 pole
        xfil=signal.filtfilt(b,a,xraw) #filter to avoid MEMS resonance       
        b,a=signal.butter(3,LPFcutoff/fDAQ/2,'low')  #note since using filtfilt a 3 pole is actually a 6 pole
        xDAQfil=signal.filtfilt(b,a,xDAQraw) #filter to avoid MEMS resonance         
        tlower=np.mean(tDAQ)-1/faxis/2 #lower time value for center cycle
        thigher=np.mean(tDAQ)+1/faxis/2 #higher time value for center cycle
        tcycleindexs=np.where(((ts>=tlower)*(ts<=thigher))>0) #find indicies coresponding to center cycle
        tcycleindexDAQ=np.where(((tDAQ>=tlower)*(tDAQ<=thigher))>0) #find indicies coresponding to center cycle
            
        xfilcycle=xfil[tcycleindexs[0][0]:tcycleindexs[0][-1]] /np.max(xfil[tcycleindexs[0][0]:tcycleindexs[0][-1]]) #filtered cycle of sawtooth for generating scan waveform
        xDAQfilcycle=xDAQfil[tcycleindexDAQ[0][0]:tcycleindexDAQ[0][-1]]/np.max(xDAQfil[tcycleindexDAQ[0][1]:tcycleindexDAQ[0][-1]]) #filtered cycle of sawtooth for generating scan waveform
    
        tqcycleindexs=np.where(((ts>=(tlower+1/faxis/4))*(ts<=(thigher-1/faxis/2)))>0) #find indicies for quarter cycle
        tqcycleindexDAQ=np.where(((tDAQ>=(tlower+1/faxis/4))*(tDAQ<=(thigher-1/faxis/2)))>0) #find indicies for quarter cycle
    
        xqfilcycle=xfil[tqcycleindexs[0][0]:tqcycleindexs[0][-1]]/np.max(xfil[tqcycleindexs[0][1]:tqcycleindexs[0][-1]]) #filtered quarter cycle to initiate scant at 0 V
        xqDAQfilcycle=xDAQfil[tqcycleindexDAQ[0][0]:tqcycleindexDAQ[0][-1]]/np.max(xDAQfil[tqcycleindexDAQ[0][0]:tqcycleindexDAQ[0][-1]]) #filtered quarter cycle to initiate scant at 0 V        
    
        A=-xqfilcycle #quarter cycle to begin scan at 0 V
        B=-np.flipud(xqfilcycle) #quarter cycle to end scan at 0 V
    
        ADAQ=-xqDAQfilcycle #quarter cycle to begin scan at 0 V
        BDAQ=-np.flipud(xqDAQfilcycle) #quarter cycle to end scan at 0 V
    
        xN=np.concatenate((A,np.tile(xfilcycle,n),B),2) #build normalized, from -1 to 1, x-waveform
        yN=np.concatenate((A,np.linspace(-1,1,xfilcycle.size*n),-B),2) #build linear ramp for y-waveform
    
        xDAQ=np.concatenate((ADAQ,np.tile(xDAQfilcycle,n),BDAQ),2) #build normalized, from -1 to 1, x-waveform
        yDAQ=np.concatenate((ADAQ,np.linspace(-1,1,xDAQfilcycle.size*n),-BDAQ),2) #build linear ramp for y-waveform
         
        Vx=FOV*Vpmmx/2
        x=Vx*xDAQ        
        Vy=FOV*Vpmmy/2
        y=Vy*yDAQ
        scanDetails.mirrOut=[x,y]
        self.plot_mirrorCmd.clear()
        self.plot_mirrorCmd.plot(x,y)
       
        # calculate the x,y pixel positions for each sampletime
        xNorm1=plotParam.xPixel*((xN/2)+0.5)  
        yNorm1=plotParam.yPixel*((yN/2)+0.5)
        scanDetails.numTrigs=xN.size
        xNorm1[0:int(xN.size/4-1)]=1e10
        xNorm1[int(xN.size/4*3):-1]=1e10
        yNorm1[0:int(yN.size/4-1)]=1e10
        yNorm1[int(yN.size/4*3):-1]=1e10
        
        # create the 3D array c2, of true/false to give which alines to average for each pixel, and then reformat this into scanDetails.c
        scanDetails.c3=np.zeros((plotParam.xPixel,plotParam.yPixel))
        for i1 in range(0,plotParam.xPixel):  
            for i2 in range(0,plotParam.yPixel):
                distance=np.sqrt((xNorm1-i1)**2+(yNorm1-i2)**2)
                if np.min(distance)<np.sqrt(2):
                    scanDetails.c3[i1,i2]=np.argmin(distance)                                   
        scanDetails.c3=np.uint(np.reshape(scanDetails.c3,(plotParam.xPixel*plotParam.yPixel)))


def setupSpiralScan(scanParams, mirrorDriver, scanDetails, plotParam, OCTtrigRate):
    """ 
    This function builds the spiral scan voltages (scanDetails.mirrOut) 
    and the variables that are used when reformatting the collected A-lines
    into the proper 3D format(scanDetails.c,scanDetails.cTile, and scanDetails.c3)
    """
    # make mirror output signals
    xAdjust = 1    
    yAdjust = scanParams.skew
    phaseShift = scanParams.phaseAdjust
    fr = scanParams.angularScanFreq  # angular scan rate (frequency of one rotation)
    fv = scanParams.volScanFreq     # plotParam scan frequency, which scans in and then out, which is actually two volumes
    DebugLog.log("VolumeScan.setupSpiralScan(): freq of one rotation (fr)= %d; scan frequency (fv)= %d" % (fr, fv))
    diameter = scanParams.length
    plotParam.xPixelsize=(diameter/plotParam.xPixel)*1000  # size on one pixel in the x dimension in microns
    plotParam.yPixelsize=(diameter/plotParam.yPixel)*1000  # size on one pixel in the y dimension in microns         
    voltsPerMM = mirrorDriver.voltsPerMillimeter
    A=voltsPerMM*diameter/2     
    fs=mirrorDriver.DAQoutputRate   # galvo output sampling rate
    t=np.arange(0,np.around(fs/fv))*1/fs  # t is the array of times for the DAQ output to the mirrors
    r=1/2*(1-np.cos(2*np.pi*fv*t))            
    x=xAdjust*A*r*np.cos(2*np.pi*fr*t) # x and y are the coordinates of the laser at each point in time
    y=yAdjust*A*r*np.sin(2*np.pi*fr*t+phaseShift*np.pi/180)
    scanDetails.mirrOut= np.vstack((x,y))
    DebugLog.log("VolumeScan.setupSpiralScan(): Number of time points to send to galvos for one volume = %d" % (len(t)))
  
    # reconstruct the image from the array of A-lines
    scanDetails.numTrigs=np.int32(OCTtrigRate/fv)     # calculate how many laser sweeps will occur during one plotParam scan               
    tReconstruct=np.linspace(0,np.max(t),scanDetails.numTrigs)
    rReconstruct=1/2*(1-np.cos(2*np.pi*fv*tReconstruct))
    theta=2*np.pi*fr*tReconstruct
    zReconstruct=rReconstruct*np.exp(-1j*theta)
    xReconstruct=np.real(zReconstruct)   # x values of the sampled A-Lines
    yReconstruct=np.imag(zReconstruct)   # y values of the sampled A-Lines
    DebugLog.log("VolumeScan.setupSpiralScan(): Number of laser sweeps that will occur in one volume = %d" % (scanDetails.numTrigs))

    # convert the x-y coordinates of each A-line to an integer between 0 and xPixel/yPixel
    xNorm1=plotParam.xPixel*((xReconstruct/2)+0.5)  
    yNorm1=plotParam.yPixel*((yReconstruct/2)+0.5)
    xNorm=np.round(xNorm1)  
    yNorm=np.round(yNorm1)
        
    # For reformat mode 0: create the 3D array c2, of true/false to give which alines to average for each pixel, and then reformat this into scanDetails.c.
    # For reformat mode 1: find the closest Aline that is within the pixel of interest, and put this index into scanDetail.c3.
    c2=np.zeros((plotParam.xPixel,plotParam.yPixel,scanDetails.numTrigs))   
    scanDetails.c3=np.zeros((plotParam.xPixel,plotParam.yPixel))
    for i1 in range(0,plotParam.xPixel):  
        a=xNorm==i1    # vector of which Alines are at this xPixel
        for i2 in range(0,plotParam.yPixel):
            b=yNorm==i2  # vector of which Alines are at this yPixel
            c1=np.logical_and(a, b)  # c1 is a vector that contains trues for each Aline that is at the pixel of interest
            c2[i1,i2,:]=c1  # c2 is a 3D array of x,y, true/falses for each Aline
            if np.count_nonzero(c1)>0:
                scanDetails.c3[i1,i2]=np.nonzero(c1)[0][0]
    cSum=np.sum(c2,axis=2)
    scanDetails.cTile=np.transpose(np.tile(cSum,[plotParam.zPixel,1,1]),(1, 2, 0))          
    scanDetails.c=np.reshape(c2,(plotParam.xPixel*plotParam.yPixel,scanDetails.numTrigs))            
    scanDetails.c3=np.uint(np.reshape(scanDetails.c3,(plotParam.xPixel*plotParam.yPixel)))
    totalNumPtstoPlot=np.count_nonzero(scanDetails.c3)
    print('plot size=', scanDetails.c3.shape,' totalNumPtstoPlot=',totalNumPtstoPlot)
    
def makeImgSliceFromVolume(volData, widthStep, norms=(0, 120), autoNorm=True, correctAspectRatio=True):
    imgData = volData.volumeImg[widthStep, :, :]
    imgData = imgData.transpose()
            
    xRes = volData.xRes
    zRes = volData.zRes
    
    # remap to 8 bit
    if(autoNorm):
        nL = np.min(imgData)
        nH = np.max(imgData)
    else:
        nL = norms[0]
        nH = norms[1]
        
    imgData = np.floor(255*(imgData - nL)/(nH - nL))
    imgData = np.clip(imgData, 0, 255)
    # imgData = np.reshape((imgData.shape[1], imgData.shape[2]))
    
    # correct the aspect ratio
    if xRes != zRes and correctAspectRatio:    
        imgData = BScan.correctImageAspectRatio(imgData, xRes, zRes)

    return imgData
    # imgData = np.transpose(imgData)
    

class EnFaceProjType(Enum):
    AVERAGE = 0
    MAX = 1
        
def makeEnfaceImgSliceFromVolume(volData, zStep, zDepth, projType=EnFaceProjType.AVERAGE, norms=(0, 120), autoNorm=True, correctAspectRatio=True):
    #imgData = volData.volumeImg[widthStep, :, :]
    shp = volData.volumeImg.shape
    
    zStart = max(zStep - zDepth // 2, 0)
    zEnd = min(zStep + zDepth // 2, shp[1]-1)
    
    if DebugLog.isLogging:
        DebugLog.log("makeEnfaceImgSliceFromVolume shp=%s zStart= %d zEnd= %d" % (repr(shp), zStart, zEnd))
        
    imgData = volData.volumeImg[:, :, zStart:zEnd+1]
    if zStart != zEnd:
        if projType == EnFaceProjType.AVERAGE:
            imgData = np.mean(imgData, 2)
        elif projType == EnFaceProjType.MAX:
            imgData = np.max(imgData, 2)
    else:
        imgData = imgData[:, :, 0]
        
    if (autoNorm):
        nL = np.min(imgData)
        nH = np.max(imgData)
    else:
        nL = norms[0]
        nH = norms[1]
        
    imgData = np.floor(255*(imgData - nL)/(nH - nL))
    imgData = np.clip(imgData, 0, 255)
            
    xRes = volData.xRes
    yRes = volData.yRes

    DebugLog.log("makeEnfaceImgSliceFromVolume imgData.shape=%s xRes= %f yRes= %f" % (repr(imgData.shape), xRes, yRes))
    
    # correct the aspect ratio
    if xRes != yRes and correctAspectRatio:    
        imgData = correctImageAspectRatio(imgData, xRes, yRes)
    else:
        imgData = np.require(imgData, np.uint8)

    return imgData
    # imgData = np.transpose(imgData)
    
    
# return command for raster volume scan of given frame number
def makeVolumeScanCommand(scanParams, frameNum, mirrorDriver, OCTtriggerRate):
    bscansPerFrame = scanParams.volBscansPerFrame
    framesPerScan = scanParams.widthSteps // bscansPerFrame
    daqOutput = None
    frameNum = np.mod(frameNum, framesPerScan)
    for n in range(0, bscansPerFrame):
        frameOffset = frameNum*bscansPerFrame + n
        if DebugLog.isLogging:
            DebugLog.log("makeVolumeScanCommand: frameNum= %d widthOffset= %g" % (frameNum // bscansPerFrame + n, scanParams.widthOffset))
        
        (cmd_x, cmd_y) = BScan.makeBscanCommand(scanParams, mirrorDriver, OCTtriggerRate, frameOffset)
        daqOutputTmp = np.vstack((cmd_x, cmd_y))
        if daqOutput is None:
            daqOutput = copy.copy(daqOutputTmp)
        else:
            daqOutput = np.hstack((daqOutput, daqOutputTmp))

    return daqOutput    
    
def getNumTrigs(scanParams, scanDetails, OCTtrigRate, mirrorDriver):
    scanP = scanParams
    
    if scanP.pattern == ScanPattern.spiral or scanP.pattern == ScanPattern.wagonWheel:
        return scanDetails.numTrigs
    
    bscansPerFrame = scanP.volBscansPerFrame
    galv = mirrorDriver

    reverseTrigs = galv.reverseTime * OCTtrigRate
    flybackTrigs = galv.flybackTime * OCTtrigRate
    startupTrigs = galv.settleTime * OCTtrigRate
    if bscansPerFrame > 1:
        if scanP.pattern == ScanPattern.rasterFast:
            trigsPerBscan = int(np.ceil(scanP.lengthSteps + reverseTrigs + flybackTrigs + startupTrigs))
        elif scanP.pattern == ScanPattern.rasterSlow:
            trigsPerPt = np.ceil(galv.settleTime * OCTtrigRate)
            trigsPerBscan = int(np.ceil(scanP.lengthSteps * trigsPerPt + reverseTrigs + flybackTrigs))
        elif scanP.pattern == ScanPattern.bidirectional:
            trigsPerBscan = int(np.ceil(scanP.lengthSteps + reverseTrigs))
    else:
        if scanP.pattern == ScanPattern.rasterFast or scanP.pattern == ScanPattern.bidirectional:
            trigsPerBscan = scanP.lengthSteps
        elif scanP.pattern == ScanPattern.rasterSlow:
            trigsPerPt = galv.settleTime * OCTtrigRate
            trigsPerBscan = int(np.ceil(scanP.lengthSteps * trigsPerPt))

    numTrigs = trigsPerBscan * bscansPerFrame
    
    return numTrigs
    
def calcuateAngledBScans(plotParam):         
    """
    This function calculates the indices for the desired B-scans
    (plotParam.xSliceList,plotParam.ySliceList)
    and makes a blank array ready for the image (plotParam.bScanArray)
    """
#        # calculate x and y arrays for each b-scan slice and build a blank Bscan array
#        angles=[0,30,60,90,120]
#        plotParam.spacer=4
#        plotParam.sliceWidth=np.round(np.average([plotParam.xPixel,plotParam.yPixel]))
#        plotParam.xSliceList=np.zeros((len(angles),plotParam.sliceWidth))        
#        plotParam.ySliceList=np.zeros((len(angles),plotParam.sliceWidth)) 
#        plotParam.bScanArray=np.zeros((len(angles)*(plotParam.sliceWidth+plotParam.spacer),plotParam.zPixel))
#
#        for i in range(len(angles)):                       
#            # define the endpoints along the ellipse
#            xEnd=np.real(np.exp(-1j*angles[i]*(2*np.pi/360)))
#            yEnd=np.imag(np.exp(-1j*angles[i]*(2*np.pi/360)))
#            xEndNorm=[np.round(plotParam.xPixel*(xEnd/2+0.5)),np.round(plotParam.xPixel*(-1*xEnd/2+0.5))]  #convert to an integer between 0 and xPixel
#            yEndNorm=[np.round(plotParam.yPixel*(yEnd/2+0.5)),np.round(plotParam.yPixel*(-1*yEnd/2+0.5))]  #convert to an integer between 0 and xPixel
#            plotParam.xSliceList[i,:]=np.linspace(xEndNorm[0],xEndNorm[1],plotParam.sliceWidth)-1
#            plotParam.ySliceList[i,:]=np.linspace(yEndNorm[0],yEndNorm[1],plotParam.sliceWidth)-1
    #resample volume at arbitrary angles        
#        angles=np.array([0, 30, 60, 90])
    angles=np.array([0,30,60,90,120,150])
    angleRad=angles*np.pi/180
    maxpixel=np.ceil(np.sqrt(plotParam.xPixel**2+plotParam.yPixel**2)/2)        
    plotParam.spacer=4
    plotParam.sliceWidth=maxpixel*2
    plotParam.xSliceList=np.zeros((len(angles),plotParam.sliceWidth))        
    plotParam.ySliceList=np.zeros((len(angles),plotParam.sliceWidth)) 
    
    r=np.arange(-maxpixel,maxpixel)
    for i in range(len(angleRad)):
        xr=r*np.cos(angleRad[i])+plotParam.xCenter
        yr=r*np.sin(angleRad[i])+plotParam.yCenter
        xrint=np.rint(xr)
        yrint=np.rint(yr)
        
        condition1=np.logical_and(xrint>=0, yrint>=0)
        condition2=np.logical_and(xrint<=(plotParam.xPixel-1), yrint<=(plotParam.yPixel-1))
        condition=np.logical_and(condition1,condition2)
        xIndex=np.int16(np.extract(condition, xrint))
        yIndex=np.int16(np.extract(condition, yrint))
#            print(xIndex.shape)
        plotParam.xSliceList[i,:xIndex.shape[0]]=xIndex
        plotParam.ySliceList[i,:yIndex.shape[0]]=yIndex      
       
def reformatScan(scanDetails,plotParam,oct_dataMag):   
    """
    This function takes the incoming block of Alines and reformats them into a 3D array (data3D).
    reformatMode 0 uses the dot product, which averages overlapping Alines, but takes longer
        scanDetails.c is an array of size numTrigs by xPixels*yPixels containing 0's and 1's.
        scanDetails.cTile is an array of size xPixels*yPixels by zPixels. 
        Each zPixel column contains the number of Alines that are being averaged for that column.
    reformatMode 1 uses an index which is a lot faster, but doesn't do averaging.
        scanDetails.c3 is an array of size xPixels*yPixels by 1, that indexes each Aline to an Aline in the collected OCT data array
    """
    reformatMode=1
    if reformatMode==0:
        datasum21=np.dot(scanDetails.c,oct_dataMag)            
        datasum22=np.reshape(datasum21,(plotParam.xPixel,plotParam.yPixel,plotParam.zPixel))
        data3D=np.nan_to_num(np.divide(datasum22,scanDetails.cTile))+1     
    else:
        plotParam.zPixel=oct_dataMag.shape[1]
        datasum23=oct_dataMag[scanDetails.c3,:]
        datasum24=np.reshape(datasum23,(plotParam.xPixel,plotParam.yPixel,plotParam.zPixel))
        data3D=np.nan_to_num(datasum24)+1     
    return data3D
        
def plotScan(plotParam,data3D, procOpts):

    """
    This simply takes the 3D data set and creates two images (surfacePlot,bScanPlot).               
    """
    # create surface plot
    v21_log=20*np.log10(np.clip(data3D, 1, np.inf))
    threshold=(procOpts.thresholdEnFace/100)*20*np.log10(2**16)
#    print('threshold',threshold,np.min(v21_log),np.max(v21_log),np.mean(v21_log))        
    v21Diff1=v21_log-threshold
    v21Diff2=scipy.stats.threshold(v21Diff1,threshmin=0, newval=2**63)
    v3=np.argmin(v21Diff2,axis=2)
    v3_sumAll=np.sum(v21_log,axis=2)           
    
    # create summed voxel projection
    centerDepth=np.mean(v3[plotParam.xCenter-plotParam.rangeCenter:plotParam.xCenter+plotParam.rangeCenter,plotParam.yCenter-plotParam.rangeCenter:plotParam.yCenter+plotParam.rangeCenter])
    if centerDepth+procOpts.enFace_avgDepth>plotParam.zPixel:
        cutoffDepth=plotParam.zPixel
    else:
        cutoffDepth=centerDepth+procOpts.enFace_avgDepth
    v22=scipy.stats.threshold(v21Diff1,threshmin=0, newval=0)
    
    v4=np.sum(v22[:,:,0:cutoffDepth],axis=2)

    # shape the array to make a proportional image
    matrix=np.zeros((2,2))
    if plotParam.xPixelsize<=plotParam.yPixelsize:
        matrix[0,0]=1
        matrix[1,1]=plotParam.xPixelsize/plotParam.yPixelsize     
        xOut=data3D.shape[0]
        yOut=int(data3D.shape[1]*plotParam.yPixelsize/plotParam.xPixelsize)             
        pixelSize=plotParam.xPixelsize    
    else:
        matrix[0,0]=plotParam.yPixelsize/plotParam.xPixelsize                  
        matrix[1,1]=1
        yOut=data3D.shape[1]
        xOut=int(data3D.shape[0]*plotParam.xPixelsize/plotParam.yPixelsize)
        pixelSize=plotParam.yPixelsize    
    v3A=scipy.ndimage.interpolation.affine_transform(v3,matrix,output_shape=(xOut,yOut))        
    v4A=scipy.ndimage.interpolation.affine_transform(v4,matrix,output_shape=(xOut,yOut))        
    v3_sumAllA=scipy.ndimage.interpolation.affine_transform(v3_sumAll,matrix,output_shape=(xOut,yOut))        
    # create a scale bar
    plotParam.pixelSize=pixelSize
#        print(pixelSize,xOut,yOut)
    scaleLength=(xOut/8)*pixelSize  # scale bar should be roughly 1/10th of the x dimension size
#        print('scalelength',scaleLength,xOut/10)
    base=50    #round to nearest 10 um
    scaleLength=int(base*round(float(scaleLength)/base))
    scaleLengthNumPix=int(scaleLength/pixelSize)
    scaleWidthNumPix=int(yOut/60)+1
    scaleXstart=int(xOut*(7/8))
    scaleYstart=int(yOut*(19/20))     
#    print('Scale Bar = ',scaleLength, ' um')
#    print(scaleLengthNumPix,scaleWidthNumPix,scaleXstart,scaleYstart)
 
    # normalize the surface and summed voxel projection plots            
    surfacePlot=np.uint8((v3A-np.min(v3A))/(np.max(v3A)-np.min(v3A))*255) # normalize the range of the surface
    summedVoxelProjectionPlot=np.uint8((v4A-np.min(v4A))/(np.max(v4A)-np.min(v4A))*255) # normalize the range of the plot
    simpleSumPlot=np.uint8((v3_sumAllA-np.min(v3_sumAllA))/(np.max(v3_sumAllA)-np.min(v3_sumAllA))*255) 
    surfacePlot=np.concatenate((surfacePlot,summedVoxelProjectionPlot,simpleSumPlot))

    # plot the scale bar on the image
    surfacePlot[scaleYstart:(scaleYstart+scaleWidthNumPix),scaleXstart:(scaleXstart+scaleLengthNumPix)]=255                    
   
    # now create bscans at different rotation angles
    plotParam.bScanArray=np.zeros((plotParam.xSliceList.shape[0]*(plotParam.sliceWidth+plotParam.spacer),plotParam.zPixel))
    for i in range(plotParam.xSliceList.shape[0]):
        xSlice=plotParam.xSliceList[i,:]
        ySlice=plotParam.ySliceList[i,:]
        plotParam.bScanArray[i*(plotParam.sliceWidth+plotParam.spacer):i*(plotParam.sliceWidth+plotParam.spacer)+plotParam.sliceWidth,0:plotParam.zPixel]=v21_log[xSlice.astype(int),ySlice.astype(int),0:plotParam.zPixel]       
    
    # shape the array to make a proportional image
    matrix=np.zeros((2,2))
    if plotParam.xPixelsize<=plotParam.zPixelsize:
        matrix[0,0]=1
        matrix[1,1]=plotParam.xPixelsize/plotParam.zPixelsize     
        xOut=plotParam.bScanArray.shape[0]
        zOut=int(plotParam.bScanArray.shape[1]*plotParam.zPixelsize/plotParam.xPixelsize)             
    else:
        matrix[0,0]=plotParam.zPixelsize/plotParam.xPixelsize                  
        matrix[1,1]=1
        zOut=plotParam.bScanArray.shape[1]
        xOut=int(plotParam.bScanArray.shape[0]*plotParam.xPixelsize/plotParam.zPixelsize)
    plotParam.bScanArray1=scipy.ndimage.interpolation.affine_transform(plotParam.bScanArray,matrix,output_shape=(xOut,zOut))        
        
#        plotParam.bScanArray1=plotParam.bScanArray
 
    nL = procOpts.normLow
    nH = procOpts.normHigh
    
    plotParam.bScanArray1 = (plotParam.bScanArray1 - nL) / (nH- nL)
    plotParam.bScanArray1 = np.clip(plotParam.bScanArray1, 0, 1)
    bscan1= plotParam.bScanArray1*255
    bScanPlot=np.transpose(np.uint8(bscan1))
    bScanPlot16b=plotParam.bScanArray1*65535
    bScanPlot16b = np.transpose(np.uint16(bScanPlot16b))
    
    return surfacePlot,bScanPlot,bScanPlot16b


def processDataSpiralScan(oct_data_mag, procOpts, scanDetails, plotParam):
    DebugLog.log("VolumeScan.processDataSpiralScan(): oct_data_mag.shape=(%d, %d)" % (oct_data_mag.shape))

    data3D=reformatScan(scanDetails,plotParam,oct_data_mag) # convert 2D array of A-lines in to 3D dataset with the proper orientation
    volDataIn = VolumeData()
#    volDataIn.scanParams = scanParams
    volDataIn.volumeImg = np.uint16(data3D)                
    volDataIn.zPixSize = procOpts.zRes     # z pixel size in um
    volDataIn.xPixSize = plotParam.xPixelsize    # x pixel size in um
    volDataIn.yPixSize = plotParam.yPixelsize    # y pixel size in um 
        
    [surfacePlot,bScanPlot,bScanPlot16b]= plotScan(plotParam,data3D, procOpts)  # generate the surface views and b-scan slice images           

    # paste the different plots together
    heightDiff=surfacePlot.shape[0]-bScanPlot.shape[0]
    if heightDiff<0:
        addZeros=np.uint8(np.zeros((-heightDiff,surfacePlot.shape[1])))
        surfacePlot1=np.vstack((surfacePlot,addZeros))
        bothPlots=np.hstack((surfacePlot1,bScanPlot))                    
    else:
        addZeros=np.uint8(np.zeros((heightDiff,bScanPlot.shape[1])))
        bScanPlot1=np.vstack((bScanPlot,addZeros))         
        bothPlots=np.hstack((surfacePlot,bScanPlot1))                    
    
    SpiralData = SpiralScanData()   
    SpiralData.bothPlots = bothPlots
    SpiralData.surfacePlot = surfacePlot
    SpiralData.bscanPlot = bScanPlot
    SpiralData.bscanPlot_16b = bScanPlot16b
    volDataIn.spiralScanData = SpiralData
    # procData.data3D = data3D
    
    data3D = 20*np.log10(data3D[:,:,::-1])
    nL = procOpts.normLow
    nH = procOpts.normHigh
    data3D = np.clip(data3D, nL, nH)
    data3D = np.uint16(65535*(data3D - nL)/(nH - nL))
    volDataIn.volumeImg = data3D
    
    return volDataIn

def processDataWagonWheelScan(rawData, procOpts, scanDetails, plotParam):
    pass
    
def processData(oct_data_mag, scanParams, mirrorDriver, OCTtrigRate, procOpts, volDataIn, frameNum, scanDetails=None, plotParam=None):
    scanP = scanParams
    bscansPerFrame = scanP.volBscansPerFrame
    if scanP.pattern == ScanPattern.spiral:
        volDataIn = processDataSpiralScan(oct_data_mag, procOpts, scanDetails, plotParam)
    elif scanP.pattern == ScanPattern.wagonWheel:
        volDataIn = processDataWagonWheelScan(oct_data_mag, procOpts, scanDetails, plotParam)
    else:
        if bscansPerFrame > 1:
            oct_data_tmp = copy.copy(oct_data_mag)
            shp = oct_data_tmp.shape
            framesPerScan = scanParams.widthSteps // bscansPerFrame
            frameInScan = np.mod(frameNum, framesPerScan)
            
            # trigsPerBscan = shp[0] // bscansPerFrame
            mirr = mirrorDriver
            reverseTrigs = np.ceil(mirr.reverseTime * OCTtrigRate)
            flybackTrigs = mirr.flybackTime * OCTtrigRate
            startupTrigs = mirr.settleTime * OCTtrigRate
            if scanP.pattern == ScanPattern.rasterFast:
                extraTrigs = flybackTrigs + reverseTrigs + startupTrigs
            elif scanP.pattern == ScanPattern.rasterSlow:
                extraTrigs = flybackTrigs + reverseTrigs 
            elif scanP.pattern == ScanPattern.bidirectional: 
                extraTrigs = reverseTrigs
                
            trigsPerBscan = scanP.lengthSteps + extraTrigs + procOpts.biDirTrigVolAdj
            
            if DebugLog.isLogging:
                DebugLog.log("VolumeScan processData(): shp= %s trigsPerBscan= %d exraTrigs= %d" % (repr(shp), trigsPerBscan, extraTrigs))
                
            isLogging = DebugLog.isLogging
            DebugLog.isLogging = False  # turn off loggig for for loop
                
            for n in range(0, bscansPerFrame):
                if scanP.pattern == ScanPattern.bidirectional and (np.mod(n, 2) != 0):
                    idx1 = np.floor(trigsPerBscan*n) + procOpts.biDirTrigVolFix
                    idx2 = idx1 + scanP.lengthSteps
                    rawData.oct_data_mag = oct_data_tmp[idx2:idx1:-1, :]
                else:
                    idx1 = trigsPerBscan*n
                    idx2 = idx1 + trigsPerBscan - extraTrigs
                    
                    oct_data_mag = oct_data_tmp[idx1:idx2, :]
#                if DebugLog.isLogging:
#                    DebugLog.log("VolumeScan processData(): n= %d idx1= %d idx2= %d" % (n, idx1, idx2))
                
                img16b, img8b = BScan.makeBscanImage(oct_data_mag, scanParams, procOpts.zRes, procOpts.normLow, procOpts.normHigh, correctAspectRatio=True)
                if volDataIn is None:
                    volDataIn = VolumeData()
                    volDataIn.scanParams = scanParams
                    volDataIn.volumeImg = np.zeros((scanParams.widthSteps, img16b.shape[0], img16b.shape[1]), dtype=np.uint16)                
                    volDataIn.zPixSize = procOpts.zRes     # z pixel size in um
                    volDataIn.xPixSize = img16b.shape[1]/scanParams.length     # x pixel size in um
                    volDataIn.yPixSize = scanParams.widthSteps/scanParams.width     # y pixel size in um 
                
                # DebugLog.log("VolumeScan processData(): step= %d img16b max= %d min=%d " % (n + bscansPerFrame*frameInScan, np.max(img16b), np.min(img16b)))
                volDataIn.volumeImg[n + bscansPerFrame*frameInScan, :, :] = img16b
                
            DebugLog.isLogging = isLogging

        else:
            img16b, img8b = BScan.makeBscanImage(oct_data_mag, scanParams, procOpts.zRes, procOpts.normLow, procOpts.normHigh, correctAspectRatio=True)
            
            if volDataIn is None:
                volDataIn = VolumeData()
                volDataIn.scanParams = scanParams
                volDataIn.volumeImg = np.zeros((scanParams.widthSteps, img16b.shape[0], img16b.shape[1]), dtype=np.uint16)                
                volDataIn.zPixSize = procOpts.zRes     # z pixel size in um
                volDataIn.xPixSize = img16b.shape[1]/scanParams.length     # x pixel size in um
                volDataIn.yPixSize = scanParams.widthSteps/scanParams.width     # y pixel size in um 
                    
            
            DebugLog.log("VolumeScan processData(): frameNum= %d img16b max= %d min=%d " % (frameNum, np.max(img16b), np.min(img16b)))
            volDataIn.volumeImg[frameNum, :, :] = img16b
            
        if DebugLog.isLogging:
            DebugLog.log("VolumeScan processData(): volumeImg.shape= " + repr(volDataIn.volumeImg.shape))
        
    return volDataIn
    
# save the processed data of this protocol
def saveVolumeData(volData, saveDir, saveOpts, scanNum):
    volImg = volData.volumeImg
    volImg = np.require(volImg, np.uint16)
    fileName = 'Volume'
    if saveOpts.subject != '':
        fileName = saveOpts.subject + 'Volume'
        
    if not saveOpts.saveOnlyMostRecentFrame:
        fileName = fileName + "_%.5d" % scanNum
            
    fileName = fileName + '.tiff'
    filePath = os.path.join(saveDir, fileName)
    
    tifffile.imsave(filePath, volImg)    


"""
    VolScanCollectFcn - function that is intended to be called by raw data collection loop (see OCTFPGAProcessingIntreface.LV_DLL_BGProcess_Adaptor)
    oct_hw is a LV_DLL_Interface
"""    

def VolScanCollectFcn(oct_hw, frameNum, extraArgs):
    t1 = time.time()
    scanParams = extraArgs[0]
    mirrorDriver = extraArgs[1]
    zROI = extraArgs[2]
    testDataDir =  extraArgs[3]
    scanDetails =  extraArgs[4]
    
    OCTtrigRate = oct_hw.GetTriggerRate()
    
    bscansPerFrame = scanParams.volBscansPerFrame
    framesPerScan = scanParams.widthSteps // bscansPerFrame   
    if frameNum >= framesPerScan and not scanParams.continuousScan:
        return None, None
        
    chanNames = [mirrorDriver.X_daqChan, mirrorDriver.Y_daqChan]
    outputRate = mirrorDriver.DAQoutputRate
    trigChan = mirrorDriver.trig_daqChan
    if not oct_hw.IsOCTTestingMode():
        from DAQHardware import DAQHardware
        daq = DAQHardware()
    
    if scanParams.pattern == ScanPattern.spiral or scanParams.pattern == ScanPattern.wagonWheel:
        mirrorOut = scanDetails.mirrOut
        startTrigOffset = 0
    else:
        isLogging = DebugLog.isLogging
        DebugLog.isLogging = False
        mirrorOut = makeVolumeScanCommand(scanParams, frameNum, mirrorDriver, OCTtrigRate)
        startTrigOffset = int(np.round(OCTtrigRate*mirrorDriver.settleTime))
        DebugLog.isLogging = isLogging

    DebugLog.log("VolScanCollectFcn: setup time= %0.1f ms" % (1000*(time.time() - t1)))
    
    t2 = time.time()
    if not oct_hw.IsDAQTestingMode():
        # setup the analog output DAQ device
        daq.setupAnalogOutput(chanNames, trigChan, outputRate, mirrorOut.transpose())        
        daq.startAnalogOutput()
    
    DebugLog.log("VolScanCollectFcn: analog output time= %0.1f ms" % (1000*(time.time() - t2)))
    
    # setup and grab the OCT data
    numTrigs = getNumTrigs(scanParams, scanDetails, OCTtrigRate, mirrorDriver)
    if oct_hw.IsOCTTestingMode():
        packedData = OCTCommon.loadRawData(testDataDir, frameNum, dataType=0)
    else:
        t3 = time.time()
        if oct_hw.setupNum == 4:
            err, packedData = oct_hw.AcquireOCTDataMagOnly(numTrigs, zROI, startTrigOffset)
            rawData = VolumeRawData(frameNum, packedData=packedData)
        else:
            err, oct_data = oct_hw.AcquireOCTDataFFT(numTrigs, zROI, startTrigOffset)
            rawData = VolumeRawData(frameNum, oct_data=oct_data)
            
        DebugLog.log("VolScanCollectFcn: acquire time= %0.1f ms" % (1000*(time.time() - t3)))
    
    t4 = time.time()
    if not oct_hw.IsDAQTestingMode():
        daq.waitDoneOutput()
        daq.stopAnalogOutput()
        daq.clearAnalogOutput()
        
    DebugLog.log("VolScanCollectFcn: analog wait, stop and clear time= %0.1f ms" % (1000*(time.time() - t4)))
        
    rawData.collectTime = time.time() - t1
    return rawData, mirrorOut

"""
    VolScanProcessingProcess - function intended to be run as background process that processes volume scan data 
"""    
def VolScanProcessingProcess(scanParams, zROI, procOpts, mirrorDriver, OCTtrigRate, scanDetails, plotParam, rawDataQ, procDataQ, procRawDataQ, msgQ, statusQ):
    shutdown = False
    volData = None
    
    frameNum = 0
    putTimeout = False  # indicates there was a timeout attempting to send data to queue
    while not shutdown:
        try:
            if not putTimeout:
                gotData = False
                if not rawDataQ.empty():
                    OCTdataMag = None
                    rawData = rawDataQ.get(timeout=0.25)
                    
                    if rawData is not None and isinstance(rawData, VolumeRawData):
                        DebugLog.log("VolScanProcessingProcess(): got raw data")
                        # process the data
                        t1 = time.time()
        
                        frameNum = rawData.frameNum
                        if rawData.packedData is not None:
                            (oct_data, mag, phase) = octfpga.unpackData(rawData.packedData)
                            DebugLog.log("VolScanProcessingProcess(): unpack time = %0.1f ms " % (1000*(time.time() - t1)))
                            mag = mag*(2**8)
                        else:
                            oct_data = rawData.oct_data
                            mag = np.abs(oct_data)
                            phase = None
                        
                        volData = processData(mag, scanParams, mirrorDriver, OCTtrigRate, procOpts, volData, frameNum, scanDetails, plotParam)
                        volData.frameNum = rawData.frameNum
                        volData.collectTime = rawData.collectTime
                        volData.processTime = time.time() - t1
                        gotData = True
                        
            # send processsed data to main program
            if gotData:
                try:
                    if procRawDataQ is not None and rawData is not None:
                        procRawDataQ.put(rawData, timeout=0.25)
                        rawData = None
                        
                    procDataQ.put(volData, timeout=0.25)   # send this data last because it contains frame number which client uss to detect whether acquisition is complete                                
                    putTimeout = False
                except queue.Full:
                    DebugLog.log("VolScanProcessingProcess(): queue.Full exception")
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
        while not msgQ.empty():
            msg = msgQ.get()
            msgType = msg[0]
            if msgType == 'shutdown':
                shutdown = True
            elif msgType == 'procOpts':
                procOpts = msg[1]
                
        sys.stdout.flush()  # flush all debugging output to console
        time.sleep(0.005)  # sleep for 5 ms to reduce CPU load

def handleStatusMessage(statusMsg):
    err = False
    # if statusMsg.msgSrc == OCTCommon.StatusMsgSource.COLLECTION:
    if statusMsg.msgType == OCTCommon.StatusMsgType.ERROR:
        err = True
    elif statusMsg.msgType == OCTCommon.StatusMsgType.DAQ_OUTPUT:
        pass
        
    return err
    
"""
    runVolScanMultiProcess - run a volume 
"""    
def runVolScanMultiProcess(appObj, testDataDir, scanParams, zROI, plotParam, scanDetails, procOpts, saveOpts, numFrames, framesPerScan):
    mirrorDriver = appObj.mirrorDriver
    rset = True
    isSaveDirInit = False
    oct_hw = appObj.oct_hw
    
    #procDataQ = mproc.Queue(4)        
#    procMsgQ = mproc.Queue(4)        
 #   dataQ = oct_hw.
    # start up the processing process
  #  procProc = mproc.Process(target=MscanProcessingProcess, args=(audioParams, scanParams, zROI, regionMscan, procOpts, dataQ, procDataQ, procMsgQ), daemon=True)
   # procProc.start()
    startTime = time.time()    
    DebugLog.log("runVolScanMultiProcess: new acquisiiton")
    oct_hw.NewAcquisition()
    oct_hw.SetSendExtraInfo(False)   # do not send mirror output
    DebugLog.log("runVolScanMultiProcess: setting acquire function")
    oct_hw.SetAcqFunction(VolScanCollectFcn)
    extraArgs = [scanParams, mirrorDriver, zROI, testDataDir, scanDetails]
    DebugLog.log("runVolScanMultiProcess: setting acquire functiona args")
    oct_hw.SetAcqFunctionArgs(extraArgs)
    
    volData = None
    startAcq = True
    bscansPerFrame = scanParams.volBscansPerFrame

    DebugLog.log("runVolScanMultiProcess: cleaning status message log")
    statusMsg = oct_hw.GetStatus()
    while statusMsg is not None:
        DebugLog.log("runVolScanMultiProcess: got status message type=" + repr(statusMsg.msgType))
        err = handleStatusMessage(statusMsg)
        statusMsg = oct_hw.GetStatus()        
    
    procDataQ = mproc.Queue(3)
    procRawDataQ = None
    if saveOpts.saveRaw:   # if saving raw data, create a raw data queue so that volume scan processing process will resned raw data to this function
        procRawDataQ = mproc.Queue(10)
    msgQ = mproc.Queue(10)
    rawDataQ = oct_hw.rawDataQ
    statusQ = oct_hw.statusQ
    # VolScanProcessingProcess(scanParams, zROI, procOpts, mirrorDriver, OCTtrigRate, scanDetails, plotParam, rawDataQ, procDataQ, procRawDataQ, msgQ, statusQ):
    OCTtrigRate = oct_hw.GetTriggerRate()
    procProcess = mproc.Process(target=VolScanProcessingProcess, args=[scanParams, zROI, procOpts, mirrorDriver, OCTtrigRate, scanDetails, plotParam, rawDataQ, procDataQ, procRawDataQ, msgQ, statusQ], daemon=True)
    DebugLog.log("runVolScanMultiProcess(): starting processing process")
    procProcess.start()
    frameNum = 0
        
    while not appObj.doneFlag and frameNum < (numFrames-1):
        # update parameters in background process
        # start the acquisitio on first loop iteration
        # we don't just do this outside the loop because we haven't loaded function args yet
        if startAcq:  
            DebugLog.log("runVolScanMultiProcess: starting acquisition")
            oct_hw.StartAcquisition() 
            startAcq = False
        
        if not procDataQ.empty():
            DebugLog.log("runVolScanMultiProcess: grabbing data")
            
            data = procDataQ.get()      
            volData = None
            if isinstance(data, VolumeData):
                DebugLog.log("runVolScanMultiProcess: received volume data")
                frameNum = data.frameNum
                volData = data
                appObj.acquisition_progressBar.setValue(round(100*(frameNum+1)/framesPerScan))                
                appObj.volCollectionTime_label.setText("%0.1f ms" % (volData.collectTime * 1000))
                appObj.volProcessTime_label.setText("%0.1f ms" % (volData.processTime * 1000))
                    
                if scanParams.pattern == ScanPattern.spiral:
                    spiralData = volData.spiralScanData
                    img8b = spiralData.bscanPlot
                    img8b = np.require(img8b, dtype=np.uint8)
                    appObj.vol_bscan_gv.setImage(img8b, ROIImageGraphicsView.COLORMAP_HOT, rset)
                    img8b = spiralData.surfacePlot
                    img8b = np.require(img8b, dtype=np.uint8)
                    appObj.vol_plane_proj_gv.setImage(img8b, ROIImageGraphicsView.COLORMAP_HOT, rset)
                    if rset:
                        gvs = (appObj.vol_bscan_gv, appObj.vol_plane_proj_gv)
                        for gv in gvs:
                            hscroll = gv.horizontalScrollBar()
                            hscroll.setSliderPosition(-500)
                            vscroll = gv.verticalScrollBar()
                            vscroll.setSliderPosition(-500)
                            
                    rset = False
                elif scanParams.pattern == ScanPattern.wagonWheel:
                    pass
                else:
                    framesPerScan = scanParams.widthSteps // bscansPerFrame
                    frameInScan = np.mod(frameNum, framesPerScan)
                    img16b = volData.volumeImg[frameInScan*bscansPerFrame, :, :]
                    img8b = np.round(255.0*img16b/65335.0)  # remap image range
                    DebugLog.log("VolumeScan runVolScan(): frame= %d img16b max= %d min=%d " % (frameNum*bscansPerFrame, np.max(img8b), np.min(img8b)))
    
                    img8b = np.require(img8b, dtype=np.uint8)
                    appObj.vol_bscan_gv.setImage(img8b, ROIImageGraphicsView.COLORMAP_HOT, rset)
                    rset = False                
                    
                if frameNum % framesPerScan == 0:
                    appObj.volDataLast = volData
                    appObj.displayVolumeImg3D(volData.volumeImg)  # update the volume image    

            # save the mscan tuning curve
            if appObj.getSaveState():
                if not isSaveDirInit:
                    saveDir = OCTCommon.initSaveDir(saveOpts, 'MScan', scanParams=scanParams, audioParams=audioParams)
                    isSaveDirInit = True
                        
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
                
        procOpts.normLow = appObj.normLow_spinBox.value()
        procOpts.normHigh = appObj.normHigh_spinBox.value()
        procOpts.biDirTrigVolFix = appObj.volBidirTrigFix_spinBox.value()
        procOpts.thresholdEnFace = appObj.thresholdEnFace_verticalSlider.value()
        procOpts.enFace_avgDepth = appObj.enFace_avgDepth_verticalSlider.value()

        msgQ.put(('procOpts', procOpts))  # update processing options

        tElapsed = time.time() - startTime
        tMins = int(np.floor(tElapsed / 60))
        tSecs = int(tElapsed - 60*tMins)
        appObj.timeElapsed_label.setText("%d mins %d secs" % (tMins, tSecs))

        # check for GUI events, particularly the "done" flag
        QtGui.QApplication.processEvents() 
        time.sleep(0.005)
    
    msgQ.put(('shutdown'))  # tell processing process to stop
    DebugLog.log("runBScanMultiProcess: finishd acquiring data")        
    oct_hw.PauseAcquisition()        
    appObj.isCollecting = False
    QtGui.QApplication.processEvents() # check for GUI events
    appObj.finishCollection()    
    
def runVolScan(appObj):
    DebugLog.log("runVolScan")
    appObj.tabWidget.setCurrentIndex(2)
   
    if not appObj.oct_hw.IsOCTTestingMode():
        from DAQHardware import DAQHardware
        daq = DAQHardware()
        
    appObj.doneFlag = False
    appObj.isCollecting = True
    OCTtrigRate = appObj.oct_hw.GetTriggerRate()

    OCTtrigRate = 10000  # make this slower for testing purposes

    mirrorDriver = appObj.mirrorDriver
    rset = True
    
    chanNames = [mirrorDriver.X_daqChan, mirrorDriver.Y_daqChan]
    trigChan = mirrorDriver.trig_daqChan
    outputRate = mirrorDriver.DAQoutputRate
    processMode = OCTCommon.ProcessMode(appObj.processMode_comboBox.currentIndex())

    # Get scan paramaeters either from GUI (running or software process testmode) or from a file (hardware process testmode)
    testDataDir = '' 
    if appObj.oct_hw.IsOCTTestingMode():
        testDataDir = os.path.join(appObj.basePath, 'exampledata', 'VolumeScan')
        if processMode == OCTCommon.ProcessMode.FPGA:
            filePath = os.path.join(testDataDir, 'ScanParams.pickle')
            f = open(filePath, 'rb')
            scanParams = pickle.load(f)
            f.close()
        elif processMode == OCTCommon.ProcessMode.SOFTWARE:        
            appObj.savedDataBuffer = JSOraw.SavedDataBuffer()     # This class holds data imported from a disk file, and loads a test data set (for software processing)
            appObj.savedDataBuffer.loadData(appObj, testDataDir, 'testData.npz')
            scanParams = appObj.getScanParams()
    else:
        # get the scan paramters that the user has entred 
        scanParams = appObj.getScanParams()

    scanDetails = None
    plotParam = None
    zROI = appObj.getZROI()    
    bscansPerFrame = scanParams.volBscansPerFrame
    numFrames = scanParams.widthSteps // bscansPerFrame            
    framesPerScan = scanParams.widthSteps // bscansPerFrame        
    if scanParams.continuousScan:
        numFrames = np.inf
    if scanParams.pattern == ScanPattern.spiral or scanParams.pattern == ScanPattern.wagonWheel:
        numFrames = 1
        framesPerScan = 1
        
    procOpts = ProcOpts()
    procOpts.normLow = appObj.normLow_spinBox.value()
    procOpts.normHigh = appObj.normHigh_spinBox.value()
    procOpts.zRes = appObj.octSetupInfo.zRes
    procOpts.biDirTrigVolFix = appObj.volBidirTrigFix_spinBox.value()
    procOpts.thresholdEnFace=appObj.thresholdEnFace_verticalSlider.value()
    procOpts.enFace_avgDepth=appObj.enFace_avgDepth_verticalSlider.value()

    timePoint1 = time.time()
    if scanParams.pattern == ScanPattern.spiral or scanParams.pattern == ScanPattern.wagonWheel:
        plotParam, scanDetails = setupScan(scanParams, mirrorDriver, zROI, OCTtrigRate, procOpts)
   
    saveOpts = appObj.getSaveOpts()
    if(appObj.multiProcess):
        runVolScanMultiProcess(appObj, testDataDir, scanParams, zROI, plotParam, scanDetails, procOpts, saveOpts, numFrames, framesPerScan)
        return

    isSaveDirInit = False
    timePoint2 = time.time()

    try: 
        frameNum = 0
        scanNum = 0
        
        while not appObj.doneFlag and frameNum < numFrames:
            # reinitialize volume data on first frame
            if frameNum % framesPerScan == 0:
                volData = None
            
            procOpts.normLow = appObj.normLow_spinBox.value()
            procOpts.normHigh = appObj.normHigh_spinBox.value()
            procOpts.thresholdEnFace=appObj.thresholdEnFace_verticalSlider.value()
            procOpts.enFace_avgDepth=appObj.enFace_avgDepth_verticalSlider.value()
                
            if scanParams.pattern == ScanPattern.spiral or scanParams.pattern == ScanPattern.wagonWheel:
                mirrorOut = scanDetails.mirrOut
                startTrigOffset = 0
            else:
                mirrorOut = makeVolumeScanCommand(scanParams, frameNum, mirrorDriver, OCTtrigRate)
                startTrigOffset = int(np.round(OCTtrigRate*mirrorDriver.settleTime))
                
            # plot command to GUI 
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
            pl2.plot(mirrorOut[0, :],mirrorOut[1, :], pen='k')
            labelStyle = appObj.xLblStyle
            pl2.setLabel('bottom', 'X galvo', 'V', **labelStyle)
            labelStyle = appObj.yLblStyle
            pl2.setLabel('left', 'Y galvo', 'V', **labelStyle)
            
            timePoint3 = time.time()

            if not appObj.oct_hw.IsDAQTestingMode():
                # setup the analog output DAQ device
                daq.setupAnalogOutput(chanNames, trigChan, outputRate, mirrorOut.transpose())        
                daq.startAnalogOutput()
            
            # setup and grab the OCT data
            numTrigs = getNumTrigs(scanParams, scanDetails, OCTtrigRate, mirrorDriver)
            if processMode == OCTCommon.ProcessMode.FPGA:
                if appObj.oct_hw.IsOCTTestingMode():
                    oct_data = OCTCommon.loadRawData(testDataDir, frameNum, dataType=0)
                else:
                    err, oct_data = appObj.oct_hw.AcquireOCTDataFFT(numTrigs, zROI, startTrigOffset)
            elif processMode == OCTCommon.ProcessMode.SOFTWARE:
                if appObj.oct_hw.IsOCTTestingMode():
                    ch0_data,ch1_data=JSOraw.getSavedRawData(numTrigs,appObj.dispData.requestedSamplesPerTrig,appObj.savedDataBuffer)
                else:
                    # def AcquireOCTDataRaw(self, numTriggers, samplesPerTrig=-1, Ch0Shift=-1, startTrigOffset=0):
                    samplesPerTrig = appObj.oct_hw.fpgaOpts.SamplesPerTrig*2
                    err, ch0_data,ch1_data = appObj.oct_hw.AcquireOCTDataRaw(numTrigs, samplesPerTrig, startTrigOffset=startTrigOffset)
                
                timePoint4 = time.time()

                oct_data, klin = JSOraw.softwareProcessing(ch0_data,ch1_data,zROI,appObj)
            else:
                QtGui.QMessageBox.critical (appObj, "Error", "Unsuppoted processing mode for current hardware")
            
            timePoint5 = time.time()

            # process the data
            oct_data_mag = np.abs(oct_data)
            volData = processData(oct_data_mag, scanParams, mirrorDriver, OCTtrigRate, procOpts, volData, frameNum, scanDetails, plotParam)
            timePoint6 = time.time()
    
            if scanParams.pattern == ScanPattern.spiral :
                spiralData = volData.spiralScanData
                img8b = spiralData.bscanPlot
                img8b = np.require(img8b, dtype=np.uint8)
                appObj.vol_bscan_gv.setImage(img8b, ROIImageGraphicsView.COLORMAP_HOT, rset)
                img8b = spiralData.surfacePlot
                img8b = np.require(img8b, dtype=np.uint8)
                appObj.vol_plane_proj_gv.setImage(img8b, ROIImageGraphicsView.COLORMAP_HOT, rset)
                if rset:
                    gvs = (appObj.vol_bscan_gv, appObj.vol_plane_proj_gv)
                    for gv in gvs:
                        hscroll = gv.horizontalScrollBar()
                        hscroll.setSliderPosition(-500)
                        vscroll = gv.verticalScrollBar()
                        vscroll.setSliderPosition(-500)
                        
                rset = False
            else:
                img16b = volData.volumeImg[frameNum*bscansPerFrame, :, :]
                img8b = np.round(255.0*img16b/65335.0)  # remap image range
                DebugLog.log("VolumeScan runVolScan(): frame= %d img16b max= %d min=%d " % (frameNum*bscansPerFrame, np.max(img8b), np.min(img8b)))

                img8b = np.require(img8b, dtype=np.uint8)
                appObj.vol_bscan_gv.setImage(img8b, ROIImageGraphicsView.COLORMAP_HOT, rset)
                rset = False                
                
            if not appObj.oct_hw.IsDAQTestingMode():
                daq.stopAnalogOutput()
                daq.clearAnalogOutput()
            timePoint7 = time.time()

            print('timePoint2-timePoint1 = ',timePoint2-timePoint1)
            print('timePoint3-timePoint2 = ',timePoint3-timePoint2)
            print('timePoint4-timePoint3 = ',timePoint4-timePoint3)
            print('timePoint5-timePoint4 = ',timePoint5-timePoint4)
            print('timePoint6-timePoint5 = ',timePoint6-timePoint5)
            print('timePoint7-timePoint6 = ',timePoint7-timePoint6)            

            frameNum += 1
            appObj.acquisition_progressBar.setValue(round(100*frameNum/framesPerScan))
            if appObj.getSaveState():
                if not isSaveDirInit:
                    saveDir = OCTCommon.initSaveDir(saveOpts, 'Volume', scanParams)
                    isSaveDirInit = True
                if saveOpts.saveRaw:
                    if processMode == OCTCommon.ProcessMode.FPGA:
                        OCTCommon.saveRawData(oct_data, saveDir, frameNum-1, dataType=0)
                    elif processMode == OCTCommon.ProcessMode.SOFTWARE:
                        outfile = os.path.join(saveDir, 'testData %d.npz' % (frameNum-1))
                        np.savez_compressed(outfile, ch0_data=ch0_data, ch1_data=ch1_data)
            
            if frameNum % framesPerScan == 0:
                if appObj.getSaveState():
                    saveVolumeData(volData, saveDir, saveOpts, scanNum)

                appObj.volDataLast = volData
                appObj.displayVolumeImg3D(volData.volumeImg)  # update the volume image
                    
                scanNum += 1
                if scanParams.continuousScan:
                    frameNum = 0
                
            # check for GUI events, particularly the "done" flag
            QtGui.QApplication.processEvents() 

        
    except Exception as ex:
#        raise ex
        traceback.print_exc(file=sys.stdout)
        QtGui.QMessageBox.critical (appObj, "Error", "Error during scan. Check command line output for details")
    finally:
        appObj.isCollecting = False
        QtGui.QApplication.processEvents() # check for GUI events
        appObj.finishCollection()