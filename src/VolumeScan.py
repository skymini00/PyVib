# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 12:22:44 2015

@author: OHNS
"""

import OCTCommon
import OCTFPGAProcessingInterface as octfpga
from DebugLog import DebugLog
from scipy import stats
from scipy import signal

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

class SpiralScanData():
    def __init__(self):
        self.bothPlots = None
        self.surfacePlot = None
        self.bscanPlot = None
        self.bscanPlot_16b = None
        
def setupScan(scanParams, mirrorDriver, zROI, OCTtrigRate, procOpts):
    """
    This function should be called to setup the output and reconstructur parameters prior to starting the protocol,
    so that the getMirrorOutput(), processData() functions will return correct values
    """
    yAdjust_val = scanParams.xskew
    phaseAdjust = scanParams.phaseAdjust
    DAQoutputRate = mirrorDriver.DAQoutputRate
    
    #scanMode = self.scanMode_comboBox.currentIndex() + 1
    #scanMode=1 
    DebugLog.log("VolumeScan.setupScan() zROI=%s OCTtrigRate=%d DAQoutputRate=%d" % (repr(zROI), OCTtrigRate, DAQoutputRate))
    # define plotting parameters and then setup plotting algorithms
    plotParam=blankRecord() # create an empty record to pass the plotting parameters in
    plotParam.zROI = zROI
    # plotParam.zROI=(250,275)
    plotParam.zPixel=plotParam.zROI[1]-plotParam.zROI[0]+1              
    plotParam.zPixelsize=procOpts.zRes    
    plotParam.depth=30    # depth of pixels to sum for summed voxel projection
    plotParam.xPixel=scanParams.lengthSteps   # number of pixels in x dimension  
    plotParam.yPixel=scanParams.widthSteps   # number of pixels in y dimension        
    plotParam.xCenter=np.int(plotParam.xPixel/2)
    plotParam.yCenter=np.int(plotParam.yPixel/2)            
    plotParam.rangeCenter=((plotParam.xCenter+plotParam.yCenter)//2)//4
    
    DebugLog.log("EndoSpiralScanProtocol.setupScan() plotParam xPixel= %d yPixel= %d zPixel= %d xCenter= %d yCenter= %d rangeCenter= %d" % (plotParam.xPixel, plotParam.yPixel, plotParam.zPixel, plotParam.xCenter, plotParam.yCenter, plotParam.rangeCenter))
    #  setup galvo voltages and calculate which samples fit into which pixels                 
    scanDetails=blankRecord()  #create an empty record to get the spiral scan parameters        if scanMode==0:
    calcuateAngledBScans(plotParam)
    if scanParams.pattern == ScanPattern.spiral:
        setupSpiralScan(mirrorDriver, scanDetails, plotParam, OCTtrigRate, DAQoutputRate)
    elif scanParams.pattern == ScanPattern.wagonWheel:
        setupWagonWheelScan(mirrorDriver, scanDetails, plotParam, OCTtrigRate, DAQoutputRate)

    return plotParam, scanDetails
    
# TODO: need to implement this function
def setupWagonWheelScan(scanParams, mirrorDriver, scanDetails,plotParam, OCTtrigRate, DAQoutputRate):
    pass

def setupSpiralScan(scanParams, mirrorDriver, scanDetails, plotParam, OCTtrigRate, DAQoutputRate):
    """ 
    This function builds the spiral scan voltages (scanDetails.mirrOut) 
    and the variables that are used when reformatting the collected A-lines
    into the proper 3D format(scanDetails.c,scanDetails.cTile, and scanDetails.c3)
    """
    # make mirror output signals
    yAdjust = scanParams.xskew
    phaseShift = scanParams.phaseAdjust
    fr = scanParams.angularScanFreq
    
#            yAdjust=self.yAdjust_spinBox.value()
#            phaseShift=self.eccentricity_spinBox.value() 
    diameter = scanParams.length
    plotParam.xPixelsize=(diameter/plotParam.xPixel)*1000  # size on one pixel in the x dimension in microns
    plotParam.yPixelsize=(diameter/plotParam.yPixel)*1000  # size on one pixel in the y dimension in microns         
        
    fv=5     # plotParam scan frequency, which scans in and then out, which is actually two volumes
#        fr=280  # angular scan rate (frequency of one rotation)
    rate=2
    voltsPerMM = mirrorDriver.voltsPerMillimeter
    xAdjust = 1
    
    A=voltsPerMM*diameter/2
       
    fs=DAQoutputRate   # galvo output sampling rate
    trigRate = OCTtrigRate
    scanDetails.numTrigs=np.int32(trigRate/fv)     # calculate how many laser sweeps will occur during one plotParam scan               
    t=np.arange(0,np.around(fs/fv))*1/fs
    r=1/2*(1-np.cos(2*np.pi*fv*t))            
    x=xAdjust*A*r*np.cos(2*np.pi*fr*t)
    y=yAdjust*A*r*np.sin(2*np.pi*fr*t+phaseShift*np.pi/180)
    DebugLog.log("EndoSpiralScanProtocol.setupSpiralScan() len(t)= %d len(x)= %d" % (len(t), len(x)))

    scanDetails.mirrOut= np.vstack((x,y))
  
    # reconstruct the image from the array of A-lines
    tReconstruct=np.linspace(0,np.max(t),scanDetails.numTrigs)
    rReconstruct=1/2*(1-np.cos(2*np.pi*fv*tReconstruct))
    theta=2*np.pi*fr*tReconstruct
    zReconstruct=rReconstruct*np.exp(-1j*theta)
    xReconstruct=np.real(zReconstruct)   # x values of the sampled A-Lines
    yReconstruct=np.imag(zReconstruct)   # y values of the sampled A-Lines

    xNorm1=plotParam.xPixel*((xReconstruct/2)+0.5)  
    yNorm1=plotParam.yPixel*((yReconstruct/2)+0.5)
    xNorm=np.round(xNorm1)  #convert to an integer between 0 and xPixel
    yNorm=np.round(yNorm1)
     
    # create the 3D array c2, of true/false to give which alines to average for each pixel, and then reformat this into scanDetails.c
    c2=np.zeros((plotParam.xPixel,plotParam.yPixel,scanDetails.numTrigs))   
    scanDetails.c3=np.zeros((plotParam.xPixel,plotParam.yPixel))
    for i1 in range(0,plotParam.xPixel):  
        a=xNorm==i1    
        for i2 in range(0,plotParam.yPixel):
            b=yNorm==i2
            c1=np.logical_and(a, b)  # c is a vector contain trues for all of the data points in the pixel of interest
            c2[i1,i2,:]=c1
#                indexList=np.nonzero(c1)
#                if len(indexList[0])>=1:
#                    scanDetails.c3[i1,i2]=indexList[0][0]
            distance=(xNorm1-i1)**2+(yNorm1-i2)**2
            if np.min(distance)<np.sqrt(2):
                scanDetails.c3[i1,i2]=np.argmin(distance)
    cSum=np.sum(c2,axis=2)
    scanDetails.cTile=np.transpose(np.tile(cSum,[plotParam.zPixel,1,1]),(1, 2, 0))          
    scanDetails.c=np.reshape(c2,(plotParam.xPixel*plotParam.yPixel,scanDetails.numTrigs))            
    scanDetails.c3=np.uint(np.reshape(scanDetails.c3,(plotParam.xPixel*plotParam.yPixel)))

def makeScanXYcoords(scanParams, galvMirror, volWidthStep=-1, xskew=1.0):
    hsl = scanParams.length/2   # half scan length
    vpm = galvMirror.voltsPerMillimeter
    hsl = vpm * hsl
    
    # rotate the scan
    rot = scanParams.rotation_Z*np.pi/180
    cos_rot = np.cos(rot)
    sin_rot = np.sin(rot)
    x2 = hsl * cos_rot
    x1 = -x2 
    
    y2 = hsl * sin_rot
    y1 = -y2
    
    # offset scan
    lo_vpm = vpm*scanParams.lengthOffset
    xOffset = lo_vpm*cos_rot
    yOffset = lo_vpm*sin_rot
    
    rot_wo = rot + np.pi/2
    width = scanParams.width
    widthOffset = scanParams.widthOffset
    if volWidthStep >=0:
        widthOffset = widthOffset + -width/2 + (volWidthStep*width)/(scanParams.widthSteps + 1)
        
    wo_vpm = vpm*widthOffset
    
    xOffset = xOffset + wo_vpm*np.cos(rot_wo)
    yOffset = yOffset + wo_vpm*np.sin(rot_wo)
    
    x1 = x1*xskew
    x2 = x2*xskew
    
    
    x1 = x1 + xOffset
    x2 = x2 + xOffset
    y1 = y1 + yOffset
    y2 = y2 + yOffset
    
    return (x1, y1, x2, y2)
    
# return the mirror B scam (slice) command
def makeBscanCommand(scanParams, mirrorDriver, OCTtriggerRate, volWidthStep=-1):
    cmd_x = None
    cmd_y = None
    galv = mirrorDriver
    DAQoutputRate = galv.DAQoutputRate
    DebugLog.log("makeBscanCommand(): outputRate= %0.1f trigRate= %0.1f" % (DAQoutputRate, OCTtriggerRate))
    DAQptsPerTrig = DAQoutputRate/OCTtriggerRate
    if scanParams.pattern == ScanPattern.rasterFast or (scanParams.pattern == ScanPattern.bidirectional):
        scanTime = scanParams.lengthSteps / OCTtriggerRate 
        scanPts = np.floor(scanTime * DAQoutputRate)
        
        startupTrigs = np.round(galv.settleTime * OCTtriggerRate)
        startupPts = np.floor(DAQoutputRate*startupTrigs/OCTtriggerRate)
        # startupPts = np.floor(galv.settleTime * DAQoutputRate)
        
        reversePts = np.floor(galv.reverseTime * DAQoutputRate)
        bidirOddFrame = False
        biDirScan = False
        if scanParams.pattern == ScanPattern.bidirectional:
            flybackPts = 0
            bidirOddFrame = (np.mod(volWidthStep, 2) != 0)
            biDirScan = True
            if volWidthStep > 0:
                startupPts = 0
                
        else:
            #flybackPts = min(scanPts, np.floor(galv.flybackTime * DAQoutputRate))
            flybackPts = np.floor(galv.flybackTime * DAQoutputRate)
        DebugLog.log("makeBscanCommand(): scanPts= %d startupPts= %d reversePts= %d flybackPts= %d" % (scanPts, startupPts, reversePts, flybackPts))
        
        # get reversal/flyback to align with OCT triggering
        extraPts = 0
        endPts = reversePts + flybackPts
        extraTrigs = endPts / DAQptsPerTrig  - np.floor(endPts/DAQptsPerTrig)
        if extraTrigs > 0:
            extraPts = int(np.round((1-extraTrigs) *  DAQoutputRate / OCTtriggerRate))
        DebugLog.log("makeBscanCommand(): extraTrigs= %0.3f extraPts= %d" % (extraTrigs, extraPts))
        reversePts = reversePts + extraPts

        xskew = 1.0
        if hasattr(scanParams, 'xskew'):
            xskew = scanParams.xskew
            
        (x1, y1, x2, y2) = makeScanXYcoords(scanParams, galv, volWidthStep, xskew)
        
        if bidirOddFrame:
            (x1, y1, x2, y2) = (x2, y2, x1, y1)
        
        # startup command
        mx = (x2 - x1)/scanPts
        my = (y2 - y1)/scanPts
        x0 = -startupPts*mx + x1
        y0 = -startupPts*my + y1
        
        startup_x = []
        startup_y = []
        if startupPts > 0:
            startup_x = np.linspace(x0, x1, startupPts)
            startup_y = np.linspace(y0, y1, startupPts)
        
        # scan command
        scan_x = np.linspace(x1, x2, scanPts)
        scan_y = np.linspace(y1, y2, scanPts)
        
        # reversa command
        x3 = x2 + mx*(reversePts/2)
        y3 = y2 + my*(reversePts/2)
        x4 = x2
        y4 = y2
        if volWidthStep >= 0 and biDirScan:
            x3 = x2 + mx*(reversePts / 4)
            y3 = y2 + my*(reversePts / 4)
            (x4, y4, x5, y5) = makeScanXYcoords(scanParams, galv, volWidthStep + 1)
            t1 = np.linspace(0, np.pi/2, reversePts // 2)
            if bidirOddFrame:
                (x4, y4) = (x5, y5)
        else:                
            t1 = np.linspace(0, np.pi/2, reversePts)
            
        reverse_x = (x3 - x2)*np.sin(t1) + x2
        reverse_y = (y3 - y2)*np.sin(t1) + y2

        if volWidthStep >= 0 and biDirScan:
            t2 = np.linspace(np.pi/2, np.pi, reversePts // 2)
            reverse_x = np.concatenate((reverse_x, (x3 - x2)*np.sin(t2) + x2))
            reverse_y = np.concatenate((reverse_y, (y3 - y2)*np.sin(t2) + y2))
        
        # flyback command
        flyback_x = []
        flyback_y = []
        if flybackPts > 0:
            t = np.linspace(0, np.pi, flybackPts)
            s = (np.cos(t)+1)/2
            flyback_x = (x3 - x0)*s + x0           
            flyback_y = (y3 - y0)*s + y0
        
        cmd_x = np.concatenate((startup_x, scan_x, reverse_x, flyback_x))
        cmd_y = np.concatenate((startup_y, scan_y, reverse_y, flyback_y))
        
    elif scanParams.pattern == ScanPattern.rasterSlow:
        scanTime = scanParams.lengthSteps / OCTtriggerRate 
        # trigsPerPt = np.ceil(galv.settleTime * OCTtriggerRate)
        DAQptsPerTrig = galvoMirror.settleTime * DAQoutputRate 
        scanPts = np.floor(DAQptsPerTrig * scanParams.lengthSteps)
        DebugLog.log("makeBscanCommand(): DAQptsPerTrig= %0.1f scanPts= %0.1f" % (DAQptsPerTrig, scanPts))

        flybackPts = np.floor(galvoMirror.flybackTime * DAQoutputRate)
        
        (x1, y1, x2, y2) = makeScanXYcoords(galvoMirror)
        
        scan_x = np.zeros(scanPts)
        scan_y = np.zeros(scanPts)
        # scan command
        for n in range(0, scanParams.lengthSteps):
            pt1 = np.floor(n*DAQptsPerTrig)
            pt2 = np.floor((n+1)*DAQptsPerTrig)
            scan_x[pt1:pt2] = x1 + n*(x2 - x1)/scanParams.lengthSteps
            scan_y[pt1:pt2] = y1 + n*(y2 - y1)/scanParams.lengthSteps
        
        # flyback command
        t = np.linspace(0, np.pi, flybackPts)
        s = (np.cos(t)+1)/2
        flyback_x = (x2 - x1)*s + x1           
        flyback_y = (y2 - y1)*s + y1
        
        cmd_x = np.concatenate((scan_x, flyback_x))
        cmd_y = np.concatenate((scan_y, flyback_y))
        
    else:
        DebugLog.log("makeBscanCommand(): unsupported scan pattern '%s'"  % repr(scanParams.pattern))
    
    if cmd_x is not None:
        return np.vstack((cmd_x, cmd_y))
    else:
        return None
    
# return command for entire volume 
def makeVolumeScanCommand(scanParams, frameNum, mirrorDriver, OCTtriggerRate):
    bscansPerFrame = scanParams.volBscansPerFrame
    framesPerScan = scanParams.widthSteps // bscansPerFrame
    daqOutput = None
    frameNum = np.mod(frameNum, framesPerScan)
    for n in range(0, bscansPerFrame):
        frameOffset = frameNum*bscansPerFrame + n
        if DebugLog.isLogging:
            DebugLog.log("makeVolumeScanCommand: frameNum= %d widthOffset= %g" % (frameNum // bscansPerFrame + n, scanParams.widthOffset))
        
        (cmd_x, cmd_y) = makeBscanCommand(scanParams, mirrorDriver, OCTtriggerRate, frameOffset)
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
    reformatMode 1 uses the dot product, which averages overlapping Alines, but takes longer
        scanDetails.c is an array of size numTrigs by xPixels*yPixels containing 0's and 1's.
        scanDetails.cTile is an array of size xPixels*yPixels by zPixels. 
        Each zPixel column contains the number of Alines that are being averaged for that column.
    reformatMode 2 uses an index which is a lot faster, but doesn't do averaging.
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
        
def plotScan(plotParam,data3D, spiralScanThresholdVal, normLow, normHigh):
    """
    This simply takes the 3D data set and creates two images (surfacePlot,bScanPlot).               
    """
    # create surface plot
    v21_log=np.log10(data3D)
#            threshold=np.mean(v21_log,axis=2) 
#            v21Mean1=np.transpose(np.tile(threshold,(oct_data.shape[1],1,1)),(1,2,0))             
#            v21Diff1=v21_log-v21Mean1
    tSlider= spiralScanThresholdVal/100
    print('tSlider',tSlider)       
    threshold=tSlider*np.log10(2**16)
    print('threshold',threshold,np.min(v21_log),np.max(v21_log),np.mean(v21_log))        
    v21Diff1=v21_log-threshold
    v21Diff2=stats.threshold(v21Diff1,threshmin=0, newval=2**63)
    v3=np.argmin(v21Diff2,axis=2)
    v3_sumAll=np.sum(v21_log,axis=2)           
    
    # create summed voxel projection
    centerDepth=np.mean(v3[plotParam.xCenter-plotParam.rangeCenter:plotParam.xCenter+plotParam.rangeCenter,plotParam.yCenter-plotParam.rangeCenter:plotParam.yCenter+plotParam.rangeCenter])
    if centerDepth+plotParam.depth>plotParam.zPixel:
        cutoffDepth=plotParam.zPixel
    else:
        cutoffDepth=centerDepth+plotParam.depth
    v22=stats.threshold(v21Diff1,threshmin=0, newval=0)
    
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
    print('Scale Bar = ',scaleLength, ' um')
    print(scaleLengthNumPix,scaleWidthNumPix,scaleXstart,scaleYstart)
 
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
 
    nL = normLow
    nH = normHigh
    
    plotParam.bScanArray1 = 20 * plotParam.bScanArray1
    plotParam.bScanArray1 = (plotParam.bScanArray1 - nL) / (nH- nL)
    plotParam.bScanArray1 = np.clip(plotParam.bScanArray1, 0, 1)
    bscan1= plotParam.bScanArray1*255
    bScanPlot=np.transpose(np.uint8(bscan1))
    bScanPlot16b=plotParam.bScanArray1*65535
    bScanPlot16b = np.transpose(np.uint16(bScanPlot16b))
    
    return surfacePlot,bScanPlot,bScanPlot16b

def makeBscanImage(oct_data_mag, scanParams, zRes, normLow, normHigh, correctAspectRatio=True):
    # numTrigs = self.scanParams.lengthSteps
    numTrigs = oct_data_mag.shape[0]
    
    # for a slwo scan, reshape and exclude unwanted triggers
    if scanParams.pattern == ScanPattern.rasterSlow: 
        # rawdata = trigsPerPt = np.ceil(galv.settleTime * OCTtriggerRate)
        shp = rawdata.shape
        
        trigsPerStep = numTrigs // scanParams.lengthSteps
        DebugLog.log("VolumeScan makeBscanImage(): rawdata.shape= %s trigsPerStep= %d" %  (repr(rawdata.shape), trigsPerStep))
        # rawdata = np.reshape(rawdata, (trigsPerStep, scanParams.lengthSteps, shp[1], shp[2]))
        
        i1 = trigsPerStep - 1
        i2 = scanParams.lengthSteps * trigsPerStep
        rng = list(range(i1, i2, trigsPerStep))
        oct_data_mag = rawdata[rng, :]
            
    mag = np.clip(oct_data_mag, 1, np.inf)
    mag = 20*np.log10(mag)
        
    DebugLog.log("VolumeScan makeBscanImage(): (before averaging) mag.shape= " + repr(mag.shape))          
    
    # average all trials
    if len(mag.shape) > 2:
        mag = np.mean(mag, 2)
    
    DebugLog.log("VolumeScan makeBscanImage(): (after averaging) mag.shape= " + repr(mag.shape))  
    mag = mag.transpose()
    
    #DebugLog.log("processBScanData: mag min= %g max= %g" % (np.min(mag), np.max(mag)))
     # map the data to 0 ... 1
    nL = normLow
    nH = normHigh
    nRng = nH - nL
    
    DebugLog.log("VolumeScan makeBscanImage(): nL= %g nH=%g nRng= %g" % (nL, nH, nRng))
        
    # remap range to 0...1
    mag = (mag - nL)/nRng  

    # remap range to 0 ... to 2^16 - 1
    mag16b = mag*65535
    mag16b = np.clip(mag16b, 0, 65535)
    
    DebugLog.log("VolumeScan makeBscanImage():  mag16b min= %g max= %g" % (np.min(mag16b), np.max(mag16b)))
    
    # produce image            
    mag16b = np.require(mag16b, 'uint32')

    s = mag16b.shape
    
    sizeX = s[1]
    sizeY = s[0]
    
    sl = scanParams.length
    xRes = sl*1e3/sizeX
    
    DebugLog.log("VolumeScan makeBscanImage(): sizeX = %d sizeY= %d xRes= %g" % (sizeX, sizeY, xRes))
    
    img16b = PIL.Image.fromarray(mag16b, 'I')     # 'I' for 32-bit images, 'P' for 8-bit
    # correct aspect ratio
    if(correctAspectRatio):
        sizeX = int(np.round(sizeX*xRes/zRes))
        sizeX = min(sizeX, 4000)    # guard against images that are too big
        sizeX = max(2, sizeX)
        DebugLog.log("VolumeScan makeBscanImage():  new sizeX = %d" % (sizeX))
        img16b = img16b.resize((sizeX, sizeY), PIL.Image.BILINEAR) 

    imgdata_16b = np.array(img16b)
    imgdata_8b = np.floor(255*imgdata_16b/65535)

    imgdata_8b = np.require(imgdata_8b, np.uint8, 'C')
    imgdata_16b = np.require(imgdata_16b, np.uint16, 'C')
    DebugLog.log("VolumeScan makeBscanImage():  img8b min= %g max= %g" % (np.min(imgdata_8b), np.max(imgdata_8b)))
    DebugLog.log("VolumeScan makeBscanImage():  img16b min= %g max= %g" % (np.min(imgdata_16b), np.max(imgdata_16b)))

    return imgdata_16b, imgdata_8b
    
def processDataSpiralScan(oct_data_mag, procOpts, scanDetails, plotParam):
    DebugLog.log("SpiralScanProtocol.processData()")

    data3D=self.reformatScan(scanDetails,plotParam,oct_data_mag) # convert 2D array of A-lines in to 3D dataset with the proper orientation
   
    [surfacePlot,bScanPlot, bScanPlot16b]= plotScan(plotParam,data3D, procOpts)  # generate the surface views and b-scan slice images           

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
    
    procData = ProcessedData()
    SpiralData = SpiralScanData()
    
    SpiralData.bothPlots = bothPlots
    SpiralData.surfacePlot = surfacePlot
    SpiralData.bscanPlot = bScanPlot
    SpiralData.bscanPlot_16b = bScanPlot16b
    procData.SpiralScanData = SpiralData
    # procData.data3D = data3D
    
    data3D = 20*np.log10(data3D[:,:,::-1])
    nL = procOpts.normLow
    nH = procOpts.normHigh
    data3D = np.clip(data3D, nL, nH)
    data3D = np.uint16(65535*(data3D - nL)/(nH - nL))
    
    return procData

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
            for n in range(0, bscansPerFrame):
                if scanP.pattern == ScanPattern.bidirectional and (np.mod(n, 2) != 0):
                    idx1 = np.floor(trigsPerBscan*n) + procOpts.biDirTrigVolFix
                    idx2 = idx1 + scanP.lengthSteps
                    rawData.oct_data_mag = oct_data_tmp[idx2:idx1:-1, :]
                else:
                    idx1 = trigsPerBscan*n
                    idx2 = idx1 + trigsPerBscan - extraTrigs
                    
                    oct_data_mag = oct_data_tmp[idx1:idx2, :]
                if DebugLog.isLogging:
                    DebugLog.log("VolumeScan processData(): n= %d idx1= %d idx2= %d" % (n, idx1, idx2))
                
                img16b, img8b = makeBscanImage(oct_data_mag, scanParams, procOpts.zRes, procOpts.normLow, procOpts.normHigh, correctAspectRatio=True)
                if volDataIn is None:
                    volDataIn = VolumeData()
                    volDataIn.scanParams = scanParams
                    volDataIn.volumeImg = np.zeros((scanParams.widthSteps, img16b.shape[0], img16b.shape[1]), dtype=np.uint16)                
                    volDataIn.zPixSize = procOpts.zRes     # z pixel size in um
                    volDataIn.xPixSize = img16b.shape[1]/scanParams.length     # x pixel size in um
                    volDataIn.yPixSize = scanParams.widthSteps/scanParams.width     # y pixel size in um 
                
                
                DebugLog.log("VolumeScan processData(): frame= %d img16b max= %d min=%d " % (n + bscanPerFrame*frameNum, np.max(img16b), np.min(img16b)))
                volDataIn.volumeImg[n + bscanPerFrame*frameNum, :, :] = img16
        else:
            img16b, img8b = makeBscanImage(oct_data_mag, scanParams, procOpts.zRes, procOpts.normLow, procOpts.normHigh, correctAspectRatio=True)
            
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
    
def runVolScan(appObj):
    DebugLog.log("runVolScan")
    appObj.tabWidget.setCurrentIndex(2)
    
    if not appObj.oct_hw.IsOCTTestingMode():
        from DAQHardware import DAQHardware
        daq = DAQHardware()
        
    appObj.doneFlag = False
    appObj.isCollecting = True
    OCTtrigRate = appObj.oct_hw.GetTriggerRate()
    mirrorDriver = appObj.mirrorDriver
    rset = True
    
    chanNames = [mirrorDriver.X_daqChan, mirrorDriver.Y_daqChan]
    trigChan = mirrorDriver.trig_daqChan
    outputRate = mirrorDriver.DAQoutputRate

    # if in testing mode, load proper paramaeters instead of getting them from GUI
    if appObj.oct_hw.IsOCTTestingMode():
        testDataDir = os.path.join(appObj.basePath, 'exampledata\\VolumeScan')
        filePath = os.path.join(testDataDir, 'ScanParams.pickle')
        f = open(filePath, 'rb')
        scanParams = pickle.load(f)
        f.close()
    else:
        # get the scan paramters tht user has entred 
        scanParams = appObj.getScanParams()
        
    zROI = appObj.getZROI()
    scanDetails = None
    if scanParams.pattern == ScanPattern.spiral or scanParams.pattern == ScanPattern.wagonWheel:
        plotParam, scanDetails = setupScan(scanParams, mirrorDriver, zROI, OCTtrigRate, procOpts)
    
    bscansPerFrame = scanParams.volBscansPerFrame
    numFrames = scanParams.widthSteps // bscansPerFrame            
    framesPerScan = scanParams.widthSteps // bscansPerFrame        
    if scanParams.continuousScan:
        numFrames = np.inf
    elif scanParams.pattern == ScanPattern.spiral or scanParams.pattern == ScanPattern.wagonWheel:
        numFrames = 1
        framesPerScan = 1
    
    procOpts = ProcOpts()
    procOpts.normLow = appObj.normLow_spinBox.value()
    procOpts.normHigh = appObj.normHigh_spinBox.value()
    procOpts.zRes = appObj.octSetupInfo.zRes
    procOpts.biDirTrigVolFix = appObj.volBidirTrigFix_spinBox.value()
    biDirTrigVolAdj = appObj.volBidirTrigAdj_spinBox.value()
    
    saveOpts = appObj.getSaveOpts()
    isSaveDirInit = False

    try: 
        
        frameNum = 0
        scanNum = 0
        while not appObj.doneFlag and frameNum < numFrames:
            # reinitialize volume data on first frame
            if frameNum % framesPerScan == 0:
                volData = None
                
            if scanParams.pattern == ScanPattern.spiral or scanParams.pattern == ScanPattern.wagonWheel:
                mirrorOut = scanDetails.mirrOut
                startTrigOffset = 0
            else:
                mirrorOut = makeVolumeScanCommand(scanParams, frameNum, mirrorDriver, OCTtrigRate)
                startTrigOffset = int(np.round(OCTtrigRate*mirrorDriver.settleTime))
                
            # plot command to GUI 
            pl = appObj.plot_mirrorCmd
            npts = mirrorOut.shape[1]
            t = np.linspace(0, npts/outputRate, npts)
            pl.clear()
            pl.plot(t, mirrorOut[0, :], pen='b')  
            pl.plot(t, mirrorOut[1, :], pen='r')  
            labelStyle = appObj.xLblStyle
            pl.setLabel('bottom', 'Time', 's', **labelStyle)
            labelStyle = appObj.yLblStyle
            pl.setLabel('left', 'Output', 'V', **labelStyle)

            if not appObj.oct_hw.IsDAQTestingMode():
                # setup the analog output DAQ device
                daq.setupAnalogOutput(chanNames, trigChan, outputRate, mirrorOut.transpose())        
                daq.startAnalogOutput()
            
            # setup and grab the OCT data
            numTrigs = getNumTrigs(scanParams, scanDetails, OCTtrigRate, mirrorDriver)
            if appObj.oct_hw.IsOCTTestingMode():
                oct_data = OCTCommon.loadRawData(testDataDir, frameNum, dataType=0)
            else:
                err, oct_data = appObj.oct_hw.AcquireOCTDataFFT(numTrigs, zROI, startTrigOffset)
            
            # process the data
            oct_data_mag = np.abs(oct_data)
            volData = processData(oct_data_mag, scanParams, mirrorDriver, OCTtrigRate, procOpts, volData, frameNum)
                
            if scanParams.pattern == ScanPattern.spiral or scanParams.pattern == ScanPattern.wagonWheel:
                pass
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
            
            frameNum += 1
            appObj.acquisition_progressBar.setValue(round(100*frameNum/framesPerScan))
            if appObj.getSaveState():
                if not isSaveDirInit:
                    saveDir = OCTCommon.initSaveDir(saveOpts, 'Volume', scanParams)
                    isSaveDirInit = True
                if saveOpts.saveRaw:
                    OCTCommon.saveRawData(oct_data, saveDir, frameNum-1, dataType=0)
                
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