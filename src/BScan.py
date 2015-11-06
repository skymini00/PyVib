# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 12:22:44 2015

@author: OHNS
"""


import numpy as np

import OCTCommon
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
import pickle
import time
import JSOraw

def correctImageAspectRatio(imgData, xRes, zRes, interpMode=PIL.Image.BILINEAR):
    if not np.isfinite(xRes) or not np.isfinite(zRes):
        return None
        
    imgData = np.require(imgData, np.uint8, 'C')
    s = imgData.shape
    img = PIL.Image.fromarray(imgData, 'P')     # 'I' for 32-bit images, 'P' for 8-bit
    # correct aspect ratio
    
    sizeX = s[1]
    sizeY = s[0]
    newSizeX = int(np.round(sizeX*xRes/zRes))
    newSizeX = min(newSizeX, 2000)    # guard against images that are too big
    newSizeX = max(2, newSizeX)
    if DebugLog.isLogging:
        DebugLog.log("correctImageAspectRatio() xRes= %g zRes= %g sizeX= %d sizeY=%d newSizeX = %d" % (xRes, zRes, sizeX, sizeY, newSizeX))
    img = img.resize((newSizeX, sizeY), interpMode) 
    imgData = np.array(img)
    
    imgData = np.require(imgData, np.uint8, 'C')
    
    return imgData


def makeBscanImage(oct_data, scanParams, zRes, normLow, normHigh, correctAspectRatio=True):
    # numTrigs = self.scanParams.lengthSteps
    rawdata = oct_data
    numTrigs = rawdata.shape[0]
    
    # for a slwo scan, reshape and exclude unwanted triggers
    if scanParams.pattern == ScanPattern.rasterSlow: 
        # rawdata = trigsPerPt = np.ceil(galv.settleTime * OCTtriggerRate)
        shp = rawdata.shape
        
        trigsPerStep = numTrigs // scanParams.lengthSteps
        DebugLog.log("makeBscanImg: rawdata.shape= %s trigsPerStep= %d" %  (repr(rawdata.shape), trigsPerStep))
        # rawdata = np.reshape(rawdata, (trigsPerStep, scanParams.lengthSteps, shp[1], shp[2]))
        
        i1 = trigsPerStep - 1
        i2 = scanParams.lengthSteps * trigsPerStep
        rng = list(range(i1, i2, trigsPerStep))
        rawdata = rawdata[rng, :]
            
    mag = np.abs(rawdata)
    mag = np.clip(mag, 1, np.inf)
    mag = 20*np.log10(mag)
        
    DebugLog.log("makeBscanImg: (before averaging) mag.shape= " + repr(mag.shape))          
    
    # average all trials
    if len(mag.shape) > 2:
        mag = np.mean(mag, 2)
    
    DebugLog.log("makeBscanImg: (after averaging) mag.shape= " + repr(mag.shape))  
    mag = mag.transpose()
    
    #DebugLog.log("processBScanData: mag min= %g max= %g" % (np.min(mag), np.max(mag)))
     # map the data to 0 ... 1
    nL = normLow
    nH = normHigh
    nRng = nH - nL
    
    DebugLog.log("makeBscanImg: nL= %g nH=%g nRng= %g" % (nL, nH, nRng))
        
    # remap range to 0...1
    mag = (mag - nL)/nRng  

    # remap range to 0 ... to 2^16 - 1
    mag16b = mag*65535
    mag16b = np.clip(mag16b, 0, 65535)
    
    DebugLog.log("makeBscanImg:  mag16b min= %g max= %g" % (np.min(mag16b), np.max(mag16b)))
    
    # produce image            
    mag16b = np.require(mag16b, 'uint32')

    s = mag16b.shape
    
    sizeX = s[1]
    sizeY = s[0]
    
    sl = scanParams.length
    xRes = sl*1e3/sizeX
    
    DebugLog.log("makeBscanImg: sizeX = %d sizeY= %d xRes= %g" % (sizeX, sizeY, xRes))
    
    img16b = PIL.Image.fromarray(mag16b, 'I')     # 'I' for 32-bit images, 'P' for 8-bit
    # correct aspect ratio
    if(correctAspectRatio):
        sizeX = int(np.round(sizeX*xRes/zRes))
        sizeX = min(sizeX, 4000)    # guard against images that are too big
        sizeX = max(2, sizeX)
        DebugLog.log("makeBscanImg:  new sizeX = %d" % (sizeX))
        img16b = img16b.resize((sizeX, sizeY), PIL.Image.BILINEAR) 

    imgdata_16b = np.array(img16b)
    imgdata_8b = np.floor(255*imgdata_16b/65535)

    imgdata_8b = np.require(imgdata_8b, np.uint8, 'C')
    imgdata_16b = np.require(imgdata_16b, np.uint16, 'C')
    DebugLog.log("ScanningProtocol.processData:  img8b min= %g max= %g" % (np.min(imgdata_8b), np.max(imgdata_8b)))
    DebugLog.log("ScanningProtocol.processData:  img16b min= %g max= %g" % (np.min(imgdata_16b), np.max(imgdata_16b)))

    return imgdata_16b, imgdata_8b
    
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
def makeBscanCommand(scanParams, galvoMirror, OCTtriggerRate, volWidthStep=-1):
    cmd_x = None
    cmd_y = None
    galv = galvoMirror
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
        
        (x1, y1, x2, y2) = makeScanXYcoords(scanParams, galvoMirror)
        
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
    
    
# makeBscanMirrorOutput
    
#   generates the mirror output function  and returns the expetected number of OCT triggers
#   parameters    
#      scanParams:   ScanParams object
#      widthStep:  interger represenging the width position of in the grid, should b e 0...scanParams.widthStpes
#      lenStep:  interger represenging the length position of in the grid, should b e 0...scanParams.lenghStpes
#       galvoMirror: GalvoMirror object        
#   returns (tuple)
#       mirr_out as a tuple (cmd_x, cmd_y), both 1D numpy arrays 
#       numTrigs
def makeBscanMirrorOutput(scanParams, galvoMirror, OCTtrigRate):
    mirrOut = makeBscanCommand(scanParams, galvoMirror, OCTtrigRate)
    numTrigs = scanParams.lengthSteps
    if scanParams.pattern == ScanPattern.rasterSlow:
        numTrigs = np.ceil(scanParams.lengthSteps * galvoMirror.settleTime * OCTtrigRate)
    
    return (mirrOut, numTrigs)
    
def saveImage(bscanImg, saveDir, saveOpts, frameNum):
    bscanImg = np.require(bscanImg, np.uint16)
    fileName = 'Bscan'
    if not saveOpts.saveOnlyMostRecentFrame:
        fileName = fileName + ("_%4d" % frameNum)
    
    fileName = fileName + '.tiff'
    filePath = os.path.join(saveDir, fileName)
    
    tifffile.imsave(filePath, bscanImg)    
    
class BScanRawData():
    def __init__(self, oct_data=None, frameNum=-1):
        self.oct_data = oct_data
        self.frameNum = frameNum
        
def BscanCollectFunction(oct_hw, frameNum, extraArgs)    :
    scanParams = extraArgs[0]
    mirrorDriver = extraArgs[1]
    zROI = extraArgs[2]
    dispCorr = extraArgs[3]
    testDataDir = extraArgs[4]
    trigRate = oct_hw.GetTriggerRate()
    # testDataDir = os.path.join(basePath, 'exampledata\\Bscan')
    
    # generate the mirrrour output function
    (mirrorOutput, numTrigs) = makeBscanMirrorOutput(scanParams, mirrorDriver, trigRate)
    numTrigs = int(np.round(numTrigs))  # ensure numTrigs is an integer
    
    # setup the analog output DAQ device
    chanNames = [mirrorDriver.X_daqChan, mirrorDriver.Y_daqChan]
    trigChan = mirrorDriver.trig_daqChan
    outputRate = mirrorDriver.DAQoutputRate
    
    if not oct_hw.IsDAQTestingMode():
        from DAQHardware import DAQHardware
        daq = DAQHardware()
        daq.setupAnalogOutput(chanNames, trigChan, outputRate, mirrorOutput.transpose())        
        daq.startAnalogOutput()
    
    # setup and grab the OCT data
    startTrigOffset = int(np.round(trigRate*mirrorDriver.settleTime))
    
    if oct_hw.IsOCTTestingMode():
        oct_data = OCTCommon.loadRawData(testDataDir, frameNum % 15, dataType=0)
        DebugLog.log('BscanCollectFunction floading test data frameNum= ' + repr(frameNum))
    else:
        err, oct_data = oct_hw.AcquireOCTDataFFT(numTrigs, zROI, startTrigOffset, dispCorr)
        DebugLog.log('BscanCollectFunction AcquireOCTDataFFT() err= ' + repr(err))
    if not oct_hw.IsDAQTestingMode():
        daq.waitDoneOutput()
        daq.stopAnalogOutput()
        daq.clearAnalogOutput()
        
    rawData = BScanRawData(oct_data, frameNum)
    return rawData, mirrorOutput

def processAndDisplayBscanData(appObj, oct_data, scanParams, rset, zROI):
    normLow = appObj.normLow_spinBox.value()
    normHigh = appObj.normHigh_spinBox.value()
    zRes = appObj.octSetupInfo.zRes
    
    img16b, img8b = makeBscanImage(oct_data, scanParams, zRes, normLow, normHigh, correctAspectRatio=True)
    appObj.imgView.setImage(img8b.transpose())    
    lut = appObj.HOT_LUT
    appObj.imgView.getImageItem().setLookupTable(lut) 
    appObj.bscan_img_gv.setImage(img8b, ROIImageGraphicsView.COLORMAP_HOT, rset)
    rset = False
    appObj.imgDataScanParams = copy.copy(scanParams)
    appObj.imgdata_8b = img8b
    appObj.imgdata_zROI = zROI
    
    return img16b

def saveBScanData(appObj, oct_data, img16b, frameNum, scanParams, saveOpts, isSaveDirInit):
    if appObj.getSaveState():
        if not isSaveDirInit:
            saveDir = OCTCommon.initSaveDir(saveOpts, 'BScan', scanParams=scanParams)
            isSaveDirInit = True
    
        saveImage(img16b, saveDir, saveOpts, frameNum)
        if saveOpts.saveRaw:
            OCTCommon.saveRawData(oct_data, saveDir, frameNum, dataType=0)
            
    return isSaveDirInit

def loadScanParams(testDataDir):
    filePath = os.path.join(testDataDir, 'ScanParams.pickle')
    f = open(filePath, 'rb')
    scanParams = pickle.load(f)
    f.close()
    
    return scanParams

def handleStatusMessage(statusMsg):
    err = False
    if statusMsg.msgSrc == OCTCommon.StatusMsgSource.COLLECTION:
        if statusMsg.msgType == OCTCommon.StatusMsgType.ERROR:
            err = True
        elif statusMsg.msgType == OCTCommon.StatusMsgType.DAQ_OUTPUT:
            pass
        
    return err
        
def runBScanMultiProcess(appObj, testDataDir):
    DebugLog.log("runBScanMultiProcess")
    oct_hw = appObj.oct_hw
    trigRate = appObj.oct_hw.GetTriggerRate()
    mirrorDriver = appObj.mirrorDriver
    rset = True
    saveOpts = appObj.getSaveOpts()
    isSaveDirInit = False
    try: 
        oct_hw.NewAcquisition()
        oct_hw.SetSendExtraInfo(False)   # do not send mirror output
        
        if oct_hw.IsOCTTestingMode():
            scanParams = loadScanParams(testDataDir)

        # tell collection background process to use this function to acquire data
        oct_hw.SetAcqFunction(BscanCollectFunction)
        startAcq = True
        
        # clear the status message queue
        statusMsg = oct_hw.GetStatus()
        while statusMsg is not None:
            err = handleStatusMessage(statusMsg)
            statusMsg = oct_hw.GetStatus()
        
        while not appObj.doneFlag:
            DebugLog.log("runBScanMultiProcess: acquiring")
            # get parameters from GUI
            if not oct_hw.IsOCTTestingMode():
                # get the scan paramters tht user has entred 
                scanParams = appObj.getScanParams()
            zROI = appObj.getZROI()
            dispCorr = appObj.dispCorr_pushButton.isChecked()  # whehter or not to do dispersion correction
            
            # update parameters in background process
            extraArgs = [scanParams, mirrorDriver, zROI, dispCorr, testDataDir]
            oct_hw.SetAcqFunctionArgs(extraArgs)

            # start the acquisitio on first loop iteration
            # we don't just do this outside the loop because we haven't loaded function args yet
            if startAcq:  
                oct_hw.StartAcquisition() 
                startAcq = False
            
            rawData = oct_hw.GetData()
            if rawData is not None and isinstance(rawData, BScanRawData):
                # DebugLog.log("runBScanMultiProcess: got data dict= " + rawData.__dict__)
                oct_data = rawData.oct_data
                frameNum = rawData.frameNum
                img16b = processAndDisplayBscanData(appObj, oct_data, scanParams, rset, zROI)
                isSaveDirInit = saveBScanData(appObj, oct_data, img16b, frameNum, scanParams, saveOpts, isSaveDirInit)
                rset = False
            else:
                DebugLog.log("runBScanMultiProcess: data is None or is not BScanRawData")
                
            statusMsg = oct_hw.GetStatus()
            while statusMsg is not None:
                DebugLog.log("runBScanMultiProcess: got status message type=" + repr(statusMsg.msgType))
                err = handleStatusMessage(statusMsg)
                if err:
                    appObj.doneFlag = True  # if error occured, stop pcollecting
                statusMsg = oct_hw.GetStatus()
            
            # check for GUI events, particularly the "done" flag
            QtGui.QApplication.processEvents() 
            
            time.sleep(0.005)  # sleep fr 5 ms to free up CPU
    except Exception as ex:
        # raise ex
        traceback.print_exc(file=sys.stdout)
        QtGui.QMessageBox.critical (appObj, "Error", "Error during scan. Check command line output for details")
    finally:
        appObj.doneFlag = True
        oct_hw.PauseAcquisition()
        appObj.isCollecting = False
        QtGui.QApplication.processEvents() # check for GUI events
        appObj.finishCollection()

"""
    starts running a Bscan
    
    appObj is of OCTWindowClass defined in PyOCT.py
"""
def runBScan(appObj):
    DebugLog.log("runBScan")
    appObj.tabWidget.setCurrentIndex(0)
    appObj.doneFlag = False
    appObj.isCollecting = True
    testDataDir = os.path.join(appObj.basePath, 'exampledata', 'Bscan')
    
    if(appObj.multiProcess):
        runBScanMultiProcess(appObj, testDataDir)
        return

    if not appObj.oct_hw.IsOCTTestingMode():
        from DAQHardware import DAQHardware
        daq = DAQHardware()
        
    trigRate = appObj.oct_hw.GetTriggerRate()
    mirrorDriver = appObj.mirrorDriver
    rset = True
    saveOpts = appObj.getSaveOpts()
    isSaveDirInit = False
    frameNum = 0
    # softwareProcess = softwareProcessing_pushButton.isChecked()
    
    try: 
        if appObj.oct_hw.IsOCTTestingMode():
            scanParams = loadScanParams(testDataDir)
                    
        while not appObj.doneFlag:
            startTime = time.time()
            processMode = OCTCommon.ProcessMode(appObj.processMode_comboBox.currentIndex())

            if not appObj.oct_hw.IsOCTTestingMode():
                # get the scan paramters tht user has entred 
                scanParams = appObj.getScanParams()
            
            # generate the mirrrour output function
            (mirrorOutput, numTrigs) = makeBscanMirrorOutput(scanParams, mirrorDriver, trigRate)
            numTrigs = int(np.round(numTrigs))  # ensure numTrigs is an integer
            
            # setup the analog output DAQ device
            chanNames = [mirrorDriver.X_daqChan, mirrorDriver.Y_daqChan]
            trigChan = mirrorDriver.trig_daqChan
            outputRate = mirrorDriver.DAQoutputRate
            
            if not appObj.oct_hw.IsDAQTestingMode():
                daq.setupAnalogOutput(chanNames, trigChan, outputRate, mirrorOutput.transpose())        
                daq.startAnalogOutput()
            
            # setup and grab the OCT data
            zROI = appObj.getZROI()
            startTrigOffset = int(np.round(trigRate*mirrorDriver.settleTime))
            dispCorr = appObj.dispCorr_pushButton.isChecked() # whehter or not to do dispersion correction
            if appObj.oct_hw.IsOCTTestingMode():
                oct_data = OCTCommon.loadRawData(testDataDir, frameNum % 15, dataType=0)
            else:
                if processMode == OCTCommon.ProcessMode.FPGA:
                    err, oct_data = appObj.oct_hw.AcquireOCTDataFFT(numTrigs, zROI, startTrigOffset, dispCorr)
                elif processMode == OCTCommon.ProcessMode.SOFTWARE:
                    fpgaOpts = appObj.oct_hw.fpgaOpts
                    err, pd_data, mzi_data = appObj.oct_hw.AcquireOCTDataRaw(numTrigs, startTrigOffset=startTrigOffset)
                    # calcOCTDataFFT(pdData, mziData, MZI_PD_shift, klinROI_idx, numklinpts, klin, dispCorr_mag, dispCorr_ph):
                    
                    klinROI = [fpgaOpts.klinRoiBegin, fpgaOpts.klinRoiEnd]
                    oct_data, appObj.klin = JSOraw.calcOCTDataFFT(pd_data, mzi_data, zROI, appObj) 
#                    fpgaOpts.Ch0Shift*2, klinROI, fpgaOpts.numKlinPts, appObj.klin, appObj.dispCorr.magWin, appObj.dispCorr.phaseCorr, zROI)
                else:
                    QtGui.QMessageBox.critical (appObj, "Error", "Unsuppoted processing mode for current hardware")
            
            img16b = processAndDisplayBscanData(appObj, oct_data, scanParams, rset, zROI)
            rset = False  # turn off reset after first loop
            
            if not appObj.oct_hw.IsDAQTestingMode():
                daq.waitDoneOutput()
                daq.stopAnalogOutput()
                daq.clearAnalogOutput()

            isSaveDirInit = saveBScanData(appObj, oct_data, img16b, frameNum, scanParams, saveOpts, isSaveDirInit)

            frameNum += 1
            # check for GUI events, particularly the "done" flag
            QtGui.QApplication.processEvents() 
            updateTime = (time.time() - startTime)*1000
            appObj.updateRate_label.setText("%0.1f ms" % updateTime)
    except Exception as ex:
        # raise ex
        traceback.print_exc(file=sys.stdout)
        QtGui.QMessageBox.critical (appObj, "Error", "Error during scan. Check command line output for details")
    finally:
        appObj.isCollecting = False
        QtGui.QApplication.processEvents() # check for GUI events
        appObj.finishCollection()
    