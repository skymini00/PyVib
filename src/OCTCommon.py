# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 14:31:38 2015

@author: OHNS
"""
import numpy as np
from enum import Enum
from DebugLog import *
import Dispersion
import re
import os
import datetime
import struct
import pickle

class StatusMsgType(Enum):
    BLANK = 0  # blank, default message which does nothing
    ERROR = 1
    DAQ_OUTPUT = 2
    ACQUISITION_COMPLETE = 3
    OCT_TRIG_RATE = 4
    KLIN = 5

class StatusMsgSource(Enum):
    UNKNOWN = 0  # unknown source
    COLLECTION = 1
    PROCESSING = 2
    LOGGING = 3

class ProcessMode(Enum):
    FPGA = 0
    SOFTWARE = 1
    GPU = 2
    
class StatusMsg:
    def __init__(self, msgSrc=StatusMsgSource.UNKNOWN, msgType=StatusMsgType.BLANK, param=None):
        self.msgSrc = msgSrc
        self.msgType = msgType
        self.param = param

class OCTSetupInfo:
    def __init__(self):
        self.setupNum = -1
        self.dispFilename = 'dispComp-initial.pickle'
        self.zRes = 8.37
        self.imgNorms = (50, 120)
        self.ZROI = (100, 100 + 4*120 - 1 )
        self.defaultScanQuickSet = ''
        self.defaultAudioQuickSet = ''
        self.mirrorConfigFile = 'GalvoMirror.txt'
        self.audioConfigFile = 'AudioHardware.txt'
        self.FPGAOptsFile = 'FPGA Opts.txt'
        self.Laser = 1       
        self.processMode = ProcessMode.SOFTWARE
        
    def encodeToString(self):
        s = ''
        s = s + 'OCTsetup=' +  repr(self.setupNum)
        s = s + '\nDispersionFilename=' +  self.dispFilename
        s = s + '\nzResolution=' +  repr(self.zRes)
        s = s + '\nimgNorms=' +  repr(self.imgNorms)
        s = s + '\nZROI=' +  repr(self.ZROI)
        s = s + '\nDefaultScanQuickSet=' + self.defaultScanQuickSet
        s = s + '\nDefaultAudioQuickSet=' + self.defaultAudioQuickSet
        s = s + '\nMirrorConfigFile=' + self.mirrorConfigFile
        s = s + '\nAudioConfigFile=' + self.mirrorConfigFile
        s = s + '\nFPGAOptsFile=' + self.FPGAOptsFile
        s = s + '\nLaser=' + self.Laser
        s = s + '\nProcessMode=' + self.processMode.name
        
        return s
        
    def getTriggerRate(self):
        trigRate=0
        if self.Laser == 1:
            trigRate = 200e3
        elif self.Laser==2:
            trigRate = 49.9598e3
        elif self.Laser==3:
            trigRate = 100e3
        elif self.Laser==4:     # test mode with a really low laser rate
            trigRate = 2000
        print('self.Laser, trigRate', self.Laser, trigRate)            
        return trigRate
        
#class DispCorr:
#   def __init__(self):
#        self.magWin = np.zeros(2048)
#        self.phaseCorr = np.zeros(2048)
        


# read in dispersino correction file, which contains magnitude and phase correction
# as complex numbers
def readDispersionFile(filePath):
    f = open(filePath, "rb")
    data = f.read()
    f.close()
    
    s = re.split('\.', filePath)
    
    if len(s) > 0 and s[1] == 'pickle':
        dispData = pickle.loads(data)
    else:
        dispData = Dispersion.DispersionData()
        
        data_len = struct.unpack("I", data[0:4])
        print("readDispersionFile() data_len = ", repr(data_len[0]))
        # need to read in 2*data_len because 
        fmt_str = "%dd" % (2*data_len[0])
        # unpack the data from binary 
        b = struct.unpack_from(fmt_str, data, 4)
        c = np.array(b)
        c = c.reshape((data_len[0], 2))
        d = c[:, 0] + 1j * c[:, 1]
        dispData.magWin = np.abs(d)
        dispData.phaseCorr = np.angle(d)
        
    return dispData


from ast import literal_eval as make_tuple
def readOCTSetupInfo(filepath):    
    setupInfo = OCTSetupInfo()

    f = open(filepath, "r")
    txt = f.read()
    f.close()
    
    #print("txt=" + txt)
    # break data up into lines
    lines = re.split('\n', txt)
    #print(lines)
    for l in lines:
        try:
            # remove comments
            #print("l= " + l)
            lnc = re.split('#', l)
            #print("lnc = " + repr(lnc[0]))   # remove commens
            s = re.split('=', lnc[0])
            #print("s = " + repr(s))
            fld = s[0].rstrip()
            val = s[1]
            # print("fld: %s val: %s" % (fld, val))
            if(fld == 'OCTsetup'):
                setupInfo.setupNum = int(val)
            elif(fld == 'DispersionFilename'):
                setupInfo.dispFilename = val
            elif(fld == 'zResolution'):            
                setupInfo.zRes = float(val)
            elif(fld == 'imgNorms'):            
                setupInfo.imgNorms = make_tuple(val)
            elif(fld == 'ZROI'):            
                setupInfo.ZROI = make_tuple(val)
            elif(fld == 'DefaultScanQuickSet'):            
                setupInfo.defaultScanQuickSet = val
            elif(fld == 'DefaultAudioQuickSet'):            
                setupInfo.defaultAudioQuickSet = val
            elif(fld == 'MirrorConfigFile'):            
                setupInfo.mirrorConfigFile = val
            elif(fld == 'AudioConfigFile'):            
                setupInfo.audioConfigFile = val
            elif(fld == 'FPGAOptsFile'):
                setupInfo.FPGAOptsFile = val
            elif(fld == 'Laser'):
                setupInfo.Laser = int(val)
            elif(fld == 'ProcessMode'):
                setupInfo.processMode = ProcessMode[val.upper()]
                    
        except Exception as ex:
            pass
            
    return setupInfo
    
def writeOCTsetupInfo(basepath, octSetupInfo):    
    filepath = os.path.join(basepath, 'oct.txt')
    f = open(filepath, "w")
    setupInfoStr = octSetupInfo.encodeToString()
    print("writeOCTsetupInfo(): setupInfoStr = %s" % setupInfoStr)
    f.write(setupInfoStr)
    f.close()
    
class blankRecord:
    pass

class SaveDirNameScheme(Enum):
    TIMESTAMP_PROTOCOL_SUBJECT = 0
    USE_BASE_DIR = 1

class SaveOpts:
    def __init__(self):
        self.notes = ''
        self.subject = ''
        self.saveBaseDir = ''
        
        self.saveRaw = False
        self.saveImages = True
        self.saveMscanTuningCurves = True
        self.saveMscanFFT = False
        self.saveMscanTDPhase = True
        
        # whether to save most recent frame or to save all frames
        self.saveOnlyMostRecentFrame = True
        
        self.dirNameScheme = SaveDirNameScheme.TIMESTAMP_PROTOCOL_SUBJECT
        
    def __repr__(self):
        s = 'saveBaseDir= %s subject= %s notes= %s dirNameScheme= %s' % (self.notes, self.subject, self.saveBaseDir, repr(self.dirNameScheme))
        s = s + '\nsaveRaw = %s saveImages= %s ' % (self.saveRaw, self.saveImages)
        s = s + '\nsave mscan tuningcurves= %s FFT= %s TD phase= %s' % (self.saveMscanTuningCurves, self.saveMscanFFT, self.saveMscanTDPhase)
        
        return s


def initSaveDir(saveOpts, protocolName, scanParams=None, audioParams=None, mirrorDriver=None, OCTtrigRate=None, processMode=None, plotParam=None, dispData=None):
    baseDir = saveOpts.saveBaseDir
        
    nameScheme = saveOpts.dirNameScheme
    d = datetime.datetime.now()
    timeStr = d.strftime('%Y-%m-%d %H_%M_%S')
    if nameScheme == SaveDirNameScheme.TIMESTAMP_PROTOCOL_SUBJECT:
        saveDir = timeStr + ' ' + protocolName
        if saveOpts.subject != '':
            saveDir = saveDir + ' ' + saveOpts.subject
    elif nameScheme == SaveDirNameScheme.USE_BASE_DIR:
        saveDir = ''
    
    if saveDir != '':
        saveDir = os.path.join(baseDir, saveDir)
    else: 
        saveDir = baseDir
        
    if not os.path.exists(saveDir):   # create directory if it does not exist
        os.makedirs(saveDir)

        
    # make a text file with general information and notes
    filepath = os.path.join(saveDir, 'miscellaneous information.txt')
    f = open(filepath, 'w')
    str1='processMode= ' + repr(processMode) + '\n'
    f.write(str1)
    str1='OCTtrigRate=' + repr(OCTtrigRate) + '\n'
    f.write(str1)    
    if plotParam is not None:      # write scan parameters to text file    
        str1='plotParam.xPixelSize= ' + repr(plotParam.xPixelSize) + '\n'
        f.write(str1)    
        str1='plotParam.yPixelSize= ' + repr(plotParam.yPixelSize) + '\n'
        f.write(str1)    
        str1='plotParam.zPixelSize= ' + repr(plotParam.zPixelSize) + '\n'
        f.write(str1)            
    if saveOpts.notes != '':     # write notes     
        f.write('Notes: \n')
        f.write(saveOpts.notes)               
    f.close()
    if scanParams is not None:      # write scan parameters to text file
        scanParamsStr = repr(scanParams)
        fileName = 'ScanParams.txt'
        filePath = os.path.join(saveDir, fileName)
        f = open(filePath, 'w')
        f.write(scanParamsStr)
        f.close()
        
        fileName = 'ScanParams.pickle'
        filePath = os.path.join(saveDir, fileName)
        f = open(filePath, 'wb')
        pickle.dump(scanParams, f)
        f.close()
    if audioParams is not None:     # write audio parameters to text file
        audioParamsStr = repr(audioParams)
        fileName = 'AudioParams.txt'
        filePath = os.path.join(saveDir, fileName)
        f = open(filePath, 'w')
        f.write(audioParamsStr)
        f.close()
        
        fileName = 'AudioParams.pickle'
        filePath = os.path.join(saveDir, fileName)
        f = open(filePath, 'wb')
        pickle.dump(audioParams, f)
        f.close()
    if mirrorDriver is not None:      # write scan parameters to text file
        mirrorDriverStr = repr(mirrorDriver)
        fileName = 'mirrorDriver.txt'
        filePath = os.path.join(saveDir, fileName)
        f = open(filePath, 'w')
        f.write(mirrorDriverStr)
        f.close()
        
        fileName = 'mirrorDriver.pickle'
        filePath = os.path.join(saveDir, fileName)
        f = open(filePath, 'wb')
        pickle.dump(mirrorDriver, f)
        f.close()
    if dispData is not None:      # write scan parameters to text file       
        fileName = 'dispData.pickle'
        filePath = os.path.join(saveDir, fileName)
        f = open(filePath, 'wb')
        pickle.dump(dispData, f)
        f.close()

    return saveDir




dataTypeFileNames = [ 'OCT complex FFT ',  'OCT interp PD ', 'OCT raw PD/MZI ', 'Mic raw ']
# save teh raw data
# dataType = 0 for complex FFT, 1 for interpolated photodioe, 2 for raw PD/MZI    
def saveRawData(octData, saveDir, frameNum,  dataType=0, trialNum=None):
    fileName = dataTypeFileNames[dataType]
        
    fileName = fileName + 'frame ' + repr(frameNum)
    if trialNum is not None:
        fileName + ' trial ' + repr(trialNum)
    
    filePath = os.path.join(saveDir, fileName)
    f = open(filePath, 'wb')
    np.save(f, octData)
    f.close()

def loadRawData(saveDir, frameNum, dataType=0, trialNum=None):
    fileName = dataTypeFileNames[dataType]
    fileName = fileName + 'frame ' + repr(frameNum)
    if trialNum is not None:
        fileName + ' trial ' + repr(trialNum)
        
    filePath = os.path.join(saveDir, fileName)
    f = open(filePath, 'rb')
    oct_data = np.load(f)
    f.close()
    return oct_data

def saveRawDataSoftwareProcessing(ch0_data, ch1_data, saveDir, frameNum): 
    outfile = os.path.join(saveDir, 'RawData_%.5d.npz' % (frameNum))
    np.savez_compressed(outfile, ch0_data=ch0_data, ch1_data=ch1_data)


# tests    
if __name__ == "__main__":   
    filepath = 'C:\\PyOCT\\oct.txt'
    setupInfo = readOCTSetupInfo(filepath)
    setupInfo.zres = 6
    basepath = 'C:\\PyOCT'
    writeOCTsetupInfo(basepath, setupInfo)
    
    print(setupInfo.__dict__)
