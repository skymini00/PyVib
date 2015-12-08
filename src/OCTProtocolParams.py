# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 21:32:28 2015

@author: OHNS
"""

import numpy as np

class ScanParams:
    def __init__(self):
        self.length = 3              # scan length in millimeters
        self.lengthOffset = 0        # scan length offset in millimeers
        self.lengthSteps = int(1200) # number of length steps (# of length  steps (size of image, without downsampling))   
        self.width = 3               # scan width in millimeters
        self.widthOffset = 0         # width offet in millimeters
        self.widthSteps = int(200)   # number of width steps (size of image)
        self.rotation_Z = 0           # rotation in degrees
        self.pattern = ScanPattern.rasterFast   # pattern of the scan
        self.ptsList = [[]]          # a list of points, fist dimension is pt #, second is the x/y 
        self.boxROIMaskX = [[]]       # a 2D grid of True/False values, used for scans that are based off of a restriction region of a volume
        self.numAverages = 1         # number of times to repeat scan
        self.downsample = 0          # used for FPGA downsampling  
        
        # the following parameters are used for scans which are reprojected along an arbitrary plane
        # although the angle of the scan is always normal to the XY plane,  the Z ROI can depend on the frame, so it is mimicing a 
        self.frameDepZROI = False    # whether this is a reprojected scan, i.e., the ZROI depends on frame,
        self.frameDepZROIArray = [[]]   # array that maps the frame # to a tuple of zROI (begin, end)
        self.rotationXZ = 0           # the rotation of th plane with respect to the XZ plane
        self.rotationYZ = 0           # the rotation of the plane with respect to the YZ plane
        
        self.volBscansPerFrame = 1    # how many bscans to acquire in one frame for a volume scan
        self.xskew = 1.0              # multiplier for X mirror command to account for differences in volts per mm during fast scan
        self.phaseAdjust = 0.0
        self.angularScanFreq = 280
        self.volScanFreq = 5
        self.continuousScan = False
        
    def __repr__(self):
        s = ''
        s = s + '\n\tlength= %0.3f lengthOffset= %0.3f lengthSteps= %d ' % (self.length, self.lengthOffset, self.lengthSteps)
        s = s + '\n\twidth= %0.3f widthOffset= %0.3f widthSteps= %d ' % (self.width, self.widthOffset, self.widthSteps)
        s = s + '\n\trotation_Z= %0.1f pattern= %s numAverages= %d downsample= %d' % (self.rotation_Z, self.pattern, self.numAverages, self.downsample)
        s = s + '\n\trotationXZ= %0.2f rotationXZ= %0.2f' % (self.rotationXZ, self.rotationYZ)
        s = s + '\n\tbscansPerFrame= %d xskew= %0.3f phaseAdjust= %0.3f' % (self.volBscansPerFrame, self.xskew, self.phaseAdjust)
        s = s + '\n\tcontinuousScan= %s angularScanFreq= %0.3f volScanFreq= %0.3f' % (repr(self.continuousScan), self.angularScanFreq, self.volScanFreq)
        ptsListStr = ''
        for n in range(0, len(self.ptsList)):
            ptsListStr = ptsListStr + repr(self.ptsList[n])
        s = s + '\nptList= %s' % ptsListStr
        return s
        
    def sameDim(self, scanP):
        return self.length == scanP.length and self.width == scanP.width
        
from enum import Enum

class ScanPattern(Enum):
    rasterFast = 1               # a scan based on a fast linear sweep in scan direction and sinuisoidal flyback
    rasterSlow = 2               # A scan based on a stairstep sweep and sinuoisoidal flyback
    bidirectional = 3             # used for volumes, a linear scan in both directions with no flyback
    spiral = 4
    wagonWheel = 5
    zigZag=6
    boxROIMask = 100              # used for mscan, a scan based off of a X/Y region like a volume, but with a mask indicatingwheter
    ptsList = 101                 # used for mscan, a discrete list o points

class AudioStimType(Enum):
    PURE_TONE = 1
    TWO_TONE_DP = 2
    NOISE = 3
    CUSTOM = 4     

class Speaker(Enum):
    LEFT = 1
    RIGHT = 2
    BOTH = 3
    
class FreqSpacing(Enum):
    LINEAR = 1
    LOGARITHMIC = 2
    CUSTOM = 3
    
class AudioOutputParams:
    def __init__(self):
        self.freq = np.array([ [1e3], [1e3] ])     # array frequencies in Hz, first index is speaker num, second  is frequency index,
        self.amp = [ 80 ]           # array amplitudes in decibels
        self.stimDuration = 50      # duration of stim in milliseconds
        self.trialDuration = 50     # duration of entire trialin milliseconds 
        self.speakerSel = Speaker.LEFT    # which speaker or speakers to use
        self.stimType = AudioStimType.PURE_TONE
        self.stimOffset = 0         # offset ofthe stim in milliseconds
        self.stimEnvelope = 1       # envelope of the stim stim milliseconds
        self.useAmpDepTrialDuration = False
        self.ampDepTrialDurationArray = []
        self.numTrials = 1
        self.downsample = 0         # used for FPGA downsampling  
        self.freqSpacing = FreqSpacing.LINEAR
        self.customSoundir = ''
        
    def getTrialDuration(self, amp):
        if(self.useAmpDepTrialDuration):
            aIdx = self.amp.index(amp)
            self.ampDepTrialDurationArray[aIdx]    
        else:
            return self.trialDuration
    
    def getNumSpeakers(self):
        numSpk = 1
        if self.speakerSel == Speaker.BOTH:
            numSpk = 2
       
        return numSpk
       
    def getNumFrequencies(self):
        return self.freq.shape[1]
        
    def __repr__(self):
        freqStr = ''
        for n in range(0, len(self.freq[0, :])):
            freqStr = freqStr + ' %0.3f' % self.freq[0, n]
        s = 'freq= %s' % freqStr 
        
        ampStr = ''
        for n in range(0, len(self.amp)):
            ampStr = ampStr + ' %0.3f' % self.amp[n]
        s = s + '\namp= %s' % ampStr 
        s = s + '\nstimDuration= %0.3f trialDuration= %0.3f ' % (self.stimDuration, self.trialDuration)
        s = s + '\nspeakerSel= %s stimType= %s ' % (repr(self.speakerSel), repr(self.stimType))
        s = s + '\nstimOffset= %0.3f stimEnvelope = %0.3f' % (self.stimOffset, self.stimEnvelope)
        s = s + '\nuseAmpDepTrialDuration= %s ' % (repr(self.useAmpDepTrialDuration))
        s = s + '\nampDepTrialDurationArray= %s ' % (repr(self.ampDepTrialDurationArray))
        s = s + '\nnumTrials= %d downsample= %d' % (self.numTrials, self.downsample)
        
        return s
            

class DisplayLayoutType(Enum):
    THREE = 1
    FOUR_2_2 = 2
    FOUR_1_3 = 3
    FIVE_3_2 = 4
    FIVE_2_3 = 5
    SIX_3_3 = 6
    
class DisplaySetup:
    def __init__(self):
        self.layout = DisplayLayoutType.THREE
        self.graphArray = [0, 1, 2]
        
class AdvancedSetup:
    def __init__(self):
        self.fpgaDownsampling = 5
        
if __name__ == "__main__":
    audioParams = AudioOutputParams()
    freqStart = 4
    freqEnd = 10
    freqSteps = 7
    freq = np.linspace(freqStart, freqEnd, freqSteps)
    audioParams.freq = np.tile(freq, (2, 1))
    audioParams.amp = np.linspace(10, 80, 9)
    print(repr(audioParams))