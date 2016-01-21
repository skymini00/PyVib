# -*- coding: utf-8 -*-
"""
Created on Thu Oct  1 09:47:06 2015

@author: OHNS
"""

import OCTCommon
from ctypes import *
import numpy.ctypeslib as np_ctypes

from enum import Enum
import numpy as np
import matplotlib.pyplot as plt
import re   # regular expressions
from DebugLog import DebugLog
import multiprocessing as mproc
import traceback
import os
import sys
import psutil
import queue # for Queue.ull exception
import copy
import time

class FPGAOpts_t(Structure):
    _fields_ = [("Ch0Shift", c_uint16), 
                ("numKlinPts", c_uint16),
                ("StartTriggerOfffset", c_uint16),
                ("PDScaling", c_int16),
                ("BGSubtract", c_uint8),   # boolean
                ("DispCorr", c_uint8),      # boolean
                ("InterpData", c_uint8),    # boolean
                ("FFT", c_uint8),           # boolean
                ("procRoiBegin", c_uint16),
                ("procRoiEnd", c_uint16),
                ("InterpFilter", c_uint8),  # boolean
                ("InterpDownsample", c_uint8),  # boolean
                ("klinRoiBegin", c_uint16),
                ("klinRoiEnd", c_uint16),
                ("SampleOffset", c_int16),
                ("MZIScaleFactor", c_int16),
                ("MZIHPFilter", c_uint8),       # boolean
                ("DownsampleFactor", c_int16),
                ("SamplesPerTrig", c_int16),
                ("Ch0OffsetTweak", c_int16),
                ("SyncTrig", c_uint8),
                ("ProcessMZI", c_uint8),
                ("FFT16b", c_uint8),
                ("Polar", c_uint8),
                ("MagOnly", c_uint8),
                ("FFTReImgScale", c_int8),
                ("InterpDSAvg", c_uint8),
                ("OnDutyCycleTrigLast", c_uint16),
                ("OffDutyCycleTrigLast", c_uint16),
                ("UseDutyCycle", c_uint8) ]

    def __init__(self):
        self.Ch0Shift = c_uint16(20)
        self.numKlinPts = c_uint16(4000)
        self.StartTriggerOfffset = c_uint16(0)
        self.PDScaling = c_int16(0)
        self.BGSubtract = c_uint8(0)
        self.DispCorr = c_uint8(55)
        self.InterpData = c_uint8(55)
        self.FFT = c_uint8(55)
        self.procRoiBegin = c_uint16(0)
        self.procRoiEnd = c_uint16(2047)
        self.InterpFilter = c_uint8(55)
        self.InterpDownsample = c_uint8(55)
        self.klinRoiBegin = c_uint16(5)
        self.klinRoiEnd = c_uint16(1165)
        self.SampleOffset = c_int16(6)
        self.MZIScaleFactor = c_int16(0)
        self.MZIHPFilter = c_uint8(55)
        self.DownsampleFactor = c_int16(0)
        self.SamplesPerTrig = c_int16(600)
        self.Ch0OffsetTweak = c_int16(0)
        self.SyncTrig = c_uint8(55)
        self.FFT16b = c_uint8(0)
        self.Polar = c_uint8(0)
        self.MagOnly = c_uint8(0)
        self.FFTReImgScale = c_int8(-8)
        self.InterpDSAvg = c_uint8(55)
        self.OnDutyCycleTrigLast = c_uint16(65535)
        self.OffDutyCycleTrigLast = c_uint16(65535)
        self.ProcessMZI = c_uint8(55)
        self.UseDutyCycle = c_uint8(0)
        
        
    def encodeToString(self):
        s = ""
        s = s + "\nCh0Shift= %d\nnumKlinPts= %d\nPDScaling= %d\nBGSubtract= %s" % (self.Ch0Shift, self.numKlinPts, self.PDScaling, repr(self.BGSubtract > 0)) 
        s = s + "\nDisp Corr= %s\nInterp= %s\nFFT= %s\nprocROIBegin= %d\nprocROIEnd= %d" % (repr(self.DispCorr > 0), repr(self.InterpData > 0), repr(self.FFT > 0), self.procRoiBegin, self.procRoiEnd)
        s = s + "\nInterp Filter= %s\nInterp Downsample= %s\nklinRoiBegin= %d\nklinRoiEnd= %d\nsampleOffset= %d" % (repr(self.InterpFilter > 0), repr(self.InterpDownsample > 0), self.klinRoiBegin, self.klinRoiEnd, self.SampleOffset)
        s = s + "\nMZIScaleFactor= %d\nMZIHPFilter=%s \nDownsampleFactor= %d" % (self.MZIScaleFactor, repr(self.MZIHPFilter > 0), self.DownsampleFactor)
        s = s + "\nSamplesPerTrig= %d \nCh0OffsetTweak=%d \nSyncTrig= %s" % (self.SamplesPerTrig, self.Ch0OffsetTweak, repr(self.SyncTrig > 0))
        s = s + "\nFFT16b= %s \nPolar= %s \nMagOnly= %s \nProcessMZI= %s" % (repr(self.FFT16b > 0), repr(self.Polar > 0), repr(self.MagOnly > 0), repr(self.ProcessMZI > 0))
        s = s + "\nFFTReImgScale= %d \nInterpDSAvg= %s\nOnDutyCycleTrigLast= %d \nOffDutyCycleTrigLast= %d" % (self.FFTReImgScale, repr(self.InterpDSAvg > 0), self.OnDutyCycleTrigLast, self.OffDutyCycleTrigLast)
        s = s + "\nUseDutyCycle= %s" % repr(self.UseDutyCycle > 0)
        
        return s

    def convertFromBool(self, val):
        a = c_uint8(0)
        if val == ' True':
            a = c_uint8(55)
            
        return a
        
    def decodeFromString(self, opts_str):
        lines = re.split('\n', opts_str)  # break up lines into array
        for s in lines:
            x = re.split('=', s)
            if(len(x) < 2):
                continue
            fld = x[0]
            val = x[1]
            if(fld == 'Ch0Shift'):
                #val2 = re.split(' ', val)
                self.Ch0Shift = c_uint16(int(val))
            elif(fld == 'numKlinPts'):
                self.numKlinPts = c_uint16(int(val))
            elif(fld == 'PDScaling'):
                self.PDScaling = c_int16(int(val))                
            elif(fld == 'BGSubtract'):
                self.BGSubtract = self.convertFromBool(val)
            elif(fld == 'Disp Corr'):
                self.DispCorr = self.convertFromBool(val)
            elif(fld == "ProcessMZI"):
                self.ProcessMZI = self.convertFromBool(val)
            elif(fld == 'Interp'):
                self.InterpData = self.convertFromBool(val)
            elif(fld == 'FFT'):
                self.FFT = self.convertFromBool(val)
            elif(fld == 'procROIBegin'):
                self.procRoiBegin = c_uint16(int(val))
            elif(fld == 'procROIEnd'):
                self.procRoiEnd = c_uint16(int(val))                
            elif(fld == 'Interp Filter'):
                self.InterpFilter = self.convertFromBool(val)
            elif(fld == 'Interp Downsample'):
                self.InterpDownsample = self.convertFromBool(val)
            elif(fld == 'klinRoiBegin'):
                self.klinRoiBegin = c_uint16(int(val))
            elif(fld == 'klinRoiEnd'):
                self.klinRoiEnd = c_uint16(int(val))
            elif(fld == 'sampleOffset'):
                self.SampleOffset = c_int16(int(val))
            elif(fld == 'MZIScaleFactor'):
                self.MZIScaleFactor = c_int16(int(val))
            elif(fld == 'MZIHPFilter'):
                self.MZIHPFilter = self.convertFromBool(val)
            elif(fld == 'DownsampleFactor'):
                self.DownsampleFactor = c_int16(int(val))
            elif(fld == 'SamplesPerTrig'):
                self.SamplesPerTrig = c_int16(int(val))
            elif(fld == 'Ch0OffsetTweak'):
                self.Ch0OffsetTweak = c_int16(int(val))
            elif(fld == 'DownsampleFactor'):
                self.SyncTrig = self.convertFromBool(val)
            elif(fld == "FFT16b"):
                self.FFT16b = self.convertFromBool(val)
            elif(fld == "Polar"):
                self.Polar = self.convertFromBool(val)
            elif(fld == "MagOnly"):
                self.MagOnly = self.convertFromBool(val)
            elif(fld == "FFTReImgScale"):
                self.FFTReImgScale = c_int8(int(val))
            elif(fld == "InterpDSAvg"):
                self.InterpDSAvg = self.convertFromBool(val)
            elif(fld == "OnDutyCycleTrigLast"):
                self.OnDutyCycleTrigLast = c_uint16(int(val))
            elif(fld == "OffDutyCycleTrigLast"):
                self.OffDutyCycleTrigLast = c_uint16(int(val))
            elif(fld == "UseDutyCycle"):
                self.UseDutyCycle = self.convertFromBool(val)
                
class LV_DLLInterface:
    def __init__(self):
        self.isInitialized = False
        self.fpgaOpts = FPGAOpts_t()        
    
    """ InitInterface: loads the LabVIEW OCT FPGA interface DLL and maps C routines to Python equivalents
    """ 
    def InitInterface(self):
        # load in DLL
        # see header file LV_OCT_Raw_Data.h for C function definitions
        print("loading OCT processed data interface DLL")
        #lvlib = CDLL('Z:\\dev\\LV OCT Interface\\DLL\\LV_OCT_Interface')
        lvlib = CDLL('..\\dll\\LV_OCT_Interface')
        self.lvlib = lvlib
        #lvlib = CDLL('LV_OCT_Interface')
        
        # setup prototypes to functions in library
        # see the header file LV_OCT_Interface.h
        klin_data_t = np_ctypes.ndpointer(c_int32, flags="C_CONTIGUOUS")
        len_ptr_t = POINTER(c_int32)
        
        #fpga_ref_t = POINTER(c_ulong)
        #imaq_session_t = POINTER(c_ulong)
        init_fpga = lvlib.InitFPGA
        init_fpga.argtypes = [c_uint16, POINTER(FPGAOpts_t), klin_data_t, len_ptr_t]
        init_fpga.resttype = c_int32
        self.init_fpga = init_fpga
        
        config_fpga_acq = lvlib.ConfigureFPGAAcquisition
        # oct setup num, num triggers, FPGAopts, ROI_out, num trigs out, fpga processing opts out
        config_fpga_acq.argtypes = [c_uint16, c_uint32, POINTER(FPGAOpts_t), POINTER(c_int32), POINTER(c_int32), POINTER(FPGAOpts_t)]
        config_fpga_acq.restypes = c_int32
        self.config_fpga_acq = config_fpga_acq
        
        acq_fpga_data = lvlib.AcquireFPGAData
        
        data_arr_t = np_ctypes.ndpointer(c_float, flags="C_CONTIGUOUS")
        packed_data_arr_t = np_ctypes.ndpointer(c_uint64, flags="C_CONTIGUOUS")
        time_ptr_t = POINTER(c_uint32)
        #int32_t __cdecl AcquireFPGAData(OCTSetup OCTSetup, uint32_t NumTriggers, 
        #	uint32_t NumSamples, uint32_t TriggerOffset, LVBoolean UnpackData, 
        #	float data_re_out[], float data_im_out[], int32_t *data_re_im_len, 
        #	uint64_t PackedDataOut[], int32_t *packed_data_len, uint32_t *TimeElapsedMs, 
        #     uint32_t *TransferTimeMs, uint32_t *UnpackTimeMs);
        
        acq_fpga_data.argtypes = [c_uint16, c_uint32, c_uint32, c_uint32, c_uint8, data_arr_t, data_arr_t, len_ptr_t, packed_data_arr_t, len_ptr_t, time_ptr_t, time_ptr_t, time_ptr_t]
        acq_fpga_data.restypes = c_int32
        self.acq_fpga_data = acq_fpga_data
        
        daqio_xyscan = lvlib.DAQIOXYScan
        scan_command_t =  np_ctypes.ndpointer(c_double, flags="C_CONTIGUOUS")
        daq_task_t = POINTER(c_ulong)
        daqio_xyscan.argtypes = [c_char_p, c_char_p, c_char_p, c_double, scan_command_t, scan_command_t, c_int32, c_int32, daq_task_t]
        daqio_xyscan.resttypes = c_int32
        self.daqio_xyscan = daqio_xyscan
        
        daqio_waitdone = lvlib.DAQIOWaitUntilDone
        daqio_waitdone.argtypes = [daq_task_t, c_double]
        daqio_waitdone.resttypes = c_int32
        self.daqio_waitdone = daqio_waitdone
        
        load_dispersion = lvlib.LoadDispersionFPGA
        disp_arr_t = np_ctypes.ndpointer(c_double, flags="C_CONTIGUOUS")
        load_dispersion.argtypes = [c_uint16, disp_arr_t, disp_arr_t, c_int32]
        load_dispersion.resttypes = c_int32
        self.load_dispersion = load_dispersion
        
        reset_fpga = lvlib.ResetFPGA
        reset_fpga.argtypes = [c_uint16, POINTER(FPGAOpts_t), klin_data_t, len_ptr_t]
        reset_fpga.resttypes = c_int32
        self.reset_fpga = reset_fpga
        
        int32_ptr = POINTER(c_int32)
        
        daqio_audio_twospk_mic_input = lvlib.DAQIOAudioTwoSpeakersMicInput
        # # args outputRate, inputRate, mic channel, trig src in, trig src out, right spk chan, left spk chan, right output, left output, num mic samples out, output task out, mic task out, right output length, left output length
        daqio_audio_twospk_mic_input.argtypes = [c_double, c_double, c_char_p, c_char_p, c_char_p, c_char_p, c_char_p, scan_command_t, scan_command_t, int32_ptr, daq_task_t, daq_task_t, c_int32, c_int32]
        daqio_audio_twospk_mic_input.resttypes = c_int32
        self.daqio_audio_twospk_mic_input = daqio_audio_twospk_mic_input
        
        daqio_audio_start_output_grab_mic_data = lvlib.DAQIOAudioStartOutputAndGrabMicData
        daqio_audio_start_output_grab_mic_data.argtypes = [ daq_task_t, c_char_p, daq_task_t, c_int32, scan_command_t] 
        daqio_audio_start_output_grab_mic_data.resttypes = c_int32
        self.daqio_audio_start_output_grab_mic_data = daqio_audio_start_output_grab_mic_data
        #int32_t __cdecl DAQIOAudioGrabMicData(uintptr_t *micTaskIn, double timeout, 
        #                  int32_t numMicSamples, double MicDataOut[], int32_t len);
        
        #int32_t __cdecl DAQIOAudioGrabMicData(uintptr_t *micTaskIn, double timeout, 
        #                  int32_t numMicSamples, double MicDataOut[], int32_t *len);
         
        daqio_audio_grab_mic_data = lvlib.DAQIOAudioGrabMicData
        daqio_audio_grab_mic_data.argtypes = [ daq_task_t, c_double, c_int32, scan_command_t, int32_ptr] 
        daqio_audio_grab_mic_data.resttypes = c_int32
        self.daqio_audio_grab_mic_data = daqio_audio_grab_mic_data
        
        daqio_set_mirrorpos = lvlib.DAQIOSetMirrorPostion
        daqio_set_mirrorpos.argtypes = [ c_double, c_double, c_char_p, c_char_p ]
        daqio_set_mirrorpos.resttypes = c_int32
        self.daqio_set_mirrorpos = daqio_set_mirrorpos
        
        daqio_set_attenlevel = lvlib.DAQIOSetAttenuatorLevel
        daqio_set_attenlevel.argtypes = [ c_uint8, c_char_p ]
        daqio_set_attenlevel.resttypes = c_int32
        self.daqio_set_attenlevel = daqio_set_attenlevel
        
        daqio_cam_trig = lvlib.DAQIOCamTriggerOnly
        # args cam trig command, trig command length, trig out line, cam trigr output channel, trig input line
        daqio_cam_trig.argtypes = [ scan_command_t, c_int32, c_double, c_char_p, c_char_p, c_char_p  ]
        daqio_cam_trig.resttypes = c_int32
        self.daqio_cam_trig = daqio_cam_trig        
        
        daqio_xyscan_cam_trig = lvlib.DAQIOXYScanWithCamTrig
        daqio_xyscan_cam_trig.argtypes = [c_char_p, c_char_p, c_char_p, scan_command_t, scan_command_t, scan_command_t, c_char_p, c_double, c_int32, daq_task_t]
        daqio_xyscan_cam_trig.resttypes = c_int32
        self.daqio_xyscan_cam_trig = daqio_xyscan_cam_trig        
        
        daqio_send_trig = lvlib.DAQIOSendTrigger
        daqio_send_trig.resttypes = c_int32
        daqio_send_trig.argtypes = [c_char_p]
        self.daqio_send_trig = daqio_send_trig        
        
        daq_reset = lvlib.DAQReset
        daq_reset.resttypes = c_int32
        daq_reset.argtypes = [c_char_p]
        self.daq_reset = daq_reset        
        
        recalc_klin = lvlib.RecalcKlin
        recalc_klin.argtypes = [c_uint16, POINTER(FPGAOpts_t), klin_data_t, len_ptr_t]
        recalc_klin.resttypes = c_int32
        self.recalc_klin = recalc_klin        
        
        close_fpga = lvlib.CloseFPGA
        close_fpga.resttypes = c_int32
        close_fpga.argtypes = [c_uint16]
        self.close_fpga = close_fpga
        
        self.isInitialized = True

    # setup = integer that indicates setup number (0 for 200 kHz, 1 for 50 kHz
    # sampleOffset
    # samplesPerTrig
    # Ch0Shift
    # klinNumPts
    # klinROI - array-like with two elements
    def InitFPGA(self, setup, sampleOffset, samplesPerTrig, Ch0Shift, klinNumPts, klinROI):
        if not self.isInitialized and setup >= 0:
            self.InitInterface()
            
        fpgaOpts = self.fpgaOpts
        fpgaOpts.SampleOffset = c_int16(sampleOffset)
        fpgaOpts.SamplesPerTrig = c_int16(samplesPerTrig)
        fpgaOpts.Ch0Shift = c_uint16(Ch0Shift)
        fpgaOpts.numKlinPts = c_uint16(klinNumPts)
        fpgaOpts.klinRoiBegin = c_uint16(klinROI[0])
        fpgaOpts.klinRoiEnd = c_uint16(klinROI[1])
        
        self.fpgaOpts = fpgaOpts
        
        numpts = c_int32(4096)
        klin = np.zeros((4096), np.int32)
        klin = np.require(klin, np.int32, ['C', 'W'])
        numpts = c_int32(4096)
        err = 0
        if setup >= 0:
            err = self.init_fpga(c_uint16(setup), byref(fpgaOpts), klin, byref(numpts))
            
        klin = klin[0:klinNumPts-1]
        self.setupNum = setup
        return err, klin
        
    # IsOCTTestingMode: returns True if OCT hardware is NOT present, otherwise returns False
    # this does not do any detection, but returns based off setup number, to determine if in testing mode
    # useful when hardware is not present for testing othe program features
    def IsOCTTestingMode(self):
        return self.setupNum<0

    # IsDAQTestingMode: returns True if DAQ hardware is NOT present, otherwise returns False
    # this does not do any detection, but returns based off setup number, to determine if in testing mode
    # useful when hardware is not present for testing othe program features
    def IsDAQTestingMode(self):
        return self.setupNum < 0
            
    # grab the output of the FFT 
    def AcquireOCTDataFFT(self, numTriggers, zROI, startTrigOffset=0, dispCorr=True, downsample=0):
        setupNum = c_uint16(self.setupNum)
        
        fpgaOpts = self.fpgaOpts
        dsFactor = downsample + 1
        fpgaOpts.procRoiBegin = c_uint16(zROI[0])
        fpgaOpts.procRoiEnd = c_uint16(zROI[1])
        fpgaOpts.StartTriggerOfffset = c_uint16(startTrigOffset*dsFactor)
        fpgaOpts.ProcessMZI = c_uint8(55) # True
        fpgaOpts.InterpData = c_uint8(55) # True
        if dispCorr:
            fpgaOpts.DispCorr = c_uint8(55) # True
        else:
            fpgaOpts.DispCorr = c_uint8(0) # False
        fpgaOpts.FFT16b = c_uint8(0)
        fpgaOpts.Polar = c_uint8(0)
        fpgaOpts.MagOnly = c_uint8(0)
        fpgaOpts.DownsampleFactor = c_int16(downsample) 
        fpgaOpts.FFT = c_uint8(55)  # True
        samplesPerTrig = fpgaOpts.SamplesPerTrig
        
        fpgaOptsOut =  FPGAOpts_t()
        numTrigsOut = c_int32(numTriggers)
        roiSize = zROI[1]-zROI[0]+1
        roiSizeOut = c_int32(roiSize)
        err = self.config_fpga_acq(setupNum, c_uint32(numTriggers*dsFactor), byref(fpgaOpts), byref(roiSizeOut), byref(numTrigsOut), byref(fpgaOptsOut))
        DebugLog.log("AcquireOCTDataFFT  roiSizeout= %d numTrigsOut=  %d" % (roiSizeOut.value,  numTrigsOut.value))
        if err < 0:
            return err, None
        
        numSamples = roiSizeOut.value * numTrigsOut.value
        numDataPts = numSamples
        numPackedDataPts = 0
                
        len_data = c_int32(numDataPts)
        len_packed_data = c_int32(numPackedDataPts)
        
        d_re = np.zeros((numDataPts), np.float32)
        d_im = np.zeros((numDataPts), np.float32)        
        d_re = np.require(d_re, np.float32, ['C', 'W'])
        d_im = np.require(d_im, np.float32, ['C', 'W'])
        
        packedData = np.zeros((numPackedDataPts), np.uint64)
        packedData = np.require(packedData, np.uint64, ['C', 'W'])
        
        #trigOffset = fpgaOpts.StartTriggerOfffset
        timeElapsed = c_uint32(0)
        transferTime = c_uint32(0)
        unpackTime = c_uint32(0)
        trigOffset = c_uint32(0)
        unpackData = True
        
        # DebugLog.log("OCTDataCollector.startFrameGetData(): isSynchOCT = " + repr(self.protocol.isSynchOCT()))
        err = self.acq_fpga_data(setupNum, c_uint32(numTrigsOut.value), c_uint32(numSamples), trigOffset, unpackData, d_re, d_im, byref(len_data), packedData, byref(len_packed_data), byref(timeElapsed), byref(transferTime), byref(unpackTime))
        DebugLog.log("AcquireOCTDataFFT  len_data= " + repr(len_data.value))
        oct_data = np.zeros(len_data.value, np.complex)
        oct_data[:] = d_re + 1j * d_im;
        oct_data = oct_data.reshape((numTrigsOut.value, roiSizeOut.value))
    
        return err, oct_data
    
    # grab the interpolated photodiode data
    def AcquireOCTDataInterpPD(self, numTrigs, downsample=0):
        setupNum = c_uint16(self.setupNum)
        
        dsFactor = downsample + 1
        fpgaOpts = self.fpgaOpts
        fpgaOpts.DispCorr = c_uint8(0) # False
        fpgaOpts.FFT = c_uint8(0)   # False
        fpgaOpts.ProcessMZI = c_uint8(55)  # True
        fpgaOpts.InterpData = c_uint8(55)  # True
    #    fpgaOpts.StartTriggerOfffset = c_uint16(startTrigOffset)
        fpgaOpts.FFT16b = c_uint8(0)
        fpgaOpts.Polar = c_uint8(0)
        fpgaOpts.MagOnly = c_uint8(0)
        fpgaOpts.DownsampleFactor = c_int16(downsample)
        
        fpgaOptsOut =  FPGAOpts_t()
        numTrigsOut = c_int32(numTrigs)
        
        roiSize = 2048
        roiSizeOut = c_int32(roiSize)
        err = self.config_fpga_acq(setupNum, c_uint32(numTrigs*dsFactor), byref(fpgaOpts), byref(roiSizeOut), byref(numTrigsOut), byref(fpgaOptsOut))
        numSamples = roiSizeOut.value * numTrigsOut.value
        if err < 0:
            return err, None
            
        numDataPts = numSamples
        numPackedDataPts = 0
                
        len_data = c_int32(numDataPts)
        len_packed_data = c_int32(numPackedDataPts)
        
        d_re = np.zeros((numDataPts), np.float32)
        d_im = np.zeros((numDataPts), np.float32)        
        d_re = np.require(d_re, np.float32, ['C', 'W'])
        d_im = np.require(d_im, np.float32, ['C', 'W'])
        
        packedData = np.zeros((numPackedDataPts), np.uint64)
        packedData = np.require(packedData, np.uint64, ['C', 'W'])
        
        #trigOffset = fpgaOpts.StartTriggerOfffset
        timeElapsed = c_uint32(0)
        transferTime = c_uint32(0)
        unpackTime = c_uint32(0)
        trigOffset = c_uint32(0)
        unpackData = c_uint8(55)
        
        # DebugLog.log("OCTDataCollector.startFrameGetData(): isSynchOCT = " + repr(self.protocol.isSynchOCT()))
        err = self.acq_fpga_data(setupNum, c_uint32(numTrigsOut.value), c_uint32(numSamples), trigOffset, unpackData, d_re, d_im, byref(len_data), packedData, byref(len_packed_data), byref(timeElapsed), byref(transferTime), byref(unpackTime))
        
        interp_pd = np.zeros(len_data.value, np.double)
        interp_pd[:] = d_re 
        interp_pd = interp_pd.reshape((numTrigsOut.value, roiSizeOut.value))
    
        return err, interp_pd
    
    # return only the magnitude in a packed 64-bit integer
    # this gives max speed for imaging only acquisiiton in multiprocess mode
    def AcquireOCTDataMagOnly(self, numTriggers, zROI, startTrigOffset=0, dispCorr=True, downsample=0):
        setupNum = c_uint16(self.setupNum)
        
        dsFactor = downsample + 1
        fpgaOpts = copy.copy(self.fpgaOpts)
        fpgaOpts.procRoiBegin = c_uint16(zROI[0])
        fpgaOpts.procRoiEnd = c_uint16(zROI[1])
        fpgaOpts.StartTriggerOfffset = c_uint16(startTrigOffset*dsFactor)
        fpgaOpts.ProcessMZI = c_uint8(55) # True
        fpgaOpts.InterpData = c_uint8(55) # True
        if dispCorr:
            fpgaOpts.DispCorr = c_uint8(55) # True
        else:
            fpgaOpts.DispCorr = c_uint8(0) # False
        fpgaOpts.FFT = c_uint8(55)  # True
        fpgaOpts.FFT16b = c_uint8(55)
        fpgaOpts.Polar = c_uint8(55)
        fpgaOpts.MagOnly = c_uint8(55)
        fpgaOpts.DownsampleFactor = c_int16(downsample)
        
        samplesPerTrig = fpgaOpts.SamplesPerTrig
        
        fpgaOptsOut =  FPGAOpts_t()
        numTrigsOut = c_int32(numTriggers)
        roiSize = zROI[1]-zROI[0]+1
        roiSizeOut = c_int32(roiSize)
        err = self.config_fpga_acq(setupNum, c_uint32(numTriggers*dsFactor), byref(fpgaOpts), byref(roiSizeOut), byref(numTrigsOut), byref(fpgaOptsOut))
        DebugLog.log("AcquireOCTDataMagOnly: roiSizeout= %d numTrigsOut=  %d" % (roiSizeOut.value,  numTrigsOut.value))
        if err < 0:
            return err, None
        
        numSamples = roiSizeOut.value * numTrigsOut.value
        numDataPts = 0
        numPackedDataPts = numSamples 
                
        len_data = c_int32(numDataPts)
        len_packed_data = c_int32(numPackedDataPts)
        
        t1 = time.time()
        
        d_re = np.zeros((numDataPts), np.float32)
        d_im = np.zeros((numDataPts), np.float32)        
        d_re = np.require(d_re, np.float32, ['C', 'W'])
        d_im = np.require(d_im, np.float32, ['C', 'W'])
        
        packedData = np.zeros((numPackedDataPts), np.uint64)
        packedData = np.require(packedData, np.uint64, ['C', 'W'])

        DebugLog.log("AcquireOCTDataMagOnly: mem alloc time= %0.1f ms " % (1000*(time.time() - t1)))

        
        #trigOffset = fpgaOpts.StartTriggerOfffset
        timeElapsed = c_uint32(0)
        transferTime = c_uint32(0)
        unpackTime = c_uint32(0)
        trigOffset = c_uint32(0)
        unpackData = c_uint8(0)
        
        # DebugLog.log("OCTDataCollector.startFrameGetData(): isSynchOCT = " + repr(self.protocol.isSynchOCT()))
        t1 = time.time()
        err = self.acq_fpga_data(setupNum, c_uint32(numTrigsOut.value), c_uint32(numSamples), trigOffset, unpackData, d_re, d_im, byref(len_data), packedData, byref(len_packed_data), byref(timeElapsed), byref(transferTime), byref(unpackTime))
        DebugLog.log("AcquireOCTDataMagOnly: grab time= %0.1f ms len_packed_data= %d" % (1000*(time.time() - t1),  len_packed_data.value))
        
        packedData = packedData.reshape((numTrigsOut.value, roiSizeOut.value))
        return err, packedData
    
    # grabb the output of the FFT 
    def AcquireOCTDataRaw(self, numTriggers, samplesPerTrig=-1, Ch0Shift=-1, startTrigOffset=0, downsample=0):
        setupNum = c_uint16(self.setupNum)
        
        dsFactor = downsample + 1
        
        fpgaOpts = copy.copy(self.fpgaOpts)
        fpgaOpts.StartTriggerOfffset = c_uint16(startTrigOffset*dsFactor)
        fpgaOpts.ProcessMZI = c_uint8(0) # False
        fpgaOpts.InterpData = c_uint8(0) # False
        fpgaOpts.DispCorr = c_uint8(0) # False
        fpgaOpts.FFT = c_uint8(0)  # False
        fpgaOpts.klinRoiBegin = c_uint16(0) 
        fpgaOpts.SampleOffset = c_int16(1)   # sample offset must be at least 1 for 200 khz system
        fpgaOpts.Ch0Shift = c_uint16(0)
        fpgaOpts.DownsampleFactor = c_int16(downsample)
        
        if Ch0Shift < 0:
            Ch0Shift = fpgaOpts.Ch0Shift
        else:
            Ch0Shift = Ch0Shift // 2
        
        fpgaOpts.Ch0Shift = c_uint16(0) 
        
        if samplesPerTrig < 0:
            samplesPerTrig = fpgaOpts.SamplesPerTrig
        else:
            samplesPerTrig = samplesPerTrig // 2
        
        samplesPerTrig += Ch0Shift 
        fpgaOpts.SamplesPerTrig = samplesPerTrig
        DebugLog.log("AcquireOCTDataRaw  samplesPerTrig= %d Ch0Shift=  %d" % (samplesPerTrig, Ch0Shift))
        fpgaOpts.klinRoiEnd = c_uint16(samplesPerTrig * 2) 
        
        fpgaOptsOut =  FPGAOpts_t()
        numTrigsOut = c_int32(numTriggers)
        roiSize = samplesPerTrig * 2
        roiSizeOut = c_int32(roiSize)
        err = self.config_fpga_acq(setupNum, c_uint32(numTriggers*dsFactor), byref(fpgaOpts), byref(roiSizeOut), byref(numTrigsOut), byref(fpgaOptsOut))
        DebugLog.log("AcquireOCTDataRaw  config_fpga_acq err= %d" % err)
        DebugLog.log("AcquireOCTDataRaw  roiSizeout= %d numTrigsOut=  %d" % (roiSizeOut.value,  numTrigsOut.value))
        if err < 0:
            return err, None, None
        
        numSamples = roiSizeOut.value * numTrigsOut.value
        numDataPts = numSamples
        numPackedDataPts = 0
                
        len_data = c_int32(numDataPts)
        len_packed_data = c_int32(numPackedDataPts)
        
        d_re = np.zeros((numDataPts), np.float32)
        d_im = np.zeros((numDataPts), np.float32)        
        d_re = np.require(d_re, np.float32, ['C', 'W'])
        d_im = np.require(d_im, np.float32, ['C', 'W'])
        
        packedData = np.zeros((numPackedDataPts), np.uint64)
        packedData = np.require(packedData, np.uint64, ['C', 'W'])
        
        #trigOffset = fpgaOpts.StartTriggerOfffset
        timeElapsed = c_uint32(0)
        transferTime = c_uint32(0)
        unpackTime = c_uint32(0)
        trigOffset = c_uint32(0)
        unpackData = c_uint8(55)
        
        # DebugLog.log("OCTDataCollector.startFrameGetData(): isSynchOCT = " + repr(self.protocol.isSynchOCT()))
        DebugLog.log("AcquireOCTDataRaw numSamples= " + repr(numSamples))
        err = self.acq_fpga_data(setupNum, c_uint32(numTrigsOut.value), c_uint32(numSamples), trigOffset, unpackData, d_re, d_im, byref(len_data), packedData, byref(len_packed_data), byref(timeElapsed), byref(transferTime), byref(unpackTime))
        DebugLog.log("AcquireOCTDataRaw len_data= " + repr(len_data.value) + " transferTime= " + repr(transferTime.value))
        pd_data = np.zeros(len_data.value, np.float32)
        mzi_data = np.zeros(len_data.value, np.float32)
        pd_data[:] = d_re
        pd_data = pd_data.reshape((numTrigsOut.value, roiSizeOut.value))
        pd_data = pd_data[:, 0:samplesPerTrig*2]
        DebugLog.log("AcquireOCTDataRaw pd_data.shape= %s" % repr(pd_data.shape)) 
        mzi_data[:] = d_im
        mzi_data = mzi_data.reshape((numTrigsOut.value, roiSizeOut.value))
        mzi_data = mzi_data[:, 0:samplesPerTrig*2]
        
        return err, pd_data, mzi_data
    
    def LoadOCTDispersion(self, magWin, phaseCorr):
        err = 0
        if self.setupNum >= 0:
            len_data = c_int32(len(magWin))
            magWin_d = np.require(magWin, np.double, ['C'])
            phaseCorr_d = np.require(phaseCorr, np.double, ['C'])
            setupNum = c_uint16(self.setupNum)
            err = self.load_dispersion(setupNum, magWin_d, phaseCorr_d, len_data)
            
        return err
    
    def SetAttenLevel(self, attenLvl, outLines):
        err = self.daqio_set_attenlevel(c_uint8(attenLvl), outLines.encode("ASCII"))
        return err
        #err1 = self._octhw.daqio_set_attenlevel(c_uint8(attenLevels[0]), self.audioHW.attenL_daqChan.encode("ASCII"))
    
    def CloseFPGA(self):
        err = 0
        if self.isInitialized:
            err = self.close_fpga(self.setupNum)
            self.isInitialized = False
            
        return err
    
    def Shutdown(self):
        pass
        
                    
def StartOCTInterfaceBGProcess(basePath, OCTtrigRate, rawDataQSize=3):
    # check if process exists 
    filePath = os.path.join(basePath, 'collector_PID_tmp.txt')
    if(os.path.isfile(filePath)):
        f = open(filePath, "r")
        pidStr = f.readline()
        f.close()
        collPID = int(pidStr)
        DebugLog.log("StartOCTInterfaceBGProcess(): killing collection PID= " + repr(collPID))
        try:
            p = psutil.Process(collPID)
            p.kill()
        except Exception as ex:
            DebugLog.log("StartOCTInterfaceBGProcess():  Exception attempting to kill collection process")
            traceback.print_exc(file=sys.stdout)

    rawDataQ = mproc.Queue(rawDataQSize)
    msgQ = mproc.Queue(10)
    statusQ = mproc.Queue(10)
    collectorProcess = mproc.Process(target=_OCTBGLoop, args=[rawDataQ, msgQ, statusQ, OCTtrigRate], daemon=True)
    DebugLog.log("StartOCTInterfaceBGProcess(): starting collector process")
    collectorProcess.start()
    collPID = collectorProcess.pid
    f = open(filePath, "w")
    f.write(repr(collPID))
    f.close()
        
    adaptor = LV_DLLInterface_BGProcess_Adaptor(rawDataQ, msgQ, statusQ)
    return adaptor
    
def _OCTBGLoop(rawDataQ, collMsgQ, statusQ, OCTtrigRate):
    bgProcess = LV_DLLInterface_BGProcess(rawDataQ, collMsgQ, statusQ, OCTtrigRate)
    bgProcess.loop()
    

# interface to the as a background process
class LV_DLLInterface_BGProcess:
    def __init__(self, rawDataQ, collMsgQ, statusQ, OCTtrigRate):
        self.oct_hw = LV_DLLInterface()
        self.shutdown = False
        self.rawDataQ = rawDataQ
        self.msgQ = collMsgQ
        self.statusQ = statusQ 
        self.OCTtrigRate = OCTtrigRate
        
    def handleMessage(self, msg):
        msgType = msg[0]
        param = msg[1]
        DebugLog.log('LV_DLLInterface_BGProcess: got message ' + repr(msgType) + " param= " + repr(param))
        if msgType == 'shutdown':
            self.shutdown = True
        elif msgType == 'initFPGA':
            err, klin = self.oct_hw.InitFPGA(param[0], param[1], param[2], param[3], param[4], param[5])
            statusMsg = OCTCommon.StatusMsg(OCTCommon.StatusMsgSource.COLLECTION, OCTCommon.StatusMsgType.KLIN)
            statusMsg.param = [err, klin]
            self.statusQ.put(statusMsg)
        elif msgType == 'closeFPGA':
            self.acquireData = False
            self.oct_hw.CloseFPGA()
        elif msgType == 'acquire':
            self.acquireData = True
        elif msgType == 'pause':
            self.acquireData = False
        elif msgType == 'newAcquisition':
            self.frameNum = 0
        elif msgType == 'setAcqFunction':
            self.acqFunction = param
        elif msgType == 'setAcqFunctionArgs':
            self.acqFunctionArgs = param
        elif msgType == 'setAttenLvl':
            err = self.oct_hw.SetAttenLevel(param[0], param[1])
            #statusMsg = OCTCommon.StatusMsg(OCTCommon.StatusMsgSource.COLLECTION, OCTCommon.StatusMsgType.BLANK)
            #statusMsg.param = err
            self.statusQ.put(err)
        elif msgType == 'loadDispersion':
            magWin = param[0]
            phaseCorr = param[1]
            err = self.oct_hw.LoadOCTDispersion(magWin, phaseCorr)
            #statusMsg = OCTCommon.StatusMsg(OCTCommon.StatusMsgSource.COLLECTION, OCTCommon.StatusMsgType.BLANK)
            #statusMsg.param = err
            self.statusQ.put(err)
        elif msgType == 'getFPGAOpts':
            self.statusQ.put(self.oct_hw.fpgaOpts)
        elif msgType == 'setFPGAOpts':
            fpgaOpts = param
            self.oct_hw.fpgaOpts = fpgaOpts
        elif msgType == 'setSendExtraInfo':
            self.sendExtraInfo = param
        elif msgType == 'acquireFFT':
            err, oct_data = self.oct_hw.AcquireOCTDataFFT(param[0], param[1], param[2], param[3])
            self.rawDataQ.put((err, oct_data))
        elif msgType == 'acquireInterpPD':
            err, interp_pd = self.oct_hw.AcquireOCTDataInterpPD(param)
            self.rawDataQ.put((err, interp_pd))
        elif msgType == 'acquireRaw':
            err, pd_data, mzi_data = self.oct_hw.AcquireOCTDataRaw(param[0], param[1], param[2], param[3])
            self.rawDataQ.put((err, pd_data, mzi_data))
        else:
            DebugLog.log('LV_DLLInterface_BGProcess: unknown message ' + repr(msgType))
            
            
    # this function once called will loop until a shutdown nmessage is received
    # it will call the acquireFunction as set by the user, when acquireData is true - which is changed b
    # by sending a 'acquire' message through the message quee
    def loop(self):
        self.acquireData = False
        self.acqFunction = None
        self.acqFunctionArgs = None
        self.sendExtraInfo = False
        self.frameNum = 0
        putTimeout = False
        
        # keep looing until we get a shutdown message
        while not self.shutdown:
            if self.acquireData and self.acqFunction is not None:
                try:
                    # call the acquire function to get OCT data - this function is set by the user
                    if not putTimeout:
                        DebugLog.log("LV_DLLInterface_BGProcess.loop(): acquiring frame " + repr(self.frameNum))
                        data, extraOutput = self.acqFunction(self.oct_hw, self.frameNum, self.OCTtrigRate, self.acqFunctionArgs) 
                        numTimeouts = 0
                    if data is not None:
                        try:
                            # try send raw data to consumer - a queue.Full exception may be thrown 
                            # in which case we will try to send it again next loop iteration
                            self.rawDataQ.put(data, True, 0.25)   
                            putTimeout = False
                            self.frameNum += 1   # increment frame number
                        except queue.Full as ex:
                            putTimeout = True
                            numTimeouts += 1 
                            if numTimeouts == 10:   # if too many timeouts, abandon collection
                                self.acquireData = False
                                putTimeout = False
                        
                        if extraOutput is not None and self.sendExtraInfo:
                            statusMsg = OCTCommon.StatusMsg(OCTCommon.StatusMsgSource.COLLECTION, OCTCommon.StatusMsgType.DAQ_OUTPUT)
                            statusMsg.param = extraOutput
                            try:
                                self.statusQ.put(statusMsg, False)   # send extraOutput to consumer
                            except queue.Full as ex:
                                pass
                    else:
                        self.acquireData = False
                    time.sleep(0.005)  # sleep for 5 ms to reduce CPU load
                except Exception as ex:
                    traceback.print_exc(file=sys.stdout)
                    statusMsg = OCTCommon.StatusMsg(OCTCommon.StatusMsgSource.COLLECTION, OCTCommon.StatusMsgType.ERROR)
                    statusMsg.param = ex
                    try:
                        self.statusQ.put(statusMsg, False)
                    except queue.Full as ex:
                        pass
                    self.acquireData = False
            else:
                time.sleep(0.02)  # sleep for 20 ms, to reduce CPU load when not acquiring
                
            while not self.msgQ.empty():
                try:
                    msg = self.msgQ.get()
                    self.handleMessage(msg)
                except Exception as ex:
                    statusMsg = OCTCommon.StatusMsg(OCTCommon.StatusMsgSource.COLLECTION, OCTCommon.StatusMsgType.ERROR)
                    statusMsg.param = ex
                    traceback.print_exc(file=sys.stdout)
                    try:
                        self.statusQ.put(statusMsg, timeout=0.25)
                    except Exception as ex:
                        traceback.print_exc(file=sys.stdout)
                    
            time.sleep(0.005)  # sleep for 5 ms, to avoid saturing CPU
            sys.stdout.flush()
            
# LV_DLLInterface_BGProcess_Adaptor: wraps most methods in LV_DLLInterface, but uses BG process interface through queues
# this can be sued on the GUI/controller/consumer side so that we can interact with the LV_DLLInterface as if it were running in same memory space
# Only functions are not wrapped are the AcquireXXX functions            
class LV_DLLInterface_BGProcess_Adaptor:
    def __init__(self, rawDataQ, collMsgQ, statusQ):
        self.rawDataQ = rawDataQ
        self.collMsgQ = collMsgQ
        self.statusQ = statusQ
        self.setupNum = None
        self.fpgaOpts = None
        self.qTimeout = 2
        
    def InitFPGA(self, setup, sampleOffset, samplesPerTrig, Ch0Shift, klinNumPts, klinROI):
        msg = ['initFPGA', (setup, sampleOffset, samplesPerTrig, Ch0Shift, klinNumPts, klinROI)]
        self.setupNum = setup
        self.collMsgQ.put(msg)
        respValid = False
        t1 = time.time()
        err = -1
        klin = []
        while not respValid and (time.time() - t1) < 5:
            if not self.statusQ.empty():
                resp = self.statusQ.get()
                if isinstance(resp, OCTCommon.StatusMsg):
                    if resp.msgType == OCTCommon.StatusMsgType.KLIN:
                        param = resp.param
                        err = param[0]
                        klin = param[1]
                        respValid = True
                
        return err, klin
        
    @property        
    def fpgaOpts(self):
        msg = ['getFPGAOpts', 0]
        self.collMsgQ.put(msg)
        resp = None
        while not isinstance(resp, FPGAOpts_t):
            resp = self.statusQ.get(self.qTimeout)
            
        fpgaOpts = resp
        return fpgaOpts
        
    @fpgaOpts.setter
    def fpgaOpts(self, fpgaOpts):
        if fpgaOpts is not None and isinstance(fpgaOpts, FPGAOpts_t):
            msg = ['setFPGAOpts', fpgaOpts]
            self.collMsgQ.put(msg)
        
    def IsOCTTestingMode(self):
        return self.setupNum < 0

    def IsDAQTestingMode(self):
        return self.setupNum == -1
        
    def LoadOCTDispersion(self, magWin, phaseCorr):
        msg = ['loadDispersion', (magWin, phaseCorr)]
        self.collMsgQ.put(msg, timeout=self.qTimeout)
        err = -1
        if not self.statusQ.empty():
            err = self.statusQ.get(timeout=self.qTimeout)
        
        if err is None:
            err = -1
            
        return err
        
    def CloseFPGA(self):
        msg = ['closeFPGA', 0]
        self.collMsgQ.put(msg, timeout=self.qTimeout)
        resp = 0
        if not self.statusQ.empty():
            resp = self.statusQ.get(timeout=self.qTimeout)
            if resp is None:
                resp = -1
        return resp

            
    def StartAcquisition(self):
        msg = ['acquire', 0]            
        self.collMsgQ.put(msg, timeout=self.qTimeout)
         
    def PauseAcquisition(self):
        msg = ['pause', 0]
        self.collMsgQ.put(msg, timeout=self.qTimeout)
        
    def NewAcquisition(self):
        msg = ['newAcquisition', 0]
        self.collMsgQ.put(msg, timeout=self.qTimeout)
        self.ClearDataQ()
        
    # sets the acquisition function, that will be called when collecting data
    # the acquisition function should have the form acqFunction(oct_hw, frameNum, extraArgs)
    # where oct_hw is an LV_DLL_Interface, frameNum is the frameNum, and extraArgs is an array-like (list or tuple), collection of arguments
    def SetAcqFunction(self, acqFunction):
        msg = ['setAcqFunction', acqFunction]
        self.collMsgQ.put(msg, timeout=self.qTimeout)
        
    # sets the "extra" arguments that will be passed to the acquisiton function when it is called
    # this is useful for passing inforation such as scan parameters, mirror Driver, etc.
    def SetAcqFunctionArgs(self, acqFunctionArgs):
        msg = ['setAcqFunctionArgs', acqFunctionArgs]
        self.collMsgQ.put(msg, timeout=self.qTimeout)
        
    # whether to enable sending of extra info generated by acqwuisition function (second output param)
    # setting to False could improve acquisition speed
    def SetSendExtraInfo(self, sendInfo):
        msg = ['setSendExtraInfo', sendInfo]
        self.collMsgQ.put(msg, timeout=self.qTimeout)
        
    def GetData(self):
        data = None
        if not self.rawDataQ.empty():
            data = self.rawDataQ.get(timeout=self.qTimeout)
        return data
        
    def GetStatus(self):
        status = None
        if not self.statusQ.empty():
            status = self.statusQ.get(timeout=self.qTimeout)
        return status
        
    def ClearDataQ(self):
        while not self.rawDataQ.empty():
            data = self.rawDataQ.get(timeout=self.qTimeout)
        
    def Shutdown(self):
        msg = ('shutdown', 0)
        self.collMsgQ.put(msg, timeout=self.qTimeout)
        
    def AcquireOCTDataFFT(self, numTriggers, zROI, startTrigOffset=0, dispCorr=True, downsample=0):
        self.ClearDataQ()
        msg = ('acquireFFT', (numTriggers, zROI, startTrigOffset, dispCorr, downsample))
        self.collMsgQ.put(msg, timeout=self.qTimeout)
        t = min((numTrigs*1e-4 + 500e-3, 5))
        (err, oct_data) = self.rawDataQ.get(timeout=t)
    
        return err, oct_data
    
    # grab the interpolated photodiode data
    def AcquireOCTDataInterpPD(self, numTrigs, downsample=0):
        self.ClearDataQ()
        msg = ('acquireInterpPD', (numTrigs, downsample))
        self.collMsgQ.put(msg, timeout=self.qTimeout)
        t = min((numTrigs*1e-4 + 500e-3, 5))
        (err, interp_pd) = self.rawDataQ.get(timeout=t)
        
        return err, interp_pd
    
    # grabb the output of the FFT 
    def AcquireOCTDataMagOnly(self, numTriggers, zROI):
        pass
    
    # grabb the output of the FFT 
    def AcquireOCTDataRaw(self, numTriggers, samplesPerTrig=-1, Ch0Shift=-1, startTrigOffset=0, downsample=0):
        self.ClearDataQ()
        msg = ('acquireRaw', (numTriggers, samplesPerTrig, Ch0Shift, startTrigOffset, downsample))
        t = min((numTriggers*1e-4 + 500e-3, 5))
        self.collMsgQ.put(msg, timeout=t)
        (err, pd_data, mzi_data) = self.rawDataQ.get(timeout=self.qTimeout)
        
        return err, pd_data, mzi_data        
        
    def SetAttenLevel(self, attenLvl, outLines):
        msg = ('setAttenLvl', (attenLvl, outLines))
        self.collMsgQ.put(msg, timeout=self.qTimeout)
        err = -1
        time.sleep(1e-3)
        if not self.statusQ.empty():
            err = self.statusQ.get(timeout=self.qTimeout)
            
        return err
        
def readFPGAOptsConfig(filepath):    
    opts = FPGAOpts_t()
    
    f = open(filepath, "r")
    txt = f.read()
    f.close()
    
    opts.decodeFromString(txt)
    return opts

# unpack data that has been packing into 64-bit values that the FPGA FIFO returns
# this has advantage of moving unpacking to the processig loop, which for a multiprocess implementaiton potentially increases speed
def unpackData(packed_data, isMagOnly=True, is16b=True, isPolar=True):
    shp = packed_data.shape
    numSmp = shp[1]
    mag = None
    phase = None
    oct_data = None
    if DebugLog.isLogging:
        DebugLog.log("OCTFPGAProcessingInterface unpackData(): packed_data.shape= %s" % repr(shp))
        
    if isMagOnly:
        mask0 = 2**16 - 1 
        mask1 = mask0 << 16
        mask2 = mask0 << 32
        mask3 = mask0 << 48
        shp1 = shp[1]
        shp1 = numSmp * 4  
        idx0 = range(0, shp1, 4)
        idx1 = range(1, shp1, 4)
        idx2 = range(2, shp1, 4)
        idx3 = range(3, shp1, 4)
        shp =  [shp[0], shp1]
        if DebugLog.isLogging:
            DebugLog.log("OCTProcessing unpackData(): shp= %s" % repr(shp))
            
        mag = np.zeros(shp, dtype=np.int16)
        
        mag[:, idx0] = np.require((packed_data & mask3) >> 48, dtype=np.int16)
        mag[:, idx1] = np.require((packed_data & mask2) >> 32, dtype=np.int16)
        mag[:, idx2] = np.require((packed_data & mask1) >> 16, dtype=np.int16)
        mag[:, idx3] = np.require((packed_data & mask0), dtype=np.int16)
        if DebugLog.isLogging:
            DebugLog.log("OCTProcessing unpackData():  mag min= %g max= %g" % (np.min(mag), np.max(mag)))
    elif is16b:
        mask0 = 2**16 - 1 
        mask1 = mask0 << 16
        mask2 = mask0 << 32
        mask3 = mask0 << 48
        
        shp1 = shp[1]
        shp1 = numSmp * 2
        shp =  [shp[0], shp1, shp[2]]
        idx0 = range(0, shp1, 2)
        idx1 = range(1, shp1, 2)
        if DebugLog.isLogging:
            DebugLog.log("OCTProcessing unpackData(): shp= %s" % repr(shp))
        
        if isPolar:
            mag = np.zeros(shp, dtype=np.int16)
            mag[:, idx0] = np.require((packed_data & mask3) >> 48, dtype=np.int16)
            mag[:, idx1] = np.require((packed_data & mask1) >> 16, dtype=np.int16)
            
            phase = np.zeros(shp, dtype=np.int16)
            phase[:, idx0] = np.require((packed_data & mask2) >> 32, dtype=np.int16)
            phase[:, idx1] = np.require((packed_data & mask0), dtype=np.int16)
        else:
            oct_data = np.zeros(shp, dtype=np.complex)
            d_re = np.require((packed_data & mask3) >> 48, dtype=np.int16)
            d_im = np.require((packed_data & mask2) >> 32, dtype=np.int16)
            oct_data[:, idx0] = d_re + 1j * d_im
            d_re = np.require((packed_data & mask1) >> 16, dtype=np.int16)
            d_im = np.require((packed_data & mask0), dtype=np.int16)
            oct_data[:, idx1] = d_re + 1j * d_im
    else:
        mask0 = 2**32 - 1 
        mask1 = mask0 << 32
        
        oct_data = np.zeros(shp, dtype=np.complex)
        d_re = np.require((packed_data & mask1) >> 32, dtype=np.int32)
        d_im = np.require((packed_data & mask0), dtype=np.int32)
        oct_data[:, :] = d_re + 1j * d_im
        
    return (oct_data, mag, phase)
    
import time

if __name__ == "__main__":
    import msvcrt

    # 50 kHz
    setup = 4
    samplesPerTrig = 1024
    klinROI = [15, 1600]

    # 200 kHz
    setup = 3
    samplesPerTrig = 650
    klinROI = [15, 1175]
    
    sampleOffset = 13
    Ch0Shift = 20
    klinNumPts = 4000

    oct_hw = LV_DLLInterface()
    
    print("Calling InitFPGA()")
    err, klin = oct_hw.InitFPGA(setup, sampleOffset, samplesPerTrig, Ch0Shift, klinNumPts, klinROI)
    print("InitFPGA() err=", err)
    
    # dataToGet = 'interpPD'
    # dataToGet = 'complexFFT'
    dataToGet = 'raw'
    
    try:
        plt.figure(1)
        plt.cla()
        plt.plot(klin)
        magWin = np.hanning(2048)
        phaseCorr = np.zeros(2048)
        err = oct_hw.LoadOCTDispersion(magWin, phaseCorr)
        print("LoadOCTDispersion err=", err)
        for n in range(0, 1):
            if dataToGet == 'interpPD': # interpolated PD
                err, interpPD = oct_hw.AcquireOCTDataInterpPD(10)
                print("AcquireOCTDataInterpPD err=", err)
                if interpPD is not None:
                    interpPDavg = np.mean(interpPD, 0)
                    plt.figure(2)
                    plt.cla()
                    plt.plot(interpPDavg)
                    plt.show()
                    
            elif dataToGet == 'raw': # interpolated PD
                for n in range(0, 5):
                    err, pd_data, mzi_data = oct_hw.AcquireOCTDataRaw(10, samplesPerTrig*2)
                    print("AcquireOCTDataRaw err=", err)
                
                if pd_data is not None:
                    #pd_data = np.mean(pd_data, 0)
                    #mzi_data = np.mean(mzi_data, 0)
                    #mzi_dataavg = np.mean(interpPD, 0)
                    plt.figure(2)
                    plt.clf()
                    plt.subplot(2, 1, 1)
                    plt.plot(pd_data[5, :])
                    plt.subplot(2, 1, 2)
                    plt.plot(mzi_data[1, :])
                    plt.show()
            else:
                err, oct_data = oct_hw.AcquireOCTDataFFT(10, [0, 1023])
                print("AcquireOCTDataFFT err=", err)
                if oct_data is not None:
                    aline = np.abs(oct_data[1, :])
                    aline = 20*np.log10(aline + 1)
                    plt.figure(2)
                    plt.cla()
                    plt.plot(aline)
                    plt.show()
                    
            time.sleep(0.1)
    except Exception as ex:
        raise ex
    finally:
        oct_hw.CloseFPGA()
        