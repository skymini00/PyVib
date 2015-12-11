# PyOCT.py
#  main file and  GUI interface

from ctypes import *
from matplotlib import *

import sys
import time
import datetime
import traceback
import os
import copy
import platform  # for differentiating platforms
import shutil # for copytree

from PyQt4 import QtCore, QtGui, uic
import pyqtgraph as pg

import OCTFPGAProcessingInterface as octfpga
import AudioHardware 
import MirrorDriver
from ROIImageGraphicsView import *
from OCTCommon import *
from OCTProtocolParams import *
from DebugLog import DebugLog
from scipy import stats
from scipy import signal

import BScan
import MScan
import VolumeScan
import Dispersion
import SpeakerCalibration
import JSOraw
import RawDataTest
import SpeakerCalTest

#from EndoSpiralScanProtocol import *
# from ORmicroscopeScanProtocol import *

#from PIL import Image
#from PIL import ImageFont
#from PIL import ImageDraw

form_class = uic.loadUiType(os.path.join("..", "ui", "PyOCT.ui"))[0]                 # Load the UI

class OCTWindowClass(QtGui.QMainWindow, form_class):
    def __init__(self, parent=None):
        print("clrmap length = %d " % (len(ROIImageGraphicsView.COLORMAP_HOT)))
        
        QtGui.QMainWindow.__init__(self, parent)
        try:
            self.setupUi(self)
        except Exception as ex:
            print(format(ex))

        # this is a flag which indicated the initialized failed in some critical wy to prevent program from being run
        # if this is true, application will exit without being shown
        self.initFailed = False
        #QtGui.QMessageBox.critical (self, "Test", "This is a test")

        self.protocolButtons = [self.Mscan_pushButton, self.Bscan_pushButton, self.Volume_pushButton, self.Dispersion_pushButton,
                                self.SpeakerCal_pushButton, self.SpeakerCalTest_pushButton, self.specialScan_pushButton]
        
        for btn in self.protocolButtons:
            btn.setEnabled(False)

        self.enableVolViewer = False
        try:
            import pyqtgraph.opengl as gl
            from OCTGLViewWidget import OCTGLViewWidget
            self.enableVolViewer = True
        except:
            DebugLog.log("OCTWindowClass.__init__: could not import pyqtgraph.opengl: volume viewer is disabled ")
            
        self.isShutdown = False            
        items = ['Single process (slower, simple code)', 'Multiprocess (faster, complex code)']
        itemSelected, okPressed = QtGui.QInputDialog.getItem(self, "QInputDialog.getItem()", "Season:", items, 0, False);
        if not okPressed:
            self.initFailed = True
            return
        else:
            self.multiProcess = (itemSelected == items[1])
        DebugLog.log("OCTWindowClass.__init__: itemSelected= %s okPressed= %s multiProcess= %s" % (repr(itemSelected), repr(okPressed), repr(self.multiProcess)))

        self.saveRaw = False
        self.saveProcessed = False
        self.saveOpts = SaveOpts()
        self.isCollecting = False
        
        sysName = platform.system()
#        if sysName == 'Windows':
#            basePath = "C:\\PyOCT\\"
#        elif sysname == 'Darwin':
#            basePath = "/Applications\PyOCT\\"
#        else:
        basePath = os.path.abspath("..")
            
        self.basePath = basePath
        defaultConfigPath = os.path.join(basePath, 'config', 'defaults')
        configBasePath = os.path.join(basePath, 'config', 'local')
        if not os.path.isdir(configBasePath):
            DebugLog.log("'%s' does not exist, copying from '%s'" % (configBasePath, defaultConfigPath))
            shutil.copytree(defaultConfigPath, configBasePath)
        
        self.configPath = configBasePath
        self.settingsPath = configBasePath

        fpgaOpts = self._readHardwareConfig(configBasePath)
        DebugLog.log(self.octSetupInfo.__dict__)


        self.JSOsaveDispersion_pushButton.setEnabled(False)
        
#        self.dispData = Dispersion.DispersionData(fpgaOpts)        
        # load up stuff for software processing (JSO) routines
        dispLoaded = False
        try:
            JSOraw.loadDispersion_onStartup(self)
            dispLoaded = True
        except Exception as ex:
            traceback.print_exc(file=sys.stdout)
            DebugLog.log('OCTWindowClass.__init__: Error loading dispersion file %s' % self.octSetupInfo.dispFilename)
            # fall back to default dispersion data
            

        self._initOCTHardware(fpgaOpts)
        if not dispLoaded:  # initialize disperino to default given the FPGA options
            fpgaOpts = self.oct_hw.fpgaOpts
            self.dispData = Dispersion.DispersionData(fpgaOpts)  

        imgNorms = self.octSetupInfo.imgNorms
        self.normLow_spinBox.setValue(imgNorms[0])
        self.normLow_slider.setValue(imgNorms[0])
        self.normHigh_spinBox.setValue(imgNorms[1])
        self.normHigh_slider.setValue(imgNorms[1])
        zROI = self.octSetupInfo.ZROI
        self.roiBeginSlider.setValue(zROI[0])
        self.roiEndSlider.setValue(zROI[1])
        self.ZROI_size_spinBox.setValue(zROI[1] - zROI[0] + 1)

        fpgaOpts = self.oct_hw.fpgaOpts
        self.ch0shift_spinBox.setValue(fpgaOpts.Ch0Shift)
        self.sampleOffset_spinBox.setValue(fpgaOpts.SampleOffset)
        
        self.FPGA_FFT16b_checkBox.setChecked(fpgaOpts.FFT16b > 0)
        self.FPGA_polar_checkBox.setChecked(fpgaOpts.Polar > 0)
        self.FPGA_magOnly_checkBox.setChecked(fpgaOpts.MagOnly > 0)
        self.FPGA_postFFT_rescale_spinBox.setValue(fpgaOpts.FFTReImgScale)

        self.isUpdating = False
        
        magWin = self.dispData.magWin
        phCorr = self.dispData.phaseCorr

        if magWin is not None:
            pl = self.plot_disp_magwin
            pl.plot(magWin, pen='b')
        if phCorr is not None:
            pl = self.plot_disp_phasecorr
            pl.plot(phCorr, pen='b')
        
        self.imgDataScanParams = None
        self.imgdata_zROI = None
        self.connect(self, QtCore.SIGNAL('triggered()'), self.closeEvent)
        self._initGraphVars()
        self.widthStepLast = -1
        self.scanParamsQuickSets = []
        self.audioParamsQuickSets = []
        self.scanParamsQuickSetsDict = {}   # mapping from set name to index in scanParamsQuickSet
        self.audioParamsQuickSetsDict = {}  # mapping from set name to index in audioParamsQuickSet
        
        self.ignoreValueChanged = False
        self.loadQuickSets()
        self.bscan_img_gv.mainOCTObj = self

        defaultSaveDir = 'D:\\Data\\OCT'
        self.saveDir_lineEdit.setText(defaultSaveDir)
        
        self.vol_bscan_gv.mainOCTObj = self
        self.vol_plane_proj_pts = [ [ 0, [0, 0], [0, 0] ], [ 0, [0, 0], [0, 0] ] ] 
        
        self.bscan_img_gv.setTransformationAnchor(QtGui.QGraphicsView.AnchorViewCenter)

        self.nextProtocol = None
        self.volDataLast = None
        # Bind the event handlers
        self._bindEventHandlers()
        
        fpgaChkBoxes = [ self.FPGA_keepDataPacked_checkBox, self.FPGA_FFT16b_checkBox, self.FPGA_polar_checkBox, self.FPGA_magOnly_checkBox]
        if self.octSetupInfo.setupNum == 3:
            for chkbox in fpgaChkBoxes:
                chkbox.setChecked(False)
                chkbox.setEnabled(False)
        
        # mag_plt = self.plot_mscan_mag_tuning
        if self.mirrorDriver.mirrorType == MirrorDriver.MirrorType.MEMS_MICROSCOPE:
            from ORuscopeStageControl import Leica
            self.leica = Leica()
            self.focalPlaneAdj = Leica()        
            
        defaultScanQuickSetName = self.octSetupInfo.defaultScanQuickSet
        if defaultScanQuickSetName is not None:
            try:
                idx = self.scanParamsQuickSetsDict[defaultScanQuickSetName]
                scanP = self.scanParamsQuickSets[idx]
                self.loadScanParams(scanP)               
                self.scan_quickSet_comboBox.setCurrentIndex(idx)
            except:
                DebugLog.log("OCTWindowClass: __init__() could not load quickset '%s'" % defaultScanQuickSetName)
            
        
        for btn in self.protocolButtons:
            btn.setEnabled(True)
            
        self.spCal = None
        try:
            filepath = os.path.join(configBasePath, 'speaker_cal_last.pickle')
            spCal = SpeakerCalibration.loadSpeakerCal(filepath)
            self.spCal = spCal
            self.audioHW.loadSpeakerCalFromProcData(spCal)
        except Exception as ex:
            traceback.print_exc(file=sys.stdout)
            DebugLog.log("OCTWindowClass: __init__() failed to load speaekr calibration '%s'" % filepath)
            
        self.klin = None
        DebugLog.log("OCTWindowClass: __init__() done")
        
        
    def _readHardwareConfig(self, configBasePath):
        self.octSetupInfo = OCTSetupInfo()
        if not os.path.exists(configBasePath):
            DebugLog.log("OCTMultiProcInterface.readHardwareConfig(): Config path does not exist, creating")
            try:
                os.makedirs(configBasePath)
                octSetupInfo = OCTSetupInfo()
                writeOCTsetupInfo(configBasePath, octSetupInfo)
            except Exception as ex:
                DebugLog.log("OCTMultiProcInterface.readHardwareConfig(): Critical error - could not create config base path '%s'" % configBasePath)
                traceback.print_exc(file=sys.stdout)
        else:
            try:
                filepath = os.path.join(configBasePath, 'oct.txt')
                self.octSetupInfo = readOCTSetupInfo(filepath)
                self.mirror_label.setText('Mirror driver: ' + self.octSetupInfo.mirrorConfigFile)
            except Exception as ex:
                DebugLog.log("OCTMultiProcInterface.readHardwareConfig(): Could not read OCT setup '%s'" + repr(filepath))
                traceback.print_exc(file=sys.stdout)
                self.octSetupInfo = OCTSetupInfo()

            self.audioHW = AudioHardware.AudioHardware()                
            try:
                audioFile = self.octSetupInfo.audioConfigFile
                filepath = os.path.join(configBasePath, audioFile)
                DebugLog.log("reading audio file '%s'" % filepath)
                self.audioHW = AudioHardware.readAudioHWConfig(filepath)
                DebugLog.log("AudioHardware= " + self.audioHW.encodeToString('\t'))
            except Exception as ex:
                DebugLog.log("OCTMultiProcInterface.readHardwareConfig(): Could not read audio hardwaare")
                traceback.print_exc(file=sys.stdout)
                self.audioHW = AudioHardware.AudioHardware()                
    
            self.mirrorDriver = MirrorDriver.MirrorDriver()
            try:
                mirrorFile = self.octSetupInfo.mirrorConfigFile
                filepath = os.path.join(configBasePath, mirrorFile)
                DebugLog.log("reading mirror file '%s'" % filepath)
                self.mirrorDriver = MirrorDriver.readMirrorDriverConfig(filepath)
                DebugLog.log("OCTMultiProcInterface.readHardwareConfig(): mirrorDriver=\n %s" % repr(self.mirrorDriver))
            except Exception as ex:
                DebugLog.log("OCTMultiProcInterface.readHardwareConfig(): Could not read mirror settings")
                traceback.print_exc(file=sys.stdout)
                self.mirrorDriver = MirrorDriver.MirrorDriver()
            
            fpgaOpts = octfpga.FPGAOpts_t()
            try:
                fpgaOptsFile = self.octSetupInfo.FPGAOptsFile
                filepath = os.path.join(configBasePath, fpgaOptsFile)
                fpgaOpts = octfpga.readFPGAOptsConfig(filepath)
#                octfpga.LV_DLLInterface.fpgaOpts = fpgaOpts
            except Exception as ex:
                DebugLog.log("OCTMultiProcInterface.readHardwareConfig(): Could not read FPGA Opts")
                traceback.print_exc(file=sys.stdout)
                fpgaOpts = octfpga.FPGAOpts_t()
                
            return fpgaOpts
        
    def _initOCTHardware(self, fpgaOpts):
        setup = self.octSetupInfo.setupNum
        
        if self.multiProcess:
            # returns LV_DLLInterface_BGProcess_Adaptor
            oct_hw = octfpga.StartOCTInterfaceBGProcess(self.basePath)  
        else:
            oct_hw = octfpga.LV_DLLInterface()
        
        samplesPerTrig = fpgaOpts.SamplesPerTrig
        sampleOffset = fpgaOpts.SampleOffset
        Ch0Shift = fpgaOpts.Ch0Shift
        klinNumPts =  fpgaOpts.numKlinPts
        klinROI = [fpgaOpts.klinRoiBegin, fpgaOpts.klinRoiEnd]
        DebugLog.log("_initOCTHardware: initializing FPGA interface")
        
        err, klin = oct_hw.InitFPGA(setup, sampleOffset, samplesPerTrig, Ch0Shift, klinNumPts, klinROI)
        DebugLog.log("_initOCTHardware: InitFPGA err= %d" % err)
        dispFilePath = os.path.join(self.settingsPath, "Dispersion")
        dispFilePath = os.path.join(dispFilePath, self.octSetupInfo.dispFilename)
        self.oct_hw = oct_hw
        
        # read in dispersion data
        try:
            # dispCorr = readDispersionFile(dispFilePath)
#            dispData = Dispersion.loadDispData(self, dispFilePath)
#            oct_hw.LoadOCTDispersion(dispData.magWin, -dispData.phaseCorr)
            oct_hw.LoadOCTDispersion(self.dispData.magWin, -self.dispData.phaseCorr)
#            self.dispData = dispData
        except Exception as ex:
            print("Could not load dispersion file '%s'" % dispFilePath)
            traceback.print_exc(file=sys.stdout)
            QtGui.QMessageBox.critical (self, "Could not load dispersion", "Could not load dispersion file '%s' file may be missing of incorrect format " % dispFilePath)
        
        
    def _bindEventHandlers(self):
        # protocol button event handlers
        self.Bscan_pushButton.clicked.connect(self.BScan_clicked)  
        self.Volume_pushButton.clicked.connect(self.VolScan_clicked) 
        self.Mscan_pushButton.clicked.connect(self.MScan_clicked)  
        self.SpeakerCal_pushButton.clicked.connect(self.SpeakerCal_clicked) 
        self.Dispersion_pushButton.clicked.connect(self.Dispersion_clicked) 
        self.SpeakerCalTest_pushButton.clicked.connect(self.SpeakerCalTest_clicked)
        self.specialScan_pushButton.clicked.connect(self.SpecialScan_clicked)  
        self.JSOraw_pushButton.clicked.connect(self.JSOraw_clicked)  
        self.rawDataTest_pushButton.clicked.connect(self.RawDataTest_clicked)
        
        self.focalPlaneAdj_spinbox.valueChanged.connect(self.focalPlaneChanged)
        self.rotZ_dial.valueChanged.connect(self.rotationZDialChanged)

        # scan quick set event handlers
        self.scan_quickSet_save_pushButton.clicked.connect(self.scan_quickSet_save_clicked)
        self.scan_quickSet_reload_pushButton.clicked.connect(self.scan_quickSet_reload_clicked)
        self.scan_quickSet_comboBox.currentIndexChanged.connect(self.scan_quickSet_comboBox_currentIndexChanged)
        
        # sound quick set event handlers
        self.sound_quickSet_save_pushButton.clicked.connect(self.sound_quickSet_save_clicked)
        self.sound_quickSet_reload_pushButton.clicked.connect(self.sound_quickSet_reload_clicked)
        self.sound_quickSet_comboBox.currentIndexChanged.connect(self.sound_quickSet_comboBox_currentIndexChanged)

        self.freqStart_dblSpinBox.valueChanged.connect(self.freqStartEndChanged)
        self.freqEnd_dblSpinBox.valueChanged.connect(self.freqStartEndChanged)
        self.freqSteps_spinBox.valueChanged.connect(self.freqStepsChanged)
        self.freqDelta_dblSpinBox.valueChanged.connect(self.freqDeltaChanged)
        
        self.ampStart_spinBox.valueChanged.connect(self.ampStartEndChanged)
        self.ampEnd_spinBox.valueChanged.connect(self.ampStartEndChanged)
        self.ampSteps_spinBox.valueChanged.connect(self.ampStepsChanged)
        self.ampDelta_spinBox.valueChanged.connect(self.ampDeltaChanged)
        
        # mscan event handlers
        self.mscan_single_pt_button.clicked.connect(self.mscan_single_pt_clicked)
        self.bmscan_box_region_button.clicked.connect(self.bmscan_box_regionn_clicked)
        self.BMscan_resolution_dblSpinBox.valueChanged.connect(self.BMscan_resolution_changed)
        self.BMscan_numSteps_spinBox.valueChanged.connect(self.BMscan_numSteps_changed)
        
        self.saveZROIandNorms_pushButton.clicked.connect(self.saveZROIandNorms)
        self.saveDir_pushButton.clicked.connect(self.saveDir_clicked)
        
        self.audio_loadFlatSpkCal_pushButton.clicked.connect(self.loadFlatSpeakerCal)
        self.JSOsaveDispersion_pushButton.clicked.connect(self.JSOsaveDispersion_pushButton_clicked)         
        self.JSOloadDispersion_pushButton.clicked.connect(self.JSOloadDispersion_pushButton_clicked)         
        
        self.roiBeginSlider.valueChanged.connect(self.ZROIChanged)
        self.roiEndSlider.valueChanged.connect(self.ZROIChanged)
        
        sliders = (self.length_horizontalSlider, self.lengthOffset_horizontalSlider, self.width_horizontalSlider, self.widthOffset_horizontalSlider)
        spinBoxes = (self.length_dblSpinBox, self.lengthOffset_dblSpinBox, self.width_dblSpinBox, self.widthOffset_dblSpinBox)
        
        for i in range(0, 4):
            slider = sliders[i]
            spinBox = spinBoxes[i]
            cb = self.scanParamSliderChanged(spinBox, slider)
            slider.valueChanged.connect(cb)
            cb = self.scanParamSpinBoxChanged(spinBox, slider)
            spinBox.valueChanged.connect(cb)

    def _initGraphVars(self):
        layout = QtGui.QHBoxLayout()
        if self.enableVolViewer:
            import pyqtgraph.opengl as gl
            from OCTGLViewWidget import OCTGLViewWidget
            self.graphicsView = OCTGLViewWidget()
            layout.addWidget(self.graphicsView)
            self.frame.setLayout(layout)
            
            w = gl.GLViewWidget()
            self.glview = w
            
            self.glvolitem = None
            self.glaxitem = None
            
        r = np.linspace(0, 255, 256)
        g = np.zeros(256)
        g[128:] = np.linspace(0, 255, 128)
        b = np.zeros(256)
        b[192:] = np.linspace(0, 255, 64)
        lut = np.zeros((256, 3))
        lut[:, 0] = r
        lut[:, 1] = g
        lut[:, 2] = b
        self.HOT_LUT = lut
            
            
        self.volumeImg = None
        # self.vol_bscan_graphics_scene = QtGui.QGraphicsScene(0, 0, 300, 500)
        
        self.xLblStyle = {'color': '#000', 'font-size': '16pt'}
        self.yLblStyle = {'color': '#000', 'font-size': '16pt'}
        # self.titleLblStyle = {'color': '#000', 'font-size': '24pt'}

        clrArray = []
        clrArray.append(QtGui.QColor(255, 0, 0))     # bright blue
        clrArray.append(QtGui.QColor(0, 0, 255))     # brigh tblue
        clrArray.append(QtGui.QColor(0, 196, 0))     # green
        clrArray.append(QtGui.QColor(255, 128, 0))   # orange 
        clrArray.append(QtGui.QColor(96, 0, 160))    # purple
        clrArray.append(QtGui.QColor(180, 96, 40))   # brown
        clrArray.append(QtGui.QColor(0, 196, 255))   # light blue
        clrArray.append(QtGui.QColor(128, 0, 0))     # dark red
        clrArray.append(QtGui.QColor(196, 196, 0))   # dark yellow
        clrArray.append(QtGui.QColor(0, 128, 0))     # dark green
        clrArray.append(QtGui.QColor(0, 0, 160))     # dark blue
        
        penArray = []
        brushArray = []
        for n in range(0, len(clrArray)):
            penArray.append(QtGui.QPen(clrArray[n]))
            brushArray.append(QtGui.QBrush(clrArray[n]))
        
        self.brushArray = brushArray
        self.penArray = penArray
        
    """
        Call back for scan paratmer sliders
        Outer function returns a closure because callbacks all have same form
        We can just write one and reuse for all parameter types (length, width, offsets)
    """    
    def scanParamSliderChanged(self, spinBox, slider):
        def dimChanged(self):
            spinBox.setValue(slider.value()/1000)
            
        return dimChanged
    
    """
        Call back for scan paratmer spinboxes
        Outer function returns a closure because callbacks all have same form
        We can just write one and reuse for all parameter types (length, width, offsets)
    """    
    def scanParamSpinBoxChanged(self, spinBox, slider):
        def dimChanged(self):
            slider.blockSignals(True)   # block the signals so the scanParam
            slider.setValue(spinBox.value*1000)
            slider.blockSignals(False)
            
        return dimChanged
        
    def focalPlaneChanged(self):
        """
        If the user changes the focal plane of the OCT beam on the OR microscope
        """
        NewPosition = self.focalPlaneAdj_spinbox.value()
        if self.mirrorDriver.mirrorType == MirrorType.OR_MICROSCOPE:
            self.focalPlaneAdj.setPosition(NewPosition)
        
    def loadQuickSets(self):
        """
        load the quick sets         
        """
        
        self.scan_quickSet_comboBox.clear()
        self.sound_quickSet_comboBox.clear()
        
        scanParamsDir = os.path.join(self.settingsPath, 'Scan Params')
        mirrorFile = self.octSetupInfo.mirrorConfigFile
        s = re.split('\.', mirrorFile)
        mirrorFileNoExt = s[0]
        scanParamsDir = os.path.join(scanParamsDir, mirrorFileNoExt)
        self.scanParamsQuickSetsDict = {}
        try:
            fileList = os.listdir(scanParamsDir)
            for fName in fileList:
                filename9, file_extension = os.path.splitext(fName)
                if file_extension=='.pickle':
                    filePath = os.path.join(scanParamsDir, fName)
                    try:
                        f = open(filePath, 'rb')
                        scanP = pickle.load(f)
                        f.close()
                        
                        fParts = re.split('\.', fName)
                        setName = fParts[0]
                        self.scanParamsQuickSets.append(scanP)
                        self.scan_quickSet_comboBox.addItem(setName)
                        idx = len(self.scanParamsQuickSets) - 1
                        self.scanParamsQuickSetsDict[setName] = idx
        
                    except Exception as ex:
                        traceback.print_exc(file=sys.stdout)
                        print("loadQuickSets: could not load scan quick set '%s'" % setName)
        except Exception as e:
            traceback.print_exc(file=sys.stdout)
            print("loadQuickSets: could not load scan quick sets from directory'%s'" % scanParamsDir)
        
                
        audioParamsDir =  os.path.join(self.settingsPath, 'Audio Output Params')
        try:
            fileList = os.listdir(audioParamsDir)
            for fName in fileList:
                filePath = os.path.join(audioParamsDir, fName)
                try:
                    f = open(filePath, 'rb')
                    audioP = pickle.load(f)
                    f.close()
                    
                    fParts = re.split('\.', fName)
                    setName = fParts[0]
                    self.audioParamsQuickSets.append(audioP)
                    self.sound_quickSet_comboBox.addItem(setName)
                    
                    idx = len(self.audioParamsQuickSets) - 1
                    self.audioParamsQuickSetsDict[setName] = idx
                except Exception as ex:
                    traceback.print_exc(file=sys.stdout)
                    print("loadQuickSets: could not load scan quick set '%s'", setName)                
        except Exception as e:
            traceback.print_exc(file=sys.stdout)
            print("loadQuickSets: could not load audio quick sets from directory'%s'" % scanParamsDir)
            
    def closeEvent(self, event): 
        print("OCTWindowClass: closeEvent() entered function") 
        
        if self.mirrorDriver.mirrorType == MirrorDriver.MirrorType.OR_MICROSCOPE:
            self.focalPlaneAdj.closeCOMport()
            
        self.shutdown()
        super().closeEvent(event)
       
    def stopCollection(self):
        self.doneFlag = True
            
    def setProtocolButtonSignalsBlocked(self, blk):
        for btn in self.protocolButtons:
            btn.blockSignals(blk)
            
    # this function should be called when the collection is finished (for single process)
    def finishCollection(self, sendAggDataRequest=True):
        self.isCollecting = False
        nextProtocol = self.nextProtocol
        
        # send one more aggregate data trequest to be sure the collection is complete
        self.setProtocolButtonSignalsBlocked(True)
        
        for btn in self.protocolButtons:
            btn.setChecked(False)
        
        self.doneFlag = False
        
        if nextProtocol is not None:
            if nextProtocol == 'Bscan':
                self.Bscan_pushButton.setChecked(True)
                self.setProtocolButtonSignalsBlocked(False)
                self.protocol = self.bscanPrtcl
                self.tabWidget.setCurrentIndex(0)
                BScan.runBScan(self)
                
            elif nextProtocol == 'Mscan':
                self.Mscan_pushButton.setChecked(True)
                self.setProtocolButtonSignalsBlocked(False)
                self.mscanTuningCurveLast = None
                MScan.runMScan(self)
          
            elif nextProtocol == 'VolScan':
                self.volumeData = None
                self.Volume_pushButton.setChecked(True)
                self.setProtocolButtonSignalsBlocked(False)
                VolumeScan.runVolScan(self)
            elif nextProtocol == 'SpeakerCal':
                self.SpeakerCal_pushButton.setChecked(True)
                self.setProtocolButtonSignalsBlocked(False)
                SpeakerCalibration.runSpeakerCal(self)
            elif nextProtocol == 'Dispersion':
                self.Dispersion_pushButton.setChecked(True)
                self.setProtocolButtonSignalsBlocked(False)
                Dispersion.runDispersion(self)
            elif nextProtocol == 'SpecialScan':
                self.specialScan_pushButton.setChecked(True)
                self.setProtocolButtonSignalsBlocked(False)
                # TODO include speaial scan code here
            elif nextProtocol == 'SpeakerCalTest':
                pass # TODO include speaial scan code here
            elif nextProtocol == 'JSOraw':
                self.JSOraw_pushButton.setChecked(True)
                self.setProtocolButtonSignalsBlocked(False)
                self.protocol = self.JSOrawPrtcl
                self.tabWidget.setCurrentIndex(7)
                JSOraw.runJSOraw(self)
            
        else:
            self.setProtocolButtonSignalsBlocked(False)
            

    def getZROI(self):
        roiBegin = self.roiBeginSlider.value()
        roiEnd = self.roiEndSlider.value()
        if(roiEnd <= roiBegin):
            roiEnd = roiBegin + 1
        roiSize = roiEnd - roiBegin 
        # ensure ROI is aligned to multiple of 4
        # this is required if data is 16-bit, such as for magnitude only collection
        if (roiSize % 4) != 0:
            roiEnd += 4 - roiSize % 4
        zROI = (roiBegin, roiEnd)
        return zROI
   
    def BScan_clicked(self):  # CtoF button event handler
        if self.Bscan_pushButton.isChecked():
            if self.isCollecting:
                self.nextProtocol = 'BScan'
                self.stopCollection()
            else:
                BScan.runBScan(self)
        else:
            self.nextProtocol = None
            self.stopCollection()
    
    def VolScan_clicked(self):  # CtoF button event handler
        if self.Volume_pushButton.isChecked():
            if self.isCollecting:
                self.nextProtocol = 'VolScan'
                self.stopCollection()
            else:
                VolumeScan.runVolScan(self)
        else:
            self.nextProtocol = None
            self.stopCollection()

    def MScan_clicked(self):  # CtoF button event handler
        if self.Mscan_pushButton.isChecked():
            if self.isCollecting:
                self.nextProtocol = 'MScan'
                self.stopCollection()
            else:
                multiProc = self.multiProcess
                MScan.runMScan(self, multiProc)
        else:
            self.nextProtocol = None
            self.stopCollection()
            
    def Dispersion_clicked(self):  # CtoF button event handler
        if self.Dispersion_pushButton.isChecked():
            if self.isCollecting:
                self.nextProtocol = 'Dispersion'
                self.stopCollection()
            else:
                Dispersion.runDispersion(self)
        else:
            self.nextProtocol = None
            self.stopCollection()            
            
    def SpecialScan_clicked(self):  # CtoF button event handler
        if self.SpecialScan_pushButton.isChecked():
            if self.isCollecting:
                self.nextProtocol = 'SpecialScan'
                self.stopCollection()
            else:
                pass  # TODO include special scan code here
                
        else:
            self.nextProtocol = None
            self.stopCollection()    
            
    def SpeakerCal_clicked(self):  # CtoF button event handler
        if self.SpeakerCal_pushButton.isChecked():
            if self.isCollecting:
                self.nextProtocol = 'SpeakerCal'
                self.stopCollection()
            else:
                SpeakerCalibration.runSpeakerCal(self, self.oct_hw.IsDAQTestingMode())
        else:
            self.nextProtocol = None
            self.stopCollection()   

    def SpeakerCalTest_clicked(self):  # CtoF button event handler
        if self.SpeakerCalTest_pushButton.isChecked():
            if self.isCollecting:
                self.nextProtocol = 'SpeakerCalTest'
                self.stopCollection()
            else:
                SpeakerCalTest.runSpeakerCalTest(self)
        else:
            self.nextProtocol = None
            self.stopCollection()   

    def JSOraw_clicked(self):  # CtoF button event handler
        if self.JSOraw_pushButton.isChecked():
            if self.isCollecting:
                self.nextProtocol = 'JSOraw'
                self.stopCollection()
            else:
                JSOraw.runJSOraw(self)
        else:
            self.nextProtocol = None
            self.stopCollection()
 
    def RawDataTest_clicked(self):  # CtoF button event handler
        if self.rawDataTest_pushButton.isChecked():
            if self.isCollecting:
                self.nextProtocol = 'rawDataTest'
                self.stopCollection()
            else:
                RawDataTest.runRawDataTest(self)
        else:
            self.nextProtocol = None
            self.stopCollection()
            
    def JSOsaveDispersion_pushButton_clicked(self):  # CtoF button event handler
        JSOraw.saveDispersion_pushButton_clicked(self)
        
    def JSOloadDispersion_pushButton_clicked(self):  # CtoF button event handler
        self.dispData = JSOraw.DispersionData()             # this should clear the previously-loaded dispersion data         
        JSOraw.loadDispersion_pushButton_clicked(self)
            
    def rotationZDialChanged(self):
        self.rotation_spinBox.setValue(self.rotZ_dial.value())
        
    def getAudioParams(self):
        audioParams = AudioOutputParams()
        freqStart = self.freqStart_dblSpinBox.value()
        freqEnd = self.freqEnd_dblSpinBox.value()
        freqSteps = self.freqSteps_spinBox.value()
        freqSpacing = self.freqSpacing_comboBox.currentIndex() + 1
        ampStart = self.ampStart_spinBox.value()
        ampEnd = self.ampEnd_spinBox.value()
        ampSteps = self.ampSteps_spinBox.value()
        numTrials = self.numTrials_spinBox.value()
        
        audioParams.trialDuration = self.trialDuration_dblSpinBox.value()
        audioParams.stimDuration = self.stimDuration_dblSpinBox.value()
        audioParams.stimOffset = self.stimOffset_dblSpinBox.value()
        audioParams.stimEnvelope = self.stimEnvelope_dblSpinBox.value()
        audioParams.speakerSel = Speaker(self.speakerSel_comboBox.currentIndex() + 1)
        audioParams.stimType = AudioStimType(self.soundStimType_comboBox.currentIndex() + 1)
        audioParams.numTrials = self.numTrials_spinBox.value()
        audioParams.freqSpacing = FreqSpacing(freqSpacing)
        audioParams.numTrials = numTrials
        audioParams.downsample = self.audio_downsample_spinBox.value() - 1
        
        DebugLog.log("getAudioParams(): freqStart= %f freqEnd= %f freqSteps= %d speakerSel= %s" % (freqStart, freqEnd, freqSteps, repr(audioParams.speakerSel)))
        freq = np.linspace(freqStart, freqEnd, freqSteps)
        s = ""
        for n in range(0, freq.shape[0]):
            s = s + " " + repr(freq[n])
        
        DebugLog.log("getAudioParams(): freq== " + s)
        audioParams.freq = np.tile(freq, (2, 1))
        audioParams.amp = np.linspace(ampStart, ampEnd, ampSteps)        
        audioParams.customSoundDir = self.customSoundDir_lineEdit.text()        
        return audioParams
        
        
    def getScanParams(self):
        sl = self.length_dblSpinBox.value()
        sl_off = self.lengthOffset_dblSpinBox.value()
        sw = self.width_dblSpinBox.value()
        sw_off = self.widthOffset_dblSpinBox.value()
        rot_Z = self.rotation_spinBox.value()
        lengthSteps = self.lengthSteps_spinBox.value()
        widthSteps = self.widthSteps_spinBox.value()
        numAvgs = self.averages_spinBox.value()
        scan_pattern = self.scan_pattern_comboBox.currentIndex() + 1
        print("getScanParams() numAvgs= %d pattern= %d" % (numAvgs, scan_pattern))
        
        scanParams = ScanParams()
        scanParams.scan_pattern=scan_pattern
        scanParams.length = sl
        scanParams.lengthOffset = sl_off
        scanParams.lengthSteps = lengthSteps
        scanParams.width = sw
        scanParams.widthOffset = sw_off
        scanParams.rotation_Z = rot_Z
        scanParams.widthSteps = widthSteps
        scanParams.numAverages = numAvgs
        scanParams.pattern = ScanPattern(scan_pattern)
        scanParams.downsample = self.scan_downsample_spinBox.value() - 1
        scanParams.volBscansPerFrame = self.scan_volBscansPerFrame_spinBox.value()
        scanParams.continuousScan = self.continuousVolume_checkBox.isChecked()

        if self.octSetupInfo.mirrorConfigFile[0:4]=='MEMS':
            scanParams.skew = self.mirrorDriver.skew
            scanParams.phaseAdjust = self.mirrorDriver.phaseAdjust
            scanParams.resonantFreq = self.mirrorDriver.resonantFreq
            scanParams.volScanFreq = self.mirrorDriver.volScanFreq
            useGUI=1
            if useGUI==1:
                scanParams.skew = self.skew_dblSpinBox.value()
                scanParams.phaseAdjust = self.scanPhaseAdjust_spinBox.value()
                scanParams.volScanFreq = self.volScanFreq_spinBox.value()
          
        return scanParams

        
    def blockScanParamsSignals(self, blocked):
        self.length_dblSpinBox.blockSignals(blocked)
        self.width_dblSpinBox.blockSignals(blocked)
        self.lengthSteps_spinBox.blockSignals(blocked)
        self.widthSteps_spinBox.blockSignals(blocked)
        self.lengthRes_dblSpinBox.blockSignals(blocked)
        self.widthRes_dblSpinBox.blockSignals(blocked)
        self.skew_dblSpinBox.blockSignals(blocked)
        
    def loadScanParams(self, scanParams):
        self.blockScanParamsSignals(True)
        self.length_dblSpinBox.setValue(scanParams.length)
        self.lengthOffset_dblSpinBox.setValue(scanParams.lengthOffset)
        self.width_dblSpinBox.setValue(scanParams.width)
        self.widthOffset_dblSpinBox.setValue(scanParams.widthOffset)
        self.rotation_spinBox.setValue(scanParams.rotation_Z)
        self.lengthSteps_spinBox.setValue(scanParams.lengthSteps)
        self.widthSteps_spinBox.setValue(scanParams.widthSteps)
        self.averages_spinBox.setValue(scanParams.numAverages)
        self.scan_pattern_comboBox.setCurrentIndex(scanParams.pattern.value - 1)
        self.scan_downsample_spinBox.setValue(scanParams.downsample)
        self.lengthRes_dblSpinBox.setValue(1e3*scanParams.length / scanParams.lengthSteps)
        self.widthRes_dblSpinBox.setValue(1e3*scanParams.width / scanParams.widthSteps)
        
        self.length_horizontalSlider.setValue(int(1000*scanParams.length))
        self.width_horizontalSlider.setValue(int(1000*scanParams.width))
        if hasattr(scanParams, 'skew'):
            self.skew_dblSpinBox.setValue(scanParams.skew)

        if hasattr(scanParams, 'volBscansPerFrame'):
            self.scan_volBscansPerFrame_spinBox.setValue(scanParams.volBscansPerFrame)

        if hasattr(scanParams, 'continuousScan'):            
            self.continuousVolume_checkBox.setChecked(scanParams.continuousScan)
            
        if hasattr(scanParams, 'phaseAdjust'):            
            self.scanPhaseAdjust_spinBox.setValue(scanParams.phaseAdjust)
            
        if hasattr(scanParams, 'volScanFreq'):            
            self.volScanFreq_spinBox.setValue(scanParams.volScanFreq)
            
        self.blockScanParamsSignals(False)
        
    def blockAudioParamsSignals(self, blocked):
        self.freqStart_dblSpinBox.blockSignals(blocked)
        self.freqEnd_dblSpinBox.blockSignals(blocked)
        self.freqSteps_spinBox.blockSignals(blocked)
        self.freqDelta_dblSpinBox.blockSignals(blocked)
        self.ampStart_spinBox.blockSignals(blocked)
        self.ampEnd_spinBox.blockSignals(blocked)
        self.ampSteps_spinBox.blockSignals(blocked)
        self.ampDelta_spinBox.blockSignals(blocked)
        
        
    def loadAudioParams(self, audioParams):
        self.blockAudioParamsSignals(True)
        freq = audioParams.freq[0, :]
        
        fStart = freq[0]
        fEnd = freq[len(freq)-1]
        numF = len(freq)
        fDelta = 0
        if numF > 1:
            fDelta = (fEnd - fStart) / (numF - 1)
        self.freqStart_dblSpinBox.setValue(fStart)
        self.freqEnd_dblSpinBox.setValue(fEnd)
        self.freqSteps_spinBox.setValue(numF)
        self.freqSpacing_comboBox.setCurrentIndex(audioParams.freqSpacing.value-1)
        self.freqDelta_dblSpinBox.setValue(fDelta)
        
        amp = audioParams.amp
        aStart = amp[0]
        aEnd = amp[len(amp)-1]
        numA = len(amp)
        aDelta = 0
        if numA > 1:
            aDelta = (aEnd - aStart) / (numA - 1)
        
        self.ampStart_spinBox.setValue(aStart)
        self.ampEnd_spinBox.setValue(aEnd)
        self.ampSteps_spinBox.setValue(numA)
        self.ampDelta_spinBox.setValue(aDelta)

        self.numTrials_spinBox.setValue(audioParams.numTrials)
        self.trialDuration_dblSpinBox.setValue(audioParams.trialDuration)
        self.stimDuration_dblSpinBox.setValue(audioParams.stimDuration)
        self.stimOffset_dblSpinBox.setValue(audioParams.stimOffset)
        self.stimEnvelope_dblSpinBox.setValue(audioParams.stimEnvelope)
        self.speakerSel_comboBox.setCurrentIndex(audioParams.speakerSel.value - 1)
        self.soundStimType_comboBox.setCurrentIndex(audioParams.stimType.value - 1)
        self.numTrials_spinBox.setValue(audioParams.numTrials)
        self.audio_downsample_spinBox.setValue(audioParams.downsample)

        if hasattr(audioParams, 'customSoundDir'):
            self.customSoundDir_lineEdit.setText(audioParams.customSoundDir)
            
        self.blockAudioParamsSignals(False)

    def scan_quickSet_save_clicked(self):
        scanP = self.getScanParams()
        saveDir = os.path.join(self.settingsPath, 'Scan Params')
        mirrorFile = self.octSetupInfo.mirrorConfigFile
        s = re.split('\.', mirrorFile)
        mirrorFileNoExt = s[0]
        saveDir = os.path.join(saveDir, mirrorFileNoExt)
        curSetName = self.scan_quickSet_comboBox.currentText() 
        setName, ok = QtGui.QInputDialog.getText(self, 'Quickset Name', 'Name:', text=curSetName)
        
        saveSet = False
        newSet = True
        if ok:
            if self.scan_quickSet_comboBox.findText(setName) >= 0:
                newSet = False
                reply = QtGui.QMessageBox.question(self, 'Set already exsists', "Overwrite existing set %s?" % (setName), QtGui.QMessageBox.Yes | QtGui.QMessageBox.No, QtGui.QMessageBox.No)
                if reply == QtGui.QMessageBox.Yes:
                    saveSet = True
            else:
                saveSet = True
                    
        if saveSet:
            filepath = os.path.join(saveDir, '%s.pickle' % setName)
            f = open(filepath, 'wb')
            pickle.dump(scanP, f)
            f.close()
            if(newSet):
                self.scanParamsQuickSets.append(scanP)
                self.scan_quickSet_comboBox.addItem(setName)
                idx = len(self.scanParamsQuickSets) - 1
                self.scanParamsQuickSetsDict[setName] = idx
            else:
                idx = self.scanParamsQuickSetsDict[setName]
                self.scanParamsQuickSets[idx] = scanP
            
        
    def sound_quickSet_save_clicked(self):
        audioP = self.getAudioParams()
        saveDir = os.path.join(self.settingsPath, 'Audio Output Params')
        curSetName = self.sound_quickSet_comboBox.currentText() 
        setName, ok = QtGui.QInputDialog.getText(self, 'Quickset Name', 'Name:', text=curSetName)

        saveSet = False
        newSet = True
        if ok:
            if self.sound_quickSet_comboBox.findText(setName) >= 0:
                newSet = False
                reply = QtGui.QMessageBox.question(self, 'Set already exsists', "Overwrite existing set %s?" % (setName), QtGui.QMessageBox.Yes | QtGui.QMessageBox.No, QtGui.QMessageBox.No)
                if reply == QtGui.QMessageBox.Yes:
                    saveSet = True
            else:
                saveSet = True
            
        if saveSet:
            filepath = os.path.join(saveDir, '%s.pickle' % setName)
            f = open(filepath, 'wb')
            pickle.dump(audioP, f)
            f.close()
            if(newSet):
                self.audioParamsQuickSets.append(audioP)
                self.sound_quickSet_comboBox.addItem(setName)
                idx = len(self.audioParamsQuickSets) - 1
                self.audioParamsQuickSetsDict[setName] = idx
            else:
                idx = self.audioParamsQuickSetsDict[setName]
                self.audioParamsQuickSets[idx] = audioP
        
        
    def scan_quickSet_reload_clicked(self):
        idx = self.scan_quickSet_comboBox.currentIndex()
        scanP = self.scanParamsQuickSets[idx]
        self.loadScanParams(scanP)
    
    def sound_quickSet_reload_clicked(self):
        idx = self.sound_quickSet_comboBox.currentIndex()
        audioP = self.audioParamsQuickSets[idx]
        self.loadAudioParams(audioP)
    
    def scan_quickSet_comboBox_currentIndexChanged(self, idx):
        scanP = self.scanParamsQuickSets[idx]
        print("scan_quickSet_comboBox_currentIndexChanged idx= %d" % (idx))
        self.loadScanParams(scanP)
    
    def sound_quickSet_comboBox_currentIndexChanged(self, idx):
        audioP = self.audioParamsQuickSets[idx]
        self.loadAudioParams(audioP)

    def saveDir_clicked(self):
        caption = "Choose save directory"
        directory = self.saveDir_lineEdit.text()
        newDir = QtGui.QFileDialog.getExistingDirectory (self, caption, directory)
        self.saveDir_lineEdit.setText(newDir)
        
    def mscan_single_pt_clicked(self):
        self.bmscan_box_region_button.setChecked(False)
        drawType = ROIImgViewROIDrawType.NONE
        if self.mscan_single_pt_button.isChecked():
            drawType = ROIImgViewROIDrawType.SINGLE_PT
            
        self.bscan_img_gv.setROIDrawType(drawType)
    
    def bmscan_box_regionn_clicked(self):
        self.mscan_single_pt_button.setChecked(False)
        drawType = ROIImgViewROIDrawType.NONE
        if self.bmscan_box_region_button.isChecked():
            drawType = ROIImgViewROIDrawType.BOX
        
        self.bscan_img_gv.setROIDrawType(drawType)

    def freqStartEndChanged(self):
        freqStart = self.freqStart_dblSpinBox.value()
        freqEnd = self.freqEnd_dblSpinBox.value()
        freqDelta = self.freqDelta_dblSpinBox.value()
            
        if freqDelta > 0:
            freqSteps = int(np.round((freqEnd - freqStart)/freqDelta)) + 1
            self.freqSteps_spinBox.blockSignals(True)
            self.freqSteps_spinBox.setValue(freqSteps)
            self.freqSteps_spinBox.blockSignals(False)
    
    def freqStepsChanged(self):
        freqSteps = self.freqSteps_spinBox.value()
        freqStart = self.freqStart_dblSpinBox.value()
        freqEnd = self.freqEnd_dblSpinBox.value()
        freqDelta = 0 
        if freqSteps > 1:
            freqDelta = (freqEnd - freqStart)/(freqSteps-1)
        if freqDelta > 0:
            self.freqDelta_dblSpinBox.blockSignals(True)
            self.freqDelta_dblSpinBox.setValue(freqDelta)
            self.freqDelta_dblSpinBox.blockSignals(False)
            
    def freqDeltaChanged(self):
        freqStart = self.freqStart_dblSpinBox.value()
        freqEnd = self.freqEnd_dblSpinBox.value()
        freqDelta = self.freqDelta_dblSpinBox.value()
        if freqDelta > 0:
            freqSteps = int(np.round((freqEnd - freqStart)/freqDelta)) + 1
            self.freqSteps_spinBox.blockSignals(True)
            self.freqSteps_spinBox.setValue(freqSteps)
            self.freqSteps_spinBox.blockSignals(False)
    
    def ampStartEndChanged(self):
        print("ampStartEndChanged")
        ampStart = self.ampStart_spinBox.value()
        ampEnd = self.ampEnd_spinBox.value()
        ampDelta = self.ampDelta_spinBox.value()
        if ampDelta > 0:
            ampSteps = int(np.round((ampEnd - ampStart)/ampDelta)) + 1
            self.ampSteps_spinBox.blockSignals(True)
            self.ampSteps_spinBox.setValue(ampSteps)
            self.ampSteps_spinBox.blockSignals(False)
    
    def ampStepsChanged(self):
        print("ampStempsChanged")
        ampStart = self.ampStart_spinBox.value()
        ampEnd = self.ampEnd_spinBox.value()
        ampSteps = self.ampSteps_spinBox.value()
        na = max(1, ampSteps - 1)
        ampDelta = int((ampEnd - ampStart)/na)
        if ampDelta > 0:
            self.ampDelta_spinBox.blockSignals(True)
            self.ampDelta_spinBox.setValue(ampDelta)
            self.ampDelta_spinBox.blockSignals(False)
    
    def ampDeltaChanged(self):
        print("ampDeltaChanged")
        ampStart = self.ampStart_spinBox.value()
        ampEnd = self.ampEnd_spinBox.value()
        ampDelta = self.ampDelta_spinBox.value()
        if ampDelta > 0:
            ampSteps = int(np.round((ampEnd - ampStart)/ampDelta)) + 1
            self.ampSteps_spinBox.blockSignals(True)
            self.ampSteps_spinBox.setValue(ampSteps)
            self.ampSteps_spinBox.blockSignals(False)
            
    def BMscanBoxRegionSet(self, pt1, pt2):
        if self.imgDataScanParams is not None:
            # calculate the # of steps based on desired resolution
            scanP = self.imgDataScanParams
            (w, h) = self.bscan_img_gv.getROIBoxWidthHeight()    
            (img_w, img_h) = self.bscan_img_gv.getImageWidthHeight()
            scanLen = scanP.length*w/img_w
            
            lenRes = 1e-3*self.BMscan_resolution_dblSpinBox.value()
            numSteps = int(np.floor(scanLen/lenRes))
            self.bmscan_box_rgn = ((pt1.x(), pt1.y()), (pt2.x(), pt2.y()))
            self.BMscan_numSteps_spinBox.blockSignals(True)
            self.BMscan_numSteps_spinBox.setValue(numSteps)
            self.BMscan_numSteps_spinBox.blockSignals(False)

    def BMscan_resolution_changed(self):
        lenRes = 1e-3*self.BMscan_resolution_dblSpinBox.value()
        
        if self.imgDataScanParams is not None:
            scanP = self.imgDataScanParams
            (w, h) = self.bscan_img_gv.getROIBoxWidthHeight()    
            (img_w, img_h) = self.bscan_img_gv.getImageWidthHeight()
            scanLen = scanP.length*w/img_w
            
            print("BMscan_resolution_changed lenRes=%0.3f w= %d img_w= %d scanLen=%0.3f" % (lenRes, w, img_w, scanLen))
            numSteps = int(np.round(scanLen/lenRes))
            
            self.BMscan_numSteps_spinBox.blockSignals(True)
            self.BMscan_numSteps_spinBox.setValue(numSteps)
            self.BMscan_numSteps_spinBox.blockSignals(False)
    
    def BMscan_numSteps_changed(self):
        numSteps = self.BMscan_numSteps_spinBox.value()
        
        if self.imgDataScanParams is not None:
            #(scanP, roiBegin, roiEnd, zRoi, zRoiSpread) = self.makeMscanScanParamsAndZROI()
            scanP = self.imgDataScanParams
            (w, h) = self.bscan_img_gv.getROIBoxWidthHeight()    
            (img_w, img_h) = self.bscan_img_gv.getImageWidthHeight()
            scanLen = scanP.length*w/img_w
            lenRes = 1e3*scanLen/numSteps
            
            self.BMscan_resolution_dblSpinBox.blockSignals(True)
            self.BMscan_resolution_dblSpinBox.setValue(lenRes)
            self.BMscan_resolution_dblSpinBox.blockSignals(False)

    def saveZROIandNorms(self):
        nL = self.normLow_slider.value()
        nH = self.normHigh_slider.value()
        
        z1 = self.roiBeginSlider.value()
        z2 = self.roiEndSlider.value()
        
        self.octSetupInfo.ZROI = [z1, z2]
        self.octSetupInfo.imgNorms = [nL, nH]
        writeOCTsetupInfo(self.configPath, self.octSetupInfo)
        
    # display volume inage as 3D in OpenGL context
    def displayVolumeImg3D(self, volImg):
        if not self.enableVolViewer:
            return

#        import pyqtgraph.opengl as gl
#        from OCTGLViewWidget import OCTGLViewWidget            
        view = self.graphicsView
        
        nL = 65535*self.vol_3dnormlow_slider.value() / 100
        nH = 65535*self.vol_3dnormhigh_slider.value() / 100
        
        # convert to 8-bit
        volImg = np.clip(volImg, nL, nH)
        volImg = 65535*(volImg - nL)/(nH - nL)
        volImg = np.floor(255*volImg / 65535)
        volImg = np.require(volImg, dtype=np.uint8)

        # remap image values to RGB colormap ("HOT" colormap)
        shp = volImg.shape
        d2 = numpy.empty(shp + (4,), dtype=numpy.ubyte)
        volImgFlat = np.reshape(volImg, np.prod(shp))    
        clrMap = ROIImageGraphicsView.COLORMAP_VAL_HOT
        r = clrMap[volImgFlat, 0]
        g = clrMap[volImgFlat, 1]
        b = clrMap[volImgFlat, 2]

        d2[..., 0] = r.reshape(shp)
        d2[..., 1] = g.reshape(shp)
        d2[..., 2] = b.reshape(shp)
      #            d2[..., 3] = numpy.exp(numpy.log(255)/(output.max()-output.min())*output) #log
        d2[..., 3] = volImg/volImg.max()*255
      #            d2[..., 3] = log(output)
    
        # color one prt o the image with an "axes
        d2[:, 0, 0] = [255,0,0,100]
        d2[0, :, 0] = [0,255,0,100]
        d2[0, 0, :] = [0,0,255,100]
    
        dataitem = pg.opengl.GLVolumeItem(d2, smooth=False)
        dx = -shp[0]/2
        dy = -shp[1]/2
        dz = -shp[2]/2
        dataitem.translate(dx, dy, dz)
        while view.items:
            view.removeItem(view.items[0])
    
        view.addItem(dataitem)          
        
    def getSaveState(self):
        return self.save_pushButton.isChecked()
        
    def getSaveOpts(self):
        saveOpts = SaveOpts()
        
        saveDir = self.saveDir_lineEdit.text()
        saveOpts.saveBaseDir = saveDir
        saveOpts.subject = self.saveSubject_lineEdit.text()
        saveOpts.dirNameScheme = SaveDirNameScheme.TIMESTAMP_PROTOCOL_SUBJECT        
        saveOpts.saveOnlyMostRecentFrame = self.save_mostRecentFrameOnly_checkBox.isChecked()
        saveOpts.saveRaw = self.save_raw_checkBox.isChecked()
        saveOpts.notes = self.save_notes_plainTextEdit.toPlainText()
        
        return saveOpts

    def mscanSinglePtSet(self, pt):   # this is called by ROIImageGraphicsView
        pass
    
    # load a flat speaker calibration, where 1V output defined to be 80 dB (60 dB = 0.1V)
    def loadFlatSpeakerCal(self):
        self.audioHW.speakerCalFreq = np.array([[20, 100e3], [20, 100e3]]) 
        self.audioHW.speakerCal = np.array([[60, 60], [60, 60]]) 
        
    def ZROIChanged(self):
        roiBegin = self.roiBeginSlider.value()
        roiEnd = self.roiEndSlider.value()
        if(roiEnd <= roiBegin):
            roiEnd = roiBegin + 1
        zroi = (roiBegin, roiEnd)            
        self.zROILast = zroi
        self.ZROI_size_spinBox.setValue(roiEnd - roiBegin + 1)
        self.bscan_img_gv.centerOn(0, 0)
        
#        if not self.singleProcess:
#            self.collMsgQ.put(CollProcMsg(CollProcMsgType.CHANGE_ZROI, zroi))        
            
    def shutdown(self):
        DebugLog.log("OCTWindowClass.shutdown(): closing FPGA")
        try:
            self.oct_hw.CloseFPGA()
            DebugLog.log("OCTWindowClass.shutdown(): shutting down")
            self.oct_hw.Shutdown()
        except Exception as ex:
            DebugLog.log("OCTWindowClass.shutdown(): exception while attempting to close/shutdown FPGA")
            traceback.print_exc(file=sys.stdout)            
            
        self.isShutdown = True
        DebugLog.log("OCTWindowClass.shutdown(): done")

            
    def __del__(self):
        if not self.isShutdown:
            self.shutdown()     
        
if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    myWindow = OCTWindowClass(None)
    if not myWindow.initFailed:
        myWindow.show()
        app.exec_()
