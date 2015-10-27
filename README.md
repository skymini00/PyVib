# PyVib

Software for special experiments

Requirements: 
  pyqt5
  pyqtgraph
  
Optional:
  pydaqmx (if actually acquiring data)
  pyopengl 3.1 (for volume viewer)
  
  
## Directory organization
```
    /config 
      oct.txt  main config file, the line 'setupNum=' determines which setup to use.  A value of -1 indicates no hardware present
        a value of -2 indicates DAQ hardware, but no FPGA
        a value of 3 is for 200 kHz setup, and 4 for 50 kHz setup.  0 through 2 are for older configs currently untested
      *.txt  config files for audio, mirror, and FPGA, as indicated in oct.txt
      /Audio Output Params   audio params quick sets are stored here
      /Scan Params  scan params quick sets are stored here
      /Dispersion  stores saved disperiosn
      /Speaker Calibration  stores saved speaker calibration files
      
    /src  directory for all python source. 
      PyOCT.py - main program file and GUI interface - calls runXXXXX.py routines for BScan, VolumeScan, etc. when user clics on a button
      OCTFPGAProcessingInterface - interface to the FPGA
      OCTCommon.py - various functions used in many parts of program
      OCTProtocolParams.py - contains ScanParams and AudioParams classes, which are frequently used in protocols
      BScan.py - Bscan protocol
      MScan.py - Mscan Protocol
      VolumeScan.py - Volume scan protocol
      Dispersion.py - Dispersion
      SpeakerCalibration.py - speaker calibraion
      DAQHardware.py  - calss provides conveient wrapper for NIDaqMX functions
      AudioHardware - class for audio hardware settings
      MirrorDriver - class for mirror hardware settings
      tiffile - used to write out TIFF images
      
    /dll - directory containig compiled LabVIEW DLLs which interface with FPGA to aquire data
    /exampledata - directory containing data used in test mode
    /ui - .ui files built in QtDesigner, loaded by program
```
  
  



