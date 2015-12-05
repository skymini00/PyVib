# -*- coding: utf-8 -*-
"""
Created on Sat Dec  5 10:56:15 2015

@author: jso
"""

def setupZigZagScan(scanParams, mirrorDriver, scanDetails, plotParam, OCTtrigRate):
    """
    This function builds the zigzag scan voltages (scanDetails.mirrOut) 
    and the variables that are used when reformatting the collected A-lines
    into the proper 3D format(scanDetails.c,scanDetails.cTile, and scanDetails.c3)

    Import MEMS parameters from text file and fill variables to be used in creating scan patterns.
    Note: x and y axis definitions reference the output of the MEMS amplifier, not labels in the 
    main python program
    """
    fname='C:\PyOCT\S4128_MEMS_param.txt'
    MEMSparams=np.loadtxt(fname, delimiter=':', skiprows = 5, usecols= (1,)) #read in parameters, see text file for labels
    Vmaxx=MEMSparams[2,] #maximum voltage for MEMS mirror for x-axis
    Vmaxy=MEMSparams[3,] #minimum voltage for MEMS mirror for y-axis
    MaxVinx = (Vmaxx / 7.5) - 9.333 #maximum input voltage to amplifier assuming
    MaxViny = (Vmaxy / 7.5) - 9.333 #maximum input voltage
    LPFcutoff=MEMSparams[10,] #recommended low pass cutoff frequency
    RFx=MEMSparams[6,] #resonant frequency for x-axis 
    RFy=MEMSparams[7,] #resonant frequency for y-axis 

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
 
 


    # System parameters
    Bscan=1 #0 for false and 1 for true
    angle=0 #angle in degrees for the Bscan 
    fs= OCTtrigRate  #laser sweep rate
    fDAQ= mirrorDriver.DAQoutputRate  #DAQ sampling rate
    
    FOV = 10 #field of view in mm
    spatialsampling= 10  #spatial sampling in microns
    Vpmmx=1 #volts per millimeter for x-axis
    Vpmmy=1 #volts per millimeter for y-axis
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
    if Bscan>0:
        angleRad=angle*np.pi/180

        tcycleindexs=np.where(((ts>=(tlower+1/faxis/4))*(ts<=(thigher+1/faxis/4)))>0) #find indicies for cycle starting at 0
        tcycleindexDAQ=np.where(((tDAQ>=(tlower+1/faxis/4))*(tDAQ<=(thigher+1/faxis/4)))>0) #find indicies for cycle starting at 0

        filcycle=xfil[tcycleindexs[0][0]:tcycleindexs[0][-1]]/np.max(xfil[tcycleindexs[0][1]:tcycleindexs[0][-1]]) #filtered quarter cycle to initiate scant at 0 V
        DAQfilcycle=xDAQfil[tcycleindexDAQ[0][0]:tcycleindexDAQ[0][-1]]/np.max(xDAQfil[tcycleindexDAQ[0][0]:tcycleindexDAQ[0][-1]]) #filtered quarter cycle to initiate scant at 0 V
        
        xN=np.cos(angleRad)*filcycle
        yN=np.sin(angleRad)*filcycle
        
        xDAQ=np.cos(angleRad)*DAQfilcycle
        yDAQ=np.sin(angleRad)*DAQfilcycle
        
    else:
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
    
    # calculate the x,y pixel positions for each sampletime
    xNorm1=plotParam.xPixel*((xN/2)+0.5)  
    yNorm1=plotParam.yPixel*((yN/2)+0.5)
    scanDetails.numTrigs=xN.size
    
    # create the 3D array c2, of true/false to give which alines to average for each pixel, and then reformat this into scanDetails.c
    scanDetails.c3=np.zeros((plotParam.xPixel,plotParam.yPixel))
    for i1 in range(0,plotParam.xPixel):  
        for i2 in range(0,plotParam.yPixel):
            distance=(xNorm1-i1)**2+(yNorm1-i2)**2
            scanDetails.c3[i1,i2]=np.argmin(distance)
    scanDetails.c3=np.uint(np.reshape(scanDetails.c3,(plotParam.xPixel*plotParam.yPixel)))
