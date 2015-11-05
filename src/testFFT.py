# -*- coding: utf-8 -*-
"""
Created on Wed Nov  4 20:39:28 2015

@author: jso
"""

import numpy as np
import matplotlib.pyplot as plt
samplingrate=1000
t = np.arange(0,1.5,1/samplingrate)
x=3*np.sin(2*np.pi*150*t)
sp=np.abs(np.fft.rfft(x))/(t.shape[-1]/2)
freq = np.fft.rfftfreq(t.shape[-1])*samplingrate

plt.subplot(2,1,1)
plt.plot(t,x)
plt.subplot(2,1,2)
plt.plot(freq,sp)
print(freq.shape,sp.shape)
plt.show()