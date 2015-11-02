# -*- coding: utf-8 -*-
"""
Created on Sun Nov  1 16:00:33 2015

@author: jso
"""

# test of Hilbert transform
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal

f=6000
fs=100000

t=np.arange(0,0.1,1/fs)
x=np.sin(2*np.pi*f*t)

x_h=scipy.signal.hilbert(x)

plt.figure(1)
plt.subplot(4,1,1)
plt.plot(t,x,'b')
plt.subplot(4,1,2)
plt.plot(t,np.abs(x_h))
plt.subplot(4,1,3)
plt.plot(t,np.unwrap(np.angle(x_h)))

instfreq = fs/(2*np.pi)*np.diff(np.unwrap(np.angle(x_h)))
#plt.subplot(4,1,4)
#plt.figure(2)
#plt.plot(instfreq)
print(instfreq[0:20])