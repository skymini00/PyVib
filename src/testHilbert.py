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

t=np.arange(0,0.1,1/100000)
x=np.sin(2*np.pi*f*t)

x_h=scipy.signal.hilbert(x)

plt.figure(1)
plt.subplot(3,1,1)
plt.plot(t,x,'b')
plt.subplot(3,1,2)
plt.plot(t,np.abs(x_h))
plt.subplot(3,1,3)
plt.plot(t,np.unwrap(np.angle(x_h)))
