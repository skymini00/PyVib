# -*- coding: utf-8 -*-
"""
Created on Tue Jul 14 03:18:23 2015

@author: OHNS
"""

class DebugLog:
    isLogging = True
    
    def log(msg):
        if (DebugLog.isLogging):
            print(msg)
        
if __name__ == "__main__":
    print("isLogging = %s" % repr(DebugLog.isLogging))
