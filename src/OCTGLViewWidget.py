from pyqtgraph.opengl import GLViewWidget
import numpy as np
try:
    from PyQt5 import uic, QtGui, QtWidgets, QtCore
except (RuntimeError, ImportError):
    from PyQt4 import QtGui, QtCore
    from PyQt4.QtCore import pyqtSlot
from OpenGL.GL import *
from OpenGL.GLU import *
import logging
from pyqtgraph import ViewBox
import time


# left mouse and middle mouse buttons do orbit and pan
# right mouse button zooms near clip in and out if ROI is not drawing
# right mouse button draws ROI if ROI is not drawing

class OCTGLViewWidget(GLViewWidget):
    positionSignal = QtCore.Signal(tuple)
    fullSignal = QtCore.Signal(tuple)
    nearClipSignal = QtCore.Signal(float)
    distanceSignal = QtCore.Signal(float)
    elevationSignal = QtCore.Signal(float)
    azimuthSignal = QtCore.Signal(float)
    
    
    def __init__(self, parent=None):
        super().__init__()
        self.drawROI = False
        # self.nearClip = self.opts['distance']*0.5
        self.opts['distance'] = 1000
        self.nearClip = 0.00001
        self._movementFwd = 0
        self._movementRight = 0
        self._flyMode = False
        
        self._camPos = QtGui.QVector3D(-200, -200, -200)
        # camera rotation in pitch/yaw
        self._camRot = (45, 45)
        self._flySpeed = 5
        self._tickLast = np.nan
        # self._camRot = 
        
        
    @pyqtSlot(float)
    def changeNearClip(self, value):
        self.nearClip = value
        self.update()

    @pyqtSlot(float)
    def changeAzimuth(self, value):
        self.opts['azimuth'] = value
        self.update()        

    @pyqtSlot(float)
    def changeElevation(self, value):
        self.opts['elevation'] = value
        self.update()        

    @pyqtSlot(float)
    def changeDistance(self, value):
        self.opts['distance'] = value
        self.update()        
        
    # converts screen coordinates to world coordinates
    def worldCoordinate(self, pos):
        x0, y0, w, h = self.getViewport()
        z = glReadPixels(pos.x(), h-pos.y(), 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT)
        worldCoordinate = gluUnProject(pos.x(), h-pos.y(), z)
        return worldCoordinate

    # can use this to emit all the points on the selection plane
    def mouseReleaseEvent(self, ev):
        worldCoordinate = self.worldCoordinate(ev)
        self.end = ev.pos()
        if self.drawROI:
            for i in range(self.begin.x(), self.end.x()+1):
                for j in range(self.begin.y(), self.end.y()+1):
                    if (i == self.begin.x() or i == self.end.x()
                        or j == self.begin.y() or j == self.end.y()):
                        self.positionSignal.emit(self.worldCoordinate(QtCore.QPoint(i,j)))
                    self.fullSignal.emit(self.worldCoordinate(QtCore.QPoint(i,j)))
                    
    def tick():
        t1 = time.time()

        if self._tickLast == np.nan:
            self._tickLast = t1
            
        t2 = self._tickLast
        if self.flyMode:
            # make camera look direction vector
            a1 = self._camRot[0] * np.pi/180
            a2 = self._camRot[1] * np.pi/180
            x = np.cos(a1) * np.sin(a2)
            y = np.sin(a1)
            z = np.cos(a1) * np.cos(a2)
            
            dt = t2 - t1
            dr = self._flySpeed * dt * self._movementFwd
            
            # move in the direction of the camera look direction
            self._camPos.x += x * dr
            self._camPos.y += y * dr
            self._camPos.z += z * dr

            # the vector to the right of the look vector            
            a = a2- np.pi()/2
            x = np.sin(a)
            y = 0
            z = np.cos(a)
            
            dr = self._flySpeed * dt * self._movementRight
            self._camPos.x += x * dr
            self._camPos.y += y * dr
            self._camPos.z += z * dr
            
        self.tickLast = time.time()
        
    def mousePressEvent(self, ev):
        super().mousePressEvent(ev)
        if self.drawROI and ev.buttons() == QtCore.Qt.RightButton:
            self.begin = ev.pos()
        
    def mouseMoveEvent(self, ev):
        if self._flyMode:
            # make vector that ponts to where camera is looking 
            a1 = self._camRot[0] * np.pi/180
            a2 = self._camRot[1]* np.pi/180
            x = np.cos(a1) * np.sin(a2)
            y = np.sin(a1)
            z = np.cos(a1) * np.cos(a2)
            v1 = QtGui.QVector3D(x, y, z)
            
            a = a2- np.pi()/2
            x = np.sin(a)
            y = 0
            z = np.cos(a)
            rgt = QtGui.QVector3D(x, y, z)
            
            up = QtGui.QVector3D.crossProduct(rgt, v1)
            
            tr = QtGui.QMatrix4x4()
            tr.lookAt(self._camPos, v1, up)
            # v QtGui.QVector3D(x, y, z)
            
#            tr.setToIdentity()
#            
#            tr.rotate(self._camPos[0], 0, 0, 1)
#            tr.rotate(self._camPos[0], 0, 0, 1)
#            tr.translate()
        else:
            if (ev.buttons() == QtCore.Qt.LeftButton
                or ev.buttons() == QtCore.Qt.MidButton):
                super().mouseMoveEvent(ev)
                self.distanceSignal.emit(self.opts['distance'])
                self.elevationSignal.emit(self.opts['elevation'])
                self.azimuthSignal.emit(self.opts['azimuth'])
                
            elif ev.buttons() == QtCore.Qt.RightButton:
                diff = ev.pos() - self.mousePos
                self.mousePos = ev.pos()
                if not self.drawROI:
                    delta = (diff.x()+diff.y())*1
                    self.nearClip += delta
                    self.nearClip = max(0.01, self.nearClip)
                    self.setProjection()
                    self.nearClipSignal.emit(self.nearClip)
        self.update()

    def wheelEvent(self, ev):
        before = self.opts['distance'] 
        self.setProjection()
        super().wheelEvent(ev)
        after = self.opts['distance']
        self.nearClip += after-before
        self.nearClipSignal.emit(self.nearClip)
        self.distanceSignal.emit(after)
        
    def setFlyMode(self, flyMode):
        self._flyMode = flyMode
        if flyMode:
            self.grabMouse()
            self.grabKeyboard()
        else:
            self.releaseMouse()
            self.releaseKeyboard()
                
    def keyPressEvent(keyEvt):
        k = keyEvt.key()
        if self._flyMode:
            if k == QtCore.Qt.Key_W:  
                self._movementFwd = 1
            elif k == QtCore.Qt.Key_S:
                self._movementFwd = -1
            elif k == QtCore.Qt.Key_A:
                self._movementRight = -1
            elif k == QtCore.Qt.Key_D:
                self._movementRight = 1
            elif k == QtCore.Qt.Key_Escape:
                self.setFlyMode(False)

            
    def keyReleaseEvent(keyEvt):
        k = keyEvt.key()
        if self._flyMode:
            if k == QtCore.Qt.Key_W:  
                self._movementFwd = 0
            elif k == QtCore.Qt.Key_S:
                self._movementFwd = 0
            elif k == QtCore.Qt.Key_A:
                self._movementRight = 0
            elif k == QtCore.Qt.Key_D:
                self._movementRight = 0


#    def projectionMatrix(self, region=None):
#        # Xw = (Xnd + 1) * width/2 + X
#        if region is None:
#            region = (0, 0, self.width(), self.height())
#        
#        x0, y0, w, h = self.getViewport()
#        dist = self.opts['distance']
#        fov = self.opts['fov']
#        nearClip = self.nearClip
#        farClip = dist * 1000.
#
#        r = nearClip * numpy.tan(fov * 0.5 * numpy.pi / 180.)
#        t = r * h / w
#
#        # convert screen coordinates (region) to normalized device coordinates
#        # Xnd = (Xw - X0) * 2/width - 1
#        ## Note that X0 and width in these equations must be the values used in viewport
#        left  = r * ((region[0]-x0) * (2.0/w) - 1)
#        right = r * ((region[0]+region[2]-x0) * (2.0/w) - 1)
#        bottom = t * ((region[1]-y0) * (2.0/h) - 1)
#        top    = t * ((region[1]+region[3]-y0) * (2.0/h) - 1)
#
#        tr = QtGui.QMatrix4x4()
#        tr.frustum(left, right, bottom, top, nearClip, farClip)
#        return tr




