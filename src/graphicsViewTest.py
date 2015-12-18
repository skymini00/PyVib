# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 16:04:33 2015

@author: OHNS
"""

from PyQt4 import QtCore, QtGui, uic
import sys
import scipy.misc
from ROIImageGraphicsView import *
import numpy as np

form_class = uic.loadUiType("..\\ui\\graphicsViewTest.ui")[0]                 # Load the UI

imgFile = "..\\exampledata\\mouse organ of corti LV.png"
imgData = scipy.misc.imread(imgFile)
print("imgData.shape= %s dtype= %s" % (repr(imgData.shape), repr(imgData.dtype)))
print("imgData max= %g min= %g" % (np.max(imgData), np.min(imgData)))

img2 = (255 << 24) + (imgData[:, :, 0] << 16)  + (imgData[:, :, 1] << 8) + imgData[:, :, 2]

# img2 = np.transpose(img2)
img2 = np.require(img2, np.uint32, 'C')
print("img2.shape= " + repr(img2.shape))

# shp - imgData.shape[0:1]
img = np.mean(imgData, 2)
img = np.round((np.clip(img, 0, 255)))
img = np.require(img, np.uint8, 'C')
#img = np.require(img, np.uint32, 'C')

imgData = np.require(imgData, np.uint8, 'C')

class GraphicsViewTestWindowClass(QtGui.QMainWindow, form_class):
    def __init__(self, parent=None):
        QtGui.QMainWindow.__init__(self, parent)
        
        try:
            self.setupUi(self)
        except Exception as ex:
            print(format(ex)) 

        print("img min= %g max= %g" % (np.min(img), np.max(img)))
        print("img.shape= " + repr(img.shape))
        
        # qimg = QtGui.QImage(w, h, QtGui.QImage.Format_Indexed8)
        # pilImg = PIL.Image.fromarray(img, 'P')
        # ch = ctypes.c_char.from_buffer(img2, 0)
        #qImg = QtGui.QImage(img2.data, h, w, 4*h, QtGui.QImage.Format_ARGB32)
            
        #imgQ = ImageQt.ImageQt(imgData)
        #pixMap = QtGui.QPixmap.fromImage(imgQ)
        #QtGui.QColor()
        #QtGui.QPen()
        # qpixMap = QtGui.QPixmap(imgFile)
            
        #clrMap = ROIImageGraphicsView.COLORMAP_HOT
        clrMap = ROIImageGraphicsView.COLORMAP_HOT
        self.ROIgraphicsView.setImage(img, clrMap)

        

        
        #scene = self.ROIgraphicsView.scene()
        #pixMapItem = scene.addPixmap(qPixMap)
        #pixMapItem.setZValue(-10)

        self.BoxROI_button.clicked.connect(self.BoxROI_button_clicked)
        self.SinglePt_button.clicked.connect(self.SinglePt_button_clicked)
        self.multiPt_button.clicked.connect(self.multiPt_button_clicked)
        self.polygon_button.clicked.connect(self.polygon_button_clicked)
        self.free_button.clicked.connect(self.free_button_clicked)
        self.clear_polygon_button.clicked.connect(self.clear_polygon_button_clicked)
        self.delete_last_button.clicked.connect(self.delete_last_button_clicked)
        self.delete_all_button.clicked.connect(self.delete_all_button_clicked)
        #self.view.fitInView(QRectF(02/        0, 0, w, h), Qt.KeepAspectRatio)
        self.multiPt_new_button.clicked.connect(self.multiPt_new_button_clicked)
        self.multiPt_move_button.clicked.connect(self.multiPt_move_button_clicked)
        
        self.clear_free_button.clicked.connect(self.clear_free_button_clicked)
        #self.graphicsView.setInteractive(True)
        #self.graphicsView.setDragMode(QtGui.QGraphicsView.RubberBandDrag)
        
        self.ROIgraphicsView.setCursorsVisible(True)
        
#        scene = self.ROIgraphicsView.scene()
#        poly =  QtGui.QPolygonF()
#        poly_item = scene.addPolygon(poly, pen=self.ROIgraphicsView.linePen)
#        poly.append(QtCore.QPointF(-100, -100))
#        poly.append(QtCore.QPointF(100, -100)) 
#        poly.append(QtCore.QPointF(100, 100))
#        poly.append(QtCore.QPointF(-100, 100))
#        scene.removeItem(poly_item)
#        poly_item = scene.addPolygon(poly, pen=self.ROIgraphicsView.linePen)        
#        poly_item.setVisible(True)

    
    def BoxROI_button_clicked(self):
        if self.BoxROI_button.isChecked():
            self.SinglePt_button.setChecked(False)
            self.multiPt_button.setChecked(False)
            self.polygon_button.setChecked(False)
            self.free_button.setChecked(False)
            self.ROIgraphicsView.setROIDrawType(ROIImgViewROIDrawType.BOX)
        else:
            self.ROIgraphicsView.setROIDrawType(ROIImgViewROIDrawType.NONE)
        
    def SinglePt_button_clicked(self):
        if self.BoxROI_button.isChecked():
            self.BoxROI_button.setChecked(False)
            self.multiPt_button.setChecked(False)
            self.polygon_button.setChecked(False)
            self.free_button.setChecked(False)
            self.ROIgraphicsView.setROIDrawType(ROIImgViewROIDrawType.SINGLE_PT)
        else:
            self.ROIgraphicsView.setROIDrawType(ROIImgViewROIDrawType.NONE)
    
    def multiPt_button_clicked(self):
        if self.multiPt_button.isChecked():
            self.SinglePt_button.setChecked(False)
            self.BoxROI_button.setChecked(False)
            self.polygon_button.setChecked(False)
            self.free_button.setChecked(False)
            self.ROIgraphicsView.setROIDrawType(ROIImgViewROIDrawType.MULTI_PT)
        else:
            self.ROIgraphicsView.setROIDrawType(ROIImgViewROIDrawType.NONE)
    
    def polygon_button_clicked(self):
        if self.polygon_button.isChecked():
            self.SinglePt_button.setChecked(False)
            self.BoxROI_button.setChecked(False)
            self.multiPt_button.setChecked(False)
            self.free_button.setChecked(False)
            self.ROIgraphicsView.setROIDrawType(ROIImgViewROIDrawType.POLYGON)
        else:
            self.ROIgraphicsView.setROIDrawType(ROIImgViewROIDrawType.NONE)
            
    def free_button_clicked(self):
        if self.free_button.isChecked():
            self.SinglePt_button.setChecked(False)
            self.BoxROI_button.setChecked(False)
            self.multiPt_button.setChecked(False)
            self.polygon_button.setChecked(False)
            self.ROIgraphicsView.setROIDrawType(ROIImgViewROIDrawType.FREE)
        else:
            self.ROIgraphicsView.setROIDrawType(ROIImgViewROIDrawType.NONE)
            
    def clear_polygon_button_clicked(self):
        self.ROIgraphicsView.clearPolygon()
        
    def delete_last_button_clicked(self):
        self.ROIgraphicsView.deleteLastMultiPt()
    
    def delete_all_button_clicked(self):
        self.ROIgraphicsView.clearMultiPt()
        
    def multiPt_new_button_clicked(self):
        if self.multiPt_new_button.isChecked():
            self.ROIgraphicsView.multPtMode = ROIImgViewMultPtMode.NEW
            self.multiPt_move_button.setChecked(False)
        else:
            self.multiPt_move_button.setChecked(True)
            self.ROIgraphicsView.multPtMode = ROIImgViewMultPtMode.MOVE
        
    def multiPt_move_button_clicked(self):
        if self.multiPt_move_button.isChecked():
            self.ROIgraphicsView.multPtMode = ROIImgViewMultPtMode.MOVE
            self.multiPt_new_button.setChecked(False)
        else:
            self.ROIgraphicsView.multPtMode = ROIImgViewMultPtMode.NEW
            self.multiPt_new_button.setChecked(True)
                    
    def clear_free_button_clicked(self):
        self.ROIgraphicsView.clearFreeDrawROI()
        
if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    myWindow = GraphicsViewTestWindowClass(None)
    myWindow.show()
    app.exec_()