# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 15:15:30 2015

@author: OHNS
"""
from PyQt4 import QtCore, QtGui, uic
import numpy as np
import pylab
import copy
from DebugLog import DebugLog

class ROIImageGraphicsView(QtGui.QGraphicsView):
    clrMap=[]
    r = np.linspace(0, 255, 256)
    g = np.zeros(256)
    g[128:] = np.linspace(0, 255, 128)
    b = np.zeros(256)
    b[192:] = np.linspace(0, 255, 64)
    clrMapVal = np.zeros((256, 4))
    clrMapVal[:, 0] = r
    clrMapVal[:, 1] = g
    clrMapVal[:, 2] = b
    clrMapVal[:, 3] = 255
    for i in range(256): 
        clrMap.append(QtGui.qRgb(r[i], g[i], b[i]))
        
    COLORMAP_HOT = clrMap    
    COLORMAP_VAL_HOT = clrMapVal
    
    clrMap=[]
    r = np.linspace(0, 255, 256)
    g = np.linspace(0, 255, 256)
    b = np.linspace(0, 255, 256)
    clrMapVal = np.zeros((256, 4))
    clrMapVal[:, 0] = r
    clrMapVal[:, 1] = g
    clrMapVal[:, 2] = b
    clrMapVal[:, 3] = 255
    for i in range(256): 
        clrMap.append(QtGui.qRgb(r[i], g[i], b[i]))
        
    COLORMAP_GRAY = clrMap
    COLORMAP_VAL_GRAY = clrMapVal
    
    clrMap=[]
    clrMapVal = np.zeros((256, 4))    
    clrMapVal[:, 3] = 255
    for i in range(256): 
        clr = pylab.cm.jet(i)
        clrMapVal[i, :] = clr
        clrMap.append(QtGui.qRgb(255*clr[0], 255*clr[1], 255*clr[2]))
        
    COLORMAP_JET = clrMap
    COLORMAP_VAL_JET = clrMapVal
    
    clrMap=[]
    for i in range(256): 
        clr = 255*pylab.cm.hsv(i)
        clrMap.append(QtGui.qRgb(255*clr[0], 255*clr[1], 255*clr[2]))
        
    COLORMAP_HSV = clrMap
    
    
    clrMap=[]
    for i in range(256): 
        if i < 128:# black and white between 0 and 127
            clr = (i/127, i/127, i/127)
        else:  # HSV between 128 and 255
            clr = pylab.cm.jet((i-128)*2)
        clrMap.append(QtGui.qRgb(255*clr[0], 255*clr[1], 255*clr[2]))
        
    COLORMAP_JET_BW = clrMap
    
    clrMap=[]
    for i in range(256): 
        if i < 128:  # black and white between 0 and 127
            clr = (i/127, i/127, i/127)
        else:  # HSV between 128 and 255
            clr = pylab.cm.hsv((i-128)*2)
        clrMap.append(QtGui.qRgb(255*clr[0], 255*clr[1], 255*clr[2]))
        
    COLORMAP_HSV_BW = clrMap
    
            
    def __init__(self, parent=None):
        QtGui.QGraphicsView.__init__(self, parent)
        self.isMouseDrag = False
        self.ROIBox_pt1 = None   # ROI box upper left
        self.ROIBox_pt2 = None   # ROI box lower right

        self.ROIbox_upperLeft_y = -1
        scene = QtGui.QGraphicsScene(self)
        
        clr = QtGui.QColor(0, 255, 64)
        pen = QtGui.QPen(clr)
        self.txtBrush = QtGui.QBrush(clr)
        self.txtClr = clr

        
        rectItem = scene.addRect(20, 20, 100, 100)
        rectItem.setPen(pen)
        rectItem.setZValue(-5)
        rectItem.setVisible(False)
        self.rectItem = rectItem
        
        self.imgData = None
        self.qImage = None
        self.qPixmap = None
        self.pixMapItem = None
        
        d = 5
        self.xhair_d = d
        x1 = 150 - d
        x2 = 150 + d
        y1 = x1
        y2 = x2
        self.xhairItem1 = scene.addLine(x1, y1, x2, y2, pen)
        self.xhairItem2 = scene.addLine(x1, y2, x2, y1, pen)
        
        self.xhairItem1.setVisible(False)
        self.xhairItem2.setVisible(False)
        
        self.ptPen = QtGui.QPen(pen)

        
        # loat x, float y, float w, float h, QPen pen = QPen(), QBrush brush = QBrush()
        d = 20
        r = d / 2
        x1 = 0
        x2 = 200
        y1 = 0
        y2 = 100
        self.circR = r
        pen.setWidth(4)
        self.linePen = QtGui.QPen(pen)
        
        self.line1_cap_start = CursorCapItem()
        self.line1_cap_start.setRect(x1-r, y1-r, d, d)
        self.line1_cap_start.setPen(pen)
        self.line1_cap_start.ROIgraphicsView = self
        
        self.line1_cap_end = CursorCapItem()
        self.line1_cap_end.isStart = False
        self.line1_cap_end.setRect(x2-r, y2-r, d, d)
        self.line1_cap_end.setPen(pen)
        self.line1_cap_end.ROIgraphicsView = self
        
        scene.addItem(self.line1_cap_start)
        scene.addItem(self.line1_cap_end)
        
        
        #self.line1 = scene.addLine(x1, y1, x2, y2, pen)
        self.line1 = CursorItem()
        self.line1.setLine(x1, y1, x2, y2)
        self.line1.setPen(pen)
        self.line1.ROIgraphicsView = self
        self.line1.setCapItems(self.line1_cap_start, self.line1_cap_end)
        self.line1_cap_start.setCursorItem(self.line1)
        self.line1_cap_end.setCursorItem(self.line1)
        
        scene.addItem(self.line1)
        
        
        clr = QtGui.QColor(255, 255, 0)
        pen.setColor(clr)
        self.line2_cap_start = CursorCapItem()
        self.line2_cap_start.setRect(x1-r, y1-r, d, d)
        self.line2_cap_start.setPen(pen)
        self.line2_cap_start.ROIgraphicsView = self
        
        self.line2_cap_end = CursorCapItem()
        self.line2_cap_end.isStart = False
        self.line2_cap_end.setRect(x2-r, y2-r, d, d)
        self.line2_cap_end.setPen(pen)
        self.line2_cap_end.ROIgraphicsView = self
        
        scene.addItem(self.line2_cap_start)
        scene.addItem(self.line2_cap_end)
        
        #self.line1 = scene.addLine(x1, y1, x2, y2, pen)
        self.line2 = CursorItem()
        self.line2.setLine(x1, y1, x2, y2)
        self.line2.setPen(pen)
        self.line2.ROIgraphicsView = self
        self.line2.setCapItems(self.line2_cap_start, self.line2_cap_end)
        self.line2_cap_start.setCursorItem(self.line2)
        self.line2_cap_end.setCursorItem(self.line2)

        
        scene.addItem(self.line2)
        
        # link the two lines so that moving the caps will move the other cap
        self.line1.linkedLine = self.line2
        self.line1.cursorNum = 0       
        self.line1_cap_start.cursorNum = 0
        self.line1_cap_end.cursorNum = 0
        self.line2.linkedLine = self.line1
        self.line2.cursorNum = 1
        self.line2_cap_start.cursorNum = 1
        self.line2_cap_end.cursorNum = 1
        
        self.setCursorsVisible(False)
        
        self.setScene(scene)
        

        self._ROIdragMode = ROIImgViewdDragMode.NONE
        self._ROIdrawType = ROIImgViewROIDrawType.NONE
        self.ptsList = []
        self.dragLastPt = (-1, -1)
        self.dragLastPtF = None
        self.cursorItemSelected = None
        
        # array of lines and text tiems for multi point 
        self.multiPt_grItems = []
        self.multPtMode = ROIImgViewMultPtMode.NEW
        
        self.ptFont = QtGui.QFont("Arial", 10)
        
        self.singlePt = (-1, -1)
        
        
        self.ROI_poly =  QtGui.QPolygonF()
        self.ROI_poly_item = scene.addPolygon(self.ROI_poly, pen=self.linePen)
        
        self.freedraw_poly =  QtGui.QPolygonF()
        self.freedraw_poly_item = scene.addPolygon(self.freedraw_poly, pen=self.linePen)
        self.freedraw_poly_item.setVisible(True)

        # reference to the main OCT object, used for calling funcitons when certain events occur (such as the defining box region)
        self.mainOCTObj = None
        
        #self.setAlignment(QtCore.Qt.AlignHCenter, QtCore.Qt.AlignVCenter)
        self.setAlignment(QtCore.Qt.AlignHCenter | QtCore.Qt.AlignVCenter)
        
    def mouseMoveEvent(self, mouseEvent):
        pt = (mouseEvent.x(), mouseEvent.y())
        qpt = QtCore.QPoint(pt[0], pt[1])
        ptf = self.mapToScene(qpt)
        pt = (ptf.x(), ptf.y())
        print("ROIImageGraphicsView.mouseMoveEvent x= %d y= %d button= %s ROIdragMode+ %s" % (pt[0], pt[1], repr(mouseEvent.button()), repr(self._ROIdragMode)))

        if self.isMouseDrag:
            if self._ROIdragMode == ROIImgViewdDragMode.DRAWROI:
                if self._ROIdrawType == ROIImgViewROIDrawType.FREE:
                    qptf = QtCore.QPointF(pt[0], pt[1])
                    if qptf != self.freedraw_poly.last():
                        self.freedraw_poly.append(qptf)
                        scene = self.scene()
                        scene.removeItem(self.freedraw_poly_item)
                        self.freedraw_poly_item = scene.addPolygon(self.freedraw_poly, pen=self.linePen)
                    print("ROIImageGraphicsView.mouseMoveEvent: FREE draw numpts= %d " % self.freedraw_poly.count())
                else:
                    self.ROIBox_pt2 = pt 
                    ul = self.ROIBox_pt1
                    lr = self.ROIBox_pt2
                    x = min((ul[0], lr[0]))
                    y = min((ul[1], lr[1]))
                    w = abs(lr[0] - ul[0])
                    h = abs(lr[1] - ul[1])
                    self.rectItem.setRect(x, y, w, h)
                

                
            elif self._ROIdragMode == ROIImgViewdDragMode.MOVE_ITEM:
                # pos = self.itemSelected.pos()
                #pos.setY(ptf.y()) 
                self.cursorItemSelected.moveEvent(ptf, self.dragLastPtF) 
                
                
            elif self._ROIdragMode == ROIImgViewdDragMode.MOVE:
                (dx, dy) = (self.dragLastPt[0] - pt[0], self.dragLastPt[1] - pt[1])
                print("\tdx= %d dy= %d " % (dx, dy))
                #self.translate(dx, dy)
                self.scrollContentsBy(dx, dy)
                # print("(dx, dy) = (%f, %f)" % (dx, dy))

                self.resetCachedContent()            
                self.update()
                self.repaint()
                # self.scale(1.0, 1.0)
                
            self.dragLastPtF = ptf
            self.dragLastPt = pt

            
    def setCursorsVisible(self, vis):
        self.line1.setVisible(vis)
        self.line2.setVisible(vis)
        
        self.line1_cap_start.setVisible(vis)
        self.line1_cap_end.setVisible(vis)
        self.line2_cap_start.setVisible(vis)
        self.line2_cap_end.setVisible(vis)
        
        
    def mousePressEvent(self, mouseEvent): 
        # call super to ensure that mouse events are correctly propogated to items
        QtGui.QGraphicsView.mousePressEvent(self, mouseEvent)
        
        if self._ROIdragMode != ROIImgViewdDragMode.NONE:  
            return  # mouse press event already handled by items

        pt = (mouseEvent.x(), mouseEvent.y())
        qpt = QtCore.QPoint(pt[0], pt[1])
        ptf = self.mapToScene(qpt)
        img_ptf = self.pixMapItem.mapFromScene(ptf)
        img_pt = (int(round(img_ptf.x())), int(round(img_ptf.y())))
        pt = (ptf.x(), ptf.y())
        print("ROIImageGraphicsView.mousePressEvent x= %d y= %d img_pt=%s  button= %s" % (pt[0], pt[1], repr(img_pt), repr(mouseEvent.button())))
        self.buttonPressed = mouseEvent.button()
        if mouseEvent.button() == QtCore.Qt.LeftButton:
            print("    left button");
            
            # hit test on circusors
            self._ROIdragMode = ROIImgViewdDragMode.DRAWROI
        elif mouseEvent.button() ==  QtCore.Qt.RightButton:
            print("    right button");
            self._ROIdragMode = ROIImgViewdDragMode.MOVE
            self.dragLastPt = pt
            self.dragLastPtF = ptf
            self.isMouseDrag = True
            return

        DebugLog.log("ROIImageGraphicsView.mouseMoveEvent: drawtype = %s" % self._ROIdrawType.name)
            
        if self._ROIdrawType == ROIImgViewROIDrawType.BOX:
            self.isMouseDrag = True
            self.ROIBox_pt1 = pt
        elif self._ROIdrawType == ROIImgViewROIDrawType.SINGLE_PT:
            d = self.xhair_d
            x1 = pt[0] - d
            x2 = pt[0] + d
            y1 = pt[1] - d
            y2 = pt[1] + d
            self.xhairItem1.setLine(x1, y1, x2, y2)
            self.xhairItem2.setLine(x1, y2, x2, y1)
            self.singlePt = img_pt
            
            if self.mainOCTObj is not None:
                self.mainOCTObj.mscanSinglePtSet(self.singlePt)
                
        elif self._ROIdrawType == ROIImgViewROIDrawType.MULTI_PT:
            d = self.xhair_d
            
            x1 = pt[0] - d
            x2 = pt[0] + d
            y1 = pt[1] - d
            y2 = pt[1] + d
            
            txtX = pt[0] + 2*d
            txtY = pt[1]
            scene = self.scene()
            pen = self.ptPen
            
            if self.multPtMode == ROIImgViewMultPtMode.NEW:
                xhairItem1 = scene.addLine(x1, y1, x2, y2, pen)
                xhairItem2 = scene.addLine(x1, y2, x2, y1, pen)
                ptNum = len(self.ptsList)
                
                #
                #textItem.textCursor().charFormat().setBackground(brush)
                #textItem.textCursor().blockCharFormat().setForeground(brush)
                #textItem.textCursor().blockCharFormat().setBackground(brush)
        
                #textItem.setZValue(-2)
                #textItem.setDefaultTextColor(clr)

                txtItem = scene.addSimpleText(repr(ptNum + 1), self.ptFont)
                txtItem.setPos(txtX, txtY)
                #txtItem.textCursor().charFormat().setForeground(self.txtBrush)
                #txtItem.textCursor().blockCharFormat().setForeground(self.txtBrush)
                #txtItem.textCursor().charFormat().setBackground(self.txtBrush)
                #txtItem.textCursor().blockCharFormat().setBackground(self.txtBrush)
                txtItem.setBrush(self.txtBrush)
                # txtItem.setDefaultColor(self.txtClr)
                
                self.ptsList.append(img_pt)
                self.multiPt_grItems.append((xhairItem1, xhairItem2, txtItem))
            elif self.multPtMode == ROIImgViewMultPtMode.MOVE:
                ptNum = len(self.ptsList) - 1
                (xhairItem1, xhairItem2, txtItem) = self.multiPt_grItems[ptNum]
                self.ptsList[ptNum] = img_pt
                
                xhairItem1.setLine(x1, y1, x2, y2)
                xhairItem2.setLine(x1, y2, x2, y1)
                txtItem.setPos(txtX, txtY)
                
                self.ptsList[ptNum] = img_pt
                
        elif self._ROIdrawType == ROIImgViewROIDrawType.POLYGON:
            qptf = QtCore.QPointF(pt[0], pt[1])
            if qptf != self.ROI_poly.last():
                scene = self.scene()
                self.ROI_poly.append(qptf)
                scene.removeItem(self.ROI_poly_item)
                self.ROI_poly_item = scene.addPolygon(self.ROI_poly, pen=self.linePen)
                
        elif self._ROIdrawType == ROIImgViewROIDrawType.FREE:
            self.isMouseDrag = True
        
    def mouseReleaseEvent(self, mouseEvent):
        pt = (mouseEvent.x(), mouseEvent.y())
        qpt = QtCore.QPoint(pt[0], pt[1])
        ptf = self.mapToScene(qpt)
        pt = (ptf.x(), ptf.y())
        
        if self._ROIdrawType == ROIImgViewROIDrawType.BOX:
            self.ROIBox_pt2 = pt
            width = np.abs(self.ROIBox_pt2[0] - self.ROIBox_pt1[0])
            height = np.abs(self.ROIBox_pt2[1] - self.ROIBox_pt1[1])
            print("ROIImageGraphicsView.mouseReleaseEvent ROI ul= %s lr = %s width= %d height= %d" % (repr(self.ROIBox_pt1), repr(self.ROIBox_pt2), width, height))
            self.resetCachedContent()
            if self.mainOCTObj is not None:
                pt1 = self.pixMapItem.mapFromScene(QtCore.QPoint(self.ROIBox_pt1[0], self.ROIBox_pt1[1]))
                pt2 = self.pixMapItem.mapFromScene(QtCore.QPoint(self.ROIBox_pt2[0], self.ROIBox_pt2[1]))
                self.mainOCTObj.BMscanBoxRegionSet(pt1, pt2)
                
        elif (self.mainOCTObj is not None) and (self.cursorItemSelected is not None):
            self.mainOCTObj.ROIImgViewCursorsMoved(self)
        
        self._ROIdragMode = ROIImgViewdDragMode.NONE
        self.isMouseDrag = False
        self.cursorItemSelected = None
        
        print("ROIImageGraphicsView.mouseReleaseEvent x= %d y= %d button= %s ROIdragMode= %s" % (pt[0], pt[1], repr(mouseEvent.button()), repr(self._ROIdragMode)))

    def getROIBoxWidthHeight(self):
        width = 0
        height = 0
        if (self.ROIBox_pt1 is not None) and (self.ROIBox_pt2 is not None):
            width = np.abs(self.ROIBox_pt2[0] - self.ROIBox_pt1[0])
            height = np.abs(self.ROIBox_pt2[1] - self.ROIBox_pt1[1])
            
        return (width, height)
        
    def getImageWidthHeight(self):
        return (self.imgData.shape[0], self.imgData.shape[1])
        
    def wheelEvent(self, evt):
        # print("wheelEvent x= %d y= %d delta= %d" % (evt.x(), evt.y(), evt.delta()))
        zoomOutF = 1.25
        zoomInF = 1/zoomOutF
        tform = self.transform()
        tf = np.zeros((3, 3))
        tf[0][0] = tform.m11()
        tf[0][1] = tform.m12()
        tf[0][2] = tform.m13()
        tf[1][0] = tform.m21()
        tf[1][1] = tform.m22()
        tf[1][2] = tform.m23()
        tf[2][0] = tform.m31()
        tf[2][1] = tform.m32()
        tf[2][2] = tform.m33()
        
        if evt.delta() > 0:
            self.scale(zoomOutF, zoomOutF)
            # print("wheelEvent zoomOut trasnform= \n%s" % (tf))
        else:
            self.scale(zoomInF, zoomInF)
            # print("wheelEvent zoomIn trasnform= \n%s" % (tf))
            
        hscroll = self.horizontalScrollBar()
        vscroll = self.verticalScrollBar()
        # print("\nwheelEvent hscroll value= %d vscroll value= %d" % (hscroll.value(), vscroll.value()))
        
    def setImage(self, imgData, clrMap=None, resetTransform=True):
        
        if imgData is None:
            self.qImage = None
            self.qPixmap = None
            if self.pixMapItem != None:
                self.scene().removeItem(self.pixMapItem)

            self.pixMapItem = None
            self.imgData = None
        else:
            if not isinstance(imgData, np.ndarray):
                raise UnsupportedImageFormatException("Image data must be numpy ndarray")

            old_shp = (-1, -1)            
            if self.imgData is not None:
                old_shp = self.imgData.shape
                
            shp = imgData.shape
            w = shp[0]
            h = shp[1]
            
#            print('ROIImageGaphicsView.setImage() imgData shape= %s dtype= %s' % (repr(shp), repr(imgData.dtype)))
            if (len(shp) == 2) and (imgData.dtype == 'uint8'):
                fmt = QtGui.QImage.Format_Indexed8
                imgData = np.require(imgData, np.uint8, 'C')
            else:
                raise UnsupportedImageFormatException("Unsupported image data format shape= %s dtype= %s" % (repr(imgData.shape), repr(imgData.dtype)))
                
            # if no change in shape, then just update image
#            if (old_shp[0] == shp[0]) and (old_shp[1] == shp[1]):
#                self.imgData[:, :] = imgData[:, :]
#                qImg = QtGui.QImage(imgData.data, h, w, h, fmt)
#                self.qImage = qImg
#                self.qPixmap.
                
            self.imgData = copy.copy(imgData)
            
            qImg = QtGui.QImage(imgData.data, h, w, h, fmt)
            
            if clrMap is not None:
                qImg.setColorTable(clrMap)
                
            qPixMap = QtGui.QPixmap.fromImage(qImg)
            
            if self.pixMapItem is None:
                pixMapItem = self.scene().addPixmap(qPixMap)
                pixMapItem.setZValue(-10)
                # center item on center of view
                self.pixMapItem = pixMapItem
            else:
                self.pixMapItem.setPixmap(qPixMap)

            # ensure item is centered in the coordinate system            
            self.pixMapItem.setPos(-h/2, -w/2)
            
            self.qImage = qImg
            self.qPixmap = qPixMap
            
            if resetTransform:
                

                # reset transform to identity with scale factor
                tForm = QtGui.QTransform()
                w_g = self.width()
                h_g = self.height()
                s = min(w_g/w, h_g/h)
                tForm.scale(s, s)
                self.setTransform(tForm)
                
                # reset line positions
                qpt1 = QtCore.QPoint(0, h // 2)
                qpt2 = QtCore.QPoint(w, h // 2)
                pt1f = self.mapToScene(qpt1)
                pt2f = self.mapToScene(qpt2)
                
                l = self.line1.line()
                l.setP1(pt1f)
                l.setP2(pt2f)
                self.line1.setLine(l)
                r = self.circR
                
                self.line1_cap_start.setPos(pt1f)
                rct = self.line1_cap_end.rect()
                pt2f.setX(pt2f.x() - r)
                pt2f.setY(pt2f.y() - r)
                rct.moveTo(pt2f)
                self.line1_cap_end.setRect(rct)
                
                self.line2.setLine(l)
                self.line2_cap_start.setPos(pt1f)
                self.line2_cap_end.setRect(rct)
                
                # center the image in the view
                self.centerOn(0, 0)
            
            # self.line1_cap_end.setPos(pt2f)
            
        self.resetCachedContent()            
        self.update()
        self.repaint()
    
    def setROIDrawType(self, drawType):
        self._ROIdrawType = drawType
        self.xhairItem1.setVisible(False)
        self.xhairItem2.setVisible(False)
        self.rectItem.setVisible(False)
        self.freedraw_poly_item.setVisible(False)
        for grItems in self.multiPt_grItems:
            for itm in grItems:
                itm.setVisible(False)
            
        self.ROI_poly_item.setVisible(False)
            
        if drawType == ROIImgViewROIDrawType.NONE:
            pass
        elif drawType == ROIImgViewROIDrawType.BOX: 
            self.rectItem.setVisible(True)
        elif drawType == ROIImgViewROIDrawType.SINGLE_PT:
            self.xhairItem1.setVisible(True)
            self.xhairItem2.setVisible(True)
        elif drawType == ROIImgViewROIDrawType.MULTI_PT: 
            for grItems in self.multiPt_grItems:
                for itm in grItems:
                    itm.setVisible(True)
        elif drawType == ROIImgViewROIDrawType.FREE:
            self.freedraw_poly_item.setVisible(True)
        elif drawType == ROIImgViewROIDrawType.POLYGON:
            self.ROI_poly_item.setVisible(True)
        
        
    def clearPolygon(self):
        self.ROI_poly.clear()
        scene = self.scene()
        scene.removeItem(self.ROI_poly_item)
        
    def clearMultiPt(self):
        scene = self.scene()
        self.ptsList.clear()
        for grItems in self.multiPt_grItems:
            for itm in grItems:
                scene.removeItem(itm)

    def deleteLastMultiPt(self):
        scene = self.scene() 
        if len(self.multiPt_grItems) > 0:
            self.ptsList.pop()
            grItems = self.multiPt_grItems.pop()
            for itm in grItems:
                scene.removeItem(itm)

    def clearFreeDrawROI(self):
        self.freedraw_poly.clear()
        scene = self.scene()
        scene.removeItem(self.freedraw_poly_item)
        # self.freedraw_poly_item = scene.addPolygon(self.freedraw_poly, pen=self.linePen)
        
        
    #def dragMoveEvent(self, dragEvent):
    #    print("dragMoveEvent")
        
class CursorItem(QtGui.QGraphicsLineItem):
    def __init__(self, parent=None):
        QtGui.QGraphicsLineItem.__init__(self, parent)
        self.ROIgraphicsView = None
        self.capStartItem = None
        self.capEndItem = None
        self.dragLastPtF = None
        self.linkedLine = None
        self.cursorNum = 0
    
    def setCapItems(self, capStart, capEnd):
        self.capStartItem = capStart
        self.capEndItem = capEnd
        
    def mousePressEvent(self, mouseEvent): 
        qpt = mouseEvent.pos()
        ptf = self.mapToScene(qpt)
        pt = (ptf.x(), ptf.y())
        #print("CursorItem.mousePressEvent x= %d y= %d button= %s" % (pt[0], pt[1], repr(mouseEvent.button())))
        
        if self.ROIgraphicsView is not None:
            self.ROIgraphicsView._ROIdragMode = ROIImgViewdDragMode.MOVE_ITEM
            self.ROIgraphicsView.isMouseDrag = True
            self.ROIgraphicsView.cursorItemSelected = self
            self.ROIgraphicsView.dragLastPt = pt
            self.ROIgraphicsView.dragLastPtF = ptf
            self.dragLastPtF = ptf
            
        mouseEvent.ignore()
        
    def moveEvent(self, ptf, ptfLast):
        lastPtF = self.ROIgraphicsView.dragLastPtF
        dy = ptf.y() - lastPtF.y()
        dx = ptf.x() - lastPtF.x()
        #print("CursorItem.moveEvent ptf.y()= %f lastPtF.y()= %f dy= %f" % (ptf.y(), self.dragLastPtF.y(), dy))
        self.moveBy(dx, dy)
        self.capStartItem.moveBy(dx, dy)
        self.capEndItem.moveBy(dx, dy)
        
            
class CursorCapItem(QtGui.QGraphicsEllipseItem):
    def __init__(self, parent=None):
        QtGui.QGraphicsLineItem.__init__(self, parent)
        self.ROIgraphicsView = None
        self.cursorItem = None
        self.dragLastPtF = None
        self.isStart = True
        
    def setCursorItem(self, cursItem):
        self.cursorItem = cursItem
        
    def mousePressEvent(self, mouseEvent): 
        qpt = mouseEvent.pos()
        ptf = self.mapToScene(qpt)
        pt = (ptf.x(), ptf.y())
        #print("CursorCapItem.mousePressEvent x= %d y= %d button= %s" % (pt[0], pt[1], repr(mouseEvent.button())))
        
        if self.ROIgraphicsView is not None:
            self.ROIgraphicsView._ROIdragMode = ROIImgViewdDragMode.MOVE_ITEM
            self.ROIgraphicsView.isMouseDrag = True
            self.ROIgraphicsView.cursorItemSelected = self
            self.ROIgraphicsView.dragLastPt = pt
            self.ROIgraphicsView.dragLastPtF = ptf
            self.dragLastPtF = ptf
            
        mouseEvent.ignore()
            
    def moveEvent(self, ptf, ptfLast):
        #print("CursorCapItem.MoveEvent ptf= %s" % (repr(ptf)))
        lastPtF = self.ROIgraphicsView.dragLastPtF
        dy = ptf.y() - lastPtF.y()
        dx = ptf.x() - lastPtF.x()
        self.moveBy(dx, dy)
        l = self.cursorItem.line()
        if(self.isStart):
            p = l.p1()
            p.setY(p.y() + dy)
            p.setX(p.x() + dx)
            l.setP1(p)
        else:
            p = l.p2()
            p.setY(p.y() + dy)
            p.setX(p.x() + dx)
            l.setP2(p)
            
        self.cursorItem.setLine(l)
        
        linkedLine = self.cursorItem.linkedLine
        if linkedLine is not None:
            l = linkedLine.line()
            if self.isStart:
                linkedLine.capStartItem.moveBy(dx, dy)
                p = l.p1()
                p.setY(p.y() + dy)
                p.setX(p.x() + dx)
                l.setP1(p)
            else:
                linkedLine.capEndItem.moveBy(dx, dy)
                p = l.p2()
                p.setY(p.y() + dy)
                p.setX(p.x() + dx)
                l.setP2(p)
                
            linkedLine.setLine(l)
        
        
from enum import Enum
class ROIImgViewROIDrawType(Enum):
    NONE = 0     # does nothing
    BOX = 1      # make ROI box  
    SINGLE_PT = 2  # draw a crosshair at single point when mouse is pressed
    MULTI_PT = 3   # draw multiple points with crosshairs 
    POLYGON = 4   # draw series of lines
    FREE = 6      # draw arbitrary closed path 

class ROIImgViewdDragMode(Enum):
    NONE = 0       # does nothing
    MOVE = 1       # move the viewport
    DRAWROI = 2    # draw ROI mode
    MOVE_ITEM = 3  
    
class ROIImgViewMultPtMode(Enum):
    NEW = 1       # add new point
    MOVE = 2      # move existing point
    
        
class UnsupportedImageFormatException(Exception):
    pass