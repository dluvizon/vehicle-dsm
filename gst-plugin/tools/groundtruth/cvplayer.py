#!/usr/bin/python2.7

import sys, os
import shutil
import signal
import numpy
import glob
import time

# we have to do this to tells python to find cv2.so here
sys.path.insert(0, '/usr/local/lib/python2.7/site-packages')
import cv2
from cv2 import cv

from threading import Thread, Lock

import xml.dom.minidom as minidom
import xml.etree.cElementTree as cet
import xml.etree.ElementTree as et

LANE1_LIMIT_POS1 = 535
LANE2_LIMIT_POS1 = 1288
LANE1_LIMIT_POS2 = 660
LANE2_LIMIT_POS2 = 1340

def strtobool(s):
    return  s in ['TRUE', 'True', 'true', '1', 'T', 't', 'Y', 'y', 'YES', 'Yes', 'yes']

class PlateDetectGt:
    def __init__(self, frame, lane, speed, plate, radar, sema, moto, x, y, w, h):
        self.frame = frame
        self.lane = lane
        self.speed = speed
        self.plate = plate
        self.radar = radar
        self.sema = sema
        self.moto = moto
        self.x = x
        self.y = y
        self.w = w
        self.h = h

class XmlWriter:
    def __init__(self, fname):
        self.fname = fname
        self.iframe = 0
        self.vfPlates = {}
        if os.path.isfile(fname):
            print ('File "' + str(fname) + '" already exist! Creating a backup')
            bkpfile = fname + '.bkp'
            os.rename(fname, bkpfile)
            # parse the xml read and write it to a new file
            self.reparse(bkpfile)
            self.write()
        else:
            self.root = cet.Element('GroundTruthRoot')
            self.doc = cet.SubElement(self.root, 'gtruth')

    def reparse(self, fname):
        self.root = cet.Element('GroundTruthRoot')
        tree = et.parse(fname)
        root = tree.getroot()
        for idoc in root:
            if idoc.tag == 'gtruth':
                self.doc = cet.SubElement(self.root, 'gtruth')
                for vehicle in idoc:
                    if vehicle.tag == 'vehicle':
                        self.reparseVehicle(vehicle)
            elif idoc.tag == 'videoframes':
                self.iframe = int(idoc.attrib.get('total'))

    def reparseVehicle(self, vehicle):
        frame = int(vehicle.attrib.get('iframe')) + 1
        lane = int(vehicle.attrib.get('lane'))
        plate = strtobool(vehicle.attrib.get('plate'))
        radar = strtobool(vehicle.attrib.get('radar'))
        sema = strtobool(vehicle.attrib.get('sema'))
        moto = strtobool(vehicle.attrib.get('moto'))
        x = 0
        y = 0
        w = 0
        h = 0
        speed = 0.0
        for child in vehicle:
            if child.tag == 'region':
                x = int(child.attrib.get('x'))
                y = int(child.attrib.get('y'))
                w = int(child.attrib.get('w'))
                h = int(child.attrib.get('h'))
            elif (child.tag == 'radar') and radar:
                speed = float(child.attrib.get('speed'))
        self.vfPlates[frame] = PlateDetectGt(frame, lane, speed, plate, radar, sema, moto, x, y, w, h)
        self.newVehicle(frame, plate, sema, moto, x, y, w, h, radar, speed, lane)
        return

    def write(self):
        xmlstr = et.tostring(self.root, encoding='utf8', method='xml')
        reparsed = minidom.parseString(xmlstr)
        file = open(self.fname, "w")
        file.write(reparsed.toprettyxml(indent="  "))
        file.close()

    def newVehicle(self, iframe, plate, sema, moto, x, y, w, h, radar, speed, lane):
        # in this software, the iframe starts from 1, but it is more usually
        # (and more logical) to start from 0. So subtract one before write
        iframe = int(iframe) - 1
        lane = int(lane)
        vehicle = cet.SubElement(self.doc, 'vehicle')
        vehicle.set('iframe', str(iframe))
        vehicle.set('plate', str(plate))
        vehicle.set('sema', str(sema))
        vehicle.set('moto', str(moto))
        vehicle.set('radar', str(radar))
        vehicle.set('lane', str(lane))

        reg = cet.SubElement(vehicle, 'region')
        reg.set('x', str(x))
        reg.set('y', str(y))
        reg.set('w', str(w))
        reg.set('h', str(h))

        if radar:
            radar = cet.SubElement(vehicle, 'radar')
            radar.set('frame_start', str(iframe))
            radar.set('frame_end', str(iframe + 40))
            radar.set('speed', str(speed))

    def saveTotalFrames(self, iframes):
        frames = cet.SubElement(self.root, 'videoframes')
        frames.set('total', str(iframes))

class Vehicle:
    def __init__(self, frame, lane, speed):
        self.frame = int(frame)
        self.lane = int(lane)
        self.speed = float(speed)
        self.valid = True
        self.used = False
        if (self.lane < 1) or (self.lane > 3):
            self.valid = False
        if (self.speed < 0.0) or (self.speed > 199.0):
            self.valid = False

class XmlGTruthReader:
    def __init__(self, fname):
        if os.path.isfile(fname) is False:
            raise Exception('Error: file ' + fname + "' does not exist!")
        self.vehicles = []
        self.parseXML(fname)
        print ('Parse of ground-truth.xml done!')

    def parseXML(self, fname):
        self.tree = cet.parse(fname)
        self.root = self.tree.getroot()
        for child in self.root:
            if child.tag == 'GroundTruth':
                self.parseGroundTruth(child)

    def parseGroundTruth(self, gtruth):
        for child in gtruth:
            if child.tag == 'vehicle':
                self.parseVehicle(child)

    def parseVehicle(self, vehicle):
        # subtract 10 frames to adjust from the loops to the bottom side of the
        # frame
        self.vehicles.append(Vehicle(int(vehicle.attrib.get('frame')) - 10,
            vehicle.attrib.get('lane'), vehicle.attrib.get('speed')))

    def showXML(self):
        for v in self.vehicles:
            if v.valid:
                print ('vehicle: ' + str(v.frame) + '\t' + str(v.lane) + '\t' + str(v.speed))
        for v in self.vehicles:
            if not v.valid:
                print ('INV.: ' + str(v.frame) + '\t' + str(v.lane) + '\t' + str(v.speed))

def enum(**enums):
    return type('Enum', (), enums)
# place selection state
pss = enum(IDLE='idle', CLICKED='clicked', READY='ready')

class PlaceSelection:
    def __init__(self):
        self.startX = 0
        self.startY = 0
        self.endX = 0
        self.endY = 0
        self.state = pss.IDLE
        self.refreshDraw = False

    def setStart(self, x, y):
        self.startX = x
        self.startY = y
        self.refreshDraw = False

    def setEnd(self, x, y):
        self.endX = x
        self.endY = y
        self.refreshDraw = True

class MouseEvent:
    def __init__(self):
        # moving event
        self.currX = 0
        self.currY = 0
        self.lastX = 0
        self.lastY = 0
        self.moved = False
        # clicked event
        self.clickX = 0
        self.clickY = 0
        self.click = False
        # released event
        self.releaseX = 0
        self.releaseY = 0
        self.release = False
        # cancellation event
        self.cancel = False
        # mutex
        self.mutex = Lock()

    def lock(self):
        self.mutex.acquire()

    def unlock(self):
        self.mutex.release()

    def testLock(self):
        return self.mutex.acquire(False)

    def newEventMove(self, x, y):
        self.currX = x
        self.currY = x
        dx = abs(self.lastX - x)
        dy = abs(self.lastY - y)
        if (dx < 2) and (dy < 2):
            # this motion is not relevant, inside a grid of 2x2
            return
        self.mutex.acquire()
        self.lastX = x
        self.lastY = y
        self.moved = True
        self.mutex.release()

    def newEventClick(self, x, y):
        self.mutex.acquire()
        self.clickX = x
        self.clickY = y
        self.click = True
        # for safety, clear the release flag
        self.release = False
        self.mutex.release()

    def newEventRelease(self, x, y):
        self.mutex.acquire()
        if self.click:
            # if is not clicked, it will be ignored
            self.releaseX = x
            self.releaseY = y
            self.release = True
        self.mutex.release()

    def newEventCalcel(self):
        self.mutex.acquire()
        self.cancel = True
        self.moved = False
        self.click = False
        self.release = False
        self.mutex.release()

class PlayVideo:
    def __init__(self, videoFile, gTruthXmlFile, outputXml, platesDir):
        self.xml = XmlWriter(outputXml)
        self.gtruth = XmlGTruthReader(gTruthXmlFile)
        self.capture = cv.CaptureFromFile(videoFile)
        cv.NamedWindow("Video", cv.CV_WINDOW_AUTOSIZE)
        self.w = cv.GetCaptureProperty(self.capture, cv.CV_CAP_PROP_FRAME_WIDTH)
        self.h = cv.GetCaptureProperty(self.capture, cv.CV_CAP_PROP_FRAME_HEIGHT)
        self.nFrames = cv.GetCaptureProperty(self.capture,
                cv.CV_CAP_PROP_FRAME_HEIGHT)
        print ('Total number of Frames in video: ' + str(self.nFrames))
        self.font = cv.InitFont(cv.CV_FONT_HERSHEY_PLAIN, 1, 1, 0, 1, 1)
        print ("resolution: " + str(int(self.w)) + "x" + str(int(self.h)))
        self.fx = 0.85
        self.fy = 0.85
        self.iframe = self.xml.iframe
        self.plate = True
        self.radar = False
        self.semaphore = False
        self.moto = False
        self.currv = None
        self.mouse = MouseEvent()
        self.select = PlaceSelection()
        ' create platesDir if needed '
        self.plateCarRadar = platesDir + "/car-radar"
        self.plateCarNoRad = platesDir + "/car-noradar"
        self.plateMotoRadar = platesDir + "/moto-radar"
        self.plateMotoNoRad = platesDir + "/moto-noradar"
        if not os.path.isdir(platesDir):
            os.mkdir(platesDir)
        if not os.path.isdir(self.plateCarRadar):
            os.mkdir(self.plateCarRadar)
        if not os.path.isdir(self.plateCarNoRad):
            os.mkdir(self.plateCarNoRad)
        if not os.path.isdir(self.plateMotoRadar):
            os.mkdir(self.plateMotoRadar)
        if not os.path.isdir(self.plateMotoNoRad):
            os.mkdir(self.plateMotoNoRad)

    def showSaveMenu(self):
        hp = 22
        pos1 = (int(self.select.startX), int(self.select.startY))
        pos2 = (int(self.select.endX), int(self.select.endY))
        pxy = (int(self.select.startX / self.fx), int(self.select.startY / self.fy))
        pax = (int(self.select.endX / self.fx), int(self.select.endY / self.fy))
        pwh = (pax[0] - pxy[0], pax[1] - pxy[1])
        numpyImg = cv2.resize(numpy.asarray(self.src[:,:]), (0, 0), fx=self.fx, fy=self.fy)
        cimg = cv.fromarray(numpyImg)
        # put info
        if self.radar:
            cv.Rectangle(cimg, (4, 4), (312, 270), cv.Scalar(220, 200, 196, 0), thickness=-1)
            cv.Rectangle(cimg, (4, 4), (312, 270), cv.Scalar(100, 80, 70, 0), thickness=2)
        else:
            cv.Rectangle(cimg, (4, 4), (312, 210), cv.Scalar(220, 200, 196, 0), thickness=-1)
            cv.Rectangle(cimg, (4, 4), (312, 210), cv.Scalar(100, 80, 70, 0), thickness=2)
        cv.PutText(cimg, '>> PRESS ENTER / BACKSPACE <<', (12 , 24),
                self.font, cv.Scalar(64, 64, 244, 0))
        cv.PutText(cimg, '(x, y) = ' + str(pxy), (12 , 2 * hp),
                self.font, cv.Scalar(16, 16, 200, 0))
        cv.PutText(cimg, '(w, h) = ' + str(pwh), (12 , 3 * hp),
                self.font, cv.Scalar(16, 16, 200, 0))
        cv.PutText(cimg, 'frame = ' + str(int(self.iframe)), (12 , 4 * hp),
                self.font, cv.Scalar(16, 16, 200, 0))
        cv.PutText(cimg, 'lane = ' + str(int(self.currlane)), (12 , 5 * hp),
                self.font, cv.Scalar(16, 16, 200, 0))
        cv.PutText(cimg, '<P>late =       ' + str(self.plate), (12 , 6 * hp),
                self.font, cv.Scalar(16, 16, 200, 0))
        cv.PutText(cimg, '<S>emaphore = ' + str(self.semaphore), (12 , 7 * hp),
                self.font, cv.Scalar(16, 16, 200, 0))
        cv.PutText(cimg, '<M>oto =      ' + str(self.moto), (12 , 8 * hp),
                self.font, cv.Scalar(16, 16, 200, 0))
        cv.PutText(cimg, '<R>adar =       ' + str(self.radar), (12 , 9 * hp),
                self.font, cv.Scalar(16, 16, 200, 0))
        if self.radar:
            cv.PutText(cimg, '    frame = ' + str(self.currv.frame), (12 , 10 * hp),
                    self.font, cv.Scalar(96, 96, 96, 0))
            cv.PutText(cimg, '    lane = ' + str(self.currv.lane), (12 , 11 * hp),
                    self.font, cv.Scalar(96, 96, 96, 0))
            cv.PutText(cimg, '    speed = ' + str(self.currv.speed), (12 , 12 * hp),
                    self.font, cv.Scalar(64, 64, 244, 0))

        # draw plate
        cv.Rectangle(cimg, pos1, pos2, cv.Scalar(36, 186, 255, 0), thickness=4)
        cv.Rectangle(cimg, pos1, pos2, cv.Scalar(8, 44, 72, 0), thickness=2)
        cv.ShowImage("Camera", cimg)

    def findNearestGTruth(self, maxDist):
        self.currv = None
        dist = sys.maxint
        for v in self.gtruth.vehicles:
            if not v.valid or v.used:
                continue
            if v.lane != self.currlane:
                continue
            d = abs(v.frame - self.iframe)
            if int(d) < int(dist):
                dist = d
                if dist < maxDist:
                    self.currv = v

    def chooseUp(self):
        found = False
        for v in self.gtruth.vehicles:
            if self.currlane != v.lane:
                continue
            if found:
                self.currv = v
                return
            if v == self.currv:
                found = True

    def chooseDown(self):
        last = None
        for v in self.gtruth.vehicles:
            if self.currlane != v.lane:
                continue
            if v == self.currv:
                if last != None:
                    self.currv = last
                return
            last = v

    def getLaneByPosition(self, x1, x2):
        mx = (x1 + x2) / 2
        if mx < LANE1_LIMIT_POS2:
            return 1
        elif mx < LANE2_LIMIT_POS2:
            return 2
        else:
            return 3

    def saveMenu(self):
        x1 = self.select.startX / self.fx
        y1 = self.select.startY / self.fy
        x2 = self.select.endX / self.fx
        y2 = self.select.endY / self.fy

        self.currlane = self.getLaneByPosition(x1, x2)
        self.findNearestGTruth(45)
        self.plate = True
        self.semaphore = False
        self.moto = False
        if self.currv != None:
            self.radar = True
        else:
            self.radar = False
        self.showSaveMenu()
        while True:
            c = cv.WaitKey(30)
            if c == -1:
                continue
            if (c == 112) or (c == 65616): # p
                self.plate = not self.plate
            elif (c == 115) or (c == 65619): # s
                self.semaphore = not self.semaphore
            elif (c == 109) or (c == 65613): # m
                self.moto = not self.moto
            elif (c == 114) or (c == 65618): # r
                self.radar = not self.radar
                if self.radar and self.currv == None:
                    self.findNearestGTruth(sys.maxint)
            elif ((c == 1113938) or (c == 65362)) and self.radar: # Up
                self.chooseUp()
            elif ((c == 1113940) or (c == 65364)) and self.radar: # Down
                self.chooseDown()
            self.showSaveMenu()
            if (c & 0xFF) == 10:
                print ('Confirm')
                x = int(x1)
                y = int(y1)
                w = int(x2 - x1)
                h = int(y2 - y1)
                if self.radar and self.currv != None:
                    self.xml.newVehicle(self.iframe, self.plate,
                            self.semaphore, self.moto, x, y, w, h,
                            self.radar, self.currv.speed, self.currlane)
                    self.currv.used = True
                else:
                    self.xml.newVehicle(self.iframe, self.plate,
                            self.semaphore, self.moto, x, y, w, h,
                            self.radar, 0.0, self.currlane)
                self.xml.write()
                return True
            elif (c & 0xFF) == 8:
                print ('Calcel')
                return False

    def on_dummy(self, e, x, y, flag, param):
        return

    def on_mouse(self, e, x, y, flag, param):
        if (e == cv.CV_EVENT_MOUSEMOVE):
            self.mouse.newEventMove(x, y)
        elif (e == cv.CV_EVENT_LBUTTONDOWN):
            self.mouse.newEventClick(x, y)
        elif (e == cv.CV_EVENT_LBUTTONUP):
            self.mouse.newEventRelease(x, y)
        elif (e == cv.CV_EVENT_RBUTTONUP) or (e == cv.CV_EVENT_RBUTTONDOWN):
            self.mouse.newEventCalcel()

    def readMouseEvents(self):
        self.mouse.lock()
        if self.mouse.cancel:
            self.select.state = pss.IDLE
            self.mouse.cancel = False
            self.drawPlate(False)
        if self.select.state == pss.IDLE:
            if self.mouse.click:
                self.select.setStart(self.mouse.lastX, self.mouse.lastY)
                self.select.setEnd(self.mouse.lastX, self.mouse.lastY)
                self.select.state = pss.CLICKED
        elif self.select.state == pss.CLICKED:
            if self.mouse.moved:
                self.select.setEnd(self.mouse.lastX, self.mouse.lastY)
                self.mouse.moved = False
            elif self.mouse.release:
                self.select.setEnd(self.mouse.lastX, self.mouse.lastY)
                self.mouse.click = False
                self.mouse.release = False
                self.select.state = pss.READY
        elif self.select.state == pss.READY:
            'do nothing here'
        else:
            self.select.state = pss.IDLE
        self.mouse.moved = False
        self.mouse.unlock()

    def handleMouse(self):
        self.readMouseEvents()
        ret = False
        if self.select.state == pss.IDLE:
            return False
        # redraw image if needed
        if self.select.refreshDraw:
            self.select.refreshDraw = False
            self.drawPlate()
        if self.select.state == pss.READY:
            # disable mouse while handling save_menu
            cv.SetMouseCallback("Camera", self.on_dummy, param=0)
            # enter in save_menu
            ret = self.saveMenu()
            self.select.state = pss.IDLE
            cv.SetMouseCallback("Camera", self.on_mouse, param=0)
            self.drawPlate(False)
        return ret

    def drawPlate(self, plate=True):
        pos1 = (self.select.startX, self.select.startY)
        pos2 = (self.select.endX, self.select.endY)
        if (pos1[0] > pos2[0]) or (pos1[1] > pos2[1]):
            cor1 = cv.Scalar(0, 0, 255, 0)
            cor2 = cv.Scalar(255, 255, 255, 0)
            thick = 8
        else:
            cor1 = cv.Scalar(36, 186, 255, 0)
            cor2 = cv.Scalar(8, 44, 72, 0)
            thick = 4
        numpyImg = cv2.resize(numpy.asarray(self.src[:,:]), (0, 0),
                fx=self.fx, fy=self.fy)
        cimg = cv.fromarray(numpyImg)
        # draw plate if setted
        if plate:
            cv.Rectangle(cimg, pos1, pos2, cor1, thickness=thick)
            cv.Rectangle(cimg, pos1, pos2, cor2, thickness=2)
        cv.ShowImage("Camera", cimg)

    def showImg(self):
        small = cv2.resize(numpy.asarray(self.src[:,:]), (0, 0), fx=self.fx,
                fy=self.fy)
        cv.ShowImage("Camera", cv.fromarray(small))

    def draw(self):
        p1 = (0, 840)
        p2 = (int(self.w), 840)
        p3 = (0, 1040)
        p4 = (int(self.w), 1040)
        shift = int((self.w - 16 * len(self.text)) / 2)
        cv.PutText(self.src, self.text, (shift , 36), self.font,
                cv.Scalar(32, 0, 220, 0))
        cv.PutText(self.src, str(self.iframe), (880 , 52), self.font,
                cv.Scalar(32, 0, 220, 0))
        cv.Line(self.src, p1, p2, cv.Scalar(0, 64, 255, 0), thickness=7)
        cv.Line(self.src, p1, p2, cv.Scalar(160, 0, 0, 0), thickness=2)
        cv.Line(self.src, p3, p4, cv.Scalar(0, 64, 255, 0), thickness=7)
        cv.Line(self.src, p3, p4, cv.Scalar(160, 0, 0, 0), thickness=2)
        self.showImg()

    def printf(self, text):
        self.text = text
        self.draw()
    
    def savePlateDetectLog(self, iframe, pdg, path):
        xmlfile = str(path) + "%05d" % iframe + ".xml"
        root = cet.Element('PlateInfo')
        root.set('x', str(pdg.x))
        root.set('y', str(pdg.y))
        root.set('w', str(pdg.w))
        root.set('h', str(pdg.h))
        xmlstr = et.tostring(root, encoding='utf8', method='xml')
        reparsed = minidom.parseString(xmlstr)
        file = open(xmlfile, "w")
        file.write(reparsed.toprettyxml(indent="  "))
        file.close()

    def checkSaveImage(self, iframe):
        ' verifica se existe no dicionario '
        imgpath = None
        if iframe in self.xml.vfPlates:
            pdg = self.xml.vfPlates[iframe]
            if pdg.plate:
                if pdg.moto:
                    if pdg.radar:
                        imgpath = self.plateMotoRadar + "/"
                    else:
                        imgpath = self.plateMotoNoRad + "/"
                else:
                    if pdg.radar:
                        imgpath = self.plateCarRadar + "/"
                    else:
                        imgpath = self.plateCarNoRad + "/"
                imgfile = str(imgpath) + "%05d" % iframe + ".png"
                print ("Save png image in: " + imgfile)
                cv.SaveImage(imgfile, self.src)
                self.savePlateDetectLog(iframe, pdg, imgpath)

    def run(self):
        self.text = ""
        play = True
        delay = 200
        pos = 0
        newPos = 0
        skipping = False
        print ("")
        print ("start running!")
        print ("")
        while True:
            if play or (newPos != pos):
                while int(self.iframe) >= int(pos):
                    if int(self.iframe) > int(pos):
                        print ('Skip frame --' + str(int(pos)) + ' [%d %%]' %
                            (100.0 * pos / self.iframe))
                        skipping = True
                    if cv.GrabFrame(self.capture):
                        pos = cv.GetCaptureProperty(self.capture,
                                cv.CV_CAP_PROP_POS_FRAMES)
                        newPos = pos
                        self.src = cv.RetrieveFrame(self.capture)
                        self.checkSaveImage(pos)
                        self.draw();
                    else:
                        print ("Grab frame failed! It must be the end of the video.")
                        self.finalize()
                        exit(0)
                self.iframe = int(pos)
            if play:
                # if playing, wait the necessary time to maintain the fps
                if skipping:
                    play = False
                    skipping = False
                c = cv.WaitKey(delay)
            else:
                # if paused, do is as fast as possible (without screw)
                c = cv.WaitKey(10)
            if not play:
                if self.handleMouse():
                    # if saved a new plate, step to the next frame
                    newPos = pos + 1
            if c == -1:
                continue
            c &= 0xFF
            if c == 32:  # backspace
                play = not play
                if play:
                    self.printf("    PLAY")
                    cv.SetMouseCallback("Camera", self.on_dummy, param=pos)
                else:
                    self.printf("              -> PAUSE")
                    cv.SetMouseCallback("Camera", self.on_mouse, param=0)
            elif c == 82:  # arrow up
                delay = 25
            elif c == 84:  # arrow down
                delay = 200
            elif c == 27:
                #print 'Are you sure about closing this program?'
                #opt = input('If the answer is "yes", type the result of 2x2: ')
                c = cv.WaitKey(0) & 0xFF
                if c == 10:
                    print ('Good bye, cruel world!')
                    self.finalize()
                    break
                else:
                    print ('To exit, press ESC followed by ENTER')

    def finalize(self):
        self.xml.saveTotalFrames(int(self.iframe))
        self.xml.write()

def _signal_handler(signal, frame):
    print ('exiting by _signal_handler')
    sys.exit(0)

if __name__ == '__main__':
    signal.signal(signal.SIGABRT, _signal_handler)
    signal.signal(signal.SIGINT, _signal_handler)
    signal.signal(signal.SIGTERM, _signal_handler)
    if len(sys.argv) is not 2:
        print ('')
        print ('Usage: ' + str(sys.argv[0]) + ' << path-to-video-dir >>')
        exit(0)
    try:
        videoin = sys.argv[1] + '/video.h264'
        gtruth = sys.argv[1] + '/olddata/ground-truth.xml'
        output = sys.argv[1] + '/vehicles.xml'
        platesdir = sys.argv[1] + '/plates'
        player = PlayVideo(videoin, gtruth, output, platesdir)
        player.run()
    except (NameError,), e:
        print str(e)
        exit(1)

