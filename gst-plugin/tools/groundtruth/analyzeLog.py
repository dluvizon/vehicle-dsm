#!/usr/bin/python2.7

'''
@file: analyzeLog.py
@author: Diogo Luvizon <diogo@luvizon.com>
'''

import os, sys
import shutil
import signal
import glob
import numpy

import matplotlib.pyplot as plt

from collections import Counter

from shapely.geometry import Polygon
from shapely.geometry import Point

import xml.dom.minidom as minidom
import xml.etree.cElementTree as cet
import xml.etree.ElementTree as et

DETECT_MOTORCYCLE = True
SIMUL_PERFECT_DETECTOR = False

''' Utilities functions:
======================== '''
def _ch(text):
    if text[-1] is '/':
        return text[:-1]
    return text

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'

def printColor(color, vmsg):
    print (color + vmsg + bcolors.ENDC)

'''
@class: lineLog
================ '''
class lineLog():
    '''
    @note: This class store the data from one line of log.
    '''
    def __init__(self, par):
        self.frame = int(par[0])
        self.speed = float(par[1])
        self.err = float(par[2])
        self.x = int(par[3])
        self.y = int(par[4])
        self.dx = float(par[5])
        self.dy = float(par[6])
        self.ipm_x = float(par[7])
        self.ipm_y = float(par[8])
        self.nfeat = int(par[9])

'''
@class: laneInfo
================ '''
class laneInfo():
    '''
    @note: This class store the information about lane's boundary.
    '''
    def __init__(self, filename):
        self.tree = cet.parse(filename)
        self.root = self.tree.getroot()
        self.poly = {}
        for child in self.root:
            if child.tag == 'lanes':
                for l in child:
                    if l.tag == 'measure':
                        p = ((int(l.find('pTopLeft').attrib.get('x')), int(l.find('pTopLeft').attrib.get('y'))),
                             (int(l.find('pTopRight').attrib.get('x')), int(l.find('pTopRight').attrib.get('y'))),
                             (int(l.find('pBotRight').attrib.get('x')), int(l.find('pBotRight').attrib.get('y'))),
                             (int(l.find('pBotLeft').attrib.get('x')), int(l.find('pBotLeft').attrib.get('y'))))
                        self.poly[int(l.attrib.get('id'))] = Polygon(p)
            if child.tag == 'nLanes':
                self.nLanes = child.attrib.get('total')
    
    def contains(self, id, x, y):
        return self.poly[id].contains(Point(x, y))

'''
@class: Vehicle
================ '''
class Vehicle():
    '''
    @note: This class store the data from one file *.txt.
    '''
    def __init__(self, filename):
        with open(filename) as f:
            content = f.readlines()
        self.filename = filename
        self.id = int(content[3].split(" : ")[-1])
        self.lane = int(content[4].split(" : ")[-1])
        self.plate = int(content[5].split(" : ")[-1])
        self.radar = int(content[6].split(" : ")[-1])
        self.sema = int(content[7].split(" : ")[-1])
        self.moto = int(content[8].split(" : ")[-1])
        self.radar_speed = float(content[9].split(" : ")[-1])
        self.lp_width = int(content[10].split(" : ")[-1])
        self.lp_height = int(content[11].split(" : ")[-1])
        self.lines = []
        for l in content[16:]:
            self.lines.append(lineLog(l.split()))
        self.haveFinalSpeed = False
        self.err = 0.0
    
    def computeResult(self, lInfo):
        speed = []
        for l in self.lines:
            if lInfo.contains(self.lane, l.x, l.y):
                speed.append(l.speed)
        if len(speed) > 0:
            self.haveFinalSpeed = True
            speed.sort()
            if len(speed) > 4:
                speed = speed[2:-2]
            self.finalSpeed = numpy.mean(speed)
        else:
            printColor(bcolors.FAIL, "Don't have final speed: %s" % self.filename)

def isMatching(a, b):
    if a.lines[0].frame > b.lines[-1].frame:
        return False
    if b.lines[0].frame > a.lines[-1].frame:
        return False
    minc = max(a.lines[0].frame, b.lines[0].frame)
    maxc = min(a.lines[-1].frame, b.lines[-1].frame)
    for f in range(minc, maxc + 1):
        aindex = f - a.lines[0].frame
        bindex = f - b.lines[0].frame
        if (aindex >= len(a.lines)) or (bindex >= len(b.lines)):
            continue;
        ax1 = a.lines[aindex].x - (a.lp_width / 2)
        ax2 = a.lines[aindex].x + (a.lp_width / 2)
        ay1 = a.lines[aindex].y - (a.lp_height / 2)
        ay2 = a.lines[aindex].y + (a.lp_height / 2)
        bx1 = b.lines[bindex].x - (b.lp_width / 2)
        bx2 = b.lines[bindex].x + (b.lp_width / 2)
        by1 = b.lines[bindex].y - (b.lp_height / 2)
        by2 = b.lines[bindex].y + (b.lp_height / 2)
        pa = ((ax1, ay1), (ax1, ay2), (ax2, ay2), (ax2, ay1))
        pb = ((bx1, by1), (bx1, by2), (bx2, by2), (bx2, by1))
        polyA = Polygon(pa)
        polyB = Polygon(pb)
        polyIntersec = polyA.intersection(polyB)
        polyUnion = polyA.union(polyB)
        if (polyIntersec.area / polyUnion.area) > 0.5:
            return True
    return False

def searchMatching(ref, searchlist):
    for s in searchlist:
        if isMatching(ref, s):
            return s
    return None

def saveLogFile(vehicles, filename):
    file = open(filename, 'w')
    file.write(str(len(vehicles)) + '\n\n')
    for v in vehicles:
        if v.haveFinalSpeed:
            file.write(v.filename + ' %05d' % v.lines[0].frame + ' %.1f\n' % v.finalSpeed)
        else:
            file.write(v.filename + ' %05d' % v.lines[0].frame + ' N.A.\n')
        if hasattr(v, 'gtlink'):
            file.write('#\t' + v.gtlink + '\n')
    file.close()

def computeDetectionPASCAL(G, D, logPath, writeLog=False, auxText=""):
    P = []
    R = []
    for d in D:
        m = searchMatching(d, G)
        if m is not None:
            P.append(d)
    for g in G:
        m = searchMatching(g, D)
        if m is not None:
            R.append(m)
    printColor(bcolors.HEADER,  '\n  -- Metrics (PASCAL) ' + str(auxText) + ' --')
    printColor(bcolors.OKBLUE,  '  Length G:               %d' % len(G))
    printColor(bcolors.OKBLUE,  '  Length D:               %d' % len(D))
    printColor(bcolors.OKBLUE,  '  Length P:               %d' % len(P))
    printColor(bcolors.OKBLUE,  '  Length R:               %d' % len(R))
    if len(D) > 0:
        v = float(len(P)) / float(len(D))
    else:
        v = float('NaN')
    printColor(bcolors.OKGREEN, '  Precision:              %.3f' % (v))
    if len(G) > 0:
        v = float(len(R)) / float(len(G))
    else:
        v = float('NaN')
    printColor(bcolors.OKGREEN, '  Recall:                 %.3f\n' % (v))
    if writeLog:
        saveLogFile(G, logPath + '/PASCAL_G.txt')
        saveLogFile(D, logPath + '/PASCAL_D.txt')
        saveLogFile(P, logPath + '/PASCAL_P.txt')
        saveLogFile(R, logPath + '/PASCAL_R.txt')
    return [P, R]

def classifySpeedError_USA(speed, logPath):
    err_high = []
    err_low = []
    err_ok = []
    for r in speed:
        if r.err > 2.0:
            err_high.append(r)
        elif r.err < -3.0:
            err_low.append(r)
        else:
            err_ok.append(r)
    saveLogFile(err_high, logPath + '/speed_high.txt')
    saveLogFile(err_low, logPath + '/speed_low.txt')
    saveLogFile(err_ok, logPath + '/speed_ok.txt')
    return [err_high, err_low, err_ok]

def computeSpeed(radar, detection, logPath):
    speed = []
    no_speed = []
    no_match = []
    for r in radar:
        match = searchMatching(r, detection)
        if match is not None:
            match.gtlink = r.filename
            if not match.radar:
                match.radar_speed = r.radar_speed
            #elif match.radar_speed != r.radar_speed:
                #printColor(bcolors.FAIL, '%s\n\t' % match.filename + 'speed not matching! %f ' % match.radar_speed + 'vs %f' % r.radar_speed)
                #match.radar_speed = r.radar_speed
            if match.haveFinalSpeed:
                match.err = match.finalSpeed - match.radar_speed
                speed.append(match)
            else:
                no_speed.append(match)
        else:
            no_match.append(r)
    saveLogFile(radar, logPath + '/speed_groundtruth.txt')
    saveLogFile(speed, logPath + '/speed.txt')
    saveLogFile(no_speed, logPath + '/no_speed.txt')
    saveLogFile(no_match, logPath + '/no_match.txt')
    return [speed, no_speed, no_match]

def parseLogFiles(path, lInfo):
    vehicles = []
    flist = glob.glob(path + "/*.txt")
    flist.sort()
    for f in flist:
        v = Vehicle(f)
        v.computeResult(lInfo)
        vehicles.append(v)
    return vehicles

def saveSpeedCurve(speed, filename):
    file = open(filename, 'w')
    for s in speed:
        estimated = s.radar_speed + s.err
        file.write('%f' % s.radar_speed + '\t%f' % estimated + '\n')
    file.close()

def saveSpeedHostigram(speed, filename):
    err = []
    for s in speed:
        if hasattr(s, 'err'):
            err.append(s.err)
        else:
            print 'God!!!!'
    plt.xlabel('Error [km/h]')
    plt.ylabel('Number of Vehicles')
    limitMin = -7.0
    limitMax = 8.0
    speedMin = -3.0
    speedMax = 2.0
    yRange = 1300
    counts, bins, patches = plt.hist(err, bins=numpy.arange(limitMin, limitMax, 0.5), facecolor='lightgray', edgecolor='black')
    plt.xticks(numpy.arange(int(float(limitMin)), int(float(limitMax)), 1.0))
    plt.ylim(0.0, yRange)
    for patch, rightside, leftside in zip(patches, bins[1:], bins[:-1]):
        if leftside < speedMin:
            patch.set_facecolor([0.2,0.2,0.2])
        elif rightside > speedMax:
            patch.set_facecolor([0.2,0.2,0.2])
        #elif (rightside > 0.0) and (leftside < 0.0):
        #    patch.set_facecolor('blue')
    tmpMin = (speedMin - limitMin + 1) / (limitMax - limitMin + 1)
    tmpMax = (speedMax - limitMin + 1) / (limitMax - limitMin + 1)
    plt.axhline(y=(1.02 * max(counts)), xmin=tmpMin, xmax=tmpMax, color='b', linestyle='dashed')
    plt.axvline(x=speedMin, ymin=(0.97 * max(counts) / yRange), ymax=(1.07 * max(counts) / yRange), color='b')
    plt.axvline(x=speedMax, ymin=(0.97 * max(counts) / yRange), ymax=(1.07 * max(counts) / yRange), color='b')
    sTotal = len(err)
    a = numpy.array(err)
    b = numpy.where(numpy.logical_and(a>=speedMin, a<=speedMax))
    sOk = len(b[0])
    plt.text(x=-1.4, y=(1.04 * max(counts)), s='%.1f %%' % (100.0 * sOk / sTotal))
    #plt.grid()
    plt.savefig(filename, format='eps')
    #plt.show()

'''
@function: main
@brief: Main program function
'''
def main(pathList):
    tplates_moto = []
    tplates_car = []
    tdetection = []
    tdetection_moto = []
    tdetection_car = []
    tpre = []
    trec = []
    tpre_moto = []
    trec_moto = []
    tpre_car = []
    trec_car = []
    tradar = []
    tehigh = []
    telow = []
    teok = []
    tspeed = []
    tnospeed = []
    tnomatch = []
    for p in pathList:
        path = _ch(p) + "/output"
        groundtruth = []
        plates_moto = []
        plates_car = []
        radar = []
        detection = []
        detection_moto = []
        detection_car = []
        lInfo = laneInfo(_ch(p) + "/lanes.xml")
        groundtruth += parseLogFiles(path + "/groundtruth/f1", lInfo)
        groundtruth += parseLogFiles(path + "/groundtruth/f2", lInfo)
        groundtruth += parseLogFiles(path + "/groundtruth/f3", lInfo)
        detection += parseLogFiles(path + "/detection/f1", lInfo)
        detection += parseLogFiles(path + "/detection/f2", lInfo)
        detection += parseLogFiles(path + "/detection/f3", lInfo)
        for r in groundtruth:
            if r.plate:
                if r.moto:
                    plates_moto.append(r)
                else:
                    plates_car.append(r)
            if r.plate and r.radar and not r.sema:
               radar.append(r)
        for r in detection:
            if r.moto:
                detection_moto.append(r)
            else:
                detection_car.append(r)
        printColor(bcolors.HEADER, 'Processing video %s' % _ch(p))
        print ('  Length of ground truth: %d' % len(groundtruth))
        print ('  Length of detection:    %d' % len(detection))
        print ('  Total license plates:   %d' % len(plates_moto + plates_car))
        preA,recA = computeDetectionPASCAL(plates_moto + plates_car, detection, path, writeLog=True, auxText="Total")
        if DETECT_MOTORCYCLE:
            preB,recB = computeDetectionPASCAL(plates_moto, detection_moto, path, auxText="Moto")
            preC,recC = computeDetectionPASCAL(plates_car, detection_car, path, auxText="Car")
        tplates_moto += plates_moto
        tplates_car += plates_car
        tdetection += detection
        tdetection_moto += detection_moto
        tdetection_car += detection_car
        tpre += preA
        trec += recA
        if DETECT_MOTORCYCLE:
            tpre_moto += preB
            trec_moto += recB
            tpre_car += recC
            trec_car += recC
        print ('  Radar ground truth:     %d' % len(radar))
        if SIMUL_PERFECT_DETECTOR:
            speed,nospeed,nomatch = computeSpeed(radar, radar, path)
        else:
            speed,nospeed,nomatch = computeSpeed(radar, detection, path)
        tspeed += speed
        ehigh,elow,eok = classifySpeedError_USA(speed, path)
        tradar += radar
        tehigh += ehigh
        telow += elow
        teok += eok
        tnospeed += nospeed
        tnomatch += nomatch
        print ('  Radar OK:               %d' % len(eok))
        print ('  Radar HIGH:             %d' % len(ehigh))
        print ('  Radar LOW:              %d\n' % len(elow))
        print ('  Radar no Speed:         %d' % len(nospeed))
        print ('  Radar no Match:         %d\n' % len(nomatch))
    print ('# --------------------------\n')
    printColor(bcolors.OKBLUE,   'Total license plates:   %d' % len(tplates_moto + tplates_car))
    if len(tdetection) > 0:
        v = float(len(tpre)) / float(len(tdetection))
    else:
        v = float('NaN')
    printColor(bcolors.OKGREEN,  'Precision (Total):      %.3f' % (v))
    if len(tplates_moto + tplates_car) > 0:
        v = float(len(trec)) / float(len(tplates_moto + tplates_car))
    else:
        v = float('NaN')
    printColor(bcolors.OKGREEN,  'Recall    (Total):      %.3f' % (v))
    if DETECT_MOTORCYCLE:
        if len(tdetection_moto) > 0:
            v = float(len(tpre_moto)) / float(len(tdetection_moto))
        else:
            v = float('NaN')
        printColor(bcolors.OKGREEN,  'Precision (Moto):       %.3f' % (v))
        if len(tplates_moto) > 0:
            v = float(len(trec_moto)) / float(len(tplates_moto))
        else:
            v = float('NaN')
        printColor(bcolors.OKGREEN,  'Recall    (Moto):       %.3f' % (v))
        if len(tdetection_car) > 0:
            v = float(len(tpre_car)) / float(len(tdetection_car))
        else:
            v = float('NaN')
        printColor(bcolors.OKGREEN,  'Precision (Car):        %.3f' % (v))
        if len(tplates_car) > 0:
            v = float(len(trec_car)) / float(len(tplates_car))
        else:
            v = float('NaN')
        printColor(bcolors.OKGREEN,  'Recall    (Car):        %.3f' % (v))
    print ('')
    printColor(bcolors.OKBLUE,   'Total Radar ground t.:  %d' % len(tradar))
    printColor(bcolors.OKBLUE,   'Total Speed OK:         %d' % len(teok))
    printColor(bcolors.OKBLUE,   'Total Speed HIGH:       %d' % len(tehigh))
    printColor(bcolors.OKBLUE,   'Total Speed LOW:        %d\n' % len(telow))
    printColor(bcolors.WARNING,  'Total no Speed:         %d' % len(tnospeed))
    printColor(bcolors.WARNING,  'Total no Match:         %d\n' % len(tnomatch))
    print ('')
    printColor(bcolors.OKBLUE,   'Total Speed: %d' % len(tspeed))
    saveSpeedHostigram(tspeed, 'speed-hist.eps')
    #tspeed.sort(key=lambda x: x.radar_speed, reverse=False)
    #saveSpeedCurve(tspeed, 'speed-graph.dat')

'''free -h

@function: _signal_handler
@brief: Exit if reach here
'''
def _signal_handler(signal, frame):
    print (bcolors.WARNING + 'Reached signal handler, exiting' + bcolors.ENDC)
    sys.exit(0)

'''
@function: _usage
@brief: Helper usage function
'''
def _usage(progname):
    printColor(bcolors.WARNING, "Usage: " + progname + " <video-output>")

''' Entry point
=============== '''
if __name__ == "__main__":
    signal.signal(signal.SIGABRT, _signal_handler)
    signal.signal(signal.SIGINT, _signal_handler)
    signal.signal(signal.SIGTERM, _signal_handler)
    try:
        if len(sys.argv) < 2:
            _usage(sys.argv[0])
            raise Exception(bcolors.FAIL + "Bad arguments: " + str(sys.argv))
        main(sys.argv[1:])
    except Exception as e:
       print (str(e) + bcolors.ENDC)
       sys.exit(1)
