#!/usr/bin/python2.7

import os, sys
import shutil
import signal
import glob
import numpy

import matplotlib.pyplot as plt

from math import sqrt
from collections import Counter

import xml.dom.minidom as minidom
import xml.etree.cElementTree as cet
import xml.etree.ElementTree as et

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

class ComputeSample():
    def __init__(self, xmlSample):
        self.index = int(xmlSample.attrib.get('_0'))
        self.id = int(xmlSample.attrib.get('_1'))
        self.frame = int(xmlSample.attrib.get('_2'))
        self.nlines = int(xmlSample.attrib.get('_3'))
        self.radar = float(xmlSample.attrib.get('_4'))
        self.speed = float(xmlSample.attrib.get('_5'))
        self.err = float(xmlSample.attrib.get('_6'))
        self.stdspeed = float(xmlSample.attrib.get('_7'))

class ComputeLane():
    def __init__(self, xmlLane):
        self.xml = xmlLane
        self.samples = []
        self.parseLane()
        self.computeErrList()
    
    def parseLane(self):
        for child in self.xml:
            if child.tag == 'InputFiles':
                self.nInputFiles = int(child.attrib.get('num'))
                self.parseInputFiles(child)
            if child.tag == 'Samples':
                for spl in child:
                    if spl.tag == 'Spl':
                        self.samples.append(ComputeSample(spl))
        printColor(bcolors.HEADER, '  found ' + str(len(self.samples)) + ' valid samples')
        print ('   nInputFiles = ' + str(self.nInputFiles))
        print ('   nEmpty = ' + str(self.nEmpty))
        print ('   nInvalidReal = ' + str(self.nInvalidReal))
        print ('   nNotTracked = ' + str(self.nNotTracked))
    
    def parseInputFiles(self, xml_input_files):
        for child in xml_input_files:
            if child.tag == 'Empty':
                self.nEmpty = int(child.attrib.get('num'))
            elif child.tag == 'InvalidReal':
                self.nInvalidReal = int(child.attrib.get('num'))
            elif child.tag == 'Valid':
                self.parseValid(child)
                
    def parseValid(self, xml_valid):
        self.nNotTracked = 0
        for parent in xml_valid:
            if (parent.tag == 'Plate') or (parent.tag == 'NoPlate'):
                for child in parent:
                    if child.tag == 'NotTracked':
                        self.nNotTracked += int(child.attrib.get('num'))
    
    def computeErrList(self):
        self.errlist = []
        for s in self.samples:
            self.errlist.append(s.err)

class ComputeAnalysis():
    def __init__(self, fpath):
        if os.path.isfile(fpath) is False:
            raise Exception(bcolors.FAIL, "Error: file '" + fpath + "' does not exist!")
        printColor(bcolors.OKBLUE, "Info: processing file " + fpath)
        self.xmlLanes = {}
        self.objLanes = {}
        self.tree = cet.parse(fpath)
        self.root = self.tree.getroot()
        self.parseXML()
    
    def parseXML(self):
        for child in self.root:
            if child.tag == "Video":
                self.vnum = child.attrib.get('num')
                for l in child:
                    if str(l.tag) == 'Lane':
                        lname = int(l.attrib.get('name').split()[1])
                        self.xmlLanes[lname] = l
        for l in self.xmlLanes:
            printColor(bcolors.HEADER, 'Parsing Lane ' + str(l))
            self.objLanes[l] = ComputeLane(self.xmlLanes[l])
    
    def getErrList(self, lnum):
        aux = []
        for i in lnum:
            aux += self.objLanes[i].errlist
        return aux
    
    def getNumInputFiles(self, lnum):
        aux = 0
        for i in lnum:
            aux += self.objLanes[i].nInputFiles
        return aux
    
    def getNumEmpty(self, lnum):
        aux = 0
        for i in lnum:
            aux += self.objLanes[i].nEmpty
        return aux
    
    def getNumInvalidReal(self, lnum):
        aux = 0
        for i in lnum:
            aux += self.objLanes[i].nInvalidReal
        return aux
    
    def getNumNotTracked(self, lnum):
        aux = 0
        for i in lnum:
            aux += self.objLanes[i].nNotTracked
        return aux

def qmean(num):
    return sqrt(sum(n*n for n in num)/len(num))

class ComputeGeneral():
    def __init__(self):
        self.comps = []
    
    def newComputeAnalysis(self, comp):
        self.comps.append(comp)
    
    def reportGeneral(self, filename):
        errlist = []
        nInputFiles = 0
        nEmpty = 0
        nInvalidReal = 0
        nNotTracked = 0
        for cp in self.comps:
            errlist += cp.getErrList([1, 2, 3])
            nInputFiles += cp.getNumInputFiles([1, 2, 3])
            nEmpty += cp.getNumEmpty([1, 2, 3])
            nInvalidReal += cp.getNumInvalidReal([1, 2, 3])
            nNotTracked += cp.getNumNotTracked([1, 2, 3])
        logfile = open(filename, 'w')
        logfile.write('General Error Report:\n')
        logfile.write('\n----------------------------------------------\n')
        logfile.write('  Average error:             %.2f\n' % numpy.mean(errlist))
        logfile.write('  Average err (POSITIVE):    %.2f\n' % numpy.mean(filter(lambda a: a>=0, errlist)))
        logfile.write('  Average err (NEGATIVE):    %.2f\n' % numpy.mean(filter(lambda a: a<0, errlist)))
        logfile.write('  Root mean square error:    %.2f\n' % qmean(errlist))
        logfile.write('  Standard deviation:        %.2f\n' % numpy.std(errlist))
        logfile.write('  Maximum error:             %.2f\n' % max(errlist))
        logfile.write('  Minimum error:             %.2f\n' % min(errlist))
        logfile.write('\n----------------------------------------------\n')
        logfile.write('  Total Input Files:         %d' % nInputFiles + ' \t%.2f%%\n' % 100.0)
        TotalKLTError = nEmpty + nNotTracked
        logfile.write('  * Total KLT Error:         %d' % TotalKLTError + ' \t%.2f%% *=\n' % (100.0 * TotalKLTError / nInputFiles))
        logfile.write('    + Total Empty:           %d' % nEmpty + ' \t%.2f%%\n' % (100.0 * nEmpty / nInputFiles))
        logfile.write('    + Total Not Tracked:     %d' % nNotTracked + ' \t%.2f%%\n' % (100.0 * nNotTracked / nInputFiles))
        logfile.write('  Total Invalid Real:        %d' % nInvalidReal + ' \t%.2f%%\n' % (100.0 * nInvalidReal / nInputFiles))
        logfile.write('\nError counting according to INMETRO/BR\n')
        logfile.write('-----------------------------------------------\n')
        errLen = len(errlist)
        nGreater = sum(i > 2.0 for i in errlist)
        nSmaller = sum(i < -3.0 for i in errlist)
        nOk = errLen - nGreater - nSmaller
        percGreater = 100.0 * float(nGreater) / float(errLen)
        percSmaller = 100.0 * float(nSmaller) / float(errLen)
        percOn = 100.0 * float(nOk) / float(errLen)
        logfile.write('  **Only tracked samples:\n')
        logfile.write('  Error greater than +2km/h: %d' % nGreater + '\t%.2f%%\n' % percGreater)
        logfile.write('  Error smaller than -3km/h: %d' % nSmaller + '\t%.2f%%\n' % percSmaller)
        logfile.write('  Inside error limits:       %d' % nOk + '\t%.2f%%\n' % percOn)
        TotalKLT = errLen + TotalKLTError
        percGreater = 100.0 * float(nGreater) / float(TotalKLT)
        percSmaller = 100.0 * float(nSmaller) / float(TotalKLT)
        percOn = 100.0 * float(nOk) / float(TotalKLT)
        percKlt = 100.0 * float(TotalKLTError) / float(TotalKLT)
        logfile.write('\n')
        logfile.write('  **With KLT Errors:\n')
        logfile.write('  Error greater than +2km/h: %d' % nGreater + '\t%.2f%%\n' % percGreater)
        logfile.write('  Error smaller than -3km/h: %d' % nSmaller + '\t%.2f%%\n' % percSmaller)
        logfile.write('  KLT error:                 %d' % TotalKLTError + '\t%.2f%%\n' % percKlt)
        logfile.write('  Inside error limits:       %d' % nOk + '\t%.2f%%\n' % percOn)
        logfile.write('  Total Vehicles:            %d' % TotalKLT + '\t%.2f%%\n' % 100.0)
        logfile.close()
        print ('\n    Final results written in "%s"\n' % filename)
        neg = numpy.arange(-0.5, min(errlist) - 0.5, -0.5)
        pos = numpy.arange(0, max(errlist) + 0.5, 0.5)
        bins = neg.tolist()[::-1] + pos.tolist()
        plt.hist(errlist, bins)
        plt.axis([-7, 7, 0, 1200])
        plt.grid(True)
        plt.title('Error histogram')
        plt.xlabel('Error [km/h]')
        plt.ylabel('Number of vehicles')
        plt.show()

"""----------------------------------------------------------------------------
# Entry point """
def main(listOfFiles):
    cg = ComputeGeneral()
    for fpath in listOfFiles:
        cg.newComputeAnalysis(ComputeAnalysis(_ch(fpath)))
    cg.reportGeneral('results.txt')

def _signal_handler(signal, frame):
        sys.exit(0)

def _usage(progname):
    printColor(bcolors.WARNING, "Usage: " + progname + " <list-of analyses.xml files>")

if __name__ == "__main__":
    signal.signal(signal.SIGABRT, _signal_handler)
    signal.signal(signal.SIGINT, _signal_handler)
    signal.signal(signal.SIGTERM, _signal_handler)
    try:
        print 'DEBUG: ' + str(len(sys.argv))
        if len(sys.argv) < 2:
            _usage(sys.argv[0])
            raise Exception(bcolors.FAIL + "Bad arguments: " + str(sys.argv))
        main(sys.argv[1:])
    except Exception as e:
       print (str(e) + bcolors.ENDC)
       sys.exit(1)
