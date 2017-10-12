#!/usr/bin/env python
"""
delayFinder.py
===========

This file calculates the delay and phase offsets for SWARM from a MIR-format data file.
"""

import os, sys, struct, mmap, argparse
import makevis, pydsm, numpy
from math import sin, cos, pi, sqrt, atan2, pi, log

neededFiles = ('bl_read', 'in_read', 'sch_read', 'sp_read')

def csv(value):
    return map(int, value.split(","))

targetRx = -1
targetRxName = '230'
rxClass = 'RxA'
rxList = []
firstScan = -1
lastScan = 100000000
trialOnly = False
verbose = False
noiseSource = False
antList = range(1,9)
nChannels = 16384
designBitcodeSpeed = 4.576
fixRx = False # Set True to force the script to ignore receiver
receiverName = 230
bandList = []
spSmallDictL = {}
spSmallDictU = {}
spBigDict = {}
antennas = {}
inDict = {}
blDictU = {}
blDictL = {}
numberOfBaselines = 0
antennaList = []
sidebandList = [0, 1]
silent = False
totalPoints = 0
totalZeros = 0
twoPi = 2.0*pi
aveDelay = {}
aveOffset = {}
#avePhase = {}
refAnt = 2
gotPhases = False
padDict = {}
loadFromTable = False

def dMod180(ang):
    while ang > 180.0:
        ang -= 360.0
    while ang < -180.0:
        ang += 180.0
    return ang

def calculateDelays():
    global twoPi, refAnt

    antennaList = []
    for key in delayData:
        if not key[0] in antennaList:
            antennaList.append(key[0])
        if not key[1] in antennaList:
            antennaList.append(key[1])
    antennaList = sorted(antennaList)
    if refAnt not in antennaList:
        print 'No data for reference antenna ', refAnt
        sys.exit(-1)
    sBList = []
    for key in delayData:
        if not key[3] in sBList:
            sBList.append(key[3])
    sBList = sorted(sBList)
    bandList = []
    for key in delayData:
        if not key[2] in bandList:
            bandList.append(key[2])
    bandList = sorted(bandList)
    for band in bandList:
        for sB in sBList:
            if sB == 0:
                sBName = 'LSB'
            else:
                sBName = 'USB'
            if not silent:
                print 'For band %d, %s: (reference ant %d)' % (band, sBName, refAnt)
            for ant1 in antennaList:
                if (ant1 != refAnt) and (ant1 in antList):
                    if (ant1, rxClass, band) not in aveDelay:
                        aveDelay[(ant1, rxClass, band)] = 0.0
                    if (ant1, rxClass, band) not in aveOffset:
                        aveOffset[(ant1, rxClass, band)] = 0.0
                    bestRMS = 1.0e30
                    if refAnt < ant1:
                        sign = -1.0
                        key = (refAnt, ant1, band, sB)
                    else:
                        sign = 1.0
                        key = (ant1, refAnt, band, sB)
                    thisData = delayData[key]
                    iData = [cos(thisData[i]) + sin(thisData[i])*1j for i in xrange(len(thisData))]
                    iTrans=numpy.fft.fft(iData)
                    delA=[abs(iTrans[i]) for i in xrange(len(iTrans))]
                    guess = numpy.argmax(delA)
                    snr = delA[numpy.argmax(delA)]/numpy.std(delA)
                    if guess > len(delA)/2:
                        guess = len(delA)-guess
                    else:
                        guess = -guess
                    if guess == 0:
                        smallSize = 32
                    else:
                        smallSize = 2**(int(log(abs(guess))/log(2.0))+4)
                    if smallSize > nChannels:
                        smallSize = nChannels
                    nToAverage = nChannels/smallSize
                    smallIData = [sum(iData[nToAverage*i:nToAverage*(i+1)]) for i in range(nChannels/nToAverage)]
                    smallPData = [atan2(smallIData[i].imag, smallIData[i].real) for i in range(len(smallIData))]
                    smallN = len(smallPData)
                    inc = twoPi/float(smallN)
                    trim = max(1, 1024/nToAverage)
                    thisRange = range(trim, smallN-trim)
                    n = len(thisRange)
                    for iDelay in xrange(4*(guess-1), 4*(guess+1)):
                        delay = float(iDelay)/4.0
                        delInc = delay*inc
                        for offset in (0.0, pi/2.0):
                            if offset == 0.0:
                                thisDataO = smallPData
                            else:
                                thisDataO = [x + offset for x in smallPData]
                            ave = ((sum(thisDataO) + delInc*sum(thisRange)) % twoPi)/float(n)
                            RMS = sum([(((thisDataO[i] + delInc*float(i)) % twoPi) - ave)**2 for i in thisRange])
                            if RMS < bestRMS:
                                bestRMS = RMS
                                bestDelay = delay
                    bestRMS = 1.0e30
                    for iDelay in xrange(int((bestDelay - 1.0)*300), int((bestDelay + 1.0)*300)):
                        delay = float(iDelay)/300.0
                        delInc = delay*inc
                        for offset in (0.0, pi/2.0):
                            if offset == 0.0:
                                thisDataO = smallPData
                            else:
                                thisDataO = [x + offset for x in smallPData]
                            ave = sum([(thisDataO[i] + delInc*float(i)) % twoPi for i in thisRange])/float(n)
                            RMS = sum([(((thisDataO[i] + delInc*float(i)) % twoPi) - ave)**2 for i in thisRange])
                            if RMS < bestRMS:
                                bestRMS = RMS
                                bestDelay = delay
                                bestOffset = offset
                                realSum = sum([cos(thisDataO[i] + delInc*float(i)) for i in thisRange])
                                imagSum = sum([sin(thisDataO[i] + delInc*float(i)) for i in thisRange])
                                testTheta = atan2(imagSum, realSum)
                    thisDataO = [x + bestOffset for x in thisData]
#                    print 'Ant %d band %d %s %f (%f)' % (ant1, band, sBName, dMod180(-(testTheta-bestOffset)*180.0/pi), bestOffset*180.0/pi)
                    if aveDelay[(ant1, rxClass, band)] == 0.0:
                        aveDelay[(ant1, rxClass, band)] += sign*bestDelay*2.0/GSamplesPerSecond
                    else:
                        aveDelay[(ant1, rxClass, band)] += sign*bestDelay*2.0/GSamplesPerSecond
                        aveDelay[(ant1, rxClass, band)] /= 2.0
                    if sBName == 'USB':
                        sBSign = 1.0
                    else:
                        sBSign = -1.0
                    if aveOffset[(ant1, rxClass, band)] == 0.0:
                        aveOffset[(ant1, rxClass, band)] = sBSign*dMod180(-(testTheta-bestOffset)*180.0/pi)
                    else:
                        oldSin = sin(aveOffset[(ant1, rxClass, band)]*pi/180.0)
                        oldCos = cos(aveOffset[(ant1, rxClass, band)]*pi/180.0)
                        newSin = sBSign*sin(-(testTheta-bestOffset))
                        newCos = sBSign*cos(-(testTheta-bestOffset))
                        aveOffset[(ant1, rxClass, band)] = atan2(oldSin+newSin, oldCos+newCos)*180.0/pi
                    if not silent:
                        n = float(len(thisRange))
                        print 'Best for ant %d: %7.2f wraps %7.2f nsec (RMS = %f)' % (ant1, bestDelay, sign*bestDelay*2.0/GSamplesPerSecond, sqrt(bestRMS/n)*180.0/pi)

def makeInt(data, size):
    tInt = sum([ord(data[i])<<(i<<3) for i in range(size)])
    return tInt

def makeFloat(data):
    return (struct.unpack('f', data[:4]))[0]

def makeDouble(data):
    return (struct.unpack('d', data[:8]))[0]

def read(dataDir):
    global bandList, antennas, inDict, blDictL, blDictU, targetRx, targetRxName, receiverName, fixRx
    global spSmallDictL, spSmallDictU, spBigDict, sourceDict, numberOfBaselines
    global antennaList

    # Check that the directory contains all the required files
    if verbose:
        print 'checking that the needed files exist in ', dataDir
    dirContents = os.listdir(dataDir)
    for needed in neededFiles:
        if not needed in dirContents:
            print "The directory %s does not have a %s file - aborting" % (dataDir, needed)
            sys.exit(-1)

    ###
    ### Read in bl_read
    ###
    if verbose:
        print 'Reading bl_read'
    f = open(dataDir+'/bl_read', 'rb')
    fSize = os.path.getsize(dataDir+'/bl_read')
    blRecLen = 158
    nBlRecords = fSize/blRecLen
    baselineList = []
    receiverList = []
    blSidebandDict = {}
    for rec in range(nBlRecords):
        if ((rec % 10000) == 0) and verbose:
            print '\t processing record %d of %d (%2.0f%% done)' % (rec, nBlRecords, 100.0*float(rec)/float(nBlRecords))
            sys.stdout.flush()
        data = f.read(blRecLen)
        if len(data) == blRecLen:
            blhid     =    makeInt(data[  0:], 4)
            inhid     =    makeInt(data[  4:], 4)
            isb       =    makeInt(data[  8:], 2)
            ipol      =    makeInt(data[ 10:], 2)
            ant1rx    =    makeInt(data[ 12:], 2)
            ant2rx    =    makeInt(data[ 14:], 2)
            pointing  =    makeInt(data[ 16:], 2)
            irec      =    makeInt(data[ 18:], 2)
            u         =  makeFloat(data[ 20:])
            v         =  makeFloat(data[ 24:])
            w         =  makeFloat(data[ 28:])
            prbl      =  makeFloat(data[ 32:])
            coh       =  makeFloat(data[ 36:])
            avedhrs   = makeDouble(data[ 40:])
            ampave    =  makeFloat(data[ 48:])
            phaave    =  makeFloat(data[ 52:])
            blsid     =    makeInt(data[ 56:], 4)
            ant1      =    makeInt(data[ 60:], 2)
            ant2      =    makeInt(data[ 62:], 2)
            if irec not in receiverList:
                receiverList.append(irec)
            if (ant1, ant2) not in baselineList:
                baselineList.append((ant1, ant2))
                numberOfBaselines +=1
            if ant1 not in antennaList:
                antennaList.append(ant1)
            if ant2 not in antennaList:
                antennaList.append(ant2)
            blSidebandDict[blhid] = isb
            if isb == 0:
                blDictL[blhid] = (inhid, ipol, ant1rx, ant2rx, pointing, irec, u, v, w, prbl, coh, avedhrs, ampave, phaave,
                                  blsid, ant1, ant2)
            else:
                blDictU[blhid] = (inhid, ipol, ant1rx, ant2rx, pointing, irec, u, v, w, prbl, coh, avedhrs, ampave, phaave,
                                  blsid, ant1, ant2)
    if ((targetRx >= 0) and (targetRx not in receiverList)) and (not fixRx):
        print 'No data was found for the receiver specified (%s) - aborting.' % (targetRxName)
        sys.exit(-1)
    if targetRx >= 0:
        rxList.append(rxClass)
    else:
        if 0 in receiverList:
            rxList.append('RxA')
        if 1 in receiverList:
            rxList.append('RxA')
        if 2 in receiverList:
            rxList.append('RxB')
        if 3 in receiverList:
            rxList.append('RxB')
    nAnts = len(antennaList)
    if verbose:
        print '%d antennas and %d baselines seen in this data set.' % (nAnts, numberOfBaselines)
    if numberOfBaselines != (nAnts*(nAnts-1)/2):
        print 'Number of baselines seen inconsistant with the number of antennas - aborting'
        sys.exit(-1)

    ###
    ### Count the number of spectral bands
    ###
    weightDict = {}
    if verbose:
        print 'Counting the number of bands'
    f = open(dataDir+'/sp_read', 'rb')
    fSize = os.path.getsize(dataDir+'/sp_read')
    spRecLen = 188
    data = f.read(spRecLen)
    nSpRecords = fSize/spRecLen
    for rec in range(nSpRecords):
        if ((rec % 100000) == 0) and verbose:
            print '\t processing record %d of %d (%2.0f%% done)' % (rec, nSpRecords, 100.0*float(rec)/float(nSpRecords))
            sys.stdout.flush()
        data = f.read(spRecLen)
        if len(data) == spRecLen:
            sphid     =    makeInt(data[  0:], 4)
            blhid     =    makeInt(data[  4:], 4)
            inhid     =    makeInt(data[  8:], 4)
            igq       =    makeInt(data[ 12:], 2)
            ipq       =    makeInt(data[ 14:], 2)
            iband     =    makeInt(data[ 16:], 2)
            fsky      = makeDouble(data[ 36:])
            fres      =  makeFloat(data[ 44:])
            wt        =  makeFloat(data[ 84:])
            flags     =    makeInt(data[ 88:], 4)
            nch       =    makeInt(data[ 96:], 2)
            dataoff   =    makeInt(data[100:], 4)
            rfreq     = makeDouble(data[104:])
            try:
                ibandRx = blDictU[blhid][5]
            except:
                ibandRx = blDictL[blhid][5]
            if ((targetRx < 0) or (ibandRx == targetRx)) or fixRx:
                rightRx = True
            else:
                rightRx = False
            if (flags != 0) and (wt > 0):
                wt *= -1.0
            if rightRx:
                spBigDict[(iband, blhid)] = (dataoff, wt)
            try:
                if weightDict[inhid] < wt:
                    weightDict[inhid] = wt
            except KeyError:
                weightDict[inhid] = wt
            if rightRx:
                if iband not in bandList:
                    bandList.append(iband)
                try:
                    if blSidebandDict[blhid] == 0:
                        if iband not in spSmallDictL:
                            spSmallDictL[iband] = (nch, fres*1.0e6, fsky*1.0e9, rfreq*1.0e9)
                    else:
                        if iband not in spSmallDictU:
                            spSmallDictU[iband] = (nch, fres*1.0e6, fsky*1.0e9, rfreq*1.0e9)
                except KeyError:
                    # This exception can occur if the last record to bl_read was not written, because dataCatcher was interrupted
                    wt = -1.0
    bandList.sort()

    ###
    ### Read in the integration header information (in_read)
    ### Produces a disctionary with integration numbers for the
    ### keys, and tuples holding the entry values
    ###
    if verbose:
        print 'Reading in_read'
    f = open(dataDir+'/in_read', 'rb')
    data = f.read()
    inRecLen = 188
    nInRecords = len(data)/inRecLen
    for rec in range(nInRecords):
        traid     =    makeInt(data[rec*inRecLen +   0:], 4)
        inhid     =    makeInt(data[rec*inRecLen +   4:], 4)
        az        =  makeFloat(data[rec*inRecLen +  12:])
        el        =  makeFloat(data[rec*inRecLen +  16:])
        hA        =  makeFloat(data[rec*inRecLen +  20:])
        iut       =    makeInt(data[rec*inRecLen +  24:], 2)
        iref_time =    makeInt(data[rec*inRecLen +  26:], 2)
        dhrs      = makeDouble(data[rec*inRecLen +  28:])
        vc        =  makeFloat(data[rec*inRecLen +  36:])
        sx        = makeDouble(data[rec*inRecLen +  40:])
        sy        = makeDouble(data[rec*inRecLen +  48:])
        sz        = makeDouble(data[rec*inRecLen +  56:])
        rinteg    =  makeFloat(data[rec*inRecLen +  64:])
        proid     =    makeInt(data[rec*inRecLen +  68:], 4)
        souid     =    makeInt(data[rec*inRecLen +  72:], 4)
        isource   =    makeInt(data[rec*inRecLen +  76:], 2)
        ivrad     =    makeInt(data[rec*inRecLen +  78:], 2)
        offx      =  makeFloat(data[rec*inRecLen +  80:])
        offy      =  makeFloat(data[rec*inRecLen +  84:])
        ira       =    makeInt(data[rec*inRecLen +  88:], 2)
        idec      =    makeInt(data[rec*inRecLen +  90:], 2)
        rar       = makeDouble(data[rec*inRecLen +  92:])
        decr      = makeDouble(data[rec*inRecLen + 100:])
        epoch     =  makeFloat(data[rec*inRecLen + 108:])
        size      =  makeFloat(data[rec*inRecLen + 112:])
        inDict[inhid] = (traid, inhid, az, el, hA, iut, iref_time, dhrs, vc, sx, sy, sz,
                         rinteg, proid, souid, isource, ivrad, offx, offy, ira, idec, rar, decr, epoch, size)

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='Derive SWARM delays from a data set.   By default the delays are propagated to the online system.')
group1 = parser.add_mutually_exclusive_group()
group2 = parser.add_mutually_exclusive_group()
parser.add_argument('dataset',           help='The data directory to process (i.e. /sma/rtdata/...)', nargs='?')
parser.add_argument('-a', '--antenna',   help='Antenna(s) to process (comma separated list)',         type=csv)
parser.add_argument('-c', '--chunk',     help='Chunk(s) to process (comma separated list)',           type=csv)
parser.add_argument('-C', '--Clock',     help='Clock frequency (GHz)',                                type=float)
parser.add_argument('-d', '--delays',    help='Print the current delay values and exit',              action='store_true')
parser.add_argument('-F', '--FirstScan', help='First scan to process (default = 1)',                  type=int)
group1.add_argument('-l', '--lower',     help='Process lower sideband only',                          action='store_true')
parser.add_argument('-L', '--LastScan',  help='Last scan to process',                                 type=int)
parser.add_argument('-n', '--noise',     help='Process data from the SWARM noise source',             action='store_true')
parser.add_argument('-r', '--receiver',  help='Receiver to use in dual Rx files',                     type=int,choices=[230,345,400,240])
parser.add_argument('-R', '--RefAnt',    help='Referance antenna (default = 2)',                      type=int,choices=antList)
group2.add_argument('-s', '--silent',    help='Run silently unless an error occurs',                  action='store_true')
group2.add_argument('-t', '--trial',     help="Trial run - don't change the delays",                  action='store_true')
parser.add_argument('-T', '--Table',     help='Load the delays from a file')
group1.add_argument('-u', '--upper',     help='Process upper sideband only',                          action='store_true')
group1.add_argument('-v', '--verbose',   help='Be extra chatty',                                      action='store_true')
parser.add_argument('-x', '--cross',     help='Solve for cross receiver delays',                      action='store_true')
parser.add_argument('-z', '--zero',      help='Zero delays and exit',                                 action='store_true')
parser.add_argument('-Z', '--Zero',      help='Zero phases and exit',                                 action='store_true')
args = parser.parse_args()
dataSet = args.dataset
possibleChunks = []

# Figure out how many SWARM quadrants are currently active
activeQuadrants = []
for line in open('/global/projects/SWARMQuadrantsInArray'):
    if '1' in line:
        possibleChunks.append(3)
        possibleChunks.append(4)
        activeQuadrants.append(1)
    if '2' in line:
        possibleChunks.append(1)
        possibleChunks.append(2)
        activeQuadrants.append(2)
    if '3' in line:
        possibleChunks.append(4)
        possibleChunks.append(5)
        activeQuadrants.append(3)
    if '4' in line:
        possibleChunks.append(5)
        possibleChunks.append(6)
        activeQuadrants.append(4)

# Examine the command line arguments
if args.antenna:
    antList = args.antenna
    for ant in antList:
        if (ant < 1) or (ant > 8):
            print 'Illegal antenna number', ant, '- aborting'
            sys.exit(1)
if args.chunk:
    chunkList = args.chunk
    for chunk in chunkList:
        if not (chunk in possibleChunks):
            print 'Illegal chunk number', chunk,' - it must be one of', possibleChunks
            sys.exit(-1)
else:
    chunkList = []
clock = args.Clock
if clock:
    GSamplesPerSecond = clock
if args.FirstScan:
    firstScan = args.FirstScan
if args.Table:
    loadFromTable = True
if args.delays:
    for quad in activeQuadrants:
        for ant in antList:
            roachName = 'roach2-%02d' % (ant+10*(quad-1))
            delays = pydsm.read(roachName, 'SWARM_FIXED_OFFSETS_X')
            delayList = [delays['DELAY_V2_D'][0][0], delays['DELAY_V2_D'][0][1]]
            phaseList = [delays['PHASE_V2_D'][0][0], delays['PHASE_V2_D'][0][1]]
            print '%d %15.11f %15.11f %15.11f %15.11f' % (ant+10*(quad-1), delayList[0], delayList[1], phaseList[0], phaseList[1])
    sys.exit(0)
if args.cross:
    crossDelays = True
else:
    crossDelays = False
if args.lower:
    sidebandList = [0]
if args.LastScan:
    lastScan = args.LastScan
if args.receiver:
    receiverSpecified = True
    if args.receiver == 230:
        targetRx = 0
        targetRxName = '230'
        rxClass = 'RxA'
    elif args.receiver == 345:
        targetRx = 1
        targetRxName = '345'
        rxClass = 'RxA'
    elif args.receiver == 400:
        targetRx = 2
        targetRxName = '400'
        rxClass = 'RxB'
    elif args.receiver == 240:
        targetRx = 3
        targetRxName = '240'
        rxClass = 'RxB'
else:
    receiverSpecified = False
if args.RefAnt:
    refAnt = args.RefAnt
if args.silent:
    silent = True
if args.trial:
    trialOnly = True
if args.verbose:
    verbose = True
if args.noise:
    noiseSource = True
    sidebandList = [1]
    ants2pads = {}
    for ant in range(1, 9):
        ants2pads[ant] = pydsm.rmread('acc%d' % (ant), 'RM_PAD_ID_B')[0]
    noiseDelta = {}
    for line in open('/application/configFiles/noiseDeltas'):
        tok = line.split()
        pad = int(tok[0])
        deltas = [float(tok[i]) for i in range(1, 9)]
        noiseDelta[pad] = deltas
if args.upper:
    sidebandList = [1]
if args.zero:
    if len(chunkList) == 0:
        chunkList = [1, 2, 3, 4]
    newDelays = {}
    for i in activeQuadrants:
        if i in chunkList:
            for j in range(1, 9):
                if j in antList:
                    roachName = 'roach2-%02d' % (j+10*(i-1))
                    if receiverSpecified:
                        oldOffsets=pydsm.read(roachName,  'SWARM_FIXED_OFFSETS_X')
                        if rxClass == 'RxA':
                            delayList = [0.0, oldOffsets['DELAY_V2_D'][0][1]]
                        else:
                            delayList = [oldOffsets['DELAY_V2_D'][0][0], 0.0]
                    else:
                        delayList = [0.0, 0.0]
                    newDelays['DELAY_V2_D'] = delayList
                    pydsm.write(roachName, 'SWARM_FIXED_OFFSETS_X', newDelays)
if args.Zero:
    if len(chunkList) == 0:
        chunkList = [1, 2, 3, 4]
    newPhases = {}
    for i in activeQuadrants:
        if i in chunkList:
            for j in range(1, 9):
                if j in antList:
                    roachName = 'roach2-%02d' % (j+10*(i-1))
                    if receiverSpecified:
                        oldOffsets=pydsm.read(roachName,  'SWARM_FIXED_OFFSETS_X')
                        if rxClass == 'RxA':
                            phaseList = [0.0, oldOffsets['PHASE_V2_D'][0][1]]
                        else:
                            phaseList = [oldOffsets['PHASE_V2_D'][0][0], 0.0]
                    else:
                        phaseList = [0.0, 0.0]
                    newPhases['PHASE_V2_D'] = phaseList
                    pydsm.write(roachName, 'SWARM_FIXED_OFFSETS_X', newPhases)
if args.zero or args.Zero:
    sys.exit(0)

if loadFromTable:
    newOffsets = {}
    fileName = args.Table
    print 'Loading delays and phases from', fileName,
    if receiverSpecified:
        print '(%s only)' % (rxClass)
    else:
        print '(RxA and RxB)'
    for line in open(fileName):
        tok = line.split()
        roach = int(tok[0])
        roachName = 'roach2-%02d' % (roach)
        if (not receiverSpecified) or (rxClass == 'RxA'):
            s1Delay = float(tok[1])
            s1Phase = float(tok[3])
        else:
            oldOffsets=pydsm.read(roachName,  'SWARM_FIXED_OFFSETS_X')
            s1Delay = oldOffsets['DELAY_V2_D'][0][0]
            s1Phase = oldOffsets['PHASE_V2_D'][0][0]
        if (not receiverSpecified) or (rxClass == 'RxB'):
            s2Delay = float(tok[2])
            s2Phase = float(tok[4])
        else:
            oldOffsets=pydsm.read(roachName,  'SWARM_FIXED_OFFSETS_X')
            s2Delay = oldOffsets['DELAY_V2_D'][0][1]
            s2Phase = oldOffsets['PHASE_V2_D'][0][1]
        delayList = [s1Delay, s2Delay]
        newOffsets['DELAY_V2_D'] = delayList
        phaseList = [s1Phase, s2Phase]
        newOffsets['PHASE_V2_D'] = phaseList
        pydsm.write(roachName, 'SWARM_FIXED_OFFSETS_X', newOffsets)
    sys.exit(0)

if noiseSource:
    for ant in range(1, 9):
        if 1 < ants2pads[ant] < 28:
            if noiseDelta[ants2pads[ant]] == [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]:
                print 'Cannot correct the delay for antenna %d on pad %d - no noiseDeltas file entry.' % (ant, ants2pads[ant])
if args.dataset:
    read(dataSet)
else:
    print 'No data set specified - aborting'
    sys.exit(-1)

if crossDelays:
    refAnt = antList[1]
    antList = [antList[0],]
    print 'Using baseline %d-%d to derive cross receiver delays\n' % (antList[0], refAnt)

if len(chunkList) == 0:
    for line in open('/global/projects/SWARMQuadrantsInArray'):
        if '1' in line:
            chunkList.append(1)
        if '2' in line:
            chunkList.append(2)
        if '3' in line:
            chunkList.append(3)
        if '4' in line:
            chunkList.append(4)
for line in open('/global/configFiles/SWARMBitcodeSpeed'):
    tok = line.split()
    if len(tok) > 0:
        nElevenths = int(tok[0]);
GSamplesPerSecond = designBitcodeSpeed*nElevenths/11
visFileLen = makevis.open(dataSet+'/sch_read')
channelRange = range(nChannels)

# Check that the data file has at least some scans for each chunk the user specified
for chunk in chunkList:
    if not chunk in bandList:
        print 'There is no chunk', chunk, 'data in this data file - aborting'
        sys.exit(0)

for band in bandList:
    for sb in range(2):
        delayData = {}
        summedSpectra = {}
        if (band in chunkList) and (sb in sidebandList):
            if sb == 0:
                sideBand = 'Lower'
                blDict = blDictL
                spSmallDict = spSmallDictL
                lowestFSky = spSmallDict[band][2] - spSmallDict[band][1]*0.5 - 52.0e6
            else:
                sideBand = 'Upper'
                blDict = blDictU
                spSmallDict = spSmallDictU
                lowestFSky = spSmallDict[band][2] + spSmallDict[band][1]*0.5 - 52.0e6

            blKeysSorted = sorted(blDict)
            if spSmallDict[band][1] < 0:
                reverseChannelOrder = True
            else:
                reverseChannelOrder = False
            scanOffset = 0
            scanNo = 0
            blPos = 0
            blCount = 0
            nChannels = spSmallDict[band][0]
            firstGoodChannel = 0
            if band == 0:
                lastGoodChannel = 1000000
            else:
                lastGoodChannel  = nChannels-1
            while scanOffset < visFileLen:
                foundBlEntry = False
                foundSpEntry = False
                nBlFound = 0
                scanNo = makevis.scanno(scanOffset, True);
                recSize = makevis.recsize(scanOffset, True)
                if scanNo > 0:
                    if sb == 0:
                        sbName = 'LSB'
                    else:
                        sbName = 'USB'
                    if firstScan <= scanNo <= lastScan:
                        if not silent:
                            print '\rSumming scan %d for band %d %s' % (scanNo, band, sbName),
                            sys.stdout.flush()
                        for bl in blKeysSorted[blPos:]:
                            blCount += 1;
                            if (targetRx < 0) or (blDict[bl][5] == targetRx) or fixRx:
                                rightRx = True
                            else:
                                rightRx = False
                            if ((not crossDelays) and (blDict[bl][0] == scanNo) and rightRx) or \
                                    (crossDelays and (blDict[bl][0] == scanNo) and (blDict[bl][2] < blDict[bl][3])):
                                foundBlEntry = True
                                nBlFound += 1
                                if (band, bl) in spBigDict:
                                    ant1 = blDict[bl][15]
                                    ant2 = blDict[bl][16]
                                    if (ant1 == refAnt) or (ant2 == refAnt):
                                        if not (ant1, ant2, band, sb) in summedSpectra:
                                            summedSpectra[(ant1, ant2, band, sb)] = numpy.zeros(shape=2*nChannels, dtype=numpy.float32)
                                        foundSpEntry = True
                                        matrixEntry  = []
                                        dataoff = scanOffset+spBigDict[(band, bl)][0] + 8
                                        scaleExp = makevis.scaleexp(dataoff, True)
                                        scale = (2.0**scaleExp) * sqrt(twoPi)*260.0
                                        weight = 1.0
                                        matrixEntry = makevis.convert(nChannels, scale, dataoff, True, weight,
                                                                      False, firstGoodChannel, lastGoodChannel,
                                                                      reverseChannelOrder)
                                        key = (ant1, ant2, band, sb)
                                        mArray = numpy.array(matrixEntry);
                                        summedSpectra[key][::2] += mArray[::3]
                                        summedSpectra[key][1::2] += mArray[1::3]

                                    if nBlFound == numberOfBaselines:
                                        break
                    else:
                        foundSpEntry = True
                        foundBlEntry = True
                    blPos += 1
                    if (not foundBlEntry) or (not foundSpEntry):
                        print 'Something not found scanNo = %d, band = %d, foundBl = %d, foundSp = %d' % (scanNo, band, foundBlEntry, foundSpEntry)
                        print '\nThis almost always means that NFS has just not flushed all the data through.'
                        print 'Try running the command again in a few seconds, or after another scan has been stored.'
                        sys.exit(-1)
                scanOffset += recSize+8
            for key in summedSpectra:
                delayData[key] = [atan2(summedSpectra[key][2*i+1], summedSpectra[key][2*i]) for i in channelRange]
            print
            calculateDelays()
newOffsets = {}
if crossDelays:
    ii = antList[0]
    antList = [1,2,3,4,5,6,7,8]
if noiseSource:
    for quad in activeQuadrants:
        key = (refAnt, rxClass, quad)
        aveDelay[key] = 0.0
xsign = 1.0
for quad in activeQuadrants:
    for i in antList:
        if (i != refAnt) or noiseSource or crossDelays:
            roachName = 'roach2-%02d' % (i + 10*(quad-1))
            delays = pydsm.read(roachName, 'SWARM_FIXED_OFFSETS_X')
            delayList = [delays['DELAY_V2_D'][0][0], delays['DELAY_V2_D'][0][1]]
            if not crossDelays:
                ii = i
            if quad == 4:
                key = (ii, 'RxA', 4)
                try:
                    delayList[0] += xsign*aveDelay[key]
                    if noiseSource:
                        delayList[0] -= noiseDelta[ants2pads[i]][6]
                except KeyError:
                    pass
                key = (ii, 'RxB', 4)
                try:
                    delayList[1] += xsign*aveDelay[key]
                    if noiseSource:
                        delayList[1] -= noiseDelta[ants2pads[i]][7]
                except KeyError:
                    pass
            elif quad == 3:
                key = (ii, 'RxA', 3)
                try:
                    delayList[0] += xsign*aveDelay[key]
                    if noiseSource:
                        delayList[0] -= noiseDelta[ants2pads[i]][4]
                except KeyError:
                    pass
                key = (ii, 'RxB', 3)
                try:
                    delayList[1] += xsign*aveDelay[key]
                    if noiseSource:
                        delayList[1] -= noiseDelta[ants2pads[i]][5]
                except KeyError:
                    pass
            elif quad == 2:
                key = (ii, 'RxA', 2)
                try:
                    delayList[0] += xsign*aveDelay[key]
                    if noiseSource:
                        delayList[0] -= noiseDelta[ants2pads[i]][2]
                except KeyError:
                    pass
                key = (ii, 'RxB', 2)
                try:
                    delayList[1] += xsign*aveDelay[key]
                    if noiseSource:
                        delayList[1] -= noiseDelta[ants2pads[i]][3]
                except KeyError:
                    pass
            elif quad == 1:
                key = (ii, 'RxA', 1)
                try:
                    delayList[0] += xsign*aveDelay[key]
                    if noiseSource:
                        delayList[0] -= noiseDelta[ants2pads[i]][0]
                except KeyError:
                    pass
                key = (ii, 'RxB', 1)
                try:
                    delayList[1] += xsign*aveDelay[key]
                    if noiseSource:
                        delayList[1] -= noiseDelta[ants2pads[i]][1]
                except KeyError:
                    pass
            newOffsets['DELAY_V2_D'] = delayList
            if trialOnly:
                if verbose:
                    print 'Aborting without adjusting the delays, as requested - would have written to', roachName, delayList
            else:
                pydsm.write(roachName, 'SWARM_FIXED_OFFSETS_X', newOffsets)

        
