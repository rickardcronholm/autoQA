#!/usr/bin/env python

# reduce VOI

import sys
import dicom
import numpy as np
import CTCtools
import copy
from scipy import ndimage

class struct:

    def __init__(self):
        pass


def getFiltSize(RD, dist):
    rd = dicom.read_file(RD)
    return int(dist / (np.average(map(float, rd.PixelSpacing))/10))


def getReducedBool(RP, RD, RS):
    # input
    # RP - DICOM RT Plan
    # RD - DICOM RT Dose
    # RS - DICOM RT Struct

    # set some variables
    confFile = '/home/mcqa/MCQA/samc.conf'
    filtDist = 1 # cm
    filtSize = getFiltSize(RD, filtDist)
    with open(confFile) as f:
        confCont = f.readlines()

    # get the variables from confFile
    whatToGet = ['CTC_auto.extName', 'CTC_auto.suppName']
    cv = struct()
    cv = CTCtools.getConfVars(cv, whatToGet, confCont)

    # get x and y coords
    rd_xmesh, rd_ymesh = CTCtools.getDICOMcoords(RD, True)  # get DICOMcoords in cm
    # get z coords
    rd_zmesh = CTCtools.getDICOMzCoord(RD, True)


    # get orientation of CT and RD
    #orientCT = CTCtools.getOrient(ct[0])
    orientRD = CTCtools.getOrient(RD)
    orientCT = orientRD
    # not sure how to handle this. I'll need to test
    # For now. assume no rotation is needed. Will do for RD(HFS) and CT(HFS)

    # generate matrix
    mtrx = np.zeros((len(rd_zmesh), len(rd_ymesh), len(rd_xmesh)), dtype=float)

    # Read RS file
    RS = dicom.read_file(RS)  # open file
    # index ROIs
    # & replace empty strucutre types with specified string
    refROIs = []
    refContSeq = []
    replaceType = 'NONE'
    for i in range(0, len(RS.ROIContourSequence)):
	try:
        	refROIs.append(int(RS.ROIContourSequence[i].ReferencedROINumber))
        	refContSeq.append(CTCtools.getCorrContSeq(RS.ROIContourSequence, refROIs[-1]))
        	if len(RS.RTROIObservationsSequence[i].RTROIInterpretedType) == 0:
	            RS.RTROIObservationsSequence[i].RTROIInterpretedType = replaceType
	except AttributeError:
		pass


    # Correlate structures between ROIContourSequence and RTROIObservationsSequence
    refObsSeq = []
    for elem in refROIs:
        refObsSeq.append(CTCtools.getCorrContSeq(RS.RTROIObservationsSequence, elem))

    # create and init structures
    structures = []
    structureShell = struct()

    # start collecting structure data
    for elem in refObsSeq:
        try:
            structure = copy.deepcopy(structureShell)
            # structure.init_name()

            structure.name = RS.RTROIObservationsSequence[elem].ROIObservationLabel
            structure.type = RS.RTROIObservationsSequence[elem].RTROIInterpretedType
            structure.refObsSeq = elem

            structures.append(structure)
        except AttributeError:
            pass

    for i in range(0, len(structures)):
        # disregard support structures
        if not structures[i].type == cv.suppName:
            # get contour
            print structures[i].name
            try:
                structures[i].contour = CTCtools.getContour(RS.ROIContourSequence[structures[i].refObsSeq], rd_zmesh, abs(orientRD[1]), True)
                # deInterpolate contour onto dose grid and generate boolean matrix
                structures[i].logicMatrix = CTCtools.interpStructToDose(structures[i].contour, rd_xmesh, rd_ymesh, rd_zmesh, rd_xmesh, rd_ymesh, rd_zmesh)
                if structures[i].type == cv.extName:  # add && not releELec !!!
                    filt = np.ones((filtSize,filtSize,filtSize)).astype(int)
                    indx = np.where(ndimage.convolve(structures[i].logicMatrix,filt,mode='nearest') >= filtSize**3)
                    structures[i].logicMatrix = np.zeros(structures[i].logicMatrix.shape).astype(int)
                    structures[i].logicMatrix[indx] = 1
            except AttributeError:
                structures[i].logicMatrix = np.zeros(mtrx.shape, dtype = float)

        else:
            structures[i].logicMatrix = np.zeros(mtrx.shape, dtype = float)

    # combine global mtrx
    for i in range(0, len(structures)):
        mtrx = np.add(mtrx, structures[i].logicMatrix)
    mtrx = np.where(mtrx >= 1, 1, mtrx)

    #import matplotlib.pyplot as plt
    #plt.figure(1)
    #plt.imshow(mtrx[30][:][:])
    #plt.show()

    return mtrx
