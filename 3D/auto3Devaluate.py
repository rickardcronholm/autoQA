#!/usr/bin/python

# pyhton script that automatically evaluates two 3D-dose distributions 
# using 3D-gamma (npgamma) and delta-DVH


class SettingsGamma:
    def __init__(self, **kwds):
        self.__dict__.update(kwds)
    def get_settings(self, confFile):
        variables = ['doseDiff', 'dta', 'local', 'cutOff', 'cut', 'stepSize', 'timeFile', 'dirToLook', 'startsWith', 'gammaPrefix', 'planPrefix', 'dosePrefix', 'timeStampLength', 'python', 'scriptName', 'homogeneousTreshold', 'heteroDens', 'heteroCorr', 'passThreshold']
        for name in variables:
            par = samcUtil.getVariable(confFile,name)
            setattr(self, name, par)
        return self
    def castType(self):
        variables = ['timeStampLength', 'local', 'cut', 'cutOff', 'doseDiff', 'dta', 'stepSize', 'homogeneousTreshold', 'passThreshold']
        types = [int, int, int, float, float, float, float, int, float]
        for i in range(0, len(variables)):
            setattr(self, variables[i], map(types[i], [getattr(self, variables[i])])[0])
        variables = ['local', 'cut', 'homogeneousTreshold']
        types = [bool, bool, bool]
        for i in range(0, len(variables)):
            setattr(self, variables[i], map(types[i], [getattr(self, variables[i])])[0])
        return self
        
        
class SettingsDVH:
    def __init__(self, **kwds):
        self.__dict__.update(kwds)
    def get_settings(self):
        self.eval_levels = [0.97, 0.5, 0.03]
        self.eval_threshold = [0.02, 0.01, 0.02]
        self.resolution =  .01  # resoultion in Gy
        return self


class DoseItem:
    def __init__(self, filename):
        self.filename = filename
        ds = dicom.read_file(filename)
        self.xCoord, self.yCoord = CTCtools.getDICOMcoords(filename, False)
        self.zCoord = CTCtools.getDICOMzCoord(filename, False)
        pixels = np.transpose(ds.pixel_array, (1, 2, 0))
        self.dose = np.ndarray.astype(pixels, 'int') * float(ds.DoseGridScaling)
        
        
class Settings:
    def __init__(self):
        gammaConf = '/home/mcqa/MCQA/autoGamma.conf'
        self.gamma = SettingsGamma()
        self.gamma = self.gamma.get_settings(gammaConf)
        self.gamma = self.gamma.castType()
        self.DVH = SettingsDVH()
        self.DVH = self.DVH.get_settings()
        
        
class StructItem:
    def __init__(self, Name, Type, RefROInr):
        self.Name = Name
        self.Type = Type
        self.Contour = None
        self.RefROInr = RefROInr


class TimedOutExc(Exception):
    pass


def deadline(timeout, *args, **kwargs):
    def decorate(f):
        def handler(signum, frame):
            raise TimedOutExc()

        def new_f(*args, **kwargs):
            signal.signal(signal.SIGALRM, handler)
            signal.alarm(timeout)
            return f(*args, **kwargs)
            signa.alarm(0)

        new_f.__name__ = f.__name__
        return new_f
    return decorate
    

class db_settings():
    def __init__(self):
        self.database = ''
        self.user = ''
        self.host = ''
        self.password = ''
        self.table = ['patient', 'rtplan', 'field', 'controlitem', 'control3ditem', 'controlstructitem']
        self.createString = ['(patientser serial PRIMARY KEY, patientariaser integer)',
            '(patientser serial REFERENCES patient, rtplanser serial PRIMARY KEY, sopinstanceuid text, creationdate timestamp)',
            '(rtplanser serial REFERENCES rtplan, fieldser serial PRIMARY KEY, referencedbeamnumber integer, treatmentmachine text, nominalenergy integer, fluencemode text)',
            '(fieldser integer, rtplanser integer, controlser serial PRIMARY KEY, controltype integer, controlpassed boolean, status text, creationdate timestamp)',
            '(controlser serial REFERENCES controlitem, control3dser serial PRIMARY KEY, dosecrit double precision, distcrit double precision, dosethreshold double precision, localgamma boolean, passingrate double precision, gamma bytea)',
            '(control3dser serial REFERENCES control3ditem, controlstructser serial PRIMARY KEY, structname text, structtype text, controlpassed boolean, status text, passingrate double precision, gamma bytea, deltad1 double precision, deltad2 double precision, deltad3 double precision)']
        #self.table = ['calibration', 'patient', 'rtplan', 'field', 'controlitem', 'control2ditem', 'backscattercorrection']
        #self.createString = ['(calibrationser serial PRIMARY KEY, treatmentmachine text, nominalenergy integer, fluencemode text, calibrationfactor double precision, creationdate timestamp)',
        #    '(patientser serial PRIMARY KEY, patientariaser integer)',
        #    '(patientser serial REFERENCES patient, rtplanser serial PRIMARY KEY, sopinstanceuid text, creationdate timestamp)',
        #    '(rtplanser serial REFERENCES rtplan, fieldser serial PRIMARY KEY, referencedbeamnumber integer, treatmentmachine text, nominalenergy integer, fluencemode text)',
        #    '(fieldser serial REFERENCES field, controlser serial PRIMARY KEY, controltype integer, controlpassed boolean, status text, creationdate timestamp)',
        #    '(controlser serial REFERENCES controlitem, control2dser serial PRIMARY KEY, alignx double precision, aligny double precision, dosecrit double precision, distcrit double precision, dosethreshold double precision, localgamma boolean, passingrate double precision, gamma bytea, calibrationser serial REFERENCES calibration)',
        #    '(correctionser serial PRIMARY KEY, treatmentmachine text, nominalenergy integer, fluencemode text, coefficient12 double precision, coefficient11 double precision, coefficient10 double precision, coefficient22 double precision, coefficient21 double precision, coefficient20 double precision, creationdate timestamp)']

class evaluation:
    def __init__(self, filename):
        ds = dicom.read_file(filename)
        self.PatientID = ds.PatientID
        self.rtplan = ds.SOPInstanceUID
        self.timestamp = datetime.datetime.now().strftime('%Y%m%d%H%M%S')
        self.ds = ds
        
        
import reduceVOI
import CTCtools
import dicom
import numpy as np
import samcUtil
import glob
import os
import shutil
import mynpgamma
import copy
import datetime
import sys
from progressbar import ProgressBar
from skimage import draw
import pylab as P
import psycopg2
import pymssql
import matplotlib.pyplot as plt
from scipy.stats import entropy
import traceback
import signal
import struct


def get_from_aria(value, item):
    server = ""
    user = ""
    passwd = ""
    db = ""
    conn = pymssql.connect(server, user, passwd, db)
    cur = conn.cursor()
    if item == 'prescription':
        sqlStr = "SELECT TotalDose=CONVERT(float, PrescriptionAnatomyItem.ItemValue) FROM RTPlan INNER JOIN PlanSetup ON PlanSetup.PlanSetupSer=RTPlan.PlanSetupSer INNER JOIN PrescriptionAnatomy ON PrescriptionAnatomy.PrescriptionSer=PlanSetup.PrescriptionSer INNER JOIN PrescriptionAnatomyItem ON PrescriptionAnatomyItem.PrescriptionAnatomySer=PrescriptionAnatomy.PrescriptionAnatomySer WHERE RTPlan.PlanUID='{:s}' AND PrescriptionAnatomyItem.ItemType='TOTAL DOSE' ORDER BY CONVERT(float, PrescriptionAnatomyItem.ItemValue) DESC".format(value)  
    elif item == 'plan':
        sqlStr = "SELECT TDose=(SUM(RP.FieldDose) * RP.NumberOfFractions), FDose=SUM(RP.FieldDose), RP.NumberOfFractions  FROM (SELECT RadiationRefPoint.*, Prescription.NumberOfFractions FROM RTPlan INNER JOIN PlanSetup ON PlanSetup.PlanSetupSer=RTPlan.PlanSetupSer INNER JOIN RefPoint ON RefPoint.PatientVolumeSer=PlanSetup.PrimaryPTVSer INNER JOIN RadiationRefPoint ON RadiationRefPoint.RefPointSer=RefPoint.RefPointSer INNER JOIN Prescription ON Prescription.PrescriptionSer=PlanSetup.PrescriptionSer WHERE RTPlan.PlanUID='{:s}' AND DoseSpecificationFlag=1 ) AS RP INNER JOIN RTPlan ON RTPlan.RTPlanSer=RP.RTPlanSer WHERE RTPlan.PlanUID='{:s}' GROUP BY RP.NumberOfFractions".format(value, value) 
    elif item == 'PatientSer':
        sqlStr = "SELECT PatientSer FROM Patient WHERE PatientId='{:s}'".format(value)
    cur.execute(sqlStr)
    row = cur.fetchone()
    cur.close()
    conn.close()
    return row[0]

@deadline(1200)
def evaluate(directory, gammaEval=True, DVHEval=True, entropyEval=False, log=True, plot=False):
    settings = Settings()
    refRD = glob.glob(os.path.sep.join([directory, 'RD.1.*']))[0]
    evlRD = glob.glob(os.path.sep.join([directory, 'RD.MC*']))[0]
    refRS = glob.glob(os.path.sep.join([directory, 'RS.1.*']))[0]
    refRP = glob.glob(os.path.sep.join([directory, 'RP.1.*']))[0]
    evalItem = evaluation(refRP)
    structures = read_structures(refRS)
    ref = DoseItem(refRD)
    evl = DoseItem(evlRD)
    dose_level, settings.gamma.cut = get_dose_level(settings.gamma.cut, ref, directory)
    
    if gammaEval: 
        gamma = eval_gamma(ref, evl, settings.gamma, dose_level)
        # run reduceVOI to get rid of fringe
        try:
            redVOI = np.transpose(reduceVOI.getReducedBool(refRP, refRD, refRS), (1, 2, 0))
        except IndexError:
            redVOI = np.ones(gamma.dose.shape)
        try:
            gamma.dose = np.where(redVOI > 0, gamma.dose, np.nan)
        except ValueError:
            pass
        gamma.pass_per_struct, evalItem.structures = report_gamma(gamma, structures, plot)
        gamma = crop_item(gamma)
        gamma.valid_gamma = gamma.dose[~np.isnan(gamma.dose)]
        gamma.valid_gamma[gamma.valid_gamma > 2] = 2
        gamma.passing_rate = 100. * np.sum(gamma.valid_gamma <= 1) / float(len(gamma.valid_gamma))
        evalItem.passing_rate = gamma.passing_rate
        evalItem.binaryGamma = gamma2binary(gamma.valid_gamma, 2., 51)
        
        cast_gamma_as_DICOM(refRP, refRD, gamma, settings.gamma, directory)
        #return gamma
    else:
        evalItem.passing_rate = 91.
        
    if DVHEval or entropyEval:
        ds = dicom.read_file(refRP)
        prescribed_dose = get_from_aria(ds.SOPInstanceUID, 'plan')
        #print prescribed_dose
        structrues = report_DVH(ref, evl, structures, settings.DVH, prescribed_dose, DVHEval=True, entropyEval=False, plot=plot)
        
    if log:
        # log to db
        DBsettings = db_settings()
        # get patientser from ARIA
        evalItem.patariaser = get_from_aria(evalItem.PatientID, 'PatientSer')
        # transfer settings to evalItem
        evalItem.dosethres = settings.gamma.doseDiff
        evalItem.dist_threshold = settings.gamma.dta
        evalItem.ldc = settings.gamma.cutOff
        evalItem.local_dose = settings.gamma.local
        # create table patient if not exist
        check_create_table('patient', DBsettings)
        # check if patient in db, otherwise add
        if not check_if_exist('patient', evalItem.patariaser, DBsettings):
            insert('patient', evalItem, DBsettings)
        # create table rtplan if not exist
        check_create_table('rtplan', DBsettings)
        # check if rtplan in db otherwise add
        if not check_if_exist('rtplan', evalItem.rtplan, DBsettings):
            evalItem.patser = get_from_table('patient', evalItem.patariaser, DBsettings)
            insert('rtplan', evalItem, DBsettings)
        evalItem.rtplanser = get_from_table('rtplan', evalItem.rtplan, DBsettings)
        #cur.execute("INSERT INTO field (rtplanser, referencedbeamnumber, treatmentmachine, nominalenergy, fluencemode) 
        #(evalItem.rtplanser, evalItem.refbeamnr, evalItem.MachineName, evalItem.NominalEnergy, evalItem.FluenceMode))
        # create table field if not exist
        check_create_table('field', DBsettings)
        for beam in evalItem.ds.BeamSequence:
            if beam.TreatmentDeliveryType == 'TREATMENT':
                evalItem.MachineName = beam.TreatmentMachineName
                evalItem.NominalEnergy = beam.ControlPointSequence[0].NominalBeamEnergy
                evalItem.refbeamnr = beam.BeamName
                evalItem.FluenceMode = beam.PrimaryFluenceModeSequence[0].FluenceMode
                insert('field', evalItem, DBsettings)
        # create table controlitem if not exist
        check_create_table('controlitem', DBsettings)
        if hasattr(evalItem, 'passing_rate'):
            evalItem.passed = evalItem.passing_rate > settings.gamma.passThreshold
            # insert
            evalItem.controlser = insert('controlitem', evalItem, DBsettings)
            # create table control3ditem if not exist
            check_create_table('control3ditem', DBsettings)
            # insert
            evalItem.control3dser = insert('control3ditem', evalItem, DBsettings)
        # create table control3ditem if not exist
        check_create_table('controlstructitem', DBsettings)
        for structure in structures:
            if hasattr(structure, 'passing_rate'):
                evalStruct = copy.deepcopy(evalItem)
                evalStruct.structName = structure.Name
                evalStruct.structType = structure.Type
                evalStruct.passing_rate = structure.passing_rate
                evalStruct.passed = evalStruct.passing_rate > settings.gamma.passThreshold
                evalStruct.binaryGamma = structure.binaryGamma
                try:
                    evalStruct.deltad = structure.deltad
                except AttributeError:
                    evalStruct.deltad = [np.nan, np.nan, np.nan]
                insert('controlstructitem', evalStruct, DBsettings)
        

def get_dose_level(cut, ref, directory):
    dose_level = np.max(ref.dose)
    if not cut:
        try:
            ds = dicom.read_file(glob.glob(os.path.sep.join([directory, 'RP.1.*']))[0])
            dose_level = get_from_aria(ds.SOPInstanceUID, 'prescription')
        except TypeError:
            cut = True
    return dose_level, cut
            
    


def eval_gamma(ref, evl, gammaSettings, dose_level):
    try: 
        ref = DoseItem(ref)
    except AttributeError:
        pass
    try: 
        evl = DoseItem(evl)
    except AttributeError:
        pass
    
    # recompute doseDiff
    doseDiff = gammaSettings.doseDiff * 0.01  # from percent to relative
    if not gammaSettings.local:
        doseDiff *= dose_level

    # cut off relative to dose_level
    cut_here = gammaSettings.cutOff * 0.01 * dose_level
    
    # set parameters for gamma evaluation
    ref.coords = (ref.yCoord, ref.xCoord, ref.zCoord)
    evl.coords = (evl.yCoord, evl.xCoord, evl.zCoord)
    distance_step_size = gammaSettings.dta / gammaSettings.stepSize
    maximum_test_distance = gammaSettings.dta * 2
    max_concurrent_calc_points = np.inf
    num_threads = 2
    
    #print dose_level, gammaSettings.dta, doseDiff, cut_here, gammaSettings.local, gammaSettings.cutOff
    
    gammaMatrix = mynpgamma.calc_gamma(
        ref.coords, ref.dose,
        evl.coords, evl.dose,
        gammaSettings.dta, doseDiff,
        lower_dose_cutoff=cut_here, 
        distance_step_size=distance_step_size,
        maximum_test_distance=maximum_test_distance,
        local_dose=gammaSettings.local,
        max_concurrent_calc_points=max_concurrent_calc_points,
        num_threads=num_threads)

    gamma = copy.deepcopy(evl)
    gamma.dose = gammaMatrix
    return gamma


def volume_below(array, val):
    array = array[~np.isnan(array)]
    return np.sum(array <= val) / float(len(array))
    
    
def dose_at_volume(array, val):
    array = array[~np.isnan(array)]
    v = np.linspace(0., 1., num=len(array))
    try:
        return np.interp(1-val, v, np.sort(array))  # note 1-val
    except ValueError:
        return np.nan

    
def report_gamma(gammaObject, structures, plot):
    bins = np.linspace(0., 2., num=41)
    valid_types = ['ORGAN', 'PTV', 'CTV', 'GTV', 'EXTERNAL']
    print 'TOTAL', ': ', volume_below(np.ravel(gammaObject.dose), 1.0) * 100.0
    gamma_s = []
    for structure in structures:
        if structure.Type in valid_types:
            array = get_values_for_struct(gammaObject, structure.Contour)
            array = array[~np.isnan(array)]
            vol_pr = volume_below(array, 1.0) * 100.0
            structure.passing_rate = vol_pr
            structure.binaryGamma = gamma2binary(array, 2., 51)
            if plot:
                #ar = array[~np.isnan(array)]
                plt.hist(array, bins=bins, normed=True, histtype='step', label='{:s}: {:.2f}%'.format(structure.Name, vol_pr))
            # print structure.Name , ': ', vol_pr
            gamma_s.append([structure.Name, vol_pr])
    # sort: worst first
    gamma_s.sort(key=lambda x: x[1])
    gamma_s = '\r\n'.join(['{:s}: {:.2f}'.format(x[0], x[1]) for x in gamma_s])
    if plot:
        plt.legend()
        plt.show()
    return gamma_s, structures
    

def gamma2binary(array, maxGamma, nBins):
    try:
        array[array > maxGamma] = maxGamma  # crop to maxGamma
        bins = np.linspace(0, maxGamma, nBins)
        hist, _ = np.histogram(array, bins=bins)
        hist = hist.astype(float)
        # normalize histogram (so that sum=1)
        hist /= np.sum(hist)
        bin_centers = (bins[:-1] + bins[1:]) / 2
        binary_data = struct.pack('i', nBins-1)
        for bc in bin_centers:
            binary_data += struct.pack('f', bc)
        for h in hist:
            binary_data += struct.pack('f', h)
        return binary_data
    except (ValueError, TypeError):
        return ''
            
    
def cast_gamma_as_DICOM(templateRP, templateRD, gamma, gammaSettings, directory):
    gammaRP = samcUtil.rewriteRP(dicom.read_file(templateRP), 'gamma')
    # remove verification points
    gammaRP.DoseReferenceSequence = dicom.sequence.Sequence()
    # info into RTPlanDescription
    if gammaSettings.local:
        localStr = 'Local dose'
    else:
        localStr = 'Global dose'
    gammaRP.RTPlanDescription += '\r\nPassing rate: {0:.2f}\r\nCriteria: {1:.1f}%/{2:.1f}mm ({3:s})\r\n'.format(gamma.passing_rate, gammaSettings.doseDiff, gammaSettings.dta, localStr)
    gammaRP.RTPlanDescription += gamma.pass_per_struct
    gammaRD = cast_gamma_as_RD(dicom.read_file(templateRD), gamma)
    gammaRD.ReferencedRTPlanSequence[0].ReferencedSOPInstanceUID = gammaRP.SOPInstanceUID
    gammaRP.save_as(os.path.join(directory, 'gammaRP.' + directory + '.dcm'))
    gammaRD.save_as(os.path.join(directory, 'gammaRD.' + directory + '.dcm'))
    # copy files to lookForDir
    shutil.copy(os.path.join(directory, 'gammaRP.' + directory + '.dcm'), gammaSettings.dirToLook)
    shutil.copy(os.path.join(directory, 'gammaRD.' + directory + '.dcm'), gammaSettings.dirToLook)
    


def cast_gamma_as_RD(dcm, gamma):
	#print gamma.dose.shape
	gamma.dose = np.transpose(gamma.dose, (2, 0, 1))
	if int(dcm.Rows) != int(gamma.dose.shape[1]) or int(dcm.Columns) != int(gamma.dose.shape[2]) or int(dcm.NumberOfFrames) != int(gamma.dose.shape[0]):
		dcm.Rows = gamma.dose.shape[1]
		dcm.Columns = gamma.dose.shape[2]
		dcm.NumberOfFrames = gamma.dose.shape[0]
		dcm.ImagePositionPatient = map(str, [gamma.xCoord[0], gamma.yCoord[0], gamma.zCoord[0]])
		dcm.GridFrameOffsetVector = map(str, [x - gamma.zCoord[0] for x in gamma.zCoord.tolist()])

	# set dose type and pixel representation
	dcm.DoseType = 'ERROR'
	#dcm.DoseUnits = 'RELATIVE'
	# dcm.PixelRepresentation = 1
	# change the VR from US to SS for all tags in the same group as PixelRepresentation
	# prGroup = dcm.data_element('PixelRepresentation').tag.group
	# items = dcm.dir()
	#for i in range(0, len(items)):
	#	try:
	#		if dcm.data_element(items[i]).tag.group == prGroup and dcm.data_element(items[i]).VR == 'US':
	#			dcm.data_element(items[i]).VR = 'SS'
	#	except AttributeError:
	#		pass
	# convert to uint32 using dcm.DoseGridScaling
	unwVal = 0.0 # value to set unwanted elements to
	gamma.dose[gamma.dose < 0] = unwVal  # get rid of negative values
	gamma.dose[np.isnan(gamma.dose)] = unwVal # set inf to value
	gamma.dose = (gamma.dose/float(dcm.DoseGridScaling)).astype('uint32')

	# PixelData
	dcm.PixelData = gamma.dose.tostring()
	
	#update UID
	dcm.SOPInstanceUID = dicom.UID.generate_uid()

	return dcm

def read_structures(RSfile):
    RS = dicom.read_file(RSfile)
    structures = []
    for elem in RS.RTROIObservationsSequence:
        try:
            structures.append(StructItem(getattr(elem, 'ROIObservationLabel'), 
                getattr(elem, 'RTROIInterpretedType'),
                getattr(elem, 'ReferencedROINumber')))
        except AttributeError:
            pass
    # get array of RefROInr
    order_index = [x.RefROInr for x in structures]
    for elem in RS.ROIContourSequence:
        try:
            # determine which struct through ROInr
            structures[order_index.index(elem.ReferencedROINumber)].Contour = getContour(elem.ContourSequence)
        except (ValueError, AttributeError):
            pass
    
    # remove "empty" structures
    structures = [x for x in structures if x.Contour]
    
    return structures


def getContour(RSData):
    # Check that the first item is not 'POINT'
    if RSData[0].ContourGeometricType == 'POINT':
        return None
    contour = []
    for i in range(0, len(RSData)):
        rawCont = map(float, RSData[i].ContourData)
        
        rawCont = np.asarray(rawCont).reshape(len(rawCont) / 3, 3)
        data = np.transpose(rawCont)
        xCont = data[0]
        yCont = data[1]
        zCont = data[2]
        # add the first element to form closed loop
        xCont = np.append(xCont, xCont[0])
        yCont = np.append(yCont, yCont[0])
        zCont = np.append(zCont, zCont[0])
        cont = np.vstack((yCont, xCont, zCont))  # note order changed
        
        #if i == 0 or RSData[i].ContourImageSequence[0].ReferencedSOPInstanceUID != RSData[i-1].ContourImageSequence[0].ReferencedSOPInstanceUID:
        #    segments = []
        segments = []
        points = cont[:]
        segments.append(points)
        contour.append(segments)
        
    return contour


def get_values_for_struct(doseObject, Contour):

    x = doseObject.xCoord
    y = doseObject.yCoord
    z = doseObject.zCoord
    # start by generating a mask based on dose
    maskM = np.zeros(doseObject.dose.shape)
    xL = np.arange(len(x))
    yL = np.arange(len(y))
    pbar = ProgressBar()
    try:
        for j in range(0, len(Contour)):
            # locate closest slice
            sliceNr = np.argwhere(np.around(z, decimals=1) ==
                np.around(Contour[j][0][2][0], decimals=1))[0][0]
            for i in range(0, len(Contour[j])):
                # take care of unbound contours
                points_x = copy.deepcopy(Contour[j][i][1][:])  # get x-points
                points_y = copy.deepcopy(Contour[j][i][0][:])  # get y-points
                points_x[np.where(points_x < min(x))] = min(x)
                points_x[np.where(points_x > max(x))] = max(x)
                points_y[np.where(points_y < min(y))] = min(y)
                points_y[np.where(points_y > max(y))] = max(y)
                # interpolate to voxel number instead of absolute coord
                points_x = np.interp(points_x, x, xL)
                points_y = np.interp(points_y, y, yL)
                # assign mask based on polygons
                [rr, cc] = draw.polygon(np.asarray(points_y),
                np.asarray(points_x), (len(y), len(x)))
                tempMask = np.zeros((len(y), len(x)))
                tempMask[rr, cc] = -1
                # add for the current slice
                maskM[:, :, sliceNr] = np.add(maskM[:, :, sliceNr],
                tempMask[:, :])
                maskM[:, :, sliceNr] = np.abs(maskM[:, :, sliceNr])

    except (IndexError, TypeError):
        pass
    try:
        del tempMask
    except UnboundLocalError:
        pass
    maskM = np.where(maskM > 1, 1, maskM)  # remove duplicates
    
    values = np.ravel(doseObject.dose)
    mask = np.ravel(maskM)
    return values[np.where(mask == 1)[0]]
    
    
def report_DVH(ref, evl, structures, DVHsettings, prescribed_dose, DVHEval=True, entropyEval=False, plot=False):
    valid_types = ['ORGAN', 'PTV', 'CTV', 'GTV', 'EXTERNAL']
    #print np.max([ref.dose, evl.dose])
    dmax = np.max([ref.dose, evl.dose])
    bins = np.linspace(0., dmax, num=int(dmax/DVHsettings.resolution))
    c_bins = np.convolve(bins, np.ones((2,))/2, mode='valid')
    cum_bins = np.linspace(0., dmax, num=int(dmax/0.01))
    for structure in structures:
        if structure.Type in valid_types:
            ref_arr = get_values_for_struct(ref, structure.Contour)
            ref_arr = ref_arr[~np.isnan(ref_arr)]
            evl_arr = get_values_for_struct(evl, structure.Contour)
            evl_arr = evl_arr[~np.isnan(evl_arr)]
            if entropyEval:
                print structure.Name

                #ref_hist, _, _ = P.hist(ref_arr, bins=bins, normed=1, histtype='bar', cumulative=False)
                #evl_hist, _, _ = P.hist(evl_arr, bins=bins, normed=1, histtype='bar', cumulative=False)
                ref_hist, _ = np.histogram(ref_arr, bins=bins, density=True)
                evl_hist, _ = np.histogram(evl_arr, bins=bins, density=True)
                #print np.sum(ref_hist), np.sum(evl_hist)
                #print np.sum(ref_hist), np.sum(evl_hist)
                entropy_kl = kl_entropy(ref_hist, qk=evl_hist)
                #P.show()
                #plt.clf()
                print entropy_kl
                # compute homogeneity index
                hi = dose_at_volume(ref_arr, .98) / dose_at_volume(ref_arr, .02)
                if plot:
                    plt.subplot(2,1,1)
                    plt.plot(c_bins, ref_hist, label='reference')
                    plt.plot(c_bins, evl_hist, label='evaluation')
                    try:
                        plt.title('{:s}: {:.4f} (HI: {:.4f})'.format(structure.Name, entropy_kl, hi))
                    except UnicodeDecodeError:
                        plt.title('{:s}: {:.4f} (HI: {:.4f})'.format('', entropy_kl, hi))
                    plt.legend()
                    plt.subplot(2,1,2)
                    ref_hist, _, _ = P.hist(ref_arr, bins=cum_bins, normed=1, histtype='step', cumulative=-1, label='reference')
                    evl_hist, _, _ = P.hist(evl_arr, bins=cum_bins, normed=1, histtype='step', cumulative=-1, label='evaluation')
                    plt.legend()
                    plt.show()
               
            if DVHEval:
                #print structure.Name
                structure.deltad = []
                for i in range(0, len(DVHsettings.eval_levels)):
                    #print DVHsettings.eval_levels[i], DVHsettings.eval_threshold[i]
                    r = (dose_at_volume(ref_arr, DVHsettings.eval_levels[i]))
                    e = (dose_at_volume(evl_arr, DVHsettings.eval_levels[i]))
                    #print 'Ref: ', r
                    #print 'Evl: ', e
                    #print 'Passed: ', np.abs(r-e) < (DVHsettings.eval_threshold[i] * 100)
                    structure.deltad.append(np.abs(r-e)/prescribed_dose * 100.)
                    if structure.Name == 'Rectum':
                        print r, e, prescribed_dose

    return structures


def kl_entropy(pk, qk):
    # compute entropy of p
    hp = entropy(pk)
    hq = entropy(qk)
    #S = sum(pk * log(pk / qk), axis=0)
    pk /= np.sum(pk)
    qk /= np.sum(qk)
    dk = np.divide(pk, qk)
    dk[np.isinf(dk)] = np.nan
    sk = np.multiply(pk, np.log(dk))
    sk[np.isnan(sk)] = 0.
    return np.sum(sk)/(hp + hq)


def crop_item(item):
    item.dose[np.isnan(item.dose)] = -1
    rel_z = (np.max(item.dose, axis=(0, 1)) > 0)
    z_start = np.max([np.where(rel_z)[0][0],0])
    z_end = np.max([np.where(rel_z)[0][-1],0]) + 1
    rel_x = (np.max(item.dose, axis=(0, 2)) > 0)
    x_start = np.max([np.where(rel_x)[0][0],0])
    x_end = np.max([np.where(rel_x)[0][-1],0]) + 1
    rel_y = (np.max(item.dose, axis=(1, 2)) > 0)
    y_start = np.max([np.where(rel_y)[0][0],0])
    y_end = np.max([np.where(rel_y)[0][-1],0]) + 1
    item.yCoord = item.yCoord[y_start:y_end]
    item.xCoord = item.xCoord[x_start:x_end]
    item.zCoord = item.zCoord[z_start:z_end]
    item.dose = item.dose[y_start:y_end, x_start:x_end, z_start:z_end]
    item.dose[np.where(item.dose < 0)] = np.nan
    return item


def logit(ts, string):
    with open(os.path.join(os.getcwd(), ts + '.log'), 'w') as flock:
        flock.write('autoEvaluate error\n{:s}'.format(string))
    shutil.copy(os.path.join(os.getcwd(), ts + '.log'), 
        os.path.join(os.path.join(os.getcwd(), ts, ts + '.gamma')))
        
        
def check_create_table(table_name, settings):
    connString = "dbname={:s} user={:s} host={:s} password={:s}".format(settings.database, settings.user, settings.host, settings.password)
    conn = psycopg2.connect(connString)
    cur = conn.cursor()
    try:
        ind = settings.table.index(table_name)
        sqlString = "CREATE TABLE IF NOT EXISTS {:s} {:s};".format(settings.table[ind], settings.createString[ind])
        cur.execute(sqlString)
    except ValueError:
        return
    conn.commit()
    cur.close()
    conn.close()


def check_if_exist(table_name, value, settings):
    connString = "dbname={:s} user={:s} host={:s} password={:s}".format(settings.database, settings.user, settings.host, settings.password)
    conn = psycopg2.connect(connString)
    cur = conn.cursor()
    if table_name == 'calibration':
        cur.execute("SELECT * FROM calibration WHERE calibrationser={:d}".format(value))
        row = cur.fetchone()
    elif table_name == 'patient':
        cur.execute("SELECT * FROM patient WHERE patientariaser={:d}".format(value))
        row = cur.fetchone()
    elif table_name == 'rtplan':
        cur.execute("SELECT * FROM rtplan WHERE sopinstanceuid='{:s}'".format(value))
        row = cur.fetchone()
    elif table_name == 'field':
        cur.execute("SELECT * FROM field WHERE rtplanser={:d} and referencedbeamnumber={:d}".format(value.rtplanser, value.refbeamnr))
        row = cur.fetchone()
    cur.close()
    conn.close()
    if row is None:
        return False
    else:
        return True
        

def get_from_table(table, value, settings):
    connString = "dbname={:s} user={:s} host={:s} password={:s}".format(settings.database, settings.user, settings.host, settings.password)
    conn = psycopg2.connect(connString)
    cur = conn.cursor()
    if table == 'patient':
        cur.execute("SELECT patientser FROM patient WHERE patientariaser={:d}".format(value))
        row = cur.fetchone()
    elif table == 'rtplan':
        cur.execute("SELECT rtplanser FROM rtplan WHERE sopinstanceuid='{:s}'".format(value))
        row = cur.fetchone()
    elif table == 'field':
        cur.execute("SELECT field.fieldser FROM field INNER JOIN rtplan ON (field.rtplanser = rtplan.rtplanser) WHERE rtplan.sopinstanceuid='{:s}' AND field.referencedbeamnumber={:d}".format(value.rtplan, value.refbeamnr))
        row = cur.fetchone()
    return row[0]


def insert(table, evalItem, settings):
    returnValue = None
    connString = "dbname={:s} user={:s} host={:s} password={:s}".format(settings.database, settings.user, settings.host, settings.password)
    conn = psycopg2.connect(connString)
    cur = conn.cursor()
    if table == "calibration":
        cur.execute("INSERT INTO calibration (treatmentmachine, nominalenergy, fluencemode, calibrationFactor, creationdate) VALUES (%s, %s, %s, %s, to_timestamp(%s,'YYYYMMDDHH24MISSMS'))", (evalItem.MachineName, evalItem.NominalEnergy, evalItem.FluenceMode, evalItem.CalibrationFactor, evalItem.timestamp))
    elif table == "patient":
        cur.execute("INSERT INTO patient (patientariaser) VALUES (%s)", ([evalItem.patariaser]))
    elif table == "rtplan":
        cur.execute("INSERT INTO rtplan (patientser, sopinstanceuid) VALUES (%s, %s)", (evalItem.patser, evalItem.rtplan))
    elif table == "field":
        cur.execute("INSERT INTO field (rtplanser, referencedbeamnumber, treatmentmachine, nominalenergy, fluencemode) VALUES (%s, %s, %s, %s, %s)", (evalItem.rtplanser, evalItem.refbeamnr, evalItem.MachineName, evalItem.NominalEnergy, evalItem.FluenceMode))
    elif table == "controlitem":
        cur.execute("INSERT INTO controlitem (rtplanser, controltype, controlpassed, status, creationdate) VALUES (%s, %s, %s, %s, to_timestamp(%s,'YYYYMMDDHH24MISSMS')) RETURNING controlser", (evalItem.rtplanser, 0, str(evalItem.passed).lower(), 'UNAPROVVED', evalItem.timestamp))
        returnValue = cur.fetchone()[0]
    elif table == "control3ditem":
        if hasattr(evalItem, 'binaryGamma') and evalItem.binaryGamma:
            sqlStr = "INSERT INTO control3ditem (controlser, dosecrit, distcrit, dosethreshold, localgamma, passingrate, gamma) VALUES ({:d}, {:.3f}, {:.2f}, {:.3f}, {:s}, {:.2f}, {:s}) RETURNING control3dser".format(evalItem.controlser, evalItem.dosethres, evalItem.dist_threshold, evalItem.ldc, str(evalItem.local_dose).lower(), evalItem.passing_rate, psycopg2.Binary(evalItem.binaryGamma))
        else:
            sqlStr = "INSERT INTO control3ditem (controlser, dosecrit, distcrit, dosethreshold, localgamma, passingrate) VALUES ({:d}, {:.3f}, {:.2f}, {:.3f}, {:s}, {:.2f}) RETURNING control3dser".format(evalItem.controlser, evalItem.dosethres, evalItem.dist_threshold, evalItem.ldc, str(evalItem.local_dose).lower(), evalItem.passing_rate)
        cur.execute(sqlStr.replace('nan', 'NULL'))
        returnValue = cur.fetchone()[0]
    elif table == "controlstructitem":
        if hasattr(evalItem, 'binaryGamma') and evalItem.binaryGamma:
            sqlStr = "INSERT INTO controlstructitem (control3dser, structname, structtype, controlpassed, status, passingrate, gamma, deltad1, deltad2, deltad3) VALUES ({:d}, '{:s}', '{:s}', {:s}, '{:s}', {:2f}, {:s}, {:2f}, {:2f}, {:2f})".format(evalItem.control3dser, str(evalItem.structName), str(evalItem.structType), str(evalItem.passed).lower(), 'UNAPROVVED', evalItem.passing_rate, psycopg2.Binary(evalItem.binaryGamma), evalItem.deltad[0], evalItem.deltad[1], evalItem.deltad[2])
        else:
            sqlStr = "INSERT INTO controlstructitem (control3dser, structname, structtype, controlpassed, status, passingrate, deltad1, deltad2, deltad3) VALUES ({:d}, '{:s}', '{:s}', {:s}, '{:s}', {:2f}, {:2f}, {:2f}, {:2f})".format(evalItem.control3dser, evalItem.structName, evalItem.structType, str(evalItem.passed).lower(), 'UNAPROVVED', evalItem.passing_rate, evalItem.deltad[0], evalItem.deltad[1], evalItem.deltad[2])
        cur.execute(sqlStr.replace('nan', 'NULL'))
    elif table == "control2ditem":
        sqlStr = "INSERT INTO control2ditem (controlser, alignx, aligny, dosecrit, distcrit, dosethreshold, localgamma, passingrate, gamma, calibrationser) VALUES ({:d}, {:.4f}, {:.4f}, {:.3f}, {:.2f}, {:.3f}, {:s}, {:.2f}, {:s}, {:d})".format(evalItem.controlser, evalItem.align_matrix[1], evalItem.align_matrix[0], evalItem.dosethres, evalItem.dist_threshold, evalItem.ldc, str(evalItem.local_dose).lower(), evalItem.passing_rate, psycopg2.Binary(evalItem.binaryGamma), evalItem.CalibrationSer)
        cur.execute(sqlStr)
    elif table == 'backscattercorrection':
        sqlStr = "INSERT INTO backscattercorrection (treatmentmachine, nominalenergy, fluencemode, coefficient12, coefficient11, coefficient10, coefficient22, coefficient21, coefficient20, creationdate) VALUES ('{:s}', {:d}, '{:s}', {:.4e}, {:.4e}, {:.4e}, {:.4e}, {:.4e}, {:.4e}, to_timestamp('{:s}','YYYYMMDDHH24MISSMS'))".format(evalItem.MachineName, evalItem.NominalEnergy, evalItem.FluenceMode, evalItem.correctionItem.coeff1[0], evalItem.correctionItem.coeff1[1], evalItem.correctionItem.coeff1[2], evalItem.correctionItem.coeff2[0], evalItem.correctionItem.coeff2[1], evalItem.correctionItem.coeff2[2], evalItem.timestamp)
        cur.execute(sqlStr)
    conn.commit()
    cur.close()
    conn.close()
    return returnValue
    

def main():
    # get current time
    now = datetime.datetime.now()
    newTime = now.strftime("%y%m%d%H%M%S")

    # change dir to where the script is
    abspath = os.path.abspath(__file__)
    dname = os.path.dirname(abspath)
    os.chdir(dname)

    # check if lock file exists, otherwise create it
    lockFile = 'ag.lock'
    lockCheck = [lockFile, 'De.lock']
    lock = list(set([os.path.isfile(x) for x in lockCheck]))
    if len(lock) > 1 or lock[0]:
        # print 'Busy, quitting'
        sys.exit(0)
    else:
        flock = open(lockFile,'w')
        flock.write('Busy\n')
        flock.close()

    # check if DICOMoutDir exists, otherwise exit
    DICOMoutDir = '/home/mcqa/MCQA/VA_TRANSFER/MonteCarlo/simulated/'
    if not os.path.exists(DICOMoutDir):
        print 'Not able to connect to DICOMoutDir'
        # do not update timeFile
        os.remove(lockFile)  # remove lockFile
        fid = open(lockFile.replace('.lock', '.jam'), 'w')
        fid.write('Not able to connect to DICOMoutDir\n')
        fid.close()
        sys.exit(0)

    # set variables
    settings = Settings()
    
    # get list of timeStamps
    timeStamps = []
    for file in [fileList for fileList in os.listdir(settings.gamma.dirToLook) if fileList.startswith(settings.gamma.startsWith)]:
        # extract timeStamp
        timeStamp = file.lstrip(settings.gamma.startsWith)[:settings.gamma.timeStampLength]
        #print timeStamp
        # check if gamma already exists for the timeStamp in question
        if len(glob.glob(''.join([settings.gamma.dirToLook, settings.gamma.gammaPrefix, '*', timeStamp, '*']))) == 0:
            # get time of modification of file
            #fileTime = os.path.getmtime(settings.dirToLook + file)
    		#fileTime = time.localtime(fileTime)
            #fileTime = time.strftime("%y%m%d%H%M%S",fileTime)
            # if more recent than last run, append to list
            #if fileTime > oldTime:
            # make sure that the "parent" directory exists
            # and that it has not been flagged as failure
            if os.path.isdir(os.path.sep.join([dname,timeStamp])) and not os.path.isfile(os.path.join(dname, timeStamp, timeStamp + '.gamma')):
                timeStamps.append(timeStamp)
    #print timeStamps            
    if len(timeStamps) > 0:
        i = 0
        repeat = True
        while repeat:
            try:
                ts = timeStamps[i]
            except IndexError:
                break
            print 'timeStamp: ', ts
            sys.stdout.flush()

            # evaluate
            try:
                evaluate(ts, gammaEval=True, DVHEval=True, entropyEval=False)
                repeat = False
            except TimedOutExc:
                i += 1
                #print traceback.format_exc()
                logit(ts, traceback.format_exc())
            except Exception as e:
                i += 1
                #print traceback.format_exc()
                logit(ts, traceback.format_exc())
	

    # update time of last run
    #samcUtil.writeTimeToFile(settings.timeFile,newTime)

    # remove lockFile
    os.remove(lockFile)

# Function chooser
func_arg = {"-main": main}
# run specifc function as called from command line
if __name__ == "__main__":
    if sys.argv[1] == "-main":
        func_arg[sys.argv[1]]()
