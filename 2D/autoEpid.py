import epidEval
import glob
import struct
import os
import sys
import mmap
import psycopg2
import pymssql
import dicom
import epid
import numpy as np
from scipy import signal
import copy
import time


class db_settings():
    def __init__(self):
        self.table = ['calibration', 'patient', 'rtplan', 'field', 'controlitem', 'control2ditem', 'backscattercorrection']
        self.createString = ['(calibrationser serial PRIMARY KEY, treatmentmachine text, nominalenergy integer, fluencemode text, calibrationfactor double precision, creationdate timestamp)',
            '(patientser serial PRIMARY KEY, patientariaser integer)',
            '(patientser serial REFERENCES patient, rtplanser serial PRIMARY KEY, sopinstanceuid text, creationdate timestamp)',
            '(rtplanser serial REFERENCES rtplan, fieldser serial PRIMARY KEY, referencedbeamnumber integer, treatmentmachine text, nominalenergy integer, fluencemode text)',
            '(fieldser integer, rtplanser integer, controlser serial PRIMARY KEY, controltype integer, controlpassed boolean, status text, creationdate timestamp)',
            '(controlser serial REFERENCES controlitem, control2dser serial PRIMARY KEY, alignx double precision, aligny double precision, dosecrit double precision, distcrit double precision, dosethreshold double precision, localgamma boolean, passingrate double precision, gamma bytea, calibrationser serial REFERENCES calibration)',
            '(correctionser serial PRIMARY KEY, treatmentmachine text, nominalenergy integer, fluencemode text, coefficient12 double precision, coefficient11 double precision, coefficient10 double precision, coefficient22 double precision, coefficient21 double precision, coefficient20 double precision, creationdate timestamp)']


class bscorr_settings():
    def __init__(self):
        self.xPos = 0.0 # perform calibraiton at center
        self.convolveFilt = 5  # use a 5x5 filter for convolution
        self.shift = 5.  # exclude elemets closer than shift mm to field edge
        self.orders = {30: 1, 50: 1, 80: 1, 100: 1, 120: 2, 150: 2, 200: 2, 300: 2}


class backscatt():
    def __init__(self):
        self.fs = []
        self.coeff = []


def generate_lockFile(lockFile):
    with open(lockFile, 'w') as f:
        f.write('Locked\n')


def delete_lockFile(lockFile):
    os.remove(lockFile)


def check_duplicity(string, globExpr):
    fileList = glob.glob(globExpr)
    for item in fileList:
        with open(item) as f:
            s = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
        if s.find(string) != -1:
            return True
    return False


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
        cur.execute("INSERT INTO controlitem (fieldser, controltype, controlpassed, status, creationdate) VALUES (%s, %s, %s, %s, to_timestamp(%s,'YYYYMMDDHH24MISSMS')) RETURNING controlser", (evalItem.fieldser, 1, str(evalItem.passed).lower(), 'UNAPROVVED', evalItem.timestamp))
        returnValue = cur.fetchone()[0]
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


def get_from_aria(item, value):
    conn = pymssql.connect(server, user, passwd, db)
    cur = conn.cursor()
    if item == 'PatientSer':
        cur.execute("SELECT PatientSer FROM Patient WHERE PatientId='{:s}'".format(value))
        row = cur.fetchone()
        return row[0]
        
   
def check_consistency(directory, sleepTime=15):
    new = sum(os.path.getsize(directory) for f in os.listdir(directory)
    if os.path.isfile(os.path.join(directory,f)))
    while True:
        print 'sleeping'
        time.sleep(sleepTime)
        old = copy.deepcopy(new)
        new = sum(os.path.getsize(directory) for f in os.listdir(directory)
        if os.path.isfile(os.path.join(directory,f)))
        if new == old:
            break



def main():
    # define variables
    lockFile = '/home/epidQA/autoEpid.lock'
    inputDir = '/home/epidQA/incoming/'
    inputFileType = '.dat'
    settings = db_settings()
    remove = True

    # check if lock file exists
    if os.path.exists(lockFile):
        print "Busy. Quitting"
        return
    else:
        generate_lockFile(lockFile)

    # check consistency of inputDir
    check_consistency(inputDir)
    # get all dat files
    dats = glob.glob('*'.join([inputDir, inputFileType]))
    dats.sort(key=os.path.getmtime)

    # evaluate
    for dat in dats:
        try:
            evalItem = epidEval.evaluation(dat, inputDir)
            valid = True
        except IOError:
            valid = False
            continue

        # load content of dat
        with open(dat) as f:
            datCont = map(str.rstrip, f.readlines())[:2]

        if hasattr(evalItem, 'CalibrationFactor'):
            # create table if not exist
            check_create_table('calibration', settings)
            # add Calibration to database
            insert('calibration', evalItem, settings)
        elif hasattr(evalItem, 'BackscatterCorrection'):
            evalItem.correctionItem, remove, items = backscatterCorrection(evalItem, dats)
            if remove:
                for item in items:
                    try:
                        os.remove(item)
                        os.remove(item.replace('RI.', '').replace('.dcm', '.dat'))
                    except Exception:
                        pass
                check_create_table('backscattercorrection', settings)
                insert('backscattercorrection', evalItem, settings)
        elif hasattr(evalItem, 'passing_rate'):
            # get patientser from ARIA
            evalItem.patariaser = get_from_aria('PatientSer', evalItem.ID)
            # create table patient if not exist
            check_create_table('patient', settings)
            # check if patient in db, otherwise add
            if not check_if_exist('patient', evalItem.patariaser, settings):
                insert('patient', evalItem, settings)
            # create table rtplan if not exist
            check_create_table('rtplan', settings)
            # check if rtplan in db otherwise add
            if not check_if_exist('rtplan', evalItem.rtplan, settings):
                evalItem.patser = get_from_table('patient', evalItem.patariaser, settings)
                insert('rtplan', evalItem, settings)
            # create table field if not exist
            check_create_table('field', settings)
            # chekc if field in db, otherwise add
            evalItem.rtplanser = get_from_table('rtplan', evalItem.rtplan, settings)
            if not check_if_exist('field', evalItem, settings):
                insert('field', evalItem, settings)
            # create table controlitem if not exist
            check_create_table('controlitem', settings)
            # insert
            evalItem.fieldser = get_from_table('field', evalItem, settings)
            evalItem.controlser = insert('controlitem', evalItem, settings)
            # create table control2ditem if not exist
            check_create_table('control2ditem', settings)
            # insert
            insert('control2ditem', evalItem, settings)
        # remove dat
        
        if remove:
            try:
                os.remove(dat)
            except OSError:
                pass
            for item in datCont:
                if not check_duplicity(item, '*'.join([inputDir, inputFileType])):
                    try:
                        os.remove('.'.join([inputDir + 'RI', item, 'dcm']))
                    except OSError:
                        pass

    # remove lock file once done
    delete_lockFile(lockFile)


def backscatterCorrection(evalItem, datList):
    RIfileList = []
    refBeamNr = []
    items = []
    dirName = os.path.sep.join(datList[0].split(os.path.sep)[:-1])
    fileItems = []
    bs = backscatt()
    remove = False
    for dat in datList:
        try:
            with open(dat) as f:
                RIfileList.append(map(str.rstrip, f.readlines())[0])
        except IOError:
            pass
    #print RIfileList
    for RI in RIfileList:
        #print RI.replace(RI, 'RI.' + RI) + '.dcm'
        fileItem = os.path.join(dirName, RI.replace(RI, 'RI.' + RI)) + '.dcm'
        try: 
            ds = dicom.read_file(fileItem)
            #print ds.SOPInstanceUID
            if ds.ReferencedRTPlanSequence[0].ReferencedSOPInstanceUID == evalItem.rtplan and ds.RadiationMachineName == evalItem.MachineName:
                refBeamNr.append(ds.ReferencedBeamNumber)
                items.append(RI)
                fileItems.append(fileItem)
        except Exception:
            pass
    if len(list(set(refBeamNr))) == 8:  # requires 8 fields
        bs_set = bscorr_settings()
        for im in fileItems:
            #print im
            image = epid.image(im)
            fs = np.around(np.sum(np.abs(map(float, image.ds.ExposureSequence[0].BeamLimitingDeviceSequence[0].LeafJawPositions)))).astype(int)
            coeffs = getCorr(image, conv=bs_set.convolveFilt, shift=bs_set.shift, xPos=bs_set.xPos, order=bs_set.orders[fs])
            bs.fs.append(int(fs/2.))
            if bs_set.orders[fs] == 1:
                bs.coeff.append([np.nan, coeffs[0]])
            else:
                bs.coeff.append(coeffs.tolist())
        non_nans = [y[0] for y in np.argwhere(~np.isnan([x[0] for x in bs.coeff])).tolist()]
        bs.coeff2 = np.polyfit([bs.fs[x] for x in non_nans], [bs.coeff[x][0] for x in non_nans], 2)
        bs.coeff1 = np.polyfit(bs.fs, [x[1] for x in bs.coeff], 2)
        remove = True

    return bs, remove, fileItems



def getCorr(image, conv=0, xPos=0., shift=10., order=2):
    if conv > 1:
        # convolve with 5x5 filter
        filt = np.ones((conv, conv),dtype=float)
        image.dose = signal.convolve2d(image.dose, filt, mode='same', boundary='fill', fillvalue=0)
        image.dose /= np.sum(filt)

    # find y-center
    indX = np.abs(image.xCoord - 0).argmin()
    y_align = alignY(copy.deepcopy(image.dose), copy.deepcopy(image.yCoord), indX, image)
    ## use this alignment further on. flip about this rather than around
    image.yCoord -= y_align
    indY = np.abs(image.yCoord - 0).argmin()

    if image.yCoord[indY] < 0:
        indY -= 1
    mod_dose = copy.deepcopy(image.dose[:indY,:])
    #mod_dose[indY:,:] = 0.
    rev_dose = np.vstack((mod_dose, mod_dose[::-1,:]))
    # fix new y_coord for rev_dose
    yCoord = np.hstack((image.yCoord[:indY], -image.yCoord[:indY][::-1]))
    # extract profiles at xPos for comparison
    indX = np.abs(image.xCoord - xPos).argmin()
    rev_prof = rev_dose[:,indX]
    image_prof = np.interp(yCoord[::-1], image.yCoord[::-1], image.dose[::-1,indX])[::-1]
    image.yCoord = yCoord
    y, x = find_cut(image, shift)
    ratio = np.divide(image_prof, rev_prof)
    x_data = copy.deepcopy(image.yCoord[indY:y[1]])
    ratio = np.where(ratio == np.inf, 1, ratio)
    #print x, np.squeeze(diff[indY:y[1], indX])
    #x_data = x_data[:,np.newaxis]
    y_data = np.squeeze(ratio[indY:y[1]]) - 1
    if order == 1:
        x_data = x_data[:,np.newaxis]
        model, residual, _, _ = np.linalg.lstsq(x_data, y_data)
        polynomial = np.poly1d([model[0], 1.])
    elif order == 2:
        x_data = np.transpose([x_data*x_data, x_data])
        model, residual, _, _ = np.linalg.lstsq(x_data, y_data)
        polynomial = np.poly1d([model[0], model[1], 1.])
    #print x, a
    r2 = 1 - np.abs(residual / (y_data.size * y_data.var()))
    return model


def alignY(dose, coord, indX, image):
    indY = np.abs(image.yCoord - 0).argmin()
    val = dose[indY, indX]
    y = []
    y.append(np.interp(.5*val, dose[indY:,indX][::-1], coord[indY:][::-1]))
    y.append(np.interp(.5*val, dose[:indY,indX], coord[:indY]))
    return np.mean(y)


def find_cut(image, shift=10.):
    lp_x = map(float, image.ds.ExposureSequence[0].BeamLimitingDeviceSequence[0].LeafJawPositions)
    lp_y = map(float, image.ds.ExposureSequence[0].BeamLimitingDeviceSequence[1].LeafJawPositions)
    lp_x -= np.divide(lp_x, np.abs(lp_x)) * shift
    lp_y -= np.divide(lp_y, np.abs(lp_y)) * shift
    #plt.plot(image.yCoord, dose_c[:,indX])
    #plt.subplot(1,3,1)
    #plt.imshow(image.dose[281:486,409:614], cmap='magma')
    #plt.axis('off')
    #plt.subplot(1,3,2)
    #plt.imshow(dose_c[281:486,409:614], cmap='magma')
    #plt.axis('off')
    #plt.subplot(1,3,3)
    y = []
    for item in lp_y:
        y.append(np.abs(image.yCoord - item).argmin())
    y.sort()
    x = []
    for item in lp_x:
        x.append(np.abs(image.xCoord - item).argmin())
    x.sort()
    return y, x


# Function chooser
func_arg = {"-main": main}
# run specifc function as called from command line
if __name__ == "__main__":
    if sys.argv[1] == "-main":
        func_arg[sys.argv[1]]()
