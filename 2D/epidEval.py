import dicom
import numpy as np
from scipy import ndimage
import epid
import ast
import image_registration
import mynpgamma
import os
import psycopg2
import copy
import matplotlib.pyplot as plt
import struct
import pymssql


class evaluation:
	
    def __init__(self, filename, inputDir, apply_align=True):
        # read instruct_file
        with open(os.path.join(filename), 'r') as f:
            instructs = f.readlines()
        instructs = map(str.strip, instructs)
        # create epid image objects
        evl = epid.image(os.path.join(inputDir, '.'.join(['RI', instructs[0], 'dcm'])))
        ref = epid.image(os.path.join(inputDir, '.'.join(['RI', instructs[1], 'dcm'])),
            MU=float(instructs[2]))
        self.date = evl.ds.ContentDate
        self.time = evl.ds.ContentTime
        self.timestamp = ''.join([self.date, self.time.replace('.','')[:9]])
        self.MachineName = evl.ds.RadiationMachineName
        self.NominalEnergy = int(evl.ds.ExposureSequence[0].KVP)
        # self.FluenceMode = 'STANDARD'  # find a way to get this from ARIA DB
        self.FluenceMode = get_from_aria('FluenceMode', [evl.ds.ReferencedRTPlanSequence[0].ReferencedSOPInstanceUID, int(evl.ds.ReferencedBeamNumber)])
        if ref.ds.PatientID.startswith('QC_Patient'):
            # compute daily output correction factor
            output_size = 20.
            ref.output = comp_output(ref, output_size)
            evl.output = comp_output(evl, output_size)
            # print results to file
            #with open(filename, 'a') as f:
		    #    f.write('Daily correction: {0:.4f}\n'.format(ref.output / evl.output))
            self.CalibrationFactor = evl.output / ref.output
        elif ref.ds.PatientID.startswith('QC_epid_corr'):
			self.BackscatterCorrection = True
			self.rtplan = evl.ds.ReferencedRTPlanSequence[0].ReferencedSOPInstanceUID
        else:
            self.ID = evl.ds.PatientID
            self.rtplan = evl.ds.ReferencedRTPlanSequence[0].ReferencedSOPInstanceUID
            self.refbeamnr = int(evl.ds.ReferencedBeamNumber)
            # compute 2D gamma
            # get latest corresponding calibration
            self.CalibrationSer, calibrationFactor = get_calibration(self)
            
            if self.CalibrationSer is None:
				return
            
            # compute alignment
            # use orientation of coordinates to adjsut sign of alignment.
            self.align = alignment(evl, ref)
            # account for XRayImageReceptorTranslation
            try:
                XRIRT = map(float, 
                   evl.ds.XRayImageReceptorTranslation)
                self.align_matrix = (self.align[0] - XRIRT[abs(evl.yDir[1])],
                   self.align[1] - XRIRT[abs(evl.xDir[1])])
            except Exception:
                self.align_matrix = self.align
            if apply_align:
			    # apply aligment
                evl.xCoord += self.align[1]
                evl.yCoord += self.align[0]
            indYr = np.abs(ref.yCoord - 0).argmin()
            indYe = np.abs(evl.yCoord - 0).argmin()
            indXr = np.abs(ref.xCoord - 0).argmin()
            indXe = np.abs(evl.xCoord - 0).argmin()
            evl.dose /= calibrationFactor
            # get backscatter correction coefficients
            coeffs1, coeffs2 = get_correction_coefficients(self)
            # compute effective field length in EPID-arm direction
            y_length_along_arm = get_arm_length(evl.ds)
            # apply backscatter correction
            corr_array = compute_corr_array(evl.yCoord, y_length_along_arm, coeffs1, coeffs2)
            corr = np.diag(corr_array)
            #evl.dose = np.dot(corr, evl.dose)
            # plot for testing
            '''
            indYr = np.abs(ref.yCoord - 0).argmin()
            indYe = np.abs(evl.yCoord - 0).argmin()
            
            indXr = np.abs(ref.xCoord - 0).argmin()
            plt.plot(ref.yCoord, ref.dose[:, indXr])
            indXe = np.abs(evl.xCoord - 0).argmin()
            plt.plot(evl.yCoord, evl.dose[:, indXe])
            plt.plot(evl.yCoord, np.dot(corr, evl.dose)[:, indXe])
            plt.plot(evl.yCoord, corr_array)
            plt.gca().grid(True)
            print calibrationFactor
            
            plt.show()
            '''
            # extract instructions from dat file
            self.dose_threshold = float(instructs[3])
            self.dist_threshold = float(instructs[4])
            self.lower_dose_cutoff = float(instructs[5])
            self.local_dose = ast.literal_eval(instructs[6])
            self.pass_crit = float(instructs[7])
            max_reference = np.max(ref.dose)
            self.ldc = copy.deepcopy(self.lower_dose_cutoff)
            self.dosethres = copy.deepcopy(self.dose_threshold)
            self.lower_dose_cutoff *= max_reference
            if not self.local_dose:
                self.dose_threshold *= max_reference
		    # compute passing rate
            self.passing_rate, self.gamma = compute_gamma(self, ref, evl)
            self.passing_rate *= 100.
            self.passed = (self.passing_rate >= self.pass_crit)
            # construct binary data from gamma
            self.binaryGamma = gamma2binary(self.gamma, 2., 51)
            # print results to file
            #with open(filename, 'a') as f:
			#    f.write('passing rate: {0:.1f}\npassed: {1:s}\nAlignment: \n\tx: {2:.2f}\n\ty: {3:.2f}\n'.format(self.passing_rate, str(self.passed), self.align_matrix[1] / 10, self.align_matrix[0] / 10))
			

def gamma2binary(array, maxGamma, nBins):
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


def comp_output(obj, output_size):
    iso = map(float, obj.ds.IsocenterPosition)
    if sum(obj.xDir) > 0:
        x_start = np.where(obj.xCoord > iso[0] - output_size/2.)[0][0]
        x_stop = np.where(obj.xCoord < iso[0] + output_size/2.)[0][-1] + 1
    else:
        x_stop = np.where(obj.xCoord > iso[0] - output_size/2.)[0][-1] + 1
        x_start = np.where(obj.xCoord < iso[0] + output_size/2.)[0][0]
    if sum(obj.yDir) > 0:
        y_start = np.where(obj.yCoord > iso[1] - output_size/2.)[0][0]
        y_stop = np.where(obj.yCoord < iso[1] + output_size/2.)[0][-1] + 1
    else:
        y_stop = np.where(obj.yCoord > iso[1] - output_size/2.)[0][-1] + 1
        y_start = np.where(obj.yCoord < iso[1] + output_size/2.)[0][0]
    return np.mean(obj.dose[y_start:y_stop,x_start:x_stop])

		
def compute_gamma(eval_obj, ref, evl):
    # set some defaults
    dist_step_size = eval_obj.dist_threshold / 10.
    max_test_distance = eval_obj.dist_threshold * 2
    max_concurrent_calc_points = np.inf
    num_threads = 1
    # check if flipping is required
    coords_ref = (ref.yCoord, ref.xCoord)
    coords_evl = (evl.yCoord, evl.xCoord)
    ref.dose, coords_ref = flippit(ref.dose, coords_ref)
    evl.dose, coords_evl = flippit(evl.dose, coords_evl)
	# compute gamma
    gamma = mynpgamma.calc_gamma(
        coords_ref, ref.dose,
        coords_evl, evl.dose,
        eval_obj.dist_threshold, eval_obj.dose_threshold,
        lower_dose_cutoff=eval_obj.lower_dose_cutoff, 
        distance_step_size=dist_step_size,
        maximum_test_distance=max_test_distance,
        local_dose=eval_obj.local_dose,
        max_concurrent_calc_points=max_concurrent_calc_points,
        num_threads=num_threads)

    # compute passing rate   
    valid_gamma = gamma[~np.isnan(gamma)]
    return np.sum(valid_gamma <= 1) / float(len(valid_gamma)), valid_gamma


def get_from_aria(item, value):
    server = ""
    user = ""
    passwd = ""
    db = ""
    conn = pymssql.connect(server, user, passwd, db)
    cur = conn.cursor()
    if item == 'PatientSer':
        cur.execute("SELECT PatientSer FROM Patient WHERE PatientId='{:s}'".format(value))
        row = cur.fetchone()
        return row[0]
    elif item == 'FluenceMode':
		cur.execute("SELECT AddOn.AddOnId FROM PlanSetup INNER JOIN RTPlan ON RTPlan.PlanSetupSer=PlanSetup.PlanSetupSer INNER JOIN Radiation ON Radiation.PlanSetupSer=PlanSetup.PlanSetupSer INNER JOIN FieldAddOn ON FieldAddOn.RadiationSer=Radiation.RadiationSer INNER JOIN AddOn ON AddOn.AddOnSer=FieldAddOn.AddOnSer WHERE RTPlan.PlanUID='{:s}' AND Radiation.RadiationOrder={:d}".format(value[0], value[1]))
		row = cur.fetchone()
		try:
			return str(row[0]).replace('HML0299', 'STANDARD')
		except TypeError:
			return 'STANDARD'
        
        
def flippit(dose, coords):
	# listify coords
	coords = list(coords)
	flip = (flip_required(coords[0]), flip_required(coords[1]))
	for i in range(0, len(coords)):
		if flip[i]:
			#  coords = np.flip(coords[i], 0)
			#  dose = np.flip(dose, i)
			coords[i] = coords[i][::-1]
			if i == 0:
			    dose = dose[::-1, :]
			else:
			    dose = dose[:, ::-1]
	# revert list to tuples
	return dose, tuple(coords)


def pad(array, shape_diff):
    y_pad = (0, max(0, shape_diff[0]))
    x_pad = (0, max(0, shape_diff[1]))
    return np.lib.pad(array, (y_pad, x_pad),
        'constant', constant_values=(0, 0))


def imreg(ref, evl, return_error=False):
    # check shape consistency
    shape_diff = np.asarray(ref.shape) - np.asarray(evl.shape)
    q = pad(ref, -shape_diff)
    ref = pad(ref, -shape_diff)
    evl = pad(evl, shape_diff)
    dx,dy = image_registration.chi2_shifts.chi2_shift(ref, evl,
        return_error=False)
    return dx+float(shape_diff[1])/2, dy+float(shape_diff[0])/2

def alignment(ref, evl):
    evl_scale = 1
    filt = .25
    ref_dose = np.where(ref.dose > filt * np.max(ref.dose), ref.dose, 0.0)
    evl_dose = np.where(evl.dose > filt * np.max(evl.dose), evl.dose, 0.0)
    # cehck consistency of shape
    shape_divided = [float(ref.dose.shape[0])/float(evl.dose.shape[0]),
        float(ref.dose.shape[1])/float(evl.dose.shape[1])]
    if shape_divided == [1, 1]:
        dx, dy = image_registration.chi2_shifts.chi2_shift(ref_dose, evl_dose,
            return_error=False)
        
    elif ref.dose.size > evl.dose.size:
        evl_scale = int(np.round(np.sqrt(ref.dose.size / evl.dose.size)))
        evl_y = ndimage.interpolation.zoom(evl.yCoord, evl_scale, output=None)
        evl_x = ndimage.interpolation.zoom(evl.yCoord, evl_scale, output=None)
        dx, dy = imreg(ref_dose,
            ndimage.interpolation.zoom(evl_dose, evl_scale, output=None), 
                return_error=False)
    else:
        scale = int(np.round(np.sqrt(evl.dose.size / ref.dose.size)))
        ref_y = ndimage.interpolation.zoom(ref.yCoord, scale, output=None)
        ref_x = ndimage.interpolation.zoom(ref.xCoord, scale, output=None)
        dx, dy = imreg(ndimage.interpolation.zoom(ref_dose,
            scale, output=None), evl_dose, return_error=False)
   
    # convert suggested translation in pixel sizes to coordinates
    dx *= evl.dx / evl_scale
    dy *= evl.dy / evl_scale
    
    # account for flippability
    dy *= -2*(int(flip_required(evl.yCoord) == True) - .5)
    dx *= -2*(int(flip_required(evl.xCoord) == True) - .5)

    return (dy, dx)


def flip_required(x):
    dx = np.diff(x)
    return np.all(dx <= 0)


def interp(values, coords_orig, coords_new):
    # listify tuples
    coords_orig = list(coords_orig)
    coords_new = list(coords_new)
    # start by flipping if necessary
    flip = (flip_required(coords_orig[0]),
    flip_required(coords_orig[1]))
    for i in range(0, len(coords_orig)):
        if flip[i]:
            #  coords_orig[i] = np.flip(coords_orig[i], 0)
            #  values = np.flip(values, i)
            # also flip coords_new
            #  coords_new[i] = np.flip(coords_new[i], 0)
            coords_orig[i] = coords_orig[i][::-1]
            coords_new[i] = coords_new[i][::-1]
            if i == 0:
                values = values[::-1, :]
            else:
                values = values[:, ::-1]
    # revert lists to tuples
    coords_orig = tuple(coords_orig)
    coords_new = tuple(coords_new)
    rgi = RegularGridInterpolator(coords_orig, values, method='linear', bounds_error=False, fill_value = 0.0)
    YY, XX = mgrid[min(coords_new[0]):max(coords_new[0]):len(coords_new[0]) * 1j,
         min(coords_new[1]):max(coords_new[1]):len(coords_new[1]) * 1j]
    values_interp = rgi((YY, XX))
    # reflip if necessary
    for i in range(0, len(coords_orig)):
        if flip[i]:
            #  values_interp = np.flip(values_interp, i)
            if i == 0:
			    values_interp = values_interp[::-1, :]
            else:
                values_interp = values_interp[:, ::-1]
    return values_interp
    
    
def get_calibration(item):
	dbname = ''
	user = ''
	passwd = ''
	host = ''
	try:
		connString = "dbname={:s} user={:s} host={:s} password={:s}".format(dbname, user, host, passwd)
		conn = psycopg2.connect(connString)
		cur = conn.cursor()
		sqlStr =  "SELECT calibrationser, calibrationfactor FROM calibration WHERE treatmentmachine='{:s}' AND  nominalenergy={:d} AND fluencemode='{:s}' AND creationdate < to_timestamp('{:s}','YYYYMMDDHH24MISSMS') ORDER BY creationdate DESC LIMIT 1;".format(item.MachineName, item.NominalEnergy, item.FluenceMode, item.timestamp)
		cur.execute(sqlStr)
		row = cur.fetchone()
		cur.close()
		conn.close()
		return row[0], row[1]
	except Exception:
		return None, 1.0


def get_correction_coefficients(item):
	dbname = ''
	user = ''
	passwd = ''
	host = ''
	try:
		connString = "dbname={:s} user={:s} host={:s} password={:s}".format(dbname, user, host, passwd)
		conn = psycopg2.connect(connString)
		cur = conn.cursor()
		sqlStr = "SELECT coefficient12, coefficient11, coefficient10, coefficient22, coefficient21, coefficient20 FROM backscattercorrection WHERE treatmentmachine='{:s}' AND  nominalenergy={:d} AND fluencemode='{:s}' AND creationdate < to_timestamp('{:s}','YYYYMMDDHH24MISSMS') ORDER BY creationdate DESC LIMIT 1;".format(item.MachineName, item.NominalEnergy, item.FluenceMode, item.timestamp)
		cur.execute(sqlStr)
		row = cur.fetchone()
		cur.close()
		conn.close()
		return np.asarray([row[0], row[1], row[2]]), np.asarray([row[3], row[4], row[5]])
	except Exception:
		return np.asarray([0., 0., 0.]), np.asarray([0., 0., 0.])


def get_arm_length(ds):
	x = get_jaw_positions(ds, order=0)
	y = get_jaw_positions(ds, order=1, multiplier=[1, 1])
	theta = np.radians(float(ds.BeamLimitingDeviceAngle))
	# use trigonometry to compute the y-coordinate of all 4 corners
	points = []
	for ypos in y:
		points.append(ypos * np.cos(theta) - x[0] * np.sin(theta))
		points.append(ypos * np.cos(theta) + x[1] * np.sin(theta))
	# compute the mean of the two farthest point in the arm direction
	points.sort()
	return np.mean(points[2:])



def get_jaw_positions(ds, order=0, multiplier=[-1, 1]):
	return np.multiply(map(float, ds.ExposureSequence[0].BeamLimitingDeviceSequence[order].LeafJawPositions), multiplier)
	
	
def compute_corr_array(yCoord, y_length_along_arm, coeffs1, coeffs2):
	if np.abs(y_length_along_arm) > 50:
		c2 = np.polyval(coeffs2, y_length_along_arm)
	else:
		c2 = 0.
	c1 = np.polyval(coeffs1, y_length_along_arm)
	corr = np.polyval([c2, c1, 1.0], -yCoord)
	return np.where(yCoord > 0., 1., corr)
