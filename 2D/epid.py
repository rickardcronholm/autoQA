import dicom
import numpy as np



class image:

    def __init__(self, filename, MU=1):
        # read DICOM file
        self.ds = dicom.read_file(filename)
        self.xCoord, self.yCoord, self.dx, self.dy, self.xDir, self.yDir = self.getCoords()
        self.dose = self.getDose()
        if getattr(self.ds, 'RescaleType').lower() == 'gy/mu':
            self.dose *= MU

    def getDose(self):
        return self.ds.pixel_array * float(self.ds.RescaleSlope)

    def getCoords(self):
        # get directiona cosines
        if hasattr(self.ds, 'RTImageOrientation'):
            # verify this
            self.cosine = map(int, self.ds.RTImageOrientation)
        else:
            self.cosine = [1, 0, 0, 0, -1, 0]
        xDir = self.cosine[:3]
        yDir = self.cosine[3:]

        # get xCoord
        if abs(xDir[0]) == 1:
            dx = float(self.ds.ImagePlanePixelSpacing[0])
            xCoord = self.getCoord(sum(xDir), float(self.ds.RTImagePosition[0]),
            dx, self.ds.Columns)
            try:
			    xCoord += float(XRayImageReceptorTranslation[0])
            except NameError:
			    pass  # Portal Prediction does not have this tag
        else:
            dx = float(self.ds.ImagePlanePixelSpacing[1])
            xCoord = self.getCoord(sum(xDir), float(self.ds.RTImagePosition[1]),
            dx, self.ds.Rows)
            try:
			    xCoord += float(XRayImageReceptorTranslation[1])
            except NameError:
			    pass  # Portal Prediction does not have this tag
        # get yCoord
        if abs(yDir[1]) == 1:
            dy = float(self.ds.ImagePlanePixelSpacing[1])
            yCoord = self.getCoord(sum(yDir), float(self.ds.RTImagePosition[1]),
            dy, self.ds.Rows)
            try:
			    yCoord += float(XRayImageReceptorTranslation[1])
            except NameError:
		        pass  # Portal Prediction does not have this tag
        else:
            dy = float(self.ds.ImagePlanePixelSpacing[0])
            yCoord = self.getCoord(sum(yDir), float(self.ds.RTImagePosition[0]),
            dy, self.ds.Columns)
            try:
			    yCoord += float(XRayImageReceptorTranslation[0])
            except NameError:
			    pass  # Portal Prediction does not have this tag
        return xCoord, yCoord, dx, dy, xDir, yDir

    def getCoord(self, dircos, position, spacing, elements):
        coord = np.linspace(position, position + (elements-1) * dircos * spacing, num = elements)
        return coord

