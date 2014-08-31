#######################################
# Routines for processing of LMI data
# SBC - Started 8 February 2014
#######################################

import pyraf
from pyraf import iraf
import astropy.io.fits as pyfits
import astropy.coordinates
from astropy.time import Time
import numpy as np
import os, shutil
from glob import glob
import matplotlib.pylab as plt
import scipy.optimize as spo
from matplotlib.lines import Line2D
import re

# Try iqpkg
try:
	import iqpkg
	import iqutils
except:
	print "No iqpkg, iqutils.  Can't run lmi_stats."
	
# Necessary packages
iraf.images()
iraf.immatch()
#iraf.imfilter()
iraf.noao()
iraf.imred()
iraf.ccdred()
iraf.stsdas()
iraf.hst_calib()
iraf.nicmos()
iraf.imutil()

yes=iraf.yes
no=iraf.no
INDEF=iraf.INDEF
globclob=yes
globver=yes

NON_DECIMAL = re.compile(r'[^\d]+') # re set to remove non-digits
LMIFILTS = ["U", "B", "V", "R", "I", "SDSS-G", "SDSS-R", "SDSS-I", "SDSS-Z", "SDSS-U"]
LMIPIXSCALE = 0.240
ASTROMCMD = "/Users/jicapone/GitHub/RATIR-GSFC/code/reduction/astrom/vlt_autoastrometry.py"
FDICT = {"SDSS-U": "UMAG", "SDSS-G": "GMAG", "SDSS-R": "RMAG",
		 "SDSS-I": "IMAG", "SDSS-Z": "ZMAG"}
SEEDICT = {"SDSS-U": "co", "SDSS-G": "go", "SDSS-R": "ro", 
		   "SDSS-I": "bo", "SDSS-Z": "ko"}
LMIGAIN = 2.95

##############################################################################

def preproc(image, fkey="FILTER", ppre="p", bkey="BIASSEC", tkey="TRIMSEC", clobber=globclob, verbose=globver):

	'''Update header keywords, subtract overscan, and trim image.'''
	
	# Create 'FILTER' keyword
	pyim = pyfits.open(image)
	f1 = pyim[0].header["FILTER1"]; f2 = pyim[0].header["FILTER2"]
	if f1 == "OPEN":
		pyim[0].header[fkey] = f2
	elif f2 == "OPEN":
		pyim[0].header[fkey] = f1
	else:
		pyim[0].header[fkey] = "%s-%s" % (f1, f2)
		
	# Grab bias and trim sections for ccdproc
	bsec = pyim[0].header[bkey]
	tsec = pyim[0].header[tkey]
	
	# Write out results
	pyim.writeto("%s%s" % (ppre, image), clobber=clobber)

	# Configure ccdproc
	ccdproc = iraf.ccdred.ccdproc
	ccdproc.ccdtype = ""
	ccdproc.noproc = no
	ccdproc.fixpix = no
	ccdproc.overscan = yes
	ccdproc.trim = yes
	ccdproc.zerocor = no
	ccdproc.darkcor = no
	ccdproc.flatcor = no
	ccdproc.illumcor = no
	ccdproc.fringecor = no
	ccdproc.readcor = no
	ccdproc.scancor = no
	ccdproc.readaxis = 'line'
	ccdproc.fixfile = ""
	ccdproc.biassec = bsec
	ccdproc.trimsec = tsec
	ccdproc.interactive = no
	ccdproc.function = "legendre"
	ccdproc.order = 1
	ccdproc.sample = "*"
	ccdproc.naverage = 1
	ccdproc.niterate = 1
	ccdproc.low_reject = 3.0
	ccdproc.high_reject = 3.0
	ccdproc.grow = 0
	ccdproc(images="%s%s" % (ppre,image), output="")
	
	return
	
#############################################################################

def lmi_cals(images, dobias=yes, dobpm=yes, doflats=yes, btype="BIAS", ftype="SKY FLAT", fkey="FILTER", ppre="p", bkey="BIASSEC", tkey="TRIMSEC", bfile="Bias.fits", bpmfile="BPM.pl", bpre="b", flatpre="Flat", clobber=globclob, verbose=globver):

	'''Process bias and twilight flats, creating relevant calibration files for
	   nightly processing of LMI data.'''

	blist = []; flist = []
	
	# Preprocess all the relevant images
	for im in images:
		hdr = pyfits.getheader(im)
		if hdr['OBSTYPE'] == btype:
			blist.append("%s%s" % (ppre, im))
			preproc(im, fkey=fkey, ppre=ppre, bkey=bkey, tkey=tkey, clobber=clobber,
					verbose=verbose)
		elif hdr['OBSTYPE'] == ftype:
			flist.append("%s%s" % (ppre, im))
			preproc(im, fkey=fkey, ppre=ppre, bkey=bkey, tkey=tkey, clobber=clobber,
					verbose=verbose)
					
	if dobias:
		
		# Setup zerocombine
		zerocombine = iraf.ccdred.zerocombine
		zerocombine.combine = 'median'
		zerocombine.reject = 'avsigclip'
		zerocombine.ccdtype = ''
		zerocombine.process = no
		zerocombine.delete = no
		zerocombine.clobber = no
		zerocombine.scale = 'none'
		zerocombine.statsec = '*'
		zerocombine.nlow = 0
		zerocombine.nhigh = 1
		zerocombine.nkeep = 1
		zerocombine.mclip = yes
		zerocombine.lsigma = 3.0
		zerocombine.hsigma = 3.0
		zerocombine.rdnoise = hdr['RDNOISE']
		zerocombine.gain = hdr['GAIN']
		zerocombine.snoise = 0
		zerocombine.pclip = -0.5
		zerocombine.blank = 0.0
		
		# Run zerocombine
		bstr = ",".join(blist)
		if os.path.exists(bfile):
			os.remove(bfile)
		zerocombine(input=bstr, output=bfile)

	if dobpm:
	
		# Setup imcombine
		imcombine = iraf.immatch.imcombine
		imcombine.sigmas = "bsigma.fits"
		imcombine.combine = "median"
		imcombine.reject = "none"
		imcombine.project = no
		imcombine.outtype = "real"
		imcombine.outlimits = ""
		imcombine.offsets = "none"
		imcombine.masktype = "none"
		imcombine.scale = "none"
		imcombine.zero = "none"
		imcombine.weight = "none"
		
		# Run imcombine
		bstr = ",".join(blist)
		if os.path.exists("junk.fits"):
			os.remove("junk.fits")
		if os.path.exists("bsigma.fits"):
			os.remove("bsigma.fits")
		imcombine(input=bstr, output="junk.fits")

		# Run ccdmask to create BPM file
		if os.path.exists(bpmfile):
			os.remove(bpmfile)      
		iraf.ccdred.ccdmask("bsigma.fits", bpmfile, ncmed=7, nlmed=7, ncsig=15,
							nlsig=15, lsigma=10, hsigma=10, ngood=1, linterp=1, 
							cinterp=1, eqinterp=1)

		os.remove("bsigma.fits")
		os.remove("junk.fits")
		
	if doflats:
	
		# Subtract bias images from all flats
		for im in flist:
			lmi_debias(im, bfile=bfile)
			
		# Loop over filters
		for filt in LMIFILTS:
				
			# ID files in appropriate filter
			flis = []
			for im in flist:
				hdr = pyfits.getheader(im)
				if hdr[fkey] == filt:
					flis.append("%s%s" % (bpre, im))
					
			if flis == []:
				continue
					
			# Set up flatcombine
			flatcombine = iraf.ccdred.flatcombine
			flatcombine.combine = 'median'
			flatcombine.reject = 'avsigclip'
			flatcombine.ccdtype = ''
			flatcombine.process = no
			flatcombine.scale = 'median'
			flatcombine.statsec = ''
			flatcombine.nlow = 1
			flatcombine.nhigh = 1
			flatcombine.nkeep = 1
			flatcombine.mclip = yes
			flatcombine.lsigma = 3.0
			flatcombine.hsigma = 3.0
			flatcombine.rdnoise = hdr["RDNOISE"]
			flatcombine.gain = hdr["GAIN"]
			flatcombine.snoise = 0.0
			flatcombine.pclip = -0.5
			flatcombine.blank = 1.0
			
			# Run flatcombine
			fstr = ",".join(flis)
			if os.path.exists("%s-%s.fits" % (flatpre, filt)):
				os.remove("%s-%s.fits" % (flatpre, filt))
			flatcombine(input=fstr, output="%s-%s.fits" % (flatpre, filt))
			
			# Normalize
			iraf.iterstat.nsigrej = 5.0
			iraf.iterstat.maxiter = 10
			iraf.iterstat.verbose = globver
			iraf.iterstat.lower = INDEF
			iraf.iterstat.upper = INDEF
			iraf.iterstat("%s-%s.fits" % (flatpre, filt),verbose=no, prin=no)
			iraf.imarith("%s-%s.fits" % (flatpre, filt), "/", 
						 iraf.iterstat.median, "%s-%s.fits" % (flatpre, filt))
			
	return
	
###########################################################################

def lmi_detrend(images, otype="OBJECT", ppre="p", bkey="BIASSEC", tkey="TRIMSEC", bpre="b", bfile="Bias.fits", fkey="FILTER", fpre="f", flatpre="Flat", skybkg="SKYBKG", skysub="SKYSUB", skysig="SKYSIG", fpfx="F", clobber=globclob, verbose=globver):

	'''Identify science frames, pre-process, subtract bias, divide by flat, 
	   and correct for non-linearity.'''
	   
	images.sort()
	
	# Loop through all images
	for im in images:
		#if os.path.isfile( fpfx + im ):
		#	continue

		hdr = pyfits.getheader(im)
		if hdr['OBSTYPE'] == otype:
		
			# Preprocess
			preproc(im, fkey=fkey, ppre=ppre, bkey=bkey, tkey=tkey, clobber=clobber,
					verbose=verbose)
		
			# Bias subtraction
			lmi_debias("%s%s" % (ppre, im), bpre=bpre, bfile=bfile)
			
			# Flat field
			lmi_flat("%s%s%s" % (bpre, ppre, im), fpre=fpre, fkey=fkey,
					 flatpre="Flat") 
					 
			# Linearity correction
			lmi_lincor("%s%s%s%s" % (fpre, bpre, ppre, im))
			
			# Write sky values to header
			iraf.iterstat("%s%s%s%s" % (fpre, bpre, ppre, im), verbose=no, prin=no)
			fimg = pyfits.open("%s%s%s%s" % (fpre, bpre, ppre, im), mode="update")
			fimg[0].header[skybkg] = iraf.iterstat.median
			fimg[0].header[skysig] = iraf.iterstat.sigma
			fimg[0].header[skysub] = 0
			fimg.flush()
			fimg.close()
			
			# Final processed image
			shutil.copy("%s%s%s%s" % (fpre, bpre, ppre, im), 
						"%s%s" % (fpfx, im))
			
	return
	
###########################################################################

def lmi_lincor(image):

	'''Apply (pre-calculated) linearity correction for LMI data.'''
	
	fimg = pyfits.open(image, mode='update')
	fimg[0].data = 1.00305 * fimg[0].data - 1.20752e-6 * np.power(fimg[0].data, 2) + 1.267053e-11 * np.power(fimg[0].data, 3)
	fimg.flush()
	fimg.close()
	return
	
###########################################################################

def lmi_astrom(images, wpre="w"):

	'''WCS fits for images.'''
	
	images.sort()
	
	if len(images) != 0:

		# Loop through all images
		for image in images:
			
			print image

			# First need to update a bunch of keywords
			fimg = pyfits.open(image, mode="update")
			fimg[0].header["PIXSCALE"] = LMIPIXSCALE
			fimg[0].header["PIXSCAL1"] = LMIPIXSCALE
			fimg[0].header["PIXSCAL2"] = LMIPIXSCALE
			fimg[0].header["CTYPE1"] = "RA---TAN"
			fimg[0].header["CTYPE2"] = "DEC--TAN"
			fimg[0].header["WCSDIM"] = 2
			fimg[0].header["WAT0_001"] = "system=image"
			fimg[0].header["WAT1_001"] = "wtype=tan axtype=ra"
			fimg[0].header["WAT2_001"] = "wtype=tan axtype=dec"
			fimg[0].header["LTM1_1"] = 1.0
			fimg[0].header["LTM2_2"] = 1.0
		
			nax1 = fimg[0].header["NAXIS1"]; nax2 = fimg[0].header["NAXIS2"]
			fimg[0].header["CRPIX1"] = nax1 / 2
			fimg[0].header["CRPIX2"] = nax2 / 2
		
			ra = fimg[0].header["RA"]; dec = fimg[0].header["DEC"]
			coo = astropy.coordinates.ICRS("%sh%sm%ss %sd%sm%ss" % (ra.split(":")[0],
										   ra.split(":")[1], ra.split(":")[2],
										   dec.split(":")[0], dec.split(":")[1],
										   dec.split(":")[2]))
			fimg[0].header["CRVAL1"] = coo.ra.deg
			fimg[0].header["CRVAL2"] = coo.dec.deg
		
			fimg[0].header["CD1_1"] = -LMIPIXSCALE / 3600.0
			fimg[0].header["CD1_2"] = 0.0
			fimg[0].header["CD2_1"] = 0.0
			fimg[0].header["CD2_2"] = LMIPIXSCALE / 3600.0
			del fimg[0].header["CTYPE1U"]
			del fimg[0].header["CRPIX1U"]
			del fimg[0].header["CRVAL1U"]
			del fimg[0].header["CD1_1U"]
			del fimg[0].header["CFINT1"]
			del fimg[0].header["CTYPE2U"]
			del fimg[0].header["CRPIX2U"]
			del fimg[0].header["CRVAL2U"]
			del fimg[0].header["CD2_2U"]
			del fimg[0].header["CD1_2U"]
			del fimg[0].header["CD2_1U"]
			del fimg[0].header["CFINT2"]
			
			if os.path.exists("%s%s" % (wpre, image)):
				os.remove("%s%s" % (wpre, image))
			
			#fimg.writeto("%s%s" % (wpre, image))
			fimg.flush()
			fimg.close()
			
			#os.system("python %s %s%s" % (ASTROMCMD, wpre, image))  

			# Setup configuration files
			o1 = open("daofind.param", "w")
			o1.write("NUMBER\nXWIN_IMAGE\nYWIN_IMAGE\nMAG_AUTO\nFLAGS\nA_IMAGE\nB_IMAGE\n")
			o1.write("ELONGATION\nFWHM_IMAGE\nXWIN_WORLD\nYWIN_WORLD\n")
			o1.write("ERRAWIN_IMAGE\nERRBWIN_IMAGE\nERRTHETAWIN_IMAGE\nERRAWIN_WORLD\n")
			o1.write("ERRBWIN_WORLD\nERRTHETAWIN_WORLD\nFLUX_AUTO\nFLUX_RADIUS\n")
			o1.write("FLUXERR_AUTO")
			o1.close()
			
			o2 = open("default.conv", "w")
			o2.write("CONV NORM\n# 5x5 convolution mask of a gaussian PSF with FWHM = 3.0 pixels.\n0.092163 0.221178 0.296069 0.221178 0.092163\n0.221178 0.530797 0.710525 0.530797 0.221178\n0.296069 0.710525 0.951108 0.710525 0.296069\n0.221178 0.530797 0.710525 0.530797 0.221178\n0.092163 0.221178 0.296069 0.221178 0.092163")
			o2.close()
			
			# Run Sextractor
			os.system("sex -CATALOG_NAME %s.cat -CATALOG_TYPE FITS_LDAC -PARAMETERS_NAME daofind.param -DETECT_THRESH 3.0 -ANALYSIS_THRESH 3.0 -GAIN_KEY GAIN -PIXEL_SCALE 0 %s" % (image[:-5], image))
			
		# Run scamp for alignment
		imlist = " ".join(images)
		os.system("scamp -ASTREF_CATALOG 2MASS -CHECKPLOT_DEV NULL -SOLVE_PHOTOM N %s" % imlist.replace(".fits", ".cat"))
		
		# Update header info
		for image in images:
			os.system("missfits %s" % image)
			os.remove("%s.back" % image)
			os.remove("%s.head" % image[:-5])
	
	return
   
	
###########################################################################

def lmi_flat(image, fpre="f", fkey="FILTER", flatpre="Flat"):

	'''Flat-field an LMI image.'''
	
	# First need to identify correct flat
	hdr = pyfits.getheader(image)
	if not os.path.exists("%s-%s.fits" % (flatpre, hdr[fkey])):
		print "Error: No Flat for specified filter: %s" % hdr[fkey]
		return
		
	ccdproc = iraf.ccdred.ccdproc
	ccdproc.ccdtype = ""
	ccdproc.noproc = no
	ccdproc.fixpix = no
	ccdproc.overscan = no
	ccdproc.trim = no
	ccdproc.zerocor = no
	ccdproc.darkcor = no
	ccdproc.flatcor = yes
	ccdproc.illumcor = no
	ccdproc.fringecor = no
	ccdproc.readcor = no
	ccdproc.scancor = no
	ccdproc.readaxis = 'line'
	ccdproc.fixfile = ""
	ccdproc.flat = "%s-%s.fits" % (flatpre, hdr[fkey])
	ccdproc.interactive = no
	
	if os.path.exists("%s%s" % (fpre, image)):
		os.remove("%s%s" % (fpre, image))
	ccdproc(images=image, output="%s%s" % (fpre, image))
	
	return
			
###########################################################################

def lmi_debias(image, bpre="b", bfile="Bias.fits", clobber=globclob):

	'''Subtract bias from LMI image.'''
	
	ccdproc = iraf.ccdred.ccdproc
	ccdproc.ccdtype = ""
	ccdproc.noproc = no
	ccdproc.fixpix = no
	ccdproc.overscan = no
	ccdproc.trim = no
	ccdproc.zerocor = yes
	ccdproc.darkcor = no
	ccdproc.flatcor = no
	ccdproc.illumcor = no
	ccdproc.fringecor = no
	ccdproc.readcor = no
	ccdproc.scancor = no
	ccdproc.readaxis = 'line'
	ccdproc.fixfile = ""
	ccdproc.zero = bfile
	ccdproc.interactive = no
	ccdproc.function = "legendre"
	ccdproc.order = 1
	ccdproc.sample = "*"
	ccdproc.naverage = 1
	ccdproc.niterate = 1
	ccdproc.low_reject = 3.0
	ccdproc.high_reject = 3.0
	ccdproc.grow = 0
	
	if os.path.exists("%s%s" % (bpre, image)):
		os.remove("%s%s" % (bpre, image))
	ccdproc(images=image, output="%s%s" % (bpre, image))
	
	return
	
###########################################################################

def lmi_coadd(object, filt, align=no):

	'''Coadd all images of a given object in a specified filter.'''
	
	# Find the appropriate images
	ims = []
	allims = glob("Flmi.????.fits")
	for im in allims:
		h = pyfits.open(im)
		if (h[0].header["OBJECT"].replace(' ','_') == object) and (h[0].header["FILTER"] == filt):
			ims.append(im)    
	
	if len(ims) == 0:
		print "No images to coadd; Exiting!"
		return
		
	imlist = " ".join(ims)
	
	if align:
	
		# Setup configuration files
		o1 = open("daofind.param", "w")
		o1.write("NUMBER\nXWIN_IMAGE\nYWIN_IMAGE\nMAG_AUTO\nFLAGS\nA_IMAGE\nB_IMAGE\n")
		o1.write("ELONGATION\nFWHM_IMAGE\nCLASS_STAR\nXWIN_WORLD\nYWIN_WORLD\n")
		o1.write("ERRAWIN_IMAGE\nERRBWIN_IMAGE\nERRTHETAWIN_IMAGE\nERRAWIN_WORLD\n")
		o1.write("ERRBWIN_WORLD\nERRTHETAWIN_WORLD\nFLUX_AUTO\nFLUX_RADIUS\n")
		o1.write("FLUX_AUTO\nFLUXERR_AUTO")
		o1.close()
		
		o2 = open("default.conv", "w")
		o2.write("CONV NORM\n# 5x5 convolution mask of a gaussian PSF with FWHM = 3.0 pixels.\n0.092163 0.221178 0.296069 0.221178 0.092163\n0.221178 0.530797 0.710525 0.530797 0.221178\n0.296069 0.710525 0.951108 0.710525 0.296069\n0.221178 0.530797 0.710525 0.530797 0.221178\n0.092163 0.221178 0.296069 0.221178 0.092163")
		o2.close()
		
		for im in ims:
			os.system("sex -CATALOG_NAME %s.cat -CATALOG_TYPE FITS_LDAC -PARAMETERS_NAME daofind.param -DETECT_THRESH 3.0 -ANALYSIS_THRESH 3.0 -GAIN_KEY GAIN -PIXEL_SCALE 0 %s" % (im[:-5], im))
			
		# Run scamp for alignment
		os.system("scamp -CHECKPLOT_DEV NULL %s" % imlist.replace(".fits", ".cat"))
	
	# Run swarp for coaddition
	os.system("swarp -IMAGEOUT_NAME %s-%s.fits %s" % (object, filt, imlist))
	
	# Remove junk files
	os.system("rm *.head")
	
	return
	
##########################################################################

def lmi_stats(imlist, outf):

	'''Basic image statistics.'''
	
	ims = glob(imlist)
	ims.sort()
	outfile = open(outf, "w")
	
	for im in ims:
	
		hdr = pyfits.open(im)
		ra = hdr[0].header["RA"]; dec = hdr[0].header["DEC"]
		obj = hdr[0].header["OBJECT"].replace(' ','_'); filt = hdr[0].header["FILTER"]
		exptime = hdr[0].header["EXPTIME"]; tobs = hdr[0].header["DATE-OBS"]
		skybkg = hdr[0].header["SKYBKG"]; skysig = hdr[0].header["SKYSIG"]
		iqpkg.iqobjs(im, 3.0, 50000.0, wtimage="", verbose=no)
		os.system("getsdss.pl -r 10.0 -f %s.reg -p %s %s %s.txt" % (obj, ra, dec, obj))
		stars = iqutils.Starlist("%s.stars" % im)
		refstars = iqutils.Starlist("%s.reg" % obj)
		if len(refstars)==0:
			continue
		refstars.wcs2pix(im)
		if not FDICT.has_key(filt):
			continue
		refstars.set_mag(FDICT[filt])
		a,b=stars.match(refstars,maxnum=1000)
		if len(a)==0:
			continue
		fwhm = np.median(a.fwhms()) * LMIPIXSCALE
		zp, zpu = stars.zeropt(refstars,method="mean",rejout=1)
		if (zp == 0.0) or (zpu > 0.20):
			tzp = 99.0
			lmag1 = 99.0
			lmag2 = 99.0
		else:
			tzp = 25.0 + zp - 2.5 * np.log10(float(exptime))
			area = np.pi * np.power(1.2 * fwhm / LMIPIXSCALE, 2)
			lmag1 = -2.5 * np.log10(3 * np.sqrt(area) * skysig / exptime) + tzp
			lmag2 = -2.5 * np.log10(3 * np.sqrt(area) * skysig / np.sqrt(100.0 * exptime)) + tzp 

		outfile.write("%s%25s%15s%8s%8.2f%10.2f%10.2f%8.2f%10.3f%10.3f%10.3f%10.3f\n" % (im, tobs, obj, filt, exptime, skybkg, skysig,fwhm, tzp, zpu, lmag1, lmag2))
		print "Writing: "
		print "%s%25s%15s%8s%8.2f%10.2f%10.2f%8.2f%10.3f%10.3f%10.3f%10.3f\n" % (im, tobs, obj, filt, exptime, skybkg, skysig, fwhm, tzp, zpu, lmag1, lmag2)
	
	outfile.close()
	return
	
##########################################################################

def lmi_defringe(ims, filts=["SDSS-Z"], fkey="FILTER", skybkg="SKYBKG", skysub="SKYSUB", fpfx="F"):

	'''Create and apply fringe frame.'''
	
	ims.sort()
	
	for filt in filts:
		
		# Identify all the images in a given filter
		fims = []
		for im in ims:
			hdr = pyfits.open(im)
			if hdr[0].header[fkey] == filt:
				fims.append(im)

				if len(fims)==0:
					continue
		
		fimlist = ""
		for fim in fims:
		
			# Run iqobjs to get mask
			iqpkg.iqobjs(fim, 2.0, 50000.0, wtimage="", verbose=no)
			
			fimlist += "%s," % fim
		
		if fimlist == "":
			continue

		iqpkg.iqfringe(fimlist[:-1], "Fringe-%s.fits" % filt, verbose=no)
		
		iqpkg.iqdefringe(fimlist[:-1], "Fringe-%s.fits" % filt, outpfx="f",
						 verbose=no)
		
		for fim in fims:
			shutil.move("f%s" % fim, "%s%s" % (fpfx, fim[-13:]))
			
##########################################################################

# set clobber to True to rerun already processed frames
def lmi_fullnight(imlist, outf, do_all=True, ftype='SKY FLAT', otype="OBJECT", clobber=False):

	# create array of file names
	ims = glob(imlist)

	# create calibration files if they don't exist
	if (not os.path.exists('Bias.fits')) or clobber:
		print '\n\t* * * {} * * *'.format('lmi_cals - bias')
		lmi_cals(ims, dobpm=False, doflats=False)
	if (len(glob('Flat*.fits'))==0) or clobber:
		print '\n\t* * * {} * * *'.format('lmi_cals - flats')
		lmi_cals(ims, dobpm=False, dobias=False, ftype=ftype)
	
	# look for science frames
	oims = [] # files to-be-processed
	if not clobber:
		rims = glob('fbp'+imlist)
		for i in range(len(rims)): rims[i] = rims[i][3:]
		for im in ims:
			if im not in rims:
				hdr = pyfits.getheader(im)
				if hdr['OBSTYPE'] == otype:
					oims.append(im) # add unprocessed frame to list
	else:
		for im in ims:
			hdr = pyfits.getheader(im)
			if hdr['OBSTYPE'] == otype:
				oims.append(im) # add frame to list

	if len(oims):
		
		# calibrate science frames
		print '\n\t* * * {} * * *'.format('lmi_detrend')
		lmi_detrend(oims)
		
		# defringe science frames
		print '\n\t* * * {} * * *'.format('lmi_defringe')
		tims = []
		for i in range(len(oims)): tims.append('fbp'+oims[i])
		lmi_defringe(tims)
		
		# set corrected WCS keywords
		print '\n\t* * * {} * * *'.format('lmi_astrom')
		tims = []
		for i in range(len(oims)): oims[i] = 'F'+oims[i]
		lmi_astrom(tims)
	
	if do_all:

		# get nightly statistics
		print '\n\t* * * {} * * *'.format('lmi_stats')
		lmi_stats("Flmi.????.fits", outf)
		
		# plot stats get imdict for coaddition
		print '\n\t* * * {} * * *'.format('lmi_plots')
		imdict = lmi_plots("Flmi.????.fits", outf)
		
		# check for default.swarp.  this file is required to propogate keywords to coadded files.
		while not os.path.exists("./default.swarp"):
			raw_input("\nWarning: default.swarp has not been copied to the working directory.  Hit return once the file has been copied.")
		
		# coadd images
		print '\n\t* * * {} * * *'.format('lmi_coadd')
		for obj in imdict:
			for filt in imdict[obj]:
				lmi_coadd(obj, filt, align=yes)
		plt.show()

	return
	
##########################################################################

def lmi_plots(imlist, dfile):

	ims = glob(imlist)
	imdict = {}
	
	for im in ims:
		fimg = pyfits.open(im)
		obj = fimg[0].header["OBJECT"].replace(' ','_')
		filt = fimg[0].header["FILTER"]
		exptime = fimg[0].header["EXPTIME"]

		if not imdict.has_key(obj):
			imdict[obj] = {filt: float(exptime)}
		elif not imdict[obj].has_key(filt):
			imdict[obj][filt] = float(exptime)
		else:
			imdict[obj][filt] += float(exptime)
			
	for ob in imdict:
		print "%s:" % ob
		for fil in imdict[ob]:
			print "\t%s: %.2f s" % (fil, imdict[ob][fil])
	print "\n\n"
	
	data = open(dfile)
	lines = data.readlines()
	
	# Create figures
	f1 = plt.figure(1)  # Seeing
	ax1 = f1.add_subplot(1, 1, 1)
	f2 = plt.figure(2)  # Sky background
	ax2 = f2.add_subplot(1, 1, 1)
	f3 = plt.figure(3)  # Zeropoint
	ax3 = f3.add_subplot(1, 1, 1)
	
	for line in lines:
	
		temp = line.strip().split()
		t = Time(temp[1].replace("T", " "), format="iso", scale="utc")
	
		# Seeing
		ax1.plot(t.mjd, float(temp[7]), SEEDICT[temp[3]])
		
		# Sky background
		ax2.plot(t.mjd, float(temp[5]) / float(temp[4]), SEEDICT[temp[3]])
		
		# Limiting magnitude
		ax3.plot(t.mjd, float(temp[11]), SEEDICT[temp[3]])
	
	ax1.set_title("Nightly Seeing") 
	ax1.set_xlabel("MJD")
	ax1.set_ylabel("Seeing (arcseconds)")
	c1 = Line2D(range(1), range(1), color="white", marker='o', mfc="cyan")
	c2 = Line2D(range(1), range(1), color="white", marker='o', mfc="green")
	c3 = Line2D(range(1), range(1), color="white", marker='o', mfc="red")
	c4 = Line2D(range(1), range(1), color="white", marker='o', mfc="blue")
	c5 = Line2D(range(1), range(1), color="white", marker='o', mfc="black")
	ax1.legend([c1, c2, c3, c4, c5], ["u'", "g'", "r'", "i'", "z'"], numpoints=1)

	
	ax2.set_title("Nightly Sky Background")
	ax2.set_xlabel("MJD")
	ax2.set_ylabel(r"Sky Background (pixel$^{-1}$ s$^{-1}$)")
	ax2.legend([c1, c2, c3, c4, c5], ["u'", "g'", "r'", "i'", "z'"], numpoints=1)

	ax3.set_title("Nightly Limiting Magnitude (per 100 s exposure)")
	ax3.set_xlabel("MJD")
	ax3.set_ylabel("Limiting Magnitude (100 s exposure)")
	ax3.legend([c1, c2, c3, c4, c5], ["u'", "g'", "r'", "i'", "z'"], numpoints=1)
	ax3.set_ylim(25.0, 20.0)

	return imdict
	
		
		

		
	
			
					  
