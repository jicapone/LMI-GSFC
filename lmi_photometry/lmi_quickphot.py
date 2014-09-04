import numpy as np
import sys
import astropy.io.fits as pf
import photprocesslibrary as pplib
import os
from astropy import wcs
import get_SEDs
import pylab as pl

"""
Does quick photometry using sextractor and get_SEDs.py comparisons
Searches for object at GRB specified coordinates
"""
def quickphot(fitsfile, filter, sigma=1.5, grbra=None, grbdec=None, grbsep=1.0):

	if grbra == None or grbdec == None:
		print 'Please enter grb coordinates (in degrees)'
		print 'Format:'
		print 'quickphot(fitsfile, filter, sigma=3.0, grbra=104.14928, grbdec=41.78660, grbsep=1.0)'
		return -1
	
	file = pf.open(fitsfile)
	header = file[0].header
	
	#Run sextractor 
	cmd = 'sex ' + fitsfile + ' -c lmi.sex -ANALYSIS_THRESH ' + str(sigma) + ' -DETECT_THRESH ' + str(sigma)
	os.system(cmd)
	
	x,y,ra,dec,mag,magerr,e,fwhm,flags = np.loadtxt('temp.cat', unpack=True)
	
	#Parse out sextractor outputs and make some cuts (flag, mag, magerr and fwhm)
	keep 	   = (flags == 0) & (mag > 10) & (mag < 25) & (magerr < 0.4) & (fwhm < 9) 
	rakeep     = ra[keep]
	deckeep    = dec[keep]
	xkeep      = x[keep]
	ykeep	   = y[keep]
	magkeep    = mag[keep]
	magerrkeep = magerr[keep]
	
	#Calculate coordinates from pixel position and header (slightly different than ra and dec) and saves to '*.im'
	w = wcs.WCS(header)
	pixcrd = []

	for ind in range(len(xkeep)):
		pixcrd.append([xkeep[ind],ykeep[ind]])

	coord_out = w.wcs_pix2world(pixcrd,0)
	ra_out  = coord_out[:,0]
	dec_out = coord_out[:,1]
	
	imfile = 'coords.im'
	np.savetxt(imfile, np.transpose([ra_out, dec_out, magkeep]), fmt='%15.6f', header='RA\t DEC\t INST_MAG\t')
	
	#Run get_SEDs to get catalog/SED fit values for best match to the image coordinates
	catfile = 'coords.cat'
	os.system('export CDSCLIENT=http')
	os.system('python /Users/jicapone/GitHub/RATIR-GSFC/code/photometry/dependencies/get_SEDs.py '+imfile+ ' ' + filter + ' ' + catfile + ' 15 True False')
	#Unpacks magnitude values
	[refra,refdec,u_mag,g_mag,r_mag,i_mag,
		z_mag,y_mag,bigB_mag,bigV_mag,bigR_mag,
		bigI_mag,J_mag,H_mag,K_mag,u_err,g_err,
		r_err,i_err,z_err,y_err,bigB_err,bigV_err,
		bigR_err,bigI_err,J_err,H_err,K_err,mode] = np.loadtxt(catfile, unpack=True)
	
	#Can expand to include other filters
	magdict = {'u_mag': u_mag, 'g_mag': g_mag, 'r_mag': r_mag, 'i_mag': i_mag, 'z_mag': z_mag, 'y_mag': y_mag, 'J_mag': J_mag, 'H_mag': H_mag, 'K_mag': K_mag} 
	errdict = {'u_err': u_err, 'g_err': g_err, 'r_err': r_err, 'i_err': i_err, 'z_err': z_err, 'y_err': y_err, 'J_err': J_err, 'H_err': H_err, 'K_err': K_err} 	
	
	#Only use values that have mode ne -1 (can change if want to force SDSS, APASS, or USNOB)
	#and save the magnitude (of file's same filter)
	use1 = (mode != -1)
	use2 = (14 < magdict[filter+'_mag'])
	use = use1 & use2
	
	#Remove all bad fits
	refra      = refra[use]
	refdec     = refdec[use]
	refmag     = magdict[filter+'_mag'][use]
	referr     = errdict[filter+'_err'][use]
	magkeep    = magkeep[use]
	magerrkeep = magerrkeep[use]
	
	difarr = refmag-magkeep
	diferr = np.sqrt(referr**2 + magerrkeep**2)

	starnum = np.arange(len(difarr))
	zpfile = filter+'.zp'
	
	print difarr
	#Creates Absolute Magnitude file with coordinates
	#np.savetxt(zpfile, np.transpose([starnum, rakeep, deckeep, difarr, diferr, refmagkeep, referrkeep]), fmt='%15.6f', header='ID\t RA\t DEC\t DIFF\t DIFF_ERR\t CAT_MAG\t CAL_MAG_ERR\t')
	
	####FOR DIAGNOSTICS####
	pl.figure()
	pl.errorbar(refmag, difarr,xerr=referr, yerr=diferr, marker='o',linestyle='None')
		
	pl.xlabel('Catalog magnitude')
	pl.ylabel('Instrumental mag - Catalog mag')
	pl.savefig(filter+'_instzp.png')
	
	pl.clf()
	####END DIAGNOSTICS####
	
	#Do a 2 sigma iterative sigma clipping rejection (USNOB has a lot of catalog sources so
	#can reject because more than enough to get instrumental offset and outliers skewing data
	#particularly the overall error
	[difMean, difSig, difMedian, difMask, difSaveMask] = pplib.djs_iterstat(difarr, SigRej=2.0)
	
	#Removes clipped values from arrays, don't want to use their mean, use weighted mean
	#and standard deviation of the weighted mean
	newdifarr = difarr[np.where(difSaveMask == 1)]
	newdiferr = diferr[np.where(difSaveMask == 1)]
	
	wmean = sum(newdifarr/newdiferr**2)/sum(1.0/newdiferr**2)
	err   = np.sqrt(sum((newdifarr-wmean)**2/newdiferr**2)/sum(1.0/newdiferr**2))

	print wmean, err

	#Keeps all data from sextractor except non-detections (had to do good data clip above to get offset)
	good = mag < 99

	#Corrects magnitude and saves to txt file under 'absolutemag.txt'
	corrmag 	= mag + wmean
	corrmagerr 	= np.sqrt(magerr**2. + err**2.)
	
	#Creates Absolute Magnitude file with coordinates
	np.savetxt('absolutemag.txt', np.transpose([x[good], y[good], ra[good], dec[good], corrmag[good], corrmagerr[good]]), fmt='%15.6f', 
		header='X\t Y\t RA\t DEC\t CAL_MAG\t CAL_MAG_ERR\t')
	
	#Finds the closest match in position to the GRB coordinates to 1 arcsec
	#Sometimes requires 2 arcsec for missing objects
	grbmatch = pplib.nearest( ra[good]*np.cos(dec[good]*np.pi/180.),dec[good],grbra*np.cos(grbdec*np.pi/180.),grbdec, grbsep/3600.)
	
	#Updates the values to only include sources with actual magnitude detections
	potra 			= ra[good]
	potdec 			= dec[good]
	potmag 			= corrmag[good]
	potmagerr 		= magerr[good]
	potcorrmagerr 	= corrmagerr[good]
	potx 			= x[good]
	poty 			= y[good]
	
	#Returns the source(s) that matches location of GRB
	#Can copy and paste into 'photometry.txt'
	print potra[grbmatch], potdec[grbmatch], potmag[grbmatch], potmagerr[grbmatch], err
	print potx[grbmatch], poty[grbmatch]
	print fitsfile, filter, potmag[grbmatch][0], potmagerr[grbmatch][0], err
	
	return fitsfile, filter, potmag[grbmatch][0], potmagerr[grbmatch][0], err
	
	
#Can run quickphot on multiple files at once using a wildcard 
def completephot(wildcardsearch):
	fitsfiles = pplib.choosefiles(wildcardsearch)
	
	f = open('photometrymultiple.txt','w')	
	
	for file in fitsfiles:
		[fitsfile, filter, mag, magerr, correrr] = quickphot(file)
		f.write(fitsfile + '\t' + filter + '\t' + ut + '\t' + str(mag) + '\t' + str(magerr) + '\t' + str(correrr) + '\n')
	
	f.close()

#Creates lightcurve from text file ('photometry.txt')
def lightcurve(photfile, grbtime=None, filters=None, grbname=None, tprdct=120.0):

	if grbtime == None:
		print 'Please include GRB start time (s)'
		print 'Format:'
		print 'lightcurve(photfile, grbtime=11511.9, filters=["r","i"], grbname="GRB140606B", tprdct=48.0)'
		return -1
		
	if grbname == None:
		grbname='GRB'
	
	#fitsfile, filter, ut, mag, magerr, correrr = np.loadtxt(photfile,dtype=[str,str,str,str,str,str])
	asc = np.genfromtxt(photfile, dtype="str")
	
	filename= asc[:,0]	
	filter  = asc[:,1]
	ut 		= asc[:,2]
	mag 	= asc[:,3]
	magerr 	= asc[:,4]
	correrr = asc[:,5]
	
	subtime   = []
	magarr    = []
	magerrarr = []
	
	for i in range(len(ut)):
		(h,m,s) = ut[i].split(':')
		(fs, ps) = s.split('.')
 		result = int(h) * 3600 + int(m) * 60 + int(fs) + int(ps)*.01
 		
 		name = filename[i]
 		if name[0] == 'd':
 			stime = 86400 * int(name[1]) + result - grbtime
 		else:
 			stime = result-grbtime
 		subtime.append(stime)
 		
 		magarr.append(float(mag[i]))
 		magerrarr.append( np.sqrt(float(magerr[i])**2 + float(correrr[i])**2))
 		
 	subtime = np.array([subtime])[0]/60.0/60.0
 	magarr = np.array([magarr])[0]
 	magerrarr = np.array([magerrarr])[0]
 	flux = 10**(magarr/-2.5)*3.631e6
 	
 	colors = ['b','g','r','c','y']
 	cind = 0
 	linmag = []
 	linflx = []
 	fprdct = []
 	xtemp = np.arange(7000)/700.0-3
 	for filt in filters:
 		fildata = np.where(filter == filt)
 		pl.errorbar(np.log10(subtime[fildata]), magarr[fildata], yerr=magerrarr[fildata], marker='*', linestyle='None', color=colors[cind],label=filt)	
 		film = np.polyfit(np.log10(subtime[fildata]), magarr[fildata],1)
 		linmag.append(film)
 		filf = np.polyfit(np.log10(subtime[fildata]),np.log10(flux[fildata]), 1)
 		linflx.append(filf)
 		tflx = 10**(filf[0]*np.log10(tprdct)+filf[1])
 		fprdct.append(-2.5*np.log10(tflx/3.63e6))
 		pl.plot(xtemp, film[0]*xtemp+film[1], color=colors[cind])
 		cind = cind+1
 	pl.ylim(25,13)
 	pl.xlim(-2, 2.5)
 	
 	pl.legend()
 	pl.xlabel('Time since trigger [Log(hrs)]')
 	pl.ylabel('Magnitude')
 	pl.title(grbname)
 	
 	print '************************************************************'
 	print 'Linear fits of log(time) vs. magnitude ' +str(filters) + ':'
 	print linmag
 	
 	print '************************************************************'
 	print 'Linear fits of log(time) vs. log(flux) ' +str(filters) + ':'
 	print linflx
 	
 	print '************************************************************'
 	print 'Predictions for ' + str(tprdct) + ' hrs after burst ' +str(filters) + ':'
 	print fprdct
 	
 	pl.savefig('lightcurve.png')
 	pl.clf()	