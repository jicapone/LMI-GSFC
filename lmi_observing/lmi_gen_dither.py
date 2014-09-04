import numpy as np

"""
fnout		= file name out
objname		= string name
exptime 	= exposure times (same shape as filts)
filts		= array of string filter names
nexp		= number of exposures per filter
max_dith	= maximum ra or dec offset from center in arcsec
"""
def gen_dither_file( fnout, objname, exptime, filts, nexp, max_dith=20 ):
	fout = open(fnout, 'w')
	if type(exptime) is not list:
		exptime = [exptime]
	if type(filts) is not list:
		filts = [filts]
	if len(exptime) != len(filts):
		print "Error: exptime and filts must have same length."
		return 1
	nfilts = len(filts)
	dithers = np.random.rand( nexp, 2 ) * 2 * max_dith - max_dith
	for i in range(nexp):
		for j in range(nfilts):
			fout.write('"{}" {} {:.1f} {:.1f} {}\n'.format(objname, exptime[j], dithers[i,0], dithers[i,1], filts[j]))
	fout.close()
	return 0