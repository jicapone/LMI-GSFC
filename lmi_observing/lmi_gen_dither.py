import numpy as np

FILT_NAMES = ['SL-u','SL-g','SL-r','SL-i','SL-z','Yish']

# tile to increase observed area
def tile_2x2( fnout, objname, exptime, filts, nsets, max_dith, alt_filts=True ):
    
    fout = open(fnout, 'w')
    
    # write dither header
    fout.write("#title=true ra=false dec=false exposureTime=true numExposures=false filter=true muRA=false muDec=false epoch=false dRA=false dDec=false rotatorPA=false rotatorFrame=false xi=true eta=true comment=true commandOption=false\n#\n")
    
    # format and check args
    if type(exptime) is not list:
        exptime = [exptime]
    if type(filts) is not list:
        filts = [filts]
    if len(exptime) == 1:
        exptime = np.ones(np.shape(filts)) * exptime[0]
    if len(exptime) != len(filts):
        print "Error: exptime and filts must have same length."
        return 1
    nfilts = len(filts)

    # generate random dither offsets
    dithers = np.random.rand( 5*nsets, 2 ) * 2 * max_dith - max_dith
    dithers[-1] = np.array([0,0])

    """
    add in tiled field offsets:
        - field #0:  0',  0'
        - field #1: +6', +6'
        - field #2: +6', -6'
        - field #3: -6', +6'
        - field #4: -6', -6'
    """
    fields_off = np.array( [[ 6.,  6.],
                            [ 6., -6.],
                            [-6.,  6.],
                            [-6., -6.],
                            [ 0.,  0.]])
    fields_off *= 60

    i = 0
    for dither in dithers:
        dither += fields_off[i%5]
        i+=1

    # write to file
    maxexpl = len('{:.1f}'.format(int(np.ceil(np.max(exptime)))))
    if alt_filts:
        for i in range(5*nsets):
            for j in range(nfilts):
                comment = '"dXi,Eta = {:5.1f},{:5.1f} {}"'.format(dithers[i,0], dithers[i,1], filts[j])
                stemp = '"{}"    {:<' + '{}'.format(maxexpl) + '.1f}    {:>4}    {:>5.1f}    {:>5.1f}    {}\n'
                fout.write(stemp.format(objname, exptime[j], filts[j], dithers[i,0], dithers[i,1], comment))
    else:
        for j in range(nfilts):
            for i in range(5*nsets):
                comment = '"dXi,Eta = {:5.1f},{:5.1f} {}"'.format(dithers[i,0], dithers[i,1], filts[j])
                stemp = '"{}"    {:<' + '{}'.format(maxexpl) + '.1f}    {:>4}    {:>5.1f}    {:>5.1f}    {}\n'
                fout.write(stemp.format(objname, exptime[j], filts[j], dithers[i,0], dithers[i,1], comment))

    fout.close()

"""
fnout       = file name out
objname     = name of current object
exptime     = exposure times (same shape as filts)
filts       = array of string filter names
nexp        = number of exposures per filter
max_dith    = maximum ra or dec offset from center in arcsec
alt_filts   = alternate between selected filters
"""
def gen_dither_file( fnout, objname, exptime, filts, nexp, max_dith=20, alt_filts=False ):
    
    fout = open(fnout, 'w')
    
    # write dither header
    fout.write("#title=true ra=false dec=false exposureTime=true numExposures=false filter=true muRA=false muDec=false epoch=false dRA=false dDec=false rotatorPA=false rotatorFrame=false xi=true eta=true comment=true commandOption=false\n#\n")
    
    # format and check args
    if type(exptime) is not list:
        exptime = [exptime]
    if type(filts) is not list:
        filts = [filts]
    if len(exptime) == 1:
        exptime = np.ones(np.shape(filts)) * exptime[0]
    if len(exptime) != len(filts):
        print "Error: exptime and filts must have same length."
        return 1
    nfilts = len(filts)

    # generate random dither offsets
    dithers = np.random.rand( nexp, 2 ) * 2 * max_dith - max_dith
    dithers[-1] = np.array([0,0])

    # write to file
    maxexpl = len('{:.1f}'.format(int(np.ceil(np.max(exptime)))))
    if alt_filts:
        for i in range(nexp):
            for j in range(nfilts):
                comment = '"dXi,Eta = {:5.1f},{:5.1f} {}"'.format(dithers[i,0], dithers[i,1], filts[j])
                stemp = '"{}"    {:<' + '{}'.format(maxexpl) + '.1f}    {:>4}    {:>5.1f}    {:>5.1f}    {}\n'
                fout.write(stemp.format(objname, exptime[j], filts[j], dithers[i,0], dithers[i,1], comment))
    else:
        for j in range(nfilts):
            for i in range(nexp):
                comment = '"dXi,Eta = {:5.1f},{:5.1f} {}"'.format(dithers[i,0], dithers[i,1], filts[j])
                stemp = '"{}"    {:<' + '{}'.format(maxexpl) + '.1f}    {:>4}    {:>5.1f}    {:>5.1f}    {}\n'
                fout.write(stemp.format(objname, exptime[j], filts[j], dithers[i,0], dithers[i,1], comment))

    fout.close()