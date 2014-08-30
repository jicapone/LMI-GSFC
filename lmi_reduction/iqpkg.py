
####################################################################
# $Id: iqpkg.py,v 2.4 2008/12/01 01:27:44 cenko Exp $ #
####################################################################

import pyraf
from pyraf import iraf
import copy, os, shutil, glob, sys, string, re, math, operator
import pyfits
from types import *

# Access to the iqutils
from iqutils import *

# Necessary packages
iraf.images()
iraf.immatch()
iraf.imfilter()
iraf.noao()
iraf.imred()
iraf.ccdred()

yes=iraf.yes
no=iraf.no
INDEF=iraf.INDEF
hedit=iraf.hedit
imgets=iraf.imgets
imcombine=iraf.imcombine

pyrafdir="python/pyraf/"
pyrafdir_key='PYRAFPARS'

if os.environ.has_key(pyrafdir_key):
    pardir=os.environ[pyrafdir_key]
else:
    pardir=os.environ['HOME']+'/'+pyrafdir

if pardir[-1] != '/':
    pardir += '/'

globclob=yes
globver=yes

######################################################################

def iqmosaic(inlist, ccdproc, outpfx="p", biassec="!BIASSEC",
             trimsec="!TRIMSEC", joinkeys="", splitkeys="",
             clobber=globclob, verbose=globver):

    """ combines multiamp readout (mosaic) sections into a single
        image for subsequent processing """

    # Defaults / constants

    # Necessary package
    iraf.imutil()

    # Parse inputs
    infiles=iraffiles(inlist)
    joinks=joinkeys.split(',')
    splitks=splitkeys.split(',')

    # Set up for ccdproc
    ccdproc=iraf.ccdproc
    ccdproc.ccdtype=""
    ccdproc.noproc=no
    ccdproc.fixpix=no
    ccdproc.zerocor=no
    ccdproc.darkcor=no
    ccdproc.flatcor=no
    ccdproc.illumco=no
    ccdproc.fringec=no
    ccdproc.readcor=no
    ccdproc.scancor=no
    ccdproc.readaxis="line"

    if not ccdproc.overscan and not ccdproc.trim:
        print "No overscan-subtraction or data trimming will be done."

    # Parse user-specified bias section
    if ccdproc.overscan:
        biaskey=None
        re1=re.search("^\!(\w+)",biassec)
        if re1:
            biaskey=re1.group(1)
        else:
            usebias=biassec

    # Parse user-specified trim section
    if ccdproc.trim:
        trimkey=None
        re1=re.search("^\!(\w+)",trimsec)
        if re1:
            trimkey=re1.group(1)
        else:
            usetrim=trimsec

    # Process each file in turn
    for image in infiles:

        check_exist(image,"r")
        if check_head(image,"DEMOSAIC"):
            if int(get_head(image,"DEMOSAIC")):
                continue

        # Get the number of extensions
        fimg=pyfits.open(image)
        next=len(fimg)
        extfiles=[]
        
        # Process each data extension
        for i in range(1,next):

            # Extension header
            head=fimg[i].header
            
            # Find bias region
            if ccdproc.overscan:
                if biaskey:
                    if head.has_key(biaskey):
                        usebias=head[biaskey]
                    else:
                        print "Couldn't find biaskey %s in %s" % \
                              (biaskey,image)
                    ccdproc.biassec=usebias

            # Find trim region
            if ccdproc.trim:
                if trimkey:
                    if head.has_key(trimkey):
                        usetrim=head[trimkey]
                    else:
                        print "Couldn't find trimkey %s in %s" % \
                              (trimkey,image)
                    ccdproc.trimsec=usetrim

            # Close the file before mucking with it
            fimg.close()

            # Copy out the extension
            tmpimg1=iraf.mktemp("iqmsc")+".fits"
            iraf.imcopy("%s[%d]" % (image,i), tmpimg1)
            tmpimg2=iraf.mktemp("iqmsc")+".fits"

            # Run ccdproc, if necessary
            if ccdproc.overscan or ccdproc.trim:
                ccdproc(tmpimg1,output=tmpimg2,Stdout=1)
                iraf.imdel(tmpimg1,verify=no,go_ahead=yes)
            else:
                iraf.imrename(tmpimg1,tmpimg2)

            # Keep track of our various files
            extfiles.append(tmpimg2)

        # Check candidate output file name
        outfile=outpfx+image
        if len(outpfx)>0:
            check_exist(outfile,"w",clobber)
        elif not clobber:
            sys.stderr.write("Clobber is false and output is same as input\n")
            sys.exit(1)
        else:
            # Clobber original image
            iraf.imdel(image,verify=no,go_ahead=yes)

        # Generate updated keyword values
        xmin={}; xmax={}; ymin={}; ymax={}
        newks=[]; newval={}

        iext=1
        for tmpimg in extfiles:

            jkvals=get_head(tmpimg,joinks)
            for i in range(0,len(joinks)):
                kwd=joinks[i]
                re1=re.search("(\d+):(\d+),(\d+):(\d+)",jkvals[i])
                if re1:
                    if xmin.has_key(kwd):
                        xmin[kwd]=min(xmin[kwd],int(re1.group(1)))
                    else:
                        xmin[kwd]=int(re1.group(1))
                    if xmax.has_key(kwd):
                        xmax[kwd]=max(xmax[kwd],int(re1.group(2)))
                    else:
                        xmax[kwd]=int(re1.group(2))
                    if ymin.has_key(kwd):
                        ymin[kwd]=min(ymin[kwd],int(re1.group(3)))
                    else:
                        ymin[kwd]=int(re1.group(3))
                    if ymax.has_key(kwd):
                        ymax[kwd]=max(ymax[kwd],int(re1.group(4)))
                    else:
                        ymax[kwd]=int(re1.group(4))

            skvals=get_head(tmpimg,splitks)
            for i in range(0,len(splitks)):
                if len(splitks[i])>7:
                    newk=splitks[i][0:7]+("%d" % iext)
                else:
                    newk=splitks[i]+("%d" % iext)
                if len(skvals[i])>0:
                    newks.append(newk)
                    newval[newk]=skvals[i]

            iext+=1

        # List of input files
        joinstr=','.join(extfiles)
        
        # Join processed images together as one
        iraf.imjoin(joinstr,outfile,2,verbose=verbose)

        # Delete temporary files
        for tmpimg in extfiles:
            iraf.imdel(tmpimg,verify=no,go_ahead=yes)

        # Update join & split keywords
        for kwd in joinks:
            if xmin.has_key(kwd):
                jkval="[%d:%d,%d:%d]" % (xmin[kwd],xmax[kwd],
                                        ymin[kwd],ymax[kwd])
                update_head(outfile,kwd,jkval,"")

        for kwd in newks:
            update_head(outfile,kwd,newval[kwd],"")

        for kwd in splitks:
            delete_head(outfile,kwd)

        # Update FITS header
        update_head(outfile,"DEMOSAIC",1,"Has file been demosaicked?")
        update_head(outfile,"NAMPS",next-1,"Number of amplifiers combined")
        update_head(outfile,"NAMEMSC",image,"Name of preceding mosaic image")


######################################################################

def iqcals(inlist, ccdproc, dobias=yes, biasproc=no, biaskey="OBJECT",
           biasre="BIAS", biasroot="Bias", doflats=yes,
           flatproc=yes, flatkey="OBJECT", flatre="DOME",
           filtkey="FILTER", flatpre="Flat-", flatscale="mode",
           statsec="", normflat=yes,
           dobpm=no, bpmmethod="biassigma", bpmroot="BPM",
           bpmffilt="R,V,r,g,I,i", classkey="", classdef="-",
           mosaic=no, clobber=globclob, verbose=globver): 

    """ intelligent creation of calibration files """

    # Defaults/constants
    biaslis=[]
    flatlis={}
    bpmlis=[]
    biasbyclass={}
    flatbyclass={}
    useclass=len(classkey)>0

    global imcombine

    # Parse inputs
    infiles=iraffiles(inlist)

    # Load this only if we're doing mosaics
    if mosaic:
        iraf.mscred()
    
    ####################
                
    # Collect keywords
    for image in infiles:

        # Criteria for bias image
        biasval=get_head(image,biaskey)
        re_b=re.search(biasre,biasval,re.I)
        if re_b:
            if useclass and check_head(image,classkey):
                classval=get_head(image,classkey)
            else:
                classval=classdef
            if biasbyclass.has_key(classval):
                biasbyclass[classval].append(image)
            else:
                biasbyclass[classval]=[image]
        else:

            # Criteria for flatfield image
            flatval,filtval=get_head(image,[flatkey,filtkey])
            re_f=re.search(flatre,flatval,re.I)
            if re_f:
                if useclass and check_head(image,classkey):
                    classval=get_head(image,classkey)
                else:
                    classval=classdef
                if not flatbyclass.has_key(classval):
                    flatbyclass[classval]={}
                if flatbyclass[classval].has_key(filtval):
                    flatbyclass[classval][filtval].append(image)
                else:
                    flatbyclass[classval][filtval]=[image]

    # Can't make a bias image if there are no biases
    if len(biasbyclass)<1 and dobias:
        print "Found no bias images"
        print "I was looking for keyword '%s' matching '%s'" % \
              (biaskey,biasre)
        print "Among these files:  " + list2str(infiles,"%s")
        dobias=no

    # Use proper ccdred/mscred routines
    ccdmask=iraf.ccdred.ccdmask

    if mosaic:
        ccdproc=iraf.mscred.ccdproc
        darkcombine=iraf.mscred.darkcombine
        flatcombine=iraf.mscred.flatcombine
        imcombine=iraf.mscred.combine
    else:
        ccdproc=iraf.ccdred.ccdproc
        darkcombine=iraf.ccdred.darkcombine
        flatcombine=iraf.ccdred.flatcombine
        imcombine=iraf.immatch.imcombine

    ####################
                
    # Construct bias image
    if dobias:

        for classval in biasbyclass.keys():

            # Make a list-file
            if classval==classdef:
                fbroot=biasroot
            else:
                fbroot=biasroot+"-"+classval
            biasname=fbroot+".fits"
            fbname=fbroot+".lis"
            check_exist(fbname,"w",clobber)
            fb=open(fbname,"w")
            for image in biasbyclass[classval]:
                fb.write("%s\n" % image)
            fb.close()

            # Setup the darkcombine
            if biasproc:
                for image in biasbyclass[classval]:
                    if check_head(image,'CCDPROC'):
                        print "CCDPROC already run for image %s" % image
                        biasproc=no

            check_exist(biasname,"w",clobber)

            darkcombine.output=biasname
            darkcombine.combine="median"
            darkcombine.reject="none"
            darkcombine.ccdtype=""
            darkcombine.process=biasproc
            darkcombine.delete=no
            darkcombine.scale="none"

            ccdproc.zerocor=no
            ccdproc.darkcor=no
            ccdproc.flatcor=no
            ccdproc.illumco=no
            ccdproc.fringec=no

            # Run darkcombine
            iraf.darkcombine("@"+fbname,Stdout=1)

    ####################
                
    # Construct flatfield image(s)
    if doflats:

        for classval in flatbyclass.keys():

            flatlis=flatbyclass[classval]

            for filt in flatlis.keys():

                fflis=flatlis[filt]

                # Make a list-file
                if classval==classdef:
                    fbroot=biasroot
                    ffroot=flatpre+filt
                else:
                    fbroot=biasroot+"-"+classval
                    ffroot=flatpre+filt+"-"+classval

                biasname=fbroot+".fits"
                flatname=ffroot+".fits"
                ffname=ffroot+".lis"
                check_exist(ffname,"w",clobber)
                ff=open(ffname,"w")
                for image in fflis:
                    ff.write("%s\n" % image)
                ff.close()

                # Setup the flatcombine
                if flatproc:
                    for image in fflis:
                        if check_head(image,'CCDPROC'):
                            print "CCDPROC already run for image %s" % image
                            flatproc=no

                ccdproc.zerocor=yes
                ccdproc.zero=biasname

                check_exist(flatname,"w",clobber)

                flatcombine.output=flatname
                flatcombine.combine="average"
                flatcombine.reject="avsigclip"
                flatcombine.ccdtype=""
                flatcombine.process=flatproc
                flatcombine.subsets=no
                flatcombine.delete=no
                flatcombine.scale=flatscale
                flatcombine.statsec=statsec
                flatcombine.mclip=yes
                flatcombine.lsigma=3.0
                flatcombine.hsigma=3.0

                # Run flatcombine
                iraf.flatcombine("@"+ffname,Stdout=1)

                # Normalize flatfield if requested
                if normflat:
                    # Check if it has been normalized already
                    if check_head(flatname,'NORMFLAT'):
                        if int(get_head(flatname,'NORMFLAT'))>0:
                            print ("Flatfield %s has already been "+\
                                   "normalized") % flatname
                            continue
                    # Use iterstat to get median
                    iraf.iterstat(flatname+statsec,nsigrej=5,maxiter=10,
                                  prin=no,verbose=no)
                    medflat=float(iraf.iterstat.median)
                    # Check the median value before using it
                    if medflat<=0:
                        print "Median of flatfield %s < 0, skipping" % \
                              flatname
                        continue
                    if medflat<500:
                        print "Small median value for flatfield %s?" % flatname
                    if verbose:
                        print "Dividing flatfield %s by %.2f..." % \
                              (flatname,medflat)
                    # Normalize and update header keywords
                    iraf.imarith(flatname,"/",medflat,flatname,
                                 verbose=no,noact=no)
                    update_head(flatname,"NORMFLAT",1,
                                "Has flatfield been normalized?")
                    update_head(flatname,"MEDFLAT",medflat,
                                "Median value of original flatfield")

    ####################
                
    # Check that we have a proper setup for the BPM
    if dobpm:

        # Double-check BPM methods
        if bpmmethod=='flatratio':
            if len(flatlis)<2:
                dobpm=0
                print "Need flats in at least two filters "+\
                      "for BPM by flatratio"
        elif bpmmethod=='biassigma':
            # Have to make this check later
            pass
        else:
            print "Fail to recognize bpmmethod '%s'" % bpmmethod
            dobpm=0

    # Construct BPM
    if dobpm:

        if bpmmethod=="flatratio":

            for classval in flatbyclass.keys():

                flatlis=flatbyclass[classval]

                if classval==classdef:
                    bpmname=bpmroot+".pl"
                else:
                    bpmname=bpmroot+"-"+classval+".pl"

                # Check if we can make BPM name
                if os.path.exists(bpmname):
                    if clobber:
                        check_exist(bpmname,"w")
                    else:
                        print "BPM file '%s' exists and clobber is no" % \
                              bpmname
                        dobpm=0

                # Construct list of flats we have made
                bpmflat=[]
                for filt in flatlis.keys():
                    if classval==classdef:
                        ffroot=flatpre+filt
                    else:
                        ffroot=flatpre+filt+"-"+classval
                    bpmff=ffroot+".fits"
                    if os.path.exists(bpmff):
                        bpmflat.append(bpmff)

                # Choose "best" filters if these are specified
                bpmuse=[]
                if len(bpmffilt)>1:
                    bpmtopfs=bpmffilt.split(',')
                    for bpmtopf in bpmtopfs:
                        if classval==classdef:
                            ffroot=flatpre+bpmtopf
                        else:
                            ffroot=flatpre+bpmtopf+"-"+classval
                        bpmff=ffroot+".fits"
                        if bpmff in bpmflat:
                            bpmuse.append(bpmff)
                            # only need two filters
                            if len(bpmuse)==2:
                                break

                # Finalize choice of comparison filters
                if len(bpmuse)<2:
                    for bpmff in bpmflat:
                        if bpmff not in bpmuse:
                            bpmuse.append(bpmff)
                            # only need two filters
                            if len(bpmuse)==2:
                                break

                # Make the ratio image
                fratio=iraf.mktemp("iqfratio")+".fits"
                iraf.imarith(bpmuse[0],"/",bpmuse[1],fratio,
                             verbose=no,noact=no)

                # Calculate bad pixel mask from the ratio
                check_exist(bpmname,"w")
                iraf.ccdmask(fratio,bpmname,ncmed=15,nlmed=15,ncsig=31,nlsig=31,
                             lsigma=10,hsigma=10,ngood=1,linterp=1,cinterp=1,
                             eqinterp=1)

                # Delete the ratio file
                check_exist(fratio,"w",yes)

        elif bpmmethod=="biassigma":

            for classval in biasbyclass.keys():

                biaslis=biasbyclass[classval]

                if classval==classdef:
                    bpmname=bpmroot+".pl"
                    biasname=biasroot+".fits"
                else:
                    bpmname=bpmroot+"-"+classval+".pl"
                    biasname=biasroot+"-"+classval+".fits"

                if len(biaslis)<3:
                    dobpm=0
                    print "Need at least three raw bias images "+\
                          "for BPM by biassigma"
                    continue

                # Workspace
                bpmlis=biaslis
                bdummy=iraf.mktemp("iqbdummy")+".fits"
                bdummy2=iraf.mktemp("iqbdummy")+".fits"
                bsigma=iraf.mktemp("iqbsigma")+".fits"

                # Setup a dummy imcombine (interested in the sigma image)
                (fbroot,fbext)=os.path.splitext(bpmname)
                fbname=fbroot+".lis"
                check_exist(fbname,"w",yes)
                fb=open(fbname,"w")
                for image in bpmlis:
                    if not biasproc or check_head(image,'CCDPROC'):
                        fb.write("%s\n" % image)
                    else:
                        print "CCDPROC not run yet for candidate " + \
                              "BPM bias %s" % image
                fb.close()

                # Run imcombine
                iraf.imcombine("@"+fbname,bdummy,sigmas=bsigma,
                               combine="median",reject="none",project=no,
                               offsets="none",masktyp="none",scale="none",
                               zero="none",weight="none",Stdout=1)

                if os.path.exists(biasname):
                    # Divide by the actual bias
                    iraf.imarith(bdummy,'/',biasname,bdummy2)
                else:
                    # Divide by a single bias image
                    iraf.imarith(bdummy,'/',bpmlis[0],bdummy2)

                # Make BPM
                check_exist(bpmname,"w")
                iraf.ccdmask(bsigma,bpmname,ncmed=7,nlmed=7,ncsig=15,nlsig=15,
                             lsigma=6,hsigma=6,ngood=5,linterp=1,cinterp=1,
                             eqinterp=1)

                # Cleanup
                check_exist(bdummy,"w",yes)
                check_exist(bdummy2,"w",yes)
                check_exist(bsigma,"w",yes)

######################################################################

def iqflatten(inlist, flat, outpfx="f", normflat=yes, statsec="",
              subsky=yes, clipflat=yes, maskin="!BPM", maskout="-bpm",
              maxfac=2.0, vignflat=no, filtsize=15, cutoff=0.5,
              border=-100000., clobber=globclob, verbose=globver):

    """ intelligent flatfield division"""

    if (len(outpfx) == 0):
        print "Warning: this will overwrite old files"

    # Parse inputs
    infiles=iraffiles(inlist)
    flatimg=fitsfile(flat)
    fflat="_F"+flatimg
    check_exist(fflat,"w",yes)

    # Normalize the flatfield if requested
    if normflat:
        imgets(flatimg,"NORMFLAT",Stdout=1,Stderr=1)
        if (int(imgets.value)!=0):
            if verbose:
                print "Flatfield %s has already been normalized" % flatimg
        else:
            iraf.iterstat(flatimg+statsec,nsigrej=5,maxiter=10,
                          prin=no,verbose=no)
            medflat=float(iraf.iterstat.median)
            if medflat<=0:
                print "Median of flatfield %s < 0, skipping" % flatimg
            else:
                if medflat<500:
                    print "Small median value for flatfield %s?" % flatimg
                if verbose:
                    print "Dividing flatfield %s by %.2f..." % \
                          (flatimg,medflat)
                iraf.imarith(flatimg,"/",medflat,flatimg,verbose=no,noact=no)
                update_head(flatimg,"NORMFLAT",1,
                            "Has flatfield been normalized?")
                update_head(flatimg,"MEDFLAT",medflat,
                            "Median value of original flatfield")

    # Clip the flatfield if requested
    if clipflat:
        imgets(flatimg,"CLIPFLAT",Stdout=1,Stderr=1)
        if (int(imgets.value)!=0):
            if verbose:
                print "Flatfield %s has already been clipped" % flatimg
        else:
            if check_head(flatimg,'NORMFLAT') and \
                   get_head(flatimg,'NORMFLAT'):
                medflat=1.0
            else:
                iraf.iterstat(flatimg+statsec,nsigrej=5,maxiter=10,
                              prin=no,verbose=no)
                medflat=float(iraf.iterstat.median)

            iqclip(flatimg,lthresh=medflat/maxfac,
                   hthresh=medflat*maxfac,bookend=no,
                   replace=medflat,maskin=maskin,maskout=maskout,
                   maskval=1,clobber=clobber,verbose=verbose)
            
            update_head(flatimg,"CLIPFLAT",1,
                        "Has flatfield been clipped?")
            update_head(flatimg,"CLIPREPL",medflat,
                        "Replacement value for clipped pixels")

    # Median-filter flatfield for clean borders
    if vignflat:
        if (verbose):
            print "Median-filtering flatfield for vignetting clip..."

        s=iraf.fmedian(flatimg,fflat,xwindow=filtsize,ywindow=filtsize,
                       boundary="nearest",hmin=-32768,hmax=32767,
                       zmin=INDEF,zmax=INDEF,zloreject=INDEF,
                       zhireject=INDEF,unmap=yes,Stdout=1)

    # Flatfield the images
    for image in infiles:

        check_exist(image,"r")

        # check processing status
	imgets(image,"FLATDIV",Stdout=1,Stderr=1)
        if (int(imgets.value) != 0):
            print "Warning: %s has already been flatfielded!" % image
            check_exist(outpfx+image,"w",clobber)
            iraf.imcopy(image,outpfx+image)
            continue

        if verbose:
            print "Processing %s with flatfield %s" % (image,flatimg)

        # Calculate image sky value if necessary
        sky=""
        if check_head(image,"SKYBKG"):
            try:
                sky=float(get_head(image,"SKYBKG"))
            except:
                sky=""

        if sky=="":
            iraf.iterstat(image+statsec,nsigrej=5,maxiter=10,
                          prin=no,verbose=no)
            sky=float(iraf.iterstat.median)
            update_head(image,"SKYBKG","%.3f" % sky,"Value of SKY median")

        if subsky:
            skyval=sky
        else:
            skyval=0.0

        # Intelligent flatfielding
        if (len(outpfx) != 0):
            check_exist(outpfx+image,"w",clobber)

        if vignflat:
            # Vignetting method
            iraf.imexpr(expr="z < %.2f ? %.1f : x/y - %s" % \
                        (cutoff,border,skyval),
                        output=outpfx+image,refim="x",
                        x=image,y=flat,z=fflat,verbose=verbose)
        else:
            # Simple method
            iraf.imarith(image,"/",flat,outpfx+image,verbose=no,
                         noact=no)
            iraf.imarith(outpfx+image,"-",skyval,outpfx+image,
                         verbose=no,noact=no)
            
        # Update FITS header
        update_head(outpfx+image,"FLATDIV",1,
                    "Has flat-field division occurred?")
        update_head(outpfx+image,"FLATNAME",flat,
                    "Name of flatfield")
        update_head(outpfx+image,"SKYBKG","%.3f" % sky,
                    "Value of SKY background")
        if subsky:
            update_head(outpfx+image,"SKYSUB",1,
                        "Has sky-subtraction occurred?")
        else:
            update_head(outpfx+image,"SKYSUB",0,
                        "Has sky-subtraction occurred?")
        update_head(outpfx+image,"BORDER",border,
                    "Value for no-exposure border")
        update_head(outpfx+image,"NAME1",image,
                    "Previous name of file")

    # Cleanup
    check_exist(fflat,"w",yes)

    if (verbose):
        print "Flatfielding finished!"
        
######################################################################

def iqsubsky(inlist, sub=yes, statsec="", skykey='SKYBKG',
             subkey='SKYSUB', minlim=no, minval=-90000.,
             clobber=globclob, verbose=globver):

    """ marginally intelligent sky subtraction """

    # Parse inputs
    infiles=iraffiles(inlist)

    # Subtracting or not
    if sub:
        subsign='-'
    else:
        subsign='+'

    # Perform sky subtraction/addition
    for image in infiles:

        check_exist(image,"r")

        # Parse user-specified sky value / header keyword
        sky=""
        if check_head(image,skykey):
            try:
                sky=float(get_head(image,skykey))
            except:
                sky=""

        if sky=="":
            iraf.iterstat(image+statsec,nsigrej=5,maxiter=10,
                          prin=no,verbose=no)
            sky=float(iraf.iterstat.median)
            update_head(image,skykey,sky,"Value of SKY median")

        # Check status of image's sky-subtraction
        subdone=no
        if check_head(image,subkey):
            try:
                subdone=int(get_head(image,subkey))
            except:
                subdone=no

        # Do we have to take action or not?
        if sub and subdone:
            continue
        if not sub and not subdone:
            continue

        # Do the subtraction
        if minlim:
            # Complicated method
            tmpimg=iraf.mktemp("iqss")+".fits"
            iraf.imexpr(expr="x < %.3f ? x : x %s %.3f" % \
                        (minval,subsign,sky),
                        output=tmpimg,x=image,
                        verbose=verbose)
            iraf.imdel(image,verify=no,go_ahead=yes)
            iraf.imrename(tmpimg,image,verbose=no)
        else:
            # Simple method
            iraf.imarith(image,subsign,sky,image,verbose=no,
                         noact=no)
            
        # Update FITS header
        if sub:
            subval=1
        else:
            subval=0

        update_head(image,subkey,subval,"Has sky-subtraction occured?")

######################################################################

def iqborder(inlist, minlim=no, minval=-100000., bkpfx="_B",
             clobber=globclob, verbose=globver):

    """ backs up original image & zeroes border region """

    # Parse inputs
    infiles=iraffiles(inlist)

    # Big loop
    for image in infiles:

        check_exist(image,"r")
        (imgroot,imgext)=os.path.splitext(image)

        # Backup current image
        bkimage=bkpfx+image
        check_exist(bkimage,"w",yes)

        if verbose:
            print "Backing up %s to %s" % (image,bkimage)

        iraf.imcopy(image,bkimage)
            
        # Replace border with zeroes
        tmpimg=iraf.mktemp("iqbr")+".fits"
        iraf.imexpr(expr="x < %.3f ? 0.0 : x" % minval,
                    output=tmpimg,x=image,
                    verbose=verbose)
        iraf.imdel(image,verify=no,go_ahead=yes)
        iraf.imcopy(tmpimg,image)
        iraf.imdel(tmpimg,verify=no,go_ahead=yes)

        if verbose:
            print "Zeroed border regions of %s" % image

######################################################################

def iqmask(inlist,mask="!BPM",method='constant',value='0.0',
           clobber=globclob,verbose=globver):

    """ flexible image masking """

    # Necessary packages
    iraf.proto()
    
    # Parse inputs
    infiles=iraffiles(inlist)

    # Mask specification
    bpmkey=None
    re1=re.search("^\!(\w+)",mask)
    if re1:
        bpmkey=re1.group(1)
    elif os.path.exists(mask):
        bpmfile=mask
    else:
        print "Failed to parse mask specification '%s'" % mask
        return

    # Reference image specification
    refkey=None
    reffiles=None
    if method=='image':
        re1=re.search("^\!(\w+)",value)
        if re1:
            refkey=re1.group(1)
        else:
            reffiles=iraffiles(value)
            if len(reffiles)==0:
                print "Failed to parse reference image specification '%s'" % value
                return

    # Check for a good method specification
    if method != 'fixpix' and method != 'constant' and \
           method != 'image':
        print "Could not parse method specification %s" % method
        return

    # Process each file in turn
    for image in infiles:

        check_exist(image,"r")

        # Identify BPM for this image
        if bpmkey:
            usebpm=get_head(image,bpmkey)
            check_exist(usebpm,"r")
        else:
            usebpm=bpmfile

        if method=='fixpix':

            # Fixpix method
            iraf.fixpix(image,usebpm,linterp=1,cinterp=2,
                        verbose=no,pixels=no)

        elif method=='constant':

            # Constant-value method
            try:
                fval=float(value)
            except:
                print "Failed to parse value '%s' as floating-point" % value
                return
            
            useoldway=None
            if useoldway:
                tmpimg=iraf.mktemp("iqmsk")+".fits"
                iraf.imexpr(expr="y > 0 ? %.2f : x" % fval,
                            output=tmpimg,x=image,y=usebpm,refim="x",
                            dims="auto",intype="auto",outtype="auto",
                            verbose=verbose)
                iraf.imdel(image,verify=no,go_ahead=yes)
                iraf.imrename(tmpimg,image,verbose=verbose)
            else:
                tmpbpm=iraf.mktemp("iqbpm")+".fits"
                iraf.imcopy(usebpm,tmpbpm,verbose=no)
                apply_bpm(image,tmpbpm,fval,
                          clobber=clobber,verbose=verbose)
                iraf.imdel(tmpbpm,verify=no,go_ahead=yes)

        elif method=='image':

            # Replace from template image
            if refkey:
                refimg=get_head(image,refkey)
            elif len(reffiles)==len(infiles):
                refimg=reffiles[infiles.index(image)]
            else:
                refimg=reffiles[0]
            
            if not os.path.exists(refimg):
                print "Failed identify reference image '%s'" % refimg
                return
                
            tmpbpm=iraf.mktemp("iqbpm")+".fits"
            iraf.imcopy(usebpm,tmpbpm,verbose=no)
            apply_bpm(image,tmpbpm,refimg,
                      clobber=clobber,verbose=verbose)
            iraf.imdel(tmpbpm,verify=no,go_ahead=yes)

        # Update FITS header
        update_head(image,'IQMASK',1,"Masked with iqmask?")
        update_head(image,'MASKMETH',method,"Mask method")
        if method=='constant':
            update_head(image,'MASKCNST',fval,"Bad pixel replacement value")
        elif method=='image':
            update_head(image,'MASKIMG',refimg,"Bad pixel replacement image")
            
        if bpmkey!='BPM':
            update_head(image,'BPM',usebpm,"Bad pixel mask used")

######################################################################

def iqclip(inlist,lthresh=-16383.0,hthresh=65536.0,
           bookend=no,replace=0.0,
           maskin="!BPM",maskout="-bpm",maskval=1,
           clobber=globclob,verbose=globver):

    """ flexible image values clipping & BPM updating """

    # Parse inputs
    infiles=iraffiles(inlist)

    # Input mask specification
    bpmkey=None
    bpmfile=None
    re0=re.search("none",maskin,re.I)
    re1=re.search("^\!(\w+)",maskin)
    if re0 or len(maskin)==0:
        pass
    elif re1:
        bpmkey=re1.group(1)
    elif os.path.exists(maskin):
        bpmfile=maskin
    else:
        print "Failed to parse input mask specification '%s'" % maskin
        return

    # Output mask specification
    bpm2key=None
    bpm2sfx=None
    bpm2pfx=None
    bpm2file=None
    re0=re.search("none",maskout,re.I)
    re1=re.search("^\!(\w+)",maskout)
    re2=re.search("^(-\w+)",maskout)
    re3=re.search("(\w+-)$",maskout)
    if re0 or len(maskout)==0:
        pass
    elif re1:
        bpm2key=re1.group(1)
    elif re2:
        bpm2sfx=re2.group(1)
    elif re3:
        bpm2pfx=re2.group(1)
    elif len(infiles)>1:
        print "Cannot process multiple files to a single output mask"
        return
    else:
        bpm2file=maskout

    # Process each file in turn
    for image in infiles:

        check_exist(image,"r")

        # Identify BPM for this image
        usebpm=None
        if bpmkey:
            usebpm=get_head(image,bpmkey)
            check_exist(usebpm,"r")
        elif bpmfile:
            usebpm=bpmfile

        tmpimg=iraf.mktemp("iqclip")+".fits"
        if bookend:
            # Reset excess values to limits
            iraf.imexpr(expr="x < %f ? %f : (x > %f ? %f : x)" % \
                        (lthresh,lthresh,hthresh,hthresh),
                        output=tmpimg,x=image,refim="x",
                        dims="auto",intype="auto",outtype="auto",
                        verbose=verbose)
        else:
            # Reset excess values to the replacement setting
            iraf.imexpr(expr="(x > %f && x < %f ) ? x : %f" % \
                        (lthresh,hthresh,replace),
                        output=tmpimg,x=image,refim="x",
                        dims="auto",intype="auto",outtype="auto",
                        verbose=verbose)

        # Update BPM if necessary
        if usebpm:
            newbpm=""
            if bpm2file:
                newbpm=bpm2file
            elif bpm2key:
                newbpm=get_head(image,bpm2key)
            elif bpm2sfx:
                newbpm=image.replace('.fits',bpm2sfx+'.pl')
            elif bpm2pfx:
                newbpm=bpm2pfx+image.replace('.fits','.pl')
            if len(newbpm)>0:
                check_exist(newbpm,'w',clobber=clobber)
                iraf.imexpr(expr="(x > %f && x < %f ) ? y : %i" % \
                            (lthresh,hthresh,maskval),
                            output=newbpm,x=image,y=usebpm,refim="y",
                            dims="auto",intype="auto",outtype="auto",
                            verbose=verbose)

        # Replace the old image
        iraf.imdel(image,verify=no,go_ahead=yes)
        iraf.imrename(tmpimg,image,verbose=verbose)

        # Update FITS headers
        update_head(image,'IQCLIP',1,"Clipped with iqclip?")
        update_head(image,'CLIPLTH',lthresh,"Lower limit for clipping")
        update_head(image,'CLIPHTH',hthresh,"Upper limit for clipping")
        update_head(image,'CLIPREP',replace,"Replacement for clipping")

        if usebpm:
            update_head(image,'ORIGBPM',usebpm,"Original bad pixel mask")
            if bpmkey:
                update_head(image,bpmkey,newbpm)
            else:
                update_head(image,'BPM',newbpm)

######################################################################

def iqfringe(inlist, fringe, combine="average", fringeamp=0.05,
             skykey="SKYBKG", subkey="SKYSUB", bpmkey="MASKNAME",
             masktype="goodvalue", maskvalue=0, maskgrow=2.0,
             doclean=yes, clobber=globclob, verbose=globver):

    """ intelligent fringe-image creation"""

    # Packages
    iraf.crutil()

    # Defaults/constants
    scpfx="_S"
    sclist=iraf.mktemp("iqfringe")+".lis"

    # Parse inputs
    infiles=iraffiles(inlist)
    if fringe.endswith('.fits'):
        frngimg=fringe
    else:
        frngimg=fringe+'.fits'
    check_exist(frngimg,"w",clobber)
    if (combine!="average" and combine!="median"):
        combine="average"

    #if abs(fringeamp)>1:
    #    print "Please specify fringeamp as a fraction of the sky level"
    #    fringeamp=0.25

    # Rescale all images
    if verbose:
        print "Rescaling images by sky value"

    # Lists of files/masks
    scfiles=[]
    bpms=[]

    for image in infiles:

        check_exist(image,"r")
        sclimg=scpfx+image
        check_exist(sclimg,"w",yes)

        # Extract image sky value
        imgets(image,skykey,Stdout=1,Stderr=1)
        if (float(imgets.value) > 0):
            sky=float(imgets.value)
        else:
            print "No sky-value keyword %s in file %s" % (skykey,image)
            sys.exit(1)

        # Subtract sky if necessary
        if check_head(image,subkey):
            subdone=int(get_head(image,subkey))
            if not subdone:
                iqsubsky(image,sub=yes,skykey=skykey,subkey=subkey)

        # Divide by the sky value
        iraf.imarith(image,"/",sky,sclimg,verbose=no,noact=no)
        scfiles.append(sclimg)

        # Deal with the object mask
        bpmfile=get_head(sclimg,bpmkey)
        newbpm=sclimg.replace('.fits','-bpm.pl')
        check_exist(newbpm,'w',yes)
        if maskgrow>0:
            #
            tmpbpm=iraf.mktemp('iqfrg')+'.pl'
            if masktype=="goodvalue":
                iraf.imexpr("a!=%d ? 1 : 0" % maskvalue,
                            tmpbpm,a=bpmfile,verbose=no)
            elif masktype=="badvalue":
                iraf.imexpr("a==%d ? 1 : 0" % maskvalue,
                            tmpbpm,a=bpmfile,verbose=no)
            elif masktype=="none":
                continue
            # Grow the mask
            iraf.crgrow(tmpbpm,newbpm,radius=maskgrow,
                        inval=1,outval="INDEF")
            iraf.imdel(tmpbpm,verify=no,go_ahead=yes)
            # Reset the mask type and so forth
            masktype="goodvalue"
            maskvalue=0
        else:
            # Copy the mask
            iraf.imcopy(bpmfile,newbpm)

        # Update the scaled image
        update_head(sclimg,'BPM',newbpm)

        # Keep track
        bpms.append(newbpm)

    # Write the list of scaled images to a file
    check_exist(sclist,"w",yes)
    scin=open(sclist,"w")
    for sclimg in scfiles:
        scin.write("%s\n" % sclimg)
    scin.close()

    # Construct the imcombine statement
    imcombine.bpmasks=""
    imcombine.rejmask=""
    imcombine.sigma=""

    # Two possible combination approaches
    if (combine=="average"):
        imcombine.combine="average"
        imcombine.reject="avsigclip"
    else:
        imcombine.combine="median"
        imcombine.reject="none"

    imcombine.project=no
    imcombine.offsets="none"
    imcombine.masktype=masktype
    imcombine.maskvalue=maskvalue
    imcombine.blank=0.0

    imcombine.scale="none"
    imcombine.zero="none"
    imcombine.weight="!%s" % skykey
    imcombine.statsec=""
    imcombine.expname=""
    
    imcombine.lthreshold=-fringeamp
    imcombine.hthreshold=fringeamp
    imcombine.mclip=yes
    imcombine.lsigma=3.0
    imcombine.hsigma=3.0
    imcombine.sigscale=0.1

    # Execute the imcombine
    s=imcombine("@"+sclist,frngimg,Stdout=1)

    # Clean up listfile
    check_exist(sclist,"w",yes)

    # Update FITS header
    for image in infiles:
        ix=infiles.index(image)
        frgfkey="FRGFIL%02d" % (1+ix)
        update_head(frngimg,frgfkey,image,"File used to make fringe")
        if doclean:
            iraf.imdel(scpfx+image,verify=no,go_ahead=yes)

    # Delete temporary BPMs
    if doclean:
        for bpmfile in bpms:
            iraf.imdel(bpmfile,verify=no,go_ahead=yes)

    # The End
    if (verbose):
        print "Fringe creation finished, %s written" % frngimg
        
######################################################################

def iqdefringe(inlist, fringe, outpfx="f", skykey='SKYBKG',
               minlim=no, minval=-90000.,
               clobber=globclob, verbose=globver):

    """ intelligent fringe subtraction """

    if (len(outpfx) == 0):
        print "Warning: this will overwrite old files"

    # Parse inputs
    infiles=iraffiles(inlist)
    frgimg=fitsfile(fringe)

    # Defringe the images
    for image in infiles:

        check_exist(image,"r")

        # check processing status
	imgets(image,"FRGSUB",Stdout=1,Stderr=1)
        if (int(imgets.value) != 0):
            print "Warning: %s has already been fringe-subtracted!" % image
            iraf.imcopy(image,outpfx+image)
            continue

        if (verbose):
            print "Processing %s with fringe image %s" % (image,frgimg)

        # Get image sky value
        sky=""
        if check_head(image,skykey):
            try:
                sky=float(get_head(image,skykey))
            except:
                sky=""

        if sky=="":
            print "Sky value not present in image %s" % image
            return
            
        # Intelligent fringe subtraction
        if (len(outpfx) != 0):
            check_exist(outpfx+image,"w",clobber)

        if minlim:
            # Border method
            iraf.imexpr(expr="x < %.3f ? x : x - y*z" % minval,
                        output=outpfx+image,x=image,y=frgimg,z=sky,
                        verbose=verbose)
        else:
            # Simple method
            iraf.imexpr(expr="x - y*z",
                        output=outpfx+image,x=image,y=frgimg,z=sky,
                        verbose=verbose)
            
        # Update FITS header
        update_head(outpfx+image,"FRGSUB",1,"Has fringe-subtraction occured?")
        update_head(outpfx+image,"FRGNAME",frgimg,"Name of fringe image")

    if (verbose):
        print "Fringe-subtraction finished!"
        
######################################################################

def iqobjs(inlist, sigma, satval, skyval="!SKYBKG", masksfx="mask", 
           zeropt=25.0, wtimage="!FLATNAME", wtcut=0.1, fwhm=1.5, pix=0.3787,
           gain=2.3, aperture=10.0, minlim=no, minval=-100000., class_cut=0.80, 
           elong_cut=1.30, clobber=globclob, verbose=globver):

    """ intelligent object detection using Sextractor"""

    # Defaults/constants
    bkpfx="_B"
    
    try:
        reduction=os.environ['REDUCTION']
    except:
        print "Please define environment variable REDUCTION"
        print "This should point to the directory with the MSSSO tools."
        sys.exit(1)
        
    # Parse inputs
    infiles=iraffiles(inlist)

    # Parse user-specified weight image / header keyword
    wtkey=None
    re1=re.search("^\!(\w+)",wtimage)
    if re1:
        wtkey=re1.group(1)

    # Parse user-specified sky value / header keyword
    sky=0.0
    skykey=None
    re1=re.search("^\!(\w+)",skyval)
    if re1:
        skykey=re1.group(1)
    else:
        try:
            sky=float(skyval)
        except TypeError:
            sky=0.0

    # Big loop
    for image in infiles:

        if verbose:
            print "Analyzing image %s:" % image

        check_exist(image,"r")
        (imgroot,imgext)=os.path.splitext(image)
        starfile=image+".stars"
        check_exist(starfile,"w",clobber)
        mskimg=imgroot+"_"+masksfx+".fits"
        check_exist(mskimg,"w",yes)

        usewt=None
        if wtkey:
            usewt=get_head(image,wtkey)
        elif len(wtimage)>0 and wtimage!="none":
            usewt=wtimage

        if skykey:
            sky=float(get_head(image,skykey))

        # Reset border regions
        if minlim:
            if verbose:
                print "Zeroing border regions..."
            # Backup current image
            bkimage=bkpfx+image
            check_exist(bkimage,"w",yes)
            iraf.imcopy(image,bkimage)
            # Replace border with zeroes for Sextractor run
            tmpimg=iraf.mktemp("iqobj")+".fits"
            iraf.imexpr(expr="x < %.3f ? 0.0 : x" % minval,
                        output=tmpimg,x=image,
                        verbose=verbose)
            iraf.imdel(image,verify=no,go_ahead=yes)
            iraf.imrename(tmpimg,image,verbose=verbose)

        # Sextractor
        cmd="$REDUCTION/runsex.pl %s %.1f -sat %.2f -zp %f -mask %s -fwhm %.2f -pix %.3f -gain %.2f -aperture %.1f -theta -errs" % \
             (image,sigma,satval-sky,zeropt,mskimg,fwhm,pix,gain,aperture)

        if usewt:
            check_exist(usewt,"r")
            cmd=cmd+" -weight %s -wcut %.2f" % (usewt,wtcut)

        if verbose:
            print "Running sextractor with this command:"
            print "   > "+cmd

        fcmd=os.popen(cmd,'r')
        sexlines=fcmd.readlines()
        fcmd.close()

        # Read in Sextractor objects
        starfile=image+".stars"
        if os.path.exists(starfile):
            stars=Starlist(file=starfile)
            # FWHMs of objects with FLAG==0, class_star > 0.90, elong < 1.3
            flags=array(stars.flags())
            elongs=array(stars.elongs())
            fwhms=array(stars.fwhms())
            fwhmx=compress(logical_and(equal(flags,0),
                  less(elongs,elong_cut)),fwhms)
            seepix=median(fwhmx,minv=0.25,maxv=20.0)
            # Write objects to region file
            regfile=image.replace('.fits','.reg')
            stars.write(regfile,format="reg")
        else:
            print "Failed to find Sextractor output file!"
            seepix="INDEF"

        if verbose:
            print "Updating image and header keywords..."

        # Replace image with backup, if necessary
        if minlim:
            iraf.imdel(image,verify=no,go_ahead=yes)
            iraf.imrename(bkimage,image,verbose=verbose)

        # Update FITS headers
        update_head(image,"MASKNAME",mskimg,
                    "Object mask image from Sextractor")
        update_head(image,"STARFILE",starfile,
                    "Objects file from Sextractor")
        update_head(image,"ZEROPT",zeropt,
                    "Photometric zero-point used for Sextractor")
        update_head(image,"SEEPIX","%.3f" % seepix,
                    "Estimated seeing from Sextractor objects (pix)")
        update_head(image,"NSTARS",len(stars),
                    "Estimated number of objects from Sextractor")

        if verbose:
            print "done."

######################################################################

def iqseeing(inlist, stars="!STARFILE", seekey="SEEING",
             method="acorr", useflags=yes, skipstars=0, usestars=10, 
             boxsize=64, strictbox=no, fwhmmin=0.10, fwhmmax=20.0,
             pixscale=1.0, imgsec="", nbins=30, mclip=yes, lsigma=3.0,
             hsigma=3.0, update=yes, clobber=globclob, verbose=globver):

    """ calculate seeing for an image from the stars therein """

    # Input file list
    infiles=iraffiles(inlist)

    # Restricted image section?
    usesec=0
    if len(imgsec)>1:
        re1=re.search("\[(\d+):(\d+),(\d+):(\d+)\]",imgsec)
        if re1:
            usesec=1
            xmin=float(re1.group(1))
            xmax=float(re1.group(2))
            ymin=float(re1.group(3))
            ymax=float(re1.group(4))

    # Size of extraction region
    halfbox=int(boxsize/2)

    # User-specified starfile name / header keyword
    sfilekey=None
    starfile=None
    re1=re.search("^\!(\w+)",stars)
    re2=re.search("^\[\[",stars)

    if re1:
        # Starfile specified by a header keyword
        sfilekey=re1.group(1)
    elif re2:
        # Star parameters passed directly as a string
	try:
	    xstar=[]
            sxy=eval(stars)
            for i in range(0,len(sxy)):
                if len(sxy[i])==2:
                    [xx,yy]=sxy[i]
	   	    xstar.append(Star(len(xstar),xx,yy))
                elif len(sxy[i])>2:
                    [xx,yy,ww]=sxy[i]
	   	    xstar.append(Star(len(xstar),xx,yy,fwhm=ww))
        except:
            sys.stderr.write("Failed to parse STARS as 2D list.\n")
            return
    else:
	# Single starfile for all images
        starfile=stars
	check_exist(starfile,"r")
        xstar=Starlist(starfile)

    # Big loop
    for image in infiles:

        check_exist(image,"r")
        if verbose:
            print "Analyzing image %s:" % image

        # Read in stars from header-keyword specified starfiles
	if sfilekey:
	    starfile=get_head(image,sfilekey)
            xstar=Starlist(starfile)

	if method=="acorr":
            
	    # Autocorrelation analysis of regions near stars

	    # Read in the data from the image
            fimg=pyfits.open(image)
            fdata=fimg[0].data
            dims=fdata.shape

            # Loop over stars in image
            acfwhm=[]
            iskip=0
            istar=0

            while len(acfwhm)<usestars:

                # Getting or skipping too many stars...
                if istar==len(xstar):
                    break

                # The star in question
                star=xstar[istar]
                istar+=1
		[x,y]=[int(star.xval+0.5),int(star.yval+0.5)]
                [xb,yb]=boxinbox(x,y,boxsize,dims)

                # Skip flagged objects if requested
                if useflags and star.flag != 0:
                    continue

                # Select objects within region if requested
                if usesec and \
                       (star.xval<xmin or star.xval>xmax or \
                        star.yval<ymin or star.yval>ymax):
                    continue

                # Select against other stars in box if requested
                if strictbox and xstar.nstarsinbox([xb,yb])>2:
                    continue

                # Skip the number of stars requested
                if iskip<skipstars:
                    iskip += 1
                    continue

		# Extract region around star
                xd=fdata[yb[0]:yb[1],xb[0]:xb[1]]

		# Find width of autocorrelation
                acfwhm.append(azfwhm(xd,autocorr=1))

            # Take the median? of these
            if verbose:
                print "Estimating seeing as median of these values:"
                print "    ",list2str(acfwhm,"%.1f")

            seeing=median(acfwhm)
            fimg.close()
		
	else:

            # Statistical analysis of provided FWHMs
            if usesec:
                xstar=xstar.starsinbox([[xmin,xmax],[ymin,ymax]])

            fwhmx=array(xstar.fwhms())
            if useflags:
                flags=array(xstar.flags())
                fwhmx=compress(equal(flags,0),fwhmx)

            if len(fwhmx)<3:
                print "Very few acceptable stars in image, probably bad focus"
                fzero=2*fwhmmin
	    elif mclip:
		fzero=median(fwhmx,minv=fwhmmin,maxv=fwhmmax)
            else:
                fzero=mean(fwhmx,minv=fwhmmin,maxv=fwhmmax)

            # Perform a single sigma-clip iteration
            fsdev=sdev(fwhmx,zero=fzero)
            min2=fzero-lsigma*fsdev
            if min2<fwhmmin:
                min2=fwhmmin
            max2=fzero+hsigma*fsdev
            if max2>fwhmmax:
                max2=fwhmmax

            # Calculate the requested statistic
            if len(fwhmx)<3:
                seeing=median(fwhmx)
            elif method=="median":
                seeing=median(fwhmx,minv=min2,maxv=max2)
            elif method=="mode":
                seeing=mode(fwhmx,minv=min2,maxv=max2,nbins=nbins)
            else:
                print "Unrecognized method '%s' -- using median" % method
                seeing=median(fwhmx)

        # Update the image header
        if update:
            if verbose:
                print "Updating image header"

            update_head(image,seekey,seeing*pixscale,
                        "Estimated seeing (arcsec)")
            update_head(image,"IQSEEING",1,"Seeing calculated by iqseeing?")
            update_head(image,"SEEMETHD",method,"Seeing calculation method")

        # Print some information
        if verbose:
            print "Image %s has estimated seeing %.1f pixels (%.1f arcsec)" % \
                  (image,seeing,seeing*pixscale)

######################################################################

def iqfocus(infile, xvals, yvals, focusvals,
            method="acorr", boxsize=32, focuskey="BESTFOC",
            seekey="SEEING", fseekey="FOCSEE",
            update=yes, clobber=globclob, verbose=globver):

    """ calculate best-focus of a sequence of exposures on a single
        star """

    # Input filename
    image=fitsfile(infile)

    # Read in the data from the image
    fimg=pyfits.open(image)
    fdata=fimg[0].data
    dims=fdata.shape

    # Input checking
    [nx,ny,nf]=[len(xvals.split(',')),len(yvals.split(',')),
                len(focusvals.split(','))]

    if (nx>1 and ny>1 and nx!=ny) or (nf != max([nx,ny])):
        print "Inconsistent arrays xvals, yvals, focusvals"
        return

    # Convert input strings to float arrays
    if nx>1:
        xs=list2float(xvals.split(','))
    else:
        xs=nf * [float(xvals)]

    if ny>1:
        ys=list2float(yvals.split(','))
    else:
        ys=nf * [float(yvals)]

    fs=list2float(focusvals.split(','))

    # Seeing estimate for each image of the star
    seevals=[]

    # Loop over images of the star
    for i in xrange(nf):

        # Extract region around star
        [xb,yb]=boxinbox(xs[i],ys[i],boxsize,dims)
        xd=fdata[yb[0]:yb[1],xb[0]:xb[1]]

        if method=="acorr":
            # Find width of autocorrelation
            seevals.append(azfwhm(xd,autocorr=1))
        else:
            # Other methods?
            print "I don't know any other methods"
            return

    # Close the datafile
    fimg.close()

    # Work only with "good values"
    usefoc=[]
    usesee=[]
    for i in range(len(seevals)):
        seev=seevals[i]
        if seev<(boxsize/2-1):
            usesee.append(seev)
            usefoc.append(fs[i])

    # Need at least ~50% good values to get a good focus
    if len(usesee)<int(0.5*len(seevals)):
        print "Failed to get good FWHM values for most focus positions"
        return None

    # Best seeing is minimum (for now -- could do a fit, of course)
    bestsee=min(usesee)
    bestfoc=usefoc[usesee.index(bestsee)]

    # Update the image header
    if update:
        if verbose:
            print "Updating image header"

        update_head(image,focuskey,bestfoc,
                    "Best focus value from this series")
        update_head(image,seekey,bestsee,
                    "Best seeing value found in this series")
        update_head(image,fseekey,list2str(seevals,"%.2f"),"")
        update_head(image,"IQFOCUS",1,"Image processed by iqfocus?")
        update_head(image,"IQFMETHD",method,"iqfocus calculation method")

    # Print some information
    if verbose:
        print "Estimated seeing is %.1f pixels at focus value %.2f" % \
              (bestsee,bestfoc)

    # Return best-focus value
    return bestfoc

######################################################################

def iqwcs(inlist, objkey="OBJECT", rakey="RA", deckey="DEC",
          pixscl="!PIXSCALE", pixtol=0.05, nstar=40, nstarmax=40, 
          irmags=no,
          starfile='!STARFILE', catalog="web", ubhost="kronos",
          forcefit=no, forcecat=no, diffuse=no, class_cut=0.8,
          clobber=globclob, verbose=globver):

    """ Calculate WCS for an image by triangle matching """

    # Defaults/constants
    webmeth="WEB"
    irmeth="IR"
    irtext=""
    lclmeth="LOCAL"
    ubcatalog="asc_ref.cat"
    ubstars="asc_sextr.cat"
    ubhead="asc_param.dat"

    if os.environ.has_key('REDUCTION'):
        reduction=os.environ['REDUCTION']
    else:
        print "Please define environment variable REDUCTION"
        print "This should point to the directory with the MSSSO tools."
        sys.exit(1)
    
    # Target files
    infiles=iraffiles(inlist)

    # User-specified starfile / extension / header keyword
    starkey=None
    starext=None
    re1=re.search("^\!(\w+)",starfile)
    re2=re.search("^\.(\w+)",starfile)
    if re1:
        starkey=re1.group(1)
    elif re2:
        starext=re2.group(1)

    # Parse user-specified pixel scale / header keyword
    pixkey=None
    re1=re.search("^\!(\w+)",pixscl)
    if re1:
        pixkey=re1.group(1)

    # Parse user-specified scale uncertainty
    if pixtol>1.0:
        pixtol=0.10

    # User-specified astrometric catalog / retrieval method
    catmeth=None
    if catalog.upper()==webmeth:
        catmeth=webmeth
    elif catalog.upper()==irmeth:
        catmeth=irmeth
    elif catalog.upper()==lclmeth:
        catmeth=lclmeth
    else:
        check_exist(catalog,"r")

    # IR mags vs. Optical mags
    if irmags:
        irtext=" -ir"

    ####################
    
    # Big loop
    for image in infiles:

        if not forcefit and check_head(image,'IQWCS'):
            wcsstr=get_head(image,'IQWCS')
            try:
                if int(wcsstr):
                    if verbose:
                        print "Image %s has IQWCS set, skipping" % image
                    continue
            except:
                pass

        if verbose:
            print "Analyzing image %s..." % image

        check_exist(image,"r")
        (imgroot,imgext)=os.path.splitext(image)

        # Pixel scale
        if pixkey:
            pixscale=float(get_head(image,pixkey))
        else:
            pixscale=float(pixscl)

        # Find the starfile
        if starkey:
            starfile=get_head(image,starkey)
        elif starext:
            starfile=image+"."+starext
        check_exist(starfile,"r")

        # Read the starfile and write to the standard place
        stars=Starlist(starfile)
        check_exist(ubstars,'w',clobber)
        # Format for input to starfit_ub
        xcol=['objnum','xval','yval','dum','dum','mag','mag',
              'mag','mag','dum','theta','dum','dum','dum',
              'flag','elong','dum','class_star','dum','dum','dum']
        stars.write(ubstars,format="sex",cols=xcol,
                    maxnum=nstarmax,clobber=clobber)

        # Get various important parameters
        [object,ra,dec,naxis1,naxis2] = \
                get_head(image,[objkey,rakey,deckey,'naxis1','naxis2'])

        # Name of astrometric catalog file
        if forcecat:
            catfile=imgroot+"_cat.cat"
            check_exist(catfile,'w',yes)
        else:
            catfile=object.replace(" ","_")+".cat"

        # Size of retrieval request (based on image size)
        boxsize=1.1*array(imextent(image,arcmin=1))

        # Pointing coords
        objcoo=astrocoords(ra,dec)
        [rasxg,dcsxg]=objcoo.sxg()
        [radeg,dcdeg]=objcoo.deg()
        rahours=radeg/15.0

        # Find or Create the astrometric catalog file
        if not os.path.exists(catfile):
            if catmeth:
                if catmeth==webmeth:
                    cmd="getusnobn -c usnob -o %s -b 5.0 -f 19.0 -m I2 %s %s %f" \
                        % (catfile,rasxg,dcsxg,max(boxsize[0],boxsize[1]))
                elif catmeth==irmeth:
                    cmd="getusnobn -c nomad -o %s -b 5.0 -f 19.0 -m I2 %s %s %f" \
                        % (catfile,rasxg,dcsxg,max(boxsize[0],boxsize[1]))
                elif catmeth==lclmeth:
                    if len(ubhost)==0 or ubhost=="localhost":
                        catshort=catfile[0:catfile.rindex('.cat')]
                        cmd=("ubcone -F R -P %.7f -p %.7f -S %f -s %f "+ \
                             "-Z S -z M -I 5 -i 19 -L I2 -O %s") % \
                            (rahours,dcdeg,boxsize[0]/60,boxsize[1]/60,catshort)
                    else:
                        cmd=("rsh %s ubcone -F R -P %.7f -p %.7f -S %f "+ \
                             "-s %f -Z S -z M -I 5 -i 19 -L I2 -O - > %s") % \
                             (ubhost,rahours,dcdeg,boxsize[0]/60,boxsize[1]/60,
                              catfile)

                # Run the command
                if verbose:
                    print "Fetching astrometry catalog with this command:"
                    print "   > "+cmd
                fcmd=os.popen(cmd,'r')
                catlines=fcmd.readlines()
                fcmd.close()

            else:
                if verbose:
                    print "Using astrometry catalog %s" % catalog
                shutil.copy(catalog,catfile)
        else:
            if verbose:
                print "Using astrometry catalog %s" % catfile

        # Copy the astrometric catalog file to the right place
        check_exist(ubcatalog,'w',clobber)
        if os.path.exists(catfile):
            catlines=getlines(catfile)
            if len(catlines)>nstarmax:
                putlines(ubcatalog,catlines[0:nstarmax+27],clobber=clobber)
            else:
                shutil.copy(catfile,ubcatalog)
        else:
            print "Failed to retrieve astrometric catalog for %s" % image
            continue
        
        # Run starfit
        sfit_extra=""
        if diffuse:
            sfit_extra+=" -diffuse 0.0"
            
        cmd=("$REDUCTION/starfit %s -final -tan -xpix %s -ypix %s -scale " + \
             "%f -tolerance %.3f -nstar %i -edge 10") % \
             (sfit_extra,naxis1,naxis2,pixscale,pixtol,nstar)

        if verbose:
            print "Running starfit_ub with this command:"
            print "   > "+cmd

        fcmd=os.popen(cmd,'r')
        wcslines=fcmd.readlines()
        fcmd.close()

        # If that worked, add the WCS to the image header
        # (should insert some test logic here)
        success=0
        if len(wcslines)==0:
            # Bad failure
            print "starfit_ub run has failed, %s header unaltered" % \
                  image
        elif re.search('^Error',wcslines[0],re.I):
            # starfit_ub reported some error
            for line in wcslines:
                print line,
            print "starfit_ub run has failed, %s header unaltered" % \
                  image
        elif len(wcslines)>1 and re.search('^starfit',wcslines[1],re.I):
            # command line error
            print "starfit_ub command line error, %s header unaltered" % \
                  image
        else:
            # Assuming this means success...
            success=1
            if verbose:
                print "Successful fit, updating %s header" % image
            cmd="$REDUCTION/addhead %s -f %s" % (image,ubhead)
            fcmd=os.popen(cmd,'r')
            addheadlines=fcmd.readlines()
            fcmd.close()

            # Estimate seeing using catalog stars
            catstars=Starlist(catfile)
            catreg=catfile.replace('.cat','.reg')
            if not os.path.exists(catreg):
                catstars.write(catreg,format="reg")
            catstars.wcs2pix(image)
            match,catmatch=stars.match(catstars,tol=3.0,useflags=yes,
                                       image=image,maxnum=50)
            if len(match)>2:
                fwhmpix=median(match.fwhms())
                seeing="%.3f" % (fwhmpix*pixscale)
                fwhmpix="%.3f" % fwhmpix
            else:
                seeing="INDEF"
                fwhmpix="INDEF"
    
            # Update FITS headers
            update_head(image,"IQWCS",1,"Has iqwcs run successfully?")
            update_head(image,"EQUINOX",2000.0,"WCS Equinox")
            update_head(image,"CATALOG",catfile,
                        "Astrometric catalog for WCS")
            update_head(image,"SEEING",seeing,
                        "Estimated seeing in arcsec or INDEF")             
            update_head(image,"SEEPIX",fwhmpix,
                        "Estimated seeing from Sextractor objects (pix)")

        # Remove the temporary files (not during testing)
        if success:
            check_exist(ubstars,'w',yes)
            check_exist(ubcatalog,'w',yes)
            check_exist(ubhead,'w',yes)

######################################################################

def iqzeropt(inlist, catmag, starfile='!STARFILE', catalog="!CATALOG", 
             pixtol=3.0, useflags=yes, maxnum=50, method="mean", rejout=yes,
             fencelim=0.5, sigma=3.0, maxfrac=0.15, zptkey="ZEROPT",
             zpukey="ZEROPTU", clobber=globclob, verbose=globver):

    """ Determine photometric zero-point by comparison to catalog """

    # Target files
    infiles=iraffiles(inlist)

    # User-specified starfile / header keyword
    starkey=None
    re1=re.search("^\!(\w+)",starfile)
    if re1:
        starkey=re1.group(1)

    # Catalog starfile / header keyword
    catkey=None
    re2=re.search("^\!(\w+)",catalog)
    if re2:
        catkey=re2.group(1)

    ####################
    
    # Big loop
    for image in infiles:

        # Find the starfile
        if starkey:
            starfile=get_head(image,starkey)
        check_exist(starfile,"r")

        # Find the catalog file
        if catkey:
            catalog=get_head(image,catkey)
        check_exist(catalog,"r")

        # Read the two Starlists
        stars=Starlist(starfile)
        catstars=Starlist(catalog)

        # Check for presence of requested filter magnitudes
        if not catstars[0].mags.has_key(catmag):
            print "Failed to find filter '%s' in starlist" % catmag
            continue

        # Set up catalog list for matching
        catstars.wcs2pix(image)
        catstars.set_mag(catmag)
        catstars.mag_sort()

        # Get zero-point by matching the two lists
        dzpt,dzptu=stars.zeropt(catstars,tol=pixtol,useflags=useflags,
                                image=image,maxnum=maxnum,method="mean",
                                rejout=rejout,fencelim=fencelim,sigma=sigma,
                                maxfrac=maxfrac)

        # Convert delta zero-point to zero-point
        if check_head(image,zptkey):
            zeropt=float(get_head(image,zptkey))+dzpt
        else:
            zeropt=dzpt

        # Update image header
        update_head(image,"IQZEROPT",1,
                    "Has iqzeropt been run?")
        update_head(image,zptkey,"%.3f" % zeropt,
                    "Zero point by comparison to catalog")
        update_head(image,zpukey,"%.3f" % dzptu,
                    "Zero point uncertainty")

######################################################################

def iqcrzap(inlist, bpm="!BPM", outpfx="c", bpmupdate=no, 
            bpmthresh=5, clobber=globclob, verbose=globver):

    """ cosmic-ray zapping, adapted from Dave Kaplan's ircrzap"""

    if len(outpfx)==0:
        print "Warning: this will overwrite old files"

    # Defaults
    batchnum=20

    # Parse inputs
    infiles=iraffiles(inlist)

    # User-specified BPM file / header keyword
    bpmkey=None
    bpmfile=None
    if len(bpm)>1:
        re1=re.search("^\!(\w+)",bpm)
        if re1:
            bpmkey=re1.group(1)
        elif os.path.exists(bpm):
            bpmfile=bpm
        elif os.path.exists(bpm+'.pl'):
            bpmfile=bpm+'.pl'
        elif os.path.exists(bpm+'.fits'):
            check_exist(bpm+'.pl','w',clobber=yes)
            iraf.imcopy(bpm+'.fits',bpm+'.pl',verbose=verbose)
            bpmfile=bpm+'.pl'
        else:
            print "Couldn't find BPM file '%s'" % bpm
            return

    # Load package if necessary
    iraf.crutil()
    iraf.imutil()
            
    # craverage settings
    craverage=iraf.craverage
    craverage.average=""
    craverage.sigma=""
    craverage.navg=5
    craverage.nrej=0
    craverage.nbkg=5
    craverage.nsig=25
    craverage.var0=0
    craverage.var1=0
    craverage.var2=0
    craverage.crval=1
    craverage.lcrsig=10
    craverage.hcrsig=5
    craverage.crgrow=0
    craverage.objval=0

    ##############################

    allbpm=[]
    inpfiles=[]
    outfiles=[]
    crmasks=[]
    
    # Big Loop
    for image in infiles:
        # Check for file
        check_exist(image,'r')
        if len(outpfx)>0:
            check_exist(outpfx+image,'w',clobber)
        # Name for the output crmask
        (imgroot,imgext)=os.path.splitext(os.path.basename(image))
        crmask=outpfx+imgroot+"_crmask.pl"
        check_exist(crmask,'w',clobber)
        # check to see if has already been done
        if check_head(image,'IQCRZAP'):
            crstat=get_head(image,'IQCRZAP')
            if int(crstat):
                print "Warning: %s has already been zapped.  Skipping." % file
                continue
        # identify the BPM
        if bpmkey:
            bpmfile=get_head(image,bpmkey)
            if not os.path.exists(bpmfile):
                print "Failed to find BPM file '%s'" % bpmfile
                bpmfile=""
        if bpmfile and len(bpmfile)>1:
            check_exist(bpmfile,'r')
            #######
            # NOTE:  IRAF V2.12 cannot deal with input masks for craverage.
            # So we have to NOT DO THIS for now.
            #######iraf.imcopy(bpmfile,crmask,verbose=verbose)
            if bpmfile not in allbpm:
                allbpm.append(bpmfile)
        crmasks.append(crmask)
        # a little verbosity
        if verbose:
            print "Zap %s -> %s mask=%s" % (image,outpfx+image,crmask)

        # Lists of files
        inpfiles.append(image)
        check_exist(outpfx+image,'w',clobber=clobber)
        outfiles.append(outpfx+image)

    ##############################

    # Perform the craverage run
    inplis=iraf.mktemp("iqcrzapi")+".lis"
    outlis=iraf.mktemp("iqcrzapo")+".lis"
    crmlis=iraf.mktemp("iqcrzapc")+".lis"

    putlines(inplis,inpfiles)
    putlines(outlis,outfiles)
    putlines(crmlis,crmasks)
    
    if verbose:
        print "Running craverage...",

    iraf.craverage(input="@"+inplis,output="@"+outlis,crmask="@"+crmlis)

    if verbose:
        print "done."

    ######
    # IRAF v2.12 modification:  Update CR masks with BPM pixels too
    ##### NOT TESTED YET!!!
    if bpmfile or bpmkey:
        crtemp=iraf.mktemp("iqcrzapm")+".pl"

        for i in xrange(len(inpfiles)):
            image=inpfiles[i]
            crmask=crmasks[i]
            if bpmkey:
                bpmfile=get_head(image,bpmkey)
                if not os.path.exists(bpmfile):
                    print "Failed to find BPM file '%s'" % bpmfile
                    bpmfile=""
            iraf.imarith(crmask,'+',bpmfile,crtemp)
            iraf.imdel(crmask)
            iraf.imrename(crtemp,crmask)

    # update headers
    for outimg in outfiles:
        update_head(outimg,'IQCRZAP',1,"iqcrzap has been run on this image")
        update_head(outimg,'CRMASK',crmask,"Cosmic-ray + bad pixel mask")

    ##############################

    # We can only update the BPM if it's the same for all images
    if len(allbpm)>1:
        if bpmupdate:
            print "BPM updating requested, but more than one static"+ \
                  "BPM was supplied or referenced"
            bpmupdate=no
    elif len(allbpm)==1:
        bpmfile=allbpm[0]
    else:
        if bpmupdate:
            print "BPM updating requested, but no static BPMs were supplied"
            bpmupdate=no

    # BPM updating:  If a single pixel is bad in more than N files, it
    # belongs in the static BPM. 
    if bpmupdate:

        if verbose:
            print "Updating static BPM..."

        # output BPM and temporary BPM
        tmpbpm1=iraf.mktemp("iqcrzp")+".pl"
        tmpbpm2=iraf.mktemp("iqcrzp")+".pl"

	# sum up all individual masks
        tmplist=iraf.mktemp("iqcrzp")+".list"
        tmpout=open(tmplist,"w")
        tmpout.write("\n".join(crmasks)+"\n")
        tmpout.close()
        imsum=iraf.imsum
        imsum.option="sum"
        imsum.low_rej=0
	imsum.high_rej=0
	iraf.imsum(input='@'+tmplist,output=tmpbpm1,verbose=verbose)
        check_exist(tmplist,'w',yes)

        # find those pixels that have > bpmthresh flags
        iraf.imexpr("x > %d ? 1 : 0" % bpmthresh,
                    output=tmpbpm2,x=tmpbpm1,verbose=verbose)
        iraf.imdel(tmpbpm1,verify=no,go_ahead=yes)
        iraf.imrename(tmpbpm2,tmpbpm1)

        # combine with BPM
        iraf.imexpr("x+y > 1 ? 1 : 0",
                    output=tmpbpm2,x=tmpbpm1,y=bpmfile,
                    verbose=verbose)
        iraf.imdel(tmpbpm1,verify=no,go_ahead=yes)

        # Backup old BPM & replace
        check_exist(bpmfile+".old",'w',clobber)
        os.rename(bpmfile,bpmfile+".old")
        os.rename(tmpbpm2,bpmfile)
        update_head(bpmfile,'UPDATED',1,"BPM updated by iqcrzap")

######################################################################

def iqcoadd(inlist, output, medsfx="-med", xshext="xsh",
            bpm="!BPM", exppfx="E", xrgpfx="X",
            bpmsfx="-bpm", crmsfx="-crm", expkey="EXPTIME",
            xcregion="[0.2:0.8,0.2:0.8]",window=15,
            crmeth="sigma", crthresh=5.0, overres=2.0, useold=yes,
            medonly=no, trimout=yes, trimmax=no, subpixel=yes,
            doclean=yes, clobber=globclob, verbose=globver):

    # Packages
    iraf.imgeom()
    iraf.imutil()
    iraf.crutil()

    # Some variables
    skysec="SKYSEC"  # header keyword to track sky area on-image
    pixsafe=5        # extra space to allow for xregister shifts
    wtsfx=".weight"  # suffix for weight image
    xcpfx="_X"       # prefix for Xregister images

    # Input checking
    if len(exppfx)==0:
        print "Expanded images are not allowed to overwrite unexpanded images"
        print "Setting exppfx to default 'E'"
        exppfx="E"

    if len(xrgpfx)==0:
        print "Xregister images are not allowed to overwrite expanded images"
        print "Setting xrgpfx to default 'X'"
        xrgpfx="X"

    if len(bpmsfx)==0:
        print "BPM names cannot be the same as expanded image names"
        print "Setting bpmsfx to default '-bpm'"
        bpmsfx="-bpm"

    if len(crmsfx)==0:
        print "Cosmic ray masknames cannot be the same as expanded image names"
        print "Setting crmsfx to default '-crm'"
        bpmsfx="-crm"

    # Region for Xregister (cross-correlation)
    re1=re.search("\[([\d\.]+):([\d\.]+),([\d\.]+):([\d\.]+)\]",xcregion)
    if re1:
        try:
            xclx=float(re1.group(1))
            xcrx=float(re1.group(2))
            xcly=float(re1.group(3))
            xcuy=float(re1.group(4))
        except:
            print "Trouble parsing XCREGION, reverting to default"
            xclx=0.2; xcrx=0.8; xcly=0.2; xcuy=0.8
    else:
        print "Using full overlap region for Xregister"
        xclx=0.0; xcrx=1.0; xcly=0.0; xcuy=1.0
    if verbose:
        print "Xcregion:  [%.2f:%.2f,%.2f:%.2f]" % \
              (xclx,xcrx,xcly,xcuy)

    # User-specified BPM file / header keyword
    bpmkey=None
    bpmpl=None
    re1=re.search("^\!(\w+)",bpm)
    re2=re.search("none",bpm,re.I)
    if re2:
        print "No BPM available"
    elif re1:
        bpmkey=re1.group(1)
    elif os.path.exists(bpm):
        bpmpl=bpm
    elif os.path.exists(bpm+'.pl'):
        bpmpl=bpm+'.pl'
    elif os.path.exists(bpm+'.fits'):
        check_exist(bpm+'.pl','w',clobber=yes)
        iraf.imcopy(bpm+'.fits',bpm+'.pl',verbose=verbose)
        bpmpl=bpm+'.pl'
    else:
        print "Couldn't find BPM file '%s' so no static BPM will be used" % bpm

    # Do we have a static BPM?
    dobpm=(bpmkey or (bpmpl and len(bpmpl)>0))

    # Parse inputs
    infiles=iraffiles(inlist)
    nfiles=len(infiles)
    if nfiles<2:
        print "Found less than two input files in iqcoadd -- quitting"
        return

    # Create Xshifts filename from output filename
    if len(xshext)==0:
        xshfile=iraf.mktemp('iqcoadd')+'.xsh'
    else:
        xshfile=output.replace('.fits','.'+xshext)

    ##############################

    # Coordinates of lower-left corner of first image
    img0=infiles[0]
    [[ra0,dc0]]=impix2wcs(img0,0,0)
    [nsmlx,nsmly]=get_head(img0,['NAXIS1','NAXIS2'])

    # Accumulate list of pixel locations for that coordinate
    # (these are the pixel shifts, basically)
    llx=[]
    lly=[]
    for image in infiles:
        [[x0,y0]]=imwcs2pix(image,ra0,dc0)
        llx.append(-x0)
        lly.append(-y0)

    # Determine the size and "location" of the coadded image,
    # assuming only these rectilinear shifts
    [meanx,meany]=[mean(llx),mean(lly)]
    [minx,miny]=[int(round(min(llx))),int(round(min(lly)))]
    [maxx,maxy]=[int(round(max(llx))),int(round(max(lly)))]
    [nbigx,nbigy]=[nsmlx+maxx-minx+2*pixsafe,
                   nsmly+maxy-miny+2*pixsafe]

    # Find the image that is closest to the mean position
    dmin=math.hypot(meanx,meany)+10
    for i in xrange(nfiles):
        dist=math.hypot(llx[i]-meanx,lly[i]-meany)
        if dist<dmin:
            dmin=dist
            imin=i
    refimg=infiles[imin]
    refxy=[round(llx[imin])-minx+pixsafe,
           round(lly[imin])-miny+pixsafe]
    refnew=exppfx+refimg

    # Convert all images to their expanded versions
    shifts=[]
    expfiles=[]
    for i in xrange(nfiles):

        # Expanded image names
        image=infiles[i]
        expimg=exppfx+image
        expfiles.append(expimg)

        # Expand to the full size we will be using
        if useold and os.path.exists(expimg):
            print "Using preexisting image %s" % expimg
            # Location of the detector within each expanded image,
            # as encoded in previous processing
            oldx0,oldy0=1,1
            skyval=get_head(expimg,skysec)
            resky=re.search("\[(\d+):(\d+),(\d+):(\d+)\]",skyval)
            if resky:
                oldx0,oldx1,oldy0,oldy1=int(resky.group(1)), \
                                        int(resky.group(2)), \
                                        int(resky.group(3)), \
                                        int(resky.group(4))
            else:
                print "Failed to retrieve %s values from old images." % \
                      skysec.upper()
                print "Cannot use these images, please set useold=no"
                return
            # Convert from IRAF-style index to Python-style index
            shifts.append([oldx0-1,oldy0-1])
        else:
            if dobpm:
                # "Fixpix" the image (should help xregister)
                if bpmkey:
                    bpmpl=get_head(image,bpmkey)
                iqmask(image,mask=bpmpl,method='fixpix',
                       clobber=clobber,verbose=verbose)
            
            # Expand the image
            xy=[round(llx[i])-minx+pixsafe,
                round(lly[i])-miny+pixsafe]
            shifts.append(xy)
            expand_image(image,expimg,size=[nbigx,nbigy],llxy=xy,
                         border=0,bpmkey="",skysec=skysec,
                         clobber=clobber,verbose=verbose)

    ##############################

    # XREGISTER cross-correlation

    # Determine the region of full-overlap of the images
    fullxrg=[maxx-minx+pixsafe,nsmlx+pixsafe]
    fullyrg=[maxy-miny+pixsafe,nsmly+pixsafe]
    fullreg="[%d:%d,%d:%d]" % (fullxrg[0]+1,fullxrg[1],
                               fullyrg[0]+1,fullyrg[1])

    # Narrow down the region prior to Xregister
    xzero=fullxrg[0]+1; xwide=fullxrg[1]-xzero+1
    yzero=fullyrg[0]+1; ytall=fullyrg[1]-yzero+1
    xcxrg=[round(xzero+xclx*xwide),round(xzero+xcrx*xwide)]
    xcyrg=[round(yzero+xcly*ytall),round(yzero+xcuy*ytall)]
    xcreg="[%d:%d,%d:%d]" % (xcxrg[0],xcxrg[1],
                             xcyrg[0],xcyrg[1])

    # Create listfile for IRAF that excludes refimage
    expfiles2=copy.deepcopy(expfiles)
    expfiles2.remove(refnew)

    # Perform Xregister, if necessary
    if (not useold) or (not os.path.exists(xshfile)):
        if verbose:
            print ("Running Xregister over region %s " +
                   "of expanded images") % xcreg

        # Take arcsinh of the images
        for image in expfiles:
            asinh_image(image,xcpfx+image)
        xexpfiles2=pfx_list(expfiles2,xcpfx)
        xrefnew=xcpfx+refnew

        # Construct list of images
        xrgin=iraf.mktemp("iqxrgin")+".lis"
        putlines(xrgin,xexpfiles2)
        check_exist(xshfile,'w',clobber)

        # Run Xregister
        try:
            # this parameter deleted in more recent IRAF versions
            iraf.xregister.databasefmt=no
        except:
            pass
        
        iraf.xregister('@'+xrgin,xrefnew,xcreg,xshfile,output="",
                       coords="",xlag=0,ylag=0,
                       dxlag=0,dylag=0,background="none",
                       apodize=0.0,filter="none",correlation="fourier",
                       xwindow=window,ywindow=window,function="centroid",
                       interactive=no,verbose=verbose)

        # Clean up
        os.remove(xrgin)
        # Delete asinh files
        os.remove(xrefnew)
        for image in xexpfiles2:
            os.remove(image)

    # Order "xshifts" the same as expfiles
    xshifts=copy.deepcopy(shifts)
    xshifts[expfiles.index(refnew)]=[0.0,0.0]
    # Read in the xregister shifts
    # Format of the xregister shifts file:
    # <filename> <xshift> <yshift>
    xshlines=getlines(xshfile)
    for xline in xshlines:
        xels=xline.split()
        xels[0]=xels[0].replace(xcpfx,'')
        ix=expfiles.index(xels[0])
        xshifts[ix]=[float(xels[1]),float(xels[2])]
    # Delete the Xregister shifts file, if temporary
    if len(xshext)==0:
        os.remove(xshfile)

    # Create eshifts, a list of the full integer shifts (including
    # original offset to expanded size) 
    eshifts=shifts
    for i in xrange(nfiles):
        eshifts[i]=[shifts[i][0]+round(xshifts[i][0]),
                    shifts[i][1]+round(xshifts[i][1])]

    # Shift images by integral-pixel amounts to Xregister positions
    xfiles=pfx_list(expfiles,xrgpfx,check='w',clobber=clobber)
    for expfile in expfiles2:
        ix=expfiles.index(expfile)
        dx,dy=xshifts[ix]
        if round(dx)!=0 or round(dy)!=0:
            # Integral shift of expanded image
            shift_image(expfile,xfiles[ix],[round(dx),round(dy)],
                        border=0,bpmkey="",skysec=skysec)
            # Replace former expfile with this one
            iraf.imdel(expfile,verify=no,go_ahead=yes)
            iraf.imrename(xfiles[ix],expfile)
            # Check for image-specific BPM
            newbpm=expfile.replace('.fits',bpmsfx+'.pl')
            if os.path.exists(newbpm):
                tmpbpmf=iraf.mktemp("iqcoaddbpmf")+'.fits'
                newbpmf=mkbpmfits(newbpm,clobber=yes,delete=no)
                shift_image(newbpmf,tmpbpmf,[round(dx),round(dy)],
                            border=1,bpmkey="",skysec=skysec)
                iraf.imdel(newbpmf,verify=no,go_ahead=yes)
                iraf.imdel(newbpm,verify=no,go_ahead=yes)
                iraf.imcopy(tmpbpmf,newbpm,verbose=no)
                iraf.imdel(tmpbpmf,verify=no,go_ahead=yes)

    # Update information in the Xshifts file
    if len(xshext)!=0:
        os.remove(xshfile)
        xshlines=[]
        for i in xrange(nfiles):
            xx,xy=xshifts[i]
            dx,dy=xx-round(xx),xy-round(xy)
            xshlines.append("%s %+.3f %+.3f" % (expfiles[i],dx,dy))
        # Write the new Xshifts file
        putlines(xshfile,xshlines,clobber=yes)

    ##############################

    # Expand & shift the BPMs for each image
    allbpms={}
    if dobpm:
        for i in xrange(nfiles):
            image=expfiles[i]
            if bpmkey:
                bpmpl=get_head(image,bpmkey)
                if len(bpmpl)==0:
                    print "Failed to fetch BPM header keyword"
                    bpmpl=None
		    dobpm=0
            if not bpmpl:
                continue
            if not os.path.exists(bpmpl):
                print "Failed to find specified BPM file %s" % bpmpl
                dobpm=0
                continue
            # Make the BPM into a FITS file
            allbpms[bpmpl]=1
            newbpm=image.replace('.fits',bpmsfx+'.pl')
            if bpmpl==newbpm:
                print "Using preexisting BPM %s" % newbpm
                continue
            bpmfits=bpmpl.replace('.pl','.fits')
            if not os.path.exists(bpmfits):
                mkbpmfits(bpmpl,clobber=yes,delete=no)
            # Expand the BPM appropriately
            newbpmf=image.replace('.fits',bpmsfx+'.fits')
            expand_image(bpmfits,newbpmf,size=[nbigx,nbigy],
                         llxy=eshifts[i],border=1,skysec=skysec,
                         clobber=clobber,verbose=no)
            # Convert expanded BPM into a PL file
            newbpm=mkbpmpl(newbpmf,clobber=yes,delete=yes)
            iraf.imarith(newbpm,'min',1,newbpm,noact=no,verbose=no)
            iraf.imarith(newbpm,'max',0,newbpm,noact=no,verbose=no)
            # Update image header
            update_head(image,'BPM',newbpm)

    ##############################

    # Flag cosmic rays &tc by comparison to median image

    # Execute median / minimum combination
    medinp=iraf.mktemp('iqmedinp')+'.lis'
    putlines(medinp,expfiles)
    if len(medsfx)==0:
        medout=iraf.mktemp('iqmedout')+'.fits'
    else:
        medout=output.replace('.fits',medsfx+'.fits')
        check_exist(medout,'w',clobber)

    # Working with BPM or not
    if dobpm:
        imcombine.masktype="goodvalue"
        imcombine.maskvalue=0
        # Rejected pixel mask name
        medrej=iraf.mktemp('iqcoaddrj')+'.pl'
        check_exist(medrej,'w',clobber)
        imcombine.rejmasks=medrej
    else:
        imcombine.masktype="none"
        imcombine.rejmasks=""

    # Median vs. Minimum combination
    if nfiles<4:
        # NOTE: this method of "minimum combination" is not
        # extensively tested...
        imcombine.reject="minmax"
        imcombine.nlow=0
        imcombine.nhigh=nfiles-1
        imcombine.nkeep=1
    else:
        imcombine.reject="none"

    # Median / minimum combination
    iraf.imcombine('@'+medinp,medout,combine="median",
                   project=no,offsets="none",scale="exposure",zero="none",
                   weight="exposure",statsec="",expname="EXPTIME",
                   lthresh="INDEF",hthresh="INDEF")

    # Create weight image if we used BPMs at this stage
    if dobpm:
        medweight=medout.replace('.fits',wtsfx+'.fits')
        check_exist(medweight,'w',clobber=clobber)
        rejlines=[]
        for i in xrange(nfiles):
            rejlines.append(medrej+'[*,*,%d]' % (i+1))
        rejlis=iraf.mktemp('iqcoaddrj')+'.lis'
        putlines(rejlis,rejlines)
        iraf.imsum('@'+rejlis,medweight,title="",hparams="",
                   pixtype="",calctype="",option="sum",
                   low_reject=0,high_reject=0,verbose=no)
        iraf.imarith(nfiles,'-',medweight,medweight,noact=no,verbose=no)
        os.remove(rejlis)
        iraf.imdel(medrej,verify=no,go_ahead=yes)
    else:
        # This variable is used below
        medweight=nfiles

    # Delete IRAF listfile
    os.remove(medinp)

    ##################################################

    # User can skip CR-rejection if they want (not recommended)
    if crmeth != "none" and not medonly:

        crmasks=[]

        # Identify flaws / cosmic rays by comparison to median image
        for i in xrange(nfiles):

            # Name for cosmic ray pixel list
            crmask=expfiles[i].replace('.fits',crmsfx+'.pl')
            check_exist(crmask,'w',clobber)
            crmasks.append(crmask)

            # Check if we already have one
            if useold and check_head(expfiles[i],'CRMASK'):
                if get_head(expfiles[i],'CRMASK')==crmask:
                    print "Using preexisting CRMASK %s" % crmask
                    continue

            # (Image+Sky) / (Median+Sky) = Ratio
            mrat=iraf.mktemp("iqcoaddr")+'.fits'
            skybkg=float(get_head(expfiles[i],'SKYBKG'))
            iraf.imexpr("(x+%.1f)/(y+%.1f)" % (skybkg,skybkg),
                        output=mrat,x=expfiles[i],y=medout,
                        verbose=verbose)

            # Create CRMASK by specified means
            if crmeth=='ccdmask':
                # Use ccdmask to convert ratio image to BPM
                iraf.ccdmask(mrat,crmask,ncmed=7,nlmed=7,ncsig=15,nlsig=15,
                             lsigma=crthresh,hsigma=crthresh,ngood=5,
                             linterp=1,cinterp=1,eqinterp=1)
            elif crmeth=='craverage':
                # Use craverage to convert ratio image to BPM
                junk=iraf.mktemp("iqcoaddj")+'.fits'
                iraf.craverage(mrat,junk,crmask=crmask,average="",
                               sigma="",navg=5,nrej=1,nbkg=5,nsig=25,
                               var0=0,var1=0,var2=0,crval=1,lcrsig=2*crthresh,
                               hcrsig=crthresh,crgrow=0,objval=1,
                               lobjsig=2*crthresh,hobjsig=crthresh,objgrow=0)
                iraf.imdel(junk,verify=no,go_ahead=yes)
            elif crmeth=='lacos':
                # Use LA_Cosmic to convert ratio image to BPM
                junk=iraf.mktemp("iqcoaddj")+'.fits'
                iraf.lacos_im(mrat,junk,crmask,gain=0,readn=0,
                              statsec=fullreg[1:-2],skyval=0,sigclip=crthresh,
                              sigfrac=0.5,objlim=1.0,niter=4,verbose=no)
                iraf.imdel(junk,verify=no,go_ahead=yes)
            elif crmeth=='sigma':
                # Apply two-sided sigma criterion
                iraf.iterstat(mrat+fullreg,nsigrej=4.0,maxiter=3,
                              verbose=no)
                rms=float(iraf.iterstat.sigma)
                cutoff=crthresh*rms
                iraf.imexpr("abs(x-1) > %.2f ? 1 : 0" % cutoff,
                            crmask,x=mrat,verbose=no)
            elif crmeth=='sighi':
                # Apply one-sided sigma criterion
                iraf.iterstat(mrat+fullreg,nsigrej=4.0,maxiter=3,
                              verbose=no)
                rms=float(iraf.iterstat.sigma)
                cutoff=1+crthresh*rms
                iraf.imexpr("x > %.2f ? 1 : 0" % cutoff,
                            crmask,x=mrat,verbose=no)
            elif crmeth=='thresh':
                # Apply two-sided value threshold
                cutoff=crthresh
                iraf.imexpr("abs(x-1) > %.2f ? 1 : 0" % \
                            (medval,cutoff),
                            crmask,x=mrat,verbose=no)
            elif crmeth=='thrhi':
                # Apply one-sided value threshold
                cutoff=1+crthresh
                iraf.imexpr("x > %.2f ? 1 : 0" % (medval,cutoff),
                            crmask,x=mrat,verbose=no)
            else:
                print "Unrecognized crmeth '%s'" % crmeth
                return

            # Remove bad pixels from CRMASK
	    # !!! This step crashes for input files without BPMs !!!
            if dobpm:
                newbpm=get_head(expfiles[i],'BPM')
                iraf.imarith(crmask,'-',newbpm,crmask,noact=no,verbose=no)
                # Some CRMASK pixels will now be negative!  Clean up.
                iraf.imarith(crmask,'max',0,crmask,noact=no,verbose=no)

            # Header updates
            update_head(expfiles[i],'CRMASK',crmask)

            # Clean up
            iraf.imdel(mrat,verify=no,go_ahead=yes)

        ##############################

        # Distinguish bad pixels / cosmics / objects by comparison
        # among the images.

        # Stack in WCS space to identify objects
        wcsstkl=iraf.mktemp("iqcoaddwst")+".lis"
        wcsstk=iraf.mktemp("iqcoaddwst")+".fits"
        objbpm=iraf.mktemp("iqcoaddobj")+".pl"
        putlines(wcsstkl,crmasks)
        iraf.imsum("@"+wcsstkl,wcsstk,title="",hparams="",
                   pixtype="",calctype="",option="sum",
                   low_reject=0,high_reject=0,verbose=no)
        # Magic threshold for identification as object
        iraf.imexpr("x >= sqrt(y+1) ? 1 : 0",
                    objbpm,x=wcsstk,y=medweight,verbose=no)
        check_exist(wcsstkl,'w',clobber=yes)
        iraf.imdel(wcsstk,verify=no,go_ahead=yes)
        # Grow the object BPM to make a better mask
        iraf.crgrow(objbpm,objbpm,radius=1.5,
                    inval="INDEF",outval="INDEF")
        
        # Exclude objects from CRMASKs
        for i in xrange(nfiles):

            crmask=crmasks[i]
            iraf.imarith(crmask,'-',objbpm,crmask,noact=no,verbose=no)
            # Some CRMASK pixels will now be negative!  Clean up.
            iraf.imarith(crmask,'min',1,crmask,noact=no,verbose=no)
            iraf.imarith(crmask,'max',0,crmask,noact=no,verbose=no)

        # Stack CRMASKs in detector space to identify bad pixels
        tmpdimgs=[]
        detlines=[]
        for i in xrange(nfiles):
            llxy=eshifts[i]
            # IRAF demands absolutely identical image dimensions, so
            # if one of the image sections (due to Xregister
            # corrections) now  extends "outside the box" of the
            # expanded images, we have to create a new image.
            # This is not technically correct (we should pad with 0's)
            # but is forgivable since we are only after a vote on each
            # pixel, anyway.
            if llxy[0]<0 or llxy[1]<0 or \
               llxy[0]+nsmlx>nbigx or llxy[1]+nsmly>nbigy:
                print "Note: sky section has been clipped for %s" % expfiles[i]
                llxy[0]=max(llxy[0],0)
                llxy[1]=max(llxy[1],0)
                if llxy[0]+nsmlx > nbigx:
                    llxy[0]=nbigx-nsmlx
                if llxy[1]+nsmly > nbigy:
                    llxy[1]=nbigy-nsmly
            # Convert from Python-style to IRAF-style indexing
            thisreg="[%d:%d,%d:%d]" % (llxy[0]+1,llxy[0]+nsmlx,
                                       llxy[1]+1,llxy[1]+nsmly)
            detlines.append(crmasks[i]+thisreg)

        # Perform detector-space stack
        detstkl=iraf.mktemp("iqcoadddst")+".lis"
        detstk=iraf.mktemp("iqcoadddst")+".fits"
        detbpm=iraf.mktemp("iqcoaddbpm")+".pl"
        putlines(detstkl,detlines)
        iraf.imsum("@"+detstkl,detstk,title="",hparams="",
                   pixtype="",calctype="",option="sum",
                   low_reject=0,high_reject=0,verbose=no)
        # Magic number for identification
        nstack=ceil(sqrt(nfiles+1))
        iraf.imexpr("x >= %i ? 1 : 0" % nstack,
                    detbpm,x=detstk,verbose=no)
        check_exist(detstkl,'w',clobber=yes)
        iraf.imdel(detstk,verify=no,go_ahead=yes)
        
        # Produce updated BPM if we can
        newbpm=output.replace('.fits',bpmsfx+'.pl')
        if dobpm and len(allbpms)==1:
            # Replace former BPM
            oldbpm=allbpms.keys()[0]
	    tmpbpmx=iraf.mktemp("iqcoaddbpm")+".pl"
            iraf.imexpr("x+y > 0 ? 1 : 0",tmpbpmx,
                        x=oldbpm,y=detbpm,verbose=no)
            check_exist(newbpm,'w',clobber=clobber)
	    iraf.imrename(tmpbpmx,newbpm)
        
        # Update individual BPMs to include CRMASKs
        for i in xrange(nfiles):

            # Current CRMASK and active detector region
            crmask=crmasks[i]

            # Image CRMASK is cosmics only; Image BPM is bad pixels +
            # cosmics.  (If you want bad pixels only, subtract CRMASK
            # from BPM.)

            # Update image and BPM
            if dobpm:
                # Create updated BPM (no cosmics)
                newbpm=get_head(expfiles[i],'BPM')

                # Incorporate newly-identified bad pixels
                llxy=eshifts[i]
                urxy=[llxy[0]+nsmlx,llxy[1]+nsmly]
                detreg=""
                # Careful calculation for images with clipped sky sections
                if llxy[0]<0 or llxy[1]<0 or \
                       urxy[0]>nbigx or urxy[1]>nbigy:
                    print "Note: sky section has been clipped for %s" % \
                          expfiles[i]
                    ddx,ddy=0,0
                    if llxy[0]<0:
                        ddx=-llxy[0]
                        llxy[0]=0
                    elif urxy[0]>nbigx:
                        urxy[0]=nbigx
                    if llxy[1]<0:
                        ddy=-llxy[1]
                        llxy[1]=0
                    elif urxy[1]>nbigy:
                        urxy[1]=nbigy
                    szx,szy=urxy[0]-llxy[0],urxy[1]-llxy[1]
                    # IRAF-style indexing
                    detreg="[%d:%d,%d:%d]" % (ddx+1,ddx+szx,
                                              ddy+1,ddy+szy)
                # Convert from Python-style to IRAF-style indexing
                thisreg="[%d:%d,%d:%d]" % (llxy[0]+1,urxy[0],
                                           llxy[1]+1,urxy[1])
                tmpbpm=iraf.mktemp("iqcoaddb")+'.pl'
                iraf.imarith(newbpm+thisreg,'+',detbpm+detreg,tmpbpm,
                             noact=no,verbose=no)
                iraf.imcopy(tmpbpm,newbpm+thisreg,verbose=no)
                iraf.imarith(newbpm,'min',1,newbpm,noact=no,verbose=no)
                iraf.imdel(tmpbpm,verify=no,go_ahead=yes)

            else:
                # The new BPM is all we've got
                detbpmf=detbpm.replace('.pl','.fits')
                if not os.path.exists(detbpmf):
                    mkbpmfits(detbpm,clobber=yes,delete=no)
                newbpmf=expfiles[i].replace('.fits',bpmsfx+'.fits')
                expand_image(detbpmf,newbpmf,size=[nbigx,nbigy],
                             llxy=eshifts[i],border=1,skysec=skysec,
                             clobber=clobber,verbose=no)
                # Convert expanded BPM into a PL file
                newbpm=mkbpmpl(newbpmf,clobber=yes,delete=yes)
                # Update image header
                update_head(expfiles[i],'BPM',newbpm)

            # Add old BPM to CRMASK to get new BPM
            iraf.imarith(newbpm,'+',crmask,newbpm,noact=no,verbose=no)
            iraf.imarith(newbpm,'min',1,newbpm,noact=no,verbose=no)

            # Propagate new BPM back to original image
            orgbpm=infiles[i].replace('.fits',bpmsfx+'.pl')
            if not os.path.exists(orgbpm) or clobber:
                check_exist(orgbpm,'w',yes)
                skyreg=iraf.hselect(newbpm,skysec,expr=yes,Stdout=1)
                if skyreg[0].find('[')>=0:
                    iraf.imcopy(newbpm+skyreg[0],orgbpm,verbose=no)
            elif os.path.exists(orgbpm):
                print "Refusing to clobber old BPM %s" % orgbpm

        # Clean up
        iraf.imdel(objbpm,verify=no,go_ahead=yes)
        iraf.imdel(detbpm,verify=no,go_ahead=yes)

    ##################################################

    # Mask expanded images according to their new BPM's, replacing
    # with pixel values from the median image.  Then reperform the
    # Xregister alignment.

    for i in xrange(nfiles):
        expfile=expfiles[i]
        # Mask image using new BPM
        newbpm=get_head(expfile,'BPM')
        newbpmf=mkbpmfits(newbpm,clobber=yes,delete=no)
        apply_bpm(expfile,newbpmf,medout)
        
    # Delete the median image if it's a temporary one
    if len(medsfx)==0 and not medonly:
        iraf.imdel(medout,verify=no,go_ahead=yes)
        if dobpm:
            iraf.imdel(medweight,verify=no,go_ahead=yes)

    # Second pass of Xregister - no asinh!
    if verbose:
        print ("Second pass of Xregister over region %s " +
               "of expanded images") % xcreg

    # Construct list of images
    xrgin=iraf.mktemp("iqxrgin")+".lis"
    putlines(xrgin,expfiles2)
    check_exist(xshfile,'w',clobber)

    # Run Xregister
    try:
        # this parameter deleted in more recent IRAF versions
        iraf.xregister.databasefmt=no
    except:
        pass

    iraf.xregister('@'+xrgin,refnew,xcreg,xshfile,output="",
                   coords="",xlag=0,ylag=0,
                   dxlag=0,dylag=0,background="none",
                   apodize=0.0,filter="none",correlation="fourier",
                   xwindow=window,ywindow=window,function="centroid",
                   interactive=no,verbose=verbose)

    # Clean up
    os.remove(xrgin)

    # Order "xshifts" the same as expfiles
    xshifts=copy.deepcopy(shifts)
    xshifts[expfiles.index(refnew)]=[0.0,0.0]
    # Read in the xregister shifts
    # Format of the xregister shifts file:
    # <filename> <xshift> <yshift>
    xshlines=getlines(xshfile)
    for xline in xshlines:
        xels=xline.split()
        ix=expfiles.index(xels[0])
        xshifts[ix]=[float(xels[1]),float(xels[2])]
    # Delete the Xregister shifts file, if temporary
    if len(xshext)==0:
        os.remove(xshfile)

    ##################################################

    # Perform integer portion of Xregister shifts
    xfiles=pfx_list(expfiles,xrgpfx,check='w',clobber=clobber)
    for i in xrange(nfiles):
        expfile=expfiles[i]
        xfile=xfiles[i]
        # The full xregister shifts
        xx,xy=xshifts[i]
        # Integer portion
        nx,ny=round(xx),round(xy)
        if nx!=0 or ny!=0:
            # Integral shift of expanded image
            shift_image(expfile,xfile,[nx,ny],
                        border=0,bpmkey="",skysec=skysec)
            # Replace former expfile with this one
            iraf.imdel(expfile,verify=no,go_ahead=yes)
            iraf.imrename(xfile,expfile)
            # Integral shift of image-specific BPM
            newbpm=get_head(expfile,'BPM')
            tmpbpmf=iraf.mktemp("iqcoaddbpmf")+'.fits'
            newbpmf=mkbpmfits(newbpm,clobber=yes,delete=no)
            shift_image(newbpmf,tmpbpmf,[nx,ny],
                        border=1,bpmkey="",skysec=skysec)
            # Replace former BPM with this one
            iraf.imdel(newbpmf,verify=no,go_ahead=yes)
            iraf.imdel(newbpm,verify=no,go_ahead=yes)
            iraf.imcopy(tmpbpmf,newbpm,verbose=no)
            iraf.imdel(tmpbpmf,verify=no,go_ahead=yes)
            # Remove integer portion of xshifts entry
            xshifts[i]=[xx-nx,xy-ny]

    # Perform sub-pixel portion of Xregister shifts
    if subpixel and not medonly:

        # Shift images by subpixel amounts to Xregister positions,
        # and blur their old BPMs to make new BPMs
        for i in xrange(nfiles):
            expfile=expfiles[i]
            # Xregister shift "remainders" dx,dy (integer portion removed)
            dx,dy=xshifts[i]
            # Round shift remainders according to the overresolution setting
            dx=round(dx*overres)/float(overres)
            dy=round(dy*overres)/float(overres)
            # Perform the remainder shift by linear interpolation
            check_exist(xfiles[i],'w',clobber)
            iraf.imshift(expfile,xfiles[i],dx,dy,shifts_file="",
                         interp_type="linear",boundary_type="constant",
                         constant=0)

            # Associated BPM
            oldbpm=get_head(expfile,'BPM')
            newbpm=xrgpfx+oldbpm

            # Blur the BPM to account for interpolation
            xyblur=[]

            if dx != 0:
                xyblur.append([sign(dx),0])
            if dy != 0:
                xyblur.append([0,sign(dy)])
            if len(xyblur)==2:
                xyblur.append([sign(dx),sign(dy)])

            if len(xyblur)>0:
                # Convert to FITS format
                oldbpmf=mkbpmfits(oldbpm,clobber=yes,delete=no)
                newbpmf=xrgpfx+oldbpmf
                # Blur old BPM to make new BPM
                blur_bpm(oldbpmf,newbpmf,xy=xyblur,clobber=clobber,
                         verbose=verbose)
                # Clean up
                iraf.imdel(oldbpmf,verify=no,go_ahead=yes)
                # Convert to PL
                newbpm=mkbpmpl(newbpmf,clobber=yes,delete=yes)
            else:
                # Same as the previous one
                check_exist(newbpm,'w',clobber)
                iraf.imcopy(oldbpm,newbpm)

            # Update BPM filename in header
            update_head(xfiles[i],'BPM',newbpm)

    # Update information in the Xshifts file
    if len(xshext)!=0:
        os.remove(xshfile)
        xshlines=[]
        for i in xrange(nfiles):
            dx,dy=xshifts[i]
            xshlines.append("%s %+.3f %+.3f" % (expfiles[i],dx,dy))
        # Write the new Xshifts file
        putlines(xshfile,xshlines,clobber=yes)

    ##############################

    # Prepare for final combination
    check_exist(output,'w',clobber=clobber)
    
    # Create IRAF listfile
    fininp=iraf.mktemp('iqcoaddfi')+'.lis'
    if subpixel and not medonly:
        putlines(fininp,xfiles)
    else:
        putlines(fininp,expfiles)

    # Rejected pixel mask name
    outrej=iraf.mktemp('iqcoaddrj')+'.pl'
    check_exist(outrej,'w',clobber)

    # Set up for final combination
    if medonly:
        # Second-pass median combination
        imcombine.combine="median"
        imcombine.reject="none"
    else:
        # Mean combination with avsigclip rejection
        imcombine.combine="average"
        imcombine.reject="avsigclip"
    
    # Final image combination
    iraf.imcombine('@'+fininp,output,rejmasks=outrej,
                   project=no,offsets="none",masktype="goodvalue",
                   maskvalue=0,scale="exposure",zero="none",
                   weight="exposure",statsec="",expname=expkey,
                   lthresh="INDEF",hthresh="INDEF",nkeep=1,
                   mclip=yes,lsigma=3.5,hsigma=3.5)

    # Convert rejected-pixel masks to weight image
    wtimage=output.replace('.fits',wtsfx+'.fits')
    check_exist(wtimage,'w',clobber=clobber)
    rejlines=[]
    for i in xrange(nfiles):
        rejlines.append(outrej+'[*,*,%d]' % (i+1))
    rejlis=iraf.mktemp('iqcoaddrj')+'.lis'
    putlines(rejlis,rejlines)
    iraf.imsum('@'+rejlis,wtimage,title="",hparams="",
               pixtype="",calctype="",option="sum",
               low_reject=0,high_reject=0,verbose=no)
    iraf.imarith(nfiles,'-',wtimage,wtimage,noact=no,verbose=no)
    os.remove(rejlis)
    iraf.imdel(outrej,verify=no,go_ahead=yes)

    # Apply average WCS of consituent images to the new image
    wcsaverage('@'+fininp,output)
    
    # Delete IRAF listfile
    os.remove(fininp)

    # Trim final image(s) if requested
    if trimout:
        if trimmax:
            # Trim to region of full overlap of all files
            trimreg=fullreg
        else:
            # Trim to single-image size at reference image location
            trimreg="[%d:%d,%d:%d]" % (refxy[0]+1,refxy[0]+nsmlx,
                                       refxy[1]+1,refxy[1]+nsmly)
        # Final image
        trim=iraf.mktemp("iqcoaddt")+".fits"
        iraf.imcopy(output+trimreg,trim)
        iraf.imdel(output,verify=no,go_ahead=yes)
        iraf.imrename(trim,output,verbose=yes)
        # Weight image
        if not medonly:
            wtrim=iraf.mktemp("iqcoaddw")+".fits"
            iraf.imcopy(wtimage+trimreg,wtrim)
            iraf.imdel(wtimage,verify=no,go_ahead=yes)
            iraf.imrename(wtrim,wtimage,verbose=yes)

    # Clean up intermediate files if requested
    if doclean:
        os.remove(xshfile)
        for image in expfiles:
            iraf.imdel(image,verify=no,go_ahead=yes)
            tmpbpm=glob.glob(image.replace('.fits',bpmsfx)+'.*')
            for bpm in tmpbpm:
                iraf.imdel(bpm,verify=no,go_ahead=yes)
            tmpcrm=glob.glob(image.replace('.fits',crmsfx)+'.*')
            for crm in tmpcrm:
                iraf.imdel(crm,verify=no,go_ahead=yes)
        for image in xfiles:
            if os.path.exists(image):
                iraf.imdel(image,verify=no,go_ahead=yes)
                tmpbpm=glob.glob(image.replace('.fits',bpmsfx)+'.*')
                for bpm in tmpbpm:
                    iraf.imdel(bpm,verify=no,go_ahead=yes)

    # The End
    return

######################################################################

# location of parameter files

_parfile=pardir + "iqmosaic.par"
t=iraf.IrafTaskFactory(taskname="iqmosaic",
                       value=_parfile,function=iqmosaic)

_parfile=pardir + "iqcals.par"
t=iraf.IrafTaskFactory(taskname="iqcals",
                       value=_parfile,function=iqcals)

_parfile=pardir + "iqflatten.par"
t=iraf.IrafTaskFactory(taskname="iqflatten",
                       value=_parfile,function=iqflatten)

_parfile=pardir + "iqsubsky.par"
t=iraf.IrafTaskFactory(taskname="iqsubsky",
                       value=_parfile,function=iqsubsky)

_parfile=pardir + "iqborder.par"
t=iraf.IrafTaskFactory(taskname="iqborder",
                       value=_parfile,function=iqborder)

_parfile=pardir + "iqmask.par"
t=iraf.IrafTaskFactory(taskname="iqmask",
                       value=_parfile,function=iqmask)

_parfile=pardir + "iqclip.par"
t=iraf.IrafTaskFactory(taskname="iqclip",
                       value=_parfile,function=iqclip)

_parfile=pardir + "iqfringe.par"
t=iraf.IrafTaskFactory(taskname="iqfringe",
                       value=_parfile,function=iqfringe)

_parfile=pardir + "iqdefringe.par"
t=iraf.IrafTaskFactory(taskname="iqdefringe",
                       value=_parfile,function=iqdefringe)

_parfile=pardir + "iqobjs.par"
t=iraf.IrafTaskFactory(taskname="iqobjs",
                       value=_parfile,function=iqobjs)

_parfile=pardir + "iqseeing.par"
t=iraf.IrafTaskFactory(taskname="iqseeing",
                       value=_parfile,function=iqseeing)

_parfile=pardir + "iqfocus.par"
t=iraf.IrafTaskFactory(taskname="iqfocus",
                       value=_parfile,function=iqfocus)

_parfile=pardir + "iqwcs.par"
t=iraf.IrafTaskFactory(taskname="iqwcs",
                       value=_parfile,function=iqwcs)

_parfile=pardir + "iqzeropt.par"
t=iraf.IrafTaskFactory(taskname="iqzeropt",
                       value=_parfile,function=iqzeropt)

_parfile=pardir + "iqcrzap.par"
t=iraf.IrafTaskFactory(taskname="iqcrzap",
                       value=_parfile,function=iqcrzap)

_parfile=pardir + "iqcoadd.par"
t=iraf.IrafTaskFactory(taskname="iqcoadd",
                       value=_parfile,function=iqcoadd)

