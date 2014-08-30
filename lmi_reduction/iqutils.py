
# Global Packages
import pyraf
from pyraf import iraf
import copy, os, shutil, glob, sys, string, re, math, operator
import pyfits
from types import *

from numpy import *
#import convolve
from scipy.signal import correlate2d

# ***Must do this after numpy, else 'copy' gets rewritten***
import copy, random

# Global variables
yes=iraf.yes
no=iraf.no
INDEF=iraf.INDEF
# Default behaviors for all functions
globclob=yes
globver=yes

######################################################################

class astrocoords:

    def __init__(self,ra,dec,equinox=2000.0,system="fk5"):

        self.set(ra,dec,equinox=equinox,system=system)

    def set(self,ra,dec,equinox=2000.0,system="fk5"):
            
        # Parse RA
        if type(ra)==StringType:
            ra_reg1=re.search("(\d+)[:\s](\d+)[:\s]([\d\.]+)",ra)
            ra_reg2=re.search("([\d\.]+)",ra)

            if ra_reg1:
                radeg = 15*(float(ra_reg1.group(1)) + \
                            (float(ra_reg1.group(2)) + \
                             float(ra_reg1.group(3))/60.0)/60.0)
            elif ra_reg2:
                radeg = float(ra_reg2.group(1))
            else:
                raise ValueError, "Failed to parse input RA-string \"%s\" as RA"

        elif type(ra)==FloatType:
            radeg = ra

        # Check range limits on RA
        if radeg<0 or radeg>360:
            raise ValueError, "Parsed RA as illegal value %.5f" % radeg

        # Parse Dec
        if type(dec)==StringType:
            dc_reg1=re.search("([\+\-]?\d+)[:\s](\d+)[:\s]([\d\.]+)",dec)
            dc_reg2=re.search("([\+\-]?[\d\.]+)",dec)

            if dc_reg1:
                sgndec=+1
                if dec[0]=='-':
                    sgndec=-1
                absdcd=abs(float(dc_reg1.group(1)))
                dcdeg = sgndec*(absdcd + \
                                (float(dc_reg1.group(2)) + \
                                 float(dc_reg1.group(3))/60.0)/60.0)
            elif dc_reg2:
                dcdeg = float(dc_reg2.group(1))
            else:
                raise ValueError, "Failed to parse input Dec-string \"%s\" as Dec"

        elif type(dec)==FloatType:
            dcdeg = dec

        # Check range limits on Dec
        if abs(dcdeg)>90:
            raise ValueError, "Parsed Dec as illegal value %.5f" % dcdeg

        # Grab the equinox
        xeqx=float(equinox)

        # Check range limits on Equinox
        if xeqx<0 or xeqx>3000:
            raise ValueError, "Illegal Equinox value %.5f" % xeqx

        # Save everything
        self.__radeg=radeg
        self.__dcdeg=dcdeg

        # Direct properties
        self.equinox=xeqx
        self.system=system

    def sxg(self):
        radeg=self.__radeg
        absdec=abs(self.__dcdeg)
        sgndec=cmp(self.__dcdeg,0)
        
        rah=int(radeg/15.0)
        ram=int(60*(radeg/15.0-rah))
        ras=60*(60*(radeg/15.0-rah)-ram)

        dcsgn="+"
        if sgndec<0:
            dcsgn="-"

        dcd=int(absdec)
        dcm=int(60*(absdec-abs(dcd)))
        dcs=60*(60*(absdec-abs(dcd))-dcm)

        return ["%02i:%02i:%07.4f" % (rah,ram,ras),
                dcsgn+"%i:%02i:%06.3f" % (dcd,dcm,dcs)]

    def deg(self):
        return [self.__radeg,self.__dcdeg]

    def rasxg(self):
        rasxg,dcsxg=self.sxg()
        return rasxg

    def dcsxg(self):
        rasxg,dcsxg=self.sxg()
        return self.dcsxg

    def radeg(self):
        return self.__radeg

    def dcdeg(self):
        return self.__dcdeg

    def set_radeg(self,radeg):
        try:
            xradeg=float(radeg)
        except:
            print "radeg must be floating-point RA in degrees"
            return
        if xradeg<0 or xradge>=360:
            print "radeg must be floating-point RA in degrees"
            return
        self.__radeg=xradeg

    def set_dcdeg(self,dcdeg):
        try:
            xdcdeg=float(dcdeg)
        except:
            print "dcdeg must be floating-point Declination in degrees"
            return
        if xdcdeg<-90 or xdcdge>90:
            print "dcdeg must be floating-point Declination in degrees"
            return
        self.__dcdeg=xdcdeg

    def set_deg(self,radeg,dcdeg):
        self.set_radeg(radeg)
        self.set_dcdeg(dcdeg)

    def rad(self):
        dtor=180.0/math.pi
        rarad=self.__radeg/dtor
        dcrad=self.__dcdeg/dtor
        return [rarad,dcrad]

    def shift(self,shift,arcsec=0,arcmin=0):
        dtor=180.0/math.pi

        if arcsec:
            # Shift is quoted in arcsec
            ashift=[shift[0]/3600.0,shift[1]/3600.0]
        elif arcmin:
            # Shift is quoted in arcmin
            ashift=[shift[0]/60.0,shift[1]/60.0]
        else:
            # Shift is quoted in degrees (default)
            ashift=shift

        # New RA
        if abs(self.__dcdeg)<90:
            radeg=self.__radeg+ashift[0]/math.cos(self.__dcdeg/dtor)
        else:
            # Degenerate case
            radeg=self.__radeg

        # Correct for RA negative or >360
        while radeg<0:
            radeg += 360
        while radeg>=360:
            radeg -= 360

        # New Dec
        dcdeg=self.__dcdeg+ashift[1]
        if abs(dcdeg)>90:
            # Oops, over the pole
            sgndec=cmp(dcdeg,0)
            dcdeg=sgndec*(180-abs(dcdeg))

        # Save 'em
        self.__radeg=radeg
        self.__dcdeg=dcdeg

    def diff(self,other,arcsec=0,arcmin=0,degree=1):
        [ra1,dc1]=self.deg()
        [ra2,dc2]=other.deg()
        if arcsec>0 or arcmin>0:
            degree=0
        coodiff=radecdiff(ra1,dc1,ra2,dc2,arcsec=arcsec,
                          arcmin=arcmin,degree=degree)
        return coodiff

######################################################################

class Star:

    def __init__(self,
                 objnum=None,xval=None,yval=None,mag=None,magu=None,
                 flag=0,amaj=None,bmin=None,elong=None,fwhm=None,
                 class_star=1.0,theta=0.0,radeg=None,dcdeg=None,
                 amajw=None,bminw=None,fwhmw=None,thetaw=0.0,
                 epoch=2000.0,equinox=2000.0,system="fk5",
                 name="",color="",mags={},magus={},coo=None,
                 pixtype=None):

        # Null properties set to None at start
        self.objnum=None
        self.xval=None
        self.yval=None
        self.mag=None
        self.magu=None
        self.amaj=None
        self.bmin=None
        self.fwhm=None
        self.elong=None
        self.amajw=None
        self.bminw=None
        self.fwhmw=None
        self.coo=None
        self.pixtype=None

        # Non-null properties set to definition or default
        self.flag=int(flag)
        self.class_star=float(class_star)
        self.theta=float(theta)
        self.name=name
        self.color=color
        self.epoch=float(epoch)
        self.thetaw=float(thetaw)

        # Reset null properties according to the object definition
        if objnum:   self.objnum=int(objnum)
        if xval:     self.xval=float(xval)
        if yval:     self.yval=float(yval)
        if mag:      self.mag=float(mag)
        if magu:     self.magu=float(magu)
        if amaj:     self.amaj=float(amaj)
        if bmin:     self.bmin=float(bmin)
        if fwhm:     self.fwhm=float(fwhm)
        if elong:    self.elong=float(elong)
        if amajw:    self.amajw=float(amajw)
        if bminw:    self.bminw=float(bminw)
        if fwhmw:    self.fwhmw=float(fwhmw)

        # Pixel type property
        if pixtype:
            self.pixtype=pixtype
        elif self.xval or self.yval:
            self.pixtype="image"

        # Coordinates object
        if coo:
            self.coo=coo
        elif radeg and dcdeg:
            self.coo=astrocoords(radeg,dcdeg,equinox=equinox,
                                 system=system)

        # Dictionary of magnitudes
        if len(mags)>0 and type(mags) is DictType:
            self.mags=mags
        elif self.mag:
            self.mags={'MAG':self.mag}
        else:
            self.mags={}

        # Dictionary of magnitude uncertainties
        if len(magus)>0 and type(magus) is DictType:
            self.magus=magus
        elif self.magu:
            self.magus={'MAGU':self.magu}
        else:
            self.magus={}

    ####################

    def set_coo(self,radeg=None,dcdeg=None,equinox=2000.0,
                system="fk5",coo=None):

        if coo:
            self.coo=coo
        elif radeg and dcdeg:
            self.coo=astrocoords(radeg,dcdeg,equinox=equinox,
                                 system=system)

    ####################

    def deg(self):

        if self.coo is None:
            return None
        else:
            return self.coo.deg()

    def radeg(self):

        if self.coo is None:
            return None
        else:
            return self.coo.radeg()

    def dcdeg(self):

        if self.coo is None:
            return None
        else:
            return self.coo.dcdeg()

    def sxg(self):

        if self.coo is None:
            return None
        else:
            return self.coo.sxg()

    def rasxg(self):

        if self.coo is None:
            return None
        else:
            return self.coo.rasxg()

    def dcsxg(self):

        if self.coo is None:
            return None
        else:
            return self.coo.dcsxg()

    def system(self):

        if self.coo is None:
            return None
        else:
            return self.coo.system

    def equinox(self):

        if self.coo is None:
            return None
        else:
            return self.coo.equinox

    ####################

    def haswcs(self):

        haswcs=not (self.coo is None)
        return haswcs

    def haspix(self):

        haspix=(self.xval and self.yval)
        return haspix

    ####################

    def isellipse(self):

        if self.haswcs():
            isellipse=(self.amajw and self.bminw)
        else:
            isellipse=(self.amaj and self.bmin)

        return isellipse

    ####################

    def lineout(self,cols="",format="sex",usepix=0):

        """ Format star as a string with output columns as specified """

        # Parse column specification for cat/sex output
        if cols != "":
            if type(cols)==ListType:
                xcol=cols
            elif type(cols)==StringType:
                xcol=cols.split()
        elif format=="sex":
            # Default columns for sextractor output
            xcol=['objnum','xval','yval','mag','flag','amaj','bmin',
                  'elong','fwhm','class_star','theta']
        elif format=="cat":
            # Default columns for catalog output
            xcol=['name','radeg','dcdeg']
            for mkey in self.mags.keys():
                xcol.append(mkey)

        # Format as line in ds9-style region file
        if format=="reg":

            # Coordinate symbol
            if not usepix and self.haswcs():
                outstr="%s;" % self.system()
                rasxg,dcsxg=self.coo.sxg()

                if self.isellipse():
                    outstr += "ellipse(%s,%s,%0.3f\",%0.3f\",%0.1f)" % \
                              (rasxg,dcsxg,self.amajw,
                               self.bminw,self.thetaw)
                else:
                    outstr += "circle(%s,%s,%0.3f\")" % \
                              (rasxg,dcsxg,self.fwhmw)

            else:
                outstr="%s;" % self.pixtype

                if self.isellipse():
                    outstr += "ellipse(%0.3f,%0.3f,%0.3f,%0.3f,%0.1f)" % \
                              (self.xval,self.yval,self.amaj,
                               self.bmin,self.theta)
                else:
                    outstr += "circle(%0.3f,%0.3f,%0.3f)" % \
                              (self.xval,self.yval,self.fwhm)

            # Comments / magnitudes
            outstr += " #"

            # Label
            if len(self.name)>0:
                outstr += " text={%s}" % self.name

            # Color
            if len(self.color)>0:
                outstr += " color={%s}" % self.color

            # Magnitudes
            if self.mag:
                outstr += " MAG={%0.3f}" % self.mag
            if self.magu:
                outstr += " MAGU={%0.3f}" % self.magu

            # Magnitudes in multiple filters
            if len(self.mags)>0:
                for magnm in self.mags.keys():
                    outstr += " %s={%0.3f}" % (magnm,self.mags[magnm])
                    magunm=magnm+"U"
                    if magunm in self.magus.keys():
                        outstr += " %s={%0.3f}" % \
                                  (magunm,self.magus[magunm])

        # Format as line in Sextractor/Catalog output file
        elif format=="sex" or format=="cat":

            # String to hold our answer
            outstr=""

            # xcol holds our formatting specification
            for col in xcol:
                ckey=col.lower()
                if ckey=="dum":
                    outstr += "0.0"
                else:
                    # Star properties
                    if self.__dict__.has_key(ckey):
                        # Express property as a string
                        outstr += str(self.__dict__[ckey])
                    # Star methods
                    elif ckey=="radeg":
                        outstr += str(self.radeg())
                    elif ckey=="rasxg":
                        outstr += self.rasxg()
                    elif ckey=="dcdeg":
                        outstr += str(self.dcdeg())
                    elif ckey=="dcsxg":
                        outstr += self.dcsxg()
                    elif ckey=="system":
                        outstr += self.system()
                    elif ckey=="equinox":
                        outstr += self.equinox()
                    elif self.mags.has_key(col):
                        outstr += "%0.3f" % self.mags[col]
                    else:
                        print "Failed to parse column specification '%s'" %\
                              ckey
                        outstr += "0.0"
                outstr += " "

            # Trim last space off
            outstr = outstr[:-1]

        # The End
        return outstr

######################################################################
######################################################################

class Starlist:

    def __init__(self,stars=[],comments=[],file=None):

        if type(stars) is ListType:
            self.stars=stars
        elif type(stars) is StringType:
            file=stars
            
        self.comments=comments

        # Region file type codes
        self.regtype="reg"
        self.cattype="cat"
        self.sextype="sex"
        self.magtype="mag"
        self.xytype="xy"

        # Read in the starlist if provided
        if file:
            self.read(file)

    ####################

    def flags(self):

        flags=[]
        nstars=len(self.stars)
        for i in xrange(nstars):
            flags.append(self.stars[i].flag)

        return flags

    def names(self):

        names=[]
        nstars=len(self.stars)
        for i in xrange(nstars):
            names.append(self.stars[i].name)

        return names

    def xvals(self):

        xvals=[]
        nstars=len(self.stars)
        for i in xrange(nstars):
            xvals.append(self.stars[i].xval)

        return xvals

    def yvals(self):

        yvals=[]
        nstars=len(self.stars)
        for i in xrange(nstars):
            yvals.append(self.stars[i].yval)

        return yvals

    def fwhms(self):

        fwhms=[]
        nstars=len(self.stars)
        for i in xrange(nstars):
            fwhms.append(self.stars[i].fwhm)

        return fwhms

    def radegs(self):

        radegs=[]
        nstars=len(self.stars)
        for i in xrange(nstars):
            radegs.append(self.stars[i].radeg())

        return radegs

    def dcdegs(self):

        dcdegs=[]
        nstars=len(self.stars)
        for i in xrange(nstars):
            dcdegs.append(self.stars[i].dcdeg())

        return dcdegs

    def fwhmws(self):

        fwhmws=[]
        nstars=len(self.stars)
        for i in xrange(nstars):
            fwhmws.append(self.stars[i].fwhmw)

        return fwhmws

    def class_stars(self):

        class_stars=[]
        nstars=len(self.stars)
        for i in xrange(nstars):
            class_stars.append(self.stars[i].class_star)

        return class_stars

    def elongs(self):
    
        elongs=[]
        nstars=len(self.stars)
        for i in xrange(nstars):
            elongs.append(self.stars[i].elong)

        return elongs

    ####################

    def set_fwhmws(self,image,pix=None):

        if type(pix)==NoneType:
            pix=impix(image)
            if type(pix)==NoneType:
                print "Can't make conversion without a pix= setting"
                return

        nstars=len(self.stars)
        for i in xrange(nstars):
            self.stars[i].fwhmw=pix*self.stars[i].fwhm

    ####################

    def set_fwhms(self,image,pix=None):

        if type(pix)==NoneType:
            pix=impix(image)
            if type(pix)==NoneType:
                print "Can't make conversion without a pix= setting"
                return

        nstars=len(self.stars)
        for i in xrange(nstars):
            self.stars[i].fwhm=self.stars[i].fwhmw/pix

    ####################

    def mags(self,mag=None):

        mags=[]
        nstars=len(self.stars)
        for i in xrange(nstars):
            if mag:
                if self.stars[i].mags.has_key(mag):
                    mags.append(self.stars[i].mags[mag])
                else:
                    mags.append(None)
            else:
                mags.append(self.stars[i].mag)

        return mags

    def magus(self,mag=None):

        magus=[]
        nstars=len(self.stars)
        for i in xrange(nstars):
            if mag:
                if self.stars[i].magus.has_key(mag+'U'):
                    magus.append(self.stars[i].magus[mag+'U'])
                else:
                    magus.append(None)
            else:
                magus.append(self.stars[i].magu)

        return magus

    ####################

    def mag_sort(self):

        # save a copy of the unsorted starlist
        self.ustars=self.stars

        stgt={}
        for star in self.stars:
            mag=star.mag
            if mag in stgt.keys():
                stgt[mag].append(star)
            else:
                stgt[mag]=[star]

        sms=stgt.keys()
        sms.sort()

        stars2=[]
        for mag in sms:
            for star in stgt[mag]:
                stars2.append(star)

        self.stars=stars2

    def mag_unsort(self):

        self.stars=self.ustars

    ####################

    def unit2asec(self,text):

        """ Parses a string like 3.45" into a number in arcsec " """

        re1=re.search("(.+)(\"|\')$",text,re.I)
        try:
            value=float(re1.group(1))
        except:
            print "Failed to parse %s as number+unit" % text
            return None

        unit=re1.group(2)
        if unit=='"':
            pass
        elif unit=="'":
            value *= 60
        else:
            print "Unrecognized WCS unit --%s--, assuming arcsec" % unit

        return value

    ####################

    def filetype(self,file):

        """ Determine what type of starlist is contained in a file """

        # Read text file into a list
        flines=getlines(file)

        # Test for region file
        if re.search('^#+\s+Region',flines[0],re.I) or \
           re.search('^\s*global',flines[0],re.I) or \
           re.search('^\s*(image|phys|icrs|fk4|fk5|galac|eclip)',flines[0],re.I):
            ftype=self.regtype

        # Test for sextractor starlist
        elif re.search('^#+\s+1',flines[0],re.I):
            ftype=self.sextype

        # Test for ubcone-type catalog file
        elif re.search('^#+FIELD',flines[0],re.I):
            ftype=self.cattype

        # Test for IRAF Mag Type
        elif re.search('^#K\s+IRAF',flines[0],re.I):
            ftype=self.magtype

        # Default is simply a text list
        else:
            ftype=self.xytype

        return ftype

    ####################

    def read(self,file,reg=1,sex=0):

        """ Read starlist from file """

        # Type the list
        ftype=self.filetype(file)

        # Read text file into a list
        flines=getlines(file)

        # The list of stars
        stars=[]
        comments=[]

        # Parse file according to its type
        
        ##############################

        if ftype==self.regtype:

            # Read in ds9-style region file
            for line in flines:
                if line.startswith('#'):
                    # comments
                    comments.append(line)
                    continue
                if line.startswith('global'):
                    # default region properties
                    comments.append(line)
                    continue
                if re.search('^\s*(icrs|fk4|fk5|galac|eclip)',line,re.I):
                    # WCS coords
                    iswcs=1
                    re1=re.search('^\s*([^;\s]+)\s*;',line,re.I)
                    system=re1.group(1)
                    if re.search('fk4',system,re.I):
                        cooeqnx='1950.00'
                    else:
                        cooeqnx='2000.00'
                elif re.search('^\s*(image|phys)',line,re.I):
                    # Pixel coords
                    iswcs=0
                    re1=re.search('^\s*([^;\s]+)\s*;',line,re.I)
                    pixtype=re1.group(1)
                else:
                    # Unrecognized coord type? Skip it.
                    print "Failed to parse this line:\n  %s" % line
                    continue
                # Symbol type
                re2=re.search('^\s*[^;]+;\s*([^\(\s]+)\s*\(',line,re.I)
                symtype=re2.group(1)
                # Coordinate center
                re3=re.search("%s\(([^,]+),([^,]+)," % symtype,line,re.I)
                xval,yval=re3.group(1),re3.group(2)
                symtype=symtype.lower()
                # Symbol size
                if symtype=="circle":
                    # Get radius
                    isellipse=0
                    re4=re.search("circle\([^,]+,[^,]+,([^,\)]+)\)",line,re.I)
                    fwhm=re4.group(1)
                    if iswcs:
                        fwhmw=self.unit2asec(fwhm)
                elif symtype=="ellipse":
                    # Get major + minor axes + angle
                    isellipse=1
                    re4=re.search("ellipse\([^,]+,[^,]+,([^,]+),([^,]+),"+\
                                  "([^\)]+)\)",line,re.I)
                    amaj,bmin,theta=re4.group(1),re4.group(2),re4.group(3)
                    if iswcs:
                        amajw=self.unit2asec(amaj)
                        bminw=self.unit2asec(bmin)
                        fwhmw=math.sqrt(amajw*bminw)
                        thetaw=theta
                    else:
                        try:
                            fwhm=math.sqrt(float(amaj)*float(bmin))
                        except:
                            fwhm=None
                # Star name
                name=""
                re5a=re.search("text=\{([^\}]+)\}",line,re.I)
                if re5a:
                    name=re5a.group(1)
                # Symbol color
                color=""
                re5b=re.search("color=[\'\"\{](\S+)[\'\"\}]",line,re.I)
                if re5b:
                    color=re5b.group(1)
                # Magnitudes
                mag=None; mags={}
                if re.search('MAG=',line,re.I):
                    lineup=line.upper()
                    mgroups=re.findall("\S*MAG=\{[^\}]+\}",lineup)
                    for magtxt in mgroups:
                        re6=re.search("(\S*MAG)=\{([^\}]+)\}",magtxt)
                        if re6:
                            mkey,mval=re6.group(1),re6.group(2)
                            try:
                                mags[mkey]=float(mval)
                                if not mag:
                                    mag=float(mval)
                            except:
                                mags[mkey]=None
                # Magnitude Uncertainties
                magu=None; magus={}
                if re.search('MAGU=',line,re.I):
                    lineup=line.upper()
                    mgroups=re.findall("\S*MAGU=\{[^\}]+\}",lineup)
                    magukey=mgroups[0]
                    for magtxt in mgroups:
                        re7=re.search("(\S*MAGU)=\{([^\}]+)\}",magtxt)
                        if re7:
                            mkey,mval=re7.group(1),re7.group(2)
                            try:
                                magus[mkey]=float(mval)
                                if not magu:
                                    magu=float(mval)
                            except:
                                magus[key]=None

                # Define star
                if iswcs:
                    coords=astrocoords(xval,yval,system=system,equinox=cooeqnx)
                    if isellipse:
                        star=Star(coo=coords,fwhmw=fwhmw,name=name,color=color,
                                  amajw=amajw,bminw=bminw,thetaw=thetaw,
                                  mag=mag,mags=mags,magu=magu,magus=magus)
                    else:
                        star=Star(coo=coords,fwhmw=fwhmw,name=name,color=color,
                                  mag=mag,mags=mags,magu=magu,magus=magus)
                else:
                    if isellipse:
                        star=Star(xval=xval,yval=yval,pixtype=pixtype,
                                  name=name,color=color,fwhm=fwhm,amaj=amaj,
                                  bmin=bmin,theta=theta,mag=mag,mags=mags,
                                  magu=magu,magus=magus)
                    else:
                        star=Star(xval=xval,yval=yval,pixtype=pixtype,
                                  name=name,color=color,fwhm=fwhm,mag=mag,
                                  mags=mags,magu=magu,magus=magus)

                # Add to the list
                stars.append(star)

        ##############################

        elif ftype==self.sextype:

            # Read in sextractor starlist

            # Parse header
            column={}
            for line in flines:
                if line.startswith('#'):
                    re1=re.search("^#\s+(\d+)\s+(\w+)",line)
                    if re1:
                        cnum,cname=re1.groups()
                        column[cname]=int(cnum)-1
                    else:
                        comments.append(line)
                else:
                    # Found the first non-comment line
                    break

            # Out of header
            ncol=len(column)
            idata=flines.index(line)
            nstars=len(flines[idata:]) # Just an estimate
            ndigit=1+int(log(nstars)/log(10))
            nmfmt="SEX_%%0%dd" % ndigit

            # Set all properties to null/default values
            objnum,xval,yval,mag,magu=None,None,None,None,None
            amaj,bmin,elong,fwhm=None,None,None,None
            flag=0
            theta=0.0
            class_star=1.0

            # Pick a preferred magnitude (if we have options)
            if column.has_key('MAG_AUTO'):
                magkey='MAG_AUTO'
                magukey='MAGERR_AUTO'
            elif column.has_key('MAG_APER'):
                magkey='MAG_APER'
                magukey='MAGERR_APER'
            elif column.has_key('MAG_ISOCOR'):
                magkey='MAG_ISOCOR'
                magukey='MAGERR_ISOCORR'
            elif column.has_key('MAG_BEST'):
                magkey='MAG_BEST'
                magukey='MAGERR_BEST'
            elif column.has_key('MAG_ISO'):
                magkey='MAG_ISO' 
                magukey='MAGERR_ISO'
            else:
                magkey=None
                magukey=None

            # Attempt to collect all magnitude and error measurements
            mkeys=[]
            for key in column.keys():
                if re.search('^MAG_',key,re.I):
                    mkeys.append(key)
            mukeys=[]
            for key in column.keys():
                if re.search('^MAGERR_',key,re.I):
                    mukeys.append(key)

            # Parse object data
            didwarn=0
            for line in flines[idata:]:
                # Exclude comments
                if line.startswith('#'):
                    comments.append(line)
                    continue
                # Grab line elements
                a=line.split()
                # Check for legal line length
                if len(a)!=ncol:
                    if not didwarn:
                        print ("Bad number of columns in at least one line "+\
                               "of %s\n") % file
                        didwarn=1
                    continue

                # Grab each property of interest
                if column.has_key('NUMBER'):
                    objnum=a[column['NUMBER']]
                elif column.has_key('OBJNUM'):
                    objnum=a[column['OBJNUM']]
                if column.has_key('X_IMAGE'):
                    xval=a[column['X_IMAGE']]
                elif column.has_key('XVAL'):
                    xval=a[column['XVAL']]
                if column.has_key('Y_IMAGE'):
                    yval=a[column['Y_IMAGE']]
                elif column.has_key('YVAL'):
                    yval=a[column['YVAL']]
                if column.has_key('FLAGS'):
                    flag=a[column['FLAGS']]
                elif column.has_key('FLAG'):
                    flag=a[column['FLAG']]
                if column.has_key('A_IMAGE'):
                    amaj=a[column['A_IMAGE']]
                if column.has_key('B_IMAGE'):
                    bmin=a[column['B_IMAGE']]
                if column.has_key('THETA_IMAGE'):
                    theta=a[column['THETA_IMAGE']]
                if column.has_key('ELONGATION'):
                    elong=a[column['ELONGATION']]
                if column.has_key('FWHM_IMAGE'):
                    fwhm=a[column['FWHM_IMAGE']]
                if column.has_key('CLASS_STAR'):
                    class_star=a[column['CLASS_STAR']]

                # Primary magnitude and uncertainty
                if magkey:
                    mag=a[column[magkey]]
                if column.has_key(magukey):
                    magu=a[column[magukey]]

                # Full set of magnitudes
                mags={}
                for mkey in mkeys:
                    try:
                        mags[mkey]=float(a[column[mkey]])
                    except:
                        mags[mkey]=None
                magus={}
                for mukey in mukeys:
                    try:
                        magus[mukey]=float(a[column[mukey]])
                    except:
                        magus[mukey]=None

                # Construct object name from its number
                name=nmfmt % int(objnum)

                # Add to the list
                stars.append(Star(objnum=objnum,xval=xval,yval=yval,mag=mag,
                                  mags=mags,magu=magu,magus=magus,flag=flag,
                                  amaj=amaj,bmin=bmin,theta=theta,elong=elong,
                                  fwhm=fwhm,class_star=class_star,
                                  pixtype="image",name=name))

        ##############################

        elif ftype==self.cattype:

            # Need a default "size" for catalog objects
            def_fwhmw=4.0
            
            # Header of the catalog file includes epoch & equinox
            # information
            system="icrs"
            equinox,epoch=2000.0,2000.0

            # Catalog parameters / comments
            for line in flines:
                if line.startswith('#'):
                    comments.append(line)
                    if re.search("^#+EQUINOX",line,re.I):
                        #re1=re.search("=\s*(\w)([\d\.]+)\s+([\d\.]+)",
                                      #line,re.I)
                        #equinox,epoch=re1.group(2),re1.group(3)
                        #estart=re1.group(1)
                        re1=re.search("=\s*(\w)([\d\.]+)", line, re.I)
                        equinox = re1.group(2)
                        if re1.group(1).startswith('B'):
                            system="fk4"
                    if re.search("^#+EPOCH", line, re.I):
                        re1=re.search("=\s*([\d\.]+)", line, re.I)
                        epoch = re1.group(1)
                    if re.search("^#+END",line,re.I):
                        break

            # Parse header line
            ihead=1+flines.index(line)
            hline=flines[ihead]
            cols=re.findall("[^\|\s]+\|",flines[ihead])
            units=re.findall("[^\|\s]+\|",flines[1+ihead])
            ncol=len(cols)
            column={}
            cunits={}
            mkeys=[]
            for i in xrange(ncol):
                re1=re.search("([^\|]+)\|",cols[i])
                cname=re1.group(1).upper()
                re2=re.search("([^\|]+)\|",units[i])
                cunit=re2.group(1).upper()
                if cunit=='MAG':
                    cname=cname+'MAG'
                    mkeys.append(cname)
                column[cname]=i
                cunits[cname]=cunit

            # Skip the other header line(s)
            for line in flines[1+ihead:]:
                if line.startswith('#'):
                    pass
                else:
                    break

            # Out of header
            idata=flines.index(line)

            # Parse object data
            didwarn=0
            for line in flines[idata:]:
                # Exclude comments
                if line.startswith('#'):
                    continue
                # Grab line elements
                a=line.split()
                # Check for legal line length
                if len(a)!=ncol:
                    if not didwarn:
                        print ("Bad number of columns in at least one line "+\
                               "of %s\n") % file
                        didwarn=1
                    continue

                # Grab each property of interest
                if column.has_key('ID'):
                    name=a[column['ID']]
                if column.has_key('RA'):
                    radeg=a[column['RA']]
                if column.has_key('DEC'):
                    dcdeg=a[column['DEC']]

                # Have no place to put the coordinate position
                # uncertainties, proper motions, and proper motion
                # uncertainties at this point.
                fwhmw=def_fwhmw

                # Full set of magnitudes
                mags={}
                for mkey in mkeys:
                    try:
                        mags[mkey]=float(a[column[mkey]])
                    except:
                        mags[mkey]=None

                # Add to the list
                stars.append(Star(name=name,radeg=radeg,dcdeg=dcdeg,
                                  fwhmw=fwhmw,mags=mags,system=system,
                                  equinox=equinox,epoch=epoch))
        
        ##############################

        elif ftype==self.magtype:

            # Need to set fwhm arbitrarily at this point
            fwhm=4.0

            # Parse comments
            for line in flines:
                if line.startswith('#'):
                    comments.append(line)
                else:
                    break

            # Determine length of List
            idata=flines.index(line)
            nstars=len(flines[idata:])/5.0
            ndigit=1+int(log(nstars)/log(10))
            nmfmt="IRAF_%%0%dd" % ndigit

            # Grab Star Properties
            i=idata
            while (i < len(flines)):
                
                try:
                    # Set initial values
                    name,xval,yval,cflag,sflag,filter,mag,magu,mflag,flag= \
                     None,None,None,None,None,None,None,None,None,None
                    mags={}
                    magus={}
                    # Set Object Name
                    name=nmfmt % ((i-idata)/5+1)
                    # Grab pixel coordinates                 
                    i+=1
                    line2=flines[i]
                    xval=line2[:14].strip()
                    yval=line2[14:25].strip()
                    cflag=line2[64:69].strip()
                    # Grab sky error
                    i+=1
                    line3=flines[i]
                    sflag=line3[64:69].strip()
                    # Grab filter name
                    i+=1
                    line4=flines[i]
                    filter=line4[33:56].strip()
                    # Grab Magnitude
                    i+=1
                    line5=flines[i]
                    mag=line5[51:58].strip()
                    magu=line5[58:64].strip()
                    mflag=line5[64:69].strip()
                    # Update i
                    i+=1
                    # Update flag
                    flag=float(cflag) or float(sflag) or float(mflag)
                    # Update Mags and stars
                    if not ((mag == 'INDEF') or (magu == 'INDEF')):
                        mag=float(mag)
                        magu=float(magu)
                        mags[filter.upper()+'MAG']=mag
                        magus[filter.upper()+'MAGU']=magu
                        stars.append(Star(objnum=(i-idata)/5,xval=xval,
                                          yval=yval,fwhm=fwhm,mag=mag,magu=magu,
                                          mags=mags,magus=magus,name=name,
                                          flag=flag,pixtype="image"))
                except:
                    print "Error reading line: %s" % file
                    i+=1

        ##############################

        elif ftype==self.xytype:

            # Read in xy list
            for line in flines:
                if line.startswith('#'):
                    comments.append(line)
                    continue
                els=string.split(line)
                if len(els)==2:
                    star=Star(len(starlist),els[0],els[1])
                elif len(els)==3:
                    if re.search('[^\d\.e]',els[2],re.I):
                        star=Star(len(starlist),els[0],els[1],name=els[2])
                    else:
                        star=Star(len(starlist),els[0],els[1],fwhm=els[2])
                elif len(els)==4:
                    if re.search('[^\d\.e]',els[3]):
                        star=Star(len(starlist),els[0],els[1],fwhm=els[2],
                                  name=els[3])
                    else:
                        star=Star(len(starlist),els[0],els[1],fwhm=els[2],
                                  mag=els[3])
                elif len(els)==5:
                    if re.search('[^\d\.e]',els[4]):
                        star=Star(len(starlist),els[0],els[1],fwhm=els[2],
                                  mag=els[3],name=els[4])
                    else:
                        star=Star(len(starlist),els[0],els[1],fwhm=els[2],
                                  mag=els[3],magu=els[4])
                elif len(els)>5:
                    if re.search('[^\d\.e]',els[-1]):
                        star=Star(len(starlist),els[0],els[1],fwhm=els[2],
                                  mag=els[3],magu=els[4],name=els[-1])
                    else:
                        star=Star(len(starlist),els[0],els[1],fwhm=els[2],
                                  mag=els[3],magu=els[4])
                else:
                    sys.stderr.write("Need at least X,Y to define Star\n")
                    continue
                
                stars.append(star)

        # Save the list & comments
        self.comments=comments
        self.stars=stars

    ##############################

    def write(self,file,format="",cols="",maxnum=0,usepix=0,
              clobber=globclob):

        """ Write starlist to file in appropriate format """

        # Default format is same as file extension
        if len(format)==0:
            root,ext=os.path.splitext(file)
            xext=ext[1:]
            if xext==self.regtype or xext==self.sextype or \
               xext==self.cattype or xext==self.xytype:
                format=xext
            else:
                print "Please specify format of output file"
                return

        # Output text as list
        #lines=copy.deepcopy(self.comments)
        lines=[]

        # Parse column specification for cat/sex output
        if cols != "":
            if type(cols)==ListType:
                xcol=cols
            elif type(cols)==StringType:
                xcol=cols.split()
        elif format==self.sextype:
            xcol=['objnum','xval','yval','mag','flag','amaj','bmin',
                  'elong','fwhm','class_star','theta']
        elif format==self.cattype:
            xcol=['name','radeg','dcdeg']
            for mkey in self.stars[0].mags.keys():
                xcol.append(mkey)
        else:
            # Dummy variable
            xcol=[]

        # Write header for sextractor output
        if format==self.sextype:
            for i in xrange(0,len(xcol)):
                lines.append("#  %2d %-16s" % (i+1,xcol[i].upper()))
        # Write header for catalog output
        elif format==self.cattype:
            cline="##"
            uline="##"
            nline="##"
            for i in xrange(0,len(xcol)):
                nline+="  %d|" % (i+1)
                if xcol[i].lower()=='name':
                    cline+='  id|'
                    uline+='  id|'
                elif xcol[i].lower()=='radeg':
                    cline+='       RA|'
                    uline+='  ddd.ddd|'
                elif xcol[i].lower()=='dcdeg':
                    cline+='     DEC|'
                    uline+='  dd.ddd|'
                elif re.search("mag$",xcol[i],re.I):
                    re1=re.search("(\S+)mag",xcol[i],re.I)
                    if re1:
                        cline+='     %s|' % re1.group(1)
                    else:
                        cline+='    MAG|'
                    uline+='   MAG|'
                else:
                    cline+="  %s|" % xcol[i]
                    uline+="  xx|"
            # These are the header lines for the catalog
            lines.append(cline)
            lines.append(uline)
            lines.append(nline)

        # Get a line for each star
        for star in self.stars:
            if maxnum>0 and self.stars.index(star)==maxnum:
                break
            lines.append(star.lineout(cols=xcol,format=format,
                                      usepix=usepix))

        # Write the lines to disk
        check_exist(file,'w',clobber=clobber)
        putlines(file,lines)

    ##############################

    def findshift(self,refstars,window=20.0,pthresh=0.95,
                  tol=3.0,useflags=no,region=None,
                  image=None,maxnum=50,skiptop=None,doshift=no):

        """ Find best (x,y) shift to refstars """

        # Use image boundary if an image is specified
        if region is None and image:
            [naxis1,naxis2]=get_head(image,['NAXIS1','NAXIS2'])
            region=[[1,naxis1],[1,naxis2]]

        # Apply region or no?
        useregion=region and len(region)==2 and \
                   len(region[0])==2 and len(region[1])==2

        # Order Starlists by magnitude (brightest first)
        self.mag_sort()
        refstars.mag_sort()

        # Construct some working arrays
        ax=array(self.xvals())
        ay=array(self.yvals())

        # Skip brightest stars if requested
        nskip=0
        if skiptop:
            try:
                nskip=int(skiptop)
            except:
                print "Failed to parse skiptop as integer"

        # Arrays we want to populate
        dx=[]; dy=[]

        # Loop over reference stars, brightest to faintest
        istar=0
        for refstar in refstars[nskip:]:

            # Exclude stars outside of specified region
            if useregion:
                if refstar.xval+tol<region[0][0] or \
                   refstar.xval-tol>region[0][1] or \
                   refstar.yval+tol<region[1][0] or \
                   refstar.yval-tol>region[1][1]:
                    continue

            # Apply flag exclusion
            if useflags and refstar.flag>0:
                continue

            # Find stars within window of refstar
            ix=where(logical_and(less(abs(ax-refstar.xval),window), 
                                 less(abs(ay-refstar.yval),window)))[0]
            # Have candidates
            if len(ix)>0:
                # Save all dx, dy
                for ii in ix:
                    # Apply flag exclusion
                    if useflags and self.stars[ii].flag>0:
                        continue
                    dx.append(refstar.xval-ax[ii])
                    dy.append(refstar.yval-ay[ii])
                # Check if we're done
                istar += 1
                if maxnum>0 and istar>=maxnum:
                    break

        # Find peak of histogram
        nn=len(dx)
        dx=array(dx); dy=array(dy)
        idx=round(dx); idy=round(dy)
        iw=round(window)
        npeak=0
        for i in xrange(-iw,iw,1):
            for j in xrange(-iw,iw,1):
                nij=len(where(logical_and(equal(idx,i),equal(idy,j)))[0])
                if nij>npeak:
                    npeak=nij
                    ishift=[i,j]

        # Test for significance
        goodshift=1
        xmean=float(nn)/float(4*window*window)
        psig=poiprob(npeak-1,xmean)

        # -- No further action for unsuccessful shifts -- #
        if psig<=pthresh:
            print "Shift probability (%.2f) does not exceed threshold (%.2f)" % \
                  (psig,pthresh)
            goodshift=0
            return None

        # -- Further analysis for successful shifts -- #

        # Centroid the shift
        ix=where(logical_and(less(abs(dx-ishift[0]),tol), 
                             less(abs(dy-ishift[1]),tol)))[0]
        xshift=[mean(dx[ix]),mean(dy[ix])]

        # Apply the shift, if requested
        if doshift:
            for star in self.stars:
                star.xval += xshift[0]
                star.yval += xshift[1]
        
        # The End
        return xshift

    ##############################

    def match(self,refstars,tol=3.0,useflags=yes,region=None,
              image=None,maxnum=50,maxskip=None):

        """ Construct matched lists of stars between two starlists """

        # Order Starlists by magnitude (brightest first)
        self.mag_sort()
        refstars.mag_sort()

        astars=copy.deepcopy(self.stars)
        ax=array(self.xvals())
        ay=array(self.yvals())

        # Use image boundary if an image is specified
        if region is None and not (image is None):
            [naxis1,naxis2]=get_head(image,['NAXIS1','NAXIS2'])
            region=[[1,naxis1],[1,naxis2]]

        # Apply region or no?
        useregion=region and len(region)==2 and \
                   len(region[0])==2 and len(region[1])==2

        nskip=0
        mystars=[]
        matches=[]

        # Loop over reference stars, brightest to faintest
        for refstar in refstars:
            if useregion:
                # Stars outside of region are not considered
                if refstar.xval+tol<region[0][0] or \
                   refstar.xval-tol>region[0][1] or \
                   refstar.yval+tol<region[1][0] or \
                   refstar.yval-tol>region[1][1]:
                    pass
            # Find stars within tolerance of refstar
            ix=where(logical_and(less(abs(ax-refstar.xval),tol), 
                                 less(abs(ay-refstar.yval),tol)))[0]
            # No stars within tolerance of refstar
            if len(ix)==0:
                nskip += 1
                if maxskip and nskip>maxskip:
                    break
                continue
            # Have candidates
            nskip=0
            # Examine the candidates
            if len(ix)==1:
                # One candidate
                mix=ix[0]
                mystar=astars[mix]
            else:
                # Multiple candidates: pick the closest
                dist=hypot(ax[ix]-refstar.xval,ay[ix]-refstar.yval)
                iy=where(equal(dist,min(dist)))[0]
                mix=ix[iy[0]]
                mystar=astars[mix]
            # Delete the matched star from the list
            astars.pop(mix)
            # Numarrays can't be pop'd, so I do this instead
            ax=compress(not_equal(arange(len(ax)),mix),ax)
            ay=compress(not_equal(arange(len(ay)),mix),ay)
            # Check the flags, if requested
            if useflags and (mystar.flag>0 or refstar.flag>0):
                continue
            # Add to the list of matches
            mystars.append(mystar)
            matches.append(copy.deepcopy(refstar))
            # Check if we're done
            if maxnum>0 and len(mystars)==maxnum:
                break

        # Turn the lists of stars into Starlists
        matchlist=Starlist(mystars,
                           comments=copy.deepcopy(self.comments))
        refmatch=Starlist(matches,
                          comments=copy.deepcopy(refstars.comments))
        
        # The End
        return [matchlist,refmatch]

    ##############################
            
    def nomatch(self,refstars,tol=3.0,region=None,
                image=None,maxskip=None):

        """ Construct list of objects NOT in a catalog """

        # Order Starlists by magnitude (brightest first)
        self.mag_sort()
        refstars.mag_sort()

        astars=copy.deepcopy(self.stars)
        ax=array(self.xvals())
        ay=array(self.yvals())

        # Use image boundary if an image is specified
        if region is None and not (image is None):
            [naxis1,naxis2]=get_head(image,['NAXIS1','NAXIS2'])
            region=[[1,naxis1],[1,naxis2]]

        # Apply region or no?
        useregion=region and len(region)==2 and \
                   len(region[0])==2 and len(region[1])==2

        nskip=0
        mystars=[]
        matches=[]

        # Loop over reference stars, brightest to faintest
        for refstar in refstars:
            if useregion:
                # Stars outside of region are not considered
                if refstar.xval+tol<region[0][0] or \
                   refstar.xval-tol>region[0][1] or \
                   refstar.yval+tol<region[1][0] or \
                   refstar.yval-tol>region[1][1]:
                    pass
            # Find stars within tolerance of refstar
            ix=where(logical_and(less(abs(ax-refstar.xval),tol), 
                                 less(abs(ay-refstar.yval),tol)))[0]
            # No stars within tolerance of refstar
            if len(ix)==0:
                nskip += 1
                if maxskip and nskip>maxskip:
                    break
                continue
            # Have candidate matches
            nskip=0
            # Examine the matches
            if len(ix)==1:
                # One candidate
                mix=ix[0]
                mystar=astars[mix]
            else:
                # Multiple candidates: pick the closest
                dist=hypot(ax[ix]-refstar.xval,ay[ix]-refstar.yval)
                iy=where(equal(dist,min(dist)))[0]
                mix=ix[iy[0]]
                mystar=astars[mix]
            # Delete the matched star from the list
            astars.pop(mix)
            ax=compress(not_equal(arange(len(ax)),mix),ax)
            ay=compress(not_equal(arange(len(ay)),mix),ay)
            # Add to the list of matches
            mystars.append(mystar)
            matches.append(copy.deepcopy(refstar))

        # Turn the leftover stars (no matches) into a Starlist
        candlist=Starlist(astars,
                          comments=copy.deepcopy(self.comments))

        # The End
        return candlist
            
    ##############################

    def zeropt(self,refstars,tol=3.0,useflags=yes,region=None,image=None,
               maxnum=50,maxskip=None,method="median",rejout=0,fencelim=0.50,
               sigma=3.0,maxfrac=0.15):

        """ Calculates delta zero-point between two starlists """

        # Do star matching against reference catalog
        matchstars,refmatch = \
            self.match(refstars,tol=tol,useflags=useflags,
                       region=region,image=image,maxnum=maxnum,
                       maxskip=maxskip)

        # Make array of mag differences (delta zero-points)
        matchmags=array(matchstars.mags())
        refmags=array(refmatch.mags())
        dzpts=refmags-matchmags

        # Make sure there are enough stars in match
        if (len(matchstars)<2):
            return 0, 0

        # Swith based on method
        if (method=="mean"):
            # Reject outliers?
            if rejout:
                dzpts=reject_outliers(dzpts,fencelim=fencelim,sigma=sigma,
                                      maxfrac=maxfrac)
            # Calculate mean and standard deviation
            dzpt=dzpts.mean()
            dzptu=dzpts.std()
        
        else:        
            # Calculate median & uncertainty therein
            dzpt=median(dzpts)
            dzptu=medunc(dzpts)

        return dzpt,dzptu

    ####################

    def set_mag(self,mkey):

        """ Set a particular magnitude to be the key magnitude """

        if not self.stars[0].mags.has_key(mkey):
            print "Couldn't find magnitude '%s' in starlist" % mkey
            return

        ukey=mkey+'U'
        dounc=self.stars[0].magus.has_key(ukey)

        nstars=len(self.stars)
        for i in xrange(nstars):
            self.stars[i].mag=self.stars[i].mags[mkey]
            if dounc:
                self.stars[i].magu=self.stars[i].magus[ukey]

    ####################

    def pix2wcs(self,image):

        """ Convert pixel coordinates to WCS using the WCS of IMAGE """

        if not os.path.exists(image):
            print "Requested reference image '%s' does not exist" % image
            return

        nstars=len(self.stars)

        x,y=[],[]
        for i in xrange(nstars):
            x.append(self.stars[i].xval)
            y.append(self.stars[i].yval)

        isphys=0
        if self.stars[0].pixtype=="physical":
            isphys=1
            
        radec=impix2wcs(image,x,y,physical=isphys)

        for i in xrange(nstars):
            self.stars[i].set_coo(radec[i][0],radec[i][1])

        # Try converting FWHM's to FWHMW's
        self.set_fwhmws(image)

    ####################

    def wcs2pix(self,image,physical=0):

        """ Convert WCS to pixel coordinates using the WCS of IMAGE """

        if not os.path.exists(image):
            print "Requested reference image '%s' does not exist" % image
            return

        nstars=len(self.stars)

        ra,dec=[],[]
        for i in xrange(nstars):
            ra.append(self.stars[i].radeg())
            dec.append(self.stars[i].dcdeg())

        xy=imwcs2pix(image,ra,dec,physical=physical)

        if physical:
            pixtype="physical"
        else:
            pixtype="image"

        for i in xrange(nstars):
            self.stars[i].xval=xy[i][0]
            self.stars[i].yval=xy[i][1]
            self.stars[i].pixtype=pixtype

        # Try converting FWHMW's to FWHM's
        self.set_fwhms(image)

    ####################

    def starsinbox(self,box):

        """ Return the starlist for a specified box-shape region """

        outstars=[]

        for star in self.stars:
            [sx,sy]=[star.xval,star.yval]
            if sx>=box[0][0] and sx<=box[0][1] and \
                   sy>=box[1][0] and sy<=box[1][1]:
                outstars.append(copy.deepcopy(star))

        outlist=Starlist(outstars)

        return outlist

    def nstarsinbox(self,box):

        num=len(self.starsinbox(box))
        return num

    ##############################

    # Override default routines to allow the "Starlist" object to
    # function in regular Python like a list of "Stars"

    def __len__(self):
        return len(self.stars)

    def __getitem__(self,ix):
        return self.stars[ix]

    def __setitem__(self,ix,value):
        self.stars[ix]=value

    def __delitem__(self,ix):
        del self.stars[ix]

    def __getslice__(self,ix,jx):
        return self.stars[ix:jx]

    def __setslice__(self,ix,jx,values):
        self.stars[ix:jx]=values

    def __delslice__(self,ix,jx):
        self.stars[ix:jx]=[]

    def __concat__(self,s2):
        self.stars=self.stars.append(s2.stars)

    def __contains__(self,element):
        return (element in self.stars)

######################################################################
######################################################################

def iraffiles(files,nfiles=0):

    if type(files) is not StringType:
        print "Input filelist is not a string"
        exit(1)

    fout=[]
    fmult=files.split(",")
    for fcand in fmult:
        re1=re.search("^(.+//)?@(.+)(//.+)?$",fcand)
        re2=re.search("[\*\?]",fcand)
        if re1:
            # Using the IRAF "@file.lis" convention
            flist=re1.group(2)
            if os.path.exists(flist):
                fflist=getlines(flist)
                for fmem in fflist:
                    if re1.group(1):
                        fmem=re1.group(1)[:-2]+fmem
                    if re1.group(3):
                        fmem=fmem+re1.group(3)[2:]
                    if (fitsfile(fmem)!=""):
                        fout.append(fitsfile(fmem))
        elif re2:
            # Using UNIX wildcards
            flist=glob.glob(fcand)
            for fmem in flist:
                if (fitsfile(fmem)!=""):
                    fout.append(fitsfile(fmem))
        else:
            # Just plain filenames (?)
            if fitsfile(fcand)!="":
                fout.append(fitsfile(fcand))
            
    return fout

######################################################################

def fitsfile(file):
    outfile=""

    if file.endswith('.fits') or file.endswith('.imh'):
        if os.path.exists(file):
            outfile=file
        else:
            print "Can't find requested file '%s'" % file
    else:
        if os.path.exists(file+'.fits'):
            outfile=file+'.fits'
        elif os.path.exists(file):
            outfile=file
        elif os.path.exists(file+'.imh'):
            outfile=file+'.imh'
        else:
            print "Can't find requested file '%s' or variants" % file
        
    return outfile

######################################################################

def check_exist(filename, status, clobber=globclob):

    """ check_exist(filename, status, clobber=yes)
    checks to see if filename exists
    if status==r, must exist, otherwise prints error + returns False
    if status==w, if exists and clobber=no then prints error + returns False
    else deletes + returns True
    """     

    if (status == "r"):
        # check to see if it exists for reading
        # (i.e. must be present)
        if (not (os.path.exists(filename))):
            print "Couldn't open input file: %s" % filename
            return False
    else:
        # check to see if it exists for writing
        # (i.e. must not exist or clobber=yes)
        if (os.path.exists(filename)):
            if (clobber):
                os.remove(filename)
            else:
                print "File %s already exists and clobber=no" % filename
                return False

    return True

######################################################################

def getlines(filename,blanks=no):

    """ lines=getlines(filename,blanks)
    opens file filename
    reads lines and strips off whitespace
    if (blanks==1), will return fully blank lines
    otherwise excises
    """

    infile=open(filename)
    lines=infile.readlines()
    infile.close()
    lines2=[]

    if (blanks==1):
        return lines
    
    for i in range(0,len(lines)):
        if (re.search("\S+",lines[i])):
            lines2.append(lines[i].rstrip())

    return lines2

######################################################################

def putlines(filename,list,clobber=no):

    """ putlines(filename,list)
    writes the string list [list] to the file [filename]
    """

    if os.path.exists(filename):
        check_exist(filename,'w',clobber=clobber)
        
    outfile=open(filename,'w')
    for line in list:
        if type(line) != StringType:
            try:
                lineout=string(line)
            except:
                print "Failed to convert list element to string in putlines"
                return
        else:
            lineout=line
        outfile.write(lineout+"\n")
    outfile.close()

    return

######################################################################

def pfx_list(list, pfx, check='', clobber=yes):

    """ pfx_list(list,pfx):  add pfx to the start of every element in list """

    if type(pfx)!=StringType or type(list[0])!=StringType:
        print "pfx_list works only with string types"
        return list

    out=[]
    for item in list:
        if type(item)!=StringType:
            print "pfx_list works only with string types"
            return list
        if len(check)>0:
            if check=='r':
                check_exist(pfx+item,'r')
            elif check=='w':
                check_exist(pfx+item,'w',clobber=clobber)
        out.append(pfx+item)

    return out

######################################################################

def check_head(file,keys,extn=0):

    check_exist(file,"r")
    out=no

    try:
        fimg=pyfits.open(file)
        if extn>=len(fimg):
            print "Requested extension [%d] does not exist" % extn
            return no
        
        head=fimg[extn].header

        if type(keys)==StringType:
            key=keys
            out=head.has_key(key)
        elif type(keys)==ListType:
            out=yes
            for key in keys:
                if not head.has_key(key):
                    out=no
                    break
        else:
            if verbose:
                print "Bad variable type for keys"

        fimg.close()

    except:
        print "Error reading header of %s" % file
        out=no

    return out

######################################################################

def get_head(file,keys,extn=0,verbose=globver):

    """ reads one or more header keywords from a FITS file
        using PyFITS """

    vals=[]
    check_exist(file,"r")

    try:
        fimg=pyfits.open(file)
        if extn>=len(fimg):
            if verbose:
                print "Requested extension [%d] does not exist" % extn
            vals=[""]*len(keys)
            return vals

        head=fimg[extn].header

        if type(keys)==StringType:
            key=keys
            if head.has_key(key):
                vals=head[key]
            else:
                if verbose:
                    print "Error reading keyword %s from %s" % (key,file)
                vals=""
        elif type(keys)==ListType:
            for key in keys:
                if head.has_key(key):
                    vals.append(head[key])
                else:
                    if verbose:
                        print "Error reading keyword %s from %s" % (key,file)
                    vals.append("")
        else:
            if verbose:
                print "Bad variable type for keys"
        fimg.close()

    except:
        print "Error reading header of %s" % file
        vals=[""]*len(keys)

    return vals

######################################################################

def keysplit(infiles,key):

    dict={}
    for image in infiles:
        val=get_head(image,key,verbose=no)
        if dict.has_key(val):
            dict[val].append(image)
        else:
            dict[val]=[image]

    return dict

######################################################################

def update_head(inlist,inkeys,invals,comments="",extn=0):

    """ update header keywords using pyfits """

    # Images to process
    if type(inlist) is ListType:
        infiles=inlist
    elif type(inlist) is StringType:
        infiles=iraffiles(inlist)
    else:
        print "Please pass a string or list of input files"
        return

    # Input checking
    if type(inkeys)==ListType and type(invals)!=ListType:
        print "Keywords and values must both be lists"
        return

    # Header keywords to update
    if type(inkeys) is StringType:
        keys=[inkeys]
        vals=[invals]
        cmts=[comments]
    elif type(inkeys) is ListType:
        keys=inkeys
        vals=invals
        if type(comments) is ListType:
            cmts=comments
        else:
            cmts=None
    else:
        print "Please pass a string or list of header keywords to update"
        return

    # Loop over files
    nkeys=len(keys)
    for file in infiles:

        check_exist(file,"r")
        try:
            inf=pyfits.open(file,"update")

            if extn>=len(inf):
                print "Requested extension [%d] does not exist" % extn
                inf.close()
                return

            for i in xrange(nkeys):
                if cmts:
                    inf[extn].header.update(keys[i],vals[i],comment=cmts[i])
                else:
                    inf[extn].header.update(keys[i],vals[i])

            inf.close()

        except:
            print "Error updating header of %s" % file
      
######################################################################

def delete_head(inlist,inkeys,extn=0):

    """ delete header keywords using pyfits """

    # Images to process
    if type(inlist) is ListType:
        infiles=inlist
    elif type(inlist) is StringType:
        infiles=iraffiles(inlist)
    else:
        print "Please pass a string or list of input files"
        return

    # Header keywords to delete
    if type(inkeys) is StringType:
        keys=[inkeys]
    elif type(inkeys) is ListType:
        keys=inkeys
    else:
        print "Please pass a string or list of header keywords to delete"
        return

    # Loop over files
    for file in infiles:

        check_exist(file,"r")
        try:
            inf=pyfits.open(file,"update")

            if extn>=len(inf):
                print "Requested extension [%d] does not exist" % extn
                inf.close()
                return

            for key in keys:
                if inf[extn].header.has_key(key):
                    del inf[extn].header[key]

            inf.close()

        except:
            print "Error updating header of %s" % file
      
##############################

def asinh_image(input,output,clobber=globclob,verbose=globver):

    """ takes the arc-sinh of an image.  this is a useful stretch for
        many applications. 

        input   name of input image
        output  name of output image

        clobber clobber output files [yes]
        verbose print messages about actions [yes]
    """

    # Input checking
    if not os.path.exists(input) and not os.path.exists(input+".fits"):
        print "Couldn't open input image %s" % input
        return
    elif os.path.exists(input+".fits"):
        input+=".fits"
    if input==output:
        if clobber:
            print "Input image %s will be overwritten" % input
        else:
            print "Refusing to overwrite input without permission"
            return
    check_exist(output,'w',clobber)

    # Open the input image as pyfits object
    fimg=pyfits.open(input)
    hdr=fimg[0].header
    D=fimg[0].data
    DX=arcsinh(D)
    fimg[0].data=DX
    hdr.update('ARCSINH',1,'Arcsinh of data from %s' % input)

    # Write/close the output image
    fimg.writeto(output)

##############################

def sinh_image(input,output,clobber=globclob,verbose=globver):

    """ takes the sinh of an image, reversing a previous "asinh" 

        input   name of input image
        output  name of output image

        clobber clobber output files [yes]
        verbose print messages about actions [yes]
    """

    # Input checking
    if not os.path.exists(input) and not os.path.exists(input+".fits"):
        print "Couldn't open input image %s" % input
        return
    elif os.path.exists(input+".fits"):
        input+=".fits"
    if input==output:
        if clobber:
            print "Input image %s will be overwritten" % input
        else:
            print "Refusing to overwrite input without permission"
            return
    check_exist(output,'w',clobber)

    # Open the input image as pyfits object
    fimg=pyfits.open(input)
    hdr=fimg[0].header
    D=fimg[0].data
    DX=sinh(D)
    fimg[0].data=DX
    hdr.update('ARCSINH',0,'Undid arcsinh of data')

    # Write/close the output image
    fimg.writeto(output)

##############################

def shift_image(input,output,shift,border=0,bpmkey="BPM",bpmnew="",
                skysec="SKYSEC",clobber=globclob,verbose=globver):

    """ shifts input image.  WCS is preserved.  

        input   name of input image
        output  name of output image
        shift   the shift:  [dx,dy] in pixels
        border  value to set in border regions [0]
        bpmkey  bad pixel mask header keyword [BPM]
        bpmnew  new bad pixel mask for shifted image [none]

        clobber clobber output files [yes]
        verbose print messages about actions [yes]
    """

    # Input checking
    if not os.path.exists(input) and not os.path.exists(input+".fits"):
        print "Couldn't open input image %s" % input
        return
    elif os.path.exists(input+".fits"):
        input+=".fits"
    if input==output:
        print "Can't shift to same filename, sorry"
        return
    check_exist(output,'w',clobber)

    # Open the input image as pyfits object
    fimg=pyfits.open(input)
    hdr=fimg[0].header
    D=fimg[0].data
    (iny,inx)=D.shape

    # The shift
    dx,dy=int(shift[0]),int(shift[1])
    szx,szy=inx-abs(dx),iny-abs(dy)

    # Coordinate ranges for zero/positive shifts
    x0,y0=0,0
    newx0,newy0=dx,dy

    # Change for negative shifts
    if dx<0:
        x0=-dx
        newx0=0
    if dy<0:
        y0=-dy
        newy0=0

    # Make the output data as a numarray
    DX=0*D+border

    # Insert data from the input image
    DX[newy0:newy0+szy,newx0:newx0+szx]=D[y0:y0+szy,x0:x0+szx]

    fimg[0].data=DX
    hdr.update('SHIFTED','Shifted from %s' % input)

    # Correct the WCS (CRPIX) settings
    if hdr.has_key('CRPIX1'):
        crpix1=hdr['CRPIX1']
        hdr.update('CRPIX1',crpix1+dx)
    if hdr.has_key('CRPIX1'):
        crpix2=hdr['CRPIX2']
        hdr.update('CRPIX2',crpix2+dy)

    # Reset the BPM keyword if requested
    if len(bpmkey)>0:
        if len(bpmnew)>0:
            hdr.update(bpmkey,bpmnew)
        else:
            del hdr['BPM']

    # Set or adjust the SKYSEC keyword, as appropriate
    if len(skysec)>0:
        # Default sky region = Full original image (IRAF style)
        skyreg="[%d:%d,%d:%d]" % (newx0+1,newx0+szx,newy0+1,newy0+szy)
        # Check if there was a sky region defined already
        if check_head(input,skysec):
            skyval=get_head(input,skysec)
            resky=re.search("\[(\d+):(\d+),(\d+):(\d+)\]",skyval)
            if resky:
                oldx0,oldx1,oldy0,oldy1=int(resky.group(1)), \
                                        int(resky.group(2)), \
                                        int(resky.group(3)), \
                                        int(resky.group(4))
                skyreg="[%d:%d,%d:%d]" % \
                        (oldx0+dx,oldx1+dx,
                         oldy0+dy,oldy1+dy)
        hdr.update(skysec,skyreg)

    # Write/close the output image
    fimg.writeto(output)

##############################

def expand_image(input,output,size=None,scale=2,llxy=None,
                 border=0,nodata=0,bpmkey="BPM",bpmnew="",
                 skysec="SKYSEC",clobber=globclob,verbose=globver):

    """ takes input image and places it in the appropriate portion
        of the (larger) output image.  WCS is preserved.  

        input   name of input image
        output  name of output image
        size    size of output image [None]
        scale   factor by which to expand if size==None [2]
        llxy    location of lower-left corner in expanded image
                (centered in output if not set)
        border  value to set in border regions [0]
        nodata  make a blank output image? [False]
        bpmkey  bad pixel mask header keyword [BPM]
        bpmnew  new bad pixel mask for expanded image [none]
        detkey  detector section header keyword [DETSEC]

        clobber clobber output files [yes]
        verbose print messages about actions [yes]
    """

    # Input checking
    if not os.path.exists(input) and not os.path.exists(input+".fits"):
        print "Couldn't open input image %s" % input
        return
    elif os.path.exists(input+".fits"):
        input+=".fits"
    if input==output:
        print "Can't expand_image to same filename, sorry"
        return
    check_exist(output,"w",clobber)

    # Open the input image as pyfits object
    fimg=pyfits.open(input)
    hdr=fimg[0].header
    D=fimg[0].data
    (iny,inx)=D.shape

    # Determine the size of the output image
    outx,outy=None,None
    if type(size) is ListType and len(size)>0:
        if len(size)==1:
            outx,outy=int(size[0]),int(size[0])
        elif len(size)>1:
            outx,outy=int(size[0]),int(size[1])
        if outx<inx or outy<iny:
            print "Can't make a smaller image than I start with..."
            return
    elif type(size) != NoneType:
        try:
            isize=int(size)
            outx,outy=isize,isize
        except:
            print "Failed to parse input size as integer"

    # Use the scale keyword
    if type(outx) is NoneType:
        scale=int(scale)
        if (scale<1):
            print "Bad scaling in expand_image, assuming scale==2"
            scale=2
        [outx,outy]=[int(scale*inx),int(scale*iny)]

    # Determine the location of the data in the output image
    if type(llxy) is not ListType or len(llxy)!=2:
        # Default:  Center old image in expanded frame
        xstart=int(outx/2-inx/2)
        ystart=int(outy/2-iny/2)
    else:
        # Follow directions if legal
        [xstart,ystart]=[int(llxy[0]),int(llxy[1])]
    # Note Python-style indices
    xend=xstart+inx
    yend=ystart+iny

    # Talk about what we're going to do
    if verbose:
        print "Image %s is %d x %d" % (input,inx,iny)
        print "Output image %s will be %d x %d" % (output,outx,outy)
        print "With data beginning at IRAF index (%d,%d)" % (xstart+1,ystart+1)

    # Make the output data as a numarray
    DX=0*resize(D,(outy,outx))+border

    # Standard settings
    x0,y0=0,0
    szx,szy=inx,iny
    newx0,newy0=xstart,ystart

    # Check against expanded image boundaries
    if newx0<0:
        newx0=0
        x0=-xstart
        szx-=x0
    elif newx0+szx>outx:
        szx=outx-newx0
    if newy0<0:
        newy0=0
        y0=-ystart
        szy-=y0
    elif newy0+szy>outy:
        szy=outy-newy0

    # Insert data from the input image
    if nodata:
        pass
    else:
        DX[newy0:newy0+szy,newx0:newx0+szx]=D[y0:y0+szy,x0:x0+szx]

    fimg[0].data=DX
    hdr.update('EXPANDED','Expanded from %s' % input)

    # Correct the WCS (CRPIX) settings
    if hdr.has_key('CRPIX1'):
        crpix1=hdr['CRPIX1']
        hdr.update('CRPIX1',crpix1+xstart)
    if hdr.has_key('CRPIX1'):
        crpix2=hdr['CRPIX2']
        hdr.update('CRPIX2',crpix2+ystart)

    # Reset the BPM keyword if requested
    if len(bpmkey)>0:
        if len(bpmnew)>0:
            hdr.update(bpmkey,bpmnew)
        else:
            del hdr['BPM']

    # Set the SKYSEC keyword if appropriate
    if len(skysec)>0:
        # Default sky region = Full original image (IRAF style)
        skyreg="[%d:%d,%d:%d]" % (newx0+1,newx0+szx,newy0+1,newy0+szy)
        # Check if there was a sky region defined already
        if check_head(input,skysec):
            skyval=get_head(input,skysec)
            resky=re.search("\[(\d+):(\d+),(\d+):(\d+)\]",skyval)
            if resky:
                oldx0,oldx1,oldy0,oldy1=int(resky.group(1)), \
                                        int(resky.group(2)), \
                                        int(resky.group(3)), \
                                        int(resky.group(4))
                newskyx0,newskyx1=oldx0+xstart,oldx1+xstart
                newskyy0,newskyy1=oldy0+ystart,oldy1+ystart
                # Clip to new image
                if newskyx0<0:     newskyx0=0
                if newskyx1>outx:  newskyx1=outx
                if newskyy0<0:     newskyy0=0
                if newskyy1>outy:  newskyy1=outy
                skyreg="[%d:%d,%d:%d]" % \
                        (newskyx0,newskyx1,
                         newskyy0,newskyy1)
        hdr.update(skysec,skyreg)

    # Write/close the output image
    fimg.writeto(output)

##############################

def mkbpmfits(bpmpl,delete=0,clobber=globclob,verbose=globver):

    bpmfits=bpmpl.replace('.pl','.fits')
    check_exist(bpmfits,'w',clobber)
    iraf.imcopy(bpmpl,bpmfits,verbose=verbose)
    if delete:
        os.remove(bpmpl)

    return bpmfits
    
##############################

def mkbpmpl(bpmfits,delete=0,clobber=globclob,verbose=globver):

    bpmpl=bpmfits.replace('.fits','.pl')
    check_exist(bpmpl,'w',clobber)
    iraf.imcopy(bpmfits,bpmpl,verbose=verbose)
    if delete:
        os.remove(bpmfits)

    return bpmpl
    
##############################

def blur_bpm(input,output,xy=[[1,1],[0,1],[1,0]],
             clobber=globclob,verbose=globver):

    """ takes input BPM (fits form) and blurs it by taking of max of
        it and its shifted version(s)

        input   name of input BPM image
        output  name of output BPM image
        xy      list-of-lists giving coordinates for blurring
                default:  [[1,1],[0,1],[1,0]]

        clobber clobber output files [yes]
        verbose print messages about actions [yes]
    """

    # Input checking
    if not os.path.exists(input) and not os.path.exists(input+".fits"):
        print "Couldn't open input image %s" % input
        return
    elif os.path.exists(input+".fits"):
        input+=".fits"
    check_exist(output,"w",clobber)

    # Open the input image as pyfits object
    fimg=pyfits.open(input)
    hdr=fimg[0].header
    D=fimg[0].data
    (iny,inx)=D.shape

    # Determine the full size of the shifted BPM's
    xshift=[]
    yshift=[]
    if len(xy)==2 and len(xy[0])==1:
        # User provided a single (x,y) shift
        xshift.append(int(round(xy[0])))
        yshift.append(int(round(xy[1])))
    else:
        # User provided a list of (x,y) shifts
        for shift in xy:
            xshift.append(int(round(shift[0])))
            yshift.append(int(round(shift[1])))

    # Make the output data as a numarray
    DX=copy.deepcopy(D)

    # Insert data from the input image
    for i in xrange(len(xshift)):
        dx,dy = xshift[i],yshift[i]
        x0,y0 = max(0,dx),max(0,dy)
        x1,y1 = min(inx,inx+dx),min(iny,iny+dy)
        x2,y2 = max(0,-dx),max(0,-dy)
        x3,y3 = min(inx,inx-dx),min(iny,iny-dy)
        
        DX[y0:y1,x0:x1] = \
             choose(greater(DX[y0:y1,x0:x1],D[y2:y3,x2:x3]),
                    (D[y2:y3,x2:x3],DX[y0:y1,x0:x1]))
                    
    fimg[0].data=DX
    hdr.update('BLURRED','Blurred BPM from %s' % input)

    # Write/close the output image
    fimg.writeto(output)

##############################

def apply_bpm(target,bpm,value,
              clobber=globclob,verbose=globver):

    """ simple BPM masking:  for all pixels where BPM is non-zero,
        change the value of the pixel in the TARGET image to VALUE. 

        target   name of target image
        bpm      name of BPM image (FITS format): goodval=0 type
        value    substitude value masked pixels -- either floating
                 point or image filename 
 
        clobber clobber output files [yes]
        verbose print messages about actions [yes]
    """

    # Input checking
    if not os.path.exists(target) and not os.path.exists(target+".fits"):
        print "Couldn't open target image %s" % target
        return
    elif os.path.exists(target+".fits"):
        target+=".fits"

    if not os.path.exists(bpm) and not os.path.exists(bpm+".fits"):
        print "Couldn't open BPM image %s" % bpm
        return
    elif os.path.exists(bpm+".fits"):
        bpm+=".fits"

    useimage=0
    if type(value) is StringType:
        refimg=value
        if not os.path.exists(refimg) and not os.path.exists(refimg+".fits"):
            print "Couldn't open reference image %s" % refimg
            return
        elif os.path.exists(refimg+".fits"):
            refimg+=".fits"
        useimage=1
    else:
        try:
            value=float(value)
        except:
            print "Couldn't parse value as floating point"
            return

    # Open the input image as pyfits object
    fimg=pyfits.open(target,"update")
    hdr=fimg[0].header
    D=fimg[0].data
    (iny,inx)=D.shape

    # Open the reference image as pyfits object
    if useimage:
        fref=pyfits.open(refimg)
        F=fref[0].data
        (refy,refx)=F.shape

    # Open the BPM image as pyfits object
    fbpm=pyfits.open(bpm)
    B=fbpm[0].data
    (bpmy,bpmx)=B.shape

    # Sanity checks
    if inx!=bpmx or iny!=bpmy:
        print "Image and BPM are not of same size, quitting"
        return

    if useimage:
        if inx!=refx or iny!=refy:
            print "Image and reference are not of same size, quitting"
            return

    # Apply BPM
    if useimage:
        # Replace bad pixels with reference image pixel values
        fimg[0].data=choose(greater(B,0),(D,F))
    else:
        # Replace bad pixels with specified floating point value
        fimg[0].data=choose(greater(B,0),(D,value))

    # Close the images
    fimg.flush()
    fimg.close()
    fbpm.close()

    if useimage:
        fref.close()

##############################

def wcsinterp(badlist,goodlist,rapoint='RA',decpoint='DEC',
              ratrue='CRVAL1',dectrue='CRVAL2',
              avgkeys=['CD1_1','CD1_2','CD2_1','CD2_2'],
              cpkeys=['CTYPE1','CTYPE2','CRPIX1','CRPIX2']):

    """ Use offsets for images with good WCS fits to refine the WCS in
        images without good WCS fits (astrometry from raw pointing only) """

    # Images with uncorrected WCS to fix
    if type(badlist) is ListType:
        badfiles=badlist
    elif type(badlist) is StringType:
        badfiles=iraffiles(badlist)
    else:
        print "Please pass a string or list of bad files"
        return

    # Images with good WCS to use as reference
    if type(goodlist) is ListType:
        goodfiles=goodlist
    elif type(goodlist) is StringType:
        goodfiles=iraffiles(goodlist)
    else:
        print "Please pass a string or list of good files"
        return

    # Key values we will copy to the new file
    good1=goodfiles[0]
    cpvals=get_head(good1,cpkeys)

    # Dictionary for key values we will average for the new file
    avgvals={}
    for akey in avgkeys:
        avgvals[akey]=[]

    # Loop over good images, collecting offsets & keyword values...
    dras=[]
    ddcs=[]
    for image in goodfiles:
        imgvals=get_head(image,avgkeys)
        for i in xrange(len(avgkeys)):
            avgvals[avgkeys[i]].append(float(imgvals[i]))
        # Coordinate offset from header position
        coovals=get_head(image,[rapoint,decpoint,ratrue,dectrue])
        coopoint=astrocoords(coovals[0],coovals[1])
        cootrue=astrocoords(coovals[2],coovals[3])
        [ra1,dc1]=coopoint.deg()
        [ra2,dc2]=cootrue.deg()
        [dra,ddc]=radecdiff(ra2,dc2,ra1,dc1,arcsec=1)
        dras.append(dra)
        ddcs.append(ddc)

    # Calculate the average values
    ameans=[]
    for i in xrange(len(avgkeys)):
        akey=avgkeys[i]
        ameans.append(mean(avgvals[akey]))

    # Calculate the average coordinate offset
    raoff=mean(dras)
    dcoff=mean(ddcs)

    # Loop over bad images, correcting them as we go...
    for badimage in badfiles:

        for i in xrange(len(cpkeys)):
            update_head(badimage,cpkeys[i],cpvals[i])

        # Calculate the average values & update the bad image
        for i in xrange(len(avgkeys)):
            update_head(badimage,avgkeys[i],ameans[i])

        # Correct the pointing coordinates
        ptvals=get_head(badimage,[rapoint,decpoint])
        coopt=astrocoords(ptvals[0],ptvals[1])
        coopt.shift([raoff,dcoff],arcsec=1)
        [radeg,dcdeg]=coopt.deg()
        update_head(badimage,ratrue,radeg)
        update_head(badimage,dectrue,dcdeg)

        # The end
        update_head(badimage,'WCSINTRP',1,
                    'WCS has been interpolated from other images')

##############################

def wcsaverage(inlist,modlist,crpix=['CRPIX1','CRPIX2'],
               ratrue='CRVAL1',dectrue='CRVAL2',
               avgkeys=['CD1_1','CD1_2','CD2_1','CD2_2'],
               cpkeys=['CTYPE1','CTYPE2']):

    """ Average the WCS on images in inlist and apply the average WCS
        to the images in modlist """

    # Images with WCS to average
    if type(inlist) is ListType:
        infiles=inlist
    elif type(inlist) is StringType:
        infiles=iraffiles(inlist)
    else:
        print "Please pass a string or list of input files"
        return

    # Images to update
    if type(modlist) is ListType:
        modfiles=modlist
    elif type(modlist) is StringType:
        modfiles=iraffiles(modlist)
    else:
        print "Please pass a string or list of files to modify"
        return

    # Key values we will copy to the new file
    in1=infiles[0]
    cpvals=get_head(in1,cpkeys)
    crpixv=get_head(in1,crpix)

    # Dictionary for key values we will average for the new file
    avgvals={}
    for akey in avgkeys:
        avgvals[akey]=[]

    # Loop over input images, collecting RAs, Decs & keyword values...
    ravs=[]
    dcvs=[]
    for image in infiles:
        # RA and Dec at the CRPIX
        [[crval1,crval2]]=impix2wcs(image,crpixv[0],crpixv[1])
        coo=astrocoords(crval1,crval2)
        ravs.append(coo.radeg())
        dcvs.append(coo.dcdeg())
        # CD Matrix values to average
        gdakeys=[]
        imgvals=get_head(image,avgkeys)
        for i in xrange(len(avgkeys)):
            try:
                avgvals[avgkeys[i]].append(float(imgvals[i]))
                gdakeys.append(avgkeys[i])
            except:
                # ignore this keyword for this image
                pass

    # Calculate the average values
    ameans=[]
    for i in xrange(len(gdakeys)):
        akey=gdakeys[i]
        ameans.append(mean(avgvals[akey]))

    # Calculate the average coordinate at center
    raavg=mean(ravs)
    dcavg=mean(dcvs)

    # Loop over new images & correct
    for image in modfiles:

        # Update copy keys (transform type)
        for i in xrange(len(cpkeys)):
            update_head(image,cpkeys[i],cpvals[i])

        # Update average keys (CD matrix)
        for i in xrange(len(gdakeys)):
            update_head(image,gdakeys[i],ameans[i])

        # Update pointing (CRPIX and RA/Dec)
        update_head(image,crpix[0],crpixv[0])
        update_head(image,crpix[1],crpixv[1])
        update_head(image,ratrue,raavg)
        update_head(image,dectrue,dcavg)

##############################

def wcsbox(inlist,rain,dcin,width,llcoo=0,
           arcsec=0,arcmin=1,degree=0):

    """ Determine pixel coordinates of a region to cut out of one or
        more images by using their WCS coords """

    # Images to cut out
    if type(inlist) is ListType:
        infiles=inlist
    elif type(inlist) is StringType:
        infiles=iraffiles(inlist)
    else:
        print "Please pass a string or list of input files"
        return

    # Width/height to cut out
    if type(width) is ListType:
        if len(width)==1:
            rawide,dchigh=float(width[0]),float(width[0])
        elif len(width)>1:
            rawide,dchigh=float(width[0]),float(width[1])
        else:
            print "Missing width/height argument, please supply"
            return
    else:
        try:
            rawide,dchigh=float(width),float(width)
        except:
            print "Failed to parse width/height as float"
            return

    # Convert width/height into arcmin
    if arcmin and not arcsec and not degree:
        # stick with arcmin
        pass
    elif arcsec:
        rawide,dchigh=rawide/60.0,dchigh/60.0
    elif degree:
        rawide,dchigh=rawide*60.0,dchigh*60.0
    else:
        # stick with arcmin
        pass
    
    # Astro coordinates of southeast corner
    coose=astrocoords(rain,dcin)
    if not llcoo:
        coose.shift([0.5*rawide,-0.5*dchigh],arcmin=1)

    # Loop over images
    boxes=[]
    for image in infiles:
    
        # Pixel coordinates of southeast corner
        pixse=imwcs2pix(image,coose.radeg(),coose.dcdeg())
        rndsex,rndsey=round(pixse[0]),round(pixse[1])
        coose2=impix2wcs(image,rndsex,rndsey)
        
        # Coordinates of the northwest corner
        coonw=astrocoords(coose2[0],coose2[1])
        coonw.shift([-1.0*rawide,1.0*dchigh],arcmin=1)
        pixnw=imwcs2pix(image,coonw.radeg(),coonw.dcdeg())
        rndnwx,rndnwy=round(pixnw[0]),round(pixnw[1])

        # Package these nicely
        box=[[int(min(rndsex,rndnwx)),int(max(rndsex,rndnwx))],
             [int(min(rndsey,rndnwy)),int(max(rndsey,rndnwy))]]
        boxes.append(box)

    # The End
    return boxes

##############################

def imextent(inlist,arcsec=0,arcmin=1,degree=0):

    """ Determine the extent of an image in RA & Dec by taking
        the Min/Max of each value from the corners of the image """

    # Images to work with
    if type(inlist) is ListType:
        infiles=inlist
    elif type(inlist) is StringType:
        infiles=iraffiles(inlist)
    else:
        print "Please pass a string or list of input files"
        return

    # Loop over images
    widehigh=[]
    for image in infiles:

        # Image size
        [n1,n2]=get_head(image,['NAXIS1','NAXIS2'])

        # Construct list of pixel coordinates
        xvals=[0.5,0.5,n1+0.5,n1+0.5]
        yvals=[0.5,n2+0.5,0.5,n2+0.5]
    
        # Transform to WCS coordinates
        wcsvals=impix2wcs(image,xvals,yvals,physical=0)

        # Find RA/Dec differences between each pair of points
        # Save the largest in each dimension
        maxrad=0.0
        maxdcd=0.0

        for i in xrange(3):
            for j in xrange(i+1,4):
                # convert to decimal degrees
                coo1=astrocoords(wcsvals[i][0],wcsvals[i][1])
                coo2=astrocoords(wcsvals[j][0],wcsvals[j][1])
                radd=radecdiff(coo1.radeg(),coo1.dcdeg(),
                               coo2.radeg(),coo2.dcdeg(),
                               arcsec=arcsec,arcmin=arcmin,
                               degree=degree)
                if abs(radd[0])>maxrad:
                    maxrad=abs(radd[0])
                if abs(radd[1])>maxdcd:
                    maxdcd=abs(radd[1])

        # Package the result
        widehigh.append([maxrad,maxdcd])

    # Unpack by one layer if we only did one image
    if len(widehigh)==1:
        widehigh=widehigh[0]

    # The End
    return widehigh

##############################

def imwcs2pix(image,rain,dcin,physical=0):

    """ Use image WCS to transform coords to pixel values """

    check_exist(image,"r")

    iraf.imcoords()

    if type(rain)==ListType:
        npts=len(rain)
        rax=rain
        dcx=dcin
    elif type(rain)==StringType or type(rain)==FloatType:
        npts=1
        rax=[rain]
        dcx=[dcin]
    else:
        print "Unrecognized type for RA, Dec"
        return None

    # Temporary file names
    finnm=iraf.mktemp("iqut")+".wcs"
    foutnm=iraf.mktemp("iqut")+".pix"

    # Write the list of coordinates
    fin=open(finnm,mode='w')

    for i in range(0,npts):
        [ra,dec]=[rax[i],dcx[i]]
        coo=astrocoords(ra,dec)
        fin.write("%14.10f %14.9f\n" % (coo.radeg()/15,coo.dcdeg()))

    fin.close()

    # Determine if we are using image or physical pixels
    if physical:
        outwcs="physical"
    else:
        outwcs="logical"

    # Run wcsctran
    #***Need to reset task parameters b/c of IRAF wcsctran bug
    # http://iraf.net/article.php?story=7306&query=wcsctran***
    iraf.flpr()
    s=iraf.wcsctran(finnm,foutnm,image,inwcs="world",outwcs=outwcs,
                    columns="1 2",units="hours degrees",
                    formats="%10.5f %10.5f",min_sigdigits=7,
                    verbose=no,Stdout=1)
        
    # Read in the answers
    pixout=getlines(foutnm)
    ans=[]
    for line in pixout:
        if re.search('^#',line):
            continue
        els=line.split()
        try:
            ans.append([float(els[0]),float(els[1])])
        except:
            ans.append([-1.0,-1.0])

    # Reduce dimensionality if npts=1
    #if npts==1:
        #ans=ans[0]

    # Clean up
    os.remove(finnm)
    os.remove(foutnm)

    # The End
    return ans

##############################

def impix2wcs(image,xin,yin,physical=0,deg=0):

    """ Use image WCS to transform pixel values to WCS coords """

    check_exist(image,"r")

    iraf.imcoords()

    if type(xin)==ListType:
        npts=len(xin)
        xs=xin
        ys=yin
    elif type(xin)==StringType:
        npts=1
        xs=[float(xin)]
        ys=[float(yin)]
    elif type(xin)==FloatType or type(xin)==IntType:
        npts=1
        xs=[xin]
        ys=[yin]
    else:
        print "Unrecognized type for pixel coords x, y"
        return None

    # Temporary file names
    finnm=iraf.mktemp("iqut")+".pix"
    foutnm=iraf.mktemp("iqut")+".wcs"

    # Write the list of coordinates
    fin=open(finnm,mode='w')

    for i in range(0,npts):
        [xpt,ypt]=[xs[i],ys[i]]
        fin.write("%.3f %.3f\n" % (xpt,ypt))

    fin.close()

    # Determine if we are using image or physical pixels
    if physical:
        inwcs="physical"
    else:
        inwcs="logical"

    # Run wcsctran
    #***Need to reset task parameters b/c of IRAF wcsctran bug
    # http://iraf.net/article.php?story=7306&query=wcsctran***
    iraf.flpr()
    if deg:
        s=iraf.wcsctran(finnm,foutnm,image,inwcs=inwcs,outwcs="world",
                        columns="1 2",units="",formats="%11.5f %10.5f",
                        min_sigdigits=7,verbose=no,Stdout=1)
    else:
        s=iraf.wcsctran(finnm,foutnm,image,inwcs=inwcs,outwcs="world",
                        columns="1 2",units="",formats="%16.6H %16.5h",
                        min_sigdigits=7,verbose=no,Stdout=1)
        
    # Read in the answers
    wcsout=getlines(foutnm)
    ans=[]
    for line in wcsout:
        if re.search('^#',line):
            continue
        els=line.split()
        ans.append([els[0],els[1]])

    # Reduce dimensionality if npts=1
    #if npts==1:
        #ans=ans[0]

    # Clean up
    os.remove(finnm)
    os.remove(foutnm)

    # The End
    return ans

##############################

def impix(image):

    """ Get the pixel scale of an image (arcsec/pixel) """

    check_exist(image,"r")

    [cd1_1,cd1_2,cd2_1,cd2_2] = \
           get_head(image,['CD1_1','CD1_2','CD2_1','CD2_2'])

    # Try to use the official CD matrix
    scale=None
    if type(cd1_2)==FloatType and type(cd2_1)==FloatType and \
           type(cd1_1)==FloatType and type(cd2_2)==FloatType:
        scale=3600*math.sqrt(abs(cd1_2*cd2_1 - cd1_1*cd2_2))
    elif type(cd1_2)==FloatType and type(cd2_1)==FloatType:
        scale=3600*math.sqrt(abs(cd1_2*cd2_1))
    elif type(cd1_1)==FloatType and type(cd2_2)==FloatType:
        scale=3600*math.sqrt(abs(cd1_1*cd2_2))
    else:
        # Try to get the pixel scale directly
        [pix0,pix1,pix2] = \
              get_head(image,['PIXSCALE','PIXSCAL1','PIXSCAL2'])

        if type(pix0)==FloatType:
            scale=pix0
        elif type(pix1)==FloatType:
            scale=pix1
        elif type(pix2)==FloatType:
            scale=pix2
        
    # Make sure we have an answer
    if type(scale)==NoneType:
        print "Failed to get pixel scale from '%s' header" % image

    return scale

##############################

def imccmap(image,refra,refdec,fit='fit',degree=2,ccmap=None,
            verbose=no,clobber=yes):

    """ Perform a standard interactive run of IRAF ccmap """

    # Images to process
    if type(image) is ListType:
        files=image
    elif type(image) is StringType:
        files=iraffiles(image)
    else:
        print "Please pass a string or list of input files"
        return

    # Reference coordinates
    if type(refra) is not StringType or type(refdec) is not StringType:
        print "Please provide reference RA, Dec as sexagesimal hours, degrees"
        return
    
    # Degree of fit
    if degree!=2 and degree!=4:
        if degree>2:
            degree=4
        else:
            degree=2
        print "Only degree=2 or degree=4 fits are allowed."
        print "Setting degree=%d" % degree
    
    # Loop over files
    for file in files:

        if not check_exist(file,"r"):
            continue

        if not ccmap:
            froot,fext=os.path.splitext(file)
            ccmap=froot+'.ccmap'
        else:
            froot,fext=os.path.splitext(ccmap)
            
        if not check_exist(ccmap,"r"):
            print "File %s not found." % ccmap
            continue

        # Output files:
        #   .ccdb (appended to)
        #   .ccout (clobbered if necessary & appropriate)

        ccdb=froot+'.ccdb'
        ccout=froot+'-'+fit+'.ccout'
        if os.path.exists(ccout):
            if clobber:
                os.remove(ccout)
            else:
                print "CCOUT file %s exists and clobber=no, skipping" % ccout
                continue

        # Run IRAF ccmap
        iraf.ccmap(ccmap,ccdb,solutions=fit,images=file,results=ccout,
                   xcolumn=1,ycolumn=2,lngcolumn=3,latcolumn=4,
                   lngunits='hours',latunits='degrees',insystem='j2000',
                   lngref=refra,latref=refdec,refsystem='j2000',
                   lngrefunits='hours',latrefunits='degrees',
                   projection='tan',fitgeometry='general',
                   function='polynomial',xxorder=degree,xyorder=degree,
                   xxterms='half',yxorder=degree,yyorder=degree,
                   yxterms='half',maxiter=0,reject=3.0,
                   update=yes,pixsystem='logical',verbose=verbose,
                   interactive=yes)

##############################

def imgeotran(image,output,refimg,name='geomap',
              fitgeometry='rscale',degree=2,xterms='half',
              geomap=None,imatch=None,rmatch=None,
              istars='!STARLIST',rstars='!STARLIST',tol=3.0,useflags=yes,
              maxnum=50,maxskip=None,verbose=no,clobber=yes):

    """ Perform interactive run of IRAF geomap+geotran """

    # Degree of fit
    if degree>2:
        if fitgeometry != 'general':
            print "Fit order>2 requires general fit geometry... resetting"
            fitgeometry='general'

    # Images to process
    if type(image) is not StringType or type(refimg) is not StringType:
        print "Please provide image and refimg as strings"
        return

    if not check_exist(image,"r"):
        print "Failed to find target image '%s'" % image
        return

    if not check_exist(refimg,"r"):
        print "Failed to find reference image '%s'" % refimg
        return

    # Size of reference image
    [refnx,refny]=get_head(refimg,['NAXIS1','NAXIS2'])

    # Geomap database name
    if not geomap:
        froot,fext=os.path.splitext(image)
        geomap=froot+'.geomap'
    else:
        froot,fext=os.path.splitext(geomap)

    # Check if we need to run geomap
    if not check_exist(geomap,"r"):

        # Check if we need to match starlists
        if not imatch or not rmatch:

            if istars[0]=='!':
                istrf=get_head(image,istars[1:])
            else:
                istrf=istars

            if not check_exist(istrf,"r"):
                print "Failed to find target starlist '%s'" % istrf
                return

            if rstars[0]=='!':
                rstrf=get_head(image,rstars[1:])
            else:
                rstrf=rstars

            if not check_exist(rstrf,"r"):
                print "Failed to find reference starlist '%s'" % rstrf
                return

            # Read in the starlists
            istr=Starlist(istrf)
            rstr=Starlist(rstrf)

            # Remember reference image pixel coordinates
            for rs in rstr:
                rs.xorg=rs.xval
                rs.yorg=rs.yval

            # Transform to target image pixel coordinates
            rstr.pix2wcs(refimg)
            rstr.wcs2pix(image)

            # Execute the match
            imatch,rmatch=istr.match(rstr,tol=tol,useflags=useflags,
                                     region=region,maxnum=maxnum,maxskip=maxskip)

            # Restore reference image pixel coordinates
            for rs in rmatch:
                rs.xval=rs.xorg
                rs.yval=rs.yorg

        # Write matched starlists to geomap input file
        gmpinp=froot+'.gmpinp'
        check_exist(gmpinp,'w',clobber=clobber)

        fgin=open(gmpinp,mode='w')
        for i in xrange(len(imatch)):
            fgin.write("%.2f %.2f %.2f %.2f\n" % \
                       (rmatch[i].xval,rmatch[i].yval,
                        imatch[i].xval,imatch[i].yval))
        fgin.close()
                       
        # Run geomap
        gmpout=froot+'.gmpout'
        check_exist(gmpout,'w',clobber=clobber)

        iraf.geomap(gmpinp,geomap,xmin=1,xmax=refnx,ymin=1,ymax=refny,
                    transforms=name,results=gmpout,
                    fitgeometry=fitgeometry,function='polynomial',
                    xxorder=degree,xyorder=degree,yxorder=degree,yyorder=degree,
                    xxterms=xterms,yxterms=xterms,
                    maxiter=0,reject=3.0,calctype='real',verbose=verbose,
                    interactive=yes)

    # Run geotran
    iraf.geotran(image,output,geomap,name,geometry='geometric',
                 xmin=1,xmax=refnx,ymin=1,ymax=refny,
                 xscale=INDEF,yscale=INDEF,ncols=INDEF,nlines=INDEF,
                 xsample=1.0,ysample=1.0,interpolant='linear',
                 boundary='nearest',constant=0.0,fluxconserve=yes,
                 xin=INDEF,yin=INDEF,xout=INDEF,yout=INDEF,
                 xshift=INDEF,yshift=INDEF,xmag=INDEF,ymag=INDEF,
                 xrotation=INDEF,yrotation=INDEF,nxblock=1024,nyblock=1024,
                 verbose=verbose)

##############################

def sphdist(ra1,dc1,ra2,dc2,arcsec=1,arcmin=0,degree=0):

    """ Use Haversine's formula to calculate distance on a sphere """

    ddtor=0.0174532925199433

    ra1r=ra1*ddtor
    dc1r=dc1*ddtor
    ra2r=ra2*ddtor
    dc2r=dc2*ddtor
    
    dra=ra2r-ra1r
    ddc=dc2r-dc1r
    a=math.sin(ddc/2)**2 + math.cos(dc1r)*math.cos(dc2r)*math.sin(dra/2)**2
    cdeg=2*math.asin(min(1,math.sqrt(a)))/ddtor

    if arcsec>0 and arcmin==0 and degree==0:
        cans=3600*cdeg
    elif arcmin>0:
        cans=60*cdeg
    else:
        cans=cdeg

    return cans

##############################

def sphdist2(ra1,dc1,ra2,dc2,arcsec=1,arcmin=0,degree=0):

    """ Use Law of Cosines to calculate distance on a sphere """

    ddtor=0.0174532925199433
    
    dra=(ra2-ra1)*ddtor

    cd = math.sin(dc1*ddtor)*math.sin(dc2*ddtor) + \
         math.cos(dc1*ddtor)*math.cos(dc2*ddtor)*math.cos(dra)
    c=math.acos(cd)/ddtor

    if arcsec>0 and arcmin==0 and degree==0:
        cans=3600*c
    elif arcmin>0:
        cans=60*c
    else:
        cans=c

    return cans

##############################

def radecdiff(ra1,dc1,ra2,dc2,arcsec=1,arcmin=0,degree=0):

    # RA Distance
    radiff=sphdist(ra1,dc1,ra2,dc1,
                   arcsec=arcsec,arcmin=arcmin,degree=degree)
    # RA Sign
    rasign=+1
    if ra1<ra2:
        rasign=-1
    # This will be wrong if the two RA's straddle 24h...
    if abs(ra1-ra2)>180:
	# ...so correct it
	rasign=-rasign
    # RA Answer
    radiff *= rasign

    # Dec Distance
    dcdiff=sphdist(ra2,dc1,ra2,dc2,
                   arcsec=arcsec,arcmin=arcmin,degree=degree)
    # Dec Sign
    dcsign=+1
    if dc1<dc2:
        dcsign=-1
    # Dec Answer
    dcdiff *= dcsign

    # The End
    return [radiff,dcdiff]

##############################

def list2str(list,format):

    """ Convert list of Python objects to a string by individual
        formatting of each list element according to the specified
        format """

    out='['
    for i in range(0,len(list)-1):
        out += (format % list[i]) + ', '
    out += (format % list[-1]) + ']'
    
    return out

##############################

def list2int(list):

    """ Convert list of Python objects to list of ints by individual
        conversion of each list element """

    out=[]
    for i in range(0,len(list)):
        try:
            out.append(int(list[i]))
        except:
            out=[]
            print "Cannot convert this list to int: "+ \
                  list.__str__()
            break

    return out

##############################

def list2float(list):

    """ Convert list of Python objects to list of floats by individual
        conversion of each list element """

    out=[]
    for i in range(0,len(list)):
        try:
            out.append(float(list[i]))
        except:
            out=[]
            print "Cannot convert this list to float: "+ \
                  list.__str__()
            break

    return out

##############################

def boxinbox(x,y,box,dims,llcoo=0,irafsty=0,
             shrinkbox=0):

    """ Check to make sure that the box of size "box", centered at
        (x,y), will fit within the larger box of size "dims".
        (x,y) specifies the lower-left corner of the box if (irafsty).
        IRAF pixel convention applies if (irafsty).  Shrink "box" to
        fit if (shrinkbox), otherwise make the box flush with the
        border but keep it the same size. """

    # Accept "box" as 1 or 2 elements
    if type(box)==ListType:
        if len(box)>=2:
            boxx,boxxy=int(box[0]),int(box[1])
        elif len(box)==1:
            boxx,boxy=int(box[0]),int(box[0])
    else:
        boxx,boxy=int(box),int(box)

    # Sanity check
    if boxx>dims[0] or boxy>dims[1]:
        print "Can't place a bigger box in a smaller one"
        return None

    # Lower-left corner of active area
    if llcoo:
        llx,lly=int(x),int(y)
    else:
        llx,lly=int(x-int(boxx/2)),int(y-int(boxy/2))

    # Sanity check #2
    if shrinkbox and llx>dims[0] or llx+boxx<0 or \
                     lly>dims[1] or lly+boxy<0:
        print "No part of box lies within DIMS"
        return [[],[]]
    
    # Note: Calculations performed using "Python style" indices

    if not shrinkbox:
        # Make box flush with borders
        if llx<0:                 xb=[0,boxx]
        elif llx+boxx>dims[0]:    xb=[dims[0]-boxx,dims[0]]
        else:                     xb=[llx,llx+boxx]

        if lly<0:                 yb=[0,boxy]
        elif lly+boxy>dims[1]:    yb=[dims[1]-boxy,dims[1]]
        else:                     yb=[lly,lly+boxy]

    else:
        # Treat out-of-bounds regions as truncated
        if llx<0:                 xb=[0,boxx+llx]
        elif llx+boxx>dims[0]:    xb=[llx,dims[0]]
        else:                     xb=[llx,llx+boxx]

        if lly<0:                 yb=[0,boxy+lly]
        elif lly+boxy>dims[1]:    yb=[lly,dims[1]]
        else:                     yb=[lly,lly+boxy]

    # Correct to IRAF-style indices, if requested
    if irafsty:
        xb[0]+=1
        yb[0]+=1

    return [xb,yb]

######################################################################
######################################################################

def azfwhm(ain,cc=[],bin=2,maxr=0,peak=0,autocorr=0):

    """ Compute the azimuthally-averaged FWHM of an array of data,
        referred to some location in the array """

    # Do autocorrelation if requested
    if autocorr:
        aa=convolve.correlate2d(ain,ain,fft=1,mode='wrap')
    else:
        aa=ain

    # Size of the input data array
    dims=aa.shape

    # User specifies nominal center with "cc" list
    if len(cc)!=2:
        if peak:
            # Find the peak of the array values
            cc=[where(aa==aa.max())[0][0],
                where(aa==aa.max())[1][0]]
        else:
            # Use the center of the array
            cc=[dims[0]/2,dims[1]/2]

    # Maximum radial extent
    if maxr==0:
        maxr=min(dims)/2

    # Construct the "distance" array
    dd=fromfunction(lambda i,j:
                    hypot(i-cc[0],j-cc[1]),dims)

    # Azimuthal radial-distance arrays
    azr=bin*arange(maxr/bin,dtype='f4')
    nbins=len(azr)
    azy=zeros(nbins,dtype='f4')

    # First bin is special case
    azy[0]=aa[where(dd<0.5*bin)].mean()

    # Equal-sized radial bins
    for i in range(1,nbins):
        #azy[i]=compress(logical_and(greater_equal(dd,azr[i]-0.5*bin),
        #                            less(dd,azr[i]+0.5*bin)),aa).mean()
        azy[i]=aa[logical_and(greater_equal(dd,azr[i]-0.5*bin),
                              less(dd,azr[i]+0.5*bin))].mean()

    # Interpolate for FWHM
    if min(azy)>0.5*azy[0]:
        print "Warning: Didn't get the FWHM in the box, setting to MAXR"
        fwhmx=maxr
    else:
        # Linear interpolation for now
        fwhmy=0.5*azy[0]
        ix1=where(greater(azy,fwhmy))[0][-1]
        if ix1>=len(azy)-1:
            ix1=len(azy)-2
        fwhmx=azr[ix1]+bin*(azy[ix1]-fwhmy)/(azy[ix1]-azy[ix1+1])

    # Correct by factor of 2 for non-autocorrelation only
    if autocorr:
        fwhm=fwhmx
    else:
        fwhm=2*fwhmx

    return fwhm

######################################################################

def median(values,minv=None,maxv=None,pick=None):

    ans=0
    xval=[]

    for val in values:
        isgood=1
        if minv:
            isgood=(isgood and val>=minv)
        if maxv:
            isgood=(isgood and val<=maxv)
        if isgood:
            xval.append(val)

    xval.sort()
    n=len(xval)

    if n<1:
        print "Need some elements for median estimate"
        return 0.0

    if n % 2 == 1:
        ans=xval[n/2]
    else:
        middle2=xval[(n/2)-1:(n/2)+1]
        if pick:
            ans=random.choice(middle2)
        else:
            try:
                ans=mean(middle2)
            except TypeError:
                ans=random.choice(middle2)

    return ans

##############################

def medunc(values,minv=None,maxv=None):

    ans=0
    xval=[]

    for val in values:
        isgood=1
        if minv:
            isgood=(isgood and val>=minv)
        if maxv:
            isgood=(isgood and val<=maxv)
        if isgood:
            xval.append(val)

    xval.sort()
    nall=len(xval)
    delta=round(math.sqrt(0.25*nall)-0.5)

    if nall<2:
        print "Need >2 elements for medunc estimate"
        return 0.0

    ix1=int(nall/2 + delta + 0.5)
    if ix1<len(xval):
        amx = xval[ix1]
    else:
        print "Too few elements for good medunc estimate"
        amx = max(xval)

    ix2=int(nall/2 - delta + 0.5)
    if ix2>=0:
        amn = xval[ix2]
    else:
        print "Too few elements for good medunc estimate"
        amn = min(xval)
        
    ans = 0.5*(amx-amn)

    return ans

##############################

def mode(values,minv=None,maxv=None,nbins=20,bin=None):

    ans=0
    xval=[]

    if minv or maxv:
        for val in values:
            isgood=1
            if minv:
                isgood=(isgood and val>=minv)
            if maxv:
                isgood=(isgood and val<=maxv)
            if isgood:
                xval.append(val)
    else:
        xval=copy.deepcopy(values)

    xval.sort()

    if nbins==0:

        ##########
        # Find the value that occurs most often (no histogram)
        ##########

        mcand=xval[0]
        run=1
        maxrun=1

        for i in range(1,len(xval)):
            if xval[i]==xval[i-1]:
                run += 1
                if run>maxrun:
                    maxrun=run
                    mcand=xval[i]
            else:
                run=1

        if maxrun==1:
            print "Warning: No element occurred with frequency>1 for mode calc"

        ans=mcand

    else:

        ##########
        # Mode calculation from a histogram of the values
        ##########

        hx=[]
        hy=[]
        
        if not minv:
            minv=xval[0]
        if not maxv:
            maxv=xval[-1]

        if nbins>0:

            # Histogram the distribution

            if bin:
                xbin=bin
                nbins=1+int((maxv-minv)/xbin)
            else:
                xbin=(1+0.01/nbins)*(maxv-minv)/(1.0*nbins)

            for i in range(0,nbins):
                hx.append(minv+(i+0.5)*xbin)
                hy.append(0)

            for i in range(0,len(xval)):
                ii=int((xval[i]-minv)/xbin)
                hy[ii] += 1

        else:

            # Adaptive binning for nbins<0
            nperbin=abs(nbins)

            i=0
            while i+nperbin<=len(xval):
                ihi=i+nperbin
                while xval[ihi-1]==xval[i]:
                    ihi += 1
                    if ihi==len(xval):
                        break
                ninbin=ihi-i
                hx.append(mean(xval[i:ihi]))
                binsz=xval[ihi-1]-xval[i]
                if binsz==0:
                    binsz=hx[-1]-hx[-2]
                hy.append((1.0*ninbin)/binsz)
                i=ihi

        # Find the peak of the histogram
        ix0=hy.index(max(hy))
        ans=hx[ix0]

    # The End
    return ans

##############################

def mean(values,minv=None,maxv=None):

    ans=0
    xval=[]

    for val in values:
        isgood=1
        if minv:
            isgood=(isgood and val>=minv)
        if maxv:
            isgood=(isgood and val<=maxv)
        if isgood:
            xval.append(val)

    if len(xval)>0:
        ans=sum(xval)/float(len(xval))

    return ans

##############################

def sign(val):
    ans = (val>0)*1 + (val<0)*(-1)
    return ans

##############################

def round(val):

    if type(val)==ListType or type(val)==ndarray:
        if type(val)==ndarray:
            ans=val.astype(Int32)
        else:
            ans=copy.deepcopy(val)
        for i in xrange(len(val)):
            ans[i]=int(math.floor(val[i]+0.5))
    else:
        ans=int(math.floor(val+0.5))
    return ans

##############################

def sdev(values,zero=None,minv=None,maxv=None):

    ans=0
    xval=[]

    if not zero:
	zval=mean(values,minv=minv,maxv=maxv)
    else:
        zval=zero

    for val in values:
        isgood=1
        if minv:
            isgood=(isgood and val>=minv)
        if maxv:
            isgood=(isgood and val<=maxv)
        if isgood:
            xval.append(val)

    if len(xval)>0:
        var=0
        for xx in xval:
	    var += (xx-zval)**2
	ans=math.sqrt(var/(len(xval)-1))

    return ans

##############################

def sum(val):
    return reduce(operator.add,val,0)	

##############################

def factorial(val):
    xout=val
    if val<2:
        xout=1
    elif val>2:
        for x in xrange(val-1,1,-1):
            xout *= x
    return xout

##############################

def poiprob(clim,lamb):
    ptot=0.0
    for i in xrange(clim+1):
        ptot += (lamb**i)*exp(-lamb)/factorial(i)
    return ptot

##############################

def reject_outliers(vals,fencelim=0.50,sigma=3,maxfrac=0.15):

    # Sort the list and get its length
    vals.sort()
    origlen=len(vals)
    if (origlen<6):
       return vals

    # Calculate "outer fence"
    lval=vals[round(fencelim/2.0*origlen)]
    uval=vals[round((1-fencelim/2.0)*origlen)]
    diff=uval-lval
    llim=lval-sigma*diff
    ulim=uval+sigma*diff
    
    # Trim outliers from list of values
    newvals=compress(greater(vals,llim),vals)
    newvals=compress(less(newvals,ulim),newvals)

    # Check to see if removed too many items
    if (len(newvals) < (1-maxfrac)*len(vals)):
       return vals
    else:
       return newvals

################################
 
def getsdss(image, outfile, cut="p60", addjkc=yes, verbose=globver):

    """Retreive SDSS Catalog for a given region of sky"""

    # Need $REDUCTION variable set for getsdss.pl
    if os.environ.has_key('REDUCTION'):
        reduction=os.environ['REDUCTION']
    else:
        print "Please define environment variable REDUCTION"
        print "This should point to the directory with the MSSSO tools."
        sys.exit(1)

    # Local variables
    tmpfile="sdss.tmp"

    # Want RA and Dec of center of array and image size
    [lenx, leny] = get_head(image, ['NAXIS1', 'NAXIS2'])
    [[ractr,dcctr]]=impix2wcs(image,0.5*(1+float(lenx)), 0.5*(1+float(leny)))
    boxsize=1.1*array(imextent(image,arcmin=1))

    # Generate command
    cmd = "$REDUCTION/getsdss.pl -r %.1f -f %s " % \
          (max(boxsize[0],boxsize[1])/2.0, outfile)
    if (cut == "p60"):
        cmd+=" -p "
    cmd+=" %s %s %s" % (ractr, dcctr, tmpfile)

    # Run command
    if verbose:
        print "Fetching SDSS Catalog with command:\n > %s" % cmd

    fcmd=os.popen(cmd)
    sdsslines=fcmd.readlines()
    fcmd.close()

    # Do we want to add JKC magnitudes as well?
    if addjkc:

        # Make sure SDSS retrieval successful
        sdssf = open(tmpfile)
        if len(sdssf.readlines()):
            sdss2jkc(outfile)

    # Remove tmpfile
    os.remove(tmpfile)

    return

################################

def jkc2sdss(infile, outfile="", outformat="reg"):

    """Augment starlist with JKC BVRI magnitudes to add SDSS griz.  Conversions
     from Jordi, Grebel, and Ammon 2006 A&A"""

    try:
        stars=Starlist(infile)
        for star in stars:
            star.mags['GPMAG']=star.mags['VMAG']+0.630*(star.mags['BMAG']- \
                               star.mags['VMAG'])-0.124
            if (star.mags['VMAG']-star.mags['RMAG']) < 0.93:
                star.mags['RPMAG']=star.mags['RMAG']+0.267*(star.mags['VMAG'] \
                                   -star.mags['RMAG'])+0.088
            else:
                star.mags['RPMAG']=star.mags['RMAG']+0.770*(star.mags['VMAG'] \
                                   -star.mags['RMAG'])-0.370
            if star.mags.has_key('IMAG'):
                star.mags['IPMAG']=star.mags['IMAG']+0.247*(star.mags['RMAG'] \
                                   -star.mags['IMAG'])+0.329
                star.mags['ZPMAG']=star.mags['RPMAG']+1.584*(star.mags['RMAG']\
                                   -star.mags['IMAG'])-0.386
            if star.mags.has_key('UMAG'):
                star.mags['UPMAG']=0.750*(star.mags['UMAG']-star.mags['BMAG'])+0.770*(star.mags['BMAG']-star.mags['VMAG'])+0.720+star.mags['GPMAG']
            else:
                star.mags['UPMAG']=0.596*(star.mags['GPMAG']-star.mags['RPMAG'])+star.mags['GPMAG']+1.110
                star.mags['UMAG']=1.333 * (star.mags['UPMAG'] - star.mags['GPMAG']) - 1.027 * (star.mags['BMAG'] - star.mags['VMAG']) - 0.960 + star.mags['BMAG']

        if not outfile:
            outfile=infile
        stars.write(outfile,format=outformat)
        return

    except:
        print "Error adding SDSS magnitudes to star list %s" % infile
        newstars=Starlist()
        newstars.write(outfile, format=outformat)
        return

################################

def sdss2jkc(infile, outfile="", outformat="reg"):

    """Augment starlist with SDSS magnitudes to add JKC BVRI.  Conversions 
       from Jordi, Grebel, and Ammon 2006 A&A"""

    try:
        stars=Starlist(infile)
        for star in stars:
            star.mags['UPMAG'] = star.mags['UMAG']
            star.magus['UPMAGU'] = star.magus['UMAGU']
            star.mags['GPMAG'] = star.mags['GMAG']
            star.magus['GPMAGU'] = star.magus['GMAGU']
            star.mags['RPMAG'] = star.mags['RMAG']
            star.magus['RPMAGU'] = star.magus['RMAGU']
            star.mags['IPMAG'] = star.mags['IMAG']
            star.magus['IPMAGU'] = star.magus['IMAGU']
            star.mags['ZPMAG'] = star.mags['ZMAG']
            star.magus['ZPMAGU'] = star.magus['ZMAGU']
            star.mags['VMAG'] = star.mags['GPMAG'] - 0.565 * \
                                (star.mags['GPMAG'] - star.mags['RPMAG'])-0.016
            star.mags['BMAG'] = star.mags['GPMAG'] + 0.313 * \
                                (star.mags['GPMAG'] - star.mags['RPMAG'])+0.219
            star.mags['RMAG'] = star.mags['RPMAG'] - 0.153 * \
                                (star.mags['RPMAG'] - star.mags['IPMAG'])-0.117
            star.mags['IMAG'] = star.mags['IPMAG'] - 0.386 * \
                                (star.mags['IPMAG'] - star.mags['ZPMAG'])-0.397
            star.mags['UMAG'] = star.mags['BMAG'] + 0.790 * \
                                (star.mags['UPMAG'] - star.mags['GPMAG'])-0.930
            #del star.mags['UMAG']
            del star.magus['UMAGU']
            del star.mags['GMAG']
            del star.magus['GMAGU']
            del star.magus['RMAGU']
            del star.magus['IMAGU']
            del star.mags['ZMAG']
            del star.magus['ZMAGU']
        
        # Write out results
        if (outfile==""):
            outfile=infile
        stars.write(outfile,format=outformat)

    except:
        print "Error reading starlist: %s" % infile
    
    return

################################

def getnomad(image, outfile, outformat="reg", addusno=yes, usnofile="",
             addsdss=yes, verbose=globver):

    """Retreive NOMAD (BVR) catalog from web.  Match objects to USNO-B catalog.
       Create new star list with BVR from NOMAD, I(2) from USNO-B, and 
       SDSS griz (conversions from Jordi, Grebel, and Ammon, A&A 2006."""

    # Need $REDUCTION variable set for getusnobn 
    if os.environ.has_key('REDUCTION'):
        reduction=os.environ['REDUCTION']
    else:
        print "Please define environment variable REDUCTION"
        print "This should point to the directory with the MSSSO tools."
        sys.exit(1)

    # Local variables
    tmpfile="nomad.tmp"

    # Want RA and Dec of center of array and image size
    [lenx, leny] = get_head(image, ['NAXIS1', 'NAXIS2'])
    [[ractr,dcctr]]=impix2wcs(image,0.5*(1+float(lenx)), 0.5*(1+float(leny)))
    boxsize=1.1*array(imextent(image,arcmin=1))

    # Generate command
    cmd = "$REDUCTION/getusnobn -c nomad -o %s -b 5.0 -f 19.0 %s %s %.2f" % \
          (tmpfile,ractr,dcctr,max(boxsize[0],boxsize[1]))

    # Run command
    if verbose:
        print "Fetching NOMAD Catalog with command:\n > %s" % cmd

    fcmd=os.popen(cmd)
    nomadlines=fcmd.readlines()
    fcmd.close()

    try:

        # Add USNO-B I-band magnitudes?
        if addusno:
            nomadpusno(tmpfile,usnofile,image,tmpfile,outformat="reg")

        # Add SDSS conversion?
        if addsdss:
            jkc2sdss(tmpfile,outfile=tmpfile,outformat="reg")

        # Write output file
        stars = Starlist(tmpfile)
        stars.write(outfile, format=outformat)
        os.remove(tmpfile)
        return

    except:
 
        print "Error creating NOMAD file: %s" % outfile
        os.remove(tmpfile)
        return

################################

def nomadpusno(nomadfile, usnofile, image, outfile, outformat="reg"):

    """Combine NOMAD and USNO-B star lists, taking BVR from NOMAD and I from
       USNO-B."""

    try:
        nomad=Starlist(nomadfile)
        usnob=Starlist(usnofile)
        nomad.wcs2pix(image)
        usnob.wcs2pix(image)
        a,b=nomad.match(usnob,maxnum=len(nomad))
        stars=[]
        
        for i in range(len(a)):
            star=a[i]
            star.mags['IMAG'] = b[i].mags['I2MAG']
            if (star.mags['BMAG'] < 25.0) and (star.mags['VMAG'] < 25.0) and \
               (star.mags['RMAG'] < 25.0) and (star.mags['IMAG'] < 25.0):
                stars.append(star)

        newstars=Starlist(stars=stars)
        newstars.write(outfile, format=outformat)
        return

    except:
        print "Error combining NOMAD (%s) and USNO-B (%s) star lists" % \
              (nomadfile, usnofile)
        newstars=Starlist()
        newstars.write(outfile, format=outformat)
        return
 
################################

def update_usnob(infile, outfile="", outformat="reg"):

    """Take USNO-B catalog file and add V and SDSS griz magnitudes.  
       Conversions from Jordi, Grebel, and Ammon, A&A 2006."""

    try:

        # Read in catalog file
        stars = Starlist(infile)
        newstars=[]
 
        for star in stars:

            # Only want objects with good B2, R2, and I2 mags
            if (star.mags['B2MAG'] < 5.0) or (star.mags['R2MAG'] < 5.0) or \
               (star.mags['I2MAG'] < 5.0):
                continue

            # Get B, R, and I from B2, R2, and I2
            star.mags['BMAG'] = star.mags['B2MAG']
            star.mags['RMAG'] = star.mags['R2MAG']
            star.mags['IMAG'] = star.mags['I2MAG']
            del star.mags['B1MAG']
            del star.mags['B2MAG']
            del star.mags['R1MAG']
            del star.mags['R2MAG']
            del star.mags['I2MAG']

            # Get B-V from B-R fits to Landolt standards
            star.mags['VMAG'] = star.mags['BMAG'] - 0.636 * (star.mags['BMAG'] \
                                - star.mags['RMAG']) + 0.008

            # Get SDSS griz from Jordi, Grebel, and Ammon
            star.mags['GPMAG'] = star.mags['BMAG'] - 0.370 * (star.mags['BMAG']\
                                 - star.mags['VMAG']) - 0.124
            if (star.mags['VMAG'] - star.mags['RMAG'] < 0.93):
                star.mags['RPMAG'] = star.mags['RMAG'] + 0.267 * \
                                     (star.mags['VMAG'] - star.mags['RMAG']) \
                                     + 0.088
            else:
                star.mags['RPMAG'] = star.mags['RMAG'] + 0.77 * \
                                     (star.mags['VMAG'] - star.mags['RMAG']) \
                                     - 0.37
            star.mags['IPMAG'] = star.mags['IMAG'] + 0.247 * (star.mags['RMAG']\
                                 - star.mags['IMAG']) + 0.329
            star.mags['ZPMAG'] = star.mags['RPMAG'] - 1.584 * \
                                 (star.mags['RMAG'] - star.mags['IMAG']) + 0.386

            # If have U-band, can use that to get u'
            if star.mags.has_key('UMAG'):
                star.mags['UPMAG']=0.750*(star.mags['UMAG']-star.mags['BMAG'])+0.770*(star.mags['BMAG']-star.mags['VMAG'])+0.720+star.mags['GPMAG']
            # Otherwise estimate u' from Landolt standards,
            # get U from Jordi, Grebel, and Ammon
            else:
                star.mags['UPMAG']=0.596*(star.mags['GPMAG']-star.mags['RPMAG'])+star.mags['GPMAG']+1.110
                star.mags['UMAG']=1.333 * (star.mags['UPMAG'] - star.mags['GPMAG']) - 1.027 * (star.mags['BMAG'] - star.mags['VMAG']) - 0.960 + star.mags['BMAG']

            # Add to new star list
            newstars.append(star)

        # Write out results
        if (outfile==""):
            outfile=infile
        newlist=Starlist(stars=newstars) 
        newlist.write(outfile,format=outformat)

    except:
        print "Error reading starlist: %s" % infile
    
    return
        
                
