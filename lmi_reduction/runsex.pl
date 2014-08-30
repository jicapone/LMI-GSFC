#!/usr/bin/env perl

&usage if ($#ARGV <1);

$snpath = $ENV{'REDUCTION'};
($exdir,$exname)=($0=~/^(.*\/)([^\/]+)$/);

$TRUE=1;
$FALSE=0;

$masktype="NONE";
$mask="check.fits";
$wttype="NONE";
$wtimage="weight.fits";
$wtcut="0";
$sat=50000;
$zp=25;
$pix=1.0;
$fwhm=1.5;
$aperture=5.0;
$minarea=5;
$gain=0;
$nobad=0;
$dotheta=$FALSE;
$doiso=$FALSE;
$doicorr=$FALSE;
$doaper=$FALSE;
$doerrs=$FALSE;

##################################################

# Parse command line
while (@ARGV) {
    $_=shift(@ARGV);
    if (/^[-\+]/) {
	# Options
	if (/^-sat/i) {
	    $sat=shift(@ARGV);
	}elsif (/^-zp/i) {
	    $zp=shift(@ARGV);
	}elsif (/^-nobad/i) {
	    $nobad=$TRUE;
	}elsif (/^\+back/i) {
	    # generate fitted background image
	    $masktype="BACKGROUND";
	    $mask=shift(@ARGV);
	}elsif (/^-back/i) {
	    # generate background-subtracted image
	    $masktype="-BACKGROUND";
	    $mask=shift(@ARGV);
	}elsif (/^-mask/i) {
	    # sky/background blanked out
	    $masktype="OBJECTS";
	    $mask=shift(@ARGV);
	}elsif (/^\+mask/i) {
	    # objects blanked out
	    $masktype="-OBJECTS";
	    $mask=shift(@ARGV);
	}elsif (/^-theta/i) {
	    $dotheta=$TRUE;
	}elsif (/^-fwhm/i) {
	    $fwhm=shift(@ARGV);
	}elsif (/^-weight/i) {
	    $wttype="MAP_WEIGHT";
	    $wtimage=shift(@ARGV);
	}elsif (/^-wcut/i) {
	    $wtcut=shift(@ARGV);
	}elsif (/^-pix/i) {
	    $pix=shift(@ARGV);
	}elsif (/^-err/i) {
	    $doerrs=$TRUE;
        }elsif (/^-gain/i) {
            $gain=shift(@ARGV);
        }elsif (/^-aperture/i) {
            $aperture=shift(@ARGV);
        }elsif (/^-minarea/i) {
            $minarea=shift(@ARGV);
        }elsif (/^-doiso/i) {
            $doiso=$TRUE;
        }elsif (/^-doicorr/i) {
            $doicorr=$TRUE;
        }elsif (/^-doaper/i) {
            $doaper=$TRUE;
	}else{
	    print STDERR "Unrecognized option $_...skipping.\n";
	}
    }else{
	# Arguments
	push(@args,$_);
    }
}

# Parse arguments
&usage if (@args<2);
($img,$sigma)=@args[0,1];

# Allow user to pass root or full name of image
if ($img=~/^(.+)\.fits$/) {
    $imgroot=$1;
    $imgname=$img;
}elsif (-e "$img.fits") {
    $imgroot=$img;
    $imgname="$img.fits";
}

##################################################

docmd("cp $snpath/Sex/daofind.param  $snpath/Sex/default.conv  $snpath/Sex/default.nnw ."); 

if ($dotheta) {
    system(qq|echo "THETA_IMAGE" >> daofind.param|);
}
if ($doiso) {
    system(qq|echo "MAG_ISO" >> daofind.param|);
}
if ($doicorr) {
    system(qq|echo "MAG_ISOCOR" >> daofind.param|);
}
if ($doaper) {
    system(qq|echo "MAG_APER" >> daofind.param|);
}
if ($doerrs) {
    system(qq|echo "MAGERR_AUTO" >> daofind.param|);
    system(qq|echo "ERRA_IMAGE" >> daofind.param|);
    system(qq|echo "ERRB_IMAGE" >> daofind.param|);
    if ($dotheta) {
        system(qq|echo "ERRTHETA_IMAGE" >> daofind.param|);
    }
    if ($doiso) {
        system(qq|echo "MAGERR_ISO" >> daofind.param|);
    }
    if ($doicorr) {
        system(qq|echo "MAGERR_ISOCOR" >> daofind.param|);
    }
    if ($doaper) {
        system(qq|echo "MAGERR_APER" >> daofind.param|);
    }
}

print "Writing default.sex file\n";
open (F,">default.sex");
print  F <<EOT;

# Default configuration file for SExtractor V1.2
# EB 18/08/97
# (*) indicates parameters which can be omitted from this config file.

#-------------------------------- Catalog ------------------------------------

CATALOG_NAME	test.cat	# name of the output catalog
CATALOG_TYPE	ASCII_HEAD	# "ASCII_HEAD","ASCII","FITS_1.0" or "FITS_LDAC"

PARAMETERS_NAME	daofind.param	# name of the file containing catalog contents

#------------------------------- Extraction ----------------------------------

DETECT_TYPE	CCD		# "CCD" or "PHOTO" (*)
DETECT_MINAREA	$minarea	# minimum number of pixels above threshold
DETECT_THRESH	$sigma		# <sigmas> or <threshold>,<ZP> in mag.arcsec-2
ANALYSIS_THRESH	$sigma		# <sigmas> or <threshold>,<ZP> in mag.arcsec-2

FILTER		Y		# apply filter for detection ("Y" or "N")?
FILTER_NAME	default.conv	# name of the file containing the filter

DEBLEND_NTHRESH	32		# Number of deblending sub-thresholds
DEBLEND_MINCONT	0.0001		# Minimum contrast parameter for deblending

CLEAN		Y		# Clean spurious detections? (Y or N)?
CLEAN_PARAM	1.0		# Cleaning efficiency

MASK_TYPE       CORRECT         # type of detection MASKing; can be one of 
                                # NONE, BLANK, or CORRECT

#------------------------------ Photometry -----------------------------------

PHOT_APERTURES	$aperture	# MAG_APER aperture diameter(s) in pixels
PHOT_AUTOPARAMS	2.5, 3.5	# MAG_AUTO parameters: <Kron_fact>,<min_radius>

SATUR_LEVEL	$sat		# level (in ADUs) at which arises saturation

MAG_ZEROPOINT	$zp		# magnitude zero-point
MAG_GAMMA	4.0		# gamma of emulsion (for photographic scans)
GAIN		$gain		# detector gain in e-/ADU.
PIXEL_SCALE	$pix		# size of pixel in arcsec (0=use FITS WCS info).

#------------------------- Star/Galaxy Separation ----------------------------

SEEING_FWHM	$fwhm		# stellar FWHM in arcsec
STARNNW_NAME	default.nnw	# Neural-Network_Weight table filename

#------------------------------ Background -----------------------------------

BACK_SIZE	64		# Background mesh: <size> or <width>,<height>
BACK_FILTERSIZE	3		# Background filter: <size> or <width>,<height>

BACKPHOTO_TYPE	GLOBAL		# can be "GLOBAL" or "LOCAL" (*)

#------------------------------ Check Image ----------------------------------

CHECKIMAGE_TYPE	$masktype	# can be one of NONE, BACKGROUND, 
                                # BACKGROUND_RMS, MINIBACKGROUND, 
                                # MINIBACK_RMS, -BACKGROUND, OBJECTS,
				# -OBJECTS, SEGMENTATION, or APERTURES

CHECKIMAGE_NAME	$mask	        # Filename for the check-image (*)

#--------------------- Memory (change with caution!) -------------------------

MEMORY_OBJSTACK	3000		# number of objects in stack
MEMORY_PIXSTACK	300000		# number of pixels in stack
MEMORY_BUFSIZE	1024		# number of lines in buffer

#----------------------------- Weight Images ---------------------------------

WEIGHT_TYPE     $wttype         # NONE, BACKGROUND, MAP_RMS, MAP_VAR, MAP_WEIGHT
WEIGHT_IMAGE    $wtimage        # Name of the weights image file
WEIGHT_THRESH   $wtcut          # Cutoff value for bad pixels

#----------------------------- Miscellaneous ---------------------------------

VERBOSE_TYPE	NORMAL		# can be "QUIET", "NORMAL" or "FULL" (*)

#------------------------------- New Stuff -----------------------------------

# Surprise!!

EOT

    close(F);

##################################################

docmd("sex $imgname");

docmd("egrep '^#' test.cat > $imgname.stars");
if ($nobad) { 
    docmd("egrep -v '^#' test.cat | sort -n -k 4 | awk \'\{ if (\$5 == 0) print\}\' >> $imgname.stars"); 
}else{
    docmd("egrep -v '^#' test.cat | sort -n -k 4 >> $imgname.stars");
}

docmd("rm -f daofind.param default.conv default.nnw default.sex test.cat");

exit(0);

######################################################################

sub docmd {
    my $cmd;
    foreach $cmd (@_) {
	print "$cmd\n";
	system($cmd);
    }
}

sub usage {
    print STDERR "Usage: $exname image sigma -sat (saturated pixels 30000) -zp (zeropoint) [-/+]mask -nobad -theta -err -pix <arcsec/pix> -fwhm <arcsec> -weight <wtimage> -wcut <wtcut> -gain <e-/ADU> -aperture <pix> -minarea <pix> <-doiso> <-doicorr> <-doaper> <-errs>\n ";
    exit(1);
}

__END__


