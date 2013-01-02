#!/usr/bin/env python

"""
Code to test DES Science Requirement R-19 and R-20:

R-19:

The mean PSF whisker length for stars per exposure must be below 0.2" in the
r, i, and z bands for the wide-area survey.

R-20:

For the wide-area survey, the residual mean whisker length for stars on scales
of 10 arcmin to 1 degree, after removal of a static component (i.e., the same
for all exposures) and a bilinear fit in (x,y) per exposure, should be below
0.06" in r, i, and z bands.  Residual mean PSF whisker length on scales of
10' - 1deg in r, i, z < 0.06"
"""

import pyfits
import optparse
import sys
import logging
import math
import numpy
import glob

# Various hard-coded values here:
nchips = 62

# File_names are assumed to be DECam_$(expnum)_$(chipnum)_cat.fits
# expnum is given on the command line
# chipnum ranges from 1..nchips
filename_prefix = 'DECam_'
filename_suffix = '_cat.fits'

script_name = 'whisker.py'

# The input catalog is in the 3rd hdu (numbered extension 2) out of a total of 3.
# If this changes, change these values appropriately.
tot_hdu = 3
cat_hdu = 2

# The column names in the input catalog that we use:
ra_col = 'ALPHA_J2000'
dec_col = 'DELTA_J2000'
mag_col = 'MAG_AUTO'
r50_col = 'FLUX_RADIUS'
flags_col = 'FLAGS'
ixx_col = 'X2_WORLD'
ixy_col = 'XY_WORLD'
iyy_col = 'Y2_WORLD'
sigixx_col = 'ERRX2_WORLD'
sigixy_col = 'ERRXY_WORLD'
sigiyy_col = 'ERRY2_WORLD'
#sg_col = 'CLASS_STAR'

# Paremeters to use for star selection:
#sg_minval = 0.9
mag_minval = 11
mag_maxval = 14

# The maximum uncertainty on our whisker length to accept an object.
sigma_maxval = 0.01

# All outlier rejection steps are done in terms of multiples of quartiles deviation
# from the median value.  e.g. we do outlier rejection of the log(ixx+iyy) size values
# and also the whisker length for the star selection.  We also do outlier rejection after
# the bilinear fit step based on the residual ixx,ixy,iyy values.
# This parameter specifies how many quartiles from the median is declared an outlier.
nquart = 3

# Skip the following chip_num's
skip = [ 61 ]

# Minimum number of stars to use a given chip
min_nstars = 80

def parse_command_line(argv):

    # Build the parser and add arguments
    description = """
This script computes a number of whisker length statistics from the stars in
FirstCut catalogs for a single exposure.
The results are output in table form to the specified output file.
A more human-readable form of the output is by default also output to stdout.
You can turn this off by running with '-v 0'.  
"""
    usage = """
%s root_dir exp_num out_file
    root_dir is the directory with the FirstCut catalogs.
    exp_num is the exposure number to process.
    out_file is the output file to produce."""%(script_name)
    parser = optparse.OptionParser(usage=usage, description=description)

    # optparse only allows string choices, so take verbosity as a string and make it int later
    parser.add_option(
            '-v', '--verbosity', type="choice", action='store', choices=('0', '1', '2', '3'),
            default='1', help='integer verbosity level: min=0, max=3 [default=1]')
    parser.add_option(
            '-l', '--log_file', type=str, action='store', default=None,
            help='filename for storing logging output [default is to stream to stdout]')
    parser.add_option(
            '-a', '--append', action='store_true', default=False,
            help='append values to out_file [default is overwrite any existing file]')
    (args, posargs) = parser.parse_args()

    # Parse the positional arguments by hand
    nposargs = 3
    if len(posargs) == nposargs:
        root_dir = posargs[0]
        exp_num = int(posargs[1])
        out_file = posargs[2]
    else:
        parser.print_help()
        if len(posargs) < nposargs:
            sys.exit('\n%s: error: too few arguments\n'%(script_name))
        else:
            argstring = posargs[nposargs]
            for addme in posargs[nposargs+1:]:
                argstring = argstring+' '+addme
            sys.exit('\n%s: error: unrecognised arguments: %s\n'%(script_name,argstring))

    # Parse the integer verbosity level from the command line args into a logging_level string
    logging_levels = { 0: logging.CRITICAL,
                       1: logging.WARNING,
                       2: logging.INFO,
                       3: logging.DEBUG }
    logging_level = logging_levels[int(args.verbosity)]

    # Setup logging to go to sys.stdout or (if requested) to an output file
    if args.log_file is None:
        logging.basicConfig(format="%(message)s", level=logging_level, stream=sys.stdout)
    else:
        logging.basicConfig(format="%(message)s", level=logging_level, filename=args.log_file)
    logger = logging.getLogger(script_name)

    logger.info('root_dir = %s',root_dir)
    logger.info('exp_num = %s',exp_num)
    logger.info('out_file = %s',out_file)

    if args.append:
        out = open(out_file,'a')
    else:
        out = open(out_file,'w')
        out.write('# expnum  chip    nstar     <ixx>       <ixy>       <iyy>        WL        RMS WL  RMS WL after fit\n')
        out.write('# chip 0 = Full exposure\n')
        out.write('# chip -1 = Use mean moments for each chip as 62 "stars"\n')
        out.write('# chip -2 = Full exposure after subtracting chip-wise bilinear fits\n')

    return root_dir, exp_num, out, logger
    
def get_stars(cat, logger):

    cols = cat.columns
    data = cat.data
    logger.debug('cols = %s',str(cols))
    logger.debug('len(data) = %d',len(data))

    # Star selection using CLASS_STAR doesn't really work.
    #sg = numpy.array(data.field(sg_col))
    #ok = (sg<sg_minval)

    # Star selection for first cut recommended by Jiangang
    mag = numpy.array(data.field(mag_col))
    flags = numpy.array(data.field(flags_col))
    ok = (mag>=mag_minval)*(mag<=mag_maxval)*(flags==0)

    r50 = numpy.array(data.field(r50_col))
    r50median = numpy.median(r50[ok])
    idx = (mag>=11)*(mag<=14)*(flags==0)*(abs(r50-r50median)<=0.5)

    # Input ra,dec values are in degrees.  Convert to radians.
    ra = numpy.array(data.field(ra_col))[idx] * numpy.pi/180.
    dec = numpy.array(data.field(dec_col))[idx] * numpy.pi/180.

    # Input moment values are in degrees^2.  Convert to arcsec^2.
    ixx = numpy.array(data.field(ixx_col))[idx] * 3600.**2
    ixy = numpy.array(data.field(ixy_col))[idx] * 3600.**2
    iyy = numpy.array(data.field(iyy_col))[idx] * 3600.**2

    # Also reject if the moment errors are too large.
    sigixx = numpy.array(data.field(sigixx_col))[idx] * 3600.**2
    sigixy = numpy.array(data.field(sigixy_col))[idx] * 3600.**2
    sigiyy = numpy.array(data.field(sigiyy_col))[idx] * 3600.**2

    # WL = ( (ixx-iyy)^2 + (2ixy)^2 )^1/4
    # sigma_WL = 1/2 WL^-3 
    #            [ (ixx-iyy)^2 (sigma_ixx^2 + sigma_iyy^2) + 16 ixy^2 sigma_ixy^2 ]^1/2
    wl = ( (ixx-iyy)**2 + (2.*ixy)**2 )**0.25
    sigma_wl = 0.5 * ( (ixx-iyy)**2 * sigixx**2 + 16. * ixy**2 * sigixy**2 )**0.5 / wl**3

    ok = (sigma_wl <= sigma_maxval)
    ra = ra[ok]
    dec = dec[ok]
    ixx = ixx[ok]
    ixy = ixy[ok]
    iyy = iyy[ok]

    # Now remove outlier stars according to their log(size) and log(wl)
    # I use the median value and a 3 * the quartile deviation as my clipping value.
    size = numpy.log(ixx+iyy)
    wl = ( (ixx-iyy)**2 + (2.*ixy)**2 )**0.25 
    n = len(ixx)
    logger.info('Found %d stars on first pass.',n)
    nclip = n
    while nclip > 0:
        sorted_size = sorted(size)
        size_median = sorted_size[n/2]
        size_quartile = 0.5*(sorted_size[n*3/4] - sorted_size[n/4])
        logger.info('median ixx+iyy = %f, 1-quartile deviation in log = %f',
                    numpy.exp(size_median),size_quartile)
    
        sorted_wl = sorted(wl)
        wl_median = sorted_wl[n/2]
        wl_quartile = 0.5*(sorted_wl[n*3/4] - sorted_wl[n/4])
        logger.info('median WL = %f, 1-quartile deviation = %f',wl_median,wl_quartile)
 
        ok = ( (numpy.fabs(size - size_median) < nquart*size_quartile) *
               (numpy.fabs(wl - wl_median) < nquart*wl_quartile) ) 

        ra = ra[ok]
        dec = dec[ok]
        ixx = ixx[ok]
        ixy = ixy[ok]
        iyy = iyy[ok]
        size = size[ok]
        wl = wl[ok]

        # Update the number of objects
        nclip = n-len(ixx)
        n = len(ixx)
        logger.info('clipped out %d objects.  Now n = %d',nclip,n)

    return ra, dec, ixx, ixy, iyy

def project(ra, dec, logger):

    # First construct the position of each object on the unit sphere.
    x = numpy.cos(dec) * numpy.cos(ra)
    y = numpy.cos(dec) * numpy.sin(ra)
    z = numpy.sin(dec)

    # Find the center.  This center avoids any problems with ra wrapping from 360 to 0.
    xcen = x.mean()
    ycen = y.mean()
    zcen = z.mean()
    # Renormalize back to the surface of the sphere.
    r = (xcen**2 + ycen**2 + zcen**2)**0.5
    xcen /= r
    ycen /= r
    zcen /= r
    logger.info('center of image is at (x,y,z) = %f,%f,%f',xcen,ycen,zcen)

    # Convert the positions to stereographic projections around the center point.
    # The equations are given at:
    #     http://mathworld.wolfram.com/StereographicProjection.html
    # u = k cos(dec) sin(ra-ra0)
    # v = k ( cos(dec0) sin(dec) - sin(dec0) cos(dec) cos(ra-ra0) )
    # k = 2 ( 1 + sin(dec0) sin(dec) + cos(dec0) cos(dec) cos(ra-ra0) )^-1
    # 
    # Using our sphere coords:
    #   x = cos(dec) * cos(ra))
    #   y = cos(dec) * sin(ra))
    #   z = sin(dec))
    # this becomes:
    # k = 2 ( 1 + x x0 + y y0 + z z0 )^-1
    # u = k (y x0 - x y0) / sqrt(1-z0^2)
    # v = k ( z - z0 (x x0 + y y0 + z z0) ) / sqrt(1-z0^2)
    cosdec0 = numpy.sqrt(1.-zcen**2)

    dot = x*xcen + y*ycen + z*zcen
    k = 2. / (1. + dot)
    u = k * (y*xcen - x*ycen) / cosdec0
    v = k * (z - zcen*dot) / cosdec0

    return u, v


def print_stats(ra, dec, ixx, ixy, iyy, exp_num, chip_num, out, logger):

    wl = ((ixx-iyy)**2 + (2.*ixy)**2 )**0.25
    theta = numpy.arctan2( 2.*ixy, ixx-iyy )
    wlx = wl * numpy.cos(theta)
    wly = wl * numpy.sin(theta)
    wl = (wlx**2 + wly**2)**0.5

    mean_ixx = ixx.mean()
    mean_ixy = ixy.mean()
    mean_iyy = iyy.mean()
    wl_meanmom = ((mean_ixx-mean_iyy)**2 + (2.*mean_ixy)**2 )**0.25
    theta_meanmom = numpy.arctan2( 2.*mean_ixy, mean_ixx-mean_iyy )
    wlx_meanmom = wl_meanmom * numpy.cos(theta_meanmom)
    wly_meanmom = wl_meanmom * numpy.sin(theta_meanmom)
    rms_wl_meanmom = ( ( (ixx-mean_ixx-iyy+mean_iyy)**2 + 4.*(ixy-mean_ixy)**2 ).mean() )**0.25
    mean_wlx = wlx.mean()
    mean_wly = wly.mean()
    mean_wl = wl.mean()
    rms_wl = numpy.sqrt( ((wlx-mean_wlx)**2 + (wly-mean_wly)**2).mean() )
    logger.warn('    Number of stars = %d',len(ixx))
    logger.warn('    Mean moments: <ixx> = %f, <ixy> = %f, <iyy> = %f',
                mean_ixx,mean_ixy,mean_iyy)
    logger.warn('    WL from mean moments = %f, theta = %f rad',wl_meanmom,theta_meanmom)
    logger.warn('    In cartesian coordinates: (%f,%f)',wlx_meanmom,wly_meanmom)
    logger.warn('    RMS WL from mean moments = %f',rms_wl_meanmom)
              
    logger.warn('    Mean WL = (%f,%f)',mean_wlx,mean_wly)
    logger.warn('    |Mean WL| = %f',(mean_wlx**2+mean_wly**2)**0.5)
    logger.warn('    Mean |WL| = %f',mean_wl)
    logger.warn('    RMS WL = %f',rms_wl)

    # Remove a bilinear fit:
    # wlx = a + bx + cy
    # wly = d + ex + fy
    # This can be expressed as a matrix:
    #
    # ( 1  x  y ) ( a  d ) = ( wlx wly )
    #             ( b  e )  
    #             ( c  f )
    # 
    # For the maximum likelihood fit, we just make many rows of (1 x y) and (wlx wly)
    # and solve for the fit matrix using QR decomposition.
    # We calculate the fit values using numpy for the matrix calculations.
    # We use A = ( 1  x_0  y_0 )
    #            ( 1  x_1  y_1 )
    #            ( ...         )
    #        M = ( a  d )
    #            ( b  e )
    #            ( c  f )
    #        W = ( wlx_0  wly_0 )
    #            ( wlx_1  wly_1 )
    #            ( ...          )
    x, y = project(ra, dec, logger)
    n = len(x)
    A = numpy.ones( shape=(n,3) )
    A[:,1] = x
    A[:,2] = y
    W = numpy.ones( shape=(n,5) )
    W[:,0] = wlx
    W[:,1] = wly
    W[:,2] = ixx
    W[:,3] = ixy
    W[:,4] = iyy

    # The version with all values, including outliers.
    Afull = A
    Wfull = W 

    nclip = n
    while nclip > 0:
        (M, resids, rank, s) = numpy.linalg.lstsq(A,W)

        logger.debug('M = %s',str(M))
        logger.debug('resids = %s',str(resids))
        logger.debug('rank = %d',rank)
        logger.debug('singular values = %s',str(s))
        logger.info('Bilinear fits:')
        logger.info('  WL_x = %f + (%f/deg) x + (%f/deg) y',
                    M[0,0],M[1,0]*numpy.pi/180.,M[2,0]*numpy.pi/180.)
        logger.info('  WL_y = %f + (%f/deg) x + (%f/deg) y',
                    M[0,1],M[1,1]*numpy.pi/180.,M[2,1]*numpy.pi/180.)
        logger.info('  Ixx = %f + (%f/deg) x + (%f/deg) y',
                    M[0,2],M[1,2]*numpy.pi/180.,M[2,2]*numpy.pi/180.)
        logger.info('  Ixy = %f + (%f/deg) x + (%f/deg) y',
                    M[0,3],M[1,3]*numpy.pi/180.,M[2,3]*numpy.pi/180.)
        logger.info('  Iyy = %f + (%f/deg) x + (%f/deg) y',
                    M[0,4],M[1,4]*numpy.pi/180.,M[2,4]*numpy.pi/180.)
        dW = W - numpy.dot(A,M)
        
        dwlx = dW[:,0]
        dwly = dW[:,1]
        dixx = dW[:,2]
        dixy = dW[:,3]
        diyy = dW[:,4]

        # Look for outliers in the residual moment values:
        sorted_dixx = sorted(dixx)
        sorted_dixy = sorted(dixy)
        sorted_diyy = sorted(diyy)
        dixx_median = sorted_dixx[n/2]
        dixx_quartile = 0.5*(sorted_dixx[n*3/4] - sorted_dixx[n/4])
        logger.info('median dixx = %f, 1-quartile deviation = %f',dixx_median,dixx_quartile)
        dixy_median = sorted_dixy[n/2]
        dixy_quartile = 0.5*(sorted_dixy[n*3/4] - sorted_dixy[n/4])
        logger.info('median dixy = %f, 1-quartile deviation = %f',dixy_median,dixy_quartile)
        diyy_median = sorted_diyy[n/2]
        diyy_quartile = 0.5*(sorted_diyy[n*3/4] - sorted_diyy[n/4])
        logger.info('median diyy = %f, 1-quartile deviation = %f',diyy_median,diyy_quartile)
        ok = ( (numpy.fabs(dixx - dixx_median) < nquart*dixx_quartile) *
               (numpy.fabs(dixy - dixy_median) < nquart*dixy_quartile) *
               (numpy.fabs(diyy - diyy_median) < nquart*diyy_quartile) ) 

        W = W[ok]
        A = A[ok]

        # Update the number of objects
        nclip = n-len(dwlx)
        n = len(dwlx)
        logger.info('clipped %d objects with large residuals.  Now n = %d',nclip,n)

    dwl = (dwlx**2 + dwly**2)**0.5

    mean_dwlx = dwlx.mean()
    mean_dwly = dwly.mean()
    mean_dwl = dwl.mean()
    mean_dixx = dixx.mean()
    mean_dixy = dixy.mean()
    mean_diyy = diyy.mean()
    rms_dwl = numpy.sqrt(((dwlx-mean_dwlx)**2 + (dwly-mean_dwly)**2).mean())
    rms_wl_dmeanmom = (((dixx-mean_dixx-diyy+mean_diyy)**2 + 4.*(dixy-mean_dixy)**2).mean())**0.25
    logger.debug('mean_dwlx = %f  (should = 0)',mean_dwlx)
    logger.debug('mean_dwly = %f  (should = 0)',mean_dwly)
    logger.debug('mean_dixx = %f  (should = 0)',mean_dixx)
    logger.debug('mean_dixy = %f  (should = 0)',mean_dixy)
    logger.debug('mean_diyy = %f  (should = 0)',mean_diyy)
    logger.warn('  After removing bilinear fit:')
    logger.warn('    Mean |WL| = %f',mean_dwl)
    logger.warn('    RMS WL = %f',rms_dwl)
    logger.warn('    RMS WL from mean moments = %f',rms_wl_dmeanmom)

    out.write('%8d   %2d   %6d  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f\n'%(
        exp_num, chip_num, len(ixx),
        mean_ixx, mean_ixy, mean_iyy,
        wl_meanmom, rms_wl_meanmom, rms_wl_dmeanmom))

    dW = Wfull - numpy.dot(Afull,M)
        
    dwlx = dW[:,0]
    dwly = dW[:,1]
    dixx = dW[:,2]
    dixy = dW[:,3]
    diyy = dW[:,4]

    return (dwlx, dwly, dixx, dixy, diyy)

def main(argv):

    (root_dir, exp_num, out, logger) = parse_command_line(argv)

    all_ra = numpy.array([], dtype=float)
    all_dec = numpy.array([], dtype=float)
    all_ixx = numpy.array([], dtype=float)
    all_ixy = numpy.array([], dtype=float)
    all_iyy = numpy.array([], dtype=float)

    # These will be lists of numpy arrays.
    chip_ra = []
    chip_dec = []
    chip_ixx = []
    chip_ixy = []
    chip_iyy = []

    # These are the values after removing a bilinear fit
    all_dwlx = numpy.array([], dtype=float)
    all_dwly = numpy.array([], dtype=float)
    all_dixx = numpy.array([], dtype=float)
    all_dixy = numpy.array([], dtype=float)
    all_diyy = numpy.array([], dtype=float)

    for chip_num in range(1,nchips+1):
        if chip_num in skip: 
            logger.info('Skipping chip_num %d because in skip list',chip_num)
            continue;
        filename_pattern = "%s/*%08d*%02d*%s"%(root_dir, exp_num, chip_num, filename_suffix)
        filename = glob.glob(filename_pattern)
        if (len(filename) == 0):
            logger.warn('Unable to find appropriate file for exp %d, chip %d',exp_num,chip_num)
            logger.warn('Expected something of the form: %s',filename_pattern)
            raise RuntimeError('Missing input file')
        if (len(filename) != 1):
            logger.warn('Filename pattern for exp %d, chip %d is not unique',exp_num,chip_num)
            logger.warn('Found: ')
            for name in filename:
                logger.warn('    %s',name)
            raise RuntimeError('Ambiguous input filename')
        filename = filename[0]
        logger.info('filename = %s',filename)
        hdulist = pyfits.open(filename)
        if not hdulist:
            logger.warn('Error opening input file %s',filename)
            raise RuntimeError('Invalid input file')
        if len(hdulist) != tot_hdu:
            logger.warn('Expecting %d hdus.  Found %d.',tot_hdu,len(hdulist))
            raise RuntimeError('Invalid input file')
        cat = hdulist[cat_hdu]
        logger.debug('  %s',str(cat.header))
        (ra, dec, ixx, ixy, iyy) = get_stars(cat, logger)
        logger.info('Got %d stars for chip %d',len(ra),chip_num)

        if len(ra) < min_nstars:
            logger.info('Skipping chip_num %d because too few stars',chip_num)

        all_ra = numpy.append(all_ra,ra)
        all_dec = numpy.append(all_dec,dec)
        all_ixx = numpy.append(all_ixx,ixx)
        all_ixy = numpy.append(all_ixy,ixy)
        all_iyy = numpy.append(all_iyy,iyy)

        chip_ra += [ ra ]
        chip_dec += [ dec ]
        chip_ixx += [ ixx ]
        chip_ixy += [ ixy ]
        chip_iyy += [ iyy ]

        logger.warn('Whisker stats for chip %d:',chip_num)
        dwlx, dwly, dixx, dixy, diyy = print_stats(
                ra, dec, ixx, ixy, iyy, exp_num, chip_num, out, logger)
        all_dwlx = numpy.append(all_dwlx,dwlx)
        all_dwly = numpy.append(all_dwly,dwly)
        all_dixx = numpy.append(all_dixx,dixx)
        all_dixy = numpy.append(all_dixy,dixy)
        all_diyy = numpy.append(all_diyy,diyy)

    logger.warn('Overall whisker stats:')
    print_stats(all_ra, all_dec, all_ixx, all_ixy, all_iyy, exp_num, 0, out, logger)

    logger.warn('Exposure whisker stats using chip-wise averages:')
    chipmean_ra = numpy.array([ c.mean() for c in chip_ra ])
    chipmean_dec = numpy.array([ c.mean() for c in chip_dec ])
    chipmean_ixx = numpy.array([ c.mean() for c in chip_ixx ])
    chipmean_ixy = numpy.array([ c.mean() for c in chip_ixy ])
    chipmean_iyy = numpy.array([ c.mean() for c in chip_iyy ])
    print_stats(chipmean_ra, chipmean_dec, chipmean_ixx, chipmean_ixy, chipmean_iyy, 
                exp_num, -1, out, logger)

    logger.warn('Overall whisker stats after chip-wise bilinear moment fits:')
    print_stats(all_ra, all_dec, all_dixx, all_dixy, all_diyy, exp_num, -2, out, logger)

if __name__ == "__main__":
    main(sys.argv)