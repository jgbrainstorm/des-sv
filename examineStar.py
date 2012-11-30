#! /usr/bin/env python
#------------------------------------------------
# This code is examine the stars of a given exposure and extension
# J. Hao, 11/20/2012 @ FNAL 
#------------------------------------------------
import sys
import glob as gl
import cPickle as p
sys.path.append('/usr/remote/user/sispi/jiangang/des-sv')

from decamImgAnalyzer_def import *



def star_viewer(img_name=None,ext=None):
    pl.ion()
    catname = img_name[0:-5]+'_star_catalog.fits'
    expid = img_name[6:14]
    img = pf.getdata(img_name,ext)
    cat = pf.getdata(catname,ext)
    detector = pf.getheader(img_name,ext)['detpos']
    exptime = pf.getheader(img_name,0)['exptime']
    x = cat.XWIN_IMAGE
    y = cat.YWIN_IMAGE
    rad = cat.FLUX_RADIUS
    mag = cat.MAG_AUTO
    flag = cat.FLAGS
    bkg = cat.BACKGROUND
    Mcc = cat.X2WIN_IMAGE
    Mrr = cat.Y2WIN_IMAGE
    Mrc = cat.XYWIN_IMAGE
    fwhm_sex = cat.FWHM_IMAGE   
    starFwhm = selectStar(mag,fwhm_sex)
    ok = (np.abs(fwhm_sex - starFwhm) < 0.4)*(x>100)*(x<2050)*(y>100)*(y<4100)*(flag == 0)*(mag<=-11.5)*(mag>-14.5)
    nstar = len(mag[ok])
    if ok.any():
        print '--- Nstars selected: '+str(nstar)+'---'
        bkg = bkg[ok]
        Mrr = Mrr[ok]
        Mcc = Mcc[ok]
        Mrc = Mrc[ok]
        x=x[ok]
        y=y[ok]
        stamplist = getStamp(data=img,xcoord=x,ycoord=y,Npix=25)
        bkglist = list(bkg)
        dispStampList(stamplist,bkglist,2.,mag,rad,ok,expid,detector,exptime)
    else:
        print '--- No Nstars selected ---'
        pl.ion()
        pl.figure(figsize=(5,5))
        pl.plot(mag,rad,'b.')
        pl.xlabel('mag')
        pl.ylabel('radius')
        pl.ylim(0,13)
        pl.xlim(-20,0)
        pl.title('Exposure: '+expid+'   CCD: '+detector)
    return '----finished one image ----'
    

if __name__ == "__main__":
    from examineStar import *
    import sys,time,glob
    startTime=time.time()
    if len(sys.argv) == 1:
        print 'syntax: '
        print 'examineStar expid ext'
        print 'Note: The image need to be reduced (bias subtraction, flat fielding'
    else:
        expid = sys.argv[1]
        ext = int(sys.argv[2])
        img_name = 'DECam_'+expid+'_reduced.fits'
        t=star_viewer(img_name,ext)
    

