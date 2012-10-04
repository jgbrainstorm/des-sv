#! /usr/bin/env python
#------------------------------------------------
# This code analyze a given DES image PSF 
#------------------------------------------------

import sys, glob
sys.path.append('/usr/remote/user/sispi/jiangang/des-sv')
sys.path.append('/usr/remote/user/sispi/jiangang/decam-fermi')
from psfFocus import *


if len(sys.argv) == 1:
    print 'syntax: '
    print 'desImgAnalysis expid'
    print 'The image need to be reduced'
else:
    expid = sys.argv[1]
    img_name = 'DECam_'+expid+'.fits'
    os.system('getstar.py '+image_name)
    catname = img_name[0:-5]+'_star_catalog.fits'
    imghdu = pf.open(image_name)
    cathdu = pf.open(catname)

    data=[]
    stamplist=[]
    dataSex=[]
    fwhmSex = np.array([])
    whiskerSex = np.array([])
    starFwhm = selectStarFwhm(catname)
    for i in range(1,63):
        print i
        img = imghdu[i].data
        cat = cathdu[i].data
        x = cat.XWIN_IMAGE
        y = cat.YWIN_IMAGE
        rad = cat.FLUX_RADIUS
        mag = cat.MAG_AUTO
        flag = cat.FLAGS
        classStar = cat.CLASS_STAR
        bkg = cat.BACKGROUND
        Mcc = cat.X2WIN_IMAGE
        Mrr = cat.Y2WIN_IMAGE
        Mrc = cat.XYWIN_IMAGE
        fwhm_sex = cat.FWHM_IMAGE
        ok = np.abs(fwhm_sex - starFwhm) < 0.2
        x = x[ok]
        y = y[ok]
        bkg = bkg[ok]
        Mrr = Mrr[ok]
        Mcc = Mcc[ok]
        Mrc = Mrc[ok]
        Nobj = len(Mrr)
        M20=np.zeros(Nobj)
        M22=np.zeros(Nobj).astype(complex)
        for k in range(Nobj):
            M20[k] = Mrr[k] + Mcc[k]
            M22[k] = np.complex(Mcc[k] - Mrr[k],2*Mrc[k])
        whisker_sex = np.sqrt(np.abs(M22))
        fwhm_sex = fwhm_sex[ok]
        M20 = np.median(M20)
        M22 = np.median(M22)
        stamp = getStamp(data=img,xcoord=x,ycoord=y,Npix=25)
        stamplist = stamplist+stamp
        xccd = eval(imghdu[i].header['extname'])[1]
        yccd = eval(imghdu[i].header['extname'])[2]
        moms = measure_stamp_moments(stamp,bkg,4.)
        data.append([xccd,yccd]+ list(moms))
        dataSex.append([xccd,yccd,M20,M22])
        fwhmSex = np.concatenate((fwhmSex,fwhm_sex))
        whiskerSex = np.concatenate((whiskerSex,whisker_sex))


    data = np.array(data)
    dataSex = np.array(dataSex)
    display_2nd_moments(dataSex)
    pl.savefig('moments_sextractor_'+expid+'.png')
    pl.close()
    display_moments(data)
    pl.savefig('moments_measurement_'+expid+'.png')
    pl.close()
    fwhm_whisker_des_plot(stamplist,whiskerSex*0.27,fwhmSex*0.27)
    pl.savefig('fwhm_whisker_'+expid+'.png')
    pl.close()
