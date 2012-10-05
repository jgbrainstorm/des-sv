#! /usr/bin/env python
#------------------------------------------------
# This code analyze a given DES image PSF 
#------------------------------------------------

import sys, glob
sys.path.append('/usr/remote/user/sispi/jiangang/des-sv')
sys.path.append('/usr/remote/user/sispi/jiangang/decam-fermi')
from psfFocus import *

def hexapodAdj(beta):
    """
    This codes give the suggested hexapod adjustment relative to the position the image is taken. The input is the zernike coefficients correpsonding to M20
    """
    knn = p.load(open('finerGridKnnObj.cp','r'))
    tmean,tstd = p.load(open('finerGridStdConst.cp','r'))
    beta = (bata - tmean)/tstd
    hexapodParameter = knn.predict(beta) #this gives the current optics status
    hexapodParameter[0] = hexapodParameter[0]*1000.
    hexapodParameter[1] = hexapodParameter[1]*1000.
    hexapodParameter[2] = hexapodParameter[2]*1000.
    return hexapodParameter*(-1)

if len(sys.argv) == 1:
    print 'syntax: '
    print 'desImgAnalysis expid'
    print 'The image need to be reduced'
else:
    expid = sys.argv[1]
    img_name = 'DECam_'+expid+'_reduced.fits'
    os.system('getstar.py '+img_name)
    catname = img_name[0:-5]+'_star_catalog.fits'
    imghdu = pf.open(image_name)
    cathdu = pf.open(catname)

    data=[]
    stamplist=[]
    bkglist=[]
    dataSex=[]
    fwhmSex = np.array([])
    whiskerSex = np.array([])
    #starFwhm = selectStarFwhm(catname)
    for i in range(1,63):
        print i
        img = imghdu[i].data
        cat = cathdu[i].data
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
        ok = (np.abs(fwhm_sex - starFwhm) < 0.3)*(x>100)*(x<2050)*(y>100)*(y<4100)*(flag == 0)*(mag<-12)*(mag>-14)
        if ok.any():
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
            bkglist = bkglist+list(bkg)
            xccd = eval(imghdu[i].header['extname'])[1]
            yccd = eval(imghdu[i].header['extname'])[2]
            moms = measure_stamp_moments(stamp,bkg,4.)
            data.append([xccd,yccd]+ list(moms))
            dataSex.append([xccd,yccd,M20,M22])
            fwhmSex = np.concatenate((fwhmSex,fwhm_sex))
            whiskerSex = np.concatenate((whiskerSex,whisker_sex))
        else:
            continue
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

    #---the hexapod adjustment ---
    beta,betaErr,R2_adj = zernikeFit(data[:,0].real,data[:,1].real,data[:,2].real,max_order=20)
    hexHao = hexapodAdj(beta)
    betaSex,betaErrSex,R2_adjSex = zernikeFit(dataSex[:,0].real,dataSex[:,1].real,dataSex[:,2].real,max_order=20)
    hexSex = hexapodAdj(betaSex)
    print '--------the suggested relative hexapod adjustment -----'
    print '        ------based on moments from Jiangang measurement --------'
    print ' -- xShift[micron], yShift[micron], zShift[micron], xTilt[arcsec], yTilt[arcsec] --'
    print hexHao
    print '        ------based on moments from sextractor --------'
    print ' -- xShift[micron], yShift[micron], zShift[micron], xTilt[arcsec], yTilt[arcsec] --'
    print hexSex
