#! /usr/bin/env python
#------------------------------------------------
# This code analyze a given DES image PSF 
#------------------------------------------------

import sys, glob
sys.path.append('/usr/remote/user/sispi/jiangang/des-sv')
sys.path.append('/usr/remote/user/sispi/jiangang/decam-fermi/pyRaytrace')
sys.path.append('/usr/remote/user/sispi/jiangang/lib/python')

from psfFocus import *

def CRAYposition(beta,removeMean=True):
    """
    This codes give the suggested hexapod adjustment relative to the position the image is taken. The input is the zernike coefficients correpsonding to M20
    """
    if removeMean == True:
        knn = p.load(open('/usr/remote/user/sispi/jiangang/des-sv/finerGridKnnObj_remMean.cp','r'))
        tmean,tstd = p.load(open('/usr/remote/user/sispi/jiangang/des-sv/finerGridStdConst_remMean.cp','r'))
        beta = beta[1:]
    else:
        knn = p.load(open('/usr/remote/user/sispi/jiangang/des-sv/finerGridKnnObj.cp','r'))
        tmean,tstd = p.load(open('/usr/remote/user/sispi/jiangang/des-sv/finerGridStdConst.cp','r'))
        beta = (beta - tmean)/tstd
    CRAYParameter = knn.predict(beta)[0] #this gives the current optics status
    CRAYParameter[0] = CRAYParameter[0]*1000.
    CRAYParameter[1] = CRAY[1]*1000.
    CRAYParameter[2] = CRAY[2]*1000.
    return CRAYParameter


def hexapodPosition(beta,removeMean=True):
    """
    This code convert the CRAY position to the hexapod position parameters
    The 75 degree rotation between the two coordinate is measured by Steve. 
    """
    x,y,z,thetax,thetay = CRAYposition(beta,removeMean=removMean)
    ang = np.deg2rad(75.)
    xh = x*np.cos(ang) - y*np.sin(ang)
    yh = x*np.sin(ang) - y*np.cos(ang)
    phi = np.arctan(thetay,thetax)
    theta = np.sqrt(thetax**2+thetay**2)
    thetaxh = theta*np.cos(phi - ang)
    thetayh = theta*np.sin(phi - ang)
    zh = z
    return np.array([xh,yh,zh,thetaxh,thetayh])



if len(sys.argv) == 1:
    print 'syntax: '
    print 'desImgAnalysis expid'
    print 'The image need to be reduced (bias subtraction, flat fielding'
else:
    expid = sys.argv[1]
    img_name = 'DECam_'+expid+'_reduced.fits'
    catname = img_name[0:-5]+'_star_catalog.fits'
    if not os.path.isfile(catname):
        os.system('getstar.py '+img_name)
    imghdu = pf.open(img_name)
    cathdu = pf.open(catname)
    dimmfwhm = pf.getheader(img_name,0)['dimmsee']
    kernelSigma = np.sqrt(dimmfwhm**2+0.55**2)/2.35482
    hexposhdr = pf.getheader(img_name,0)['telfocus']
    data=[]
    stamplist=[]
    bkglist=[]
    dataSex=[]
    dataAmom=[]  # for the adaptive moments, only M20 and M22
    fwhmSex = np.array([])
    whiskerSex = np.array([])
    magall = []
    radall = []
    okall = []
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
        magall.append(mag)
        radall.append(rad)
        okall.append(ok)
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
            moms = measure_stamp_moments(stamp,bkg,kernelSigma/scale)
            momsA = measure_stamp_moments(stamp,bkg,kernelSigma/scale,adaptive=True)
            data.append([xccd,yccd]+ list(moms))
            dataAmom.append([xccd,yccd]+ list(momsA))
            dataSex.append([xccd,yccd,M20,M22])
            fwhmSex = np.concatenate((fwhmSex,fwhm_sex))
            whiskerSex = np.concatenate((whiskerSex,whisker_sex))
        else:
            continue
    data = np.array(data)
    dataSex = np.array(dataSex)
    dataAmom = np.array(dataAmom)
    magall = np.array(magall)
    radall = np.array(radall)
    okall = np.array(okall)
    display_2nd_moments(dataSex)
    pl.savefig('moments_sextractor_'+expid+'.png')
    pl.close()
    display_2nd_moments(dataAmom)
    pl.savefig('moments_adaptive_'+expid+'.png')
    pl.close()
    display_moments(data)
    pl.savefig('moments_measurement_'+expid+'.png')
    pl.close()
    fwhm_whisker_des_plot(stampImgList=stamplist,bkgList=bkglist,whkSex=whiskerSex*0.27,fwhmSex=fwhmSex*0.27,sigma=kernelSigma/scale,dimmfwhm=dimmfwhm)
    #fwhm_whisker_des_plot(stamplist,whiskerSex*0.27,fwhmSex*0.27,dimmfwhm)
    pl.savefig('fwhm_whisker_'+expid+'.png')
    pl.close()
    pl.plot(magall,radall,'b,')
    pl.plot(magall[ok],radall[ok],'r,')
    pl.ylim(0,10)
    pl.savefig('mag_radius_'+expid+'.png')
    pl.close()
    

    #---the hexapod adjustment ---
    beta,betaErr,R2_adj = zernikeFit(data[:,0].real,data[:,1].real,data[:,2].real,max_order=20)
    hexHao = hexapodPosition(beta)
    hexHaoCRAY = CRAYposition(beta)
    betaSex,betaErrSex,R2_adjSex = zernikeFit(dataSex[:,0].real,dataSex[:,1].real,dataSex[:,2].real,max_order=20)
    hexSex = hexapodPosition(betaSex)
    hexSexCRAY = CRAYposition(betaSex)
    betaA,betaErrA,R2_adjA = zernikeFit(dataAmom[:,0].real,dataAmom[:,1].real,dataAmom[:,2].real,max_order=20)
    hexA = hexapodPosition(betaA)
    hexACRAY = CRAYposition(betaA)

    print '----hexpod configuration from header -----'
    print hexposhdr

    print '--------the suggested relative hexapod adjustment -----'
    print '        ------based on weighted moments --------'
    print ' -- xShift[micron], yShift[micron], zShift[micron], xTilt[arcsec], yTilt[arcsec] --'
    print hexHao
    print '        ------based on Adaptive moments  --------'
    print ' -- xShift[micron], yShift[micron], zShift[micron], xTilt[arcsec], yTilt[arcsec] --'
    print hexA
    print '        ------based on moments from sextractor --------'
    print ' -- xShift[micron], yShift[micron], zShift[micron], xTilt[arcsec], yTilt[arcsec] --'
    print hexSex

    hexposhdr = np.array(hexposhdr.split(',')).astype(float)[0:5]
    p.dump([hexposhdr,hexHao,hexA,hexSex],open(img_name[0:-5]+'_hexpod.p','w'))
    p.dump([hexposhdr,hexHaoCRAY,hexACRAY,hexSexCRAY],open(img_name[0:-5]+'_CRAY.p','w'))
    
