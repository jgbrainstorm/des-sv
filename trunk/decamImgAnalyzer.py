#! /usr/bin/env python
#------------------------------------------------
# This code is the replacement of desImgAnalysis.py
# J. Hao, 11/13/2012 @ FNAL 
#------------------------------------------------
import sys
sys.path.append('/usr/remote/user/sispi/jiangang/des-sv')

from decamImgAnalyzer_def import *

def hexapodPosition(beta):
    """
    the CRAY position to the hexapod position parameters. There is a 15 deg rotation between the two coordinate. However, this is accounted in the sispi. So, the hexapod position in the header is acutally rotated to its designed coordiante, which relate to the CRAY coordinate by the last page of des-docdb #6551
    """
    x,y,z,thetax,thetay = CRAYposLinearModel(b=beta)
    xh = x
    yh = -y
    zh = -z
    thetaxh = -thetay
    thetayh = -thetax
    return np.array([xh,yh,zh,thetaxh,thetayh])

def CRAYposLinearModel(b=None):
    """
    here, the convention for b is: M20 (0-19), M22real(20 - 39),M22imag (40-59) 
    """
    M22realTrefoil2 = b[29] 
    M22imagTrefoil1 = b[48] 
    M22TrefoilXshift = 0.5*(M22realTrefoil2+M22imagTrefoil1) # for x decenter
    M22realTrefoil1 = b[26] 
    M22imagTrefoil2 = b[49] 
    M22TrefoilYshift = 0.5*(M22realTrefoil1 - M22imagTrefoil2) # for y decenter
    M20defocus = b[4]    # for defocus
    M22realComa2 = b[28] 
    M22imagComa1 = b[47]
    M22ComaXtilt = 0.5*(M22realComa2+M22imagComa1) # for x-tilt
    M22realComa1 = b[27] 
    M22imagComa2 = b[48]
    M22ComaYtilt = 0.5*(M22realComa1 - M22imagComa2) # for y-tilt
    x = -3.0063 * M22TrefoilXshift -0.0053
    y = -2.9318 * M22TrefoilYshift - 0.0005
    z = 0.4046 * M20defocus - 0.0705
    xtilt = 1075.8892 * M22ComaXtilt - 0.4876
    ytilt = -1064.6332 * M22ComaYtilt - 0.1234
    return x*1000,y*1000,z*1000,xtilt,ytilt


def runanalysis(img_name=None):
    catname = img_name[0:-5]+'_star_catalog.fits'
    if not os.path.isfile(catname):
        os.system('getstar.py '+img_name)
    imghdu = pf.open(img_name)
    cathdu = pf.open(catname)
    expid = img_name[6:14]
    dimmfwhm = pf.getheader(img_name,0)['dimmsee']
    kernelSigma = np.sqrt(dimmfwhm**2+0.55**2)/2.35482
    hexposhdr = pf.getheader(img_name,0)['telfocus']
    bcamDX = pf.getheader(img_name,0)['BECAMDX']
    bcamDY = pf.getheader(img_name,0)['BECAMDY']
    bcamAX = pf.getheader(img_name,0)['BECAMAX']
    bcamAY = pf.getheader(img_name,0)['BECAMAY']
    data=[]
    stamplist=[]
    bkglist=[]
    dataSex=[]
    fwhmSex = np.array([])
    whiskerSex = np.array([])
    magall = []
    radall = []
    okall = []
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
        ok = (np.abs(fwhm_sex - starFwhm) < 0.4)*(x>100)*(x<2050)*(y>100)*(y<4100)*(flag == 0)*(mag<-11.5)*(mag>-15)*(flag ==0)
        nstar = len(mag[ok])
        print '--- Nstars selected: '+str(nstar)+'---'
        magall.append(mag)
        radall.append(rad)
        okall.append(ok)
        if ok.any():
            bkg = bkg[ok]
            Mrr = Mrr[ok]
            Mcc = Mcc[ok]
            Mrc = Mrc[ok]
            x=x[ok]
            y=y[ok]
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
            stamp = getStamp(data=img,xcoord=x,ycoord=y,Npix=30)
            stamplist = stamplist+stamp
            bkglist = bkglist+list(bkg)
            xccd = eval(imghdu[i].header['extname'])[1]
            yccd = eval(imghdu[i].header['extname'])[2]
            moms = measure_stamp_moments(stamp,bkg,2.)
            data.append([xccd,yccd]+ list(moms))
            dataSex.append([xccd,yccd,M20,M22])
            fwhmSex = np.concatenate((fwhmSex,fwhm_sex))
            whiskerSex = np.concatenate((whiskerSex,whisker_sex))
        else:
            continue
    data = np.array(data)
    dataSex = np.array(dataSex)
    magall = np.array(magall)
    radall = np.array(radall)
    okall = np.array(okall)
    #display_2nd_moments(dataSex)
    #pl.savefig('moments_sextractor_'+expid+'.png')
    #pl.close()
    #display_moments(data)
    display_2nd_moments(data)
    pl.savefig('moments_whisker_'+expid+'.png')
    pl.close()
    fwhm_whisker_des_plot(stampImgList=stamplist,bkgList=bkglist,whkSex=whiskerSex*0.27,fwhmSex=fwhmSex*0.27,sigma=kernelSigma/scale,dimmfwhm=dimmfwhm)
    pl.savefig('fwhm_whisker_'+expid+'.png')
    pl.close()
    #---check the fitted value of the moments ---
    #datafitted = data.copy()
    #datafitted[:,2].real = zernikeFit(data[:,0].real,data[:,1].real,data[:,2].real,max_order=20)[3]
    #datafitted[:,3].real = zernikeFit(data[:,0].real,data[:,1].real,data[:,3].real,max_order=20)[3]
    #datafitted[:,3].imag = zernikeFit(data[:,0].real,data[:,1].real,data[:,3].imag,max_order=20)[3]
    #datafitted[:,4].real = zernikeFit(data[:,0].real,data[:,1].real,data[:,4].real,max_order=20)[3]
    #datafitted[:,4].imag = zernikeFit(data[:,0].real,data[:,1].real,data[:,4].imag,max_order=20)[3]
    #datafitted[:,5].real = zernikeFit(data[:,0].real,data[:,1].real,data[:,5].real,max_order=20)[3]
    #datafitted[:,5].imag = zernikeFit(data[:,0].real,data[:,1].real,data[:,5].imag,max_order=20)[3]
    #display_moments(datafitted)
    #pl.savefig('fitted_moments_measurement_'+expid+'.png')
    #pl.close()
    #---the hexapod adjustment using M20,M22---
    beta=[]
    betaErr=[]
    betaM20 = zernikeFit(data[:,0].real,data[:,1].real,data[:,2].real,max_order=20)
    beta.append(betaM20[0])
    betaErr.append(betaM20[1])
    betaM22real = zernikeFit(data[:,0].real,data[:,1].real,data[:,2].real,max_order=20)
    beta.append(betaM22real[0])
    betaErr.append(betaM22real[1])

    betaM22imag = zernikeFit(data[:,0].real,data[:,1].real,data[:,3].imag,max_order=20)
    beta.append(betaM22imag[0])
    betaErr.append(betaM22imag[1])
    betaforplot = beta
    beta=np.array(beta)
    beta=beta.flatten()
    posCRAY = CRAYposLinearModel(beta)
    hexHao = hexapodPosition(beta)
    dispM202Coeff(betaAll = betaforplot, betaErrAll = betaErr,hexinfo=hexHao)
    pl.savefig('zernike_coeff_'+expid+'.png')
    pl.close()
 
    #betaSex=[]
    #betaSex.append(zernikeFit(dataSex[:,0].real,dataSex[:,1].real,dataSex[:,2].real,max_order=20)[0])
    #betaSex.append(zernikeFit(dataSex[:,0].real,dataSex[:,1].real,dataSex[:,3].real,max_order=20)[0])
    #betaSex.append(zernikeFit(dataSex[:,0].real,dataSex[:,1].real,dataSex[:,3].imag,max_order=20)[0])
    #dispM202Coeff(betaSex)
    #pl.savefig('coeff_Sex_'+expid+'.png')
    #pl.close()
    #betaSex = np.array(betaSex)
    #betaSex=betaSex.flatten()
    #posCRAYsex = CRAYposLinearModel(betaSex)
    #hexSex = hexapodPosition(betaSex)
    #print '----hexpod configuration from header -----'
    #print hexposhdr
    #print '--------the hexapod positions -----'
    #print '        ------based on weighted moments --------'
    #print ' -- xShift[micron], yShift[micron], zShift[micron], xTilt[arcsec], yTilt[arcsec] --'
    #print hexHao
    #print '        ------based on moments from sextractor --------'
    #print ' -- xShift[micron], yShift[micron], zShift[micron], xTilt[arcsec], yTilt[arcsec] --'
    #print hexSex
    #print ' '
    #print ' '
    #print ' '
    #print ' '
    #print '--------the CRAY positions -----'
    #print '        ------based on weighted moments --------'
    #print ' -- xShift[micron], yShift[micron], zShift[micron], xTilt[arcsec], yTilt[arcsec] --'
    #print posCRAY
    #print '        ------based on moments from sextractor --------'
    #print ' -- xShift[micron], yShift[micron], zShift[micron], xTilt[arcsec], yTilt[arcsec] --'
    #print posCRAYsex
    #---save files---
    hexposhdr = np.array(hexposhdr.split(',')).astype(float)[0:5]
    #np.savetxt('hexapod_cray_position_'+expid+'.txt',[hexposhdr,hexHao,hexSex,posCRAY,posCRAYsex],fmt='%10.5f')
    hexBCAM = np.array([bcamDX,bcamDY,-999,bcamAX,bcamAY])
    np.savetxt('hexapod_position_'+expid+'.txt',[hexposhdr,hexHao,hexBCAM],fmt='%10.5f')
    
    return '----finished one image ----'
    

if __name__ == "__main__":
    from decamImgAnalyzer import *
    import sys,time,glob
    startTime=time.time()
    if len(sys.argv) == 1:
        print 'syntax: '
        print 'desImgAnalysis expid'
        print 'or'
        print 'desImgAnalysis all'
        print 'Note: The image need to be reduced (bias subtraction, flat fielding'
    elif sys.argv[1] == 'all':
        img_nameList = glob.glob('*_reduced.fits')
        nimg = len(img_nameList)
        for i in range(nimg):
            t=runanalysis(img_nameList[i])
    else:   
        expid = sys.argv[1]
        img_name = 'DECam_'+expid+'_reduced.fits'
        t=runanalysis(img_name)
    endTime=time.time()
    elapseTime=endTime-startTime
    print '---elapsed time: ' + str(elapseTime)

