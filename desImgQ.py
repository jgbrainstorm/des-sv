#! /usr/bin/env python
#------------------------------------------------
# This code is the replacement of desImgAnalysis.py
# J. Hao, 11/13/2012 @ FNAL 
#------------------------------------------------
import sys
import glob as gl
import cPickle as p
sys.path.append('/usr/remote/user/sispi/jiangang/des-sv')

from decamImgAnalyzer_def import *

def averageN30(data=None):
    """
    this function use the mean moments of N29,N31 to replace the moments of N30 because of the failure of N30. 
    """
    idxN29 = 59
    idxN30 = 60
    idxN31 = 61
    datanew = data.copy()
    datanew[idxN30,2:] = 0.5*(datanew[idxN29,2:]+datanew[idxN31,2:])
    return datanew

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
    bcamDX = pf.getheader(img_name,0)['BCAMDX']
    bcamDY = pf.getheader(img_name,0)['BCAMDY']
    bcamAX = pf.getheader(img_name,0)['BCAMAX']
    bcamAY = pf.getheader(img_name,0)['BCAMAY']
    data=[]
    stamplist=[]
    bkglist=[]
    dataSex=[]
    fwhmSex = np.array([])
    whiskerSex = np.array([])
    r50Sex = np.array([])
    nstarall = 0
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
        Mcc = cat.X2_IMAGE
        Mrr = cat.Y2_IMAGE
        Mrc = cat.XY_IMAGE
        fwhm_sex = cat.FWHM_IMAGE
        starFwhm = selectStar(mag,fwhm_sex)
        ok = (np.abs(fwhm_sex - starFwhm) < 0.4)*(x>100)*(x<2050)*(y>100)*(y<4100)*(flag == 0)*(mag<=-11.5)*(mag>-14.5)
        nstar = len(mag[ok])
        nstarall = nstarall + nstar
        print '--- Nstars selected: '+str(nstar)+'---'
        xccd = eval(imghdu[i].header['extname'])[1]
        yccd = eval(imghdu[i].header['extname'])[2]
        if ok.any():
            bkg = bkg[ok]
            Mrr = robust_mean(Mrr[ok])
            Mcc = robust_mean(Mcc[ok])
            Mrc = robust_mean(Mrc[ok])
            r50Sex = robust_mean(rad[ok])
            fwhmSex = robust_mean(fwhm_sex[ok])
            x=x[ok]
            y=y[ok]
            stamp = getStamp(data=img,xcoord=x,ycoord=y,Npix=25)
            moms = measureIQstamp(stamp,bkg,2.)
            data.append(moms)
            dataSex.append([Mcc,Mrr,Mrc,r50Sex,fwhmSex])
    data = np.array(data)
    dataSex = np.array(dataSex)
    datamean = np.array([robust_mean(data[:,0]),robust_mean(data[:,1]),robust_mean(data[:,2])])
    dataSexmean = np.array([robust_mean(dataSex[:,0]),robust_mean(dataSex[:,1]),robust_mean(dataSex[:,2]),robust_mean(dataSex[:,3]),robust_mean(dataSex[:,4])])
    datasubmean = data - datamean
    dataSexsubmean = dataSex - dataSexmean
    whk = ((datamean[0]-datamean[1])**2 + (2.*datamean[2])**2)**(0.25)*0.27
    phi = np.rad2deg(0.5*np.arctan2(2.*datamean[2],(datamean[0]-datamean[1])))
    whkrms = (robust_mean((datasubmean[:,0] - datasubmean[:,1])**2 + 4.*datasubmean[:,2]**2))**(0.25)*0.27
    r50=0.5*2.35482*np.sqrt((datamean[0]+datamean[1])/2.)*0.27
    whkSex = ((dataSexmean[0]-dataSexmean[1])**2 + (2.*dataSexmean[2])**2)**(0.25)*0.27
    phiSex = np.rad2deg(0.5*np.arctan2(2.*dataSexmean[2],(dataSexmean[0]-dataSexmean[1])))
    whkrmsSex = (robust_mean((dataSexsubmean[:,0] - dataSexsubmean[:,1])**2 + 4.*dataSexsubmean[:,2]**2))**(0.25)*0.27
    r50Sex = dataSexmean[3]*0.27
    fwhmSex = dataSexmean[4]*0.27
    p.dump([int(expid),whk,phi,whkrms,r50,whkSex,phiSex,whkrmsSex,r50Sex,fwhmSex],open('desIQ_measures_'+expid+'.p','w'))
    return '----finished one image ----'
    

if __name__ == "__main__":
    from desImgQ import *
    import sys,time,glob
    startTime=time.time()
    if len(sys.argv) == 1:
        print 'syntax: '
        print 'desImgQ.py expid'
        print 'or'
        print 'desImgQ.py all'
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

