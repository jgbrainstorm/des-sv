#! /usr/bin/env python
#----------------------------------------------
# this code compare the various IQ measurement
# J. Hao, 12/4/2012
#----------------------------------------------
import sys,os
sys.path.append('/usr/remote/user/sispi/jiangang/des-sv')
sys.path.append('/usr/remote/user/sispi/jiangang/decam-fermi')
sys.path.append('/usr/remote/user/sispi/jiangang/lib/python')
from decamImgAnalyzer_def import *
import cPickle as p

def measureIQ(img_name=None):
    catname = img_name[0:-5]+'_star_catalog_selected.fits'
    imghdu = pf.open(img_name)
    cathdu = pf.open(catname)
    expid = img_name[6:14]
    dimmfwhm = pf.getheader(img_name,0)['dimmsee']
    kernelSigma = np.sqrt(dimmfwhm**2+0.55**2)/2.35482
    bkglist=[]
    dataSex=[]
    fwhm = []
    whk = []
    r50 = []
    nstartotal = 0.
    for i in range(1,63):
        print i
        img = imghdu[i].data
        cat = cathdu[i].data
        x = cat.XWIN_IMAGE
        y = cat.YWIN_IMAGE
        rad = cat.FLUX_RADIUS
        bkg = cat.BACKGROUND
        Mcc = cat.X2WIN_IMAGE
        Mrr = cat.Y2WIN_IMAGE
        Mrc = cat.XYWIN_IMAGE
        fwhm_sex = cat.FWHM_IMAGE
        radall.append(rad)
        Nobj = len(Mrr)
        nstartotal=nstartotal+Nobj
        M20=np.zeros(Nobj)
        M22=np.zeros(Nobj).astype(complex)
        stamp = getStamp(data=img,xcoord=x,ycoord=y,Npix=25)
        for k in range(Nobj):
            M20[k] = Mrr[k] + Mcc[k]
            M22[k] = np.complex(Mcc[k] - Mrr[k],2*Mrc[k])
            whisker_sex = np.sqrt(np.abs(M22[k]))
            res = get_fwhm_whisker(stampImg=stamp[k],bkg =bkg[k],sigma=2.)
            fwhm.append(np.append(res[1],fwhm_sex[k]*0.27))
            whk.append(np.append(res[0],whisker_sex*0.27)) 
            r50.append(np.append(res[2],rad[k]*0.27))
    r50 = np.array(r50)
    fwhm = np.array(fwhm)
    whk = np.array(whk)
    #fwhm: wmoment,amoment,moffat,gauss1d,sech2,sex
    #whk: wmoments,amoments,sex
    #r50: sech2, moffat,gauss1d,sex
    if nstartotal > 100:
        p.dump([fwhm,whk,r50],open('compare_fwhm_whisker_data_'+expid+'.p','w')) # save the fwhm and whisker data.
        pl.figure(figsize=(12,8))
        pl.subplot(2,1,1)
        pl.boxplot(list(fwhm.T))
        pl.ylabel('fwhm')
        pl.ylim(0.5,3.5)
        pl.grid()
        pl.xticks(np.arange(1,7),['wmoments','amoments','moffat','gauss','sech2','sextractor'])
        pl.title('exposure id: '+expid)
        pl.subplot(2,1,2)
        pl.boxplot(list(r50.T))
        pl.ylabel('R50')
        pl.xticks(np.arange(1,5),['sech2','moffat','gauss1d','sextractor'])
        pl.grid()
        pl.ylim(0,1.5)
        pl.savefig('fwhm_r50_'+expid+'.png')
        pl.close()
    return '---done!---'

if __name__ == "__main__":
    from compare_IQ import *
    import sys,time,glob
    startTime=time.time()
    if len(sys.argv) == 1:
        print 'syntax: '
        print 'compare_IQ expid'
        print 'or'
        print 'compare_IQ all'
        print 'Note: The image need to be reduced (bias subtraction, flat fielding'
    elif sys.argv[1] == 'all':
        img_nameList = glob.glob('*_reduced.fits')
        nimg = len(img_nameList)
        for i in range(nimg):
            t=measureIQ(img_nameList[i])
    else:   
        expid = sys.argv[1]
        img_name = 'DECam_'+expid+'_reduced.fits'
        t=measureIQ(img_name)
    endTime=time.time()
    elapseTime=endTime-startTime
    print '---elapsed time: ' + str(elapseTime)


