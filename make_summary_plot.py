#! /usr/bin/env python
# ---- make the plot for Gary ----

import numpy as np, pylab as pl, cPickle as p,glob as gl
from DECamCCD import *

def getMeanStd(filelist,ext):
    nfile = len(filelist)
    resMean = []
    resStd=[]
    for i in range(nfile):
        t = p.load(open(filelist[i],'r'))    
        resMean.append(robust_mean(t[ext][~np.isnan(t[ext])]))
        resStd.append(t[ext][~np.isnan(t[ext])].std())
    return np.array(resMean)

def getMedian(filelist,ext):
    nfile = len(filelist)
    resMedian = []
    resStd=[]
    for i in range(nfile):
        t = p.load(open(filelist[i],'r'))    
        resMedian.append(np.median(t[ext][~np.isnan(t[ext])]))
    return np.array(resMedian)



f1 = gl.glob('/home/jghao/research/desSV/12022012/fwhm*.p')
f1.sort()

f2 = gl.glob('/home/jghao/research/desSV/11232012/fwhm*.p')
f2.sort()
       
f3 = gl.glob('/home/jghao/research/desSV/11242012/rband/fwhm*.p') +gl.glob('/home/jghao/research/desSV/11242012/iband/fwhm*.p')+ gl.glob('/home/jghao/research/desSV/11242012/zband/fwhm*.p')

# --- whisker ----
whk1202w = getMeanStd(f1,6)
whk1123w = getMeanStd(f2,6)
whk1124w = getMeanStd(f3,6)

whk1202a = getMeanStd(f1,7)
whk1123a = getMeanStd(f2,7)
whk1124a = getMeanStd(f3,7)

whk1202s = getMeanStd(f1,8)
whk1123s = getMeanStd(f2,8)
whk1124s = getMeanStd(f3,8)

#----fwhm-- weighted, moffat,sextractor ---
fwhm1202w = getMedian(f1,0)
fwhm1123w = getMedian(f2,0)
fwhm1124w = getMedian(f3,0)

fwhm1202m = getMedian(f1,2)
fwhm1123m = getMedian(f2,2)
fwhm1124m = getMedian(f3,2)


fwhm1202s = getMedian(f1,5)
fwhm1123s = getMedian(f2,5)
fwhm1124s = getMedian(f3,5)


#---R50, moffat, gauss, sech2 ---
r501202g = getMedian(f1,3)/2.
r501123g = getMedian(f2,3)/2.
r501124g = getMedian(f3,3)/2.


r501202sech2 =getMedian(f1,4)*0.592938
r501123sech2 = getMedian(f2,4)*0.592938
r501124sech2 = getMedian(f3,4)*0.592938





pl.figure(figsize=(15,15))
pl.subplot(3,3,1)
pl.hist(whk1123w,bins=20,alpha=0.5)
pl.xlabel('Whiker Length [weighed moments]')
pl.vlines(0.2,0,10,color='r',lw=2)
pl.ylim(0,10)
pl.title('11/23/2012')
pl.grid()

pl.subplot(3,3,2)
pl.hist(whk1124w,bins=20,alpha=0.5)
pl.xlabel('Whiker Length [weighed moments]')
pl.vlines(0.2,0,10,color='r',lw=2)
pl.ylim(0,10)
pl.title('11/24/2012')
pl.grid()

pl.subplot(3,3,3)
pl.hist(whk1202w,bins=20,alpha=0.5)
pl.xlabel('Whiker Length [weighed moments]')
pl.vlines(0.2,0,10,color='r',lw=2)
pl.ylim(0,10)
pl.title('12/02/2012')
pl.grid()

pl.subplot(3,3,4)
pl.hist(whk1123s,bins=20,alpha=0.5)
pl.xlabel('Whiker Length [sextractor]')
pl.vlines(0.2,0,10,color='r',lw=2)
pl.title('11/23/2012')
pl.grid()

pl.subplot(3,3,5)
pl.hist(whk1124s,bins=20,alpha=0.5)
pl.xlabel('Whiker Length [sextractor]')
pl.vlines(0.2,0,10,color='r',lw=2)
pl.title('11/24/2012')
pl.grid()

pl.subplot(3,3,6)
pl.hist(whk1202s,bins=20,alpha=0.5)
pl.xlabel('Whiker Length [sextractor]')
pl.vlines(0.2,0,10,color='r',lw=2)
pl.ylim(0,10)
pl.title('12/02/2012')
pl.grid()

pl.subplot(3,3,7)
pl.hist(whk1123a,bins=20,alpha=0.5)
pl.xlabel('Whiker Length [adaptive moments]')
pl.vlines(0.2,0,10,color='r',lw=2)
pl.title('11/23/2012')
pl.subplot(3,3,8)
pl.grid()

pl.hist(whk1124a,bins=20,alpha=0.5)
pl.xlabel('Whiker Length [adaptive moments]')
pl.vlines(0.2,0,10,color='r',lw=2)
pl.title('11/24/2012')
pl.grid()

pl.subplot(3,3,9)
pl.hist(whk1202a,bins=20,alpha=0.5)
pl.xlabel('Whiker Length [adaptive moments]')
pl.vlines(0.2,0,10,color='r',lw=2)
pl.ylim(0,10)
pl.title('12/02/2012')
pl.grid()
pl.figtext(0.35,0.95,'Distribution of Mean Whisker Length',color='b',fontsize=18)

pl.savefig('whisker_robust_mean_summary.png')


#----median fwhm distribution  ----
pl.figure(figsize=(15,15))
pl.subplot(3,3,1)
pl.hist(fwhm1123w,bins=20,alpha=0.5)
pl.xlabel('Whiker Length [weighed moments]')
pl.vlines(0.9,0,10,color='r',lw=2)
pl.ylim(0,10)
pl.xlim(0.6,2.0)
pl.title('11/23/2012')
pl.grid()

pl.subplot(3,3,2)
pl.hist(fwhm1124w,bins=20,alpha=0.5)
pl.xlabel('Whiker Length [weighed moments]')
pl.vlines(0.9,0,10,color='r',lw=2)
pl.ylim(0,10)
pl.xlim(0.6,2.0)
pl.title('11/24/2012')
pl.grid()

pl.subplot(3,3,3)
pl.hist(fwhm1202w,bins=20,alpha=0.5)
pl.xlabel('Whiker Length [weighed moments]')
pl.vlines(0.9,0,10,color='r',lw=2)
pl.ylim(0,10)
pl.xlim(0.6,2.0)
pl.grid()
pl.title('12/02/2012')

pl.subplot(3,3,4)
pl.hist(fwhm1123s,bins=20,alpha=0.5)
pl.xlabel('Whiker Length [sextractor]')
pl.vlines(0.9,0,10,color='r',lw=2)
pl.title('11/23/2012')
pl.ylim(0,10)
pl.xlim(0.6,2.0)
pl.grid()

pl.subplot(3,3,5)
pl.hist(fwhm1124s,bins=20,alpha=0.5)
pl.xlabel('Whiker Length [sextractor]')
pl.vlines(0.9,0,10,color='r',lw=2)
pl.title('11/24/2012')
pl.ylim(0,10)
pl.xlim(0.6,2.0)
pl.grid()

pl.subplot(3,3,6)
pl.hist(fwhm1202s,bins=20,alpha=0.5)
pl.xlabel('Whiker Length [sextractor]')
pl.vlines(0.9,0,10,color='r',lw=2)
pl.ylim(0,10)
pl.xlim(0.6,2.0)
pl.grid()
pl.title('12/02/2012')
pl.savefig('whisker_summary.png')

pl.subplot(3,3,7)
pl.hist(fwhm1123m,bins=20,alpha=0.5)
pl.xlabel('Whiker Length [moffat]')
pl.vlines(0.9,0,10,color='r',lw=2)
pl.title('11/23/2012')
pl.ylim(0,10)
pl.xlim(0.6,2.0)
pl.grid()

pl.subplot(3,3,8)
pl.hist(fwhm1124m,bins=20,alpha=0.5)
pl.xlabel('Whiker Length [moffat]')
pl.vlines(0.9,0,10,color='r',lw=2)
pl.title('11/24/2012')
pl.ylim(0,10)
pl.xlim(0.6,2.0)
pl.grid()

pl.subplot(3,3,9)
pl.hist(fwhm1202m,bins=20,alpha=0.5)
pl.xlabel('Whiker Length [moffat]')
pl.vlines(0.9,0,10,color='r',lw=2)
pl.ylim(0,10)
pl.xlim(0.6,2.0)
pl.grid()
pl.title('12/02/2012')
pl.figtext(0.35,0.95,'Distribution of Median FWHM',color='b',fontsize=18)
pl.savefig('fwhm_median_summary.png')



#----median R50 distribution  ----
pl.figure(figsize=(15,10))
pl.subplot(2,3,1)
pl.hist(r501123g,bins=20,alpha=0.5)
pl.xlabel('R50 [Gaussian Profle]')
pl.vlines(0.522,0,10,color='r',lw=2)
pl.ylim(0,10)
pl.xlim(0.2,1.0)
pl.title('11/23/2012')
pl.grid()

pl.subplot(2,3,2)
pl.hist(r501124g,bins=20,alpha=0.5)
pl.xlabel('R50 [Gaussian Profile]')
pl.vlines(0.522,0,10,color='r',lw=2)
pl.ylim(0,10)
pl.xlim(0.2,1.0)
pl.title('11/24/2012')
pl.grid()

pl.subplot(2,3,3)
pl.hist(r501202g,bins=20,alpha=0.5)
pl.xlabel('R50 [Gaussian Profile]')
pl.vlines(0.522,0,10,color='r',lw=2)
pl.ylim(0,10)
pl.xlim(0.2,1.0)
pl.grid()
pl.title('11.02/2012')

pl.subplot(2,3,4)
pl.hist(r501123sech2,bins=20,alpha=0.5)
pl.xlabel('R50 [Sech2 profile]')
pl.vlines(0.522,0,10,color='r',lw=2)
pl.title('11/23/2012')
pl.ylim(0,10)
pl.xlim(0.2,1.0)
pl.grid()

pl.subplot(2,3,5)
pl.hist(r501124sech2,bins=20,alpha=0.5)
pl.xlabel('R50 [Sech2 profile]')
pl.vlines(0.522,0,10,color='r',lw=2)
pl.title('11/24/2012')
pl.ylim(0,10)
pl.xlim(0.2,1.0)
pl.grid()

pl.subplot(2,3,6)
pl.hist(r501202sech2,bins=20,alpha=0.5)
pl.xlabel('R50 [Sech2 Profile]')
pl.vlines(0.522,0,10,color='r',lw=2)
pl.ylim(0,10)
pl.xlim(0.2,1.0)
pl.grid()
pl.title('12/02/2012')
pl.figtext(0.35,0.95,'Distribution of Median R50',color='b',fontsize=18)
pl.savefig('r50_median_summary.png')


#---first cut plot ----

b=np.genfromtxt('firstcut_before_163831.txt',delimiter=',')
#b=np.genfromtxt('firstcut_163831.txt',delimiter=',')
ok = b[:,2] != -999.
b = b[ok,:]
pl.figure(figsize=(10,5))
pl.subplot(1,2,1)
pl.plot(b[:,2],b[:,6],'b.')
pl.plot([0,1],[0,1],'r-')
pl.xlim(0,0.5)
pl.ylim(0,0.5)
pl.xlabel('whisker_hao')
pl.ylabel('whisker_mike')
pl.grid()

pl.subplot(1,2,2)
pl.plot(b[:,3],b[:,7],'b.')
pl.plot([0,1],[0,1],'r-')
pl.xlim(0,0.5)
pl.ylim(0,0.5)
pl.xlabel('whisker_rms_hao')
pl.ylabel('whisker_rms_mike')
pl.grid()

pl.savefig('firstcut_before_163831_whisker_hao_mike.png')
