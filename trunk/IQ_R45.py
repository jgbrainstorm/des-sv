#---------------------------------
# This is the code for IQ-R4,5.
# Initiated by J. Hao, 7/13/2012
# Comments sent to: jghao@fnal.gov 
#---------------------------------
import numpy as np
import pyfits as pf
import scipy.ndimage as nd
import pylab as pl
import sys
import binplot as bp
from decamRay import *

scale = 0.27

def sech(x,width,height):
    """
    hyperbolic secant function
    """
    z = x/width
    res = height*2./(np.exp(z)+np.exp(-z))
    return res


def findbstr(data=None, hdr=None):
    """
    find the bright stars on the image
    """
    saturate = hdr['saturate']
    bsIDX = (data >= 0.3*saturate)* (data <= 0.8*saturate)
    good=nd.binary_opening(bsIDX,structure=np.ones((3,3)))  
    objData = data*good
    seg,nseg=nd.label(good,structure=np.ones((3,3)))  
    coords=nd.center_of_mass(objData,seg,range(1,nseg+1))
    xcoord=np.array([x[1] for x in coords])
    ycoord=np.array([x[0] for x in coords])
    return xcoord, ycoord

def getStamp(data=None,xcoord=None,ycoord=None,Npix = None):
    """
    Input: CCD image in maxtrix, x, y centroid of stars,the stamp npix
    Output: a list of stamp image around the x,y centroid
    """
    Nstar = len(xcoord)
    rowcen = ycoord
    colcen = xcoord
    stampImg=[]
    for i in range(Nstar):
        Img = data[int(rowcen[i]-Npix/2):int(rowcen[i]+Npix/2),int(colcen[i]-Npix/2):int(colcen[i]+Npix/2)]
        stampImg.append(Img)
    return stampImg

def wr(x=None,y=None,xcen=None,ycen=None,sigma=None):
    """
    Returns a gaussian weight function with the given parameters
    """
    res=np.exp(-((x-xcen)**2+(y-ycen)**2)/(2.*sigma**2))/(2.*np.pi*sigma**2) 
    return res

def adaptiveCentroid(data=None,sigma=None):
    """
    calculate the centroid using the adaptive approach
    """
    nrow,ncol=data.shape
    Isum = data.sum()
    Icol = data.sum(axis=0) # sum over all rows
    Irow = data.sum(axis=1) # sum over all cols
    IcolSum = np.sum(Icol)
    IrowSum = np.sum(Irow)
    colgrid = np.arange(ncol)
    rowgrid = np.arange(nrow)
    ROW,COL=np.indices((nrow,ncol))
    rowmean=np.sum(rowgrid*Irow)/IrowSum
    colmean=np.sum(colgrid*Icol)/IcolSum
    maxItr = 50
    EP = 0.0001
    for i in range(maxItr):
        wrmat = wr(ROW,COL,rowmean,colmean,sigma)
        IWmat = data*wrmat
        wrcol = wrmat.sum(axis=0)
        wrrow = wrmat.sum(axis=1)
        wrcolsum = np.sum(Icol*wrcol)
        wrrowsum = np.sum(Irow*wrrow)
        drowmean = np.sum((rowgrid-rowmean)*Irow*wrrow)/wrrowsum
        dcolmean = np.sum((colgrid-colmean)*Icol*wrcol)/wrcolsum
        rowmean = rowmean+2.*drowmean
        colmean = colmean+2.*dcolmean
        print drowmean**2+dcolmean**2
        if drowmean**2+dcolmean**2 <= EP:
            break
    return rowmean,colmean


def complexMoments(image=None,sigma=None):
    """
    This one calcualte the 3 2nd moments and 4 thrid moments with the Gaussian weights.
    col : x direction
    row : y direction
    the centroid is using the adpative centroid.
    sigma is the stand deviation of the measurement kernel in pixel
    The output is in pixel coordinate
    """
    nrow,ncol=image.shape
    Isum = image.sum()
    Icol = image.sum(axis=0) # sum over all rows
    Irow = image.sum(axis=1) # sum over all cols
    IcolSum = np.sum(Icol)
    IrowSum = np.sum(Irow)
    colgrid = np.arange(ncol)
    rowgrid = np.arange(nrow)
    ROW,COL=np.indices((nrow,ncol))
    rowmean=np.sum(rowgrid*Irow)/IrowSum
    colmean=np.sum(colgrid*Icol)/IcolSum
    maxItr = 50
    EP = 0.0001 # start getting the adaptive centroid
    for i in range(maxItr):  
        wrmat = wr(ROW,COL,rowmean,colmean,sigma)
        IWmat = image*wrmat
        wrcol = wrmat.sum(axis=0)
        wrrow = wrmat.sum(axis=1)
        wrcolsum = np.sum(Icol*wrcol)
        wrrowsum = np.sum(Irow*wrrow)
        drowmean = np.sum((rowgrid-rowmean)*Irow*wrrow)/wrrowsum
        dcolmean = np.sum((colgrid-colmean)*Icol*wrcol)/wrcolsum
        rowmean = rowmean+2.*drowmean
        colmean = colmean+2.*dcolmean
        if drowmean**2+dcolmean**2 <= EP:
            break
    rowgrid = rowgrid - rowmean # centered
    colgrid = colgrid - colmean
    Mr = np.sum(rowgrid*Irow*wrrow)/wrrowsum
    Mc = np.sum(colgrid*Icol*wrcol)/wrcolsum
    Mrr = np.sum(rowgrid**2*Irow*wrrow)/(wrrowsum)
    Mcc = np.sum(colgrid**2*Icol*wrcol)/(wrcolsum)
    Mrc = np.sum(np.outer(rowgrid,colgrid)*IWmat)/np.sum(IWmat)
    Mrrr = np.sum(rowgrid**3*Irow*wrrow)/(wrrowsum)
    Mccc = np.sum(colgrid**3*Icol*wrcol)/(wrcolsum)
    Mrrc = np.sum(np.outer(rowgrid**2,colgrid)*IWmat)/np.sum(IWmat)
    Mrcc = np.sum(np.outer(rowgrid,colgrid**2)*IWmat)/np.sum(IWmat)
    M20 = Mrr + Mcc
    M22 = complex(Mcc - Mrr,2*Mrc)
    M31 = complex(3*Mc - (Mccc+Mrrc)/sigma**2, 3*Mr - (Mrcc + Mrrr)/sigma**2)
    M33 = complex(Mccc-3*Mrrc, 3.*Mrcc - Mrrr)
    return M20, M22, M31, M33

def fwhm_whisker(image=None, sigma = 1.1/scale):
    """
    This code calculate the fwhm and wisker length defined as (M22.real^2 + M22.imag^2)^{1/4}
    input: 
         data: 2d stamp image
         sigma: std of the Gaussian weight Kernel in pixel
    """
    M20, M22, M31, M33 = complexMoments(image=image,sigma=sigma)
    #fwhm = np.sqrt(M20/2.)*2.35482*scale
    whisker_length = np.sqrt(np.abs(M22))*scale
    lambdap = 0.5*(M20 + abs(M22))
    lambdam = 0.5*(M20 - abs(M22))
    fwhm = np.sqrt(2.*np.log(2.))*(np.sqrt(lambdap)+np.sqrt(lambdam))*scale
    #fwhm = (1./(M20/2.) - 1./sigma**2)**(-0.5)*2.35482*scale
    return fwhm,whisker_length


def dispStamp(stampImg=None,sigma=1.1/scale):
    npix = stampImg.shape[0]
    pl.figure(figsize=(18,6))
    pl.subplot(1,3,1)
    pl.matshow(stampImg,fignum=False)
    mfit = mfwhm(stampImg)
    gfit = gfwhm(stampImg)
    s2fit = s2fwhm(stampImg)
    g2dfit = g2dfwhm(stampImg)
    pl.xlabel('Pixel')
    pl.ylabel('Pixel')
    pl.grid(color='y')
    rowCen,colCen = adaptiveCentroid(data=stampImg,sigma=sigma)
    M20, M22, M31, M33 =complexMoments(stampImg,sigma=sigma)
    e1 = M22.real/M20.real
    e2 = M22.imag/M20.real
    whiskerLength = np.sqrt(np.abs(M22))*scale
    lambdap = 0.5*(M20 + abs(M22))
    lambdam = 0.5*(M20 - abs(M22))
    fwhm = (1./(M20/2.) - 1./sigma**2)**(-0.5)*2.35482*scale
    pl.figtext(0.15,0.8, 'e1: '+str(round(e1,3)) + ',  e2: '+str(round(e2,3)), color='r')
    pl.figtext(0.15,0.75, 'rowCen: '+str(round(rowCen,4)) + ',  colCen: '+str(round(colCen,4)), color='r')
    pl.figtext(0.15,0.7, 'PSF whisker_Wmoments: '+str(round(whiskerLength,4))+' [arcsec]', color='r')
    pl.figtext(0.15,0.65, 'PSF whisker_Amoments: '+str(round(g2dfit[2]*scale,4))+' [arcsec]', color='r')
    pl.subplot(1,3,2)
    row,col = np.mgrid[0:npix,0:npix]
    row = row - rowCen
    col = col - colCen
    radius = np.sqrt(row**2+col**2)
    img = stampImg.flatten()
    radius = radius.flatten()
    idx = np.argsort(radius)
    img = img[idx]
    radius = radius[idx]
    halfmax = np.median(img[0:10])/2.
    pl.plot(radius,img,'k.')
    pl.grid(color='y')
    pl.hlines(halfmax,0,radius.max(),linestyle='solid',colors='b')
    pl.hlines(mfit[2]/2.,0,radius.max(),linestyle='solid',colors='r')
    pl.hlines(gfit[1]/2.,0,radius.max(),linestyle='solid',colors='g')
    pl.hlines(s2fit[1]/2.,0,radius.max(),linestyle='solid',colors='m')
    pl.hlines(g2dfit[0]/2.,0,radius.max(),linestyle='solid',colors='c',label='Adaptive Moments')
    pl.vlines(fwhm/scale/2.,0, halfmax*4,linestyle='solid',colors='b',label='Weighted Moments')
    pl.vlines(mfit[4]/2.,0, halfmax*4,linestyle='solid',colors='r')
    pl.vlines(gfit[3]/2.,0, halfmax*4,linestyle='solid',colors='g')
    pl.vlines(s2fit[3]/2.,0, halfmax*4,linestyle='solid',colors='m')
    pl.vlines(g2dfit[3]/2.,0, halfmax*4,linestyle='solid',colors='c')
    pl.plot(radius,mprofile(radius,mfit[0],mfit[1],mfit[2],mfit[3]),'r-',label='Moffat Fit')
    pl.plot(radius,gprofile(radius,gfit[0],gfit[1],gfit[2]),'g-',label='Gaussian Fit')
    pl.plot(radius,s2profile(radius,s2fit[0],s2fit[1],s2fit[2]),'m-',label='Sech2 Fit')
    pl.legend(loc='best')
    pl.ylim(0,halfmax*4)
    pl.xlim(0,npix/4.) 
    pl.xlabel('Radius [pixels]')
    pl.ylabel('Mean counts [ADU]')
    pl.title('Radial profile')
    pl.figtext(0.65,0.7,'Gaussian Weight '+r'$\sigma$: '+str(round(sigma*scale,3))+ ' arcsec',color='r')
    pl.figtext(0.65,0.6,'FWHM_Gaussian: '+str(round(gfit[3]*scale,3))+ ' arcsec')
    pl.figtext(0.65,0.55,'FWHM_Moffat: '+str(round(mfit[4]*scale,3))+ ' arcsec')
    pl.figtext(0.65,0.5,'FWHM_Sech2: '+str(round(s2fit[3]*scale,3))+ ' arcsec')
    pl.figtext(0.65,0.45,'FWHM_Wmoments: '+str(round(fwhm,3))+ ' arcsec') 
    pl.figtext(0.65,0.4,'FWHM_Amoments: '+str(round(g2dfit[3]*scale,3))+ ' arcsec')
    return '---- Done! ----'
    
def dispStampList(stampImgList=None,sigma=None):
    if sigma == None:
        print 'syntax: dispStampList(stampImgList,sigma)'
        sys.exit()
    Nstamp = len(stampImgList)
    for i in range(Nstamp):
        t=dispStamp(stampImg=stampImgList[i],sigma=sigma)
        raw_input('--- hit the enter key to proceed ---')
        pl.close()
    return ' ---- Done ! ----'
        

def rowcol2XY(row,col,CCD):
    """
    convert the row/col [in pixels] of a given CCD to the x, y
    of the Focal plane [in mm] by assuming a constant pixel scale 0.015mm/pix
    Input: row, col coordinate of the object, the CCD position i.e. S1, S2, etc
    Output: the x, y coordinate in the FP coordiante in mm.
    Convention: 
    1. each ccd, the origin of row and col is the south east corner.
    2. the direction row increase is West
    3. the direction col increase is North.
    4. In my Focal Plane definition file: DECam_def.py,
        positive X is South
        positive Y is East
    So, the row increase as -Y direction.
        the col increase as -X direction.
    """
    pixscale = 0.015 #mm/pix
    X = CCD[1]+1024*pixscale-(col*pixscale+pixscale/2.)
    Y = CCD[2]+2048*pixscale-(row*pixscale+pixscale/2.)
    return X,Y
    

if __name__ == "__main__":
    from IQ_R45 import *
    pl.ion()
    #dir = '/home/jghao/research/data/des_optics_psf/dc6b_image/'
    dr = '/home/jghao/research/data/des_optics_psf/dc6b_image/goodseeing/decam--28--49-r-1/'
    starfile='/home/jghao/research/data/des_optics_psf/dc6b_image/goodseeing/catfile/decam_-27.72186_-48.60000-objects.fit'
    #imgname = 'decam-34-0-r-0_30.fits'
    #bkgname = 'decam-34-0-r-0_30_bkg.fits'
    #catname = 'decam-34-0-r-0_scamp.fits'
    ext='01'
    imgname= 'decam--28--49-r-1_'+ext+'.fits.fz'
    bkgname = 'decam--28--49-r-1_'+ext+'_bkg.fits.fz'
    hdr = pf.getheader(dr+imgname,1)
    data = pf.getdata(dr+imgname) - pf.getdata(dr+bkgname)
    #xc,yc=findbstr(data=data, hdr=hdr)
    xc = pf.getdata(starfile,int(ext)).xccd
    yc = pf.getdata(starfile,int(ext)).yccd
    rmag = pf.getdata(starfile,int(ext)).mag_3
    ok = (rmag > 12)*(rmag < 16)
    xc=xc[ok]
    yc = yc[ok]
    stamp=getStamp(data=data,xcoord=xc,ycoord=yc,Npix =25)
    Nstamp = len(stamp)
    for i in range(Nstamp):
        print fwhm_whisker(image=stamp[i], sigma = 4.)
        
