#! /usr/bin/env python
#-------------------------------------------------------------
# This set of codes test the DECam image quality, IQ_R4, IQ-R5
# It measure the FWHM and the whisker using various ways based
# on imput image and a list of star positions
# Created by: Jiangang Hao @ Fermilab, 8/1/2012
#-------------------------------------------------------------

try:
    import numpy as np
    import pyfits as pf
    import scipy.ndimage as nd
    import pylab as pl
    import sys
    from scipy.optimize import leastsq
except ImportError:
    print "Error: missing one of the libraries (numpy, pyfits, scipy, matplotlib)"
    sys.exit()


scale=0.27


def findbstr(data=None, hdr=None):
    """
    find the bright stars on the image
    """
    saturate = hdr['saturate']
    bsIDX = (data >= 0.3*saturate)* (data <= 0.5*saturate)
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


def moments(data):
    """
    Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution by calculating its
    moments
    """
    total = data.sum()
    if total != 0.:
        X, Y = np.indices(data.shape)
        x = (X*data).sum()/total
        y = (Y*data).sum()/total
        if int(y) < data.shape[0] and int(x)< data.shape[0]:
            col = data[:, int(y)]
            row = data[int(x), :]
            if col.sum() != 0. and row.sum() != 0.:      
                width_x = np.sqrt(abs((np.arange(col.size)-y)**2*col).sum()/col.sum())
                width_y = np.sqrt(abs((np.arange(row.size)-x)**2*row).sum()/row.sum())
                height = data.max()
            else:
                height=0
                x=0
                y=0
                width_x=0
                width_y=0
        else:
            height=0
            x=0
            y=0
            width_x=0
            width_y=0
    else:
        height=0
        x=0
        y=0
        width_x=0
        width_y=0
    return height,np.sqrt(width_x**2 + width_y**2)


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
    colgrid = np.arange(ncol)
    rowgrid = np.arange(nrow)
    rowmean=np.sum(rowgrid*Irow)/Isum
    colmean=np.sum(colgrid*Icol)/Isum
    ROW,COL=np.indices((nrow,ncol))
    maxItr = 50
    EP = 0.0001
    for i in range(maxItr):
        wrmat = wr(ROW,COL,rowmean,colmean,sigma)
        IWmat = data*wrmat
        IWcol = IWmat.sum(axis=0)
        IWrow = IWmat.sum(axis=1)
        drowmean = np.sum((rowgrid-rowmean)*IWrow)/np.sum(IWrow)
        dcolmean = np.sum((colgrid-colmean)*IWcol)/np.sum(IWcol)
        rowmean = rowmean+2.*drowmean
        colmean = colmean+2.*dcolmean
        if drowmean**2+dcolmean**2 <= EP:
            break

    return rowmean,colmean

 
def complexMoments(data=None,sigma=None):
    """
    This one calcualte the 3 2nd moments and 4 thrid moments with the Gaussian weights.
    col : x direction
    row : y direction
    the centroid is using the adpative centroid.
    sigma is the stand deviation of the measurement kernel in pixel
    The output is in pixel coordinate
    """
    nrow,ncol=data.shape
    Isum = data.sum()
    Icol = data.sum(axis=0) # sum over all rows
    Irow = data.sum(axis=1) # sum over all cols
    colgrid = np.arange(ncol)
    rowgrid = np.arange(nrow)
    rowmean=np.sum(rowgrid*Irow)/Isum
    colmean=np.sum(colgrid*Icol)/Isum
    ROW,COL=np.indices((nrow,ncol))
    maxItr = 50
    EP = 0.0001
    for i in range(maxItr):
        wrmat = wr(ROW,COL,rowmean,colmean,sigma)
        IWmat = data*wrmat
        IWcol = IWmat.sum(axis=0)
        IWrow = IWmat.sum(axis=1)
        IWsum = IWmat.sum()
        drowmean = np.sum((rowgrid-rowmean)*IWrow)/IWsum
        dcolmean = np.sum((colgrid-colmean)*IWcol)/IWsum
        rowmean = rowmean+2.*drowmean
        colmean = colmean+2.*dcolmean
        if drowmean**2+dcolmean**2 <= EP:
            break
    rowgrid = rowgrid - rowmean # centered
    colgrid = colgrid - colmean
    Mr = np.sum(rowgrid*IWrow)/IWsum
    Mc = np.sum(colgrid*IWcol)/IWsum
    Mrr = np.sum(rowgrid**2*IWrow)/IWsum
    Mcc = np.sum(colgrid**2*IWcol)/IWsum
    Mrc = np.sum(np.outer(rowgrid,colgrid)*IWmat)/IWsum
    Mrrr = np.sum(rowgrid**3*IWrow)/IWsum
    Mccc = np.sum(colgrid**3*IWcol)/IWsum
    Mrrc = np.sum(np.outer(rowgrid**2,colgrid)*IWmat)/IWsum
    Mrcc = np.sum(np.outer(rowgrid,colgrid**2)*IWmat)/IWsum
    print Mrrr, Mccc, Mrrc, Mrcc
    M20 = Mrr + Mcc
    M22 = complex(Mcc - Mrr,2*Mrc)
    M31 = complex(3*Mc - (Mccc+Mrrc)/sigma**2, 3*Mr - (Mrcc + Mrrr)/sigma**2)
    M33 = complex(Mccc-3*Mrrc, 3.*Mrcc - Mrrr)
    return M20, M22, M31, M33


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

def s2profile(r,r0,A,B):
    """
    hyperbolic secant square function
    """
    x = r/r0
    res = A*4./(np.exp(x)+np.exp(-x))**2 + B
    return res

def gprofile(r,sig,A,B):
    """
    Fit the binned distribution to a 1D gaussian profile with a constant
    """
    res = A*np.exp(-0.5*(r/sig)**2)+B
    return res

def mprofile(r, alpha, beta,A,B):
    """
    Fit the light distribution to a Moffat profile
    """
    res = A*(1+(r/alpha)**2)**(-beta)+B
    return res

def gaussian2d(x,y,xc,yc,sigmax,sigmay,rho,A,B):
    """
    2D Gaussian profile with a constant
    """
    res = A*np.exp(-0.5/(1-rho**2)*(x**2/sigmax**2+y**2/sigmay**2-2.*rho*x*y/(sigmax*sigmay)))+B
    return res


def gfwhm(img):
    """
    measure the fwhm based on a 1D Gaussian fit
    """
    npix = img.shape[0]
    rowCen,colCen = adaptiveCentroid(img,1.1/scale)
    row,col = np.mgrid[0:npix,0:npix]
    row = row - rowCen
    col = col - colCen
    A0,sig0 = moments(img)
    radius = np.sqrt(row**2+col**2)
    img = img.flatten()
    ok = img >0
    img = img[ok]
    radius = radius.flatten()
    radius = radius[ok]
    def residualg(p,r,I):
        sig,A,B = p
        Ierr = np.sqrt(abs(I))
        res = (gprofile(radius,sig,A,B) - I)/Ierr
        return res 
    B0 = 0.
    p0=np.array([sig0,A0,B0])
    p = leastsq(residualg,p0,args=(radius,img))[0]
    sig,A,B = p
    fwhm_gauss= 2. * sig * np.sqrt(2. * np.log(2.))
    return sig,A,B,fwhm_gauss

def s2fwhm(img):
    """
    measure the fwhm by fitting a sech2 profile
    """
    npix = img.shape[0]
    rowCen,colCen = adaptiveCentroid(img,1.1/scale)
    row,col = np.mgrid[0:npix,0:npix]
    row = row - rowCen
    col = col - colCen
    A0,r0_0 = moments(img)
    radius = np.sqrt(row**2+col**2)
    img = img.flatten()
    ok = img >0
    img = img[ok]
    radius = radius.flatten()
    radius = radius[ok]
    def residuals2(p,r,I):
        r0,A,B = p
        Ierr = np.sqrt(abs(I))
        res = (s2profile(radius,r0,A,B) - I)/Ierr
        return res 
    B0 = 0.
    p0=np.array([r0_0,A0,B0])
    p = leastsq(residuals2,p0,args=(radius,img))[0]
    r0,A,B = p
    fwhm_sech2= 1.7627471*r0 # obtained by solving the equation
    return r0,A,B,fwhm_sech2


def g2dfwhm(img):
    """
    x is col, y is row
    sigmax -> sigmac, sigmay -> sigmar
    """
    npix = img.shape[0]
    rowCen,colCen = adaptiveCentroid(img,1.1/scale)
    row,col = np.mgrid[0:npix,0:npix]
    row = row - rowCen
    col = col - colCen
    A0,sigmac0 = moments(img)
    sigmar0 = sigmac0
    rho0 = 0.
    B0 = 0.
    p0=np.array([sigmac0,sigmar0,rho0,A0, B0])
    def residualg2d(p,x,y,xc,yc,I):
        sigmax,sigmay,rho,A,B = p
        Ierr = np.sqrt(abs(I))+0.00001 # to avoid those = 0, add a small number 
        res = (gaussian2d(x,y,xc,yc,sigmax,sigmay,rho,A,B) - I)/Ierr
        return res.flatten()
    p = leastsq(residualg2d,p0,args=(col,row,colCen,rowCen,img))[0]
    sigmac,sigmar,rho,A,B = p
    Mcc = sigmac**2
    Mrr = sigmar**2
    Mrc = rho**2*Mcc*Mrr
    M20 = Mrr + Mcc
    M22 = complex(Mcc - Mrr,2*Mrc)
    whisker_g2d = np.sqrt(np.abs(M22))
    lambdap = 0.5*(M20 + abs(M22))
    lambdam = 0.5*(M20 - abs(M22))
    fwhm_g2d = np.sqrt(2.*np.log(2.))*(np.sqrt(lambdap)+np.sqrt(lambdam))
    return A, B, whisker_g2d, fwhm_g2d

def wfwhm(img,sigma):
    """
    This code calculate the fwhm and wisker length defined as (M22.real^2 + M22.imag^2)^{1/4} using the weighted moments method.
    input: 
         data: 2d stamp image
         sigma: std of the Gaussian weight Kernel in pixel
    """
    nrow,ncol=img.shape
    Isum = img.sum()
    Icol = img.sum(axis=0) # sum over all rows
    Irow = img.sum(axis=1) # sum over all cols
    colgrid = np.arange(ncol)
    rowgrid = np.arange(nrow)
    rowmean=np.sum(rowgrid*Irow)/Isum
    colmean=np.sum(colgrid*Icol)/Isum
    ROW,COL=np.indices((nrow,ncol))
    maxItr = 50
    EP = 0.0001
    for i in range(maxItr):
        wrmat = wr(ROW,COL,rowmean,colmean,sigma)
        IWmat = img*wrmat
        IWcol = IWmat.sum(axis=0)
        IWrow = IWmat.sum(axis=1)
        IWsum = IWmat.sum()
        drowmean = np.sum((rowgrid-rowmean)*IWrow)/IWsum
        dcolmean = np.sum((colgrid-colmean)*IWcol)/IWsum
        rowmean = rowmean+2.*drowmean
        colmean = colmean+2.*dcolmean
        if drowmean**2+dcolmean**2 <= EP:
            break
    rowgrid = rowgrid - rowmean # centered
    colgrid = colgrid - colmean
    Mrr = np.sum(rowgrid**2*IWrow)/IWsum
    Mcc = np.sum(colgrid**2*IWcol)/IWsum
    Mrc = np.sum(np.outer(rowgrid,colgrid)*IWmat)/IWsum
    Cm = np.matrix([[Mcc,Mrc],[Mrc,Mrr]]) # cov matrix from measurement
    Cw = np.matrix([[sigma**2,0.],[0.,sigma**2]])# cov matrix from weight
    Cimg = (Cm.I - Cw.I).I #cov matrix after subtract weight
    Mcc = Cimg[0,0]
    Mrr = Cimg[1,1]
    Mrc = Cimg[0,1]
    M20 = Mrr + Mcc
    M22 = complex(Mcc - Mrr,2*Mrc)
    e1 = M22.real/M20.real
    e2 = M22.imag/M20.real
    whiskerLength = np.sqrt(np.abs(M22))
    lambdap = 0.5*(M20 + abs(M22))
    lambdam = 0.5*(M20 - abs(M22))
    fwhmw = np.sqrt(2.*np.log(2.))*(np.sqrt(lambdap)+np.sqrt(lambdam))
    return e1,e2,whiskerLength,fwhmw 


def mfwhm(img=None):
    """
    measure the fwhm using Moffat fit.
    output: 
    """
    npix = img.shape[0]
    rowCen,colCen = adaptiveCentroid(img,1.1/scale)
    row,col = np.mgrid[0:npix,0:npix]
    row = row - rowCen
    col = col - colCen
    radius = np.sqrt(row**2+col**2)
    A0,alpha0 = moments(img)
    beta0=1.5
    B0 = 0.
    p0=np.array([alpha0,beta0,A0, B0])
    img = img.flatten()
    ok = img >0
    img = img[ok]
    radius = radius.flatten()
    radius = radius[ok]
    def residualm(p,r,I):
        alpha,beta,A,B = p
        Ierr = np.sqrt(abs(I))
        res = (mprofile(radius,alpha,beta,A,B) - I)/Ierr
        return res 
    p = leastsq(residualm,p0,args=(radius,img))[0]
    alpha,beta,A,B = p
    fwhm_moffat= 2. * abs(alpha) * np.sqrt(2.**(1./beta)-1)
    return alpha,beta,A,B,fwhm_moffat

    
def get_fwhm_whisker(stampImg=None,sigma=1.1/scale):
    """
    Calcualte the fwhm, whisker using various approach. 
    return the results in arcsec. 
    output: [whisker_weighted_moments, whisker_Amoments]
            [fwhm_weighted, fwhm_Amoments,fwhm_moffat, fwhm_gauss,fwhm_sech2]
    """
    if stampImg.shape[0] == stampImg.shape[1] and stampImg.shape[1] != 0:
        stampImg = stampImg - np.min(stampImg) # this renormalize the negative pixel values
        npix = stampImg.shape[0]
        mfit = mfwhm(stampImg)
        gfit = gfwhm(stampImg)
        s2fit = s2fwhm(stampImg)
        g2dfit = g2dfwhm(stampImg)
        wfit = wfwhm(stampImg,sigma=sigma)
        fwhm = np.array([wfit[3],g2dfit[3],mfit[4],gfit[3],s2fit[3]])*scale
        whisker = np.array([wfit[2],g2dfit[2]])*scale
        fwhm[np.isnan(fwhm)]=-999
        whisker[np.isnan(whisker)]=-999
    else:
        fwhm = np.array([-999,-999,-999,-999,-999])
        whisker = np.array([-999,-999])
    return whisker, fwhm

def get_fwhm_whisker_list(stampImgList=None):
    """
    Calcualte the fwhm, whisker using various approach from the list of stamp image.
    return the results in arcsec. 
    output: [whisker_weighted_moments, whisker_Amoments]
            [fwhm_weighted, fwhm_Amoments,fwhm_moffat, fwhm_gauss,fwhm_sech2]
    """
    n=len(stampImgList)
    whisker=[]
    fwhm=[]
    for i in range(n):
        print i
        whker,fw = get_fwhm_whisker(stampImgList[i])
        whisker.append(whker)
        fwhm.append(fw)
    whisker = np.array(whisker)
    fwhm = np.array(fwhm)
    return whisker, fwhm

def fwhm_whisker_plot(stampImgList=None):
    whk,fwhm = get_fwhm_whisker_list(stampImgList)
    whk=list(whk.T)
    fwh=list(fwhm.T)
    pl.figure(figsize=(7,5))
    pl.boxplot(whk)
    pl.hlines(0.2,0,3,linestyle='solid',color='g')
    pl.ylim(0.,.4)
    pl.xticks(np.arange(1,3),['whisker_Wmoments','whisker_Amoments'])
    pl.figure(figsize=(12,5))
    pl.boxplot(fwh)
    pl.ylim(0.4,1.5)
    pl.hlines(0.9,0,6,linestyle='solid',color='g')
    pl.xticks(np.arange(1,6),['fwhm_weighted', 'fwhm_Amoments','fwhm_moffat', 'fwhm_gauss','fwhm_sech2'])
    return '-----done !----'


def fwhm_whisker_des_plot(stampImgList=None,whkSex=None,fwhmSex=None):
    whk,fwhm = get_fwhm_whisker_list(stampImgList)
    whk=list(whk.T)
    fwh=list(fwhm.T)
    fwh.append(fwhmSex)
    whk.append(whkSex)
    pl.figure(figsize=(15,10))
    pl.subplot(2,1,1)
    pl.boxplot(whk)
    pl.hlines(0.2,0,4,linestyle='solid',color='g')
    pl.ylim(np.median(whk[2])-0.3,np.median(whk[2])+0.3)
    pl.grid()
    pl.xticks(np.arange(1,4),['whisker_Wmoments','whisker_Amoments','whisker_sx'])
    pl.subplot(2,1,2)
    pl.boxplot(fwh)
    pl.ylim(np.median(fwh[5])-2,np.median(fwh[5])+2)
    pl.grid()
    pl.hlines(0.9,0,7,linestyle='solid',color='g')
    pl.xticks(np.arange(1,7),['fwhm_weighted', 'fwhm_Amoments','fwhm_moffat', 'fwhm_gauss','fwhm_sech2','fwhm_sx'])
    return '-----done !----'



def dispStamp(stampImg=None,sigma=1.08/scale):
    if stampImg.shape[0] != stampImg.shape[1]:
        sys.exit('bad stamp image')
    npix = stampImg.shape[0]
    pl.figure(figsize=(18,6))
    pl.subplot(1,3,1)
    pl.matshow(stampImg,fignum=False)
    #pl.contour(stampImg,nlevels=20)
    mfit = mfwhm(stampImg)
    gfit = gfwhm(stampImg)
    s2fit = s2fwhm(stampImg)
    g2dfit = g2dfwhm(stampImg)
    wfit = wfwhm(stampImg,sigma=sigma)
    pl.xlabel('Pixel')
    pl.ylabel('Pixel')
    pl.grid(color='y')
    rowCen,colCen = adaptiveCentroid(data=stampImg,sigma=sigma)
    M20, M22, M31, M33 =complexMoments(stampImg,sigma=sigma)
    print M20, M22, M31, M33
    e1 = M22.real/M20.real
    e2 = M22.imag/M20.real
    pl.figtext(0.15,0.8, 'e1: '+str(round(e1,3)) + ',  e2: '+str(round(e2,3)), color='r')
    pl.figtext(0.15,0.75, 'rowCen: '+str(round(rowCen,4)) + ',  colCen: '+str(round(colCen,4)), color='r')
    pl.figtext(0.15,0.7, 'PSF whisker_Wmoments: '+str(round(wfit[2]*scale,4))+' [arcsec]', color='r')
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
    pl.vlines(wfit[3]/2.,0, halfmax*4,linestyle='solid',colors='b',label='Weighted Moments')
    pl.vlines(mfit[4]/2.,0, halfmax*4,linestyle='solid',colors='r')
    pl.vlines(gfit[3]/2.,0, halfmax*4,linestyle='solid',colors='g')
    pl.vlines(s2fit[3]/2.,0, halfmax*4,linestyle='solid',colors='m')
    pl.vlines(g2dfit[3]/2.,0, halfmax*4,linestyle='solid',colors='c')
    pl.plot(radius,mprofile(radius,mfit[0],mfit[1],mfit[2],mfit[3]),'r-',label='Moffat Fit')
    pl.plot(radius,gprofile(radius,gfit[0],gfit[1],gfit[2]),'g-',label='Gaussian Fit')
    pl.plot(radius,s2profile(radius,s2fit[0],s2fit[1],s2fit[2]),'m-',label='Sech2 Fit')
    pl.legend(loc='best')
    pl.ylim(0,halfmax*4)
    pl.xlim(0,npix/2.) 
    pl.xlabel('Radius [pixels]')
    pl.ylabel('Mean counts [ADU]')
    pl.title('Radial profile')
    pl.figtext(0.65,0.7,'Gaussian Weight '+r'$\sigma$: '+str(round(sigma*scale,3))+ ' arcsec',color='r')
    pl.figtext(0.65,0.6,'FWHM_Gaussian: '+str(round(gfit[3]*scale,3))+ ' arcsec')
    pl.figtext(0.65,0.55,'FWHM_Moffat: '+str(round(mfit[4]*scale,3))+ ' arcsec')
    pl.figtext(0.65,0.5,'FWHM_Sech2: '+str(round(s2fit[3]*scale,3))+ ' arcsec')
    pl.figtext(0.65,0.45,'FWHM_Wmoments: '+str(round(wfit[3]*scale,3))+ ' arcsec') 
    pl.figtext(0.65,0.4,'FWHM_Amoments: '+str(round(g2dfit[3]*scale,3))+ ' arcsec')
    pl.figtext(0.65,0.35,'M20: '+str(round(M20,5))+ ' pix')
    pl.figtext(0.65,0.3,'M22.real: '+str(round(M22.real,5))+ ' pix')
    pl.figtext(0.8,0.3,'M22.imag: '+str(round(M22.imag,5))+ ' pix')
    pl.figtext(0.65,0.25,'M31.real: '+str(round(M31.real,5))+ ' pix')
    pl.figtext(0.8,0.25,'M31.imag: '+str(round(M31.imag,5))+ ' pix')
    pl.figtext(0.65,0.2,'M33.real: '+str(round(M33.real,5))+ ' pix')
    pl.figtext(0.8,0.2,'M33.imag: '+str(round(M33.imag,5))+ ' pix')
    return '---- Done! ----'
   
def dispStampList(stampImgList=None,sigma=1.08/scale):
    if sigma == None:
        print 'syntax: dispStampList(stampImgList,sigma)'
        sys.exit()
    Nstamp = len(stampImgList)
    for i in range(Nstamp):
        t=dispStamp(stampImg=stampImgList[i],sigma=sigma)
        raw_input('--- hit the enter key to proceed ---')
        pl.close()
    return ' ---- Done ! ----'

if __name__ == "__main__":
    from ImgQuality import *
    pl.ion()
    dr = '/home/jghao/research/data/des_optics_psf/dc6b_image/goodseeing/decam--28--49-r-1/'
    starfile='/home/jghao/research/data/des_optics_psf/dc6b_image/goodseeing/catfile/decam_-27.72186_-48.60000-objects.fit'
    extension = np.arange(1,10)
    stamp=[]
    for ext in extension:
        print ext
        if ext < 10:
            ext = '0'+str(ext)
        else:
            ext = str(ext)
        imgname= 'decam--28--49-r-1_'+ext+'.fits.fz'
        bkgname = 'decam--28--49-r-1_'+ext+'_bkg.fits.fz'
        data = pf.getdata(dr+imgname) - pf.getdata(dr+bkgname)
        xc = pf.getdata(starfile,3*(int(ext)-1)+1).xccd
        yc = pf.getdata(starfile,3*(int(ext)-1)+1).yccd
        rmag = pf.getdata(starfile,3*(int(ext)-1)+1).mag_3
        ok = (rmag > 16.5)*(rmag < 17.5)
        xc=xc[ok]
        yc = yc[ok]
        stamp=stamp+getStamp(data=data,xcoord=xc,ycoord=yc,Npix =25)
    fwhm_whisker_plot(stamp)
    pl.savefig()
