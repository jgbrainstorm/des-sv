#------------------------------------------------------
# this file contains functions used for the decamImgAnalyzer
# the functions are adopted from ImgQuality.py and psfFocus.py
# this file contains only those useful functions for the decamImgAnalyzer 
# some other functions not used are still in the ImgQuality.py and
# psfFocus.py. 
# created by Jiangang Hao, 11/13/2012

import sys,os
sys.path.append('/usr/remote/user/sispi/jiangang/des-sv')
sys.path.append('/usr/remote/user/sispi/jiangang/decam-fermi')
sys.path.append('/usr/remote/user/sispi/jiangang/lib/python')


try:
    import numpy as np
    import pyfits as pf
    import scipy.ndimage as nd
    import pylab as pl
    from scipy.optimize import leastsq
    from DECamCCD_def import *
    from scipy.misc import factorial as fac
except ImportError:
    print "Error: missing one of the libraries (numpy, pyfits, scipy, matplotlib)"
    sys.exit()

scale=0.27

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
        if np.sum(Img) > 0:
            stampImg.append(Img)
    return stampImg


def momentsold(data):
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

def moments(data):
    """
    Returns (height, and width)
    the gaussian parameters of a 2D distribution by calculating its
    moments
    """
    sigma=2. # 2 pix as weight kernel
    nrow,ncol=data.shape
    Isum = data.sum()
    Icol = data.sum(axis=0) # sum over all rows
    Irow = data.sum(axis=1) # sum over all cols
    colgrid = np.arange(ncol)
    rowgrid = np.arange(nrow)
    rowmean=np.sum(rowgrid*Irow)/Isum
    colmean=np.sum(colgrid*Icol)/Isum
    ROW,COL=np.indices((nrow,ncol))
    wrmat = wr(ROW,COL,rowmean,colmean,sigma)
    IWmat = data*wrmat
    IWcol = IWmat.sum(axis=0)
    IWrow = IWmat.sum(axis=1)
    IWsum = IWmat.sum()
    rowgrid = rowgrid - rowmean # centered
    colgrid = colgrid - colmean
    Mrr = np.sum(rowgrid**2*IWrow)/IWsum
    Mcc = np.sum(colgrid**2*IWcol)/IWsum
    width = np.sqrt((Mrr+Mcc)*0.5)
    datasrt = np.sort(data.flatten())
    height = np.median(datasrt[-20:-1])
    return height,width


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
    #print Mrrr, Mccc, Mrrc, Mrcc
    M20 = Mrr + Mcc
    M22 = complex(Mcc - Mrr,2*Mrc)
    M31 = complex(3*Mc - (Mccc+Mrrc)/sigma**2, 3*Mr - (Mrcc + Mrrr)/sigma**2)
    M33 = complex(Mccc-3*Mrrc, 3.*Mrcc - Mrrr)
    return M20, M22, M31, M33

def AcomplexMoments(img,sigma=1.1/scale):
    """
    calculate M20, M22 using the adaptive moments.
    x is col, y is row
    sigmax -> sigmac, sigmay -> sigmar
    """
    npix = img.shape[0]
    rowCen,colCen = adaptiveCentroid(img,sigma)
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
    return M20, M22



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

    
def get_fwhm_whisker(stampImg=None,bkg = None,sigma=1.1/scale):
    """
    Calcualte the fwhm, whisker using various approach. 
    return the results in arcsec. 
    output: [whisker_weighted_moments, whisker_Amoments]
            [fwhm_weighted, fwhm_Amoments,fwhm_moffat, fwhm_gauss,fwhm_sech2]
    """
    if bkg != None:
        stampImg = stampImg - bkg
    if stampImg.shape[0] == stampImg.shape[1] and stampImg.shape[1] != 0 and np.sum(stampImg) > 0:
        npix = stampImg.shape[0]
        try:
            mfit = mfwhm(stampImg)
            gfit = gfwhm(stampImg)
            s2fit = s2fwhm(stampImg)
            g2dfit = g2dfwhm(stampImg)
            wfit = wfwhm(stampImg,sigma=sigma)
            fwhm = np.array([wfit[3],g2dfit[3],mfit[4],gfit[3],s2fit[3]])*scale
            whisker = np.array([wfit[2],g2dfit[2]])*scale
        except ValueError:
            fwhm = np.array([-999,-999,-999,-999,-999])
            whisker = np.array([-999,-999])
            pass
        #fwhm[np.isnan(fwhm)]=-999
        #whisker[np.isnan(whisker)]=-999
    else:
        fwhm = np.array([-999,-999,-999,-999,-999])
        whisker = np.array([-999,-999])
    return whisker, fwhm

def get_fwhm_whisker_list(stampImgList=None,bkgList=None,sigma=1.1/scale):
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
        whker,fw = get_fwhm_whisker(stampImgList[i],bkgList[i],sigma=sigma)
        if len(whker[whker>1])>0 or len(fw[fw==-999])>0:
            continue
        else:
            whisker.append(whker)
            fwhm.append(fw)
    whisker = np.array(whisker)
    fwhm = np.array(fwhm)
    return whisker, fwhm

def fwhm_whisker_plot(stampImgList=None,bkgList=None,sigma=1.1/scale):
    whk,fwhm = get_fwhm_whisker_list(stampImgList,bkgList,sigma=sigma)
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


def fwhm_whisker_des_plot(stampImgList=None,bkgList=None,whkSex=None,fwhmSex=None,sigma=1.1/scale,dimmfwhm=None):
    whk,fwhm = get_fwhm_whisker_list(stampImgList,bkgList,sigma=sigma)
    whk=list(whk.T)
    fwh=list(fwhm.T)
    fwh.append(fwhmSex)
    whk.append(whkSex)
    pl.figure(figsize=(15,10))
    pl.subplot(2,1,1)
    pl.boxplot(whk)
    pl.hlines(0.2,0,4,linestyle='solid',color='g')
    pl.ylim(np.median(whk[2])-0.3,np.median(whk[2])+0.6)
    pl.grid()
    pl.xticks(np.arange(1,4),['whisker_Wmoments','whisker_Amoments','whisker_sx'])
    if dimmfwhm != None:
        pl.title('DIMM Seeing FWHM: '+str(round(dimmfwhm,3)) +'(arcsec)    sqrt(DIMM fwhm^2 + 0.55^2): '+str(round(np.sqrt(dimmfwhm**2 + 0.55**2),3)))
    pl.subplot(2,1,2)
    pl.boxplot(fwh)
    pl.ylim(0,np.median(fwh[5])+2)
    pl.grid()
    pl.hlines(0.9,0,7,linestyle='solid',color='g')
    pl.xticks(np.arange(1,7),['fwhm_weighted', 'fwhm_Amoments','fwhm_moffat', 'fwhm_gauss','fwhm_sech2','fwhm_sx'])
    return fwh,whk

def dispM202Coeff(betaAll=None,betaErrAll=None,hexinfo=None):
    ind = np.arange(len(betaAll[0]))
    momname = ('M20','M22.Real','M22.imag')
    fmtarr = ['bo-','ro-','go-']
    if betaErrAll == None:
        betaErrAll = np.zeros(len(ind))
    pl.figure(figsize=(17,7))
    for i in range(3):
        pl.subplot(4,1,i+1)
        pl.errorbar(ind[1:],betaAll[i][1:],yerr = betaErrAll[i][1:],fmt=fmtarr[i])
        pl.grid()
        pl.xlim(-1,len(betaAll[i])+1)
        if i==0 and hexinfo != None:
            pl.title('Hexapod deviation from perfect: x'+str(round(hexinfo[0],2))+'  y:'+str(round(hexinfo[1],2))+'  z:'+str(round(hexinfo[2],2))+'  xtilt:'+str(round(hexinfo[3],2))+'  ytilt:'+str(round(hexinfo[4],2)))
        #pl.ylim(min(betaAll[i][1:])-0.01,max(betaAll[i][1:])+0.01)
        pl.xticks(ind,('','','','','','','','','','','','','','','','','','','',''))
        pl.ylim(-0.3,0.3)
        pl.ylabel(momname[i])
    pl.xticks(ind,('Piston','Tip','Tilt','Astignism','Defocus','Astignism','Trefoil','Coma','Coma','Trefoil','Ashtray','Astigm.5th','Spherical','Astigm.5th','Ashtray','16','17','18','19','20'),rotation=90)
    pl.xlabel('Zernike Coefficients')
    return '---done!---'


def dispStamp(stampImg=None,bkg=None,sigma=1.08/scale,mag=None,rad=None,ok=None,expid=None,detector=None):
    if stampImg.shape[0] != stampImg.shape[1]:
        sys.exit('bad stamp image')
    if bkg != None:
        stampImg=stampImg - bkg
    npix = stampImg.shape[0]
    pl.figure(figsize=(10,10))
    pl.subplot(2,2,1)
    pl.plot(mag,rad,'b.')
    pl.plot(mag[ok],rad[ok],'r.')
    pl.xlabel('mag')
    pl.ylabel('radius')
    pl.ylim(0,13)
    pl.xlim(-20,0)
    pl.title('Exposure: '+expid+'   CCD: '+detector)
    pl.subplot(2,2,2)
    pl.matshow(stampImg,fignum=False)
    pl.colorbar(fraction=0.1)
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
    pl.figtext(0.15,0.85, 'e1: '+str(round(e1,3)) + ',  e2: '+str(round(e2,3)), color='r')
    pl.figtext(0.15,0.82, 'whisker_Wmoments: '+str(round(wfit[2]*scale,4))+' [arcsec]', color='r')
    pl.figtext(0.15,0.79, 'whisker_Amoments: '+str(round(g2dfit[2]*scale,4))+' [arcsec]', color='r')
    pl.subplot(2,2,3)
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
    pl.figtext(0.55,0.45,'Gaussian Weight '+r'$\sigma$: '+str(round(sigma*scale,3))+ ' arcsec',color='r')
    pl.figtext(0.55,0.42,'FWHM_Gaussian: '+str(round(gfit[3]*scale,3))+ ' arcsec')
    pl.figtext(0.55,0.39,'FWHM_Moffat: '+str(round(mfit[4]*scale,3))+ ' arcsec')
    pl.figtext(0.55,0.36,'FWHM_Sech2: '+str(round(s2fit[3]*scale,3))+ ' arcsec')
    pl.figtext(0.55,0.33,'FWHM_Wmoments: '+str(round(wfit[3]*scale,3))+ ' arcsec') 
    pl.figtext(0.55,0.30,'FWHM_Amoments: '+str(round(g2dfit[3]*scale,3))+ ' arcsec')
    pl.figtext(0.55,0.27,'M20: '+str(round(M20,5))+ ' pix^2')
    pl.figtext(0.55,0.24,'M22.real: '+str(round(M22.real,5))+ ' pix^2')
    pl.figtext(0.55,0.21,'M22.imag: '+str(round(M22.imag,5))+ ' pix^2')
    return '---- Done! ----'
   
def dispStampList(stampImgList=None,bkgList=None,sigma=1.08/scale,mag=None,rad=None,ok=None,expid=None,detector=None):
    if sigma == None:
        print 'syntax: dispStampList(stampImgList,bkgList,sigma)'
        sys.exit()
    Nstamp = len(stampImgList)
    for i in range(Nstamp):
        t=dispStamp(stampImgList[i],bkgList[i],sigma,mag,rad,ok,expid,detector)
        raw_input('--- hit the enter key to proceed ---')
        pl.close()
    return ' ---- Done ! ----'

def measure_stamp_moments(stamp,bkg=None,sigma=4.,adaptive=False):
    """
    measure the moments of stamp image on one CCD
    return the median moments
    """
    Nobj = len(stamp)
    M20=np.zeros(Nobj)
    M22=np.zeros(Nobj).astype(complex)
    M31=np.zeros(Nobj).astype(complex)
    M33=np.zeros(Nobj).astype(complex)
    for i in range(Nobj):
        if bkg == None:
            M20[i],M22[i],M31[i],M33[i]=complexMoments(data=stamp[i],sigma=sigma)
        else:
            data = stamp[i]-bkg[i]
            if data.sum > 0.:
                if adaptive == False:
                    M20[i],M22[i],M31[i],M33[i]=complexMoments(data=data,sigma=sigma)               
                else:
                    M20[i],M22[i] = AcomplexMoments(img=data,sigma = sigma)
    return [np.median(M20), np.median(M22), np.median(M31), np.median(M33)]


def zernike_rad(m, n, rho):
    """
    Calculate the radial component of Zernike polynomial (m, n) 
    given a grid of radial coordinates rho.
    """
    if (n < 0 or m < 0 or abs(m) > n):
        raise ValueError
    if ((n-m) % 2):
        return rho*0.0
    pre_fac = lambda k: (-1.0)**k * fac(n-k) / ( fac(k) * fac( (n+m)/2.0 - k ) * fac( (n-m)/2.0 - k ) )
    return sum(pre_fac(k) * rho**(n-2.0*k) for k in xrange((n-m)/2+1))

def zernike(m, n, rho, phi):
    """
    Calculate Zernike polynomial (m, n) given a grid of radial
    coordinates rho and azimuthal coordinates phi.
    """
    if (m > 0): return zernike_rad(m, n, rho) * np.cos(m * phi)
    if (m < 0): return zernike_rad(-m, n, rho) * np.sin(-m * phi)
    return zernike_rad(0, n, rho)


def zernikel(j, rho, phi):
    """
    Calculate Zernike polynomial with Noll coordinate j given a grid of radial coordinates rho and azimuthal coordinates phi.
    """
    n = 0
    while (j > n):
        n += 1
        j -= n
    m = -n+2*j
    return zernike(m, n, rho, phi)


def zernikeFit(x, y, z,max_rad=225.,cm=[0,0],max_order=20):
    """
    Fit a set of x, y, z data to a zernike polynomial with the least square fitting. Note that here x, y, z are all 1 dim array. Here the max_rad is by default equal to 225 mm, the size of the decam focal plane.
    It will return the beta and the adjusted R2
    """
    x = x - cm[0]
    y = y - cm[1]
    n = len(x)
    p = max_order
    rho = np.sqrt(x**2+y**2)/max_rad #normalize to unit circle.
    phi = np.arctan2(y,x)
    dataX = []
    ok = rho <= 1.
    for j in range(max_order):
        dataX.append(zernikel(j,rho[ok],phi[ok]))
    dataX=np.array(dataX).T
    beta,SSE,rank,sing = np.linalg.lstsq(dataX,z[ok])# SSE is the residual sum square
    sigma = np.sqrt(SSE/(n-p))
    betaErr = sigma/np.dot(dataX.T,dataX).diagonal()
    SST = np.var(z[ok])*(len(z[ok])-1)# SST is the sum((z_i - mean(z))^2)
    R2 = 1 - SSE/SST
    R2adj = 1-(1-R2)*(len(z[ok])-1)/(len(z[ok])-max_order)# adjusted R2 for quality of fit.             
    fitted = np.dot(dataX,beta) # fitted value
    return beta,betaErr, R2adj,fitted



def measure_stamp_coeff(data = None, zernike_max_order=20):
    """
    the convention of data is: x, y, M20, M22, M31, M33
    """
    betaAll=[]
    betaErrAll=[]
    R2adjAll=[]
    beta,betaErr,R2_adj = zernikeFit(data[:,0].real,data[:,1].real,data[:,2].real,max_order=zernike_max_order)
    betaAll.append(beta)
    betaErrAll.append(betaErr)
    R2adjAll.append(R2_adj)
    for i in range(3,6):
        beta,betaErr,R2_adj = zernikeFit(data[:,0].real,data[:,1].real,data[:,i].real,max_order=zernike_max_order)
        betaAll.append(beta)
        betaErrAll.append(betaErr)
        R2adjAll.append(R2_adj)
        beta,betaErr,R2_adj = zernikeFit(data[:,0].real,data[:,1].real,data[:,i].imag,max_order=zernike_max_order)
        betaAll.append(beta)
        betaErrAll.append(betaErr)
        R2adjAll.append(R2_adj)
    betaAll = np.array(betaAll)
    betaErrAll = np.array(betaErrAll)
    R2adjAll = np.array(R2adjAll)
    return betaAll,betaErrAll, R2adjAll


def display_moments(data=None):
    # remove the mean for all moments
    data = subMeanAll(data)
    pl.figure(figsize=(11,11))
    pl.subplot(2,2,1)
    phi22 = 0.5*np.arctan2(data[:,3].imag,data[:,3].real)
    x = data[:,0].real
    y = data[:,1].real
    phi22[x<0] = phi22+np.deg2rad(180)
    u = np.abs(data[:,3])*np.cos(phi22)
    v = np.abs(data[:,3])*np.sin(phi22)
    qvr = pl.quiver(x,y,u,v,width = 0.004, color='r',pivot='middle',headwidth=0.,headlength=0.,headaxislength=0.,scale_units='width')
    qk = pl.quiverkey(qvr, -150,-240,np.max(np.sqrt(u**2+v**2)),str(round(np.max(np.sqrt(u**2+v**2)),3))+' pix^2',coordinates='data',color='blue')
    pl.plot(x,y,'b,')
    pl.xlim(-250,250)
    pl.ylim(-250,250)
    pl.grid(color='g')
    pl.xlabel('X [mm] (WEST)')
    pl.ylabel('Y [mm] (NORTH)')
    pl.title('M22')
    pl.subplot(2,2,2)
    phi31 = np.arctan2(data[:,4].imag,data[:,4].real)
    u = np.abs(data[:,4])*np.cos(phi31)
    v = np.abs(data[:,4])*np.sin(phi31)
    qvr=pl.quiver(x,y,u,v,width=0.003,color='r',pivot='middle',headwidth=4)
    qk = pl.quiverkey(qvr, -150,-240,np.max(np.sqrt(u**2+v**2)),str(round(np.max(np.sqrt(u**2+v**2)),3))+' pix^3',coordinates='data',color='blue')
    pl.plot(x,y,'b,')
    pl.xlim(-250,250)
    pl.ylim(-250,250)
    pl.grid(color='g')
    pl.xlabel('X [mm] (WEST)')
    pl.ylabel('Y [mm] (NORTH)')
    pl.title('M31')
    pl.subplot(2,2,3)
    phi33 = np.arctan2(data[:,5].imag,data[:,5].real)/3.
    u = np.abs(data[:,5])*np.cos(phi33)
    v = np.abs(data[:,5])*np.sin(phi33)
    pl.quiver(x,y,u,v,width=0.003,color='r',headwidth=4)
    u = np.abs(data[:,5])*np.cos(phi33+np.deg2rad(120))
    v = np.abs(data[:,5])*np.sin(phi33+np.deg2rad(120))
    pl.quiver(x,y,u,v,width=0.003,color='r',headwidth=4)
    u = np.abs(data[:,5])*np.cos(phi33+np.deg2rad(240))
    v = np.abs(data[:,5])*np.sin(phi33+np.deg2rad(240))
    qvr=pl.quiver(x,y,u,v,width=0.003,color='r',headwidth=4)
    qk = pl.quiverkey(qvr, -150,-240,np.max(np.sqrt(u**2+v**2)),str(round(np.max(np.sqrt(u**2+v**2)),3))+' pix^3',coordinates='data',color='blue')
    pl.plot(x,y,'b,')
    pl.xlim(-250,250)
    pl.ylim(-250,250)
    pl.grid(color='g')
    pl.xlabel('X [mm] (WEST)')
    pl.ylabel('Y [mm] (NORTH)')
    pl.title('M33')
    pl.subplot(2,2,4)
    m20sqr = np.sqrt(data[:,2].real)
    x = data[:,0].real
    y = data[:,1].real
    m20sqr_med = np.median(m20sqr)
    m20sqr_diff = m20sqr - m20sqr_med
    m20sqr_diff_absmed = np.median(np.abs(m20sqr_diff))
    plotScale = 1./m20sqr_diff_absmed*100
    pos = m20sqr_diff >=0
    neg = m20sqr_diff < 0
    pl.scatter(x[pos],y[pos],s=m20sqr_diff[pos]*plotScale,c='r',alpha=0.5)
    pl.scatter(x[neg],y[neg],s=-m20sqr_diff[neg]*plotScale,c='b',alpha=0.5)
    pl.scatter(-230,-210,s=m20sqr_diff_absmed*plotScale,c='b',alpha=0.5)
    pl.text(-200,-215,'-'+str(round(m20sqr_diff_absmed,6))+' pix')
    pl.scatter(-230,-230,s=m20sqr_diff_absmed*plotScale,c='r',alpha=0.5)
    pl.text(-200,-235,str(round(m20sqr_diff_absmed,6))+' pix')
    pl.plot(x,y,'y,')
    pl.grid(color='g')
    pl.xlim(-250,250)
    pl.ylim(-250,250)
    pl.xlabel('X [mm] (WEST)')
    pl.ylabel('Y [mm] (NORTH)')
    pl.title('median '+r'$\sqrt{M20}$: '+str(round(scale*m20sqr_med,3))+' [arcsec]')
    return '---done!--'

def display_2nd_moments(data=None):
    # remove the mean for all moments
    data = subMeanAll(data)
    pl.figure(figsize=(11,5.5))
    pl.subplot(1,2,1)
    phi22 = 0.5*np.arctan2(data[:,3].imag,data[:,3].real)
    x = data[:,0].real
    y = data[:,1].real
    phi22[x<0] = phi22+np.deg2rad(180)
    u = np.abs(data[:,3])*np.cos(phi22)
    v = np.abs(data[:,3])*np.sin(phi22)
    qvr = pl.quiver(x,y,u,v,width = 0.004, color='r',pivot='middle',headwidth=0.,headlength=0.,headaxislength=0.,scale_units='width')
    #qk = pl.quiverkey(qvr, -150,-240,np.max(np.sqrt(u**2+v**2)),str(round(np.max(np.sqrt(u**2+v**2)),3))+' pix^2',coordinates='data',color='blue')
    qk = pl.quiverkey(qvr, -150,-240,0.3,str(0.3)+' pix^2',coordinates='data',color='blue')
    pl.plot(x,y,'b,')
    pl.xlim(-250,250)
    pl.ylim(-250,250)
    pl.grid(color='g')
    pl.xlabel('X [mm] (WEST)')
    pl.ylabel('Y [mm] (NORTH)')
    pl.title('M22')
    pl.subplot(1,2,2)
    m20sqr = np.sqrt(data[:,2].real)
    x = data[:,0].real
    y = data[:,1].real
    m20sqr_med = np.median(m20sqr)
    m20sqr_diff = m20sqr - m20sqr_med
    m20sqr_diff_absmed = np.median(np.abs(m20sqr_diff))
    plotScale = 1./m20sqr_diff_absmed*100
    pos = m20sqr_diff >=0
    neg = m20sqr_diff < 0
    pl.scatter(x[pos],y[pos],s=m20sqr_diff[pos]*plotScale,c='r',alpha=0.5)
    pl.scatter(x[neg],y[neg],s=-m20sqr_diff[neg]*plotScale,c='b',alpha=0.5)
    #pl.scatter(-230,-210,s=m20sqr_diff_absmed*plotScale,c='b',alpha=0.5)
    #pl.text(-200,-215,'-'+str(round(m20sqr_diff_absmed,6))+' pix')
    #pl.scatter(-230,-230,s=m20sqr_diff_absmed*plotScale,c='r',alpha=0.5)
    #pl.text(-200,-235,str(round(m20sqr_diff_absmed,6))+' pix')
    pl.scatter(-230,-210,s=0.01*plotScale,c='b',alpha=0.5)
    pl.text(-200,-215,'-'+str(0.01)+' pix')
    pl.scatter(-230,-230,s=0.01*plotScale,c='r',alpha=0.5)
    pl.text(-200,-235,str(0.01)+' pix')
    pl.plot(x,y,'y,')
    pl.grid(color='g')
    pl.xlim(-250,250)
    pl.ylim(-250,250)
    pl.xlabel('X [mm] (WEST)')
    pl.ylabel('Y [mm] (NORTH)')
    pl.title('median '+r'$\sqrt{M20}$: '+str(round(scale*m20sqr_med,3))+' [arcsec]')
    return '---done!--'


def showZernike(beta=None,betaErr=None,gridsize = 1, max_rad = 1,significance=False):
    """
    significance shows how significant the coefficients are constrained. 
    """
    x,y = np.meshgrid(np.arange(-gridsize,gridsize,0.01),np.arange(-gridsize,gridsize,0.01))
    rho = np.sqrt(x**2+y**2)
    phi = np.arctan2(y,x)
    ok = rho < max_rad
    if significance != False:
        sigIdx = np.abs(beta)/betaErr >= significance
        beta = beta[sigIdx]
    nn = len(beta)
    znk=0
    for j in range(nn):
        znk = znk + beta[j]*zernikel(j,rho,phi)*ok
    pl.imshow(znk)
    return znk


def display_zernike(data,zernike_max_order=20):
    
    colnames = ['x','y','M20','M22','M31','M33']
    data=np.array(data)
    pl.figure(figsize=(15,15))
    betaAll=[]
    betaErrAll=[]
    R2adjAll=[]
    pl.subplot(3,3,1)
    beta,betaErr,R2_adj = zernikeFit(data[:,0].real,data[:,1].real,data[:,2].real,max_order=zernike_max_order)
    betaAll.append(beta)
    betaErrAll.append(betaErr)
    R2adjAll.append(R2_adj)
    znk=showZernike(beta=beta)
    pl.colorbar()
    pl.title(colnames[2])
    for i in range(3,6):
        pl.subplot(3,3,2*i-4)
        beta,betaErr,R2_adj = zernikeFit(data[:,0].real,data[:,1].real,data[:,i].real,max_order=zernike_max_order)
        betaAll.append(beta)
        betaErrAll.append(betaErr)
        R2adjAll.append(R2_adj)
        znk=showZernike(beta=beta)
        pl.colorbar()
        pl.title(colnames[i]+'_real')
        print '--- R2_adj of the fit is: '+str(R2_adj) +'---'
        pl.subplot(3,3,2*i-3)
        beta,betaErr,R2_adj = zernikeFit(data[:,0].real,data[:,1].real,data[:,i].imag,max_order=zernike_max_order)
        betaAll.append(beta)
        betaErrAll.append(betaErr)
        R2adjAll.append(R2_adj)
        znk=showZernike(beta=beta)
        pl.colorbar()
        pl.title(colnames[i]+'_imag')
        print '--- R2_adj of the fit is: '+str(R2_adj) +'---'
    betaAll = np.array(betaAll)
    betaErrAll = np.array(betaErrAll)
    R2adjAll = np.array(R2adjAll)
    return '----done !---'


def display_coeff(data=None):
    betaAll,betaErrAll, R2adjAll = measure_stamp_coeff(data = data, zernike_max_order=20)
    ind = np.arange(len(betaAll[0]))
    momname = ('M20','M22.Real','M22.imag','M31.real','M31.imag','M33.real','M33.imag')
    fmtarr = ['bo-','ro-','go-','co-','mo-','yo-','ko-']
    pl.figure(figsize=(17,13))
    for i in range(7):
        pl.subplot(7,1,i+1)
        pl.errorbar(ind,betaAll[i],yerr = betaErrAll[i],fmt=fmtarr[i])
        pl.grid()
        pl.xlim(-1,21)
        if i ==0:
            pl.ylim(-10,65)
        elif i ==1:
            pl.ylim(-5,6)
        elif i ==2:
            pl.ylim(-5,6)
        elif i == 3:
            pl.ylim(-0.1,0.1)
        elif i == 4:
            pl.ylim(-0.1,0.1)
        elif i ==5:
            pl.ylim(-100,100)
        elif i == 6:
            pl.ylim(-100,100)
        pl.xticks(ind,('','','','','','','','','','','','','','','','','','','',''))
        pl.ylabel(momname[i])
    pl.xticks(ind,('0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19'))
    pl.xlabel('Zernike Coefficients')
    return '--- done ! ----'

def standardizeData(tdata,vdata):
    """
    This code standardize the training data and validation data by the training data.
    """
    tmean = tdata.mean(axis=0)
    tstd = tdata.std(axis=0)
    tdataNew = (tdata - tmean)/tstd
    vdataNew = (vdata - tmean)/tstd
    return tdataNew, vdataNew

def remM3xZernike(tdata):
    """
    This code remove the 0th zernike coefficient for the M31, M33 from the training and validation data object. The data has 140 columns, from 0 - 59 are the 2nd moments coefficients. 60, 80, 100, 120 are the 0th coefficients for the 3rd moments. We remove them from the data structure.    
    """
    idx = np.concatenate((np.arange(0,60),np.arange(61,80),np.arange(81,100),np.arange(101,120),np.arange(121,140)))
    datanew = tdata[:,idx]
    return datanew


def subMeanM3x(data=None):
    """
    this code subtract the mean of the 3rd moments from the data. This is to remove the tracking errors.
    """
    datamean = data.mean(axis = 0)
    data[:,4:6] = data[:,4:6] - datamean[4:6]
    return data

def subMeanAll(data=None):
    """
    this subtract the mean of all moments except M20 from the data
    """
    datamean = data.mean(axis = 0)
    data[:,3:] = data[:,3:] - datamean[3:]
    return data

def selectStarFwhm(catname):
    ext = [1,2,3,4]
    fwhm_sex=np.zeros(0)
    mag = np.zeros(0)
    for i in ext:
        cat=pf.getdata(catname,i)
        fwhm_sex=np.concatenate((fwhm_sex,cat.FWHM_IMAGE))
        mag = np.concatenate((mag,cat.MAG_AUTO))
    ok = (mag > -15)*(mag<-13)*(fwhm_sex > 0)*(fwhm_sex < 6.)
    md = np.median(fwhm_sex[ok])
    return md
    
def selectStar(mag,fwhm_sex):
    ok = (mag > -15)*(mag<-13)*(fwhm_sex > 0)*(fwhm_sex < 6.)
    md = np.median(fwhm_sex[ok])
    return md
   

