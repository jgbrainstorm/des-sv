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

scale = 0.27
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
        print i
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

def fwhm_whisker(image=None, sigma = None):
    """
    This code calculate the fwhm and wisker length defined as (M22.real^2 + M22.imag^2)^{1/4}
    input: 
         data: 2d stamp image
         sigma: std of the Gaussian weight Kernel in pixel
    """
    M20, M22, M31, M33 = complexMoments(image=image,sigma=sigma)
    fwhm = np.sqrt(M20)*2.35482*scale
    whisker_length = np.sqrt(np.abs(M22))*scale
    return fwhm,whisker_length


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
    dir = '/home/jghao/research/data/des_optics_psf/dc6b_image/'
    imgname = 'decam-34-0-r-0_01.fits'
    bkgname = 'decam-34-0-r-0_01_bkg.fits'
    hdr = pf.getheader(dir+imgname)
    data = pf.getdata(dir+imgname) - pf.getdata(dir+bkgname)
    xc,yc=findbstr(data=data, hdr=hdr)
    stamp=getStamp(data=data,xcoord=xc,ycoord=yc,Npix =20)
    Nstamp = len(stamp)
    for i in range(Nstamp):
        print fwhm_whisker(image=stamp[i], sigma = 4.)
        
