#---------------------------------
# This is the code for IQ-R4,5.
# Initiated by J. Hao, 7/13/2012
# Comments sent to: jghao@fnal.gov 
#---------------------------------
import numpy as np


def wr(x,y,xcen,ycen,sigma):
    """
    Returns a gaussian weight function with the given parameters
    """
    res=np.exp(-((x-xcen)**2+(y-ycen)**2)/(2.*sigma**2))/(2.*np.pi*sigma**2) 
    return res


def complexMoments(data=None,sigma=None):
    """
    This one calcualte the 3 2nd moments and 4 thrid moments with the Gaussian weights.
    col : x direction
    row : y direction
    sigma is the stand deviation of the measurement kernel in pixel
    input: 2d stamp image, Gaussian weight kernel size (in pixels)
    output: 
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
    wrmat = wr(ROW,COL,rowmean,colmean,sigma)
    IWmat = data*wrmat
    wrcol = wrmat.sum(axis=0)
    wrrow = wrmat.sum(axis=1)
    wrcolsum = np.sum(Icol*wrcol)
    wrrowsum = np.sum(Irow*wrrow)
    rowmean = np.sum(rowgrid*Irow*wrrow)/wrrowsum
    colmean = np.sum(colgrid*Icol*wrcol)/wrcolsum
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


def whisker(data=None, sigma = None):
    """
    This code calculate the wisker length defined as sqrt(e1^2+e2^2)
    input: 
         data: 2d stamp image
         sigma: std of the Gaussian weight Kernel in pixel
    """
    M20, M22, M31, M33 = complexMoments(data=data,sigma=sigma)
    whisker_length = np.sqrt(np.abs(M22))
    return whisker_length


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
    
