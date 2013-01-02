#--------------------------------------------------------------------------
# This codes display the 2nd moments whisker based on the moments measured 
# from the Quick Reduce
# J. Hao @ FNAL, 12/4/2012
#--------------------------------------------------------------------------
import pylab as pl, numpy as np

def subMeanAll(data=None):
    """
    this subtract the mean of all moments except M20 from the data
    """
    datamean = data.mean(axis = 0)
    data[:,3:] = data[:,3:] - datamean[3:]
    return data

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


def whisker4QReduce(X2WIN_IMAGE=None,Y2WIN_IMAGE=None,XYWIN_IMAGE=None):
    """
    This function make the whisker plot based on the sextractor export from QuickReduce
    J.Hao, 12/4/2012
    """
    xpos = np.genfromtxt('xpos_ypos_fp.txt').T[0]
    ypos = np.genfromtxt('xpos_ypos_fp.txt').T[1]
    temp = np.zeros(62).astype('complex')
    temp.real = X2WIN_IMAGE - Y2WIN_IMAGE
    temp.imag = 2*XYWIN_IMAGE
    data=np.array([xpos, ypos,X2WIN_IMAGE + Y2WIN_IMAGE,temp]).T
    data = averageN30new(data)
    data = subMeanAll(data)
    pl.figure(figsize=(11,5.5))
    pl.subplot(1,2,1)
    phi22 = 0.5*np.arctan2(data[:,3].imag,data[:,3].real)
    x = data[:,0].real
    y = data[:,1].real
    u = np.abs(data[:,3])*np.cos(phi22)
    v = np.abs(data[:,3])*np.sin(phi22)
    qvr = pl.quiver(x,y,u,v,width = 0.004, color='r',pivot='middle',headwidth=0.,headlength=0.,headaxislength=0.,scale_units='width')
    qk = pl.quiverkey(qvr, -150,-240,0.3,str(0.3)+' pix^2',coordinates='data',color='blue')
    pl.plot(x,y,'b,')
    pl.xlim(-250,250)
    pl.ylim(-250,250)
    pl.grid(color='g')
    pl.xlabel('Camera WEST [mm]')
    pl.ylabel('Camera NORTH [mm]')
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
    pl.scatter(-230,-210,s=0.01*plotScale,c='b',alpha=0.5)
    pl.text(-200,-215,'-'+str(0.01)+' pix')
    pl.scatter(-230,-230,s=0.01*plotScale,c='r',alpha=0.5)
    pl.text(-200,-235,str(0.01)+' pix')
    pl.plot(x,y,'y,')
    pl.grid(color='g')
    pl.xlim(-250,250)
    pl.ylim(-250,250)
    pl.xlabel('Camera WEST [mm]')
    pl.ylabel('Camera NORTH [mm]')
    pl.title('median '+r'$\sqrt{M20}$: '+str(round(0.27*m20sqr_med,3))+' [arcsec]')
    return '---done!--'
