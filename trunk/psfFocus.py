import sys
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
    import cPickle as p
    import sklearn.neighbors as nb
    from sklearn.svm import SVR
    from decamRay import *
    from ImgQuality import *
except ImportError:
    print "Error: missing one of the libraries (numpy, pyfits, scipy, matplotlib)"
    sys.exit()


#scale=0.27/4.
#npix=160

def regulate_img(img=None):
    """
    trim and normalize the image
    """
    #img = (img - img.min())/img.sum()
    img = img*(img>=0)
    return img/img.sum()

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
    return beta,betaErr, R2adj



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
    qk = pl.quiverkey(qvr, -150,-240,np.max(np.sqrt(u**2+v**2)),str(round(np.max(np.sqrt(u**2+v**2)),3))+' pix^2',coordinates='data',color='blue')
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

def KNeighborRegression(trainingObs,trainingParam,Obs,n_neighbors):
    """
    using the k nearest neighbor regression to get the parameter
    trainingObs: the zernike coefficients data matrix. each row is a new observation
    trainingParam: the hexapod configration and seeing, each row is a new obs.
    Obs: a given set of measured coefficients
    """
    #knn = nb.KNeighborsRegressor(algorithm='ball_tree',n_neighbors=n_neighbors,weights = 'distance')
    knn = nb.KNeighborsRegressor(algorithm='ball_tree',n_neighbors=n_neighbors)
    knn.fit(trainingObs,trainingParam)
    return knn.predict(Obs)

def remM3xZernike(tdata):
    """
    This code remove the 0th zernike coefficient for the M31, M33 from the training and validation data object. The data has 140 columns, from 0 - 59 are the 2nd moments coefficients. 60, 80, 100, 120 are the 0th coefficients for the 3rd moments. We remove them from the data structure.    
    """
    idx = np.concatenate((np.arange(0,60),np.arange(61,80),np.arange(81,100),np.arange(101,120),np.arange(121,140)))
    datanew = tdata[:,idx]
    return datanew


  
def get_hexapod_pos(data=None):
    """
    output are: xshift, yshift, z, theta_x, theta_y
    unit for shift is mm
    unit for tilt angle is arcsec
    """
    betaAll, betaErrAll, R2adjAll=measure_stamp_coeff(data)
    Tfile='/home/jghao/research/decamFocus/psf_withseeing/finerGrid_coeff_matrix/zernike_coeff_finerGrid_all.cp'
    b=p.load(open(Tfile))
    nobs = len(b)
    tdata=b[:,8:].copy()
    ttpara=b[:,0:5].copy()
    tpara=b[:,0:5].copy()
    tpara[:,3] = ttpara[:,3]*np.cos(np.deg2rad(ttpara[:,4]))
    tpara[:,4] = ttpara[:,3]*np.sin(np.deg2rad(ttpara[:,4]))
    vdata = betaAll.flatten()
    # remove the 0th coeff for 3rd moments
    #tdata = remM3xZernike(tdata)
    #vdata = remM3xZernike(vdata)
    tdata,vdata = standardizeData(tdata,vdata)
    vparaReg=KNeighborRegression(tdata,tpara,vdata,15)
    xshift = np.round(vparaReg[0][0]*1000,2)
    yshift = np.round(vparaReg[0][1]*1000,2)
    zshift = np.round(vparaReg[0][2]*1000,2)
    xtilt = np.round(vparaReg[0][3],2)
    ytilt = np.round(vparaReg[0][4],2)
    print '--- x,y,z are in microns, tilt is in arcsec ---'
    return xshift,yshift,zshift,xtilt,ytilt

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
    ok = (mag > -15)*(mag<-13)*(fwhm_sex > 0)*(fwhm_sex < 20.)
    md = np.median(fwhm_sex[ok])
    return md
    
def selectStar(mag,fwhm_sex):
    ok = (mag > -15)*(mag<-13)*(fwhm_sex > 0)*(fwhm_sex < 20.)
    md = np.median(fwhm_sex[ok])
    return md
   

def moments_uncertainty(ccd=None,mag=15,Nstar=1,seeing=[0.9,0.,0.],theta=0,phi=0,x=0,y=0,z=0,npix=npix):
    momtsW = []
    momtsA = []
    momtsWTrue = []
    momtsATrue = []
    for i in range(50):        
        img,bkg,psf = des_image(exptime=100,mag=mag, Nstar=Nstar,ccd=ccd,seeing=seeing,npix=npix,zenith=0,filter='r', theta=theta, phi=phi,corrector='corrector',x=x,y=y,z=z)
        momtsW.append(complexMoments(img-bkg,sigma = 2.))
        momtsWTrue.append(complexMoments(psf,sigma = 2.))
        momtsA.append(AcomplexMoments(img - bkg))
        momtsATrue.append(AcomplexMoments(psf))
    return np.array(momtsW),np.array(momtsWTrue),np.array(momtsA),np.array(momtsATrue)


def moments_uncertainty_new(ccd=None,Nstar=1,mag=None,seeing=[0.9,0.,0.],theta=0,phi=0,x=0,y=0,z=0,npix=npix):
    momtsW = []
    momtsA = []
    momtsWTrue = []
    momtsATrue = []
    img,bkg,psf,magg = des_image(exptime=100,mag=[16.5,18], Nstar=Nstar,ccd=ccd,seeing=seeing,npix=npix,zenith=0,filter='r', theta=theta, phi=phi,corrector='corrector',x=x,y=y,z=z)
    for i in range(len(img)):
        momtsW.append(complexMoments(img[i]-bkg,sigma = 2.))
        momtsWTrue.append(complexMoments(psf[i],sigma = 2.))
        momtsA.append(AcomplexMoments(img[i] - bkg,sigma = 2.))
        momtsATrue.append(AcomplexMoments(psf[i],sigma = 2.))
    return np.array(momtsW),np.array(momtsWTrue),np.array(momtsA),np.array(momtsATrue)



    

if __name__ == "__main__":
    from psfFocus import *
    pl.ion()
    #dr = '/home/jghao/research/data/des_optics_psf/dc6b_image/decam--18--38-r-10/'
    dr = '/home/jghao/research/data/des_optics_psf/dc6b_image/goodseeing/decam--28--49-r-1/'
    starfile='/home/jghao/research/data/des_optics_psf/dc6b_image/goodseeing/catfile/decam_-27.72186_-48.60000-objects.fit'
    extension = np.arange(1,63)
    data=[]
    ccd = ['S1', 'S2', 'S3', 'N1', 'N2', 'N3', 'S8', 'S9', 'S14', 'S15', 'S20', 'S25', 'N8', 'N9', 'N14', 'N15', 'N20', 'N25', 'S10', 'S11', 'S12', 'S13', 'S18', 'S19', 'S16', 'S17', 'S21', 'S22', 'S23', 'S24', 'S26','S27', 'S28', 'S29', 'S30', 'S31', 'N4', 'N5', 'N6', 'N7', 'S4', 'S5','S6', 'S7', 'N10', 'N11', 'N12', 'N13', 'N18', 'N19', 'N16', 'N17', 'N21', 'N22', 'N23', 'N24', 'N26', 'N27', 'N28', 'N29', 'N30', 'N31']
    for ext in extension:
        stamp=[]
        print ext
        xccd = eval(ccd[ext-1])[1]
        yccd = eval(ccd[ext-1])[2]
        if ext < 10:
            ext = '0'+str(ext)
        else:
            ext = str(ext)
        imgname= 'decam--28--49-r-1_'+ext+'.fits.fz'
        bkgname = 'decam--28--49-r-1_'+ext+'_bkg.fits.fz'
        img = pf.getdata(dr+imgname) - pf.getdata(dr+bkgname)
        xc = pf.getdata(starfile,3*(int(ext)-1)+1).xccd
        yc = pf.getdata(starfile,3*(int(ext)-1)+1).yccd
        rmag = pf.getdata(starfile,3*(int(ext)-1)+1).mag_3
        ok = (rmag > 16.5)*(rmag < 19.5)*(xc > 100)*(xc<1900)*(yc>200)*(yc<3900)
        xc=xc[ok]
        yc = yc[ok]    
        stamp=stamp+getStamp(data=img,xcoord=xc,ycoord=yc,Npix =40)
        moms = measure_stamp_moments(stamp,sigma=1.)
        data.append([xccd,yccd]+ list(moms))
    data = np.array(data)
    display_moments(data)
    display_coeff(data)
    get_hexapod_pos(data)

    #-----verify the noise effect ---
    momtsW,momtsWTrue,momtsA,momtsATrue = moments_uncertainty_new(ccd=S4,mag=[16.5,18],Nstar=50)

    pl.figure(figsize=(10,8))
    pl.subplot(2,1,1)
    pl.boxplot(momtsW.real)
    pl.xticks(np.arange(1,5),['M20.real','M22.real','M31.real','M33.real'])
    pl.grid()
    pl.title('mag: 15, CCD: S4')
    pl.subplot(2,1,2)
    pl.boxplot(momtsW.imag)
    pl.xticks(np.arange(1,5),['M20.imag','M22.imag','M31.imag','M33.imag'])
    pl.grid()
    
    pl.figure(figsize=(10,8))
    pl.subplot(2,1,1)
    pl.boxplot(momtsWTrue.real)
    pl.xticks(np.arange(1,5),['M20.real','M22.real','M31.real','M33.real'])
    pl.grid()
    pl.title('mag: 15, CCD: S4, Truth')
    pl.subplot(2,1,2)
    pl.boxplot(momtsWTrue.imag)
    pl.xticks(np.arange(1,5),['M20.imag','M22.imag','M31.imag','M33.imag'])
    pl.grid()

    #----distribution ----
    pl.figure(figsize=(15,9))
    pl.subplot(2,3,1)
    pl.hist(momtsW[:,0].real,bins=10,normed=False)
    pl.title(str(round(np.mean(momtsW[:,0].real),6)) + r'$\pm$'+str(round(np.std(momtsW[:,0].real)/np.sqrt(momtsW.shape[0]),6)))
    pl.xlabel('Wmoment M20')
    pl.subplot(2,3,2)
    pl.hist(momtsW[:,1].real,bins=10,normed=False)
    pl.title(str(round(np.mean(momtsW[:,1].real),6)) + r'$\pm$'+str(round(np.std(momtsW[:,1].real)/np.sqrt(momtsW.shape[0]),6)))
    pl.xlabel('Wmoment M22.real')
    pl.subplot(2,3,3)
    pl.hist(momtsW[:,1].imag,bins=10,normed=False)
    pl.title(str(round(np.mean(momtsW[:,1].imag),6)) + r'$\pm$'+str(round(np.std(momtsW[:,1].imag)/np.sqrt(momtsW.shape[0]),6)))
    pl.xlabel('Wmoment M22.imag')
    pl.subplot(2,3,4)
    pl.hist(momtsA[:,0].real,bins=10,normed=False)
    pl.title(str(round(np.mean(momtsA[:,0].real),6)) + r'$\pm$'+str(round(np.std(momtsA[:,0].real)/np.sqrt(momtsA.shape[0]),6)))
    pl.xlabel('Amoment M20')
    pl.subplot(2,3,5)
    pl.hist(momtsA[:,1].real,bins=10,normed=False)
    pl.title(str(round(np.mean(momtsA[:,1].real),6)) + r'$\pm$'+str(round(np.std(momtsA[:,1].real)/np.sqrt(momtsA.shape[0]),6)))
    pl.xlabel('Amoment M22.real')
    pl.subplot(2,3,6)
    pl.hist(momtsA[:,1].imag,bins=10,normed=False)
    pl.title(str(round(np.mean(momtsA[:,1].imag),6)) + r'$\pm$'+str(round(np.std(momtsA[:,1].imag)/np.sqrt(momtsA.shape[0]),6)))
    pl.xlabel('Amoment M22.imag')
  
    #--true value --
    pl.figure(figsize=(15,9))
    pl.subplot(2,3,1)
    pl.hist(momtsW[:,0].real,bins=10,normed=False)
    pl.title(str(round(np.median(momtsWTrue[:,0].real),6)) + r'$\pm$'+str(round(np.std(momtsWTrue[:,0].real)/np.sqrt(momtsWTrue.shape[0]),6)))
    pl.xlabel('Wmoment M20')
    pl.subplot(2,3,2)
    pl.hist(momtsWTrue[:,1].real,bins=10,normed=False)
    pl.title(str(round(np.median(momtsWTrue[:,1].real),6)) + r'$\pm$'+str(round(np.std(momtsWTrue[:,1].real)/np.sqrt(momtsWTrue.shape[0]),6)))
    pl.xlabel('Wmoment M22.real')
    pl.subplot(2,3,3)
    pl.hist(momtsWTrue[:,1].imag,bins=10,normed=False)
    pl.title(str(round(np.median(momtsWTrue[:,1].imag),6)) + r'$\pm$'+str(round(np.std(momtsWTrue[:,1].imag)/np.sqrt(momtsWTrue.shape[0]),6)))
    pl.xlabel('Wmoment M22.imag')
    pl.subplot(2,3,4)
    pl.hist(momtsATrue[:,0].real,bins=10,normed=False)
    pl.title(str(round(np.median(momtsATrue[:,0].real),6)) + r'$\pm$'+str(round(np.std(momtsATrue[:,0].real)/np.sqrt(momtsATrue.shape[0]),6)))
    pl.xlabel('Amoment M20')
    pl.subplot(2,3,5)
    pl.hist(momtsATrue[:,1].real,bins=10,normed=False)
    pl.title(str(round(np.median(momtsATrue[:,1].real),6)) + r'$\pm$'+str(round(np.std(momtsATrue[:,1].real)/np.sqrt(momtsATrue.shape[0]),6)))
    pl.xlabel('Amoment M22.real')
    pl.subplot(2,3,6)
    pl.hist(momtsATrue[:,1].imag,bins=10,normed=False)
    pl.title(str(round(np.median(momtsATrue[:,1].imag),6)) + r'$\pm$'+str(round(np.std(momtsATrue[:,1].imag)/np.sqrt(momtsATrue.shape[0]),6)))
    pl.xlabel('Amoment M22.imag')
  
def dispM202Coeff(betaAll=None,betaErrAll=None):
    ind = np.arange(len(betaAll[0]))
    momname = ('M20','M22.Real','M22.imag')
    fmtarr = ['bo-','ro-','go-']
    if betaErrAll == None:
        betaErrAll = np.zeros(len(ind))
    pl.figure(figsize=(17,7))
    for i in range(3):
        pl.subplot(4,1,i+1)
        pl.plot(ind[1:],betaAll[i][1:],yerr = betaErrAll[i],fmt=fmtarr[i])
        pl.grid()
        pl.xlim(-1,len(betaAll[i])+1)
        pl.ylim(min(betaAll[i][1:])-0.5,max(betaAll[i][1:])+0.5)
        pl.xticks(ind,('','','','','','','','','','','','','','','','','','','',''))
        pl.ylabel(momname[i])
    pl.xticks(ind,('Piston','Tip','Tilt','Defocus','Astignism','Astignism','Coma','Coma','Trefoil','Trefoil','Spherical','12','13','14','15','16','17','18','19','20'),rotation=90)
    pl.xlabel('Zernike Coefficients')
    return '---done!---'
