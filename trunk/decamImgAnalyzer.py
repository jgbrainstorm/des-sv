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

def hexapod_multilinear(beta):
    coeff = p.load(open('/usr/remote/user/sispi/jiangang/des-sv/multiLinearCoeff.p','r'))
    res = np.dot(beta,coeff.T)
    return res



def analyze_hex():
    f = gl.glob('hexapod_position_*.txt')
    f.sort()
    expid=[]
    data=[]
    nexp = len(f)
    if nexp == 0:
        return '-- no image to analyze, exit --'
    for fi in f:
        expid.append(int(fi[-12:-4]))
        data.append(np.genfromtxt(fi))
    expid = np.array(expid)
    xtick = expid.astype('S10')
    data = np.array(data)
    pl.figure(figsize=(16,9))
    pl.subplot(5,1,1)
    pl.plot(expid,data[:,1,0],'bo-',label='Image Analysis')
    pl.plot(expid,data[:,2,0],'ro-',label='BCAM')
    pl.plot(expid,data[:,0,0],'go-',label='Hexapod ')
    pl.ylabel('x-decenter')
    pl.xticks(expid,np.repeat('',nexp))
    pl.legend(loc='upper center',bbox_to_anchor=(0.5,1.68))
    pl.grid()
    pl.subplot(5,1,2)
    pl.plot(expid,data[:,1,1],'bo-',label='Image Analysis')
    pl.plot(expid,data[:,2,1],'ro-',label='BCAM')
    pl.plot(expid,data[:,0,1],'go-',label='Hexapod ')
    pl.ylabel('y-decenter')
    pl.xticks(expid,np.repeat('',nexp))
    pl.grid()
    pl.subplot(5,1,3)
    pl.plot(expid,data[:,1,2],'bo-',label='Image Analysis')
    pl.plot(expid,data[:,0,2],'go-',label='Hexapod ')
    #pl.plot(expid,data[:,2,2],'ro',label='BCAM') #no focus for BCAM
    pl.ylabel('z-defocus')
    pl.xticks(expid,np.repeat('',nexp))
    pl.grid()
    pl.subplot(5,1,4)
    pl.plot(expid,data[:,1,3],'bo-',label='Image Analysis')
    pl.plot(expid,data[:,2,3],'ro-',label='BCAM')
    pl.plot(expid,data[:,0,3],'go-',label='Hexapod ')
    pl.ylabel('x-tilt')
    pl.xticks(expid,np.repeat('',nexp))
    pl.grid()
    pl.subplot(5,1,5)
    pl.plot(expid,data[:,1,4],'bo-',label='Image Analysis')
    pl.plot(expid,data[:,2,4],'ro-',label='BCAM')
    pl.plot(expid,data[:,0,4],'go-',label='Hexapod ')
    pl.xlabel('exposure_id')
    pl.ylabel('y-tilt')
    pl.xticks(expid,xtick,rotation=90)
    pl.grid()
    pl.savefig('hexapod_pos_summary_all.png')
    pl.close()
    pl.figure(figsize=(16,9))
    pl.subplot(5,1,1)
    pl.plot(expid,data[:,1,0],'bo-',label='Image Analysis')
    pl.ylabel('x-decenter')
    pl.xticks(expid,np.repeat('',nexp))
    pl.legend(loc='upper center',bbox_to_anchor=(0.5,1.68))
    pl.grid()
    pl.subplot(5,1,2)
    pl.plot(expid,data[:,1,1],'bo-',label='Image Analysis')
    pl.ylabel('y-decenter')
    pl.xticks(expid,np.repeat('',nexp))
    pl.grid()
    pl.subplot(5,1,3)
    pl.plot(expid,data[:,1,2],'bo-',label='Image Analysis')
    pl.ylabel('z-defocus')
    pl.xticks(expid,np.repeat('',nexp))
    pl.grid()
    pl.subplot(5,1,4)
    pl.plot(expid,data[:,1,3],'bo-',label='Image Analysis')
    pl.ylabel('x-tilt')
    pl.xticks(expid,np.repeat('',nexp))
    pl.grid()
    pl.subplot(5,1,5)
    pl.plot(expid,data[:,1,4],'bo-',label='Image Analysis')
    pl.xlabel('exposure_id')
    pl.ylabel('y-tilt')
    pl.xticks(expid,xtick,rotation=90)
    pl.grid()
    pl.savefig('hexapod_pos_summary_image.png')
    pl.close()
    return '---done!----'


def analyze_r50_whisker():
    """
    make a summary plot for the R50 and whisker length for each exposure
    """
    f = gl.glob('fwhm_whisker_*.p')
    ff = gl.glob('*reduced.fits')
    f.sort()
    ff.sort()
    expid=[]
    r50=[]
    r50err=[]
    whkerr=[]
    whk=[]
    e1 = []
    e1err = []
    e2 = []
    e2err=[]
    flter = []
    nexp = len(f)
    if nexp == 0:
        return '-- no image to analyze, exit --'
    for fi in f:
        expid.append(int(fi[-10:-2]))
        t = p.load(open(fi,'r'))
        r50.append(np.median(t[12]))
        r50err.append(robust_mean_std(np.array(t[12]))[1])
        whk.append(robust_mean_std(np.array(t[6]))[0])
        whkerr.append(robust_mean_std(np.array(t[6]))[1])
        e1.append(robust_mean_std(np.array(t[13]))[0])
        e1err.append(robust_mean_std(np.array(t[13]))[1])
        e2.append(robust_mean_std(np.array(t[14]))[0])
        e2err.append(robust_mean_std(np.array(t[14]))[1])
        flter.append(pf.getheader('DECam_'+fi[-10:-2]+'_reduced.fits')['filter'][0])
    flter = np.array(flter)
    unqfltr = np.unique(flter)
    expid = np.array(expid)
    xtick = expid.astype('S10')
    r50 = np.array(r50)
    r50err = np.array(r50err)
    whk = np.array(whk)
    whkerr = np.array(whkerr)
    e1 = np.array(e1)
    e1err = np.array(e1err)
    e2 = np.array(e2)
    e2err = np.array(e2err)
    xidx = np.arange(len(expid))
    pl.figure(figsize=(16,16))
    fmtarray = ['go','ro','bo','ko','co','mo']
    pl.subplot(4,1,1)
    for k in range(len(unqfltr)):
        ok = flter == unqfltr[k]
        pl.errorbar(xidx[ok],whk[ok],whkerr[ok],fmt=fmtarray[k],label=unqfltr[k])
    pl.hlines(0.2,-1,len(expid),color='green')
    pl.legend(loc='best')
    pl.grid()
    pl.ylabel('whisker length (weighted momts.)')
    pl.xticks(xidx,np.repeat('',len(expid)))
    pl.ylim(0,0.5)
    pl.subplot(4,1,2)
    for k in range(len(unqfltr)):
        ok = flter == unqfltr[k]
        pl.errorbar(xidx[ok],r50[ok],r50err[ok],fmt=fmtarray[k])
    pl.grid()
    pl.ylabel('R50 (sextractor)')
    pl.hlines(0.522,-1,len(expid),color='green')
    pl.xticks(np.arange(len(expid)),np.repeat('',len(expid)))
    pl.subplot(4,1,3)
    for k in range(len(unqfltr)):
        ok = flter == unqfltr[k]
        pl.errorbar(xidx[ok],e1[ok],e1err[ok],fmt=fmtarray[k])
    pl.hlines(0.,-1,len(expid),color='green')
    pl.grid()
    pl.ylabel('e1 (weighted momts.')
    pl.ylim(-0.3,0.3)
    pl.xticks(np.arange(len(expid)),np.repeat('',len(expid)))
    pl.subplot(4,1,4)
    for k in range(len(unqfltr)):
        ok = flter == unqfltr[k]
        pl.errorbar(xidx[ok],e2[ok],e2err[ok],fmt=fmtarray[k])
    pl.hlines(0.,-1,len(expid),color='green')
    pl.grid()
    pl.ylabel('e2 (weighted momts.')
    pl.ylim(-0.3,0.3)
    pl.xticks(np.arange(len(expid)),xtick,rotation=90)
    pl.savefig('exposure_iq_summary.png')
    pl.close()
    for k in range(len(unqfltr)):
        ok = flter == unqfltr[k]
        pl.figure(figsize=(14,14))
        pl.subplot(2,2,1)
        pl.hist(r50[ok],bins=10)
        pl.grid()
        pl.xlabel('R50 (sextractor)')
        pl.subplot(2,2,2)
        pl.hist(whk[ok],bins=10)
        pl.xlabel('Whisker (weighted momts.)')
        pl.grid()
        pl.subplot(2,2,3)
        pl.hist(e1[ok],bins=10)
        pl.xlabel('e1 (weighted momts.)')
        pl.grid()
        pl.subplot(2,2,4)
        pl.hist(e2[ok],bins=10)
        pl.xlabel('e2 (weighted momts.)')
        pl.grid()
        pl.figtext(0.45,0.95,'filter: '+unqfltr[k],fontsize=20,color='blue')
        pl.savefig('iq_distribution_'+unqfltr[k]+'_.png')
        pl.close()
    return '---done!----'


def hexapodPosition(beta,betaErr,weighted=True):
    """
    the CRAY position to the hexapod position parameters. There is a 15 deg rotation between the two coordinate. However, this is accounted in the sispi. So, the hexapod position in the header is acutally rotated to its designed coordiante, which relate to the CRAY coordinate by the last page of des-docdb #6551
    """
    x,y,z,thetax,thetay = CRAYposLinearModel(beta,betaErr,weighted)
    #xh = x                
    #yh = -y
    yh = x       # change on 11/30/12 after comparing with BCAM and hexapod data on 11/22/12
    xh = -y   
    zh = -z
    thetaxh = -thetay
    thetayh = -thetax
    return np.array([xh,yh,zh,thetaxh,thetayh])

def CRAYposLinearModel(b=None,bErr=None,weighted=True):
    """
    here, the convention for b is: M20 (0-19), M22real(20 - 39),M22imag (40-59) 
    """
    M22realTrefoil2 = b[29] 
    M22imagTrefoil1 = b[48] 
    M22realTrefoil2Err = bErr[29] 
    M22imagTrefoil1Err = bErr[48] 
    if weighted == False:
        M22TrefoilXshift = 0.5*(M22realTrefoil2+M22imagTrefoil1) # for x decenter
    else:
        M22TrefoilXshift = (M22realTrefoil2*M22imagTrefoil1Err**2+M22imagTrefoil1*M22realTrefoil2Err**2)/(M22imagTrefoil1Err**2+M22realTrefoil2Err**2) # for x decenter
    M22realTrefoil1 = b[26] 
    M22imagTrefoil2 = b[49] 
    M22realTrefoil1Err = bErr[26] 
    M22imagTrefoil2Err = bErr[49] 
    if weighted == False:
        M22TrefoilYshift = 0.5*(M22realTrefoil1 - M22imagTrefoil2) # for y decenter
    else:
        M22TrefoilYshift = (M22realTrefoil1*M22imagTrefoil2Err**2 - M22imagTrefoil2*M22realTrefoil1Err**2)/(M22imagTrefoil2Err**2+M22realTrefoil1Err**2) # for y decenter
    M20defocus = b[4]    # for defocus
    M22realComa2 = b[28] 
    M22imagComa1 = b[47]
    M22realComa2Err = bErr[28] 
    M22imagComa1Err = bErr[47]
    if weighted == False:
        M22ComaXtilt = 0.5*(M22realComa2+M22imagComa1) # for x-tilt
    else:
        M22ComaXtilt = (M22realComa2*M22imagComa1Err**2+M22imagComa1*M22realComa2Err**2)/(M22imagComa1Err**2+M22realComa2Err**2) # for x-tilt
    M22realComa1 = b[27] 
    M22imagComa2 = b[48]
    M22realComa1Err = bErr[27] 
    M22imagComa2Err = bErr[48]
    if weighted == False:
        M22ComaYtilt = 0.5*(M22realComa1 - M22imagComa2) # for y-tilt
    else:
        M22ComaYtilt = (M22realComa1*M22imagComa2Err**2 - M22imagComa2*M22realComa1Err**2)/(M22imagComa2Err**2+M22realComa1Err**2) # for y-tilt
    x = -3.0063 * M22TrefoilXshift -0.0053
    y = -2.9318 * M22TrefoilYshift - 0.0005
    z = 0.4046 * M20defocus - 0.0705
    xtilt = 1075.8892 * M22ComaXtilt - 0.4876
    ytilt = -1064.6332 * M22ComaYtilt - 0.1234
    return np.array([x*1000,y*1000,z*1000,xtilt,ytilt])


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
        Mcc = cat.X2WIN_IMAGE
        Mrr = cat.Y2WIN_IMAGE
        Mrc = cat.XYWIN_IMAGE
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
            Mrr = Mrr[ok]
            Mcc = Mcc[ok]
            Mrc = Mrc[ok]
            x=x[ok]
            y=y[ok]
            Nobj = len(Mrr)
            M20=np.zeros(Nobj)
            M22=np.zeros(Nobj).astype(complex)
            for k in range(Nobj):
                M20[k] = Mrr[k] + Mcc[k]
                M22[k] = np.complex(Mcc[k] - Mrr[k],2*Mrc[k])
            whisker_sex = np.sqrt(np.abs(M22))
            fwhm_sex = fwhm_sex[ok]
            rad = rad[ok]
            M20 = np.median(M20)
            M22 = np.median(M22)
            stamp = getStamp(data=img,xcoord=x,ycoord=y,Npix=25)
            stamplist = stamplist+stamp
            bkglist = bkglist+list(bkg)
            moms = measure_stamp_moments(stamp,bkg,2.)
            data.append([xccd,yccd]+ list(moms))
            dataSex.append([xccd,yccd,M20,M22])
            fwhmSex = np.concatenate((fwhmSex,fwhm_sex))
            whiskerSex = np.concatenate((whiskerSex,whisker_sex))
            r50Sex = np.concatenate((r50Sex,rad))
        else:
            momsfake=[0.,(0.+0.j),(0.+0.j),(0.+0.j)]
            data.append([xccd,yccd]+list(momsfake))
    if nstarall <=200:
        return '-- This image is bad, skipped ---'
    data = np.array(data)
    data = averageN30(data) # use the average of N31 and N29 to replace N30
    dataSex = np.array(dataSex)
    display_2nd_moments(data)
    pl.savefig('moments_whisker_'+expid+'.png')
    pl.close()
    fwh,whk,r50,e1,e2 = fwhm_whisker_des_plot(stampImgList=stamplist,bkgList=bkglist,whkSex=whiskerSex*0.27,fwhmSex=fwhmSex*0.27,r50Sex=r50Sex*0.27,sigma=2.,dimmfwhm=dimmfwhm)
    pl.savefig('fwhm_whisker_'+expid+'.png')
    pl.close()
    p.dump(fwh+whk+r50+[e1,e2],open('fwhm_whisker_data_'+expid+'.p','w')) # save the fwhm and whisker data.
    #---the hexapod adjustment using M20,M22---
    beta=[]
    betaErr=[]
    betaM20 = zernikeFit(data[:,0].real,data[:,1].real,data[:,2].real,max_order=20)
    beta.append(betaM20[0])
    betaErr.append(betaM20[1])
    betaM22real = zernikeFit(data[:,0].real,data[:,1].real,data[:,3].real,max_order=20)
    beta.append(betaM22real[0])
    betaErr.append(betaM22real[1])
    betaM22imag = zernikeFit(data[:,0].real,data[:,1].real,data[:,3].imag,max_order=20)
    beta.append(betaM22imag[0])
    betaErr.append(betaM22imag[1])
    betaforplot = beta
    betaErrforplot = betaErr
    beta=np.array(beta)
    betaErr = np.array(betaErr)
    beta=beta.flatten()
    betaErr = betaErr.flatten()
    posCRAY = CRAYposLinearModel(beta,betaErr,weighted=False)
    hexHao = hexapodPosition(beta,betaErr,weighted=False)
    #hexHao = hexapod_multilinear(beta)
    dispM202Coeff(betaAll = betaforplot, betaErrAll = betaErrforplot,hexinfo=hexHao)
    pl.savefig('zernike_coeff_'+expid+'.png')
    pl.close()
    #---save files---
    hexposhdr = np.array(hexposhdr.split(',')).astype(float)[0:5]
    hexBCAM = np.array([bcamDX,bcamDY,-999,bcamAX,bcamAY])
    np.savetxt('hexapod_position_'+expid+'.txt',[hexposhdr,hexHao,hexBCAM,posCRAY],fmt='%10.5f')
    return '----finished one image ----'
    

if __name__ == "__main__":
    from decamImgAnalyzer import *
    import sys,time,glob
    startTime=time.time()
    if len(sys.argv) == 1:
        print 'syntax: '
        print 'desImgAnalysis expid'
        print 'or'
        print 'desImgAnalysis all'
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
    tt=analyze_hex()
    tt = analyze_r50_whisker()
    print '---elapsed time: ' + str(elapseTime)

