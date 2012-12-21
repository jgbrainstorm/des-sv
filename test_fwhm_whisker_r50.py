# this code test the various measurement for the fwhm, whisker, r50

from decamImgAnalyzer_def import *
import scipy.stats as st


def add_imageNoise(img):
    """
    add poisson noise to images
    """
    if not np.all(img >= 0):
        print 'make sure the image pixel values are positive definite'
        sys.exit()
    noise = st.poisson.rvs(1.,loc = -1.,scale=1.,size=img.shape)*np.sqrt(img)
    return noise


def moffat_seeing(npix = None, alpha=None,beta=None):
    row,col = np.mgrid[-npix/2:npix/2,-npix/2:npix/2]
    rowc = row.mean()
    colc = col.mean()
    #img = (beta - 1)/(np.pi*alpha**2)/(1+((row**2+col**2)/alpha**2))**beta
    img = 10.*(1+(row**2+col**2)/alpha**2)**(-beta)
    res = img/img.sum()
    return res



def gauss_seeing(npix = None,fwhm=None,e1=None,e2=None,scale=scale):
    """
    generate a seeing PSF of given fwhm and e1 and e2
    fwhm in the unit of arcsec
    """
    fwhm = fwhm/scale
    M20 = 2.*(fwhm/2.35482)**2
    row,col = np.mgrid[-npix/2:npix/2,-npix/2:npix/2]
    rowc = row.mean()
    colc = col.mean()
    Mcc = 0.5*M20*(1+e1)
    Mrc = 0.5*e2*M20
    Mrr = 0.5*M20*(1-e1)
    rho = Mrc/np.sqrt(Mcc*Mrr)
    img = np.exp(-0.5/(1-rho**2)*(row**2/Mrr + col**2/Mcc - 2*rho*row*col/np.sqrt(Mrr*Mcc)))
    res = img/img.sum()
    return res


def des_psf_image(exptime=100,mag=None,seeing=[0.9,0.,0.],alpha=None,beta=None,setbkg=True,seeingType='gauss'):
    
    """
    This code generate a PSF star with seeing and sky background (no optics psf)
    exptime is given in sec
    seeing is give in terms of [fwhm (arcsec),e1,e2]
    """
    gain = 0.21 # convert electrons to ADU
    npix = 40
    zeropoint = 26.794176 # r band, from Nikolay
    objectphoton = exptime*10**(0.4*(zeropoint - mag))
    if setbkg == False:
        skyphoton = 0.
    else:
        skyphoton = 8.460140*exptime #(sky level per pix per sec)
    bkg = skyphoton*gain  # background in ADU
    if seeingType == 'gauss':
        psf = gauss_seeing(25,seeing[0],seeing[1],seeing[2],scale = 0.27)
    if seeingType == 'moffat':
        psf = moffat_seeing(25,alpha,beta)
    img = (psf * objectphoton + skyphoton)*gain
    img = img + add_imageNoise(img)
    return img,bkg,psf

def fwhm_whisker_des_plot_sim(stampImgList=None,bkgList=None,whkSex=None,fwhmSex=None,r50Sex=None,whkInput=None,fwhmInput=None,r50Input=None,sigma=1.1/scale,seeingType='gauss'):
    whk,fwhm,r50 = get_fwhm_whisker_list(stampImgList,bkgList,sigma=sigma)
    whk=list(whk.T)
    fwh=list(fwhm.T)
    r50=list(r50.T)
    fwh.append(fwhmSex)
    whk.append(whkSex)
    r50.append(r50Sex)
    fwh.append(fwhmInput)
    whk.append(whkInput)
    r50.append(r50Input)
    pl.figure(figsize=(15,15))
    pl.subplot(3,1,1)
    pl.boxplot(whk)
    pl.hlines(0.2,0,5,linestyle='solid',color='g')
    pl.ylim(np.median(whk[2])-0.3,np.median(whk[2])+0.6)
    pl.grid()
    pl.xticks(np.arange(1,5),['whisker_Wmoments','whisker_Amoments','whisker_sx','whisker_input'])
    if seeingType == 'gauss':
        pl.title('PSF generated assuming Bivariate Gaussian profile')
    if seeingType == 'moffat':
        pl.title('PSF generated assuming Moffat profile')
    pl.subplot(3,1,2)
    pl.boxplot(fwh)
    pl.ylim(0,np.median(fwh[5])+2)
    pl.grid()
    pl.hlines(0.9,0,8,linestyle='solid',color='g')
    pl.xticks(np.arange(1,8),['fwhm_weighted', 'fwhm_Amoments','fwhm_moffat', 'fwhm_gauss','fwhm_sech2','fwhm_sx','fwhm_Input'])
    pl.subplot(3,1,3)
    pl.boxplot(r50)
    pl.ylim(0,np.median(r50[1])+0.5)
    pl.grid()
    pl.hlines(0.5,0,6,linestyle='solid',color='g')
    pl.xticks(np.arange(1,6),['R50_Sech2', 'R50_Moffat','R50_Gaussian', 'R50_sx','R50_input'])
    return fwh,whk,r50


if __name__== "__main__":
    from test_fwhm_whisker_r50 import *
    np.random.seed(1)
    fwhm = np.random.randint(700,1500,500)/1000.
    np.random.seed(2)
    e1 = np.random.randint(-120,120,500)/1000.
    np.random.seed(3)
    e2 = np.random.randint(-120,120,500)/1000.
    np.random.seed(4)
    alpha = np.random.randint(2500,4500,500)/1000.
    np.random.seed(5)
    beta = np.random.randint(2500,4500,500)/1000.
    """
    #---gaussian profile ----
    img = []
    bkg = []
    hduList = pf.HDUList()
    for i in range(len(fwhm)):
        sim = des_psf_image(exptime=100,mag=16.5,seeing=[fwhm[i],e1[i],e2[i]],setbkg=True)
        img.append(sim[0])
        bkg.append(sim[1])
        hdu = pf.PrimaryHDU(sim[0])
        hduList.append(hdu)
    hduList.writeto('sim_500_gaussian.fits')
    hduList.close()


    #---moffat profile ---
    imgm = []
    bkgm = []
    hduList = pf.HDUList()
    for i in range(len(alpha)):
        sim = des_psf_image(exptime=100,mag=16.5,alpha=alpha[i],beta=beta[i],setbkg=True,seeingType='moffat')
        imgm.append(sim[0])
        bkgm.append(sim[1])
        hdu = pf.PrimaryHDU(sim[0])
        hduList.append(hdu)
    hduList.writeto('sim_500_moffat.fits')
    hduList.close()
    """
    
    #---analyze gaussian profile ------
    
    imghduG = pf.open('sim_500_gaussian.fits') 
    cathduG = pf.open('sim_500_gaussian_star_catalog.fits')
    stampG = []
    bkgG = []
    r50SexG = []
    fwhmSexG = []
    whkSexG = []
    for i in range(len(imghduG)):
        stampG.append(imghduG[i].data)
        bkgG.append(cathduG[i+1].data.BACKGROUND)
        fwhmSexG.append(cathduG[i+1].data.FWHM_IMAGE*0.27)
        r50SexG.append(cathduG[i+1].data.FLUX_RADIUS*0.27)
        Mcc = cathduG[i+1].data.X2WIN_IMAGE
        Mrr = cathduG[i+1].data.Y2WIN_IMAGE
        Mrc = cathduG[i+1].data.XYWIN_IMAGE
        whkSexG.append(((Mcc- Mrr)**2+(2*Mrc)**2)**(0.25)*0.27)

    imghduG.close()
    cathduG.close()
 
    whkSexG = np.array(whkSexG)
    r50SexG = np.array(r50SexG)
    fwhmSexG = np.array(fwhmSexG)

    fwhmInputG = fwhm
    whkInputG = (fwhm/1.665)*(e1**2+e2**2)**0.25
    r50InputG = fwhm/2.
    fakemag = np.arange(10)
    fakerad = np.zeros(10)+2.

    fwG,whkG, r50G =fwhm_whisker_des_plot_sim(stampImgList=stampG,bkgList=bkgG,whkSex=whkSexG,fwhmSex=fwhmSexG,r50Sex=r50SexG,whkInput=whkInputG,fwhmInput=fwhmInputG,r50Input=r50InputG,sigma=2.,seeingType='gauss')
    pl.savefig('sim_test_gauss.png')
    pl.close()

    pl.figure(figsize=(18,6))
    pl.subplot(1,3,1)
    pl.plot(whkG[3],whkG[0],'b.',label='weighted',alpha=0.5)
    pl.plot(whkG[3],whkG[1],'r.',label='ataptive',alpha=0.5)
    pl.plot(whkG[3],whkG[2],'k.',label='sextractor',alpha=0.8)
    pl.plot([0.,0.5],[0,0.5],'r-',lw=2)
    pl.legend(loc='best')
    pl.xlabel('Input Whisker')
    pl.ylabel('Output Whisker')

    pl.subplot(1,3,2)
    pl.plot(fwG[6],fwG[0],'b.',label='weighted',alpha=0.5)
    pl.plot(fwG[6],fwG[1],'r.',label='ataptive',alpha=0.5)
    pl.plot(fwG[6],fwG[2],'k.',label='moffat',alpha=0.7)
    pl.plot(fwG[6],fwG[3],'g.',label='gauss',alpha=0.8)
    pl.plot(fwG[6],fwG[4],'c.',label='sech2',alpha=0.8)
    pl.plot(fwG[6],fwG[5],'y.',label='sextractor',alpha=0.8)
    pl.plot([0.5,1.6],[0.5,1.6],'r-',lw=2)
    pl.legend(loc='best')
    pl.xlabel('Input FWHM')
    pl.ylabel('Output FWHM')

    pl.subplot(1,3,3)
    pl.plot(r50G[4],r50G[0],'c.',label='sech2',alpha=0.5)
    pl.plot(r50G[4],r50G[1],'r.',label='moffat',alpha=0.5)
    pl.plot(r50G[4],r50G[2],'k.',label='gauss',alpha=0.7)
    pl.plot(r50G[4],r50G[3],'g.',label='sextractor',alpha=0.8)
    pl.plot(r50G[4],fwG[0]/2.,'b.',label='weighted',alpha=0.8)
    pl.plot([0.1,0.8],[0.1,.8],'r-',lw=2)
    pl.legend(loc='best')
    pl.xlabel('Input R50')
    pl.ylabel('Output R50')

    pl.figtext(0.35,0.95,'PSF from Bivariate Gaussian',fontsize=18)
    pl.savefig('sim_compare_gauss.png')
    pl.close()

    #---analyze moffat profile ------
    
    imghduM = pf.open('sim_500_moffat.fits') 
    cathduM = pf.open('sim_500_moffat_star_catalog.fits')
    stampM = []
    bkgM = []
    r50SexM = []
    fwhmSexM = []
    whkSexM = []
    for i in range(len(imghduM)):
        stampM.append(imghduM[i].data)
        bkgM.append(cathduM[i+1].data.BACKGROUND)
        fwhmSexM.append(cathduM[i+1].data.FWHM_IMAGE*0.27)
        r50SexM.append(cathduM[i+1].data.FLUX_RADIUS*0.27)
        Mcc = cathduM[i+1].data.X2WIN_IMAGE
        Mrr = cathduM[i+1].data.Y2WIN_IMAGE
        Mrc = cathduM[i+1].data.XYWIN_IMAGE
        whkSexM.append(((Mcc- Mrr)**2+(2*Mrc)**2)**(0.25)*0.27)

    imghduM.close()
    cathduM.close()
    
    whkSexM = np.array(whkSexM)
    r50SexM = np.array(r50SexM)
    fwhmSexM = np.array(fwhmSexM)

    fwhmInputM = 2.*alpha*np.sqrt(2.**(1./beta)-1)*0.27
    whkInputM = 0.
    r50InputM = alpha*np.sqrt(2.**(1./(beta-1))-1)*0.27
    fakemag = np.arange(10)
    fakerad = np.zeros(10)+2.

    fwM, whkM,r50M =fwhm_whisker_des_plot_sim(stampImgList=stampM,bkgList=bkgM,whkSex=whkSexM,fwhmSex=fwhmSexM,r50Sex=r50SexM,whkInput=whkInputM,fwhmInput=fwhmInputM,r50Input=r50InputM,sigma=2.,seeingType='moffat')

    pl.savefig('sim_test_moffat.png')
    pl.close()

    
    pl.figure(figsize=(18,6.))
    pl.subplot(1,3,1)
    pl.plot(np.zeros(len(whkM[0])),whkM[0],'b.',label='weighted',alpha=0.5)
    pl.plot(np.zeros(len(whkM[0])),whkM[1],'r.',label='ataptive',alpha=0.5)
    pl.plot(np.zeros(len(whkM[0])),whkM[2],'k.',label='sextractor',alpha=0.8)
    pl.plot([0.,0],[0,0.5],'r-',lw=2)
    pl.legend(loc='best')
    pl.xlabel('Input Whisker')
    pl.ylabel('Output Whisker')

    pl.subplot(1,3,2)
    pl.plot(fwM[6],fwM[0],'b.',label='weighted',alpha=0.5)
    pl.plot(fwM[6],fwM[1],'r.',label='ataptive',alpha=0.5)
    pl.plot(fwM[6],fwM[2],'k.',label='moffat',alpha=0.7)
    pl.plot(fwM[6],fwM[3],'g.',label='gauss',alpha=0.8)
    pl.plot(fwM[6],fwM[4],'c.',label='sech2',alpha=0.8)
    pl.plot(fwM[6],fwM[5],'y.',label='sextractor',alpha=0.8)
    pl.plot([0.5,1.6],[0.5,1.6],'r-',lw=2)
    pl.legend(loc='best')
    pl.xlabel('Input FWHM')
    pl.ylabel('Output FWHM')

    pl.subplot(1,3,3)
    pl.plot(r50M[4],r50M[0],'c.',label='sech2',alpha=0.5)
    pl.plot(r50M[4],r50M[1],'r.',label='moffat',alpha=0.5)
    pl.plot(r50M[4],r50M[2],'k.',label='gauss',alpha=0.7)
    pl.plot(r50M[4],r50M[3],'g.',label='sextractor',alpha=0.8)
    pl.plot(r50M[4],fwM[0]/2.,'b.',label='weighted',alpha=0.8)
    pl.plot([0.2,1.2],[0.2,1.2],'r-',lw=2)
    pl.legend(loc='best')
    pl.xlabel('Input R50')
    pl.ylabel('Output R50')

    pl.figtext(0.35,0.95,'PSF from Moffat profile',fontsize=18)
    pl.savefig('sim_compare_moffat.png')
    pl.close()
