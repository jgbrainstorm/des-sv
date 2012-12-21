import pylab as pl, numpy as np,cPickle as p


def compare_qr_hao_scatter(expid=None):
    qr = np.genfromtxt('/home/jghao/research/desSV/psfChallenge/'+str(expid)+'.dat',usecols=(1,2,3,4))
    hao = p.load(open('/home/jghao/research/desSV/psfChallenge/compare_fwhm_whisker_data_00'+str(expid)+'.p','r'))
    #fwhm: wmoment,amoment,moffat,gauss1d,sech2,sex
    #whk: wmoments,amoments,sex
    #r50: sech2, moffat,gauss1d,sex
    pl.figure(figsize=(15,5))
    pl.subplot(1,3,1)
    pl.plot(qr[:,0],hao[0][:,5],'b.')
    pl.plot([0,2],[0,2],'r-')
    pl.xlim(0.7,1.6)
    pl.ylim(0.7,1.6)
    pl.grid()
    pl.xlabel('fwhm: quick reduce')
    pl.ylabel('fwhm: sextractor Hao')
    pl.subplot(1,3,2)
    pl.plot(np.sqrt(qr[:,1]**2+qr[:,2]**2),hao[1][:,2]**2*2.77259/hao[0][:,5]**2,'b.')
    pl.plot([0,2],[0,2],'r-')
    pl.xlim(0.,0.1)
    pl.ylim(0.,0.1)
    pl.grid()
    pl.xlabel('ellipticity: quick reduce')
    pl.ylabel('ellipticity: sextractor Hao')
    pl.subplot(1,3,3)
    pl.plot(qr[:,3]*0.27,hao[2][:,3],'b.')
    pl.plot([0.4,0.8],[0.4,.8],'r-')
    pl.xlim(0.4,0.8)
    pl.ylim(0.4,0.8)
    pl.grid()
    pl.xlabel('R50: quick reduce')
    pl.ylabel('R50: sextractor Hao')
    pl.figtext(0.4,0.95,'exposure_id:'+str(expid))
    pl.savefig('sextractor_compare_qr_hao_'+str(expid)+'.png')


def compare_qr_hao_boxplot(expid=None):
    qr = np.genfromtxt('/home/jghao/research/desSV/psfChallenge/'+str(expid)+'.dat',usecols=(1,2,3,4))
    hao = p.load(open('/home/jghao/research/desSV/psfChallenge/compare_fwhm_whisker_data_00'+str(expid)+'.p','r'))
    #fwhm: wmoment,amoment,moffat,gauss1d,sech2,sex
    #whk: wmoments,amoments,sex
    #r50: sech2, moffat,gauss1d,sex
    pl.figure(figsize=(15,5))
    pl.subplot(1,3,1)
    pl.boxplot([qr[:,0],hao[0][:,5]])
    pl.ylim(0,2)
    pl.grid()
    pl.xticks(np.arange(1,3),('fwhm_QR','fwhm_Sx_Hao'))
    pl.subplot(1,3,2)
    pl.boxplot([np.sqrt(qr[:,1]**2+qr[:,2]**2),hao[1][:,2]**2*2.77259/hao[0][:,5]**2])
    pl.ylim(0.,0.1)          
    pl.grid()
    pl.xticks(np.arange(1,3),('ellip_QR','ellip_Sx_Hao'))          
    pl.subplot(1,3,3)
    pl.boxplot([qr[:,3]*0.27,hao[2][:,3]])
    pl.grid()
    pl.ylim(0.4,0.8)
    pl.xticks(np.arange(1,3),('R50_QR','R50_Sx_Hao'))          
    pl.figtext(0.4,0.95,'exposure_id:'+str(expid))
    pl.savefig('sextractor_compare_qr_hao_boxplot_'+str(expid)+'.png')


def median_fwhm_r50():
    """
    This code calculate the scaling between R50 and fwhm
    """
    haof = gl.glob('/home/jghao/research/desSV/psfChallenge/compare*.p')
    haof.sort()
    haof = haof[0:35]
    nf = len(haof)
    fwhmsex = np.zeros(nf)
    r50sex = np.zeros(nf)
    for i in range(nf):
        print i
        fwhmsex[i] = np.median(p.load(open(haof[i],'r'))[0][:,5])
        r50sex[i] = np.median(p.load(open(haof[i],'r'))[2][:,3])
    pl.plot(fwhmsex,r50sex,'b.')
    res = linfit(fwhmsex,r50sex,np.sqrt(r50sex))
    pl.plot(np.array([0.8,1.6]),np.array([0.8,1.6])*0.48+0.09,'r-')
    pl.figtext(0.2,0.8,'slope: 0.48')
    pl.figtext(0.2,0.77,'Intercept: 0.09')
    pl.xlim(0.7,1.6)
    pl.ylim(0.4,0.9)
    pl.xlabel('fwhm_sx')
    pl.ylabel('R50_sx')
    pl.grid()
