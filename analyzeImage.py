#---this codes analyze the images ----
from decamImgAnalyzer import *
from getDatabase import *

def getMeanStd(datafile):
    t = p.load(open(datafile,'r'))
    resMean = []
    resStd=[]
    for i in range(len(t)):
        resMean.append(t[i][~np.isnan(t[i])].mean())
        resStd.append(t[i][~np.isnan(t[i])].std())
    return np.array(resMean)

def whisker2ellip(whisker,fwhm):
    ellip = whisker**2/fwhm**2*2.77259
    return ellip

def ellip2whisker(ellip,fwhm):
    whk = np.sqrt(ellip/2.77259*fwhm**2)
    return whk

def getfwhm_ellip(expid):
    nexp = len(expid)
    fwhmQR = np.zeros(nexp)
    ellipQR = np.zeros(nexp)
    fwhmIH = np.zeros(nexp)
    ellipIH = np.zeros(nexp)
    for i in range(nexp):
        resQR = getQRbyExpid(expid=expid[i])
        resIH = getIHbyExpid(expid=expid[i])
        fwhmQR[i]=resQR[1]
        ellipQR[i]=resQR[2]
        fwhmIH[i]=resIH[1]
        ellipIH[i]=resIH[2]
    return fwhmQR,ellipQR,fwhmIH,ellipIH



rdir = 'rband/'
idir = 'iband/'
zdir = 'zband/'


hexfile_r = gl.glob(rdir+'*.txt')
hexfile_r.sort()
hexfile_i = gl.glob(idir+'*.txt')
hexfile_i.sort()
hexfile_z = gl.glob(zdir+'*.txt')
hexfile_z.sort()

Nfile = len(hexfile_r)
expidr = []
expidi = []
expidz = []
for i in range(Nfile):
    expidr.append(hexfile_r[i][-12:-4])
    expidi.append(hexfile_i[i][-12:-4])
    expidz.append(hexfile_z[i][-12:-4])
expidr = np.array(expidr).astype('i')
expidi = np.array(expidi).astype('i')
expidz = np.array(expidz).astype('i')

fwf_r = gl.glob(rdir+'*.p')
fwf_r.sort()
fwf_i = gl.glob(idir+'*.p')
fwf_i.sort()
fwf_z = gl.glob(zdir+'*.p')
fwf_z.sort()

fw_r = []
fw_i = []
fw_z = []

for i in range(Nfile):
    fw_r.append(getMeanStd(fwf_r[i]))
    fw_i.append(getMeanStd(fwf_i[i]))
    fw_z.append(getMeanStd(fwf_z[i]))

fw_r = np.array(fw_r)
fw_i = np.array(fw_i)
fw_z = np.array(fw_z)

#---get the value from QR and IH
fwhmQR_r,ellipQR_r,fwhmIH_r,ellipIH_r=getfwhm_ellip(expidr)
fwhmQR_i,ellipQR_i,fwhmIH_i,ellipIH_i=getfwhm_ellip(expidi)
fwhmQR_z,ellipQR_z,fwhmIH_z,ellipIH_z=getfwhm_ellip(expidz)
whkQR_r = ellip2whisker(ellipQR_r,fwhmQR_r)
whkQR_i = ellip2whisker(ellipQR_i,fwhmQR_i)
whkQR_z = ellip2whisker(ellipQR_z,fwhmQR_z)

whkIH_r = ellip2whisker(ellipIH_r,fwhmIH_r)
whkIH_i = ellip2whisker(ellipIH_i,fwhmIH_i)
whkIH_z = ellip2whisker(ellipIH_z,fwhmIH_z)


#------fwhm plots -----
expidr = expidr.astype('i')
pl.plot(fw_r[:,0],fw_r[:,1],'ro',label='Adapt')
pl.plot(fw_r[:,0],fw_r[:,2],'go',label='Moffat')
pl.plot(fw_r[:,0],fw_r[:,3],'bo',label='Gauss')
pl.plot(fw_r[:,0],fw_r[:,4],'co',label='Sech2')
pl.plot(fw_r[:,0],fw_r[:,5],'mo',label='Sextractor')
pl.plot(fw_r[:,0],fwhmQR_r,'r*',label='QR')
pl.plot(fw_r[:,0],fwhmIH_r,'b*',label='IH')
pl.xlabel('fwhm_weighted [arcsec]')
pl.ylabel('other fwhm [arcsec]')
pl.plot([0.7,1.3],[0.7,1.3],'k-')
pl.grid()
pl.ylim(0.7,1.3)
pl.xlim(0.7,1.3)
pl.legend(loc='best')
pl.title('r-band')
pl.savefig('fwhm_rband_comp.png')
pl.close()


pl.plot(fw_i[:,0],fw_i[:,1],'ro',label='Adapt')
pl.plot(fw_i[:,0],fw_i[:,2],'go',label='Moffat')
pl.plot(fw_i[:,0],fw_i[:,3],'bo',label='Gauss')
pl.plot(fw_i[:,0],fw_i[:,4],'co',label='Sech2')
pl.plot(fw_i[:,0],fw_i[:,5],'mo',label='Sextractor')
pl.plot(fw_r[:,0],fwhmQR_i,'r*',label='QR')
pl.plot(fw_r[:,0],fwhmIH_i,'b*',label='IH')
pl.xlabel('fwhm_weighted [arcsec]')
pl.ylabel('other fwhm [arcsec]')
pl.plot([0.7,1.3],[0.7,1.3],'k-')
pl.grid()
pl.ylim(0.7,1.3)
pl.xlim(0.7,1.3)
pl.legend(loc='best')
pl.title('i-band')
pl.savefig('fwhm_iband_comp.png')
pl.close()

pl.plot(fw_z[:,0],fw_z[:,1],'ro',label='Adapt')
pl.plot(fw_z[:,0],fw_z[:,2],'go',label='Moffat')
pl.plot(fw_z[:,0],fw_z[:,3],'bo',label='Gauss')
pl.plot(fw_z[:,0],fw_z[:,4],'co',label='Sech2')
pl.plot(fw_z[:,0],fw_z[:,5],'mo',label='Sextractor')
pl.plot(fw_r[:,0],fwhmQR_z,'r*',label='QR')
pl.plot(fw_r[:,0],fwhmIH_z,'b*',label='IH')
pl.xlabel('fwhm_weighted [arcsec]')
pl.ylabel('other fwhm [arcsec]')
pl.plot([0.7,1.3],[0.7,1.3],'k-')
pl.grid()
pl.ylim(0.7,1.3)
pl.xlim(0.7,1.3)
pl.legend(loc='best')
pl.title('z-band')
pl.savefig('fwhm_zband_comp.png')
pl.close()


#---whikser length ------
pl.plot(fw_r[:,6],fw_r[:,7],'ro',label='Adapt')
pl.plot(fw_r[:,6],fw_r[:,8],'go',label='Sextractor')
pl.plot(fw_r[:,6],whkQR_r],'r*',label='QR')
pl.plot(fw_r[:,6],whkIH_r],'b*',label='IH')

pl.xlabel('whisker_weighted [arcsec]')
pl.ylabel('other whisker [arcsec]')
pl.plot([0.1,0.4],[0.1,0.4],'k-')
pl.grid()
pl.ylim(0.05,0.4)
pl.xlim(0.05,0.4)
pl.legend(loc='best')
pl.title('r-band')
pl.savefig('whisker_rband_comp.png')
pl.close()


pl.plot(fw_i[:,6],fw_i[:,7],'ro',label='Adapt')
pl.plot(fw_i[:,6],fw_i[:,8],'go',label='Sextractor')
pl.plot(fw_i[:,6],whkQR_i],'r*',label='QR')
pl.plot(fw_i[:,6],whkIH_i],'b*',label='IH')
pl.xlabel('whisker_weighted [arcsec]')
pl.ylabel('other whisker [arcsec]')
pl.plot([0.1,0.4],[0.1,0.4],'k-')
pl.grid()
pl.ylim(0.05,0.4)
pl.xlim(0.05,0.4)
pl.legend(loc='best')
pl.title('i-band')
pl.savefig('whisker_iband_comp.png')
pl.close()


pl.plot(fw_z[:,6],fw_z[:,7],'ro',label='Adapt')
pl.plot(fw_z[:,6],fw_z[:,8],'go',label='Sextractor')
pl.plot(fw_z[:,6],whkQR_z],'r*',label='QR')
pl.plot(fw_z[:,6],whkIH_z],'b*',label='IH')

pl.xlabel('whisker_weighted [arcsec]')
pl.ylabel('other whisker [arcsec]')
pl.plot([0.1,0.4],[0.1,0.4],'k-')
pl.grid()
pl.ylim(0.05,0.4)
pl.xlim(0.05,0.4)
pl.legend(loc='best')
pl.title('z-band')
pl.savefig('whisker_zband_comp.png')
pl.close()


#-----check R-19----
# mean psf whisker lenth for each exposure must be less than 0.2" in r, i, z
pl.figure(figsize=(15,15))
pl.subplot(3,1,1)
pl.plot(expidr,fw_r[:,6],'bo',label='weighted')
pl.plot(expidr,fw_r[:,7],'ro',label='Adapt')
pl.plot(expidr,fw_r[:,8],'go',label='Sextractor')
pl.plot(expidr,whkQR_r,'r*',label='QR')
pl.plot(expidr,whkIH_r,'b*',label='IH')
pl.legend(loc='best')
pl.xticks(expidr,expidr.astype('S10'),rotation=40)
pl.grid()
pl.ylabel('mean whisker length')
pl.hlines(0.2,expidr[0]-1,expidr[-1]+1,color='k')
pl.ylim(0.,0.5)
pl.title('r-band')

pl.subplot(3,1,2)
pl.plot(expidi,fw_i[:,6],'bo',label='weighted')
pl.plot(expidi,fw_i[:,7],'ro',label='Adapt')
pl.plot(expidi,fw_i[:,8],'go',label='Sextractor')
pl.plot(expidi,whkQR_i,'r*',label='QR')
pl.plot(expidi,whkIH_i,'b*',label='IH')
pl.legend(loc='best')
pl.xticks(expidi,expidi.astype('S10'),rotation=40)
pl.grid()
pl.ylabel('mean whisker length')
pl.hlines(0.2,expidr[0]-1,expidr[-1]+1,color='k')
pl.ylim(0.,0.5)
pl.title('i-band')

pl.subplot(3,1,3)
pl.plot(expidz,fw_z[:,6],'bo',label='weighted')
pl.plot(expidz,fw_z[:,7],'ro',label='Adapt')
pl.plot(expidz,fw_z[:,8],'go',label='Sextractor')
pl.plot(expidz,whkQR_z,'r*',label='QR')
pl.plot(expidz,whkIH_z,'b*',label='IH')
pl.legend(loc='best')
pl.xticks(expidz,expidz.astype('S10'),rotation=40)
pl.grid()
pl.ylabel('mean whisker length')
pl.hlines(0.2,expidr[0]-1,expidr[-1]+1,color='k')
pl.ylim(0.,0.5)
pl.title('z-band')
pl.savefig('IQ-R5_whker.png')
pl.close()



#-----check R-18 b----
# Moreover, for 95% of the survey area, there should be at least one exposure in each of these bands for which the mean PSF FWHM is 0.9” or smaller

pl.figure(figsize=(15,15))
pl.subplot(3,1,1)
pl.plot(expidr,fw_r[:,1],'ro',label='Adapt')
pl.plot(expidr,fw_r[:,2],'go',label='Moffat')
pl.plot(expidr,fw_r[:,3],'bo',label='Gauss')
pl.plot(expidr,fw_r[:,4],'co',label='Sech2')
pl.plot(expidr,fw_r[:,5],'mo',label='Sextractor')
pl.plot(expidr,fwhmQR_r,'r*',label='QR')
pl.plot(expidr,fwhmIH_r,'b*',label='IH')
pl.ylabel('mean fwhm [arcsec]')
pl.legend(loc='upper right',bbox_to_anchor=(1.1,1.3))
pl.xticks(expidr,expidr.astype('S10'),rotation=40)
pl.grid()
pl.hlines(0.9,expidr[0]-1,expidr[-1]+1,color='k')
pl.ylim(0.5,1.4)
pl.title('r-band')

pl.subplot(3,1,2)
pl.plot(expidi,fw_i[:,1],'ro',label='Adapt')
pl.plot(expidi,fw_i[:,2],'go',label='Moffat')
pl.plot(expidi,fw_i[:,3],'bo',label='Gauss')
pl.plot(expidi,fw_i[:,4],'co',label='Sech2')
pl.plot(expidi,fw_i[:,5],'mo',label='Sextractor')
pl.plot(expidi,fwhmQR_i,'r*',label='QR')
pl.plot(expidi,fwhmIH_i,'b*',label='IH')
pl.ylabel('mean fwhm [arcsec]')
pl.xticks(expidi,expidi.astype('S10'),rotation=40)
pl.grid()
pl.hlines(0.9,expidi[0]-1,expidi[-1]+1,color='k')
pl.ylim(0.5,1.4)
pl.title('i-band')

pl.subplot(3,1,3)
pl.plot(expidz,fw_z[:,1],'ro',label='Adapt')
pl.plot(expidz,fw_z[:,2],'go',label='Moffat')
pl.plot(expidz,fw_z[:,3],'bo',label='Gauss')
pl.plot(expidz,fw_z[:,4],'co',label='Sech2')
pl.plot(expidz,fw_z[:,5],'mo',label='Sextractor')
pl.plot(expidz,fwhmQR_z,'r*',label='QR')
pl.plot(expidz,fwhmIH_z,'b*',label='IH')
pl.ylabel('mean fwhm [arcsec]')
pl.xticks(expidz,expidz.astype('S10'),rotation=40)
pl.grid()
pl.hlines(0.9,expidz[0]-1,expidz[-1]+1,color='k')
pl.ylim(0.5,1.4)
pl.title('z-band')
pl.savefig('IQ-R4_fwhm_b.png')
pl.close()
