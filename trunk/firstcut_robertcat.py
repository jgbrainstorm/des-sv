# this code calcuate the quantities from the firstcut catlogs from Robert's query. 
from func_def import *
import glob as gl
import numpy as np, pylab as pl

def check_mag_radius(filename,ext):
    b=np.genfromtxt(filename,delimiter=',')
    idx = b[:,3] == ext
    rad = b[idx,8]
    mag = b[idx,6]
    pl.plot(mag,rad,'b.')
    pl.xlabel('mag_auto')
    pl.ylabel('flux_radius')
    pl.ylim(0,10)

def firstcutStarRobert(b):
    rad = b[:,8]
    mag = b[:,6]
    flag = b[:,9]
    ok = (mag>=10.5)*(mag<=12)*(flag ==0)*(rad<5.)
    radmedian = np.median(rad[ok])
    idx = (mag>=10.5)*(mag<=13)*(flag ==0)*(abs(rad-radmedian)<=0.2)
    return b[idx,:]

def correctMoments(Mcc=None, Mrr=None, Mrc=None,measureSigma=None):
    '''
    This function correct the weighted moments to remove the weight effect
    All based on the Gaussian assumption
    '''
    for i in range(len(Mcc)):
        Cm = np.matrix([[Mcc[i],Mrc[i]],[Mrc[i],Mrr[i]]])
        Cw = np.matrix([[measureSigma[i]**2,0.],[0.,measureSigma[i]**2]])
        Cimg = (Cm.I - Cw.I).I
        Mcc[i] = Cimg[0,0]
        Mrr[i] = Cimg[1,1]
        Mrc[i] = Cimg[0,1]
    return Mcc, Mrr, Mrc
    
f =gl.glob('/home/jghao/research/data/firstcutcat/iqcat_robert_feb_2013/*.csv') 
f.sort()

n=len(f)
res = []
fltr=[]
for i in range(n):
    data = []
    print i
    b=np.genfromtxt(f[i],delimiter=',')
    bflt=np.genfromtxt(f[i],delimiter=',',dtype='str')
    if b.shape[0]<400:
        continue
    fltr.append(bflt[0,2])
    for ext in np.append(range(1,61),62):
        idx = b[:,3] == ext
        if np.any(idx):
            bb = firstcutStarRobert(b[idx,:])
            fluxrad = robust_mean(bb[:,8])
            #Mcc = robust_mean(bb[:,10])
            #Mrr = robust_mean(bb[:,11])
            #Mrc = robust_mean(bb[:,12])
            #Mcc = robust_mean(bb[:,16]) # windowed
            #Mrr = robust_mean(bb[:,17])
            #Mrc = robust_mean(bb[:,18])    
            # ----using size corrected weighted moments: 3/13/2013----
            measureSigma = bb[:,8]*2./2.35482
            Mcc,Mrr,Mrc = correctMoments(bb[:,16],bb[:,17],bb[:,18],measureSigma)
            Mcc = robust_mean(Mcc)
            Mrr = robust_mean(Mrr)
            Mrc = robust_mean(Mrc)
            #--------------------------------
            fwhmworld = robust_mean(bb[:,19])
            data.append([Mcc,Mrr,Mrc,fluxrad,fwhmworld])
        else:
            continue
    data = np.array(data)
    datamean =np.array([robust_mean(data[:,0]),robust_mean(data[:,1]),robust_mean(data[:,2]),robust_mean(data[:,3]),robust_mean(data[:,4])])
    whk = ((datamean[0]-datamean[1])**2 + (2.*datamean[2])**2)**(0.25)*3600
    phi = np.rad2deg(0.5*np.arctan2(2.*datamean[2],(datamean[0]-datamean[1])))
    r50 = datamean[3]*0.27
    fwhm = datamean[4]*3600
    datasubmean = data - datamean
    whkrms = (robust_mean((datasubmean[:,0] - datasubmean[:,1])**2 + 4.*datasubmean[:,2]**2))**(0.25)*3600
    expid = int(f[i][-10:-4])
    res.append([expid,r50,whk,whkrms,phi,fwhm])

res=np.array(res)
np.savetxt('/home/jghao/research/data/firstcutcat/iqcat_robert_feb_2013/firstcut_iq_hao_20130214_windowed_corrected_Moments_newstar.txt',res,fmt=['%i','%10.5f','%10.5f','%10.5f','%10.5f','%10.5f'])


f =gl.glob('/home/jghao/research/data/firstcutcat/iqcat_robert_feb_2013/*.csv') 
f.sort()

n=len(f)
res = []
fltr=[]
for i in range(n):
    data = []
    print i
    b=np.genfromtxt(f[i],delimiter=',',dtype='str')
    if b.shape[0]<400:
        continue
    fltr.append(b[0,2])
np.savetxt('/home/jghao/research/data/firstcutcat/iqcat_robert_feb_2013/new_iq_filter_newstar.txt',fltr,fmt='%10s')

#----analyze the results------------

# expid, r50, whk, whkrms

#b=np.genfromtxt('/home/jghao/research/data/firstcutcat/iqcat_robert_feb_2013/firstcut_iq_hao_20130214_newstar.txt')
b=np.genfromtxt('/home/jghao/research/data/firstcutcat/iqcat_robert_feb_2013/firstcut_iq_hao_20130214_windowed_Moments_newstar.txt')
#b=np.genfromtxt('/home/jghao/research/data/firstcutcat/iqcat_robert_feb_2013/firstcut_iq_hao_20130214_windowed_corrected_Moments_newstar.txt')
fltr = np.genfromtxt('/home/jghao/research/data/firstcutcat/iqcat_robert_feb_2013/new_iq_filter_newstar.txt',dtype='S10')

ok=((fltr=='r')+(fltr=='i')+(fltr=='z'))*(b[:,1]>0)*(b[:,1]<1.)
bb= b[ok,:]

pl.subplot(1,2,1)
pl.hist(bb[:,1],bins=20)
pl.xlabel('R50')
pl.title('median: '+str(np.median(bb[:,1])))
pl.subplot(1,2,2)
pl.hist(bb[:,2],bins=20)
pl.xlabel('Whisker')
pl.title('median: '+str(np.median(bb[:,2])))
#pl.figtext(0.35,0.95,'r, i, z, Feb.1 - Feb.14, 2013, windowed moments corrected')
#pl.savefig('/home/jghao/research/data/firstcutcat/iqcat_robert_feb_2013/r50_whker_riz_window_correct.png')
pl.figtext(0.35,0.95,'r, i, z, Feb.1 - Feb.14, 2013, windowed moments without correction')
pl.savefig('/home/jghao/research/data/firstcutcat/iqcat_robert_feb_2013/r50_whker_riz_window.png')




#----combine to get the final catalog ----
b=np.genfromtxt('/home/jghao/research/data/firstcutcat/iqcat_robert_feb_2013/firstcut_iq_hao_20130214_newstar.txt')
expid = np.genfromtxt('/home/jghao/research/data/firstcutcat/iqcat_robert_feb_2013/iq_feb.cat')
nrow = expid.shape[0]
data = np.zeros((nrow,6))
for i in range(nrow):
    ok = b[:,0] == expid[i]
    data[i,0] = expid[i]
    if np.any(ok):
        if b[ok,1] > 0:
            data[i,1:] = b[ok,1:]
        else:
            data[i,1:] = np.array([-999,-999,-999,-999,-999])
    else:
        data[i,1:] = np.array([-999,-999,-999,-999,-999])

np.savetxt('/home/jghao/research/data/firstcutcat/iqcat_robert_feb_2013/firstcut_iq_feb_hao.txt',data,fmt=['%i','%10.5f','%10.5f','%10.5f','%10.5f','%10.5f'])

