from psfFocus import *
#from decamRay import *

#imghdu = pf.open('/home/jghao/research/data/des_realimage/test2.fits')
#cathdu = pf.open('/home/jghao/research/data/des_realimage/test2.cat')

imghdu = pf.open('/home/jghao/research/data/des_realimage/fromGary/psf134136.fits')
#cathdu = pf.open('/home/jghao/research/data/des_realimage/fromGary/Gary_cat/psf134136.cat')
cathdu = pf.open('/home/jghao/research/data/des_realimage/fromGary/psf134136_star_catalog.fits')

imghdu = pf.open('/home/jghao/research/data/des_realimage/goodseeing/DECam_00137625_reduced.fits.fz')
cathdu = pf.open('/home/jghao/research/data/des_realimage/goodseeing/starcat.fits')

data=[]
stamplist=[]
for i in range(1,63):
    print i
    img = imghdu[i].data
    cat = cathdu[i*2].data
    x = cat.XWIN_IMAGE
    y = cat.YWIN_IMAGE
    rad = cat.FLUX_RADIUS
    mag = cat.MAG_AUTO
    flag = cat.FLAGS
    bkg = cat.BACKGROUND
    #ok = (rad > 1.7)*(rad < 2.25)*(mag > -14)*(mag < -12)
    ok = (rad > 4)*(rad <6)*(mag > -14)*(mag < -12)
    x = x[ok]
    y = y[ok]
    bkg = bkg[ok]
    stamp = getStamp(data=img,xcoord=x,ycoord=y,Npix=40)
    stamplist = stamplist+stamp
    xccd = eval(imghdu[i].header['extname'])[1]
    yccd = eval(imghdu[i].header['extname'])[2]
    moms = measure_stamp_moments(stamp,bkg,1.5)
    data.append([xccd,yccd]+ list(moms))

data = np.array(data)
display_moments(data)
display_coeff(data)
xshift,yshift,zshift,thetax,thetay = get_hexapod_pos(data)
phi = np.rad2deg(np.arctan2(thetay,thetax))
theta = np.sqrt(thetax**2+thetay**2)

moments_display(Nstar=1,npix = npix,x=xshift*0.001,y=yshift*0.001,z=zshift*0.001,theta=theta,phi=phi)
fwhm_whisker_plot(stamplist)

#---on new star catalog. 

imghdu = pf.open('/home/jghao/research/data/des_realimage/goodseeing/DECam_00137625_reduced.fits.fz')
cathdu = pf.open('/home/jghao/research/data/des_realimage/goodseeing/starcat.fits')

data=[]
datasex=[]
stamplist=[]
fwhm = np.array([])
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
    fwhmsex = cat.FWHM_IMAGE
    #ok = (rad > 1.7)*(rad < 2.25)*(mag > -14)*(mag < -12)
    ok = selectStar(rad,mag)
    x = x[ok]
    y = y[ok]
    bkg = bkg[ok]
    Mrr = np.median(Mrr[ok])
    Mcc = np.median(Mcc[ok])
    Mrc = np.median(Mrc[ok])
    fwhmsex = fwhmsex[ok]
    M20 = Mrr + Mcc
    M22 = np.complex(Mcc - Mrr,2*Mrc)
    stamp = getStamp(data=img,xcoord=x,ycoord=y,Npix=25)
    stamplist = stamplist+stamp
    xccd = eval(imghdu[i].header['extname'])[1]
    yccd = eval(imghdu[i].header['extname'])[2]
    moms = measure_stamp_moments(stamp,bkg,4)
    data.append([xccd,yccd]+ list(moms))
    datasex.append([xccd,yccd,M20,M22])
    fwhm = np.concatenate((fwhm,fwhmsex))

data = np.array(data)
datasex = np.array(datasex)
display_2nd_moments(datasex)
display_moments(data)
display_coeff(data)
xshift,yshift,zshift,thetax,thetay = get_hexapod_pos(data)
phi = np.rad2deg(np.arctan2(thetay,thetax))
theta = np.sqrt(thetax**2+thetay**2)

moments_display(Nstar=1,npix = npix,x=xshift*0.001,y=yshift*0.001,z=zshift*0.001,theta=theta,phi=phi)
fwhm_whisker_plot(stamplist)
