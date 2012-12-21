import os, numpy as np
from decamImgAnalyzer_def import *

expid = np.genfromtxt('/home/jghao/research/ggsvn/des-sv/iq.cat',dtype='S10')
expid.sort()

"""
runid = ['20121210144218_20121207','20121211104015_20121208','20121212123519_20121209','20121217091522_20121215','20121219150303_20121216']

for i in range(len(runid)):
    for j in range(len(expid)):
        print i, j
        for k in range(1,10):
            catname = 'https://desar.cosmology.illinois.edu:7443/DESFiles/desardata/OPS/red/'+runid[i]+'/red/DECam_00'+expid[j]+'/DECam_00'+expid[j]+'_0'+str(k)+'_cat.fits'
            os.system('wget --no-check-certificate --user=jghao --password="jgh70chips" --directory-prefix=/home/jghao/research/data/firstcutcat/ '+catname)
        for k in range(10,63):
            catname = 'https://desar.cosmology.illinois.edu:7443/DESFiles/desardata/OPS/red/'+runid[i]+'/red/DECam_00'+expid[j]+'/DECam_00'+expid[j]+'_'+str(k)+'_cat.fits'
            os.system('wget --no-check-certificate --user=jghao --password="jgh70chips" --directory-prefix=/home/jghao/research/data/firstcutcat/ '+catname)

"""

res = []
for eid in expid:
    print eid
    res.append(whiskerStat_firstcut(eid))

np.savetxt('r50_whk_whkrms_phi_nonwindowed.txt',np.array(res),fmt='%10.5f')
