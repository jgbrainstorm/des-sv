#! /usr/bin/env python

import os, sys,numpy as np
from decamImgAnalyzer_def import *

#expid = np.genfromtxt('/home/jghao/research/ggsvn/des-sv/downiq.cat',dtype='S10')
#expid = np.genfromtxt('/home/s1/jghao/ggsvn/des-sv/downiq.cat',dtype='S10')
expid = np.genfromtxt('/home/s1/jghao/ggsvn/des-sv/iq_163831.cat',dtype='S10')
expid.sort()

i = int(sys.argv[1])
#runid = ['20121210144218_20121207','20121211104015_20121208','20121212123519_20121209','20121217091522_20121215','20121219150303_20121216']
#runid=['20121219090047_20121217','20121219150555_20121218','20120720091430_20121219','20120720174936_20121220']
#runid=['20121219145915_20121217','20121223151438_20121220','20121222135459_20121221','20121223151000_20121222']
"""
runid=['20121227110356_20121223', '20130105122302_20121226','20130104154724_20121227','20121231093746_20121228','20130101152443_20121229','20130105113400_20130103','20130105171408_20130104']


for j in range(len(expid)):
    print i, j
    for k in range(1,10):
        catname = 'https://desar.cosmology.illinois.edu:7443/DESFiles/desardata/OPS/red/'+runid[i]+'/red/DECam_00'+expid[j]+'/DECam_00'+expid[j]+'_0'+str(k)+'_cat.fits'
        os.system('wget --no-check-certificate --user=jghao --password="jgh70chips" --directory-prefix=/data/des08.b/data/jiangang/firstcut/ '+catname)
    for k in range(10,63):
        catname = 'https://desar.cosmology.illinois.edu:7443/DESFiles/desardata/OPS/red/'+runid[i]+'/red/DECam_00'+expid[j]+'/DECam_00'+expid[j]+'_'+str(k)+'_cat.fits'
        os.system('wget --no-check-certificate --user=jghao --password="jgh70chips" --directory-prefix=/data/des08.b/data/jiangang/firstcut/ '+catname)




"""


reshao = []
resmike = []

if i ==0:
    for eid in expid:
        print eid
        reshao.append(whiskerStat_firstcut(eid))
    np.savetxt('firstcut_stat_hao_163831.txt',np.array(reshao),fmt='%10.5f')

if i == 1:
    for eid in expid:
        print eid
        resmike.append(whiskerStat_firstcut_mike(eid))
        np.savetxt('firstcut_stat_mike_163831.txt',np.array(resmike),fmt='%10.5f')

