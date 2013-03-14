#! /usr/bin/env python

import os, sys,numpy as np
from decamImgAnalyzer_def import *

#expid = np.genfromtxt('/home/jghao/research/ggsvn/des-sv/downiq.cat',dtype='S10')
#expid = np.genfromtxt('/home/s1/jghao/ggsvn/des-sv/downiq.cat',dtype='S10')
#expid = np.genfromtxt('/home/s1/jghao/ggsvn/des-sv/iq_163831_new.cat',dtype='S10')
#expid = np.genfromtxt('/home/s1/jghao/ggsvn/des-sv/aos_exposures_1_13_2013.txt',dtype='S10')
#expid = np.genfromtxt('/home/s1/jghao/ggsvn/des-sv/iq_20130131_newentry.cat',dtype='S10')

expid = np.genfromtxt('/home/s1/jghao/ggsvn/des-sv/iq_20130214.cat',dtype='S10')
expid.sort()

i = int(sys.argv[1])



#runid = ['20121210144218_20121207','20121211104015_20121208','20121212123519_20121209','20121217091522_20121215','20121219150303_20121216']
#runid=['20121219090047_20121217','20121219150555_20121218','20120720091430_20121219','20120720174936_20121220']
#runid=['20121219145915_20121217','20121223151438_20121220','20121222135459_20121221','20121223151000_20121222']
#runid=['20121227110356_20121223', '20130105122302_20121226','20130104154724_20121227','20121231093746_20121228','20130101152443_20121229','20130105113400_20130103','20130105171408_20130104']
#runid=['20130105171408_20130104','20130106142401_20130105','20130107104210_20130106','20130108095202_20130107','20130111160812_20130108','20130114122604_20130113','20130115124305_20130114','20130116101510_20130115']
#runid =['20130117092652_20130116','20130118080849_20130117','20130119114402_20130118','20130120104633_20130119','20130122093733_20130121','20130124110605_20130122','20130130115812_20130129','20130131160808_20130130']

runid=['20130201081044_20130131','20130202103930_20130201','20130203103627_20130202','20130204112807_20130203','20130205101526_20130204','20130206100153_20130205','20130207100905_20130206','20130211151905_20130207','20130209073451_20130208','20130221124016_20130220','20130222105744_20130221']

for j in range(len(expid)):
    print i, j
    for k in range(1,10):
        catname = 'https://desar.cosmology.illinois.edu:7443/DESFiles/desardata/OPS/red/'+runid[i]+'/red/DECam_00'+expid[j]+'/DECam_00'+expid[j]+'_0'+str(k)+'_cat.fits'
        os.system('wget --no-check-certificate --user=jghao --password="jgh70chips" --directory-prefix=/data/des08.b/data/jiangang/firstcut/iq_20130214/ '+catname)
        #os.system('wget --no-check-certificate --user=jghao --password="jgh70chips" --directory-prefix=/home/jghao/research/data/firstcutcat/fircutMoments/ '+catname)
    for k in range(10,63):
        catname = 'https://desar.cosmology.illinois.edu:7443/DESFiles/desardata/OPS/red/'+runid[i]+'/red/DECam_00'+expid[j]+'/DECam_00'+expid[j]+'_'+str(k)+'_cat.fits'
        os.system('wget --no-check-certificate --user=jghao --password="jgh70chips" --directory-prefix=/data/des08.b/data/jiangang/firstcut/iq_20130214/ '+catname)
        #os.system('wget --no-check-certificate --user=jghao --password="jgh70chips" --directory-prefix=/home/jghao/research/data/firstcutcat/fircutMoments/ '+catname)

"""

reshao = []
resmike = []

if i ==0:
    for eid in expid:
        print eid
        reshao.append(whiskerStat_firstcut(eid))
    np.savetxt('firstcut_stat_165131_hao.txt',np.array(reshao),fmt='%10.5f')

if i == 1:
    for eid in expid:
        print eid
        resmike.append(whiskerStat_firstcut_mike(eid))
        np.savetxt('firstcut_stat_165131_mike.txt',np.array(resmike),fmt='%10.5f')


"""
