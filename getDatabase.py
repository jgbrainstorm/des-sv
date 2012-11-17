#! /usr/bin/env python
#------------------------------------------------
# This code compare the image quality measurement
# from IH and QR.
# J. Hao, 11/16/2012 @ FNAL 

import psycopg2 as psy
import numpy as np
import pylab as pl
import sys


def getExpidByDate(start_date='2012-11-08',end_date='2012-11-09'):
    """
    get info from the IH database
    """
    conn_str="host='server2.ctio.noao.edu' port=5442 dbname='decam_prd' user='decam_reader' password='reader'"
    d = psy.connect(conn_str)
    c = d.cursor()
    cmd_str = "SELECT exposure_id FROM telemetry.image_health WHERE time_recorded BETWEEN '%s' AND '%s' GROUP BY exposure_id"%(start_date,end_date)
    c.execute(cmd_str)
    res = np.array(c.fetchall())
    return res[:,0]

def getIHbyDate(start_date='2012-11-08',end_date='2012-11-09'):
    """
    get info from the IH database
    """
    conn_str="host='server2.ctio.noao.edu' port=5442 dbname='decam_prd' user='decam_reader' password='reader'"
    d = psy.connect(conn_str)
    c = d.cursor()
    cmd_str = "SELECT exposure_id,AVG(harmonic_mean_seeing[1]), AVG(median_ellipticity[1]) FROM telemetry.image_health WHERE time_recorded BETWEEN '%s' AND '%s' GROUP BY exposure_id"%(start_date,end_date)
    c.execute(cmd_str)
    res = np.array(c.fetchall())
    expid = res[:,0].astype('i')
    fwhm = res[:,1]
    ellip = res[:,2]
    return expid,fwhm,ellip

def getIHbyExpid(expid=None):
    """
    get info from the IH database
    """
    conn_str="host='server2.ctio.noao.edu' port=5442 dbname='decam_prd' user='decam_reader' password='reader'"
    d = psy.connect(conn_str)
    c = d.cursor()
    cmd_str = "SELECT exposure_id,AVG(harmonic_mean_seeing[1]), AVG(median_ellipticity[1]) FROM telemetry.image_health WHERE exposure_id='%s' GROUP BY exposure_id"%str(expid)
    c.execute(cmd_str)
    res = c.fetchall()
    if len(res)==0:
        exposure_id = expid
        fwhm = 0.
        ellip = 0.
    else:
        exposure_id = res[0][0]
        fwhm = res[0][1]
        ellip = res[0][2]
    return exposure_id,fwhm,ellip


def getDIMMbyDate(start_date='2012-11-08',end_date='2012-11-09'):
    """
    get info from the IH database
    """
    conn_str="host='server2.ctio.noao.edu' port=5442 dbname='decam_prd' user='decam_reader' password='reader'"
    d = psy.connect(conn_str)
    c = d.cursor()
    cmd_str = "SELECT en.time_recorded, en.dimm_seeing FROM telemetry.environmental_data en WHERE en.time_recorded BETWEEN '%s' AND '%s'"%(start_date,end_date)
    c.execute(cmd_str)
    res = np.array(c.fetchall())
    return res


def getQRbyDate(start_date='2012-11-08',end_date='2012-11-09'):
    """
    get info from the Quick Reduce database
    """
    conn_str="host='server2.ctio.noao.edu' port=5442 dbname='desbrdev' user='qr_reader' password='QRreader'"
    d = psy.connect(conn_str)
    c = d.cursor()
    cmd_str = "SELECT e.expnum as exposure_id, AVG(c.fwhm), AVG(c.ellipticity) FROM exposure e, ccd c, night n where e.exposure_id = c.exposure_id AND n.date BETWEEN '%s' AND '%s' GROUP BY e.expnum"%(start_date,end_date)
    c.execute(cmd_str)
    res = np.array(c.fetchall())
    expid = res[:,0].astype('i')
    fwhm = res[:,1]
    ellip = res[:,2]
    return expid,fwhm,ellip

def getQRbyExpid(expid=None):
    """
    get info from the Quick Reduce database
    """
    conn_str="host='server2.ctio.noao.edu' port=5442 dbname='desbrdev' user='qr_reader' password='QRreader'"
    d = psy.connect(conn_str)
    c = d.cursor()
    cmd_str = "SELECT e.expnum as exposure_id, AVG(c.fwhm), AVG(c.ellipticity),MAX(e.filter) FROM exposure e, ccd c WHERE e.exposure_id = c.exposure_id AND e.expnum ='%s' GROUP BY e.expnum"%str(expid)
    c.execute(cmd_str)
    res = c.fetchall()
    if len(res)==0:
        exposure_id = expid
        fwhm = 0.
        ellip = 0.
        filter='NA'
    else:
        exposure_id = res[0][0]
        fwhm = res[0][1]
        ellip = res[0][2]
        filter = res[0][3]
    return exposure_id,fwhm,ellip,filter


if __name__ == "__main__":
    import sys
    sys.path.append('/usr/remote/user/sispi/jiangang/des-sv')
    from getDatabase import *
    if len(sys.argv) <2:
        print 'syntax: getDatabase start_date end_date'
        print 'example: getDatabase.py 2012-11-06 2012-11-07'
        sys.exit()
    start_date = sys.argv[1]
    end_date = sys.argv[2]
    expid = getExpidByDate(start_date,end_date)
    nexposure = len(expid)
    fltr=[]
    fwhmIH = np.zeros(nexposure)
    ellipIH = np.zeros(nexposure)
    fwhmQR = np.zeros(nexposure)
    ellipQR = np.zeros(nexposure)
    for i in range(nexposure)
        print i
        resQR = getQRbyExpid(expid[i])
        resIH = getIHbyExpid(expid[i])
        fwhmQR[i]=resQR[1]
        ellipQR[i] = resQR[2]
        fltr.append(resQR[3])
        fwhmIH[i] = resIH[1]
        ellipIH[i] = resIH[2]
    fltr = np.array(fltr)
    idx = np.arange(nexposure)
    idxg = idx[fltr=='g']
    idxr = idx[fltr=='r']
    idxi = idx[fltr=='i']
    idxz = idx[fltr=='z']
    idxY = idx[fltr=='Y']
    pl.figure(figsize=(10,10))
    pl.plot(fwhmQR[idxg],fwhmIH[idxg],'g.',label='g')
    pl.plot(fwhmQR[idxr],fwhmIH[idxr],'r.',label='r')
    pl.plot(fwhmQR[idxi],fwhmIH[idxi],'c.',label='i')
    pl.plot(fwhmQR[idxz],fwhmIH[idxz],'m.',label='z')
    pl.plot(fwhmQR[idxY],fwhmIH[idxY],'b.',label='Y')
    pl.grid()
    pl.xlabel('FWHM from QR')
    pl.ylabel('FWHM from IH')
    pl.title(start_date+' --  '+end_date)
    pl.legend(loc='best')
    pl.plot([0,2.5],[0,2.5],'k-')
    pl.savefig('fwhm_QR_IH.png')
    pl.close()
    pl.figure(figsize=(10,10))
    pl.plot(ellipQR[idxg],ellipIH[idxg],'g.',label='g')
    pl.plot(ellipQR[idxr],ellipIH[idxr],'r.',label='r')
    pl.plot(ellipQR[idxi],ellipIH[idxi],'c.',label='i')
    pl.plot(ellipQR[idxz],ellipIH[idxz],'m.',label='z')
    pl.plot(ellipQR[idxY],ellipIH[idxY],'b.',label='Y')
    pl.grid()
    pl.plot([0,0.5],[0,0.5],'k-')
    pl.xlabel('ELLIPTICITY from QR')
    pl.ylabel('ELLIPTICITY from IH')
    pl.title(start_date+' --  '+end_date)
    pl.legend(loc='best')
    pl.savefig('ellip_QR_IH.png')
    pl.close()
   
  
