import psycopg2 as psy
import numpy as np

def getIHseeing(start_date='2012-11-08',end_date='2012-11-09'):
    """
    get info from the IH database
    """
    conn_str="host='server2.ctio.noao.edu' port=5442 dbname='decam_prd' user='decam_reader' password='reader'"
    d = psy.connect(conn_str)
    c = d.cursor()
    cmd_str = "SELECT exposure_id,AVG(harmonic_mean_seeing[1]) FROM telemetry.image_health WHERE time_recorded BETWEEN '%s' AND '%s' GROUP BY exposure_id"%(start_date,end_date)
    #cmd_str = "SELECT exposure_id, AVG(harmonic_mean_seeing[1]) FROM telemetry.image_health GROUP BY exposure_id LIMIT 2"
    c.execute(cmd_str)
    res = np.array(c.fetchall())
    return res

    

def getDIMM(start_date='2012-11-08',end_date='2012-11-09'):
    """
    get info from the IH database
    """
    conn_str="host='server2.ctio.noao.edu' port=5442 dbname='decam_prd' user='decam_reader' password='reader'"
    d = psy.connect(conn_str)
    c = d.cursor()
    cmd_str = "SELECT en.time_recorded, en.dimm_seeing FROM telemetry.environmental_data en WHERE en.time_recorded BETWEEN '%s' AND '%s'"%(start_date,end_date)
    c.execute(cmd_str)
    res = np.array(c.fetchall())
    


def getQRseeing(start_date='2012-11-08',end_date='2012-11-09'):
    """
    get info from the Quick Reduce database
    """
    conn_str="host='server2.ctio.noao.edu' port=5442 dbname='desbrdev' user='qr_reader' password='QRreader'"
    d = psy.connect(conn_str)
    c = d.cursor()
    cmd_str = "SELECT e.expnum AS exposure_id, AVG(c.fwhm) FROM exposure e, ccd c, night n where e.exposure_id = c.exposure_id AND n.date BETWEEN '%s' AND '%s' GROUP BY e.expnum"%(start_date,end_date)
    c.execute(cmd_str)
    res = np.array(c.fetchall())
    return res

#psql -h server2.ctio.noao.edu -U qr_reader -p 5442 -d desbrdev -c "SELECT e.expnum AS exposure_id, AVG(c.fwhm) FROM exposure e, ccd c, night n where e.exposure_id = c.exposure_id AND n.date BETWEEN '2012-11-06' AND '2012-11-07' GROUP BY e.expnum"

    
