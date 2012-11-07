import psycopg2 as psy

def getIHinfo(expid=None):
    #conn_str ="host='http://system1.ctio.noao.edu' dbname='' user='' password='' "
    d = psy.connect(conn_str)
    c = d.cursor()
    cmd_str = 'SELECT exposure_id, median_ellipticity_value, harmonic_mean_seeing_value FROM image_health WHERE exposure_id =' +str(expid)
    c.execute(cmd_str)
