import psycopg2 as psy

def getIHinfo(expid=None):
    DBHOST ='server2.ctio.noao.edu'
    DBPORT ='5442'
    DBNAME='decam_prd'
    DBUSER='decam_reader'
    DBUSER_PASSWORD='reader'
    #conn_str ="host= %s port=%s dbname=%s user=%s password=%s" %(DBHOST,DBPORT,DBNAME,DBUSER,DBUSER_PASSWORD)
    conn_str="host='server2.ctio.noao.edu' port=5442 dbname='decam_prd' user='decam_reader' password='reader'"
    d = psy.connect(conn_str)
    c = d.cursor()
    cmd_str = 'SELECT exposure_id, median_ellipticity_value, harmonic_mean_seeing_value FROM image_health WHERE exposure_id =' +str(expid)
    c.execute(cmd_str)
