#!/usr/bin/env python

import os,argparse
import cx_Oracle

parser = argparse.ArgumentParser(description='Run single file')
parser.add_argument('--exposure',type=int,help='exposure number',default=-1)
parser.add_argument('--file',help='exposures from a file',
                    default="")
parser.add_argument('--start',type=int,default=-1,
                    help='starting exposure number')
parser.add_argument('--end',type=int,help='end exposure number')
parser.add_argument('--run',default='',help='run')
parser.add_argument('--outdir',default='',help='run')

args = parser.parse_args()
exp_list=[]

if args.exposure>0:
    exp_list.append(args.exposure)
elif args.file!="":
    file=open(args.file)
    for line in file:
        exp=int(line)
        exp_list.append(exp)
else:
     exp_list=range(args.start,args.end+1)


for exp in exp_list:

    print "Downloading exposure %d"%exp
    
    connection = cx_Oracle.connect('[user]/[password]@leovip148.ncsa.uiuc.edu/desoper')
    cursor = connection.cursor()
    run=args.run
    if run=="":
        query="select distinct(run) from exposure,image where image.imagetype='red' and image.exposureid=exposure.id and expnum=%d order by run desc" % exp
        
        
        cursor.execute(query)
        results=cursor.fetchall();
        if len(results)==0:
            print 'Could not find exposure %d'%exp
            continue
        if len(results)>1:
            print 'Found multitple runs using latest'
                
        run=results[0][0]

    outdir='%s/DECam_%d' % (args.outdir,exp)

        
    loc=('ftp://desar.cosmology.illinois.edu/DESFiles/desardata/OPS/red/%s/red/DECam_%08d/*' % (run,exp))
    print loc

    if not os.path.exists(outdir): os.makedirs (outdir)
    os.system('wget -q -P %s %s ' % (outdir,loc))

