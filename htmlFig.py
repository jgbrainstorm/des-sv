#! /usr/bin/env python
import glob as gl
import sys,os

if len(sys.argv) == 1:
    print 'syntax: htmlFig.py date MJD'
    print 'example: htmlFig.py 11-17-2012 56248'
    sys.exit()
date = sys.argv[1]
mjd = sys.argv[2]
htmlName = 'Image_Quality_'+date+'-MJD-'+mjd+'.html'
Fig_coeff = gl.glob('zernike_coeff*.png')
Fig_moments = gl.glob('moments*.png') 
Fig_fwhm = gl.glob('fwhm*.png') 

Fig_coeff.sort()
Fig_moments.sort()
Fig_fwhm.sort()


nfig = len(Fig_coeff)

htm=open(htmlName,'w')
htm.write('<HTML> \n')

htm.write('<HEAD> \n')
htm.write('<TITLE>Image Quality of '+date+'</TITLE>\n')
htm.write('</HEAD> \n')
htm.write('<BODY> \n')
htm.write('<p>Image Quality some r band images taken during: '+date+'</p>\n')
htm.write('<p>MJD: '+mjd+'</p>\n')
htm.write('<p>J. Hao @ FNAL</p>\n')

htm.write('<p>Hexapod Parameter Summary </p>\n')
htm.write('<img src="%s" width="1000">\n'%'hexapod_pos_summary.png')
for i in range(nfig):
    expid=Fig_fwhm[i][-12:-4]
    htm.write('exposure ID:'+expid+'\n')
    htm.write('<p>\n')
    htm.write('<img src="%s" width="1000">\n'%Fig_fwhm[i])
    htm.write('<img src="%s" width="1000">\n'%Fig_moments[i])
    htm.write('<img src="%s" width="1000">\n'%Fig_coeff[i])
    htm.write('</p>\n')
htm.write('</BODY> \n')
htm.write('</HTML> \n')
htm.close()
os.system('tar -czf '+htmlName+'.tar.gz *.png *.html')
os.system('tar -czf hexpod_parameter_'+date+'-MJD-'+mjd+'.tar.gz *.txt')
