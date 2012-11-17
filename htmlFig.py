#! /usr/bin/env python
import glob as gl
import sys,os

htmlName = sys.argv[1]
htmlName = 'Image_Quality_'+htmlName+'.html'
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
htm.write('<TITLE>Image Quality of 11/15/2012</TITLE>\n')
htm.write('</HEAD> \n')
htm.write('<BODY> \n')
htm.write('<p>Image Quality some r band images taken during 11/15/2012</p>\n')
htm.write('<p>J. Hao @ FNAL</p>\n')

for i in range(nfig):
    expid=Fig_fwhm[i][-12:-5]
    htm.write('exposure ID:'+expid+'\n')
    htm.write('<p>\n')
    htm.write('<img src="%s" width="1000">\n'%Fig_fwhm[i])
    htm.write('<img src="%s" width="1000">\n'%Fig_moments[i])
    htm.write('<img src="%s" width="1000">\n'%Fig_coeff[i])
    htm.write('</p>\n')
htm.write('</BODY> \n')
htm.write('</HTML> \n')
htm.close()
os.system('tar -czf '+htmlName+'.tar.gz *')
