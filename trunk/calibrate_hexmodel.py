from decamImgAnalyzer import *

f = gl.glob('hexapod_position_*.txt')
fzerni = gl.glob('zernike_coeff_*.p')
fzerni.sort()
f.sort()
expid=[]
data=[]
nexp = len(f)
coeff=[]
#if nexp == 0:
#    return '-- no image to analyze, exit --'
for fi in f:
    expid.append(int(fi[-12:-4]))
    data.append(np.genfromtxt(fi))
    expid = np.array(expid)
    xtick = expid.astype('S10')
data = np.array(data)
for fz in fzerni:
    coeff.append(p.load(open(fz),'r'))
coeff = np.array(coeff)

    
    pl.figure(figsize=(16,9))
    pl.subplot(5,1,1)
    pl.plot(expid,data[:,1,0],'bo-',label='Image Analysis')
    pl.plot(expid,data[:,2,0],'ro-',label='BCAM')
    pl.plot(expid,data[:,0,0],'go-',label='Hexapod ')
    pl.ylabel('x-decenter')
    pl.xticks(expid,np.repeat('',nexp))
    pl.legend(loc='upper center',bbox_to_anchor=(0.5,1.68))
    pl.grid()
    pl.subplot(5,1,2)
    pl.plot(expid,data[:,1,1],'bo-',label='Image Analysis')
    pl.plot(expid,data[:,2,1],'ro-',label='BCAM')
    pl.plot(expid,data[:,0,1],'go-',label='Hexapod ')
    pl.ylabel('y-decenter')
    pl.xticks(expid,np.repeat('',nexp))
    pl.grid()
    pl.subplot(5,1,3)
    pl.plot(expid,data[:,1,2],'bo-',label='Image Analysis')
    pl.plot(expid,data[:,0,2],'go-',label='Hexapod ')
    #pl.plot(expid,data[:,2,2],'ro',label='BCAM') #no focus for BCAM
    pl.ylabel('z-defocus')
    pl.xticks(expid,np.repeat('',nexp))
    pl.grid()
    pl.subplot(5,1,4)
    pl.plot(expid,data[:,1,3],'bo-',label='Image Analysis')
    pl.plot(expid,data[:,2,3],'ro-',label='BCAM')
    pl.plot(expid,data[:,0,3],'go-',label='Hexapod ')
    pl.ylabel('x-tilt')
    pl.xticks(expid,np.repeat('',nexp))
    pl.grid()
    pl.subplot(5,1,5)
    pl.plot(expid,data[:,1,4],'bo-',label='Image Analysis')
    pl.plot(expid,data[:,2,4],'ro-',label='BCAM')
    pl.plot(expid,data[:,0,4],'go-',label='Hexapod ')
    pl.xlabel('exposure_id')
    pl.ylabel('y-tilt')
    pl.xticks(expid,xtick,rotation=90)
    pl.grid()
 
