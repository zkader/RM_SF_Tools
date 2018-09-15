from RM_toolbox import *
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
from astropy.coordinates import SkyCoord
import astropy.units as u
import sys
from scipy.optimize import curve_fit

linfit = lambda x,a,b: a*x+b
expfit = lambda x,a,b: b*np.exp(a*x)
powfit = lambda x,a,b: b*x**a
pow2fit = lambda x,a,b,c: b + np.log(x**a + c/b)
powp2fit = lambda x,a,b,c: b*x**a + c
logfit = lambda x,a,b: a*np.log(x) + b

testfiles = glob("C:/Users/smile_000/Documents/Py_Docs/uncutRMdata/src*RMvals.txt")
#testfilesv = glob("C:/Users/smile_000/Documents/Py_Docs/testRMdata/*RMvals.v00.txt")
testfilesv0 = glob("C:/Users/smile_000/Documents/Py_Docs/RMdata/*RMvals.Ext.v03.txt")
testfilesv = glob("C:/Users/smile_000/Documents/Py_Docs/RMdata/*RMvals.v03.txt")

testhaorcfiles = "C:/Users/smile_000/Documents/Py_Docs/testhaorc/*haorc.txt"
hashfiles = "C:/Users/smile_000/Documents/Py_Docs/testsubsub/*Hash.txt"
orshfiles = "C:/Users/smile_000/Documents/Py_Docs/testsubsub/*ORsh.txt"
testsubfiles = "C:/Users/smile_000/Documents/Py_Docs/testsubsub2/*subsub.txt"
"""
testemfiles = glob("C:/Users/smile_000/Documents/Py_Docs/EMdata/*EMvals.uncor.txt")
emv,emerrv,emradecv = readRMfile(testemfiles[0])
for tf in testemfiles[1:]:
    em1,emerr1,emradec1 = readRMfile(tf)
    emv = np.append(emv,em1)
    emerrv = np.append(emerrv,emerr1)
    emradecv = np.append(emradecv,emradec1,axis=0)

testemfiles = glob("C:/Users/smile_000/Documents/Py_Docs/EMdata/*EMvals.txt")
emv1,emerrv1,emradecv1 = readRMfile(testemfiles[0])
for tf in testemfiles[1:]:
    em1,emerr1,emradec1 = readRMfile(tf)
    emv1 = np.append(emv1,em1)
    emerrv1 = np.append(emerrv1,emerr1)
    emradecv1 = np.append(emradecv1,emradec1,axis=0)

sf3,sfe3,d3,d3wid,sfn3,sfpnt3 = SFLoop(emv,emradecv,emerrv,10,15,0.01,0.001,100,True)
sf2,sfe2,d2,d2wid,sfn2,sfpnt2 = SFLoop(emv1,emradecv1,emerrv1,10,15,0.01,0.001,100,True)

plt.figure()
plt.errorbar(d2,sfpnt2,xerr=d2wid,fmt='s',color='r')
plt.errorbar(d3,sfpnt3,xerr=d3wid,fmt='.',color='k')
plt.xscale('log')
plt.yscale('log')
plt.tight_layout()
plt.show()

plt.figure()
plt.errorbar(d2*(np.pi/180.0)*50000,sf2,yerr=sfe2,xerr=d2wid*(np.pi/180.0)*50000,fmt='.',color='k',label='Calibrated Data')
plt.xlabel('$\delta \\theta $ [$pc$]',fontsize=18)
plt.ylabel('$D_{EM}(\delta \\theta) $ [$pc^2/cm^12$]',fontsize=18)
plt.yscale('log')
plt.xscale('log')
plt.legend(loc=0,numpoints=1)
plt.tight_layout()
plt.show()


plt.figure()
plt.errorbar(d3*(np.pi/180.0)*50000,sf3,yerr=sfe3,xerr=d3wid*(np.pi/180.0)*50000,fmt='.',color='k',label='Uncalibrated Data')
plt.xlabel('$\delta \\theta $ [$pc$]',fontsize=18)
plt.ylabel('$D_{EM}(\delta \\theta) $ [$???$]',fontsize=18)
plt.yscale('log')
plt.xscale('log')
plt.legend(loc=0,numpoints=1)
plt.tight_layout()
plt.show()
"""
rmcat,rmerrcat,radeccat = readcatfile('testRMdata/lmccatdat.txt')
rm3,rmerr3,radec3 = readlistfile("testRMdata/rm_list.txt")

rmv,rmerrv,radecv = readRMfile(testfilesv[0])
rmv0,rmerrv0,radecv0 = readRMfile(testfilesv0[0])
#extsrc = ['src1080', 'src1119', 'src1130', 'src1184', 'src1192', 'src1206', 'src1211', 'src1216', 'src1221', 'src1225', 'src1247', 'src1251', 'src1256', 'src1290']
#1247 try
extsrc = ['src1080', 'src1119', 'src1130', 'src1192', 'src1216', 'src1221', 'src1247', 'src1290']
extsrc = np.array(extsrc)

testsrc = '0' #src1221-src1211, src1119,src1192
for tf in testfilesv[1:]:
    #if tf.split('\\')[-1].split('.')[0] in ['src1211',testsrc]: #16,22
    #    continue
    #if tf.split('\\')[-1].split('.')[0] in extsrc:#['src1119','src1290','src1211','src1247','src1256','src1165','src1251']:
    #    continue
    rm1,rmerr1,radec1 = readRMfile(tf)
    rmv = np.append(rmv,rm1)
    rmerrv = np.append(rmerrv,rmerr1)
    radecv = np.append(radecv,radec1,axis=0)

for tf in testfilesv0[1:]:
    #if tf.split('\\')[-1].split('.')[0] in ['src1221',testsrc]:
    #    continue
    #if tf.split('\\')[-1].split('.')[0] in extsrc:#['src1119','src1290','src1211','src1247','src1256','src1165','src1251']:
    #    continue
    rm1,rmerr1,radec1 = readRMfile(tf)
    rmv0 = np.append(rmv0,rm1)
    rmerrv0 = np.append(rmerrv0,rmerr1)
    radecv0 = np.append(radecv0,radec1,axis=0)

rm = np.append(rmv,rmv0)
rmerr = np.append(rmerrv,rmerrv0)
radec = np.append(radecv,radecv0,axis=0)

#chang those noise numbers
prebin1 = 10#8
#prebin1 = 9
postbin1 = 15#19

prebin2 = 10
postbin2 = 15

nbin3 = 60

#catmask1 = PositionFilter(radecv,radeccat)
#catmask2 = PositionFilter(radecv0,radeccat)
catmask3 = PositionFilter(radec,radeccat)

"""
rmv = np.append(rmv,rmcat[catmask1])
rmerrv = np.append(rmerrv,rmerrcat[catmask1])
radecv = np.append(radecv,radeccat[catmask1,:],axis=0)

rmv0 = np.append(rmv0,rmcat[catmask2])
rmerrv0 = np.append(rmerrv0,rmerrcat[catmask2])
radecv0 = np.append(radecv0,radeccat[catmask2,:],axis=0)
"""
rm3 = np.append(rm,rmcat[catmask3])
rmerr3 = np.append(rmerr,rmerrcat[catmask3])
radec3 = np.append(radec,radeccat[catmask3,:],axis=0)

rmcat = rmcat[catmask3]
rmerr = rmerrcat[catmask3]
radeccat = radeccat[catmask3,:]

"""
catmask = PositionFilter(radecv,radeccat)
rmv = np.append(rmv,rmcat[catmask])
rmerrv = np.append(rmerrv,rmerrcat[catmask])
radecv = np.append(radecv,radeccat[catmask,:],axis=0)
"""
#sf2,sfe2,d2,d2wid,sfn2,sfpnt2 = SFLoop(rmcat,radeccat,rmerrcat,prebin2,postbin2,0.01,0.001,100)
#print "--------------------------break------------------------------"
#sf1,sfe1,d1,d1wid,sfn1,sfpnt1 = SFLoop(rm,radec,rmerr,prebin1,postbin1,0.01,0.0001,100)
#sf2,sfe2,d2,d2wid,sfn2,sfpnt2 = SFLoop(rmv0,radecv0,rmerrv0,prebin1,postbin1,0.01,0.001,100)
#sf2,sfe2,d2,d2wid,sfn2,sfpnt2 = SFLoop(rmv0,radecv0,rmerrv0,prebin1,postbin1,0.01,0.001,100,True)
#print "--------------------------break------------------------------"
#sf1,sfe1,d1,d1wid,sfn1,sfpnt1 = SFLoop(rmv,radecv,rmerrv,prebin1,postbin1,0.01,0.001,100)
#sf1,sfe1,d1,d1wid,sfn1,sfpnt1 = SFLoop(rmv,radecv,rmerrv,prebin1,postbin1,0.01,0.001,100,True)
#print "--------------------------break------------------------------"
sf3,sfe3,d3,d3wid,sfn3,sfpnt3 = SFLoop(rm3,radec3,rmerr3,prebin1,postbin1,0.01,0.001,100,True)
#sf3,sfe3,d3,d3wid,sfn3,sfpnt3 = SFLoop(rm3,radec3,rmerr3,prebin1,postbin1,10**(-1.8),0.001,100,True)
#sf3,sfe3,d3,d3wid,sfn3,sfpnt3 = SFLoop(rm3,radec3,rmerr3,prebin1,postbin1,0.01,0.001,100,True)
print "--------------------------break------------------------------"
#try 0.01
"""
##try all combined (seems to have ra dec not match) 
newrm = np.append(rm,rmcat)
newradec = np.append(radec,radeccat,axis=0)
newrmerr = np.append(rmerr,rmerrcat)
newrm = np.append(newrm,rm3)
newradec = np.append(newradec,radec3,axis=0)
newrmerr = np.append(newrmerr,rmerr3)
sf4,sfe4,d4,d4wid,sfn4,sfpnt4 = SFLoop(newrm,newradec,newrmerr,prebin2,postbin2,0,0.001,100)
"""

plt.figure()
#plt.errorbar(d1,sfpnt1,xerr=d1wid,fmt='^',color='b')
#plt.axhline(y=2,color='k')
#plt.errorbar(d2,sfpnt2,xerr=d2wid,fmt='s',color='r')
#plt.axhline(y=1,color='r')
plt.errorbar(d3,sfpnt3,xerr=d3wid,fmt='.',color='k')
#plt.axhline(y=2,color='g')
plt.axvline(x=6.0/(3600),color='k')
#plt.errorbar(d4,sfpnt4,xerr=d4wid,fmt='.',color='b')
#plt.axhline(y=9,color='b')
plt.xscale('log')
plt.yscale('log')
plt.tight_layout()
plt.show()

"""
##pnt cut
pntcut1 = np.logical_and(np.logical_and(np.greater(d1,(6.0/3600)),np.greater(sfpnt1,0)),np.greater(sf1,0))
sf1 = sf1[pntcut1]
sfe1 = sfe1[pntcut1]
d1 = d1[pntcut1]
d1wid = d1wid[pntcut1]
sfn1 = sfn1[pntcut1]

pntcut2 = np.logical_and(np.logical_and(np.greater(d2,(6.0/3600)),np.greater(sfpnt2,0)),np.greater(sf2,0))
sf2 = sf2[pntcut2]
sfe2 = sfe2[pntcut2]
d2 = d2[pntcut2]
d2wid = d2wid[pntcut2]
sfn2 = sfn2[pntcut2]
"""
pntcut3 = np.logical_and(np.logical_and(np.greater(d3,(6.0/3600)),np.greater(sfpnt3,0)),np.greater(sf3,0))
sf3 = sf3[pntcut3]
sfe3 = sfe3[pntcut3]
d3 = d3[pntcut3]
d3wid = d3wid[pntcut3]
sfn3 = sfn3[pntcut3]

#pntcut4 = np.greater(sfpnt4,9)
#sf4 = sf4[pntcut4]
#sfe4 = sfe4[pntcut4]
#d4 = d4[pntcut4]
#d4wid = d4wid[pntcut4]
#sfn4 = sfn4[pntcut4]


#preb = len(np.where(np.less(d1,0.03))[0])
#preb2 = len(np.where(np.less(d2,1.0))[0])
preb1 = np.less(d3,0.01)
posb1 = np.logical_not(preb1)
preb2 = np.less(d3,0.1)
posb2 = np.logical_not(preb2)
preb3 = np.less(d3,0.029)#0.14
posb3 = np.logical_not(preb3)
powfit = lambda x,a,b: a*np.power(x,b) 
linfit = lambda x,a,b: a + b*x
constfit = lambda x,b: b
kolvfit = lambda x,a: a + (5/3.0)*x

popt1,pcov1 = curve_fit(kolvfit,np.log10(d3[preb1]),np.log10(sf3[preb1]),sigma=(sfe3/(sf3*np.log(10)))[preb1])
kolvamp = popt1[0]
print popt1,pcov1

popt1,pcov1 = curve_fit(linfit,np.log10(d3[preb1]),np.log10(sf3[preb1]),sigma=(sfe3/(sf3*np.log(10)))[preb1])

print popt1[0],"+/-",np.sqrt(np.abs(pcov1))[0,0]
print popt1[1],"+/-",np.sqrt(np. abs(pcov1))[1,1]


sfconst1 = np.sum(np.power(sfe3[posb1],-2)*sf3[posb1])/np.sum(np.power(sfe3[posb1],-2))
sfvar1 = (np.sum(np.power(sfe3[posb1],-2)*(sf3[posb1]-sfconst1)**2)/float(len(sfe3[posb1])-1))/(np.sum(np.power(sfe3[posb1],-2))/float(len(sfe3[posb1])))

print sfconst1,'+/-',np.sqrt(sfvar1)
print (sfconst1/2)**(0.5),'+/-',0.5*((np.sqrt(sfvar1)/2)/(sfconst1/2))*(sfconst1/2)**(0.5)

"""
popt2,pcov2 = curve_fit(linfit,np.log10(d3[preb2]),np.log10(sf3[preb2]),sigma=(sfe3/(sf3*np.log(10)))[preb2])

print 10**popt2[0],"+/-",100*np.log(10)*np.sqrt(np.abs(pcov2))[0,0],"%"
print popt2[1],"+/-",100*np.sqrt(np.abs(pcov2))[1,1]/popt2[1],"%"
sfconst2 = np.sum(np.power(sfe3[posb2],-2)*sf3[posb2])/np.sum(np.power(sfe3[posb2],-2))

popt3,pcov3 = curve_fit(linfit,np.log10(d3[preb3]),np.log10(sf3[preb3]),sigma=(sfe3/(sf3*np.log(10)))[preb3])

print 10**popt3[0],"+/-",100*np.log(10)*np.sqrt(np.abs(pcov3))[0,0],"%"
print popt3[1],"+/-",100*np.sqrt(np.abs(pcov3))[1,1]/popt3[1],"%"
sfconst3 = np.sum(np.power(sfe3[posb3],-2)*sf3[posb3])/np.sum(np.power(sfe3[posb3],-2))
"""
xtrans = 10**((np.log10(sfconst1) - popt1[0])/popt1[1])

xtrans2 = 10**((np.log10(sfconst1) - kolvamp)/(5/3.0))


def Getlinval(x,y):
    A = (y[1]-y[0])/(x[1]-x[0])
    B = 10**(y[1] - A*x[1])
    return A,B

def DataRead(fname):
    f = open(fname,'r')
    x = []
    y = []
    f.readline()
    for line in f:  
        tline = map(float,line.split(','))
        x.append(tline[0])
        y.append(tline[1])
    f.close()
    return np.array(x),np.array(y)

def SlopeFunc(xvals,xint,pwr,const):
    xx,yy,zz = np.meshgrid(xint,pwr,xvals)
    ans = const*np.power(xx,-yy)*np.power(zz,yy)
    return ans

def Chi2(yval,yerr,yfvals):
    ans = np.sum(np.power((yval - yfvals)/yerr,2))
    return ans

def Model(xvals,xint,pwr):
    lpnt = np.less(xvals,xint)
    gpnt = np.logical_not(lpnt)
    ans = np.empty(shape=xvals.shape)
    avy = AVSFval(xint)
    amp = avy*xint**(-pwr)
    ans[lpnt] = amp*np.power(xvals[lpnt],pwr) 
    ans[gpnt] = avy
    return ans

def AVSFval(xint):
    lpnt = np.less(td3,xint)
    gpnt = np.logical_not(lpnt)
    avy = np.sum(np.power(tsfe3[gpnt],-2)*tsf3[gpnt])/np.sum(np.power(tsfe3[gpnt],-2))
    return avy

td3 = d3
tsf3 = sf3
tsfe3 = sfe3
pv,cov = curve_fit(Model,d3,sf3,sigma=sfe3)

trials = 5#50000
params = np.empty(shape=(trials,2))

for i in range(trials):
    xtvals = np.random.uniform(d3-d3wid,d3+d3wid)
    ytvals = np.random.normal(sf3,sfe3)
    td3 = xtvals
    tsf3 = ytvals
    tsfe3 = np.ones(shape=ytvals.shape)
    tv,tcov = curve_fit(Model,xtvals,ytvals,pv)
    params[i] = tv


sfslopeav = params[:,1].mean()
sfslopestd = params[:,1].std()

sfxav = params[:,0].mean()*(np.pi/180)*50000
sfxstd = params[:,0].std()*(np.pi/180)*50000


plt.figure()
plt.hist(params[:,1],100)
plt.title('Parameter $b$')
plt.xlabel('Parameter $b$')
plt.tight_layout()
plt.show()


plt.figure()
plt.hist(params[:,0]*(np.pi/180)*50000,100)
plt.title('Parameter $\delta \\theta_0$')
plt.xlabel('Parameter $\delta \\theta_0$ [$pc$]')
plt.tight_layout()
plt.show()


plt.figure()
plt.scatter(params[:,0]*(np.pi/180)*50000,params[:,1],s=2,edgecolor='none')
plt.xlabel('Transition Point $\delta \\theta_{0} [$pc$] $',fontsize=18)
plt.ylabel('Power $b$',fontsize=18)
plt.axhline(y=sfslopeav,color='k')
plt.axhline(y=sfslopeav-sfslopestd,color='r')
plt.axhline(y=sfslopeav+sfslopestd,color='r')
plt.axvline(x=sfxav,color='k')
plt.axvline(x=sfxav-sfxstd,color='r')
plt.axvline(x=sfxav+sfxstd,color='r')
plt.show()

td3 = d3
tsf3 = sf3
tsfe3 = sfe3
"""
SFchi2 = np.empty((len(randb),len(randx)))
xx,bb = np.meshgrid(randx,randb)
for i in range(len(randx)):
    for j in range(len(randb)):
        SFchi2[i,j] = Chi2(sf3,sfe3,Model(d3,sf3,sfe3,randx[i],randb[j]))

plt.figure()
plt.contourf(xx,bb,(SFchi2-np.amin(SFchi2)))
plt.colorbar()
plt.tight_layout()
plt.show()
"""

radc = np.pi/180.0
#Mao et al something
##disk
Mk15x1 = np.array([-2.1058,-1.3528,-1.2425,-0.8756])
Mk15y1 = np.array([3.7854,3.8567,3.9369,4.4464])
Mk15file = glob("C:/Users/smile_000/Documents/COOP/uoft/fig19maodat*.txt")[0]
Mk15x2,Mk15y2 = DataRead(Mk15file)

Mk15_ds = radc*7000000


#Still et al 2010
sk10rgn = ["NGP","SGP","Orion","Anti-Orion","Region C'","Anti region C'","Region A'","Anti region A'","Fan","Anti Fan","Gum"]
sk10A2 = [2.40,2.67,2.92,2.67,3.06,3.15,2.81,2.68,2.78,2.98,3.72]
sk10alph2 = [0.05,0.02,0.5,0.31,0.39,0.58,0.51,0.38,0.58,0.59,1.14]
sk10color = ['darkblue','crimson','g','darkblue','crimson','g','darkblue','crimson','g','darkblue','crimson']
sk10line = ['-','-','--','--','-.','-.',':',':','-','-','--']
sk10size = [1,1,1,1,1,1,1,1,2,2,3]
sk10arr = np.logspace(-1,2,10)
sk10_ds = 2000*radc

#haverkorn et al 2005
##carinaarm
hk5_ds = 2000*radc
hk5y1 = np.array([5.2744,5.2652])
hk5x1 = np.array([-0.3033,0.398])
hk51A,hk51B = Getlinval(hk5x1,hk5y1)
hk51_ds = radc*14500
##i1
hk5y2 = np.array([4.5949,4.6718])
hk5x2 = np.array([-0.1652,0.5285])
hk52A,hk52B = Getlinval(hk5x2,hk5y2)
hk52_ds = radc*11500
##i2
hk5y3 = np.array([4.7062,5.1869])
hk5x3 = np.array([-0.2814,0.487])
hk53A,hk53B = Getlinval(hk5x3,hk5y3)
hk53_ds = radc*16000
##i3
hk5y4 = np.array([4.7617,5.3372])
hk5x4 = np.array([-0.2177,0.4705])
hk54A,hk54B = Getlinval(hk5x4,hk5y4)
hk54_ds = radc*18000
##cruxarm
hk5y5 = np.array([4.5528,4.6322])
hk5x5 = np.array([-0.2477,0.6024])
hk55A,hk55B = Getlinval(hk5x5,hk5y5)
hk55_ds = radc*17000

#print hk52A,hk52B
#print hk53A,hk53B
#print hk54A,hk54B

hk5darr = np.logspace(-0.4,0.4,10)

#minter&spangler: D_RM = (340 +- 30) * x^(0.64 +- 0.06)
linearr = np.array([d3[0],xtrans,xtrans,d3[-1]])
linearr2 = np.array([d3[0],xtrans2,xtrans2,d3[-1]])

msy1 = np.log10(np.array([74.4183,152.5267]))
msx1 = np.log10(np.array([0.0902,0.2678]))
msy2 = np.log10(np.array([42.2496,28.9795]))
msx2 = np.log10(np.array([0.0555,0.0438]))
msA1,msB1 = Getlinval(msx1,msy1)
msA2,msB2 = Getlinval(msx2,msy2)
msxint = (msB1/msB2)**(1/(msA2-msA1))

thetaarr = np.append(np.linspace(0.02,msxint,10,endpoint=False),np.linspace(msxint,10,10))
msdrmarr = np.append(powfit(np.linspace(0.02,msxint,10,endpoint=False),msB2,msA2),powfit(np.linspace(msxint,10,10),msB1,msA1))

lmc_ds = 50000*radc
mintspang_ds = 2900*radc


#kolvamp = 10**(np.log10(sfconst1) - (5/3.0)*np.log10(xtrans) )

plt.figure()
"""
plt.plot(d1[:preb+1],powfit(d1[:preb+1],10**popt[0],popt[1]),color='b')
plt.axhline(y=sfconst,color='b')
plt.plot(d2[:preb2+1],powfit(d2[:preb2+1],10**popt2[0],popt2[1]),color='r')
plt.axhline(y=sfconst2,color='r')


plt.plot(linearr[:2]*lmc_ds,powfit(linearr[:2],10**popt3[0],popt3[1]),color='k')
plt.plot(linearr[2:]*lmc_ds,[sfconst3,sfconst3],color='k')

plt.plot(linearr[:2]*lmc_ds,powfit(linearr[:2],10**popt2[0],popt2[1]),color='r')
plt.plot(linearr[2:]*lmc_ds,[sfconst2,sfconst2],color='r')
"""

plt.plot(linearr[:2]*lmc_ds,powfit(linearr[:2],10**popt1[0],popt1[1]),color='k',label='Best Fit')
plt.plot(linearr[2:]*lmc_ds,[sfconst1,sfconst1],color='k')
plt.plot(np.logspace(np.log10(linearr[0]),np.log10(linearr[-1]),200)*lmc_ds,Model(np.logspace(np.log10(linearr[0]),np.log10(linearr[-1]),200),pv[0],pv[1]))

#plt.plot(linearr2[:2]*lmc_ds,powfit(linearr2[:2],10**kolvamp,5/3.0),color='r',label='Kolmogorov Fit')


#plt.plot(thetaarr*mintspang_ds,msdrmarr,color='r')
#plt.scatter(x1*lmc_ds,x2,s=10,color='black',label='SFs of src1192,src1247,src1290')
#plt.errorbar(d1*lmc_ds,sf1,yerr=sfe1,xerr=d1wid*lmc_ds,fmt='^',color='b',label='Central Source')
#plt.errorbar(d1,sfn1,xerr=d1wid,fmt='-.',color='k')
#plt.errorbar(d2*lmc_ds,sf2,yerr=sfe2,xerr=d2wid*lmc_ds,fmt='s',color='r',label='Ann\'s Data')
#plt.errorbar(d2,sfn2,xerr=d2wid,fmt='-.',color='r')
plt.errorbar(d3*lmc_ds,sf3,yerr=sfe3,xerr=d3wid*lmc_ds,fmt='.',color='k',label='Data')

#plt.errorbar(d3,sfn3,xerr=d3wid,fmt='-.',color='g')
#plt.errorbar(d4,sf4,yerr=sfe4,xerr=d4wid,fmt='.',color='b')

#plt.errorbar(d3,sf3,yerr=sfe3,xerr=d3wid,fmt='.',color='g')
plt.xlabel('$\delta \\theta $ [$pc$]',fontsize=18)
plt.ylabel('$D_{RM}(\delta \\theta) $ [$rad^2/m^4$]',fontsize=18)
plt.yscale('log')
plt.xscale('log')
plt.legend(loc=0,numpoints=1)
plt.tight_layout()
plt.show()
sys.exit()

plt.figure()
ax = plt.subplot(111)
ax.plot(linearr[:2]*lmc_ds,powfit(linearr[:2],10**popt1[0],popt1[1]),color='k',label='Best Fit')
ax.plot(linearr[2:]*lmc_ds,[sfconst1,sfconst1],color='k')
for i in range(len(sk10rgn)):
    ax.plot(sk10arr*sk10_ds,powfit(sk10arr,10**sk10A2[i],sk10alph2[i]),color=sk10color[i],linestyle=sk10line[i],linewidth=sk10size[i],label=sk10rgn[i])
plt.xlabel('$\delta \\theta $ [$pc$]',fontsize=18)
plt.ylabel('$D_{RM}(\delta \\theta) $ [$rad^2/m^4$]',fontsize=18)
plt.yscale('log')
plt.xscale('log')
plt.title('Comparison with Stil et al\'s work $[4]$')
box = ax.get_position()
ax.set_position([box.x0,box.y0,box.width*0.8,box.height])
ax.legend(bbox_to_anchor=(1.0,1.0),loc=2)
#plt.tight_layout()
plt.show()

plt.figure()
plt.plot(linearr[:2]*lmc_ds,powfit(linearr[:2],10**popt1[0],popt1[1]),color='k',label='Best Fit')
plt.plot(linearr[2:]*lmc_ds,[sfconst1,sfconst1],color='k')
plt.plot(thetaarr*mintspang_ds,msdrmarr,color='r',label='ISM')
plt.xlabel('$\delta \\theta $ [$pc$]',fontsize=18)
plt.ylabel('$D_{RM}(\delta \\theta) $ [$rad^2/m^4$]',fontsize=18)
plt.yscale('log')
plt.xscale('log')
plt.title('Comparison with Minter and Spangler\'s work $[3]$ ')
plt.legend(loc=0,numpoints=1)
plt.tight_layout()
plt.show()

plt.figure()
plt.plot(linearr[:2]*lmc_ds,powfit(linearr[:2],10**popt1[0],popt1[1]),color='k',label='Best Fit')
plt.plot(linearr[2:]*lmc_ds,[sfconst1,sfconst1],color='k')
plt.plot(10**Mk15x1*Mk15_ds,10**Mk15y1,label='M51 Disk')
plt.plot(10**Mk15x2*Mk15_ds,10**Mk15y2,label='M51 Halo')
plt.xlabel('$\delta \\theta $ [$pc$]',fontsize=18)
plt.ylabel('$D_{RM}(\delta \\theta) $ [$rad^2/m^4$]',fontsize=18)
plt.yscale('log')
plt.xscale('log')
plt.title('Comparison with Mao et al\'s work $[2]$ ')
plt.legend(loc=0,numpoints=1)
plt.tight_layout()
plt.show()

plt.figure()
plt.plot(linearr[:2]*lmc_ds,powfit(linearr[:2],10**popt1[0],popt1[1]),color='k',label='Best Fit')
plt.plot(linearr[2:]*lmc_ds,[sfconst1,sfconst1],color='k')
plt.plot(hk5darr*hk5_ds,powfit(hk5darr,hk55B,hk55A),label='Crux Arm')
plt.plot(hk5darr*hk5_ds,powfit(hk5darr,hk54B,hk54A),label='Interarm 3')
plt.plot(hk5darr*hk5_ds,powfit(hk5darr,hk53B,hk53A),label='Interarm 2')
plt.plot(hk5darr*hk5_ds,powfit(hk5darr,hk52B,hk52A),label='Interarm 1')
plt.plot(hk5darr*hk5_ds,powfit(hk5darr,hk51B,hk51A),label='Carina Arm')
plt.xlabel('$\delta \\theta $ [$pc$]',fontsize=18)
plt.ylabel('$D_{RM}(\delta \\theta) $ [$rad^2/m^4$]',fontsize=18)
plt.yscale('log')
plt.xscale('log')
plt.title('Comparison with Haverkorn et al\'s work $[1]$')
plt.legend(loc=0,numpoints=1)
plt.tight_layout()
plt.show()


plt.figure()
#plt.errorbar(d1,sfn1,xerr=d1wid,fmt='-^',color='b')
#plt.errorbar(d2,sfn2,xerr=d2wid,fmt='-s',color='r')
plt.errorbar(d3,sfn3,xerr=d3wid,fmt='-.',color='k')
plt.xlabel('$\delta \\theta $ [$\degree$]',fontsize=18)
plt.ylabel('$D_{RM,noise}(\delta \\theta) $ [$rad^2/m^4$]',fontsize=18)
plt.yscale('log')
plt.xscale('log')
plt.tight_layout()
plt.show()
