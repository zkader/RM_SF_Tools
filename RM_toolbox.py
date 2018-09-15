import numpy as np
import matplotlib.pyplot as plt
from glob import glob
from astropy.coordinates import SkyCoord
import astropy.units as u
import sys
from scipy.optimize import curve_fit

def sigfig(val,err):
    return round(val,-int(np.floor(np.log10(abs(err)))))

def latexSI(val,err,unit="",sfig=True):
    if sfig:
        val = sigfig(val,err)
        err = sigfig(err,err)
    if err == 0:
        return '\\SI{'+str(int(val))+'}{'+unit+'}'
    elif np.log10(abs(err)) < 0:
        return '\\SI{'+str(int(val))+'('+ str(err)[-1] +')}{'+unit+'}'
    else:
        return '\\SI{'+str(int(val))+'('+ str(int(err))+')}{'+unit+'}'

def readvecfile(fname):
    tf1 = open(fname)
    t1 = []
    t2 = []
    rarr = []
    darr = []
    narr = []
    i = 0
    j = 0
    phisize = 113
    str0,str1 = 0,0
    for line in tf1:
        tline = map(float,line.split())
        t1.append(tline[2])
        t2.append(tline[3])
        if str0 != tline[0] or str1 != tline[1]:
            str0 = tline[0]
            str1 = tline[1]
            narr.append(tline[4])
            rarr.append(str0)
            darr.append(str1)
            if i != 0:
                j = i
        if j == 0:
            i += 1
    if j== 0:
        j = i
    tf1.close()
    vec1 = np.array(t1)
    vec2 = np.array(t2)
    #vec1.reshape((len(rarr),vec1.size/len(rarr)))
    #vec2.reshape((len(rarr),vec1.size/len(rarr)))
    vec1 = vec1.reshape((vec1.size/j,j))
    vec2 = vec2.reshape((vec1.size/j,j))
    #worldcoord = SkyCoord(rarr,darr,unit=(u.deg,u.deg),frame='fk5')
    #worldcoord = worldcoord.to_string('hmsdms')
    srcstr = fname.split('\\')[-1].split('.')[0]
    worldcoord = []
    for i in range(len(rarr)):
        tstr = srcstr +" @  x: " + str(rarr[i]) + "  y: " + str(darr[i])
        worldcoord.append(tstr)
    #for i in range(len(worldcoord)):
    #    worldcoord[i] = srcstr + " : " + worldcoord[i]
    return vec1,vec2,worldcoord,np.array(narr)

def readfile(fname):
    tf1 = open(fname)
    d_arr = []
    sf_arr = []
    for line in tf1:
        tline = map(float,line.split(" ")[0:2])
        d_arr.append(tline[0])
        sf_arr.append(tline[1])
    tf1.close()
    return np.array(d_arr),np.array(sf_arr)

def readlistfile(fname):
    tf1 = open(fname)
    rm_arr = []
    rmerr_arr = []
    ra_arr = []
    dec_arr = []
    tf1.readline()
    for line in tf1:
        tline = map(float,line.split())
        rm_arr.append(tline[1])
        rmerr_arr.append(tline[2])
        ra_arr.append(tline[3])
        dec_arr.append(tline[4])
    tf1.close()
    worldcoord = np.array(list(zip(ra_arr,dec_arr)))
    lmcoord = SkyCoord("05:16:03 -68:41:45",unit=(u.hourangle,u.deg),frame='fk5')
    tworld = SkyCoord(worldcoord,unit=(u.deg,u.deg),frame='fk5')    
    tsep = lmcoord.separation(tworld).degree
    sepcut = np.less(tsep,3)
    rm_arr = np.array(rm_arr)[sepcut]
    rmerr_arr = np.array(rmerr_arr)[sepcut]
    worldcoord = worldcoord[sepcut,:]
    return rm_arr,rmerr_arr,worldcoord

def readcatfile(fname):
    tf1 = open(fname)
    rm_arr = []
    rmerr_arr = []
    ra_arr = []
    dec_arr = []
    #radec_strs = []
    for line in tf1:
        tline = line.split()
        rm_arr.append(float(tline[6]))
        rmerr_arr.append(float(tline[7]))
        ra_arr.append(":".join(tline[0:3]))
        dec_arr.append(":".join(tline[3:6]))
        #radec_strs.append(":".join(tline[0:3]) + " " + ":".join(tline[3:6]))
    tf1.close()
    c = SkyCoord(ra_arr,dec_arr,unit=(u.hourangle,u.deg),frame='fk5')
    worldcoord = np.array(list(zip(c.ra.degree,c.dec.degree)))
    
    return np.array(rm_arr),np.array(rmerr_arr),worldcoord#,radec_strs

def readRMfile(fname):
    tf1 = open(fname)
    rm_arr = []
    rmerr_arr = []
    ra_arr = []
    dec_arr = []
    for line in tf1:
        tline = map(float,line.split(" ")[0:4])
        rm_arr.append(tline[0])
        rmerr_arr.append(tline[1])
        ra_arr.append(tline[2])
        dec_arr.append(tline[3])
    tf1.close()
    worldcoord = np.array(list(zip(ra_arr,dec_arr)))    
    return np.array(rm_arr),np.array(rmerr_arr),worldcoord


def GenSF(aarr,warr,earr):
    x0 = warr[0,0]
    y0 = warr[0,1]
    p0 = SkyCoord(x0*u.deg,y0*u.deg,frame='fk5')
    a0 = aarr[0]
    e0 = earr[0]

    SFarr = np.power(aarr[1:] - a0,2)
    SFEarr = np.power(earr[1:] - e0,2)

    #SFErrarr = np.sqrt(e0*e0 + np.power(earr[1:],2))
    #SFErrarr = (SFErrarr/(aarr[1:] - a0))*2.0*SFarr
    parr = SkyCoord(warr[1:,0]*u.deg,warr[1:,1]*u.deg,frame='fk5')
    darr = parr.separation(p0).degree
    
    return SFarr,darr,SFEarr#,SFErrarr

def equalN(x, nbin):
    npt = len(x)
    x= np.sort(x)
    #nbinarr = np.linspace(0, npt, nbin + 1)
    #ans = np.interp(xarr,np.arange(npt),np.sort(x))
    xarr = [x[0]]
    parr = [0]
    for i in range(1,npt):
        pchk = (np.abs(xarr[-1] -x[i]))/x[i]
        if pchk >= 1.0e-3:
            if  np.abs(i - parr[-1]) >=  nbin:
                xarr.append(x[i])
                parr.append(i)
    xarr = np.array(xarr)
    if (npt - parr[-1]) < nbin:
        xarr = xarr[:-1]
    return xarr

def bootstrap_resample(X,n=None):
    if n == None:
        n = len(X)        
    resample_i = np.floor(np.random.rand(n)*len(X)).astype(int)
    X_resample = X[resample_i]
    return X_resample

def mp_bootstrap(yarr,xarr,xbarr,itnum=500,nbin=1):
    SFunval,pntvals,SFEunval = AverageSF(xarr,xbarr,yarr,xarr)
    SFsampval = np.ndarray(shape=(itnum,xbarr.shape[0]-1))
    for i in range(itnum):
        SFsampval[i] = bootAvSF(xarr,xbarr,yarr,xarr)
    SFsampmean = np.sum(SFsampval,axis=0)/float(itnum)
    SFsamperr = np.sqrt(np.sum((SFsampval-SFsampmean)**2,axis=0)/float(itnum-1))
    minpoint = np.where(np.less_equal(pntvals,nbin))[0]
    SFsamperr[minpoint] = SFEunval[minpoint]
    return SFunval,SFsamperr

def bootAvSF(xarr,eqxarr,yarr,yerr):
    avSFvals = []
    for i in range(eqxarr.size):
        if i == eqxarr.size - 2:
            logmask = np.greater_equal(xarr,eqxarr[i])
        elif i == eqxarr.size - 1:
            break
        else:
            logmask = np.logical_and(np.greater_equal(xarr,eqxarr[i]),np.less(xarr,eqxarr[i+1]))
        tvals = yarr[logmask]
        tsize = len(tvals)
        if tsize == 0:
            avSFvals += [0]
            continue
        bootarr = np.random.randint(tvals.shape[0],size=tvals.shape[0])
        tvals = tvals[bootarr]
        unitvals,unicnts = np.unique(tvals,return_counts=True)
        if tsize > 1 :
            tmean = np.sum(unicnts*unitvals)/np.sum(unicnts)
            avSFvals += [tmean]
        elif tsize == 1:
            avSFvals += [tvals[0]]
    return avSFvals

def AverageSF(xarr,eqxarr,yarr,yerr):
    avSFvals = []
    pntvals = []
    avSFEvals = []

    for i in range(eqxarr.size):
        if i == eqxarr.size - 2:
            logmask = np.greater_equal(xarr,eqxarr[i])
        if i == eqxarr.size - 1:
            break
        else:
            logmask = np.logical_and(np.greater_equal(xarr,eqxarr[i]),np.less(xarr,eqxarr[i+1]))
        tvals = yarr[logmask]
        pntvals.append(len(tvals))
        if len(tvals) > 1 :
            #tevals = np.power(yerr[logmask],-2.0)
            tmean = np.sum(tvals)/len(tvals)
            avSFvals += [tmean]
            tstd = np.sum((tvals-tmean)**2)/(len(tvals) -1.0)
            avSFEvals += [tstd]
            #twnum = np.sum(tevals*np.power(tvals - tmean,2))/(len(tevals) - 1.0)
            #avSFEvals += [np.sqrt(np.divide(twnum,np.sum(tevals)/float(len(tevals))))]
        elif len(tvals) == 1:
            avSFvals += [tvals[0]] 
            avSFEvals += [yerr[logmask][0]]
        else:
            avSFvals += [0]
            avSFEvals += [0]
    return np.array(avSFvals),np.array(pntvals),np.array(avSFEvals)

def SFLoop(aarr,warr,earr,prebin=10,postbin=10,cutval=0.1,begval=0.001,endval=100,eqprebin=False):    
    SFvals,davals,SFEvals = GenSF(aarr,warr,earr)
    #SFvals,davals,SFEvals,SFErrvals = GenSF(aarr,warr,earr)
    
    for i in range(1,aarr.size-1):
        tempsf,tempda,tempsfe = GenSF(aarr[i:],warr[i:,:],earr[i:])
        #tempsf,tempda,tempsfe,tempsferr = GenSF(aarr[i:],warr[i:,:],earr[i:])
        SFvals = np.append(SFvals,tempsf)
        davals = np.append(davals,tempda)
        SFEvals = np.append(SFEvals,tempsfe)        
        #SFErrvals = np.append(SFErrvals,tempsferr)

    dasort = np.argsort(davals)
    davals = davals[dasort]
    SFvals = SFvals[dasort]
    SFEvals = SFEvals[dasort]

    if cutval == 0:
        eqdavals = np.logspace(np.log10(begval),np.log10(endval),num=prebin+postbin)        
    else:
        if eqprebin:
            degcut = np.less(davals, cutval)
            eqda1vals = equalN(davals[degcut],prebin)
            eqda2vals = np.logspace(np.log10(cutval),np.log10(endval),num=postbin)
            dacutvals = davals[np.logical_not(degcut)]
            numarr = []
            for i in range(len(eqda2vals)-1):
                tmask = np.logical_and(np.less(davals,eqda2vals[i+1]),np.greater_equal(davals,eqda2vals[i]))
                numarr.append(len(davals[tmask]))
            badstats = np.where(np.logical_and(np.less(numarr,prebin),np.greater(numarr,0)))[0]
            
            tmask = np.ones(shape=eqda2vals.shape,dtype=np.bool_)
            if len(badstats) != 0:
                for i in range(len(badstats)):
                    if badstats[i]+1 in badstats:
                        tmask[badstats[i]] = False
                    else:
                        tmask[badstats[i]+1] = False 
            eqda2vals = eqda2vals[tmask]
            
            eqdavals = np.append(eqda1vals,eqda2vals)
            
        else:
            if prebin <= 0:
                eqdavals = np.logspace(np.log10(cutval),np.log10(endval),num=postbin)
            elif postbin <= 0:
                eqdavals = np.logspace(np.log10(begval),np.log10(cutval),num=prebin)
            else:   
                eqda1vals = np.logspace(np.log10(begval),np.log10(cutval),num=prebin,endpoint=False)
                eqda2vals = np.logspace(np.log10(cutval),np.log10(endval),num=postbin)
                eqdavals = np.append(eqda1vals,eqda2vals)

    avSFvals,avSFEvals = mp_bootstrap(SFvals,davals,eqdavals)
    avSFNoisevals,avPNTvals,avSFNEvals = AverageSF(davals,eqdavals,SFEvals,davals)

    eqdawidth = np.ediff1d(eqdavals)/2.0
    eqdavals = eqdavals[:-1] + eqdawidth
    
    avSFvals = avSFvals - avSFNoisevals
    return avSFvals,avSFEvals,eqdavals,eqdawidth,avSFNoisevals,avPNTvals      

def ReadFolder(folname):
    folname = glob(folname)
    
    t1,t2 = readfile(folname[0])
    for tf in folname[1:]:
        #if tf.split('\\')[-1].split('.')[0] in ['src1184','src1195']:
        #    continue 
        th1,th2 = readfile(tf)
        t1 = np.append(t1,th1)
        t2 = np.append(t2,th2)

    return t1,t2


def RaDecStr(radec):
    rh,rm,rs = radec.ra.hms
    dd,dm,ds = radec.dec.dms    
    rs = sigfig(rs,1e-2)
    ds = sigfig(ds,1e-2)
    if rh < 10:
        rh = "0"+str(int(rh))
    else:
        rh = str(int(rh))
    if rm < 10:
        rm = "0"+str(int(rm))
    else:
        rm = str(int(rm))
    if rs < 10:
        rs = "0"+str(rs)
    else:
        rs = str(rs)   
    if abs(dm) < 10:
        dm = "0"+str(int(abs(dm)))
    else:
        dm = str(int(abs(dm)))
    if abs(ds) < 10:
        ds = "0"+str(abs(ds))
    else:
        ds = str(abs(ds))
    tstr = rh+":"+rm+":"+rs+" "+str(int(dd))+":"+dm+":"+ds
    return tstr

def PositionFilter(rcrd,rcrd2,ratol=0.9,dectol=6.0,selfcheck=False):
    rcrd = SkyCoord(rcrd,unit=(u.deg,u.deg),frame='fk5')
    rcrd2 = SkyCoord(rcrd2,unit=(u.deg,u.deg),frame='fk5')
    newrcrd = []
    newrcrd2 = []
    amask = np.ones(shape=(len(rcrd2),),dtype=np.bool_)
    for i in range(len(rcrd2)):
        r2h,r2m,r2s = rcrd2[i].ra.hms
        d2d,d2m,d2s = rcrd2[i].dec.dms
        if selfcheck:
            tmask = np.copy(amask)
            tmask[i] = False
            r1h,r1m,r1s = rcrd.ra[tmask].hms
            d1d,d1m,d1s = rcrd.dec[tmask].dms            
        else:
            r1h,r1m,r1s = rcrd.ra.hms
            d1d,d1m,d1s = rcrd.dec.dms
        
        rhcheck = np.equal(r1h,r2h)
        rmcheck = np.equal(r1m,r2m)
        rscheck = np.less_equal(np.abs(r1s-r2s),ratol)
        rhmscheck = np.logical_and(np.logical_and(rhcheck,rmcheck),rscheck)
        ddcheck = np.equal(d1d,d2d)
        dmcheck = np.equal(d1m,d2m)
        dscheck = np.less_equal(np.abs(d1s-d2s),dectol)
        ddmscheck = np.logical_and(np.logical_and(ddcheck,dmcheck),dscheck)
        check = np.logical_and(rhmscheck,ddmscheck)
        if len(np.where(check)[0]) > 0:
            amask[i] = False
    return amask

def WriteToFile(fname,arr1,arr2):
    f = open(fname,'w')
    for i in range(len(arr1)):
        f.write('%d %d \n'%(arr1[-(i+1)],arr2[-(i+1)]))
    f.close()
    return


