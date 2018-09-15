from RM_toolbox import *


linfit = lambda x,a,b: a*x+b
expfit = lambda x,a,b: b*np.exp(a*x)
powfit = lambda x,a,b: b*x**a
pow2fit = lambda x,a,b,c: b + np.log(x**a + c/b)
powp2fit = lambda x,a,b,c: b*x**a + c
logfit = lambda x,a,b: a*np.log(x) + b

hashfiles = "C:/Users/smile_000/Documents/Py_Docs/testsubsub/*Hash.txt"
orshfiles = "C:/Users/smile_000/Documents/Py_Docs/testsubsub/*ORsh.txt"

corfiles = "C:/Users/smile_000/Documents/Py_Docs/testhaorc/*.HaOR.txt"

#a1,sha2 = readfile("C:/Users/smile_000/Documents/Py_Docs/testhaorc/orvsshassapix.txt")
#a1,sha2 = ReadFolder("C:/Users/smile_000/Documents/Py_Docs/testsubsub/*Hash.txt")
#sha2,a1 = readfile("C:/Users/smile_000/Documents/Py_Docs/testhaorc/src1247.sumsha.txt")
#a1,sha2 = readfile("C:/Users/smile_000/Documents/Py_Docs/testhaorc/testline.txt")

cha,cor = ReadFolder(corfiles)

poly,covy = np.polyfit(cor,cha,deg=1,cov=True)
polyerr = np.sqrt(np.abs(covy))
pchaa,pchab = poly[0],poly[1]

print pchaa,'+/-',polyerr[0,0]
print pchab,'+/-',polyerr[1,1]
linarr = np.linspace(np.amin(cor),np.amax(cor),100)

poly,covy = np.polyfit(np.log(cor),np.log(cha),deg=1,cov=True)
polyerr = np.sqrt(np.abs(covy))
ppowa,ppowb = poly[0],np.exp(poly[1])

print ppowa,polyerr[0,0]
print ppowb,ppowb*polyerr[1,1]

plt.figure()
plt.plot(linarr,linfit(linarr,pchaa,pchab),label='$y_{lin} = %.2E x+ %.2E$'%(pchaa,pchab))
plt.plot(linarr,powfit(linarr,ppowa,ppowb),label='$y_{pow} = %.2E x^{%.2E}$'%(ppowb,ppowa))
plt.scatter(cor,cha,s=2,color='k')
plt.ylabel('$H \\alpha $ data  [$R$]',fontsize=18)
plt.xlabel('$OR $ data [$R$]',fontsize=18)
plt.title("Continuum Matching$")
plt.tight_layout()
#plt.xscale('log')
#plt.yscale('log')
plt.legend(loc=2)
plt.show()

ha1,sha1 = ReadFolder(hashfiles)
or1,sor1 = ReadFolder(orshfiles)
#WriteToFile("C:/Users/smile_000/Documents/Py_Docs/testhaorc/hashall.txt",ha1,sha1)
#WriteToFile("C:/Users/smile_000/Documents/Py_Docs/testhaorc/orshall.txt",or1,sor1)

shasort=np.argsort(sha1)
sha1= sha1[shasort]
ha1 = ha1[shasort]

shasort=np.argsort(sor1)
sor1= sor1[shasort]
or1 = or1[shasort]

sha1 = sha1*1e-1
sor1 = sor1*1e-1
#or1 = or1/1.0e10
#ha1 = ha1/1.0e10

cut1 = np.less(sha1,400)
cut2 = np.less(sha1,10000)

#ha1 = -2.5*np.log10(ha1)
#sha1 = -2.5*np.log10(sha1)

#linear no bs
poly,covy = curve_fit(linfit,sha1,ha1)
polyerr = np.sqrt(np.abs(covy))
plina,plinb = poly[0],poly[1]

#linearcut
poly,covy = curve_fit(linfit,sha1[cut1],ha1[cut1])
polyerr = np.sqrt(np.abs(covy))
plinlrah,plinlrbh = poly[0],poly[1]

#exp all
#poly,covy = curve_fit(linfit,sha1,np.log(ha1))
poly,covy = np.polyfit(sha1,np.log(ha1),deg=1,cov=True)
polyerr = np.sqrt(np.abs(covy))
pexpa,pexpb = poly[0],np.exp(poly[1])

#pow all
#poly,covy = curve_fit(linfit,np.log(sha1),np.log(ha1))
poly,covy = np.polyfit(np.log(sha1),np.log(ha1),deg=1,cov=True)
polyerr = np.sqrt(np.abs(covy))
ppowhaa,ppowhab = poly[0],np.exp(poly[1])

print ppowhaa, polyerr[0,0]
print ppowhab, ppowhab*polyerr[1,1]

#log all
poly,covy = np.polyfit(np.log(sha1),ha1,deg=1,cov=True)
polyerr = np.sqrt(np.abs(covy))
ploga,plogb = poly[0],poly[1]

sha1arrhr = np.linspace(np.amin(sha1),np.amax(sha1),100)
sha1arrlr = np.linspace(np.amin(sha1),np.amax(sha1[cut1]),100)

plt.figure()
#plt.plot(sha1arrhr,linfit(sha1arrhr,plina,plinb),label='$y_{lin} = %.2E x+ %.2E$'%(plina,plinb))
#plt.plot(sha1arrhr,linfit(sha1arrhr,plinlrah,plinlrbh),label='$y_{lin \ low \ R} = %.2E x+ %.2E$'%(plinlrah,plinlrbh))
plt.plot(sha1arrhr,expfit(sha1arrhr,pexpa,pexpb),label='$y_{exp} = %.2E e^{%.2E x}$'%(pexpb,pexpa),color='red')
plt.plot(sha1arrhr,powfit(sha1arrhr,ppowa,ppowb),label='$y_{pow} = %.2E x^{%.2E}$'%(ppowb,ppowa),color='cyan')
plt.plot(sha1arrhr,logfit(sha1arrhr,ploga,plogb),label='$y_{ln} = %.2E ln(x) + %.2E $'%(ploga,plogb),color='magenta')
plt.scatter(sha1,ha1,s=2,color='k')
plt.ylabel('$LMC$ $H \\alpha $ data  [$counts/pixel$]',fontsize=18)
plt.xlabel('$SHASSA$ $H \\alpha $ data [$R$]',fontsize=18)
plt.xlim([np.amin(sha1)-10,np.amax(sha1)+1000])
plt.title("Fits for $H \\alpha$")
plt.tight_layout()
plt.legend(loc=2)
#plt.xscale('log')
#plt.yscale('log')
plt.show()

plt.figure()
plt.hist(100*(expfit(sha1,pexpa,pexpb)-ha1)/ha1,30,histtype='step',label='$y_{exp} = %.2E e^{%.2E x}$'%(pexpb,pexpa))
plt.hist(100*(powfit(sha1,ppowa,ppowb)-ha1)/ha1,30,histtype='step',label='$y_{pow} = %.2E x^{%.2E}$'%(ppowb,ppowa))
plt.hist(100*(logfit(sha1,ploga,plogb)-ha1)/ha1,30,histtype='step',label='$y_{ln} = %.2E ln(x) + %.2E $'%(ploga,plogb))
plt.xlabel('$H \\alpha $ Percent Residuals [%]',fontsize=18)
plt.ylabel('Counts [arb. units]',fontsize=18)
#plt.xlim([np.amin(sha1)-10,np.amax(sha1)+1000])
plt.title("Hist for $H \\alpha$")
plt.tight_layout()
plt.legend(loc=0)
plt.show()


cut1 = np.less(sor1,400)
#or1 = -2.5*np.log10(or1)
#sor1 = -2.5*np.log10(sor1)

#linear no bs
poly,covy = curve_fit(linfit,sor1,or1)
polyerr = np.sqrt(np.abs(covy))
plina,plinb = poly[0],poly[1]

#linearcut
poly,covy = curve_fit(linfit,sor1[cut1],or1[cut1])
polyerr = np.sqrt(np.abs(covy))
plinlra,plinlrb = poly[0],poly[1]

#exp all
poly,covy = np.polyfit(sor1,np.log(or1),deg=1,cov=True)
polyerr = np.sqrt(np.abs(covy))
pexpa,pexpb = poly[0],np.exp(poly[1])

#pow all
poly,covy = np.polyfit(np.log(sor1),np.log(or1),deg=1,cov=True)
polyerr = np.sqrt(np.abs(covy))
ppowora,ppoworb = poly[0],np.exp(poly[1])

print ppowora, polyerr[0,0]
print ppoworb, ppoworb*polyerr[1,1]

#log all
poly,covy = np.polyfit(np.log(sor1),or1,deg=1,cov=True)
polyerr = np.sqrt(np.abs(covy))
ploga,plogb = poly[0],poly[1]

sor1arrhr = np.linspace(0.01,np.amax(sor1),100000)
sor1arrlr = np.linspace(np.amin(sor1),np.amax(sor1[cut1]),100)


plt.figure()
#plt.plot(sor1arrhr,linfit(sor1arrhr,plina,plinb),label='$y_{lin} = %.2E x+ %.2E$'%(plina,plinb))
#plt.plot(sor1arrhr,linfit(sor1arrhr,plinlra,plinlrb),label='$y_{lin \ low \ R} = %.2E x+ %.2E$'%(plinlra,plinlrb))
#plt.plot(sha1arrhr,linfit(sha1arrhr,plinlrah,plinlrbh),label='$g_{lin \ low \ R} = %.2E x+ %.2E$'%(plinlrah,plinlrbh))
plt.plot(sor1arrhr,expfit(sor1arrhr,pexpa,pexpb),label='$y_{exp} = %.2E e^{%.2E x}$'%(pexpb,pexpa),color='red')
plt.plot(sor1arrhr,powfit(sor1arrhr,ppowa,ppowb),label='$y_{pow} = %.2E x^{%.2E}$'%(ppowb,ppowa),color='cyan')
plt.plot(sor1arrhr,logfit(sor1arrhr,ploga,plogb),label='$y_{ln} = %.2E ln(x) + %.2E$'%(ppowb,ppowa),color='magenta')
plt.scatter(sor1,or1,s=2,color='k')
plt.ylabel('$LMC$ OR data  [$counts/pixel$]',fontsize=18)
plt.xlabel('$SHASSA$ OR data [$R$]',fontsize=18)
plt.xlim([np.amin(sor1)-10,np.amax(sor1)+1000])
plt.title("Fits for OR ")
plt.legend(loc=0)
plt.tight_layout()
plt.show()

plt.figure()
#plt.hist(100*(linfit(sor1,plinlra,plinlrb)-or1)/or1,30,histtype='step',label='$y_{lin} = %.2E e^{%.2E x}$'%(pexpb,pexpa),color='g')
plt.hist(100*(expfit(sor1,pexpa,pexpb)-or1)/or1,30,histtype='step',label='$y_{exp} = %.2E e^{%.2E x}$'%(pexpb,pexpa),color='r')
plt.hist(100*(powfit(sor1,ppowa,ppowb)-or1)/or1,30,histtype='step',label='$y_{pow} = %.2E x^{%.2E}$'%(ppowb,ppowa),color='cyan')
plt.hist(100*(logfit(sor1,ploga,plogb)-or1)/or1,30,histtype='step',label='$y_{ln} = %.2E ln(x) + %.2E $'%(ploga,plogb),color='magenta')
plt.xlabel('$OR $ Percent Residuals [%]',fontsize=18)
plt.ylabel('Counts [arb. units]',fontsize=18)
#plt.xlim([np.amin(sha1)-10,np.amax(sha1)+1000])
plt.title("Hist for $ OR $")
plt.tight_layout()
plt.legend(loc=0)
plt.show()


