from RM_toolbox import *

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

testvecfiles = glob("C:/Users/smile_000/Documents/Py_Docs/testhaorc/src1247.Pol*.v02.txt")
testvecv00files = glob("C:/Users/smile_000/Documents/Py_Docs/testhaorc/*Polvec.Ext.v02.txt")

thet_arr = np.linspace(0,360,10000)
xcos = np.cos(thet_arr)
ysin = np.sin(thet_arr)

v1,v2,wc1,n1= readvecfile(testvecfiles[0])
for tf in testvecfiles[1:]:
    tv1,tv2,twc1,tn1 = readvecfile(tf)
    v1 = np.append(v1,tv1,axis=0)
    v2 = np.append(v2,tv2,axis=0)
    wc1 = np.append(wc1,twc1)
    n1 = np.append(n1,tn1) 

for i in range(v1.shape[0]):
    plt.figure()
    pmax = (v1[i]**2 + v2[i]**2)**(0.5)
    pmax = pmax == np.amax(pmax)
    pval = (v1[i][pmax]**2 + v2[i][pmax]**2)**(0.5)
    plt.scatter(v1[i]/pval,v2[i]/pval,s=1.5,label='$F(\phi)$')
    plt.plot(n1[i]*xcos/pval,n1[i]*ysin/pval,color='red',label='$5\sigma$ Noise')
    plt.plot([0,v1[i][pmax]/pval],[0,v2[i][pmax]/pval],color='b',label='Delta Function')
    plt.xlabel('Q')
    plt.ylabel('U')
    plt.title('%s'%(wc1[i]))
    plt.xlim([-1,1])
    plt.ylim([-1,1])
    plt.axhline(y=0,ls='--',color='k')
    plt.axvline(x=0,ls='--',color='k')
    #plt.grid()
    plt.legend(loc=0,numpoints=1)
    plt.tight_layout()
    plt.show()

