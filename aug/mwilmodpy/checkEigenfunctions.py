
import sys,os
sys.path.append("/afs/ipp-garching.mpg.de/home/m/mwillens/Programme/python/lib/")

import numpy as np
from IPython import embed
import matplotlib.pylab as plt
import commands

def readAllEigenfunctions(ending='_n12_axi_90',pathEig='../out/',nEig=1):
	Eigenfunctions = []
	radEigRe, radEigIm, radVelRe, radVelIm, radTorRe, radTorIm = [], [],[], [],[], []
	
	for e in np.arange(nEig):
		radEigReTmp, radEigImTmp, radVelReTmp, radVelImTmp, radTorReTmp, radTorImTmp = [], [],[], [],[], []
		for c in np.arange(3)+1:	
			radEigReTmp.append(readSingleEigenfunctions(fileIn=pathEig+'rrev'+str(c)+str(e+1).zfill(3)+ending))
			radEigImTmp.append(readSingleEigenfunctions(fileIn=pathEig+'rimv'+str(c)+str(e+1).zfill(3)+ending))
			radVelReTmp.append(readSingleEigenfunctions(fileIn=pathEig+'rrea'+str(c)+str(e+1).zfill(3)+ending))
			radVelImTmp.append(readSingleEigenfunctions(fileIn=pathEig+'rima'+str(c)+str(e+1).zfill(3)+ending))
		
		radEigRe.append(radEigReTmp)
		radEigIm.append(radEigImTmp)
		radVelRe.append(radVelReTmp)
		radVelIm.append(radVelImTmp)

		radTorRe.append(readSingleEigenfunctions(fileIn=pathEig+'rret'+str(e+1).zfill(3)+ending))
		radTorIm.append(readSingleEigenfunctions(fileIn=pathEig+'rimt'+str(e+1).zfill(3)+ending))
		
	embed()		
	rho,eigFun,velFun,torFun,mArr,nArr=[],[],[],[],[],[]
	
	for e in np.arange(nEig):
		for c in np.arange(3):
			for i in np.arange(len(radEigRe[0][0])):
				eigFun.append(radEigRe[e][c][i]['eigfun']+np.complex('j')*radEigIm[e][c][i]['eigfun'])
				velFun.append(radVelRe[e][c][i]['eigfun']+np.complex('j')*radVelIm[e][c][i]['eigfun'])
				if c==0:	
					torFun.append(radTorRe[e][i]['eigfun']+np.complex('j')*radTorIm[e][i]['eigfun'])
				
					rho.append(radEigRe[e][c][i]['rho'])
					mArr.append(radEigRe[e][c][i]['m'])
					nArr.append(radEigRe[e][c][i]['n'])
		
	rho,mArr,nArr=np.array(rho).mean(axis=0),np.array(mArr),np.array(nArr)
	newShape = (nEig,3,len(radEigRe[0][0]),rho.size)
	newShape2 = (nEig,1,len(radEigRe[0][0]),rho.size)
	eigFun,velFun,torFun = np.reshape(np.ravel(eigFun),newShape),np.reshape(np.ravel(velFun),newShape),np.reshape(np.ravel(torFun),newShape2)
	
	return rho,mArr,nArr,eigFun,velFun,torFun
	
	#plt.plot(mArr,np.abs(eigfun).max(axis=0))
	#plt.show()
	#embed()


def readEigenValues():
	print 'todo'

def checkAllnTor(end='_axi',pathEig='/draco/u/mwillens/CASTOR3D/AUGD/34424/micdu_eqb_5_ja/a5_j4/m6fAxiEven/out/',fileName='rrev1001_',plot=True):
	
	status, outputList = commands.getstatusoutput("ls "+pathEig+fileName+"n*"+end)
	fileNameList = np.array(outputList.split('\n'))
	nEigen = []
	
	for f in fileNameList:
		fIdx = f.find('_n')
		sIdx = f[fIdx+1:].find('_')
		nfile =  int(f[fIdx+2:fIdx+1+sIdx]) 
		rho,eigFun,mArr,nArr = readEigenfunctions(ending=f[fIdx:],pathEig=pathEig)
		eigArr=np.abs(eigFun).max(axis=0)
		#eigArr_low = np.exp(np.log(eigArr).max()-np.log(eigArr).std()*1.5)
		eigArr_low = np.exp(-12)
		idxValid = np.where( eigArr >=eigArr_low)[0]
		mLow=mArr[idxValid[0]]-5
		mHigh=mArr[idxValid[-1]]+8
		nEigen.append({'n':nfile,'rho':rho,'eigFun':eigFun,'mArr':mArr,'nArr':nArr,'mLow':mLow,'mHigh':mHigh})
	
	if plot:	
		co=['m','k','r','b','g','y']
		for i in np.arange(len(nEigen)):
			plt.plot(nEigen[i]['mArr'],np.abs(nEigen[i]['eigFun']).max(axis=0),label='n=%d'%nEigen[i]['n'],color=co[i])
			plt.axvline(nEigen[i]['mHigh'],color=co[i])
			plt.axvline(nEigen[i]['mLow'],color=co[i])
			print nEigen[i]['n'],nEigen[i]['mLow'],nEigen[i]['mHigh']
		plt.legend()	
		plt.show()
	
		embed()
	
	

def readEigenfunctions(ending='_n12_axi_90',pathEig='../out/',nEig=1,useEig=None):
	Eigenfunctions = []

	if useEig == None:	
		for e in np.arange(nEig):
			try:
				radEigRe = readEigenfunctionsFiles(fileIn=pathEig+'rrev1'+str(e+1).zfill(3)+ending)
				radEigIm = readEigenfunctionsFiles(fileIn=pathEig+'rimv1'+str(e+1).zfill(3)+ending)
			except:
				print 'file not found'
	else:
		radEigRe = readEigenfunctionsFiles(fileIn=pathEig+'rrev1'+str(useEig).zfill(3)+ending)
		radEigIm = readEigenfunctionsFiles(fileIn=pathEig+'rimv1'+str(useEig).zfill(3)+ending)	


	rho,eigfun,mArr,nArr=[],[],[],[]
	
	for i in np.arange(len(radEigRe)):
		eigfun.append(radEigRe[i]['eigfun']+np.complex('j')*radEigIm[i]['eigfun'])
		rho.append(radEigRe[i]['rho'])
		mArr.append(radEigRe[i]['m'])
		nArr.append(radEigRe[i]['n'])
		
	rho,eigfun,mArr,nArr=np.array(rho).mean(axis=0),np.array(eigfun).T,np.array(mArr),np.array(nArr)
	 
	return rho,eigfun,mArr,nArr
	#plt.plot(mArr,np.abs(eigfun).max(axis=0))
	#plt.show()
	#embed()


def readEigenfunctionsFiles(fileIn='../out/rrev1001_n6_axi_90'):
	print 'reading '+fileIn
	f = open(fileIn, 'r')
	#eigen={'m':}
	dimens= np.array(f.readline().split()[-2:],dtype='int')
	eigenfunctions=[]
	for c in np.arange(dimens[0]):
		ind = np.array(f.readline().split()[-1:],dtype='int')
		modeNumbers = np.array(f.readline().split()[-2:],dtype='int')
		
		eig=[]
		line=[]
		condition=True
		while condition:
			line = f.readline().split()		
			if line[0] != "&&":
				eig.append(np.array(line,dtype='double'))	
			else:	
				condition=False		
			
		eigenfunc = np.array(eig).T
		eigenfunctions.append({'i':ind[0],'m':modeNumbers[0],'n':modeNumbers[1],'eigfun':eigenfunc[1],'rho':eigenfunc[0]})
		
	f.close()
	return eigenfunctions


#checkAllnTor(pathEig='/tokp/work/mwillens/VMEC/AUGD/34424/micdu_eqb_5_ja_trunc/a5_j4/castor3DF/out/',end='_3D_m28_317_vt')

#end='_n024_m20_3D_m30_317_m10_vt'

#testEigenfunctions()	
#readAllEigenfunctions()	
	

