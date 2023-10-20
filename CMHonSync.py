#! /usr/bin/env python
from numpy import *
from scipy import stats
import re
import sys
import os

MAF=0.05
COV=10
k=2 # Set that to the column that should be used to compute the distance, starting at 0.

if ".sync" not in sys.argv[1]:
	exit("The input file for CMH test is not in *.sync format, exiting.")

Pheno_File=sys.argv[1].split('.sync')[0]

if '/' in Pheno_File:
	Pheno_Dir=sys.argv[1].rsplit('/',1)[0]
	Pheno_File=sys.argv[1].rsplit('/',1)[1]
else:
	Pheno_Dir='.'

if not os.path.isfile(sys.argv[1]):
	exit("The input files for CMH test were not found, exiting.")

def colsplit(x):
	return(x.split(":"))

def calc_cmh_chisq(counts,cont_correct=.5):
	global nStacks
	numer_sum = 0
	denom_sum = 0

	for i in range(nStacks):
		a = counts[0,:,i]
		b = counts[1,:,i]
		c = counts[2,:,i]
		d = counts[3,:,i]
		n = a + b + c + d
		numer_sum += a - (a+b)*(a+c)/(1.0*n)
		denom_sum += (a+b)*(a+c)*(b+d)*(c+d)/(1.0*n*n*(n-1))

	chi_sq = ((abs(numer_sum) - cont_correct)*(abs(numer_sum) - cont_correct))/(1.0*denom_sum)
	return((chi_sq,1-stats.chi2.cdf(chi_sq,1)))

POP=genfromtxt(sys.argv[2], delimiter=",",missing_values='NA', dtype=str)[:,0]
FEAT=genfromtxt(sys.argv[2], delimiter=",",missing_values='NA', dtype=float)[:,1:]

Case=FEAT[:,k]>50 
Ctrl=FEAT[:,k]<50
POP_pair=array(meshgrid(nonzero(Case), nonzero(Ctrl))).T.reshape(-1,2)
nStacks=shape(POP_pair)[0]

n = shape(FEAT)[0]
DIST = zeros((n,n))

for i in range(n):
    for j in range(n):
	    DIST[i, j] = sqrt(sum((FEAT[:,k][i] - FEAT[:,k][j]) ** 2))

Freq_File=open(sys.argv[1],'rb').readlines()
SNP_call=['A','T','C','G','N','D']

CMH_out=[]

for SNP in Freq_File:
	SNP=SNP.decode("utf-8").replace('\n','').split('\t')
	SNP_name=[SNP[0],SNP[1]]
	SNP_K=array(list(map(colsplit,SNP[3:])),dtype=float)
	COVER=sum(SNP_K[:,:4],axis=1)
	freq=SNP_K[:,:4]/(COVER[:,None]+1e-6) # here I have to force non-zero coverage
	freq=array([[1-x if x>.5 else x for x in i] for i in freq])
	freq_max=amax(freq,axis=1)
	#here the rep stacking starts
	if any(freq_max>MAF) and all(COVER>(COV/n)):
		for i in range(4):
			Allele_freq=mean(freq[:,i])
			if len(CMH_out)>0 and round(Allele_freq,4)==round(CMH_line[5],4): continue
			elif Allele_freq>(MAF/n):
				Allele_count=array([]).reshape(4, 0)
				for x,y in POP_pair:
					Allele=vstack((SNP_K[x,i],SNP_K[y,i],sum(delete(SNP_K[x,:4], i, 0)),sum(delete(SNP_K[y,:4], i, 0))))
					Allele_count=hstack([Allele_count,Allele])
				
				Allele_count=Allele_count.reshape((4,1,nStacks))
				CMH=calc_cmh_chisq(Allele_count)
				CMH_line=[SNP_name[0],SNP_name[1],SNP_call[i],float(CMH[0]),float(CMH[1]),Allele_freq,COVER[i]]
				CMH_out.append(CMH_line)

header="Chromosome,Position,Mutation,X2_CMH,p_CMH,MAF,COV"
filename=Pheno_Dir+"/CMH_"+Pheno_File+"_out.csv"
savetxt(filename, array(CMH_out), delimiter=",",header=header,fmt="%s")
