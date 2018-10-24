#!/usr/bin/python

import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time
import seaborn as sns
import os
from pathlib import Path
#import math


forces = [ 5.0 ,6.0 ,7.0 , 8.0, 9.0, 10.0, 11.0, 12.0]
delimiter = "   "


for i in range(0,8):
	f=forces[i]
	folder='../%s/*.dat'%(f)
	#average file
	averagefile='../%s/average.txt'%(f)
	var_averagefile='../%s/var_average.txt'%(f)
	my_file=Path(averagefile)
	if my_file.is_file():
		print ("Average file already created")
	
		#create an average file
		#create a matrix with all data points of all the files
		#number of columns = numfiles
		numfiles=0
		for filename in glob.glob(folder):
			numfiles+=1
		numfilestotal=numfiles
		#number of lines= lines of 1 file
		cont=1
		numfiles=0
		for filename in glob.glob(folder):
			g=open(filename,"r")
			lines=g.readlines()
			for y in lines:
				cont+=1
			g.close()
			numfiles+=1
			if numfiles>=1:
				break;
		numlines=cont
		print("Numfile: ",numfilestotal)
		print("Numlines: ", numlines)
		
		t=[0 for x in range(numlines)]
		x=[[0 for x in range(numlines)] for y in range(numfilestotal)]
	
		#print ("Average curve")
		print(averagefile)
		numfiles=0
		for filename in glob.glob(folder):	
			h=open(filename,"r")
			lines=h.readlines()
			cont=0
			for y in lines:
				line=y.rstrip()
				if numfiles==1:
					t[cont]+=float(line.split(delimiter)[0])
				x[numfiles][cont]+=float(line.split(delimiter)[1])
				cont+=1
			h.close()
			#adding x of each file at each line
			if numfiles==0:
				xav=np.zeros(cont+1)
			for i in range(0,cont+1):
				xav[i]+=x[numfiles][i]
			numfiles+=1
			#writing it to the average file
	#else:
		fout=open(averagefile,"w")
		for i in range(0,cont):
			xav[i]/=numfiles		
			fout.write('%.12f %.12f\n'%(t[i],xav[i]))
		fout.close()
	#calculation of the variance
	numfiles=0
	varianceav=[0 for x in range(numlines)]
	fout2=open(var_averagefile,"w")
	for i in range(0,numlines):
	#for i in range(0,9999):
		for  j in range(0,numfilestotal):
			varianceav[i]+=(x[j][i]-xav[i])**2
			#print(x[j][i])
			#print(xav[i])
			#print((x[j][i]-xav[i])**2)
		#print (t[i])
		#print (varianceav[i]/numfilestotal)
		fout2.write('%.12f %.12f \n'%(t[i],varianceav[i]/numfilestotal))
	fout2.close()
	
