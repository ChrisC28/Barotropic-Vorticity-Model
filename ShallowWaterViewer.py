import numpy as np
from pylab import *

import subprocess                 # For issuing commands to the OS.
import os
import sys      

dataFileBase = "/home/cchlod/Barotropic_Vorticity_Model/BarotropicVorticityModel/output/streamFunctionOut_"


for i in range(0,49):

	inputFileName = dataFileBase + str(i) + ".dat"	
	appendFileName = '%02d.png'%i
	outputImageName = dataFileBase + appendFileName	
	print inputFileName
	eta_data = np.loadtxt(inputFileName, dtype = float, delimiter=',')

	CS = contour(eta_data, 20, linewidths=0.5,colors='k')
	CS = contourf(eta_data, 20,cmap=plt.cm.jet)
	xlim(0,100)
	ylim(0,50)
	print "Writing  "  + outputImageName
	savefig(outputImageName, dpi=100)
	clf()
	
systemString = "rm " + fileName
os.system("convert -delay 10 *.png anim.gif")













