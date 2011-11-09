""" Module for the integration of the TUM structural clustering algorithm
	author: Tobias Girschick; tobias.girschick@in.tum.de
			TUM - I12 (wwwkramer.in.tum.de/girschic)
	dependencies: gSpan, java
		
	Please cite the following article if you use the structural clustering procedure or results produced with it in any publication:
	
@InProceedings{ seeland2010,
	title = "Online Structural Graph Clustering Using Frequent Subgraph Mining",
	booktitle = "Proceedings of the ECML/PKDD'10",
	series = "Lecture Notes in Computer Science",
	author = "M. Seeland and T. Girschick and F. Buchwald and S. Kramer",
	editor = "J.L. Balc{\'a}zar and F. Bonchi and A. Gionis and M. Sebag",
	pages = "213--228",
	volume = "6323",
	year = "2010",
	booktitle1 = "Machine Learning and Knowledge Discovery in Databases, European Conference, {ECML} {PKDD} 2010, Barcelona, Spain, September 20-24, 2010, Proceedings, Part {III}"
}

======
parameters:
1) Dataset
2) Pfad zu Verzeichnis, in dem sich der gSpan Ordner befindet
3) threshold
4) gibt an, ab welcher Clustergröße die Cluster als sdf gespeichert
werden (0, falls überhaupt keine sdf's erstellt werden sollen)
5) gibt an, ab welcher Clustergröße die Cluster im output_cluster File
gespeichert werden
6) Output Pfad
7) Anzahl Threads
8) Timeout für gSpan
"""

import os
import sys
import time

import subprocess
from subprocess import Popen, PIPE
import tempfile
from AZutilities import getStructuralDesc


def getStructuralClusters(data, gspanpath = ".", threshold, minSaveSDFsize = 0, minClusterSize, numThreads=1, timeout=20):
	""" just the clustering
	"""
	sdf_temp = getStructuralDesc.makeTempSDF(data)
	# create tempdir for usage as 6) outputpath

	# call clustering routine
	# java -jar StructuralClustering.jar ../../SAR/631_active.sdf . 0.5 5 5 test/ 2 20

	# parse output 

	sdf_mols.close()
	return None


#def ftm(temp_occ, freq, mols):
#    ftm_opt = ' -f ' + str(freq) + ' -o ' + temp_occ	+ ' ' + mols
#    cmd = FTM + ftm_opt
    #print cmd
#    p = Popen(cmd, shell=True, close_fds=True, stdout=PIPE)
#    stdout = p.communicate()



def getReferenceStructures(data,):
	""" calls the clustering? 
		gets the clustering result?
		calls the different (now 1) methods to get representatives
	"""

def getRandomRepresentative():
	""" returns list of id\tsmiles pairs for usage in descriptor calculation
	"""


