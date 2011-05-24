""" Module for structural descriptor calculation
	author: Tobias Girschick; tobias.girschick@in.tum.de
			TUM - I12 (wwwkramer.in.tum.de/girschic)
	dependencies: FreeTreeMiner http://wwwkramer.in.tum.de/research/data_mining/pattern_mining/graph_mining
	to come: gSpan', BBRC and LAST-PM 
"""

import os
import sys
import time
import random
import string

import subprocess
import tempfile
import orange
from cinfony import rdk
import AZOrangeConfig as AZOC

#FTM = os.path.join(os.environ["HOME"], "OpenTox", "serverFiles", "FTM", "ftm")
#FTM = os.path.join("/home/kgvf414/projects/SimBoostedQSAR/FTM", "ftm")
FTM = os.path.join("", "ftm")

def getFTMDescResult(data,minSup):
	""" Calculates the structural FTM descriptors for the data using FTM
		It expects a relative minimum frequency parameter and a data attribute containing smiles with a name defined in AZOrangeConfig.SMILESNAMES
		It returns a dataset with the same smiles input variable, and as many variables as the descriptors returned by the toolkit
	"""
	sdf_mols, temp_occ = makeTempFiles(data)
		
	# start FTM subprocess
	ftm(temp_occ.name, minSup, sdf_mols.name)
	#ftm(temp_occ.name, minSup, 'testftm.sdf')
	
	# parse result file and create new data
	occ_flag = 0
	atts = []
	occ_list = []
	for line in temp_occ:

		if occ_flag == 1:
			line_split = line.split()
			#
			# line_split[0] is the "name" of the feature
			# line_split[1] is the occurrence string for all instances separated by comma
			#

			# add new attributes to list
			atts.append(orange.FloatVariable(line_split[0], numberOfDecimals=1))
			#atts.append(orange.EnumVariable(line_split[0], values=['0','1']))
									
			occ = line_split[1].split(",")
			occ_list.append(occ)
					
		if line.startswith('Occurrences:'):
			occ_flag = 1
	
	newdomain = orange.Domain(data.domain.attributes + atts, data.domain.classVar)
	newdata = orange.ExampleTable(newdomain, data)

	count1 = 1
	for j in range(len(newdata)):
		instance = newdata[j]
		count = 1	
				
		for k in range(len(occ_list)):
			if str(occ_list[k][j]) == '1':
				tmp = orange.Value(atts[k],1.0)
				instance[atts[k]] = tmp
			elif str(occ_list[k][j]) == '0':
				tmp = orange.Value(atts[k],0.0)
				instance[atts[k]] = tmp
							
			count += 1
		
	temp_occ.close()
	sdf_mols.close()
	return newdata
	

def makeTempFiles(data):
	"""	create temporary files for FTM I/O
		returns to file objects that still have to be closed!
		temp_occ contains the occurrence list of substructures (result)
		sdf_mols is for the FTM input (SDF structure file)
	"""
	
	smilesName = getSMILESAttr(data)
	if not smilesName: return None
	
	temp_occ = tempfile.NamedTemporaryFile()
	sdf_mols = tempfile.NamedTemporaryFile(suffix='.sdf') 
    
    # write structures in sdf format to this file
	count = 0
	w = rdk.Chem.SDWriter(sdf_mols.name)
	#w = rdk.Chem.SDWriter('testftm.sdf')
	count_mols = 0
	for a in data:
		smile = str(a[smilesName].value)
		m = rdk.Chem.MolFromSmiles(smile)
		if m is None: 
			count += 1
			#print "ALERT"
			continue
		# set smiles as molname (just to have any name)	
		m.SetProp("_Name", "-->"+smile)
		w.write(m)

		count_mols+=1
	print "Number of molecules that could not be read: ", count
	print "Number of molecules read: ", count_mols

	return sdf_mols, temp_occ
	
def getSMARTSrecalcDesc(data,smarts):
	""" Calculates structural descriptors for test and training data.
		In other words, checks for the substructure occurrence (0/1) in the 
		test or prediction molecules. Uses RDK.
		Expects the test/prediction data and a list of SMARTS strings.
		Returns the data including the new features. 
	"""
	
	smilesName = getSMILESAttr(data)
	if not smilesName: return None
    		
	atts = []
	for attr in smarts:
		#atts.append(orange.EnumVariable(attr.name, values=['0','1']))
		atts.append(orange.FloatVariable(attr.name, numberOfDecimals=1))
		
	newdomain = orange.Domain(data.domain.attributes + atts, data.domain.classVar)
	newdata = orange.ExampleTable(newdomain, data)
		
	count = 0
	for a in newdata:
		smile = str(a[smilesName].value)
		m = rdk.Chem.MolFromSmiles(smile)
		if m is None: 
			count += 1
			continue
		
		for b in range(len(smarts)):
			#print b.name
			patt = rdk.Chem.MolFromSmarts(smarts[b].name)
			if m.HasSubstructMatch(patt):
				tmp = orange.Value(atts[b],1.0)
			else:
				tmp = orange.Value(atts[b],0.0)
			a[atts[b]] = tmp
				
	return newdata
	
	
def getSMILESAttr(data):
    # Check that the data contains a SMILES attribute
    smilesName = None
    for attr in [a.name for a in  data.domain] + [a.name for a in data.domain.getmetas().values()]:
        if attr in AZOC.SMILESNAMES:
            smilesName = attr
    if not smilesName:
        print "Warning: The data set does not contain any known smiles attribute!"
        print "No Cinfony descriptors added!"
        return None
    else:       
        return smilesName
	
	
def ftm(temp_occ, freq, mols):
	ftm_opt = ' -f ' + str(freq) + ' -o ' + temp_occ	+ ' ' + mols
	ftm_call = FTM + ftm_opt
	print ftm_call
	p = subprocess.call(ftm_call, shell=True)

