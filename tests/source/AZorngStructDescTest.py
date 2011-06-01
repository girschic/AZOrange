import unittest
import os
import time

from AZutilities import dataUtilities
from AZutilities import getCinfonyDesc
import AZOrangeConfig as AZOC

class evalUtilitiesTest(unittest.TestCase):

    def setUp(self):
        
        smiDataPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/mol.smi")
        self.smiData = dataUtilities.loadSMI(smiDataPath)
        
	def test_makeTmpFile(self):
		
		
	def test_ftm(self):
		
		
	def test_smartsRecalc(self):	
		

if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(evalUtilitiesTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
    #unittest.main()
