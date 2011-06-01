import unittest
import os
import time

from AZutilities import dataUtilities
from AZutilities import getSimBoostedQSAR
import AZOrangeConfig as AZOC

class evalUtilitiesTest(unittest.TestCase):

    def setUp(self):
        
        smiDataPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/mol.smi")
        self.smiData = dataUtilities.loadSMI(smiDataPath)


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(evalUtilitiesTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
    #unittest.main()
