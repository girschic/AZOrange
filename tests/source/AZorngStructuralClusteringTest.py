import unittest
import os

import orange
from AZutilities import structuralClustering
from AZutilities import dataUtilities
import AZOrangeConfig as AZOC

class StructuralClusteringTest(unittest.TestCase):

    def setUp(self):
        qsarDataPath = os.path.join(AZOC.AZORANGEHOME,"tests/source/data/AID677_100mols.tab")
        self.data = orange.ExampleTable(qsarDataPath)

    def test_clustering(self):
        """ Test if clustering is running properly """
	refs = structuralClustering.getReferenceStructures(data,threshold=0.5,minClusterSize=5)
        expected_length = 4
        self.assertEqual(expected_length, len(refs))
                
    def test_mkTMPfile(self):
        """Test if temporary SDF files can be written """
        import tempfile
        from subprocess import Popen, PIPE
        from cinfony import rdk
				
        sdf_mols = dataUtilities.makeTempFilesFTM(self.data)
        cmd = 'cat ' + sdf_mols.name + ' | grep \'\$\$\$\$\' | wc'
        p = Popen(cmd, shell=True, close_fds=True, stdout=PIPE)
        stdout = p.communicate()
        self.assertEqual(stdout[0].strip(),"100      100      500")
	
		
if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(StructuralClusteringTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
    #unittest.main()
