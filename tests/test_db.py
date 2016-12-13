import unittest
import tables
from Database import download
import Database.loadDB as DB
import os


class dbTest(unittest.TestCase):

    def setUp(self):
        download.loadDB()
        self.localdb = os.environ['ICDIR'] + '/Database/localdb.h5'

    def test_numberOfPMTs(self):
        """
        Check that we retrieve the correct number of PMTs
        """
        pmts = DB.DataPMT()
        self.assertEqual(pmts.shape[0], 12)

    def test_numberOfSiPMs(self):
        """
        Check that we retrieve the correct number of SiPMs
        """
        sipms = DB.DataSiPM()
        self.assertEqual(sipms.shape[0], 1792)

    def test_SiPMNoise(self):
        """
        Check we have noise for all SiPMs and energy of each bin
        """
        noise, energy, baseline = DB.SiPMNoise()
        self.assertEqual(noise.shape[0], baseline.shape[0])
        self.assertEqual(noise.shape[0], 1792)
        self.assertEqual(noise.shape[1], energy.shape[0])

    def test_DetectorGeometry(self):
        """
        Check Detector Geometry
        """
        geo = DB.DetectorGeo()
        self.assertEqual(geo['XMIN'][0], -198)
        self.assertEqual(geo['XMAX'][0], 198)
        self.assertEqual(geo['YMIN'][0], -198)
        self.assertEqual(geo['YMAX'][0], 198)
        self.assertEqual(geo['ZMIN'][0], 0)
        self.assertEqual(geo['ZMAX'][0], 532)
        self.assertEqual(geo['RMAX'][0], 198)

if __name__ == '__main__':
    unittest.main(verbosity=2)
