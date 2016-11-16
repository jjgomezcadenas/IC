import unittest
import tables
from Database import download
from Database import loadDB


class dbTest(unittest.TestCase):

    def test_numberOfPMTs(self):
        """
        Check that we retrieve the correct number of PMTs
        """
        download.loadDB()
        h5out = tables.open_file('Database/localdb.h5', 'r')
        nrows = h5out.root.Sensors.DataPMT.nrows
        h5out.close()
        self.assertEqual(nrows, 12)

    def test_numberOfSiPMs(self):
        """
        Check that we retrieve the correct number of SiPMs
        """
        download.loadDB()
        h5out = tables.open_file('Database/localdb.h5', 'r')
        nrows = h5out.root.Sensors.DataSiPM.nrows
        h5out.close()
        self.assertEqual(nrows, 1792)

    def test_SiPMNoise(self):
        """
        Check we have noise for all SiPMs and energy of each bin
        """
        download.loadDB()
        noise, energy, baseline = loadDB.SiPMNoise()
        self.assertEqual(noise.shape[0], baseline.shape[0])
        self.assertEqual(noise.shape[0], 1792)
        self.assertEqual(noise.shape[1], energy.shape[0])

if __name__ == '__main__':
    unittest.main(verbosity=2)
