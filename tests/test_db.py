import unittest
import tables
from Database import download


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

if __name__ == '__main__':
    unittest.main(verbosity=2)
