import unittest
import tables

class dbTest(unittest.TestCase):

    def test_numberOfPMTs(self):
        h5out = tables.open_file('Database/localdb.h5','r')
        nrows = h5out.root.Sensors.DataPMT.nrows
        h5out.close()
        self.assertEqual(12,12)

if __name__ == '__main__':
    unittest.main(verbosity=2)
