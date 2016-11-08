import tables as tb
import os
import MySQLdb
from base64 import b64decode as dec


class DataSensor(tb.IsDescription):
    """
    Stores metadata information for the SiPMs
    (position, gain, calibration-constant, mask)
    """
    sensorID = tb.Int32Col(pos=0)
    channel = tb.Int32Col(pos=1)
    pmtid = tb.StringCol(5,pos=2)
    active = tb.Int32Col(pos=3)
    x = tb.Float32Col(pos=4)
    y = tb.Float32Col(pos=5)
    deconvClean = tb.Float64Col(pos=6)
    deconvCoeff = tb.Float64Col(pos=7)
    noise_rms = tb.Float64Col(pos=8)


def loadPMTs(h5f, cursor):
    group = h5f.create_group(h5f.root, "Sensors")
    pmt_table = h5f.create_table(group, "DataPMT", DataSensor,
                                 "DataPMT", tb.Filters(0))

    cursor.execute("SELECT * FROM Sensors")

    for row in cursor.fetchall():
        pmt = pmt_table.row
        pmt['sensorID'] = row[0]
        pmt['channel'] = row[1]
        pmt['pmtid'] = row[2]
        pmt['active'] = row[3]
        pmt['x'] = row[4]
        pmt['y'] = row[5]
        pmt['deconvClean'] = row[6]
        pmt['deconvCoeff'] = row[7]
        pmt['noise_rms'] = row[8]
        pmt.append()

    pmt_table.flush()


def loadDB():
    fname = os.environ['ICDIR'] + '/Database/localdb.h5'
    h5out = tb.open_file(fname, 'w')
    db = MySQLdb.connect(host="neutrinos1.ific.uv.es", user=dec('am1iZW5sbG9jaA=='),
                         passwd=eval(dec('Jycuam9pbihtYXAobGFtYmRhIGM6IGNocihjLTUpLCBbNzIsIDEwMiwgMTE1LCAxMDcsIDExOSwgMTAyLCAxMTUsIDEwNF0pKQ==')), db="ICNEWDB")
    cur = db.cursor()
    loadPMTs(h5out, cur)

    db.close()
    h5out.close()


if __name__ == '__main__':
    loadDB()
