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
    active = tb.Int32Col(pos=2)
    x = tb.Float32Col(pos=3)
    y = tb.Float32Col(pos=4)
    coeff = tb.Float64Col(pos=5)
    adc_to_pes = tb.Float32Col(pos=6)
    noise_rms = tb.Float32Col(pos=7)


def loadPMTs(h5f, cursor):
    group = h5f.create_group(h5f.root, "Sensors")
    pmt_table = h5f.create_table(group, "DataPMT", DataSensor,
                                 "DataPMT", tb.Filters(0))

    cursor.execute("SELECT * FROM Sensors")

    for row in cursor.fetchall():
        pmt = pmt_table.row
        pmt['sensorID'] = row[0]
        pmt['channel'] = row[1]
        pmt['active'] = row[2]
        pmt['noise_rms'] = row[3]
        pmt['coeff'] = row[4]
        pmt['adc_to_pes'] = row[5]
        pmt['x'] = row[6]
        pmt['y'] = row[7]
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
