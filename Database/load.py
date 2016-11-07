import tables
import numpy as np
import MySQLdb

class DataSensor(tables.IsDescription):
    """
    Stores metadata information for the SiPMs
    (position, gain, calibration-constant, mask)
    """
    sensorID = tables.Int32Col(pos=0)
    channel = tables.Int32Col(pos=1)
    active = tables.Int32Col(pos=2)
    position = tables.Float32Col(shape=2, pos=3)
    coeff = tables.Float64Col(pos=4)
    adc_to_pes = tables.Float32Col(pos=5)
    noise_rms = tables.Float32Col(pos=6)

def loadPMTs(h5f, cursor):
    group = h5f.create_group(h5f.root, "Sensors")
    pmt_table = h5f.create_table(group, "DataPMT", DataSensor,"DataPMT",tables.Filters(0))

    cursor.execute("SELECT * FROM Sensors")

    for row in cursor.fetchall():
        pmt = pmt_table.row
        pmt['sensorID'] = row[0]
        pmt['channel'] = row[1]
        pmt['active'] = row[2]
        pmt['noise_rms'] = row[3]
        pmt['coeff'] = row[4]
        pmt['adc_to_pes'] = row[5] 
        pmt['position'] = np.array((row[6],row[7]))
        pmt.append()

    pmt_table.flush()


def loadDB():
    h5out = tables.open_file('Database/localdb.h5','w')
    db = MySQLdb.connect(host="localhost", user="next", passwd="Canfranc", db="ICNEW")
    cur = db.cursor()
    
    loadPMTs(h5out, cur)
    
    db.close()
    h5out.close()


if __name__ == '__main__':
    loadDB()
