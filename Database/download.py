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
    adc_to_pes = tb.Float64Col(pos=8)
    noise_rms = tb.Float64Col(pos=9)


class DataSensorSipm(tb.IsDescription):
    """
    Stores metadata information for the SiPMs
    (position, gain, calibration-constant, mask)
    """
    sensorID = tb.Int32Col(pos=0)
    active = tb.Int32Col(pos=1)
    x = tb.Float32Col(pos=2)
    y = tb.Float32Col(pos=3)
    adc_to_pes = tb.Float64Col(pos=4)
    noise_rms = tb.Float64Col(pos=5)


class DetectorGeo(tb.IsDescription):
    """
    Stores metadata information for the SiPMs
    (position, gain, calibration-constant, mask)
    """
    xmin = tb.Float32Col(pos=0)
    xmax = tb.Float32Col(pos=1)
    ymin = tb.Float32Col(pos=2)
    ymax = tb.Float32Col(pos=3)
    zmin = tb.Float32Col(pos=4)
    zmax = tb.Float32Col(pos=5)
    rmax = tb.Float32Col(pos=6)


def loadPMTs(h5f, group, cursor):
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
        pmt['adc_to_pes'] = row[8]
        pmt['noise_rms'] = row[9]
        pmt.append()

    pmt_table.flush()


def loadSiPMs(h5f, group, cursor):
    sipm_table = h5f.create_table(group, "DataSiPM", DataSensorSipm,
                                 "DataSiPM", tb.Filters(0))

    cursor.execute("SELECT * FROM Sipms")

    for row in cursor.fetchall():
        sipm = sipm_table.row
        sipm['sensorID'] = row[0]
        sipm['active'] = row[1]
        sipm['x'] = row[2]
        sipm['y'] = row[3]
        sipm['adc_to_pes'] = row[4]
        sipm['noise_rms'] = row[5]
        sipm.append()

    sipm_table.flush()


def loadDetGeo(h5f, group, cursor):
    geo_table = h5f.create_table(group, "DetectorGeometry", DetectorGeo,
                                 "DetectorGeometry", tb.Filters(0))

    cursor.execute("SELECT * FROM DetectorGeo")

    for row in cursor.fetchall():
        geo = geo_table.row
        geo['xmin'] = row[0]
        geo['xmax'] = row[1]
        geo['ymin'] = row[2]
        geo['ymax'] = row[3]
        geo['zmin'] = row[4]
        geo['zmax'] = row[5]
        geo['rmax'] = row[6]
        geo.append()

    geo_table.flush()


def loadDB():
    fname = os.environ['ICDIR'] + '/Database/localdb.h5'
    h5out = tb.open_file(fname, 'w')
    db = MySQLdb.connect(host="neutrinos1.ific.uv.es", user=dec('am1iZW5sbG9jaA=='),
                         passwd=eval(dec('Jycuam9pbihtYXAobGFtYmRhIGM6IGNocihjLTUpLCBbNzIsIDEwMiwgMTE1LCAxMDcsIDExOSwgMTAyLCAxMTUsIDEwNF0pKQ==')), db="ICNEWDB")
    cur = db.cursor()

    gSensors = h5out.create_group(h5out.root, "Sensors")
    loadPMTs(h5out, gSensors, cur)
    loadSiPMs(h5out, gSensors, cur)
    gGeo = h5out.create_group(h5out.root, "Detector")
    loadDetGeo(h5out, gGeo, cur)

    db.close()
    h5out.close()


if __name__ == '__main__':
    loadDB()
