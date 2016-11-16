import tables as tb
import pandas as pd
import os


def DataPMT():
    fname = os.environ['ICDIR'] + '/Database/localdb.h5'
    with tb.open_file(fname, 'r') as h5f:
        return pd.DataFrame.from_records(h5f.root.Sensors.DataPMT[:])


def DataSiPM():
    fname = os.environ['ICDIR'] + '/Database/localdb.h5'
    with tb.open_file(fname, 'r') as h5f:
        return pd.DataFrame.from_records(h5f.root.Sensors.DataSiPM[:])


def DetectorGeo():
    fname = os.environ['ICDIR'] + '/Database/localdb.h5'
    with tb.open_file(fname, 'r') as h5f:
        return pd.DataFrame.from_records(h5f.root.Detector.DetectorGeometry[:])
