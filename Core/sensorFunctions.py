"""
Sensor Functions
JJGC, September 2016
"""
import pandas as pd
import numpy as np


def read_data_geom(geom_t):
    """
    Reads the geom table en returns a PD Series
    """
    ga = geom_t.read()
    G = pd.Series([ga[0][0][0], ga[0][0][1], ga[0][1][0], ga[0][1][1],
                   ga[0][2][0], ga[0][2][1], ga[0][3]],
                  index=["xdet_min", "xdet_max", "ydet_min", "ydet_max",
                         "zdet_min", "zdet_max", "R"])
    return G


def read_data_FEE(fee_t):
    """
    Reads the FEE table en returns a PD Series for the simulation parameters
    and a PD series for the values of the capacitors used in the simulation
    """

    fa = fee_t.read()
    F = pd.Series([fa[0][0], fa[0][1], fa[0][2], fa[0][3], fa[0][5], fa[0][6],
                   fa[0][7], fa[0][8], fa[0][9], fa[0][10], fa[0][11],
                   fa[0][12]],
                  index=["offset", "pmt_gain", "V_gain", "R", "time_step",
                         "time_daq", "freq_LPF", "freq_HPF", "LSB",
                         "volts_to_adc", "noise_fee_rms", "noise_adc"])
    C = pd.Series([fa[0][4]], index=["C12"])
    return F, C


def get_column_(pmta, ic):
    """
    Access column ic of table pmta and returns column as an array
    """
    col = []
    for i in range(pmta.shape[0]):
        col.append(pmta[i][ic])
    return np.array(col)


def read_data_sensors(sensor_table):
    """
    Reads the sensors table and returns a data frame
    """
    pmta = sensor_table.read()
    PMT = {}
    PMT["channel"] = get_column_(pmta, 0)
    PMT["sensorID"] = get_column_(pmta, 1)
    PMT["x"], PMT["y"], dummy = get_column_(pmta, 2).T
    PMT["coeff"] = get_column_(pmta, 3)
    PMT["adc_to_pes"] = get_column_(pmta, 4)
    PMT["noise_rms"] = get_column_(pmta, 5)
    return pd.DataFrame(PMT)


def get_energy_sensors(energy_v):
    """
    Reads the sensors energy and returns a data frame
    """
    return pd.DataFrame(energy_v.read())


def find_sensor(x0, y0, df):
    """
    Finds the sensor with given coordinates.

    Parameters
    ----------
    x0, y0 : floats
        x and y coordinates of the sensors.
    df : pandas DataFrame
        Contains the sensors' information.

    Returns
    -------
    index : int
        Index in the data frame
    ID : int
        Channel ID of the sensor with the given coordinates.
    """
    data = df[ (df.x == x0) & (df.y == y0) ]['channel']
    return data.index[0], data.values[0]
