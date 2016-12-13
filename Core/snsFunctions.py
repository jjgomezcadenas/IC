"""
Sensor Functions
JJGC, September 2016
"""

import numpy as np


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
    data = df[(df.x == x0) & (df.y == y0)]['channel']
    return data.index[0], data.values[0]


def baricenter(q, sipmdf, thrs=0.25):
    """
    Compute the baricenter of a slice.

    Parameters
    ----------
    q : 1-dim np.ndarray
        Signal for each SiPM.
    sipmdf : pd.DataFrame
        Contains the sensors' information.
    thrs : float, optional
        Cut value to filter data based on the signal relative to maximum.
        Defaults to 0.25.

    Returns
    -------
    x, y: floats
        Weighted average for the x and y coordinates.
    """
    selection = q > np.max(q) * thrs
    q = q[selection]
    if q.sum() < 1.:
        return 1e4, 1e4
    x = np.average(sipmdf["X"].values[selection], weights=q)
    y = np.average(sipmdf["Y"].values[selection], weights=q)
    return x, y
