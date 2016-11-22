"""
Sensor Functions
JJGC, September 2016
"""


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
