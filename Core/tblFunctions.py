"""
Table Functions
GML, October 2016

ChangeLog
14/10: import table io functions from wfmFunctions
"""

from __future__ import print_function
import pandas as pd
import wfmFunctions as wfm
from LogConfig import *

def store_wf(event, table, WF):
    """
    Store a wavform in a table
    """
    row = table.row
    for isens,wf in WF.items():
        for t,e in zip(wf.time_mus, wf.ene_pes):
            row['event'] = event
            row['ID'] = isens
            row['time_mus'] = t
            row['ene_pes'] = e
            row.append()
    table.flush()

# def read_twf(pmttwf, event_number):
#     """
#     Reads back the TWF: old version kept for backward compatibility to
#     be deleted asap
#     """
#     PMT={}
#     for row in table.where("event == event_number"):
#         pmt = row['ID']
#         time_mus = row['time_mus']
#         ene_pes =  row['ene_pes']
#
#         #print('pmt = {},time_mus = {},ene_pes = {}'.format(pmt,time_mus,ene_pes))
#         if pmt not in PMT:
#             WF={}
#             TIME =[]
#             ENE = []
#             TIME.append(time_mus)
#             ENE.append(ene_pes)
#             WF['time_mus'] = TIME
#             WF['ene_pes'] = ENE
#             PMT[pmt] = WF
#         else:
#             WF = PMT[pmt]
#             TIME = WF['time_mus']
#             ENE  = WF['ene_pes']
#             TIME.append(time_mus)
#             ENE.append(ene_pes)
#
#     return PMT

# def read_twf(twf, event_number):
#     """
#     Reads back the TWF of the PMTs/SiPMs for event number:
#     input: the twf table of the PMTs,(SiPMs) a list with the PMT (SiPMs) indexes and the event number
#     outputs: a PMT/SiPM panel
#
#     """
#     sensors ={}
#     for isensor in sensor_list:
#         try:
#             time_mus, ene_pes = zip(*[ (row['time_mus'],row['ene_pes']) for row in twf.iterrows() if row['event']== event_number and row['ID']== isensor])
#             sensors[isensor] = wf2df(time_mus,ene_pes)
#         except ValueError:
#             logger.error('found an empty sensor')
#             exit()
#
#     return pd.Panel(sensors)

def read_wf(table,event_number,isens):
    '''
        Reads table and returns the waveform (time_mus and ene_pes) corresponding
        to sensor isens of event event_number.
    '''
    try:
        return zip(*[ (row['time_mus'],row['ene_pes']) for row in table.iterrows() if row['event']== event_number and row['ID']== isens])
    except:
        logger.error('[read_wf]: empty sensor found: {}'.format(isens))


def read_wf_table( table, event_number ):
    """
    Reads back the TWF of the PMTs/SiPMs for event number:
    input: the twf table of the PMTs,(SiPMs) a list with the PMT (SiPMs) indexes and the event number
    outputs: a PMT/SiPM panel

    """
    sensor_list = set(table.read_where('event == {}'.format(event_number),field='ID'))
    return pd.Panel({ isens : wfm.wf2df(*read_wf(table,event_number,isens)) for isens in sensor_list})
