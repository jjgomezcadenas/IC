"""
Configure running options for the cities
JJGC August 2016
"""
from LogConfig import logger
import pandas as pd
import getopt
import sys
import os


def cdf_to_dict(cdf):
    """
    transforms the configuration data frame into a dictionary
    """
    dc = {}
    for k in cdf.keys():
        value = cdf[k][0]
        if isinstance(value, str) and "$" in value:
            value = os.path.expandvars(value)
        dc[k] = value
    dc["PATH_DB"] = os.environ["ICDBDIR"]
    return dc


def usage(program_name):
    """
    Usage of program
    """
    print("""
        Usage: python (run) {} [args]


        where args are:
         -h (--help) : this text
         -i (--info) : print a text describing the invisible city of DIOMIRA
         -d (--debug) : can be set to 'DEBUG','INFO','WARNING','ERROR'
         -c (--cfile) : full path to a configuration file

        """.format(program_name))


def configure(pname, argv):
    """
    reads arguments from the command line and configures job
    """
    DEBUG = 'INFO'
    INFO = False
    cfile = ''

    try:
        opts, args = getopt.getopt(argv,
                                   "hid:c:",
                                   ["help", "info", "debug", "cfile"])
    except getopt.GetoptError:
        usage(pname)
        sys.exit(2)

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage(pname)
            sys.exit()
        elif opt in ("-d", "--debug"):
            DEBUG = arg
        elif opt in ("-i", "--info"):
            INFO = True
        elif opt in ("-c", "--cfile"):
            cfile = arg

    logger.setLevel(DEBUG)

    if cfile == '':
        print("Path to configuration file not given")
        usage(pname)
        sys.exit()

    cfp = pd.read_csv(cfile, comment="#")
    CFP = cdf_to_dict(cfp)
    return DEBUG, INFO, CFP


def define_event_loop(FIRST_EVT, LAST_EVT, NEVENTS, NEVENTS_DST, RUN_ALL):
    """
    defines the number of events to run in the loop
    """
    if RUN_ALL:
        return 0, NEVENTS_DST, NEVENTS_DST//20
    first = FIRST_EVT
    last = LAST_EVT
    if NEVENTS > NEVENTS_DST and RUN_ALL is False:
        print("""
                Refusing to run: you have requested
                FIRST_EVT = {}
                LAST_EVT  = {}
                Thus you want to run over {} events
                but the size of the DST is {} events.
                Please change your choice or select RUN_ALL = TRUE
                to run over the whole DST when this happens
                """.format(FIRST_EVT, LAST_EVT, NEVENTS, NEVENTS_DST))
        sys.exit(0)

    elif NEVENTS > NEVENTS_DST and RUN_ALL is True:
        first = 0
        last = NEVENTS_DST
    return first, last, (last-first)//20

#
#
#
#
# if __name__ == '__main__':
#     INFO, CFP = configure(sys.argv[0],sys.argv[1:])
