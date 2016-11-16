"""
Configure running options for the cities
JJGC August 2016
"""
import getopt
import sys
import os

from Core.LogConfig import logger


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

    CFP = read_config_file(cfile)
    return DEBUG, INFO, CFP


def define_event_loop(FIRST_EVT, LAST_EVT, NEVENTS, NEVENTS_DST, RUN_ALL):
    """
    defines the number of events to run in the loop
    """
    if RUN_ALL:
        return 0, NEVENTS_DST, NEVENTS_DST//20
    first = FIRST_EVT
    last = LAST_EVT
    if NEVENTS > NEVENTS_DST:
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

    print_mod = (last-first)//20 if last-first >= 20 else 1
    return first, last, print_mod


def cast(value):
    """
    Cast value from string to a python type.
    """
    if value == "True":
        return True
    if value == "False":
        return False
    if value.isdigit():
        return int(value)
    if value.replace(".", "").isdigit():
        return float(value)
    if "$" in value:
        value = os.path.expandvars(value)
    return value


def read_config_file(cfile):
    """
    Read a configuration file of the form
    PARAMETER VALUE
    """
    d = {}
    for line in open(cfile, "r"):
        if line == "\n" or line[0] == "#":
            continue
        tokens = line.rstrip().split(" ")
        key = tokens[0]

        value = map(cast, tokens[1:])
        d[key] = value[0] if len(value) == 1 else value
    return d


def filter_options(options, name):
    """
    Construct a new option dictionary with the parameters relevant to some
    module.

    Parameters
    ----------
    options : dictionary
        Dictionary of options with format "MODULE:PARAMETER": value.
    name : string
        Selected module name.

    Returns
    -------
    out : dictionary
        Filtered dictionary.
    """
    out = {}
    for key, value in options.items():
        if name in key:
            out[key.split(":")[1]] = value
    return out
