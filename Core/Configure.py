"""
Configure running options for the cities
JJGC August 2016
"""
from LogConfig import *
import pandas as pd
import getopt


def cdf_to_dict(cdf):
    """
    transforms the configuration data frame into a dictionary
    """

    dc ={}
    for k in cdf.keys():
        dc[k] = cdf[k][0]
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
         
         example of configuration file 

         # comment line  
        Names of parameters (comma separated)
        Values of parameters (comma separated)
        
        The parameters for DIOMIRA are:

        PATH_IN = path to input DST file (must be a MCRD file)
        FILE_IN = name of input DST file
        PATH_OUT = path to output DST file (RWF file)
        FILE_OUT = name of ouput DST file (RWF file)
        FIRST_EVT,LAST_EVT,RUN_ALL,

        RUN_ALL is used to decide whether to run all the events in the file
        in case that the total number of events requested (LAST_EVT-FIRST_EVT) 
        exceeds the number of events in the DST file. If RUN_ALL is set to 1 (True), 
        the script will run over all elements in the DST, 
        otherwise it will exit with a warning.


        """.format(program_name))

def configure(pname,argv):


    """
    reads arguments from the command line and configures job
    """
    
    #print("argv ={}".format(argv))
    
    DEBUG='INFO'
    INFO = False
    cfile =''
    CYTHON = False

    try:
        opts, args = getopt.getopt(argv, "hixd:c:", ["help","info","cython","debug","cfile"])

    except getopt.GetoptError:
        usage(pname)
        sys.exit(2)

    #print("opts ={}".format(opts))
    #print("args ={}".format(args))

    for opt, arg in opts:
        #print("opt ={}, arg = {}".format(opt,arg))
        if opt in ("-h", "--help"):
            usage(pname)
            sys.exit()
        elif opt in ("-d", "--debug"):
            DEBUG = arg
        elif opt in ("-i", "--info"):
            INFO = True
        elif opt in ("-c", "--cfile"):
            cfile = arg
        elif opt in ("-x", "--cython"):
            CYTHON = True
 
    lg = 'logging.'+DEBUG
    
    logger.setLevel(eval(lg))

    if cfile == '':
        print("Path to configuration file not given")
        usage(pname)
        sys.exit()

    cfp =pd.read_csv(cfile,comment="#")
    
    CFP = cdf_to_dict(cfp)
    
    logger.info("Configuration Parameters (CFP) dictionary  = {}".format(CFP))
    return DEBUG, INFO, CYTHON, CFP
    


if __name__ == '__main__':
    INFO, CFP = configure(sys.argv[0],sys.argv[1:])
