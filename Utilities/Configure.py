def usage(program_name):
    """
    Usage of program
    """
    print("""
        Usage: python (run) % [args]
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
def configure(argv):
    """
    reads arguments from the command line and configures job
    """
    
    print("argv ={}".format(argv))
    global DEBUG, PATH_IN, PATH_OUT, FILE_IN, FILE_OUT
    global  FIRST_EVT, LAST_EVT,NEVENTS, RUN_ALL, INFO
    
    DEBUG='INFO'
    INFO = False
    cfile =''
    try:
        opts, args = getopt.getopt(argv, "hid:c:", ["help","info","debug","cfile"])

    except getopt.GetoptError:
        usage()
        sys.exit(2)

    print("opts ={}".format(opts))

    for opt, arg in opts:
        print("opt ={}, arg = {}".format(opt,arg))
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-d", "--debug"):
            DEBUG = arg
        elif opt in ("-i", "--info"):
            INFO = True
        elif opt in ("-c", "--cfile"):
            cfile = arg
 
    lg = 'logging.'+DEBUG
    print('INFO = {} DEBUG={}'.format(INFO,DEBUG))
    print("lg ={}".format(lg))
    logging.basicConfig(level=eval(lg))

    if cfile == '':
        print("Path to configuration file not given")
        usage()
        sys.exit()

    CFP =pd.read_csv(cfile,comment="#")
    print("""
        Configuration parameters \n 
        {}
        """.format(CFP))

    PATH_IN=CFP['PATH_IN'][0] 
    PATH_OUT=CFP['PATH_OUT'][0]
    FILE_IN=CFP['FILE_IN'][0]
    FILE_OUT=CFP['FILE_OUT'][0]
    FIRST_EVT=CFP['FIRST_EVT'][0]
    LAST_EVT=CFP['LAST_EVT'][0]
    RUN_ALL=CFP['RUN_ALL'][0]
    NEVENTS = LAST_EVT -  FIRST_EVT
