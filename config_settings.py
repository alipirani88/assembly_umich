__author__ = 'alipirani'

import ConfigParser
Config = ConfigParser.ConfigParser()
parser = ConfigParser.SafeConfigParser()
#change the name of config file in case using different one.
Config.read("/home2/apirani/bin/assembly_umich/config")
parser.read("/home2/apirani/bin/assembly_umich/config")

###############Read Config File########################
def ConfigSectionMap(section):
    dict1 = {}
    options = Config.options(section)
    for option in options:
        try:
            dict1[option] = Config.get(section, option)
            if dict1[option] == -1:
                DebugPrint("skip: %s" % option)
        except:
            print("exception on %s!" % option)
            dict1[option] = None
    return dict1
###############Read Config File########################

def set_config_values(paired, unpaired):
    #parser.add_section('check_clean')
    parser.set('check_clean', 'paired', paired)
    parser.set('check_clean', 'unpaired', unpaired)
    # Writing our configuration file to 'example.ini'
    with open('./config', 'wb') as configfile:
        parser.write(configfile)

def create_config_section():
    #parser.add_section('check_clean')
    parser.set('check_clean', 'paired', '1')
    parser.set('check_clean', 'unpaired', '1')
    with open('./config', 'wb') as configfile:
        #config.write(configfile)
        parser.write(configfile)
