__author__ = 'alipirani'
import subprocess
import re
import os
import errno
import glob
import datetime
from config_settings import ConfigSectionMap
from modules.logging_subprocess import *
from modules.log_modules import keep_logging


def ariba_AMR(forward_paired, reverse_paired, output_folder, prefix, logger, Config):
    ariba_AMR_dir = output_folder + "/" + prefix + "_AMR"
    ariba_cmd = "%s run --force --verbose %s %s %s %s" % (ConfigSectionMap("ariba", Config)['base_cmd'], ConfigSectionMap("ariba", Config)['ariba_amr_db'], forward_paired, reverse_paired, ariba_AMR_dir)
    #print ariba_cmd
    keep_logging("Using Ariba DB path mentioned in config file: %s" % ConfigSectionMap("ariba", Config)['ariba_amr_db'], "Using Ariba DB path mentioned in config file: %s" % ConfigSectionMap("ariba", Config)['ariba_amr_db'], logger, 'info')
    try:
        keep_logging(ariba_cmd, ariba_cmd, logger, 'debug')
        call(ariba_cmd, logger)
    except sp.CalledProcessError:
        keep_logging('Error in running Ariba AMR. Exiting.', 'Error in running Ariba AMR. Exiting.', logger, 'exception')
        sys.exit(1)


def ariba_MLST(forward_paired, reverse_paired, output_folder, prefix, logger, Config):
    ariba_mlst_dir = output_folder + "/" + prefix + "_MLST"
    ariba_mlst_cmd = "%s run --force --verbose %s %s %s %s" % (ConfigSectionMap("ariba", Config)['base_cmd'], ConfigSectionMap("ariba", Config)['ariba_mlst_db'], forward_paired, reverse_paired, ariba_mlst_dir)
    #print ariba_mlst_cmd
    keep_logging("Using Ariba DB path mentioned in config file: %s" % ConfigSectionMap("ariba", Config)['ariba_mlst_db'], "Using Ariba DB path mentioned in config file: %s" % ConfigSectionMap("ariba", Config)['ariba_mlst_db'], logger, 'info')
    try:
        keep_logging(ariba_mlst_cmd, ariba_mlst_cmd, logger, 'debug')
        call(ariba_mlst_cmd, logger)
    except sp.CalledProcessError:
        keep_logging('Error in running Ariba MLST. Exiting.', 'Error in running Ariba MLST. Exiting.', logger, 'exception')
        sys.exit(1)
