_author__ = 'alipirani'
import subprocess
import re
import os
import errno
import glob
import datetime
from config_settings import ConfigSectionMap
from modules.log_modules import *
from modules.logging_subprocess import *


####################################################################### Genome annotation using Prokka #######################################################################################################
def prokka(final_ordered_contigs, out_path, first_part, logger, Config):
    keep_logging('Genome Annotation using Prokka', 'Genome Annotation using Prokka', logger, 'info')
    prokka_cmd = "%s -outdir %s/%s_prokka/ -prefix %s %s %s" % (ConfigSectionMap("prokka", Config)['base_cmd'], out_path, first_part, first_part, ConfigSectionMap("prokka", Config)['prokka_parameters'], final_ordered_contigs)
    try:
        keep_logging(prokka_cmd, prokka_cmd, logger, 'debug')
        call(prokka_cmd, logger)
        #print ""
    except sp.CalledProcessError:
        keep_logging('Error in Prokka Annotation step. Exiting.', 'Error in Prokka Annotation Step. Exiting.', logger, 'exception')
        sys.exit(1)
    final_annotation_folder = "%s/%s_prokka/" % (out_path, first_part)
    return final_annotation_folder
