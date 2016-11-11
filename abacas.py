__author__ = 'alipirani'
import subprocess
import re
import os
import errno
import glob
import datetime
from config_settings import ConfigSectionMap
from logging_subprocess import *
from log_modules import keep_logging


####################################################################### Post-Assembly processing using Abacas #######################################################################################################
def abacas(reference_genome_path, final_l500_contig, out_path, first_part, logger, Config):
    keep_logging('Contig Reordering using ABACAS', 'Contig Reordering using ABACAS', logger, 'info')
    abacas_cmd = "perl %s/%s/%s -r %s -q %s %s -o %s/%s_contigs_ordered" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("abacas", Config)['abacas_bin'], ConfigSectionMap("abacas", Config)['base_cmd'], reference_genome_path, final_l500_contig, ConfigSectionMap("abacas", Config)['abacas_parameters'], out_path, first_part)
    try:
        keep_logging(abacas_cmd, abacas_cmd, logger, 'debug')
        call(abacas_cmd, logger)
        #print ""
    except sp.CalledProcessError:
        keep_logging('Error in reordering Contigs using Abacas. Exiting.', 'Error in reordering Contigs using Abacas. Exiting.', logger, 'exception')
        sys.exit(1)
    final_ordered_contigs = "%s/%s_contigs_ordered.fasta" % (out_path, first_part)
    return final_ordered_contigs