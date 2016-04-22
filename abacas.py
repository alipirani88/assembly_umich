__author__ = 'alipirani'
import subprocess
import re
import os
import errno
import glob
import datetime
from config_settings import ConfigSectionMap
from check_subroutines import *


####################################################################### Post-Assembly processing using Abacas #######################################################################################################
def abacas(reference_genome_path, final_l500_contig, out_path, first_part):
    print "\n################## Ordering Contigs using Abacas. ##################\n"
    abacas_cmd = "perl %s/%s/%s -r %s -q %s %s -o %s/%s_contigs_ordered" % (ConfigSectionMap("bin_path")['binbase'], ConfigSectionMap("abacas")['abacas_bin'], ConfigSectionMap("abacas")['base_cmd'], reference_genome_path, final_l500_contig, ConfigSectionMap("abacas")['abacas_parameters'], out_path, first_part)
    os.system(abacas_cmd)
    final_ordered_contigs = "%s/%s_contigs_ordered.fasta" % (out_path, first_part)
    return final_ordered_contigs