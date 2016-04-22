__author__ = 'alipirani'
import subprocess
import re
import os
import errno
import glob
import datetime
from config_settings import ConfigSectionMap
from check_subroutines import *


####################################################################### Post-Assembly processing using Bioawk #######################################################################################################
def bioawk(contigs, out_path, first_part):
    print "\n################## Removing Contigs less than 500 bp. ##################\n"
    contig_l500_cmd = "%s/%s/%s -c fastx '{ if(length($seq) > 500) { print \">\"$name; print $seq }}' %s > %s/%s_l500_contigs.fasta" % (ConfigSectionMap("bin_path")['binbase'], ConfigSectionMap("bioawk")['bioawk_bin'], ConfigSectionMap("bioawk")['base_cmd'], contigs, out_path, first_part)

    os.system(contig_l500_cmd)
    final_l500_contig = "%s/%s_l500_contigs.fasta" % (out_path, first_part)
    return final_l500_contig