_author__ = 'alipirani'
import subprocess
import re
import os
import errno
import glob
import datetime
from config_settings import ConfigSectionMap
from check_subroutines import *


####################################################################### Genome annotation using Prokka #######################################################################################################
def prokka(final_ordered_contigs, out_path, first_part):
    print "\n################## Genome Annotation using Prokka. ##################\n"
    prokka_cmd = "%s/%s/%s -outdir %s/%s_prokka/ -prefix %s %s %s" % (ConfigSectionMap("bin_path")['binbase'], ConfigSectionMap("prokka")['prokka_bin'], ConfigSectionMap("prokka")['base_cmd'], out_path, first_part, first_part, ConfigSectionMap("prokka")['prokka_parameters'], final_ordered_contigs)
    os.system(prokka_cmd)
    final_annotation_folder = "%s/%s_prokka/" % (out_path, first_part)
    return final_annotation_folder
