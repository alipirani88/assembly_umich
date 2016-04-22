__author__ = 'alipirani'
import subprocess
import re
import os
import errno
import glob
import datetime
from config_settings import ConfigSectionMap
from check_subroutines import *


####################################################################### Raw data Pre-processing using Trimmomatic #######################################################################################################
def clean_reads(input1, input2, out_path, crop):
    print "\n################## Pre-processing using Trimmomatic. ##################\n"
    make_sure_path_exists(out_path)
    forward_paired = out_path + ConfigSectionMap("Trimmomatic")['f_p']
    reverse_paired = out_path + ConfigSectionMap("Trimmomatic")['r_p']
    forward_unpaired = out_path + ConfigSectionMap("Trimmomatic")['f_up']
    reverse_unpaired = out_path + ConfigSectionMap("Trimmomatic")['r_up']
    adapter_file = ConfigSectionMap("bin_path")['binbase'] + "/" + ConfigSectionMap("Trimmomatic")['trimmomatic_bin'] + "/" + ConfigSectionMap("Trimmomatic")['adaptor_filepath']
    clean_filenames = forward_paired + " " + forward_unpaired + " " + reverse_paired + " " + reverse_unpaired
    illumina_string = 'ILLUMINACLIP:' + adapter_file + ConfigSectionMap("Trimmomatic")['colon'] + ConfigSectionMap("Trimmomatic")['seed_mismatches'] + ConfigSectionMap("Trimmomatic")['colon'] + ConfigSectionMap("Trimmomatic")['palindrome_clipthreshold'] + ConfigSectionMap("Trimmomatic")['colon'] + ConfigSectionMap("Trimmomatic")['simple_clipthreshold']
    sliding_string = 'SLIDINGWINDOW:' + ConfigSectionMap("Trimmomatic")['window_size'] + ConfigSectionMap("Trimmomatic")['colon'] + ConfigSectionMap("Trimmomatic")['window_size_quality']
    minlen_string = 'MINLEN:' + ConfigSectionMap("Trimmomatic")['minlength']
    headcrop_string = 'HEADCROP:' + ConfigSectionMap("Trimmomatic")['headcrop_length']
    if not crop:
        cmdstring = "java -jar " + ConfigSectionMap("bin_path")['binbase'] + ConfigSectionMap("Trimmomatic")['trimmomatic_bin'] + "trimmomatic-0.33.jar PE " + input1 + " " + input2 + " " + clean_filenames + " " + illumina_string + " " + sliding_string + " " + minlen_string + " " + headcrop_string
        print "\nRunning:\n [%s] \n" % cmdstring
        os.system(cmdstring)
        print "\n################## End: Data Pre-processing ##################\n"
        return forward_paired, reverse_paired, forward_unpaired, reverse_unpaired
    else:
        crop_string = 'CROP:' + crop
        cmdstring = "java -jar " + ConfigSectionMap("bin_path")['binbase'] + ConfigSectionMap("Trimmomatic")['trimmomatic_bin'] + "trimmomatic-0.33.jar PE " + input1 + " " + input2 + " " + clean_filenames + " " + crop_string + " " + illumina_string + " " + sliding_string + " " + minlen_string
        print "\nRunning:\n [%s] \n" % cmdstring
        os.system(cmdstring)
        print "\n################## End: Data Pre-processing ##################\n"
        return forward_paired, reverse_paired, forward_unpaired, reverse_unpaired
################################################################################### End #################################################################################################################################
