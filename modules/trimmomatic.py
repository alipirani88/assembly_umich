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
#from check_subroutines import *


####################################################################### Raw data Pre-processing using Trimmomatic #######################################################################################################
def clean_reads(input1, input2, out_path, crop, logger, Config):
    if input2 != "None":
        forward_paired = out_path + ConfigSectionMap("Trimmomatic", Config)['f_p']
        reverse_paired = out_path + ConfigSectionMap("Trimmomatic", Config)['r_p']
        forward_unpaired = out_path + ConfigSectionMap("Trimmomatic", Config)['f_up']
        reverse_unpaired = out_path + ConfigSectionMap("Trimmomatic", Config)['r_up']
        adapter_file = ConfigSectionMap("bin_path", Config)['binbase'] + "/" + ConfigSectionMap("Trimmomatic", Config)['trimmomatic_bin'] + "/" + ConfigSectionMap("Trimmomatic", Config)['adaptor_filepath']
        clean_filenames = forward_paired + " " + forward_unpaired + " " + reverse_paired + " " + reverse_unpaired
        illumina_string = 'ILLUMINACLIP:' + adapter_file + ConfigSectionMap("Trimmomatic", Config)['colon'] + ConfigSectionMap("Trimmomatic", Config)['seed_mismatches'] + ConfigSectionMap("Trimmomatic", Config)['colon'] + ConfigSectionMap("Trimmomatic", Config)['palindrome_clipthreshold'] + ConfigSectionMap("Trimmomatic", Config)['colon'] + ConfigSectionMap("Trimmomatic", Config)['simple_clipthreshold']
        sliding_string = 'SLIDINGWINDOW:' + ConfigSectionMap("Trimmomatic", Config)['window_size'] + ConfigSectionMap("Trimmomatic", Config)['colon'] + ConfigSectionMap("Trimmomatic", Config)['window_size_quality']
        minlen_string = 'MINLEN:' + ConfigSectionMap("Trimmomatic", Config)['minlength']
        headcrop_string = 'HEADCROP:' + ConfigSectionMap("Trimmomatic", Config)['headcrop_length']
        if not crop:
            cmdstring = "java -jar " + ConfigSectionMap("bin_path", Config)['binbase'] + ConfigSectionMap("Trimmomatic", Config)['trimmomatic_bin'] + "trimmomatic-0.36.jar PE -phred33 " + input1 + " " + input2 + " " + clean_filenames + " " + illumina_string + " " + sliding_string + " " + minlen_string + " " + headcrop_string
            keep_logging(cmdstring, cmdstring, logger, 'debug')
            try:
                call(cmdstring, logger)
                print ""
            except sp.CalledProcessError:
                keep_logging('Error in Trimmomatic Pre-processing step. Exiting.', 'Error in Trimmomatic Pre-processing step. Exiting.', logger, 'exception')
                sys.exit(1)
            return forward_paired, reverse_paired, forward_unpaired, reverse_unpaired
        else:
            crop_string = 'CROP:' + crop
            cmdstring = "java -jar " + ConfigSectionMap("bin_path", Config)['binbase'] + ConfigSectionMap("Trimmomatic", Config)['trimmomatic_bin'] + "trimmomatic-0.36.jar PE " + input1 + " " + input2 + " " + clean_filenames + " " + crop_string + " " + illumina_string + " " + sliding_string + " " + minlen_string
            keep_logging(cmdstring, cmdstring, logger, 'debug')
            try:
                call(cmdstring, logger)
                print ""
            except sp.CalledProcessError:
                keep_logging('Error in Trimmomatic Pre-processing step. Exiting.', 'Error in Trimmomatic Pre-processing step. Exiting.', logger, 'exception')
                sys.exit(1)
            return forward_paired, reverse_paired, forward_unpaired, reverse_unpaired
    else:
        forward_paired = out_path + ConfigSectionMap("Trimmomatic", Config)['f_p']
        # reverse_paired = out_path + ConfigSectionMap("Trimmomatic", Config)['r_p']
        #forward_unpaired = out_path + ConfigSectionMap("Trimmomatic", Config)['f_up']
        # reverse_unpaired = out_path + ConfigSectionMap("Trimmomatic", Config)['r_up']
        adapter_file = ConfigSectionMap("bin_path", Config)['binbase'] + "/" + ConfigSectionMap("Trimmomatic", Config)['trimmomatic_bin'] + "/" + ConfigSectionMap("Trimmomatic", Config)['adaptor_filepath_se']
        clean_filenames = forward_paired
        illumina_string = 'ILLUMINACLIP:' + adapter_file + ConfigSectionMap("Trimmomatic", Config)['colon'] + ConfigSectionMap("Trimmomatic", Config)['seed_mismatches'] + ConfigSectionMap("Trimmomatic", Config)['colon'] + ConfigSectionMap("Trimmomatic", Config)['palindrome_clipthreshold'] + ConfigSectionMap("Trimmomatic", Config)['colon'] + ConfigSectionMap("Trimmomatic", Config)['simple_clipthreshold']
        sliding_string = 'SLIDINGWINDOW:' + ConfigSectionMap("Trimmomatic", Config)['window_size'] + ConfigSectionMap("Trimmomatic", Config)['colon'] + ConfigSectionMap("Trimmomatic", Config)['window_size_quality']
        minlen_string = 'MINLEN:' + ConfigSectionMap("Trimmomatic", Config)['minlength']
        headcrop_string = 'HEADCROP:' + ConfigSectionMap("Trimmomatic", Config)['headcrop_length']
        if not crop:
            cmdstring = "java -jar " + ConfigSectionMap("bin_path", Config)['binbase'] + ConfigSectionMap("Trimmomatic", Config)['trimmomatic_bin'] + "trimmomatic-0.36.jar SE -phred33 " + input1 + " " + clean_filenames + " " + illumina_string + " " + sliding_string + " " + minlen_string + " " + headcrop_string
            keep_logging(cmdstring, cmdstring, logger, 'debug')
            try:
                call(cmdstring, logger)
                print ""
            except sp.CalledProcessError:
                keep_logging('Error in Trimmomatic Pre-processing step. Exiting.', 'Error in Trimmomatic Pre-processing step. Exiting.', logger, 'exception')
                sys.exit(1)
            reverse_paired = "None"
            reverse_unpaired = "None"
            forward_unpaired = "None"
            return forward_paired, reverse_paired, forward_unpaired, reverse_unpaired
        else:
            crop_string = 'CROP:' + crop
            cmdstring = "java -jar " + ConfigSectionMap("bin_path", Config)['binbase'] + ConfigSectionMap("Trimmomatic", Config)['trimmomatic_bin'] + "trimmomatic-0.36.jar SE " + input1 + " " + clean_filenames + " " + crop_string + " " + illumina_string + " " + sliding_string + " " + minlen_string
            keep_logging(cmdstring, cmdstring, logger, 'debug')
            try:
                call(cmdstring, logger)
                print ""
            except sp.CalledProcessError:
                keep_logging('Error in Trimmomatic Pre-processing step. Exiting.', 'Error in Trimmomatic Pre-processing step. Exiting.', logger, 'exception')
                sys.exit(1)
            reverse_paired = "None"
            reverse_unpaired = "None"
            forward_unpaired = "None"
            return forward_paired, reverse_paired, forward_unpaired, reverse_unpaired
################################################################################### End #################################################################################################################################
