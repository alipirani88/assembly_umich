__author__ = 'alipirani'

import subprocess
import re
import os
import errno
import glob
from config_settings import *




# # Usage Message:
# def usage():
#     print "Usage: pipeline.py [-h] [-v] [-f1 FILE_1] [-f2 FILE_2] [-o OUTPUT_FOLDER] [-start_step START] [-end_step END] [--qa] [-A Assembler] [--err] \n"
#
# # Check Java Availability:
# def java_check():
#     print "Checking Java Availability....\n"
#     jd = subprocess.check_output(["java", "-version"], stderr=subprocess.STDOUT)
#     if len(jd) < 1:
#         print "Unable to find a java runtime environment. The pipeline requires java 6 or later."
#     else:
#         print "Java Availability Check completed ...\n\n" + jd
#
# # Make sure input raw reads files exists at given location.
# def file_exists(path1, path2):
#     if path1:
#         if not os.path.isfile(path1):
#             file_basename = os.path.basename(path1)
#             print "The input file " + file_basename + " does not exists. \nPlease provide another file.\n"
#             exit()
#     if path2:
#         if not os.path.isfile(path2):
#             file_basename = os.path.basename(path2)
#             print "The input file " + file_basename + " does not exists. \nPlease provide another file.\n"
#             exit()
#
# # Make sure the output folder exists or create at given path
# def make_sure_path_exists(out_path):
#     try:
#         os.makedirs(out_path)
#     except OSError as exception:
#         if exception.errno != errno.EEXIST:
#             print "Errors in output folder path! please change the output path or analysis name\n"
#             exit()
#
# # Prepare the paths for clean reads when all the input files exists.
# def prepare_cleanreadsinput():
#     global forward_paired, forward_unpaired, reverse_unpaired, reverse_paired
#     forward_paired = out_name + "forward_paired.fq.gz"
#     forward_unpaired = out_name + "forward_unpaired.fq.gz"
#     reverse_paired = out_name + "reverse_paired.fq.gz"
#     reverse_unpaired = out_name + "reverse_unpaired.fq.gz"
#     return forward_paired, forward_unpaired, reverse_paired, reverse_unpaired
#
# # check if the paired and unpaired clean reads exists in output directory. Set the paired and unpaired constants accordingly.
# def check_cleanreads(forward_paired, reverse_paired, forward_unpaired, reverse_unpaired):
#     # Initialize the constants
#     unpaired = "1"
#     paired = "1"
#     if not os.path.isfile(forward_paired):
#         print "The pre-processed paired reads file does not exists. This step requires clean reads namely: forward_paired.fq.gz, reverse_paired.fq.gz.\n"
#         paired = "0"
#     if not os.path.isfile(forward_unpaired):
#         print "The pre-processed unpaired reads file does not exists. This step requires (optional) unpaired clean reads named: forward_unpaired.fq.gz, reverse_unpaired.fq.gz.\n"
#         unpaired = "0"
#     return paired, unpaired
#
# # Set appropriate path for contigs/scaffolds fasta files. Needed only when quast or reapr runs individually.
# def get_contigs(out_path, assembler):
#     if assembler == "velvet":
#         contigs = out_path + ConfigSectionMap("velvet")['contigs_path']
#         scaffolds = ""
#         return contigs, scaffolds
#     elif assembler == "spades":
#         contigs = out_path + ConfigSectionMap("spades")['contigs_path']
#         scaffolds = out_path + ConfigSectionMap("spades")['scaffolds_path']
#         plasmid_contigs = out_path + ConfigSectionMap("spades")['plasmid_contigs_path']
#         plasmid_scaffolds = out_path + ConfigSectionMap("spades")['plasmid_scaffolds_path']
#         return contigs, scaffolds, plasmid_contigs, plasmid_scaffolds
# ###########################################################################################################################################################################
#
# # Check if the clean Paired and Unpaired reads exists in the output folder.
# # def check_cleanreads(out_path):
# #     paired_path = out_path + "*_paired.fq.gz"
# #     paired_files = glob.glob(paired_path)
# #     paired = ConfigSectionMap("check_clean")['paired']
# #     unpaired = ConfigSectionMap("check_clean")['unpaired']
# #     if not paired_files:
# #         #file_basename = os.path.basename(path)
# #         print "The pre-processed paired reads file does not exists. This step requires clean reads namely: forward_paired.fq.gz, reverse_paired.fq.gz.\n"
# #         paired = "0"
# #     unpaired_path = out_path + "*_unpaired.fq.gz"
# #     unpaired_files = glob.glob(unpaired_path)
# #     if not unpaired_files:
# #         #file_basename = os.path.basename(path)
# #         print "The pre-processed unpaired reads file does not exists. This step requires (optional) unpaired clean reads named: forward_unpaired.fq.gz, reverse_unpaired.fq.gz.\n"
# #         unpaired = "0"
# #     set_config_values(paired, unpaired)
# #     #return paired, unpaired
#
# ###########################################################################################################################################################################










