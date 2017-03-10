__author__ = 'alipirani'
import subprocess
import re
import os
import errno
import glob
import datetime
import sys
from config_settings import ConfigSectionMap

# Quality assessment using FatQC. Runs only when --qa option is specified.
def qa_fastqc(out_path, file_1, file_2):
    print "\n################## Running FASTQC on input files. ##################\n"
    fastqc_dir = os.path.join(out_path, "fastqc" + datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S') + "/")
    cmdstring = "mkdir " + fastqc_dir
    os.system(cmdstring)
    fastqc = ConfigSectionMap("bin_path")['binbase'] + "FastQC/fastqc"
    # check both the files
    if file_1 and file_2:
        cmd = "%s %s %s -o %s" % (fastqc, file_1, file_2, fastqc_dir)
        print "Running: %s" % cmd
        os.system(cmd)
        fastqc_zipfolder = fastqc_dir + "*.zip"
        fastqc_zipfoldernames = glob.glob(fastqc_zipfolder)
        print "\n" + "FastQC results are saved in resp. folders.\n"
        print fastqc_dir
        for x in fastqc_zipfoldernames:
            p_unzip = subprocess.Popen(["unzip", x, "-d", fastqc_dir], stdout=subprocess.PIPE)
            p_unzip.communicate()
        fastqc_results = fastqc_dir + "*_fastqc"
        fastqc_resultsfolder = glob.glob(fastqc_results)
        for y in fastqc_resultsfolder:
            print "\n" + "FastQC summary of Sample Files:\n"
            summary_file = y + "/summary.txt"
            base = os.path.basename(y)
            f = open(summary_file,"r")
            for line in f:
                print line
    else:
        cmd = "%s %s %s -o %s" % (fastqc, file_1, file_2, fastqc_dir)
        print "Running: %s" % cmd
        os.system(cmd)
        fastqc_zipfolder = fastqc_dir + "*.zip"
        print "FastQC results are saved in resp. folders.\n"
        print fastqc_zipfolder
        fastqc_zipfoldernames = glob.glob(fastqc_zipfolder)
        print fastqc_zipfoldernames
        for x in fastqc_zipfoldernames:
            p_unzip = subprocess.Popen(["unzip", x, "-d", fastqc_dir], stdout=subprocess.PIPE)
            p_unzip.communicate()
        fastqc_results = fastqc_dir + "*_fastqc"
        fastqc_resultsfolder = glob.glob(fastqc_results)
        for y in fastqc_resultsfolder:
            print "FastQC summary of Sample Files:\n"
            summary_file = y + "/summary.txt"
            f = open(summary_file,"r")
            for line in f:
                print line
    print "\n################## END: FASTQC ##################\n"
