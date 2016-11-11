__author__ = 'alipirani'
import subprocess
import re
import os
import errno
import glob
import datetime
from config_settings import ConfigSectionMap
from log_modules import *
from logging_subprocess import *


####################################################################### Post-Assembly processing using Bioawk #######################################################################################################
def bioawk(contigs, plasmid_contigs, out_path, first_part, logger, Config):
    keep_logging('Removing Contigs less than 500 bp using BIOAWK', 'Removing Contigs less than 500 bp using BIOAWK', logger, 'info')
    contig_l500_cmd = "%s/%s/%s -c fastx '{ if(length($seq) > 500) { print \">\"$name; print $seq }}' %s > %s/%s_l500_contigs.fasta" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("bioawk", Config)['bioawk_bin'], ConfigSectionMap("bioawk", Config)['base_cmd'], contigs, out_path, first_part)
    plasmid_contig_l500_cmd = "%s/%s/%s -c fastx '{ if(length($seq) > 500) { print \">\"$name; print $seq }}' %s > %s/%s_l500_plasmid_contigs.fasta" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("bioawk", Config)['bioawk_bin'], ConfigSectionMap("bioawk", Config)['base_cmd'], plasmid_contigs, out_path, first_part)
    try:
        keep_logging(contig_l500_cmd, contig_l500_cmd, logger, 'debug')
        call(contig_l500_cmd, logger)
        #print ""
        keep_logging('Contigs greater than 500 bp are: {}/{}_l500_contigs.fasta'.format(out_path, first_part), 'Contigs greater than 500 bp are: {}/{}_l500_contigs.fasta'.format(out_path, first_part), logger, 'info')
    except sp.CalledProcessError:
        keep_logging('Error removing contigs greater than 500 bp using Bioawk. Exiting.', 'Error removing contigs greater than 500 bp using Bioawk. Exiting.', logger, 'exception')
        sys.exit(1)
    try:
        keep_logging(plasmid_contig_l500_cmd, plasmid_contig_l500_cmd, logger, 'debug')
        call(plasmid_contig_l500_cmd, logger)
        #print ""
        keep_logging('Plasmid Contigs greater than 500 bp are: {}/{}_l500_plasmid_contigs.fasta'.format(out_path, first_part), 'Contigs greater than 500 bp are: {}/{}_l500_plasmid_contigs.fasta'.format(out_path, first_part), logger, 'info')
    except sp.CalledProcessError:
        keep_logging('Error removing Plasmid contigs greater than 500 bp using Bioawk. Exiting.', 'Error removing Plasmid contigs greater than 500 bp using Bioawk. Exiting.', logger, 'exception')
        sys.exit(1)

    final_l500_contig = "%s/%s_l500_contigs.fasta" % (out_path, first_part)
    final_l500_plasmid_contig = "%s/%s_l500_plasmid_contigs.fasta" % (out_path, first_part)
    return final_l500_contig, final_l500_plasmid_contig