__author__ = 'alipirani'
import subprocess
import re
import os
import errno
import glob
import datetime
from config_settings import ConfigSectionMap
from modules.log_modules import *
from modules.logging_subprocess import *


def bioawk(contigs, plasmid_contigs, out_path, first_part, logger, Config, do_assembly):
    keep_logging('Removing Contigs less than 500 bp using BIOAWK', 'Removing Contigs less than 500 bp using BIOAWK', logger, 'info')
    contig_l500_cmd = "%s -c fastx '{ if(length($seq) > 500) { print \">\"$name; print $seq }}' %s > %s/%s_l500_contigs.fasta" % (ConfigSectionMap("bioawk", Config)['base_cmd'], contigs, out_path, first_part)
    plasmid_contig_l500_cmd = "%s -c fastx '{ if(length($seq) > 500) { print \">\"$name; print $seq }}' %s > %s/%s_l500_plasmid_contigs.fasta" % (ConfigSectionMap("bioawk", Config)['base_cmd'], plasmid_contigs, out_path, first_part)
    if do_assembly == "both":
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
    if do_assembly == "wga":
        try:
            keep_logging(contig_l500_cmd, contig_l500_cmd, logger, 'debug')
            call(contig_l500_cmd, logger)
            #print ""
            keep_logging('Contigs greater than 500 bp are: {}/{}_l500_contigs.fasta'.format(out_path, first_part), 'Contigs greater than 500 bp are: {}/{}_l500_contigs.fasta'.format(out_path, first_part), logger, 'info')
        except sp.CalledProcessError:
            keep_logging('Error removing contigs greater than 500 bp using Bioawk. Exiting.', 'Error removing contigs greater than 500 bp using Bioawk. Exiting.', logger, 'exception')
            sys.exit(1)
    if do_assembly == "plasmid":
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



def bioawk_make_reference_size(reference, logger, Config):
    base_cmd = ConfigSectionMap("bioawk", Config)['base_cmd']
    command = base_cmd + " -c fastx '{ print $name, length($seq) }' < %s > %s.size" % (reference, reference)
    try:
        call(command, logger)
    except sp.CalledProcessError:
        keep_logging('Error in Bioawk step. Exiting.', 'Error in Bioawk step. Exiting.', logger, 'exception')
        sys.exit(1)
    reference_size_file = "%s.size" % reference
    return reference_size_file