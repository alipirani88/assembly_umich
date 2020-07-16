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


####################################################################### Post-Assembly processing using Abacas #######################################################################################################
def abacas(reference_genome_path, final_l500_contig, out_path, first_part, logger, Config):
    keep_logging('Contig Reordering using ABACAS', 'Contig Reordering using ABACAS', logger, 'info')
    abacas_cmd = "perl %s/%s/%s -r %s -q %s %s -o %s/%s_contigs_ordered" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("abacas", Config)['abacas_bin'], ConfigSectionMap("abacas", Config)['base_cmd'], reference_genome_path, final_l500_contig, ConfigSectionMap("abacas", Config)['abacas_parameters'], out_path, first_part)
    try:
        keep_logging(abacas_cmd, abacas_cmd, logger, 'debug')
        call(abacas_cmd, logger)
        #print ""
        fasta_header = ">%s" % first_part
        header_cmd = "echo \"%s\" > %s/fasta_header" % (fasta_header, out_path)
        print (header_cmd)
        keep_logging(abacas_cmd, abacas_cmd, logger, 'debug')
        call(header_cmd, logger)
        abacas_ordered_multifasta = "%s/%s_contigs_ordered.MULTIFASTA.fa" % (out_path, first_part)
        abacas_ordered_contigsInbin = "%s/%s_contigs_ordered.contigsInbin.fas" % (out_path, first_part)
        join_all_contigs = "cat %s %s > %s/all_contigs.fasta" % (abacas_ordered_multifasta, abacas_ordered_contigsInbin, out_path)
        #print join_all_contigs
        keep_logging(join_all_contigs, join_all_contigs, logger, 'debug')
        call(join_all_contigs, logger)
        add_linker = "sed -i 's/>.*/NNNNNCATTCCATTCATTAATTAATTAATGAATGAATGNNNNN/g' %s/all_contigs.fasta" % out_path
        #print add_linker
        keep_logging(add_linker, add_linker, logger, 'debug')
        call(add_linker, logger)
        remove_spaces = "tr -d '[:space:]' < %s/all_contigs.fasta > %s/all_contigs.fasta_changed.fasta" % (out_path, out_path)
        #print remove_spaces
        keep_logging(remove_spaces, remove_spaces, logger, 'debug')
        call(remove_spaces, logger)
        join_files = "cat %s/fasta_header %s/all_contigs.fasta_changed.fasta > %s/%s_contigs_ordered.fasta" % (out_path, out_path, out_path, first_part)
        #print join_files
        keep_logging(join_files, join_files, logger, 'debug')
        call(join_files, logger)
    except sp.CalledProcessError:
        keep_logging('Error in reordering Contigs using Abacas. Exiting.', 'Error in reordering Contigs using Abacas. Exiting.', logger, 'exception')
        sys.exit(1)
    # fasta_header = ">%s" % first_part
    # header_cmd = "echo \"%s\" > %s/fasta_header" % (fasta_header, out_path)
    # print header_cmd
    # call(header_cmd, logger)
    # abacas_ordered_multifasta = "%s/%s_contigs_ordered.MULTIFASTA.fa" % (out_path, first_part)
    # abacas_ordered_contigsInbin = "%s/%s_contigs_ordered.contigsInbin.fas" % (out_path, first_part)
    # print "here"
    # join_all_contigs = "cat %s %s > %s/all_contigs.fasta" % (abacas_ordered_multifasta, abacas_ordered_contigsInbin, out_path)
    # print join_all_contigs
    # call(join_all_contigs, logger)
    # add_linker = "sed -i 's/>.*/NNNNNCATTCCATTCATTAATTAATTAATGAATGAATGNNNNN/g' %s/all_contigs.fasta" % out_path
    # print add_linker
    # call(add_linker, logger)
    # remove_spaces = "tr -d '[:space:]' < %s/all_contigs.fasta > %s/all_contigs.fasta_changed.fasta" % (out_path, out_path)
    # print remove_spaces
    # call(remove_spaces, logger)
    # join_files = "cat %s/fasta_header %s/all_contigs.fasta_changed.fasta > %s/%s_contigs_ordered.fasta" % (out_path, out_path, out_path, first_part)
    # print join_files
    # call(join_files, logger)
    final_ordered_contigs = "%s/%s_contigs_ordered.fasta" % (out_path, first_part)
    return final_ordered_contigs