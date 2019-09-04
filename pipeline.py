__author__ = 'alipirani'
""" Declaring required python modules """
import argparse
import ConfigParser
import subprocess
import re
import os
import sys
import errno
import glob
import gzip
from modules.log_modules import  *
from modules.logging_subprocess import *
from config_settings import ConfigSectionMap
from modules.trimmomatic import *
from modules.assembly import *
from modules.qa_fastqc import qa_fastqc
from modules.quast import quast_evaluation
from modules.reapr import reapr
from modules.abacas import abacas
from modules.bioawk import *
from modules.abacas import *
from modules.qa_fastqc import *
from modules.prokka import prokka
from modules.ariba import ariba_AMR, ariba_MLST
from datetime import datetime



""" Command line argument parsing """
def parser():
    parser = argparse.ArgumentParser(description='Assembly pipeline for Illumina SE/PE data')
    required = parser.add_argument_group('Required arguments')
    optional = parser.add_argument_group('Optional arguments')
    required.add_argument('-f1', action='store', dest="file_1", help='Paired End file 1')
    optional.add_argument('-f2', action='store', dest="file_2", help='Paired End file 2', )
    required.add_argument('-config', action='store', dest="config", help='Path to Config file', required=True)
    required.add_argument('-analysis', action='store', dest="analysis_name", help='Unique analysis name to save the results', required=True)
    required.add_argument('-o', action='store', dest="output_folder", help='Output Path ending with output directory name to save the results', required=True)
    required.add_argument('-start_step', action="store", dest="start_step", type=int, help='Provide the start step. Only 1 works for now.')
    required.add_argument('-end_step', action="store", dest="end_step", type=int, help='Provide the end step. 2/3/4 All three works')
    required.add_argument('-A', action='store', dest="assembler", help='Choose the assembler to assemble the sample reads. Velvet Optimiser or Spades')
    required.add_argument('-type', action='store', dest='type', help='Type of analysis: SE or PE')
    optional.add_argument('-c', action='store', dest="crop", help='choose crop value to crop the reads')
    optional.add_argument('-reference', action='store', dest="reference", help='Provide a reference genome for Abacas Contig ordering')
    optional.add_argument('-ariba', action='store', dest="ariba", help='Run ariba AMR or MLST on clean reads. expected values: AMR/MLST/BOTH')
    optional.add_argument('-assembly', action='store', dest="assembly", help='Select one of the following assembly options:\n\"wga\":Only Spades Whole Genome Assembly or\n\"plasmid\": Only Plasmid Assembly or\n\"Both\": Perform both wga and plasmid assembly. Default:wga Options: wga/plasmid/both')
    return parser

# Main Pipeline
def pipeline(args, logger):
    keep_logging('START: Pipeline', 'START: Pipeline', logger, 'info')

    """ Check Subroutines and create logger object: Arguments, Input files, Reference Index """
    keep_logging('START: Checking Dependencies...', 'Checking Dependencies', logger, 'info')

    """ Check if the input file exists """
    if args.type != "PE":
        reverse_raw = "None"
        file_exists(args.file_1, reverse_raw)
    else:
        file_exists(args.file_1, args.file_2)

    """ Check java availability """
    java_check()

    """ Get Reference Genome Path """
    if args.reference and args.reference != "None":
        reference_genome_path = ConfigSectionMap(args.reference, Config)['ref_path'] + "/" + ConfigSectionMap(args.reference, Config)['ref_name']

    """ Set reads """
    if args.type != "PE":
        Reverse_read = "None"
        Forward_read = args.file_1
    else:
        Forward_read = args.file_1
        Reverse_read = args.file_2

    if args.assembly:
        do_assembly = args.assembly
    else:
        do_assembly = "wga"

    """ Start the pipeline with the given start and end steps: """
    if args.ariba == "AMR" and args.type == "PE":
        keep_logging('Running Ariba AMR...', 'Running Ariba AMR...', logger, 'info')
        ariba_AMR(Forward_read, Reverse_read, args.output_folder, args.analysis_name, logger, Config)
    if args.ariba == "MLST" and args.type == "PE":
        keep_logging('Running Ariba MLST...', 'Running Ariba MLST...', logger, 'info')
        ariba_MLST(Forward_read, Reverse_read, args.output_folder, args.analysis_name, logger, Config)
    if args.ariba == "both" and args.type == "PE":
        keep_logging('Running Ariba AMR and MLST...', 'Running Ariba AMR and MLST...', logger, 'info')
        ariba_AMR(Forward_read, Reverse_read, args.output_folder, args.analysis_name, logger, Config)
        ariba_MLST(Forward_read, Reverse_read, args.output_folder, args.analysis_name, logger, Config)


    if args.start_step and args.end_step:

        if args.start_step == 1 and args.end_step == 1:
            keep_logging('You chose steps 1 to 1: Pre-processing using Trimmomatic', 'You chose steps 1 to 1: Pre-processing using Trimmomatic', logger, 'info')
            keep_logging('START: Pre-processing step using Trimmomatic.', 'START: Pre-processing step using Trimmomatic.', logger, 'info')
            (forward_paired, reverse_paired, forward_unpaired, reverse_unpaired) = clean_reads(Forward_read, Reverse_read, args.output_folder, args.crop, logger, Config)
            keep_logging('END: Pre-processing step using Trimmomatic.', 'END: Pre-processing step using Trimmomatic.', logger, 'info')
            # if args.ariba == "AMR" and reverse_paired:
            #     keep_logging('Running Ariba AMR...', 'Running Ariba AMR...', logger, 'info')
            #     ariba_AMR(forward_paired, reverse_paired, args.output_folder, args.analysis_name, logger, Config)
            # if args.ariba == "MLST" and reverse_paired:
            #     keep_logging('Running Ariba MLST...', 'Running Ariba MLST...', logger, 'info')
            #     keep_logging('You chose steps 1 to 1: Pre-processing using Trimmomatic', 'You chose steps 1 to 1: Pre-processing using Trimmomatic', logger, 'info')
            #     ariba_MLST(forward_paired, reverse_paired, args.output_folder, args.analysis_name, logger, Config)
            # if args.ariba == "both" and reverse_paired:
            #     keep_logging('Running Ariba AMR and MLST...', 'Running Ariba AMR and MLST...', logger, 'info')
            #     keep_logging('You chose steps 1 to 1: Pre-processing using Trimmomatic', 'You chose steps 1 to 1: Pre-processing using Trimmomatic', logger, 'info')
            #     ariba_AMR(forward_paired, reverse_paired, args.output_folder, args.analysis_name, logger, Config)
            #     ariba_MLST(forward_paired, reverse_paired, args.output_folder, args.analysis_name, logger, Config)
        elif args.start_step == 1 and args.end_step == 2:
            keep_logging('You chose steps 1 to 2: Pre-processing and Assembly', 'You chose steps 1 to 2: Pre-processing and Assembly', logger, 'info')
            keep_logging('START: Pre-processing step using Trimmomatic.', 'START: Pre-processing step using Trimmomatic.', logger, 'info')
            (forward_paired, reverse_paired, forward_unpaired, reverse_unpaired) = clean_reads(Forward_read, Reverse_read, args.output_folder, args.crop, logger, Config)
            keep_logging('END: Pre-processing step using Trimmomatic.', 'END: Pre-processing step using Trimmomatic.', logger, 'info')
            keep_logging('START: Starting Assembly using {}'.format(args.assembler), 'START: Starting Assembly using {}'.format(args.assembler), logger, 'info')
            (contigs, scaffolds, plasmid_contigs, plasmid_scaffolds) = assembly(forward_paired, reverse_paired, forward_unpaired, reverse_unpaired, args.assembler, args.output_folder, logger, Config, do_assembly)
            keep_logging('END: Starting Assembly using {}'.format(args.assembler), 'END: Starting Assembly using {}'.format(args.assembler), logger, 'info')
        elif args.start_step == 1 and args.end_step == 3:
            keep_logging('You chose steps 1 to 3: Pre-processing, Assembly and Evaluation', 'You chose steps 1 to 3: Pre-processing, Assembly adn Evaluation', logger, 'info')
            keep_logging('START: Pre-processing step using Trimmomatic.', 'START: Pre-processing step using Trimmomatic.', logger, 'info')
            (forward_paired, reverse_paired, forward_unpaired, reverse_unpaired) = clean_reads(Forward_read, Reverse_read, args.output_folder, args.crop, logger, Config)
            keep_logging('END: Pre-processing step using Trimmomatic.', 'END: Pre-processing step using Trimmomatic.', logger, 'info')
            keep_logging('START: Starting Assembly using {}'.format(args.assembler), 'START: Starting Assembly using {}'.format(args.assembler), logger, 'info')
            (contigs, scaffolds, plasmid_contigs, plasmid_scaffolds) = assembly(forward_paired, reverse_paired, forward_unpaired, reverse_unpaired, args.assembler, args.output_folder, logger, Config, do_assembly)
            keep_logging('END: Starting Assembly using {}'.format(args.assembler), 'END: Starting Assembly using {}'.format(args.assembler), logger, 'info')
            (final_l500_contig, final_l500_plasmid_contig) = bioawk(contigs, plasmid_contigs, args.output_folder, args.analysis_name, logger, Config, do_assembly)
            keep_logging('START: Assembly Evaluation using QUAST', 'START: Assembly Evaluation using QUAST', logger, 'info')
            quast_evaluation(args.output_folder, final_l500_contig, final_l500_plasmid_contig, logger, Config)
            keep_logging('END: Assembly Evaluation using QUAST', 'END: Assembly Evaluation using QUAST', logger, 'info')
        elif args.start_step == 1 and args.end_step == 4:
            keep_logging('You chose steps 1 to 4: Pre-processing, Assembly, Evaluation and Annotation', 'You chose steps 1 to 4: Pre-processing, Assembly, Evaluation and Annotation', logger, 'info')
            keep_logging('START: Pre-processing step using Trimmomatic.', 'START: Pre-processing step using Trimmomatic.', logger, 'info')
            (forward_paired, reverse_paired, forward_unpaired, reverse_unpaired) = clean_reads(Forward_read, Reverse_read, args.output_folder, args.crop, logger, Config)
            keep_logging('END: Pre-processing step using Trimmomatic.', 'END: Pre-processing step using Trimmomatic.', logger, 'info')
            keep_logging('START: Starting Assembly using {}'.format(args.assembler), 'START: Starting Assembly using {}'.format(args.assembler), logger, 'info')
            (contigs, scaffolds, plasmid_contigs, plasmid_scaffolds) = assembly(forward_paired, reverse_paired, forward_unpaired, reverse_unpaired, args.assembler, args.output_folder, logger, Config, do_assembly)
            keep_logging('END: Starting Assembly using {}'.format(args.assembler), 'END: Starting Assembly using {}'.format(args.assembler), logger, 'info')
            (final_l500_contig, final_l500_plasmid_contig) = bioawk(contigs, plasmid_contigs, args.output_folder, args.analysis_name, logger, Config, do_assembly)
            keep_logging('START: Assembly Evaluation using QUAST', 'START: Assembly Evaluation using QUAST', logger, 'info')
            quast_evaluation(args.output_folder, final_l500_contig, final_l500_plasmid_contig, logger, Config)
            keep_logging('END: Assembly Evaluation using QUAST', 'END: Assembly Evaluation using QUAST', logger, 'info')
            #final_l500_contig = "%s/%s_l500_contigs.fasta" % (args.output_folder, args.analysis_name)
            #final_l500_plasmid_contig = "%s/%s_l500_plasmid_contigs.fasta" % (args.output_folder, args.analysis_name)
            if args.reference and args.reference != "None":
                final_ordered_contigs = abacas(reference_genome_path, final_l500_contig, args.output_folder, args.analysis_name, logger, Config)
                #final_ordered_contigs = "/nfs/esnitkin/Ali/Project_MRSA_analysis/MRSA_assembly///6154_R1.fastq.g_contigs_ordered.fasta"
                sed_cmd = "sed -i 's/>.*/>%s/g' %s" % (args.analysis_name[0:20], final_ordered_contigs)
                keep_logging(sed_cmd, sed_cmd, logger, 'debug')
                os.system(sed_cmd)
                final_annotation_folder = prokka(final_ordered_contigs, args.output_folder, args.analysis_name, logger, Config)
                keep_logging('Final Prokka Annotation files for genome are in: {}'.format(final_annotation_folder), 'Final Prokka Annotation files for genome are in: {}'.format(final_annotation_folder), logger, 'debug')
            else:
                print "\nPlease provide a path to reference genome for Abacas\n"
                final_ordered_contigs = final_l500_contig
                #final_ordered_contigs = "/nfs/esnitkin/Ali/Project_MRSA_analysis/MRSA_assembly///6154_R1.fastq.g_contigs_ordered.fasta"
                sed_cmd = "sed -i 's/>NODE/>%s/g' %s" % (args.analysis_name[0:15], final_ordered_contigs)
                sed_cmd_2 = "sed -i 's/_length_.*//g' %s" % (final_ordered_contigs)
                #sed_cmd_3 = "sed -i 's/>.*/>%s/g' %s" % (args.analysis_name[0:20], final_ordered_contigs)
                keep_logging(sed_cmd, sed_cmd, logger, 'debug')
                keep_logging(sed_cmd_2, sed_cmd_2, logger, 'debug')
                #keep_logging(sed_cmd_3, sed_cmd_3, logger, 'debug')
                os.system(sed_cmd)
                os.system(sed_cmd_2)
                #os.system(sed_cmd_3)
                final_annotation_folder = prokka(final_ordered_contigs, args.output_folder, args.analysis_name, logger, Config)
                keep_logging('Final Prokka Annotation files for genome are in: {}'.format(final_annotation_folder), 'Final Prokka Annotation files for genome are in: {}'.format(final_annotation_folder), logger, 'debug')
            sed_cmd = "sed -i 's/>NODE/>%s/g' %s" % (args.analysis_name[0:15], final_l500_plasmid_contig)
            #sed_cmd_2 = "sed -i 's/_length_.*component/_plasmid/g' %s" % (final_l500_plasmid_contig)
            sed_cmd_2 = "sed -i 's/_length_.*component//g' %s" % (final_l500_plasmid_contig)
            #sed_cmd_3 = "sed -i 's/>.*/>%s/g' %s" % (args.analysis_name[0:20], final_l500_plasmid_contig)
            keep_logging(sed_cmd, sed_cmd, logger, 'debug')
            keep_logging(sed_cmd_2, sed_cmd, logger, 'debug')
            os.system(sed_cmd)
            os.system(sed_cmd_2)
            #os.system(sed_cmd_3)
            plasmid_first_part = args.analysis_name + "_plasmid"
            final_plasmid_annotation_folder = prokka(final_l500_plasmid_contig, args.output_folder, plasmid_first_part, logger, Config)
            keep_logging('Final Prokka Annotation files for plasmid are in: {}'.format(final_plasmid_annotation_folder), 'Final Prokka Annotation files for plasmid are in: {}'.format(final_plasmid_annotation_folder), logger, 'debug')
        elif args.start_step == 2 and args.end_step == 3:
            keep_logging('You chose steps 2 to 3: Assembly and Evaluation', 'You chose steps 2 to 3: Assembly and Evaluation', logger, 'info')
            (contigs, scaffolds) = assembly(forward_paired, reverse_paired, forward_unpaired, reverse_unpaired, args.assembler, args.output_folder, do_assembly)
            quast_evaluation(args.output_folder, contigs, scaffolds, plasmid_contigs, plasmid_scaffolds)

        elif args.start_step == 2 and args.end_step == 2:
            keep_logging('START: Starting Assembly using {}'.format(args.assembler), 'START: Starting Assembly using {}'.format(args.assembler), logger, 'info')
            (contigs, scaffolds, plasmid_contigs, plasmid_scaffolds) = assembly(forward_paired, reverse_paired, forward_unpaired, reverse_unpaired, args.assembler, args.output_folder, logger, Config, do_assembly)
            keep_logging('END: Starting Assembly using {}'.format(args.assembler), 'END: Starting Assembly using {}'.format(args.assembler), logger, 'info')

        elif args.start_step == 2 and args.end_step == 4:
            keep_logging('START: Starting Assembly using {}'.format(args.assembler), 'START: Starting Assembly using {}'.format(args.assembler), logger, 'info')
            (contigs, scaffolds, plasmid_contigs, plasmid_scaffolds) = assembly(forward_paired, reverse_paired, forward_unpaired, reverse_unpaired, args.assembler, args.output_folder, logger, Config, do_assembly)
            keep_logging('END: Starting Assembly using {}'.format(args.assembler), 'END: Starting Assembly using {}'.format(args.assembler), logger, 'info')
            (final_l500_contig, final_l500_plasmid_contig) = bioawk(contigs, plasmid_contigs, args.output_folder, args.analysis_name, logger, Config, do_assembly)
            keep_logging('START: Assembly Evaluation using QUAST', 'START: Assembly Evaluation using QUAST', logger, 'info')
            quast_evaluation(args.output_folder, final_l500_contig, final_l500_plasmid_contig, logger, Config)
            keep_logging('END: Assembly Evaluation using QUAST', 'END: Assembly Evaluation using QUAST', logger, 'info')
            #final_l500_contig = "%s/%s_l500_contigs.fasta" % (args.output_folder, args.analysis_name)
            #final_l500_plasmid_contig = "%s/%s_l500_plasmid_contigs.fasta" % (args.output_folder, args.analysis_name)
            if args.reference and args.reference != "None":
                final_ordered_contigs = abacas(reference_genome_path, final_l500_contig, args.output_folder, args.analysis_name, logger, Config)
                #final_ordered_contigs = "/nfs/esnitkin/Ali/Project_MRSA_analysis/MRSA_assembly///6154_R1.fastq.g_contigs_ordered.fasta"
                sed_cmd = "sed -i 's/>.*/>%s/g' %s" % (args.analysis_name[0:20], final_ordered_contigs)
                keep_logging(sed_cmd, sed_cmd, logger, 'debug')
                os.system(sed_cmd)
                final_annotation_folder = prokka(final_ordered_contigs, args.output_folder, args.analysis_name, logger, Config)
                keep_logging('Final Prokka Annotation files for genome are in: {}'.format(final_annotation_folder), 'Final Prokka Annotation files for genome are in: {}'.format(final_annotation_folder), logger, 'debug')
            else:
                print "\nPlease provide a path to reference genome for Abacas\n"
                final_ordered_contigs = final_l500_contig
                #final_ordered_contigs = "/nfs/esnitkin/Ali/Project_MRSA_analysis/MRSA_assembly///6154_R1.fastq.g_contigs_ordered.fasta"
                sed_cmd = "sed -i 's/>NODE/>%s/g' %s" % (args.analysis_name[0:15], final_ordered_contigs)
                sed_cmd_2 = "sed -i 's/_length_.*//g' %s" % (final_ordered_contigs)
                #sed_cmd_3 = "sed -i 's/>.*/>%s/g' %s" % (args.analysis_name[0:20], final_ordered_contigs)
                keep_logging(sed_cmd, sed_cmd, logger, 'debug')
                keep_logging(sed_cmd_2, sed_cmd_2, logger, 'debug')
                #keep_logging(sed_cmd_3, sed_cmd_3, logger, 'debug')
                os.system(sed_cmd)
                os.system(sed_cmd_2)
                #os.system(sed_cmd_3)
                final_annotation_folder = prokka(final_ordered_contigs, args.output_folder, args.analysis_name, logger, Config)
                keep_logging('Final Prokka Annotation files for genome are in: {}'.format(final_annotation_folder), 'Final Prokka Annotation files for genome are in: {}'.format(final_annotation_folder), logger, 'debug')
            sed_cmd = "sed -i 's/>NODE/>%s/g' %s" % (args.analysis_name[0:15], final_l500_plasmid_contig)
            #sed_cmd_2 = "sed -i 's/_length_.*component/_plasmid/g' %s" % (final_l500_plasmid_contig)
            sed_cmd_2 = "sed -i 's/_length_.*component//g' %s" % (final_l500_plasmid_contig)
            #sed_cmd_3 = "sed -i 's/>.*/>%s/g' %s" % (args.analysis_name[0:20], final_l500_plasmid_contig)
            keep_logging(sed_cmd, sed_cmd, logger, 'debug')
            keep_logging(sed_cmd_2, sed_cmd, logger, 'debug')
            os.system(sed_cmd)
            os.system(sed_cmd_2)
            #os.system(sed_cmd_3)
            plasmid_first_part = args.analysis_name + "_plasmid"
            final_plasmid_annotation_folder = prokka(final_l500_plasmid_contig, args.output_folder, plasmid_first_part, logger, Config)
            keep_logging('Final Prokka Annotation files for plasmid are in: {}'.format(final_plasmid_annotation_folder), 'Final Prokka Annotation files for plasmid are in: {}'.format(final_plasmid_annotation_folder), logger, 'debug')



        elif args.start_step == 3 and args.end_step == 3:
            keep_logging('START: Assembly Evaluation using QUAST', 'START: Assembly Evaluation using QUAST', logger, 'info')
            quast_evaluation(args.output_folder, final_l500_contig, final_l500_plasmid_contig, logger, Config)
            keep_logging('END: Assembly Evaluation using QUAST', 'END: Assembly Evaluation using QUAST', logger, 'info')

        elif args.start_step == 3 and args.end_step == 4:
            keep_logging('START: Assembly Evaluation using QUAST', 'START: Assembly Evaluation using QUAST', logger, 'info')
            quast_evaluation(args.output_folder, final_l500_contig, final_l500_plasmid_contig, logger, Config)
            keep_logging('END: Assembly Evaluation using QUAST', 'END: Assembly Evaluation using QUAST', logger, 'info')
            #final_l500_contig = "%s/%s_l500_contigs.fasta" % (args.output_folder, args.analysis_name)
            #final_l500_plasmid_contig = "%s/%s_l500_plasmid_contigs.fasta" % (args.output_folder, args.analysis_name)
            if args.reference and args.reference != "None":
                final_ordered_contigs = abacas(reference_genome_path, final_l500_contig, args.output_folder, args.analysis_name, logger, Config)
                #final_ordered_contigs = "/nfs/esnitkin/Ali/Project_MRSA_analysis/MRSA_assembly///6154_R1.fastq.g_contigs_ordered.fasta"
                sed_cmd = "sed -i 's/>.*/>%s/g' %s" % (args.analysis_name[0:20], final_ordered_contigs)
                keep_logging(sed_cmd, sed_cmd, logger, 'debug')
                os.system(sed_cmd)
                final_annotation_folder = prokka(final_ordered_contigs, args.output_folder, args.analysis_name, logger, Config)
                keep_logging('Final Prokka Annotation files for genome are in: {}'.format(final_annotation_folder), 'Final Prokka Annotation files for genome are in: {}'.format(final_annotation_folder), logger, 'debug')
            else:
                print "\nPlease provide a path to reference genome for Abacas\n"
                final_ordered_contigs = final_l500_contig
                #final_ordered_contigs = "/nfs/esnitkin/Ali/Project_MRSA_analysis/MRSA_assembly///6154_R1.fastq.g_contigs_ordered.fasta"
                sed_cmd = "sed -i 's/>NODE/>%s/g' %s" % (args.analysis_name[0:15], final_ordered_contigs)
                sed_cmd_2 = "sed -i 's/_length_.*//g' %s" % (final_ordered_contigs)
                #sed_cmd_3 = "sed -i 's/>.*/>%s/g' %s" % (args.analysis_name[0:20], final_ordered_contigs)
                keep_logging(sed_cmd, sed_cmd, logger, 'debug')
                keep_logging(sed_cmd_2, sed_cmd_2, logger, 'debug')
                #keep_logging(sed_cmd_3, sed_cmd_3, logger, 'debug')
                os.system(sed_cmd)
                os.system(sed_cmd_2)
                #os.system(sed_cmd_3)
                final_annotation_folder = prokka(final_ordered_contigs, args.output_folder, args.analysis_name, logger, Config)
                keep_logging('Final Prokka Annotation files for genome are in: {}'.format(final_annotation_folder), 'Final Prokka Annotation files for genome are in: {}'.format(final_annotation_folder), logger, 'debug')
            sed_cmd = "sed -i 's/>NODE/>%s/g' %s" % (args.analysis_name[0:15], final_l500_plasmid_contig)
            #sed_cmd_2 = "sed -i 's/_length_.*component/_plasmid/g' %s" % (final_l500_plasmid_contig)
            sed_cmd_2 = "sed -i 's/_length_.*component//g' %s" % (final_l500_plasmid_contig)
            #sed_cmd_3 = "sed -i 's/>.*/>%s/g' %s" % (args.analysis_name[0:20], final_l500_plasmid_contig)
            keep_logging(sed_cmd, sed_cmd, logger, 'debug')
            keep_logging(sed_cmd_2, sed_cmd, logger, 'debug')
            os.system(sed_cmd)
            os.system(sed_cmd_2)
            #os.system(sed_cmd_3)
            plasmid_first_part = args.analysis_name + "_plasmid"
            final_plasmid_annotation_folder = prokka(final_l500_plasmid_contig, args.output_folder, plasmid_first_part, logger, Config)
            keep_logging('Final Prokka Annotation files for plasmid are in: {}'.format(final_plasmid_annotation_folder), 'Final Prokka Annotation files for plasmid are in: {}'.format(final_plasmid_annotation_folder), logger, 'debug')



        elif args.start_step == 4 and args.end_step == 4:
            final_l500_contig = "%s/%s_l500_contigs.fasta" % (args.output_folder, args.analysis_name)
            final_l500_plasmid_contig = "%s/%s_l500_plasmid_contigs.fasta" % (args.output_folder, args.analysis_name)
            if args.reference and args.reference != "None":
                print "here"
                final_ordered_contigs = abacas(reference_genome_path, final_l500_contig, args.output_folder, args.analysis_name, logger, Config)
                #final_ordered_contigs = "/nfs/esnitkin/Ali/Project_MRSA_analysis/MRSA_assembly///6154_R1.fastq.g_contigs_ordered.fasta"
                sed_cmd = "sed -i 's/>.*/>%s/g' %s" % (args.analysis_name[0:20], final_ordered_contigs)
                keep_logging(sed_cmd, sed_cmd, logger, 'debug')
                os.system(sed_cmd)
                final_annotation_folder = prokka(final_ordered_contigs, args.output_folder, args.analysis_name, logger, Config)
                keep_logging('Final Prokka Annotation files for genome are in: {}'.format(final_annotation_folder), 'Final Prokka Annotation files for genome are in: {}'.format(final_annotation_folder), logger, 'debug')
            else:
                print "\nPlease provide a path to reference genome for Abacas\n"
                final_ordered_contigs = final_l500_contig
                #final_ordered_contigs = "/nfs/esnitkin/Ali/Project_MRSA_analysis/MRSA_assembly///6154_R1.fastq.g_contigs_ordered.fasta"
                sed_cmd = "sed -i 's/>NODE/>%s/g' %s" % (args.analysis_name[0:15], final_ordered_contigs)
                sed_cmd_2 = "sed -i 's/_length_.*//g' %s" % (final_ordered_contigs)
                #sed_cmd_3 = "sed -i 's/>.*/>%s/g' %s" % (args.analysis_name[0:20], final_ordered_contigs)
                keep_logging(sed_cmd, sed_cmd, logger, 'debug')
                keep_logging(sed_cmd_2, sed_cmd_2, logger, 'debug')
                #keep_logging(sed_cmd_3, sed_cmd_3, logger, 'debug')
                os.system(sed_cmd)
                os.system(sed_cmd_2)
                #os.system(sed_cmd_3)
                final_annotation_folder = prokka(final_ordered_contigs, args.output_folder, args.analysis_name, logger, Config)
                keep_logging('Final Prokka Annotation files for genome are in: {}'.format(final_annotation_folder), 'Final Prokka Annotation files for genome are in: {}'.format(final_annotation_folder), logger, 'debug')
            sed_cmd = "sed -i 's/>NODE/>%s/g' %s" % (args.analysis_name[0:15], final_l500_plasmid_contig)
            #sed_cmd_2 = "sed -i 's/_length_.*component/_plasmid/g' %s" % (final_l500_plasmid_contig)
            sed_cmd_2 = "sed -i 's/_length_.*component//g' %s" % (final_l500_plasmid_contig)
            #sed_cmd_3 = "sed -i 's/>.*/>%s/g' %s" % (args.analysis_name[0:20], final_l500_plasmid_contig)
            keep_logging(sed_cmd, sed_cmd, logger, 'debug')
            keep_logging(sed_cmd_2, sed_cmd, logger, 'debug')
            os.system(sed_cmd)
            os.system(sed_cmd_2)
            #os.system(sed_cmd_3)
            plasmid_first_part = args.analysis_name + "_plasmid"
            final_plasmid_annotation_folder = prokka(final_l500_plasmid_contig, args.output_folder, plasmid_first_part, logger, Config)
            keep_logging('Final Prokka Annotation files for plasmid are in: {}'.format(final_plasmid_annotation_folder), 'Final Prokka Annotation files for plasmid are in: {}'.format(final_plasmid_annotation_folder), logger, 'debug')

        else:
            keep_logging('ERROR: Please provide start and end steps for the pipeline.', 'ERROR: Please provide start and end steps for the pipeline.', logger, 'exception')


""" Start: Check Subroutines """

""" Usage Message: """
def usage():
    print "Usage: pipeline.py [-h] [-v] [-f1 FILE_1] [-f2 FILE_2] [-o OUTPUT_FOLDER] [-start_step START] [-end_step END] [--qa] [-A Assembler] [--err] \n"


""" Check Java Availability: """
def java_check():
    print "Checking Java Availability....\n"
    jd = subprocess.check_output(["java", "-version"], stderr=subprocess.STDOUT)
    if len(jd) < 1:
        print "Unable to find a java runtime environment. The pipeline requires java 6 or later."
    else:
        print "Java Availability Check completed ...\n\n" + jd


""" Make sure input raw reads files exists at given location. """
def file_exists(path1, path2):
    if path1:
        if not os.path.isfile(path1):
            file_basename = os.path.basename(path1)
            print "The input file " + file_basename + " does not exists. \nPlease provide another file.\n"
            exit()
    if path2 != "None":
        if not os.path.isfile(path2):
            file_basename = os.path.basename(path2)
            print "The input file " + file_basename + " does not exists. \nPlease provide another file.\n"
            exit()


""" Make sure the output folder exists or create at given path """
def make_sure_path_exists(out_path):
    try:
        os.makedirs(out_path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            print "Errors in output folder path! please change the output path or analysis name\n"
            exit()


""" Prepare the paths for clean reads when all the input files exists. """
def prepare_cleanreadsinput():
    global forward_paired, forward_unpaired, reverse_unpaired, reverse_paired
    forward_paired = out_name + "forward_paired.fq.gz"
    forward_unpaired = out_name + "forward_unpaired.fq.gz"
    reverse_paired = out_name + "reverse_paired.fq.gz"
    reverse_unpaired = out_name + "reverse_unpaired.fq.gz"
    return forward_paired, forward_unpaired, reverse_paired, reverse_unpaired


""" Set appropriate path for contigs/scaffolds fasta files. Needed only when quast or reapr runs individually. """
def get_contigs(out_path, assembler):
    if assembler == "velvet":
        contigs = out_path + ConfigSectionMap("velvet")['contigs_path']
        scaffolds = ""
        return contigs, scaffolds
    elif assembler == "spades":
        contigs = out_path + ConfigSectionMap("spades")['contigs_path']
        scaffolds = out_path + ConfigSectionMap("spades")['scaffolds_path']
        plasmid_contigs = out_path + ConfigSectionMap("spades")['plasmid_contigs_path']
        plasmid_scaffolds = out_path + ConfigSectionMap("spades")['plasmid_scaffolds_path']
        return contigs, scaffolds, plasmid_contigs, plasmid_scaffolds

def prepare_readgroup(forward_read, aligner, logger):
    keep_logging('Preparing ReadGroup Info', 'Preparing ReadGroup Info', logger, 'info')
    samplename = os.path.basename(forward_read)
    if forward_read.endswith(".gz"):
        #output = gzip.open(forward_read, 'rb')
        #firstLine = output.readline()
        #split_field = re.split(r":",firstLine)
        #id_name = split_field[1]
        #id_name = id_name.strip()
        #split_field = "\"" + "@RG" + "\\tID:" + split_field[1] + "\\tSM:" + samplename + "\\tLB:1\\tPL:Illumina" + "\""
        #return split_field
        output = gzip.open(forward_read, 'rb')
        firstLine = output.readline()
        if ":" in firstLine:
            split_field = re.split(r":",firstLine)
            id_name = split_field[1].rstrip()
            id_name = id_name.rstrip()
        # if aligner == "bowtie":
        #     split_field = "--rg-id %s --rg SM:%s --rg LB:1 --rg PL:Illumina" % (split_field[1], samplename)
        # elif aligner == "bwa":
        #     split_field = "\"" + "@RG" + "\\tID:" + split_field[1] + "\\tSM:" + samplename + "\\tLB:1\\tPL:Illumina" + "\""

        ###Pending
        elif "/" in firstLine:
            split_field = re.split(r"/",firstLine)
            id_name = split_field[1].rstrip()
            id_name = id_name.rstrip()
            #id_name = split_field[1].rstrip()
            #id_name = id_name.rstrip()
            split_field = "\"" + "@RG" + "\\tID:" + id_name + "\\tSM:" + samplename + "\\tLB:1\\tPL:Illumina" + "\""
        else:
            id_name = "1"
            split_field = "\"" + "@RG" + "\\tID:" + id_name + "\\tSM:" + samplename + "\\tLB:1\\tPL:Illumina" + "\""
        if aligner == "bowtie":
            split_field = "--rg-id %s --rg SM:%s --rg LB:1 --rg PL:Illumina" % (split_field[1], samplename)
        elif aligner == "bwa":
            split_field = "\"" + "@RG" + "\\tID:" + split_field[1] + "\\tSM:" + samplename + "\\tLB:1\\tPL:Illumina" + "\""
        return split_field

    elif forward_read.endswith(".fastq"):
        output = open(forward_read, 'r')
        firstLine = output.readline()
        split_field = re.split(r":",firstLine)
        split_field = "\"" + "@RG" + "\\tID:" + split_field[1] + "\\tSM:" + samplename + "\\tLB:1\\tPL:Illumina" + "\""
        return split_field

    elif forward_read.endswith(".fq"):
        ###
        output = open(forward_read, 'r')
        firstLine = output.readline()
        split_field = re.split(r":",firstLine)
        split_field = "\"" + "@RG" + "\\tID:" + split_field[1] + "\\tSM:" + samplename + "\\tLB:1\\tPL:Illumina" + "\""
        return split_field
""" End: Check Subroutines """



"""
Main Method:
1. parse the command-line arguments.
2. Create and make sure the output folder exist.
3. Generate the logging method object logger to initiate log file.
4. start pipeline
"""
if __name__ == '__main__':
    start_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    start_time_2 = datetime.now()
    args = parser().parse_args()
    global config_file
    if args.config:
        config_file = args.config
    else:
        config_file = "./config"
    if args.output_folder != '':
        args.output_folder += '/'
        make_sure_path_exists(args.output_folder)
    global logger
    log_unique_time = datetime.now().strftime('%Y_%m_%d_%H_%M_%S')
    logger = generate_logger(args.output_folder, args.analysis_name, log_unique_time)
    global Config
    Config = ConfigParser.ConfigParser()
    Config.read(config_file)
    pipeline(args, logger)
    # final_l500_contig = "%s/%s_l500_contigs.fasta" % (args.output_folder, args.analysis_name)
    # forward_paired = args.output_folder + ConfigSectionMap("Trimmomatic", Config)['f_p']
    # reverse_paired = args.output_folder + ConfigSectionMap("Trimmomatic", Config)['r_p']
    # forward_unpaired = args.output_folder + ConfigSectionMap("Trimmomatic", Config)['f_up']
    # reverse_unpaired = args.output_folder + ConfigSectionMap("Trimmomatic", Config)['r_up']
    # split_field = prepare_readgroup(forward_paired, "bwa", logger)
    # base_cmd = ConfigSectionMap("bin_path", Config)['binbase'] + "/" + ConfigSectionMap("bowtie", Config)[
    #     'bowtie_bin'] + "/" + ConfigSectionMap("bowtie", Config)['align_cmd']
    # files_to_delete = ""
    # parameters = ConfigSectionMap("bowtie", Config)['parameters']
    #
    # align_bowtie(base_cmd, forward_paired, reverse_paired, forward_unpaired, reverse_unpaired, args.output_folder, final_l500_contig,
    #              split_field, args.analysis_name, files_to_delete, logger, Config, args.type, parameters)
    keep_logging('End: Pipeline', 'End: Pipeline', logger, 'info')
    time_taken = datetime.now() - start_time_2
    keep_logging('Total Time taken: {}'.format(time_taken), 'Total Time taken: {}'.format(time_taken), logger, 'info')





