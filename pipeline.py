__author__ = 'alipirani'
""" Declaring required python modules """
import argparse
import subprocess as sp
import re
import os
import sys
import errno
import glob
import gzip
from datetime import datetime
import configparser
from config_settings import ConfigSectionMap
from datetime import *
from modules.log_modules import  *
from modules.logging_subprocess import *
from modules.trimmomatic import *
from modules.assembly import *
from modules.qa_fastqc import qa_fastqc
from modules.quast import quast_evaluation
from modules.abacas import abacas
from modules.bioawk import *
from modules.abacas import *
from modules.qa_fastqc import *
from modules.prokka import prokka
from modules.ariba import ariba_AMR, ariba_MLST
from modules.bowtie import *
from modules.samtools import *
from modules.picard import *
from modules.bedtools import *
from modules.pilon import *

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
    optional.add_argument('-downsample', action='store', dest="downsample",
                          help='yes/no: Downsample Reads data to default depth of 100X or user specified depth')
    optional.add_argument('-coverage_depth', action='store', dest="coverage_depth",
                          help='Downsample Reads to this user specified depth')
    optional.add_argument('-genome_size', action='store', dest="genome_size",
                          help='Genome Size. If not provided, will be estimated from Mash')
<<<<<<< HEAD
    optional.add_argument('-pilon', action='store', dest="pilon", help='Run pilon yes/no. Default - yes')
=======
>>>>>>> 09565177f3a85722989ef80167ec09f603f993ae
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

    """ Downsample original data to default/ user-specified coverage"""
    if args.downsample == "yes":
        read1, read2 = downsample(args, logger)
        Forward_read = read1
        Reverse_read = read2
        print ("Using downsampled forward reads %s" % Forward_read)
        print ("Using downsampled reverse reads %s" % Reverse_read)




    if args.start_step and args.end_step:

        if args.start_step == 1 and args.end_step == 1:
            keep_logging('You chose steps 1 to 1: Pre-processing using Trimmomatic', 'You chose steps 1 to 1: Pre-processing using Trimmomatic', logger, 'info')
            keep_logging('START: Pre-processing step using Trimmomatic.', 'START: Pre-processing step using Trimmomatic.', logger, 'info')
            (forward_paired, reverse_paired, forward_unpaired, reverse_unpaired) = clean_reads(Forward_read, Reverse_read, args.output_folder, args.crop, logger, Config)
            keep_logging('END: Pre-processing step using Trimmomatic.', 'END: Pre-processing step using Trimmomatic.', logger, 'info')
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

        # Main Method to run end-to-end Assembly pipeline
        elif args.start_step == 1 and args.end_step == 4:
            keep_logging('You chose steps 1 to 4: Pre-processing, Assembly, Evaluation and Annotation', 'You chose steps 1 to 4: Pre-processing, Assembly, Evaluation and Annotation', logger, 'info')
            keep_logging('START: Pre-processing step using Trimmomatic.', 'START: Pre-processing step using Trimmomatic.', logger, 'info')

            forward_paired = args.output_folder + ConfigSectionMap("Trimmomatic", Config)['f_p']
            reverse_paired = args.output_folder + ConfigSectionMap("Trimmomatic", Config)['r_p']
            forward_unpaired = args.output_folder + ConfigSectionMap("Trimmomatic", Config)['f_up']
            reverse_unpaired = args.output_folder + ConfigSectionMap("Trimmomatic", Config)['r_up']

            (forward_paired, reverse_paired, forward_unpaired, reverse_unpaired) = clean_reads(Forward_read, Reverse_read, args.output_folder, args.crop, logger, Config)
            keep_logging('END: Pre-processing step using Trimmomatic.', 'END: Pre-processing step using Trimmomatic.', logger, 'info')
            keep_logging('START: Starting Assembly using {}'.format(args.assembler), 'START: Starting Assembly using {}'.format(args.assembler), logger, 'info')
            (contigs, scaffolds, plasmid_contigs, plasmid_scaffolds) = assembly(forward_paired, reverse_paired, forward_unpaired, reverse_unpaired, args.assembler, args.output_folder, logger, Config, do_assembly)
            keep_logging('END: Starting Assembly using {}'.format(args.assembler), 'END: Starting Assembly using {}'.format(args.assembler), logger, 'info')

            forward_paired = args.output_folder + ConfigSectionMap("Trimmomatic", Config)['f_p']
            reverse_paired = args.output_folder + ConfigSectionMap("Trimmomatic", Config)['r_p']
            forward_unpaired = args.output_folder + ConfigSectionMap("Trimmomatic", Config)['f_up']
            reverse_unpaired = args.output_folder + ConfigSectionMap("Trimmomatic", Config)['r_up']

            # Use Spades assembly as a reference to map back reads and generate BAM for Pilon
            reference = "%s/spades_results/contigs.fasta" % (args.output_folder)
            reference_plasmid = "%s/spades_plasmid_results/contigs.fasta" % (args.output_folder)

            if args.pilon and args.pilon == "no":
                # Add bioawk here
                (final_l500_contig, final_l500_plasmid_contig) = bioawk(reference, reference_plasmid,
                                                                        args.output_folder,
                                                                        args.analysis_name, logger, Config,
                                                                        args.assembly)
            else:
                # Adding Pilon to the pipeline
                (final_l500_contig, final_l500_plasmid_contig) = run_pilon(forward_paired, reverse_paired, forward_unpaired, reverse_unpaired, reference)

<<<<<<< HEAD
=======
            # Adding Pilon to the pipeline
            (final_l500_contig, final_l500_plasmid_contig) = run_pilon(forward_paired, reverse_paired, forward_unpaired, reverse_unpaired, reference)

>>>>>>> 09565177f3a85722989ef80167ec09f603f993ae
            final_l500_contig = "%s/%s_l500_contigs.fasta" % (args.output_folder, args.analysis_name)
            final_l500_plasmid_contig = "%s/%s_l500_plasmid_contigs.fasta" % (args.output_folder, args.analysis_name)

            # Order contigs with reference to the reference genome provided.
            if args.reference and args.reference != "None":
                final_ordered_contigs = abacas(reference_genome_path, final_l500_contig, args.output_folder, args.analysis_name, logger, Config)
                sed_cmd = "sed -i 's/>.*/>%s/g' %s" % (args.analysis_name[0:20], final_ordered_contigs)
                keep_logging(sed_cmd, sed_cmd, logger, 'debug')
                os.system(sed_cmd)
                final_annotation_folder = prokka(final_ordered_contigs, args.output_folder, args.analysis_name, logger, Config)
                keep_logging('Final Prokka Annotation files for genome are in: {}'.format(final_annotation_folder), 'Final Prokka Annotation files for genome are in: {}'.format(final_annotation_folder), logger, 'debug')
            else:
                print ("\nSkipping Abacas...\n")
                final_ordered_contigs = final_l500_contig
                sed_cmd = "sed -i 's/>NODE/>%s/g' %s" % (args.analysis_name[0:15], final_ordered_contigs)
                sed_cmd_2 = "sed -i 's/_length_.*//g' %s" % (final_ordered_contigs)
                keep_logging(sed_cmd, sed_cmd, logger, 'debug')
                keep_logging(sed_cmd_2, sed_cmd_2, logger, 'debug')
                os.system(sed_cmd)
                os.system(sed_cmd_2)
                final_annotation_folder = prokka(final_ordered_contigs, args.output_folder, args.analysis_name, logger, Config)
                keep_logging('Final Prokka Annotation files for genome are in: {}'.format(final_annotation_folder), 'Final Prokka Annotation files for genome are in: {}'.format(final_annotation_folder), logger, 'debug')
            sed_cmd = "sed -i 's/>NODE/>%s/g' %s" % (args.analysis_name[0:15], final_l500_plasmid_contig)
            sed_cmd_2 = "sed -i 's/_length_.*component//g' %s" % (final_l500_plasmid_contig)
            keep_logging(sed_cmd, sed_cmd, logger, 'debug')
            keep_logging(sed_cmd_2, sed_cmd, logger, 'debug')
            os.system(sed_cmd)
            os.system(sed_cmd_2)
            plasmid_first_part = args.analysis_name + "_plasmid"
            final_plasmid_annotation_folder = prokka(final_l500_plasmid_contig, args.output_folder, plasmid_first_part, logger, Config)
            keep_logging('Final Prokka Annotation files for plasmid are in: {}'.format(final_plasmid_annotation_folder), 'Final Prokka Annotation files for plasmid are in: {}'.format(final_plasmid_annotation_folder), logger, 'debug')

            keep_logging('START: Assembly Evaluation using QUAST', 'START: Assembly Evaluation using QUAST', logger,
                         'info')
            quast_evaluation(args.output_folder, final_l500_contig, final_l500_plasmid_contig, logger, Config)
            keep_logging('END: Assembly Evaluation using QUAST', 'END: Assembly Evaluation using QUAST', logger, 'info')

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

            if args.reference and args.reference != "None":
                final_ordered_contigs = abacas(reference_genome_path, final_l500_contig, args.output_folder, args.analysis_name, logger, Config)
                #final_ordered_contigs = "/nfs/esnitkin/Ali/Project_MRSA_analysis/MRSA_assembly///6154_R1.fastq.g_contigs_ordered.fasta"
                sed_cmd = "sed -i 's/>.*/>%s/g' %s" % (args.analysis_name[0:20], final_ordered_contigs)
                keep_logging(sed_cmd, sed_cmd, logger, 'debug')
                os.system(sed_cmd)
                final_annotation_folder = prokka(final_ordered_contigs, args.output_folder, args.analysis_name, logger, Config)
                keep_logging('Final Prokka Annotation files for genome are in: {}'.format(final_annotation_folder), 'Final Prokka Annotation files for genome are in: {}'.format(final_annotation_folder), logger, 'debug')
            else:
                print ("\nPlease provide a path to reference genome for Abacas\n")
                final_ordered_contigs = final_l500_contig

                sed_cmd = "sed -i 's/>NODE/>%s/g' %s" % (args.analysis_name[0:15], final_ordered_contigs)
                sed_cmd_2 = "sed -i 's/_length_.*//g' %s" % (final_ordered_contigs)

                keep_logging(sed_cmd, sed_cmd, logger, 'debug')
                keep_logging(sed_cmd_2, sed_cmd_2, logger, 'debug')

                os.system(sed_cmd)
                os.system(sed_cmd_2)

                final_annotation_folder = prokka(final_ordered_contigs, args.output_folder, args.analysis_name, logger, Config)
                keep_logging('Final Prokka Annotation files for genome are in: {}'.format(final_annotation_folder), 'Final Prokka Annotation files for genome are in: {}'.format(final_annotation_folder), logger, 'debug')
            sed_cmd = "sed -i 's/>NODE/>%s/g' %s" % (args.analysis_name[0:15], final_l500_plasmid_contig)

            sed_cmd_2 = "sed -i 's/_length_.*component//g' %s" % (final_l500_plasmid_contig)

            keep_logging(sed_cmd, sed_cmd, logger, 'debug')
            keep_logging(sed_cmd_2, sed_cmd, logger, 'debug')
            os.system(sed_cmd)
            os.system(sed_cmd_2)

            plasmid_first_part = args.analysis_name + "_plasmid"
            final_plasmid_annotation_folder = prokka(final_l500_plasmid_contig, args.output_folder, plasmid_first_part, logger, Config)
            keep_logging('Final Prokka Annotation files for plasmid are in: {}'.format(final_plasmid_annotation_folder), 'Final Prokka Annotation files for plasmid are in: {}'.format(final_plasmid_annotation_folder), logger, 'debug')

            keep_logging('START: Assembly Evaluation using QUAST', 'START: Assembly Evaluation using QUAST', logger,
                         'info')
<<<<<<< HEAD
=======
            quast_evaluation(args.output_folder, final_l500_contig, final_l500_plasmid_contig, logger, Config)
            keep_logging('END: Assembly Evaluation using QUAST', 'END: Assembly Evaluation using QUAST', logger, 'info')

        elif args.start_step == 3 and args.end_step == 3:
            keep_logging('START: Assembly Evaluation using QUAST', 'START: Assembly Evaluation using QUAST', logger, 'info')
>>>>>>> 09565177f3a85722989ef80167ec09f603f993ae
            quast_evaluation(args.output_folder, final_l500_contig, final_l500_plasmid_contig, logger, Config)
            keep_logging('END: Assembly Evaluation using QUAST', 'END: Assembly Evaluation using QUAST', logger, 'info')

        elif args.start_step == 3 and args.end_step == 3:
            keep_logging('START: Assembly Evaluation using QUAST', 'START: Assembly Evaluation using QUAST', logger, 'info')
            quast_evaluation(args.output_folder, final_l500_contig, final_l500_plasmid_contig, logger, Config)
            keep_logging('END: Assembly Evaluation using QUAST', 'END: Assembly Evaluation using QUAST', logger, 'info')
        elif args.start_step == 3 and args.end_step == 4:

            contigs = args.output_folder + "/spades_results/contigs.fasta"
            scaffolds = args.output_folder + "/spades_results/contigs.fasta"
            plasmid_contigs = args.output_folder + "/spades_plasmid_results/contigs.fasta"
            plasmid_scaffolds =  args.output_folder + "/spades_plasmid_results/contigs.fasta"

            forward_paired = args.output_folder + ConfigSectionMap("Trimmomatic", Config)['f_p']
            reverse_paired = args.output_folder + ConfigSectionMap("Trimmomatic", Config)['r_p']
            forward_unpaired = args.output_folder + ConfigSectionMap("Trimmomatic", Config)['f_up']
            reverse_unpaired = args.output_folder + ConfigSectionMap("Trimmomatic", Config)['r_up']


            # Use Spades assembly as a reference to map back reads and generate BAM for Pilon
            reference = "%s/spades_results/contigs.fasta" % (args.output_folder)
            reference_plasmid = "%s/spades_plasmid_results/contigs.fasta" % (args.output_folder)

            if args.pilon and args.pilon == "no":
                # Add bioawk here
                (final_l500_contig, final_l500_plasmid_contig) = bioawk(reference, reference_plasmid,
                                                                        args.output_folder,
                                                                        args.analysis_name, logger, Config,
                                                                        args.assembly)
            else:
                # Adding Pilon to the pipeline
                (final_l500_contig, final_l500_plasmid_contig) = run_pilon(forward_paired, reverse_paired,
                                                                           forward_unpaired, reverse_unpaired,
                                                                           reference)


            final_l500_contig = "%s/%s_l500_contigs.fasta" % (args.output_folder, args.analysis_name)
            final_l500_plasmid_contig = "%s/%s_l500_plasmid_contigs.fasta" % (args.output_folder, args.analysis_name)

            # Order contigs with reference to the reference genome provided.
            if args.reference and args.reference != "None":
                final_ordered_contigs = abacas(reference_genome_path, final_l500_contig, args.output_folder,
                                               args.analysis_name, logger, Config)
                sed_cmd = "sed -i 's/>.*/>%s/g' %s" % (args.analysis_name[0:20], final_ordered_contigs)
                keep_logging(sed_cmd, sed_cmd, logger, 'debug')
                os.system(sed_cmd)
                final_annotation_folder = prokka(final_ordered_contigs, args.output_folder, args.analysis_name, logger,
                                                 Config)
                keep_logging('Final Prokka Annotation files for genome are in: {}'.format(final_annotation_folder),
                             'Final Prokka Annotation files for genome are in: {}'.format(final_annotation_folder),
                             logger, 'debug')
            else:
<<<<<<< HEAD
                print ("\nSkipping Abacas...\n")
=======
                print ("\nPlease provide a path to reference genome for Abacas\n")
>>>>>>> 09565177f3a85722989ef80167ec09f603f993ae
                final_ordered_contigs = final_l500_contig
                print ("Replacing fasta header NODE with %s" % args.analysis_name[0:15])
                sed_cmd = "sed -i 's/>NODE/>%s/g' %s" % (args.analysis_name[0:15], final_ordered_contigs)
                sed_cmd_2 = "sed -i 's/_length_.*//g' %s" % (final_ordered_contigs)
                keep_logging(sed_cmd, sed_cmd, logger, 'debug')
                keep_logging(sed_cmd_2, sed_cmd_2, logger, 'debug')
                os.system(sed_cmd)
                os.system(sed_cmd_2)
                final_annotation_folder = prokka(final_ordered_contigs, args.output_folder, args.analysis_name, logger,
                                                 Config)
                keep_logging('Final Prokka Annotation files for genome are in: {}'.format(final_annotation_folder),
                             'Final Prokka Annotation files for genome are in: {}'.format(final_annotation_folder),
                             logger, 'debug')
            sed_cmd = "sed -i 's/>NODE/>%s/g' %s" % (args.analysis_name[0:15], final_l500_plasmid_contig)
            sed_cmd_2 = "sed -i 's/_length_.*component//g' %s" % (final_l500_plasmid_contig)
            keep_logging(sed_cmd, sed_cmd, logger, 'debug')
            keep_logging(sed_cmd_2, sed_cmd, logger, 'debug')
            os.system(sed_cmd)
            os.system(sed_cmd_2)
            plasmid_first_part = args.analysis_name + "_plasmid"
            final_plasmid_annotation_folder = prokka(final_l500_plasmid_contig, args.output_folder, plasmid_first_part,
                                                     logger, Config)
            keep_logging('Final Prokka Annotation files for plasmid are in: {}'.format(final_plasmid_annotation_folder),
                         'Final Prokka Annotation files for plasmid are in: {}'.format(final_plasmid_annotation_folder),
                         logger, 'debug')

            keep_logging('START: Assembly Evaluation using QUAST', 'START: Assembly Evaluation using QUAST', logger,
                         'info')
            quast_evaluation(args.output_folder, final_l500_contig, final_l500_plasmid_contig, logger, Config)
            keep_logging('END: Assembly Evaluation using QUAST', 'END: Assembly Evaluation using QUAST', logger, 'info')


        elif args.start_step == 4 and args.end_step == 4:
            final_l500_contig = "%s/%s_l500_contigs.fasta" % (args.output_folder, args.analysis_name)
            final_l500_plasmid_contig = "%s/%s_l500_plasmid_contigs.fasta" % (args.output_folder, args.analysis_name)
            if args.reference and args.reference != "None":
                final_ordered_contigs = abacas(reference_genome_path, final_l500_contig, args.output_folder, args.analysis_name, logger, Config)

                sed_cmd = "sed -i 's/>.*/>%s/g' %s" % (args.analysis_name[0:20], final_ordered_contigs)
                keep_logging(sed_cmd, sed_cmd, logger, 'debug')
                os.system(sed_cmd)
                final_annotation_folder = prokka(final_ordered_contigs, args.output_folder, args.analysis_name, logger, Config)
                keep_logging('Final Prokka Annotation files for genome are in: {}'.format(final_annotation_folder), 'Final Prokka Annotation files for genome are in: {}'.format(final_annotation_folder), logger, 'debug')
            else:
                print ("\nPlease provide a path to reference genome for Abacas\n")
                final_ordered_contigs = final_l500_contig

                sed_cmd = "sed -i 's/>NODE/>%s/g' %s" % (args.analysis_name[0:15], final_ordered_contigs)
                sed_cmd_2 = "sed -i 's/_length_.*//g' %s" % (final_ordered_contigs)

                keep_logging(sed_cmd, sed_cmd, logger, 'debug')
                keep_logging(sed_cmd_2, sed_cmd_2, logger, 'debug')

                os.system(sed_cmd)
                os.system(sed_cmd_2)
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


""" Start: Check Subroutines """

""" Usage Message: """
def usage():
    print ("Usage: pipeline.py [-h] [-v] [-f1 FILE_1] [-f2 FILE_2] [-o OUTPUT_FOLDER] [-start_step START] [-end_step END] [--qa] [-A Assembler] [--err] \n")

""" Check Java Availability: """
def java_check():
    print ("Checking Java Availability....\n")
    jd = subprocess.check_output(["java", "-version"], stderr=subprocess.STDOUT)
    if len(jd) < 1:
        print ("Unable to find a java runtime environment. The pipeline requires java 6 or later.")
    else:
        print ("Java Availability Check completed ...\n\n" + str(jd))

""" Make sure input raw reads files exists at given location. """
def file_exists(path1, path2):
    if path1:
        if not os.path.isfile(path1):
            file_basename = os.path.basename(path1)
            print ("The input file " + file_basename + " does not exists. \nPlease provide another file.\n")
            exit()
    if path2 != "None":
        if not os.path.isfile(path2):
            file_basename = os.path.basename(path2)
            print ("The input file " + file_basename + " does not exists. \nPlease provide another file.\n")
            exit()

""" Make sure the output folder exists or create at given path """
def make_sure_path_exists(out_path):
    try:
        os.makedirs(out_path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            print ("Errors in output folder path! please change the output path or analysis name\n")
            exit()

""" Prepare the paths for clean reads when all the input files exists. """
def prepare_cleanreadsinput():
    global forward_paired, forward_unpaired, reverse_unpaired, reverse_paired
    forward_paired = out_name + "forward_paired.fq.gz"
    forward_unpaired = out_name + "forward_unpaired.fq.gz"
    reverse_paired = out_name + "reverse_paired.fq.gz"
    reverse_unpaired = out_name + "reverse_unpaired.fq.gz"
    return forward_paired, forward_unpaired, reverse_paired, reverse_unpaired

""" Set appropriate path for contigs/scaffolds fasta files. Needed only when quast runs individually. """
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

""" End: Check Subroutines """

""" Methods """

""" Downsample Raw reads data to Default 100 X or user specified coverage """
def downsample(args, logger):
    if args.coverage_depth:
        keep_logging('Downsampling Coverage Depth to: %s' % args.coverage_depth, 'Downsampling Coverage Depth to: %s' % args.coverage_depth, logger, 'info')

    keep_logging("Running: mash sketch -o /tmp/sketch_out -k 32 -m 3 -r %s >& /tmp/sketch_stdout" % args.file_1, "Running: mash sketch -o /tmp/sketch_out -k 32 -m 3 -r %s >& /tmp/sketch_stdout" % args.file_1, logger, 'info')

    mash_cmd = "mash sketch -o /tmp/sketch_out -k 32 -m 3 -r %s >& /tmp/sketch_stdout" % args.file_1

    keep_logging('Running: %s' % mash_cmd,
                 'Running: %s' % mash_cmd, logger, 'info')

    try:
        call(mash_cmd, logger)
    except sp.CalledProcessError:
        keep_logging('Error running Mash for estimating genome size.', 'Error running Mash for estimating genome size', logger, 'exception')
        sys.exit(1)


    # with open("/tmp/sketch_stdout", 'rU') as file_open:
    #     for line in file_open:
    #         if line.startswith('Estimated genome size:'):
    #             gsize = float(line.split(': ')[1].strip())
    #         if line.startswith('Estimated coverage:'):
    #             est_cov = float(line.split(': ')[1].strip())
    # file_open.close()
    #
    # keep_logging('Estimated Genome Size from Mash Sketch: %s' % gsize,
    #              'Estimated Genome Size from Mash Sketch: %s' % gsize, logger, 'info')


    if args.genome_size:
        gsize = int(args.genome_size)
        keep_logging('Using Genome Size: %s' % gsize,
                     'Using Genome Size: %s' % gsize, logger, 'info')
    else:
        try:
            call(mash_cmd, logger)
        except sp.CalledProcessError:
            keep_logging('Error running Mash for estimating genome size.', 'Error running Mash for estimating genome size', logger, 'exception')
            sys.exit(1)

        with open("/tmp/sketch_stdout", 'rU') as file_open:
            for line in file_open:
                if line.startswith('Estimated genome size:'):
                    gsize = float(line.split(': ')[1].strip())
                if line.startswith('Estimated coverage:'):
                    est_cov = float(line.split(': ')[1].strip())
        file_open.close()


        keep_logging('Estimated Genome Size from Mash Sketch: %s' % gsize,
                     'Estimated Genome Size from Mash Sketch: %s' % gsize, logger, 'info')

    seqtk_check = "seqtk fqchk -q3 %s > /tmp/%s_fastqchk.txt" % (
    args.file_1, os.path.basename(args.file_1))
    keep_logging('Running seqtk to extract Fastq statistics: %s' % seqtk_check,
                 'Running seqtk to extract Fastq statistics: %s' % seqtk_check, logger, 'info')


    try:
        call(seqtk_check, logger)
    except sp.CalledProcessError:
        keep_logging('Error running seqtk for extracting fastq statistics.', 'Error running seqtk for extracting fastq statistics.', logger, 'exception')
        sys.exit(1)

    with open("/tmp/%s_fastqchk.txt" % os.path.basename(args.file_1), 'rU') as file_open:
        for line in file_open:
            if line.startswith('min_len'):
                line_split = line.split(';')
                min_len = line_split[0].split(': ')[1]
                max_len = line_split[1].split(': ')[1]
                avg_len = line_split[2].split(': ')[1]
            if line.startswith('ALL'):
                line_split = line.split('\t')
                total_bases = int(line_split[1]) * 2
    file_open.close()

    keep_logging('Average Read Length: %s' % avg_len,
                 'Average Read Length: %s' % avg_len, logger, 'info')

    keep_logging('Total number of bases in fastq: %s' % total_bases,
                 'Total number of bases in fastq: %s' % total_bases, logger, 'info')

    ori_coverage_depth = int(total_bases / gsize)
    keep_logging('Original Covarage Depth: %s x' % ori_coverage_depth,
                 'Original Covarage Depth: %s x' % ori_coverage_depth, logger, 'info')

    if not args.coverage_depth and ori_coverage_depth > 100:
        # Downsample to 100
        factor = float(100 / float(ori_coverage_depth))
    else:
        factor = float(float(args.coverage_depth) / float(ori_coverage_depth))

    # proc = sp.Popen(["nproc"], stdout=sp.PIPE, shell=True)
    # (nproc, err) = proc.communicate()
    # nproc = nproc.strip()

    proc = (sp.Popen(["nproc"], stdout=sp.PIPE, shell=True, universal_newlines=True)).communicate()[0]
    nproc = proc.strip()

    if ori_coverage_depth > 100:
        # Downsample to 100
        factor = float(100 / float(ori_coverage_depth))
        print(factor)
        r1_sub = "/tmp/%s" % os.path.basename(args.file_1)

        # Downsample using seqtk
        try:
            keep_logging("seqtk sample %s %s | pigz --fast -c -p %s > /tmp/%s" % (
                args.file_1, factor, nproc, os.path.basename(args.file_1)),
                         "/nfs/esnitkin/bin_group/seqtk/seqtk sample %s %s | pigz --fast -c -p %s > /tmp/%s" % (
                             args.file_1, factor, nproc, os.path.basename(args.file_1)), logger, 'info')
            call("seqtk sample %s %s | pigz --fast -c -p %s > /tmp/%s" % (
                args.file_1, factor, nproc, os.path.basename(args.file_1)), logger)
        except sp.CalledProcessError:
            keep_logging('Error running seqtk for downsampling raw fastq reads.',
                         'Error running seqtk for downsampling raw fastq reads.', logger, 'info')
            sys.exit(1)

        if args.file_2:
            r2_sub = "/tmp/%s" % os.path.basename(args.file_2)

            try:
                keep_logging("seqtk sample %s %s | pigz --fast -c -p %s > /tmp/%s" % (
                    args.file_2, factor, nproc, os.path.basename(args.file_2)),
                             "seqtk sample %s %s | pigz --fast -c -p %s > /tmp/%s" % (
                                 args.file_2, factor, nproc, os.path.basename(args.file_2)), logger, 'info')
                call("seqtk sample %s %s | pigz --fast -c -p %s > /tmp/%s" % (
                    args.file_2, factor, nproc, os.path.basename(args.file_2)), logger)
            except sp.CalledProcessError:
                keep_logging('Error running seqtk for downsampling raw fastq reads.',
                             'Error running seqtk for downsampling raw fastq reads.', logger, 'exception')
                sys.exit(1)
        else:
            r2_sub = "None"
    else:
        r1_sub = args.file_1
        if args.file_2:
            r2_sub = args.file_2
        else:
            r2_sub = "None"

    # if not args.coverage_depth and ori_coverage_depth > 100:
    #     # Downsample to 100
    #     factor = float(100 / float(ori_coverage_depth))
    #     print (factor)
    #     r1_sub = "/tmp/%s" % os.path.basename(args.file_1)
    #
    #     # Downsample using seqtk
    #     try:
    #         keep_logging("seqtk sample %s %s | pigz --fast -c -p %s > /tmp/%s" % (
    #             args.file_1, factor, nproc, os.path.basename(args.file_1)),
    #                      "seqtk sample %s %s | pigz --fast -c -p %s > /tmp/%s" % (
    #             args.file_1, factor, nproc, os.path.basename(args.file_1)), logger, 'info')
    #         call("seqtk sample %s %s | pigz --fast -c -p %s > /tmp/%s" % (
    #             args.file_1, factor, nproc, os.path.basename(args.file_1)), logger)
    #     except sp.CalledProcessError:
    #         keep_logging('Error running seqtk for downsampling raw fastq reads.',
    #                      'Error running seqtk for downsampling raw fastq reads.', logger, 'info')
    #         sys.exit(1)
    #
    #     if args.file_2:
    #         r2_sub = "/tmp/%s" % os.path.basename(args.file_2)
    #
    #         try:
    #             keep_logging("seqtk sample %s %s | pigz --fast -c -p %s > /tmp/%s" % (
    #                 args.file_2, factor, nproc, os.path.basename(args.file_2)),
    #                          "seqtk sample %s %s | pigz --fast -c -p %s > /tmp/%s" % (
    #                 args.file_2, factor, nproc, os.path.basename(args.file_2)), logger, 'info')
    #             call("seqtk sample %s %s | pigz --fast -c -p %s > /tmp/%s" % (
    #                 args.file_2, factor, nproc, os.path.basename(args.file_2)), logger)
    #         except sp.CalledProcessError:
    #             keep_logging('Error running seqtk for downsampling raw fastq reads.',
    #                          'Error running seqtk for downsampling raw fastq reads.', logger, 'exception')
    #             sys.exit(1)
    #     else:
    #         r2_sub = "None"
    # elif not args.coverage_depth and ori_coverage_depth < 100:
    #     r1_sub = args.file_1
    #     r2_sub = args.file_2
    # else:
    #     factor = float(float(args.coverage_depth) / float(ori_coverage_depth))
    #     #print round(factor, 3)
    #     r1_sub = "/tmp/%s" % os.path.basename(args.file_1)
    #
    #     # Downsample using seqtk
    #     try:
    #         keep_logging("seqtk sample %s %s | pigz --fast -c -p %s > /tmp/%s" % (
    #             args.file_1, factor, nproc, os.path.basename(args.file_1)),
    #                      "seqtk sample %s %s | pigz --fast -c -p %s > /tmp/%s" % (
    #             args.file_1, factor, nproc, os.path.basename(args.file_1)), logger, 'info')
    #         call("seqtk sample %s %s | pigz --fast -c -p %s > /tmp/%s" % (
    #             args.file_1, factor, nproc, os.path.basename(args.file_1)), logger)
    #     except sp.CalledProcessError:
    #         keep_logging('Error running seqtk for downsampling raw fastq reads.',
    #                      'Error running seqtk for downsampling raw fastq reads.', logger, 'exception')
    #         sys.exit(1)
    #
    #     if args.file_2:
    #         r2_sub = "/tmp/%s" % os.path.basename(args.file_2)
    #
    #         try:
    #             keep_logging("seqtk sample %s %s | pigz --fast -c -p %s > /tmp/%s" % (
    #                 args.file_2, factor, nproc, os.path.basename(args.file_2)),
    #                          "seqtk sample %s %s | pigz --fast -c -p %s > /tmp/%s" % (
    #                 args.file_2, factor, nproc, os.path.basename(args.file_2)), logger, 'info')
    #             call("seqtk sample %s %s | pigz --fast -c -p %s > /tmp/%s" % (
    #                 args.file_2, factor, nproc, os.path.basename(args.file_2)), logger)
    #         except sp.CalledProcessError:
    #             keep_logging('Error running seqtk for downsampling raw fastq reads.',
    #                          'Error running seqtk for downsampling raw fastq reads.', logger, 'exception')
    #             sys.exit(1)
    #     else:
    #         r2_sub = "None"
    #
    return r1_sub, r2_sub

""" Prepare ReadGroup option for BWA alignment """
def prepare_readgroup(forward_read, aligner, logger):
    keep_logging('Preparing ReadGroup Info', 'Preparing ReadGroup Info', logger, 'info')
    samplename = os.path.basename(forward_read)
    if forward_read.endswith(".gz"):
        output = gzip.open(forward_read, 'rb')
        firstLine = output.readline()
        firstLine = firstLine.decode('utf-8')
        if ":" in str(firstLine):
            split_field = re.split(r":",firstLine)
            id_name = split_field[1].rstrip()
            id_name = id_name.rstrip()

        ###Pending
        elif "/" in firstLine:
            split_field = re.split(r"/",firstLine)
            id_name = split_field[1].rstrip()
            id_name = id_name.rstrip()
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

""" Prepare Reference Genome Index for Bowtie, samtools and Picard """
def create_index(reference,ref_index_suffix1, ref_index_suffix2, ref_index_suffix3, ref_index_suffix4, ref_index_suffix5):
    aligner = ConfigSectionMap("pipeline", Config)['aligner']
    keep_logging('Creating Index of reference fasta file for {} aligner.'.format(aligner), 'Creating Index of reference fasta file for {} aligner'.format(aligner), logger, 'info')
    if aligner == "bwa":
        cmd = "%s/%s/%s %s %s" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("bwa", Config)['bwa_bin'], ConfigSectionMap("bwa", Config)['base_cmd'], ConfigSectionMap("bwa", Config)['index'], reference)
        keep_logging(cmd, cmd, logger, 'debug')
        try:
            call(cmd, logger)
        except sp.CalledProcessError:
                keep_logging('Error in {} Indexer. Exiting.'.format(aligner), 'Error in {} Indexer. Exiting.'.format(aligner), logger, 'exception')
                sys.exit(1)
        if not os.path.isfile(ref_index_suffix1):
            keep_logging('The {} reference index files were not created properly. Please try to create the index files again or manually.'.format(aligner), 'The {} reference index files were not created properly. Please try to create the index files again or manually.'.format(aligner), logger, 'exception')
    elif aligner == "bowtie":
        cmd = "%s/%s/%s %s %s" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("bowtie", Config)['bowtie_bin'], ConfigSectionMap("bowtie", Config)['build_cmd'], reference, reference)
        keep_logging(cmd, cmd, logger, 'debug')
        try:
            call(cmd, logger)
        except sp.CalledProcessError:
                keep_logging('Error in {} Indexer. Exiting.'.format(aligner), 'Error in {} Indexer. Exiting.'.format(aligner), logger, 'exception')
                sys.exit(1)
        if not os.path.isfile(ref_index_suffix1):
            keep_logging('The {} reference index files were not created properly. Please try to create the index files again or manually.'.format(aligner), 'The {} reference index files were not created properly. Please try to create the index files again or manually.'.format(aligner), logger, 'exception')

    else:
        print ("Different Aligner in config file")

def create_fai_index(reference, ref_fai_index):
    keep_logging('Creating FAI Index using Samtools.', 'Creating FAI Index using Samtools.', logger, 'info')
    cmd = "%s %s %s" % (ConfigSectionMap("samtools", Config)['base_cmd'], ConfigSectionMap("samtools", Config)['faiindex'], reference)
    keep_logging(cmd, cmd, logger, 'debug')
    try:
        call(cmd, logger)
    except sp.CalledProcessError:
        keep_logging('Error in Samtools FAI Indexing step. Exiting.', 'Error in Samtools FAI Indexing step. Exiting.', logger, 'exception')
        sys.exit(1)


    if not os.path.isfile(ref_fai_index):
        keep_logging('The reference fai index file {} was not created properly.\n Please try to create the samtools fai index files manually. \n'.format(ref_fai_index), 'The reference fai index file {} was not created properly.\n Please try to create the samtools fai index files manually. \n'.format(ref_fai_index), logger, 'exception')
    else:
        keep_logging('Samtools Fai Index file created.', 'Samtools Fai Index file created.', logger, 'info')

def picard_seqdict(dict_name, reference):
    keep_logging('Creating Sequence Dictionary using Picard.', 'Creating Sequence Dictionary using Picard.', logger, 'info')
    cmd = "%s CreateSequenceDictionary REFERENCE=%s OUTPUT=%s/%s" % (ConfigSectionMap("picard", Config)['base_cmd'], reference, os.path.dirname(reference), dict_name)
    keep_logging(cmd, cmd, logger, 'debug')
    try:
        call(cmd, logger)
    except sp.CalledProcessError:
        keep_logging('Error in Picard Sequence Dictionary creation step. Exiting.', 'Error in Picard Sequence Dictionary creation step. Exiting.', logger, 'exception')
        sys.exit(1)

""" samtools: Post-Alignment SAM/BAM conversion, sort, index """
def prepare_bam(out_sam, out_path, analysis, files_to_delete, logger, Config):
    out_bam = samtobam(out_sam, out_path, analysis, files_to_delete, logger, Config)
    out_sort_bam = sort_bam(out_bam, out_path, analysis, logger, Config)
    #files_to_delete.append(out_sort_bam)
    out_marked_bam = markduplicates(out_sort_bam, out_path, analysis, files_to_delete, logger, Config)
    out_sort_bam = sort_bam(out_marked_bam, out_path, analysis, logger, Config)
    index_bam(out_sort_bam, out_path, logger, Config)
    if not os.path.isfile(out_sort_bam):
        keep_logging('Error in SAM/BAM conversion, sort, index. Exiting.', 'Error in SAM/BAM conversion, sort, index. Exiting.', logger, 'exception')
        exit()
    else:
        return out_sort_bam

def post_align(out_sam, reference):
    keep_logging('START: Post-Alignment using SAMTOOLS, PICARD etc...', 'START: Post-Alignment using SAMTOOLS, PICARD etc...', logger, 'info')
    out_sorted_bam = prepare_bam(out_sam, args.output_folder, args.analysis_name, files_to_delete, logger, Config)
    keep_logging('END: Post-Alignment using SAMTOOLS, PICARD etc...', 'END: Post-Alignment using SAMTOOLS, PICARD etc...', logger, 'info')
    #out_sorted_bam = "%s/%s_aln_sort.bam" % (args.output_folder, args.analysis_name)
    # keep_logging('START: Creating BedGraph Coverage', 'START: Creating BedGraph Coverage', logger, 'info')
    # bedgraph_coverage(out_sorted_bam, args.output_folder, args.analysis_name, reference, logger, Config)
    # only_unmapped_positions_file = bedtools(out_sorted_bam, args.output_folder, args.analysis_name, logger, Config)
    # keep_logging('END: Creating BedGraph Coverage', 'END: Creating BedGraph Coverage', logger, 'info')
    return out_sorted_bam

def file_exists_for_alignment(path1, path2, reference):
    if not os.path.isfile(reference):
        file_basename = os.path.basename(reference)
        keep_logging('The reference fasta file {} does not exists. Please provide another with full path file with full path or check the files path.\n'.format(file_basename), 'The reference fasta file {} does not exists. Please provide another file or check the files path.\n'.format(file_basename), logger, 'exception')
        exit()
    if ConfigSectionMap("pipeline", Config)['aligner'] == "bwa":
        ref_index_suffix1 = reference + ".bwt"
        ref_index_suffix2 = reference + ".amb"
        ref_index_suffix3 = reference + ".ann"
        ref_index_suffix4 = reference + ".sa"
        ref_index_suffix5 = reference + ".pac"
    elif ConfigSectionMap("pipeline", Config)['aligner'] == "bowtie":
        ref_index_suffix1 = reference + ".1.bt2"
        ref_index_suffix2 = reference + ".2.bt2"
        ref_index_suffix3 = reference + ".3.bt2"
        ref_index_suffix4 = reference + ".4.ebwt"
        ref_index_suffix5 = reference + ".rev.1.bt2"
        ref_index_suffix6 = reference + ".rev.2.bt2"
    if not os.path.isfile(ref_index_suffix1):
        keep_logging('The reference index files given below does not exists:\n {}\n {}\n {}\n {}\n {}'.format(ref_index_suffix1, ref_index_suffix2, ref_index_suffix3, ref_index_suffix4, ref_index_suffix5), 'The reference index files given below does not exists:\n {}\n {}\n {}\n {}\n {}'.format(ref_index_suffix1, ref_index_suffix2, ref_index_suffix3, ref_index_suffix4, ref_index_suffix5), logger, 'warning')
        create_index(reference, ref_index_suffix1, ref_index_suffix2, ref_index_suffix3, ref_index_suffix4, ref_index_suffix5)
    else:
        keep_logging('Index file already exists.', 'Index file already exists.', logger, 'info')

    ref_fai_index = reference + ".fai"
    if not os.path.isfile(ref_fai_index):
        keep_logging('The reference fai index file {} required for samtools does not exists.'.format(ref_fai_index), 'The reference fai index file {} required for samtools does not exists.'.format(ref_fai_index), logger, 'warning')
        create_fai_index(reference, ref_fai_index)
    else:
        keep_logging('Samtools fai Index file already exists.', 'Samtools fai Index file already exists.', logger, 'info')

    dict_name = os.path.splitext(os.path.basename(reference))[0] + ".dict"
    if not os.path.isfile("%s" % reference + "/" + dict_name):
        keep_logging('The reference seq dict file {} required for GATK and PICARD does not exists.'.format(dict_name), 'The reference seq dict file {} required for GATK and PICARD does not exists.'.format(dict_name), logger, 'warning')
        picard_seqdict(dict_name, reference)
    else:
        keep_logging('The reference seq dict file required for GATK and PICARD exists.', 'The reference seq dict file required for GATK and PICARD exists.', logger, 'info')

def run_pilon(forward_paired, reverse_paired, forward_unpaired, reverse_unpaired, reference):
    # Post assembly improvements
    # Changes Added on August 26 2019
    split_field = prepare_readgroup(forward_paired, "bwa", logger)
    files_to_delete = []

    # Check if FASTQ files exists; Generate Reference Genome index for alignments
    if args.type != "PE" and args.type != "BAM":
        reverse_raw = "None"
        file_exists_for_alignment(args.file_1, args.file_1, reference)
    elif args.type != "PE" and args.type != "SE":
        print ("BAM type... Not Integrated... continue")
    else:
        file_exists_for_alignment(args.file_1, args.file_2, reference)

    # Align reads to reference genome

    # Prepare Aligner base command
    base_cmd = ConfigSectionMap("bin_path", Config)['binbase'] + "/" + ConfigSectionMap("bowtie", Config)[
        'bowtie_bin'] + "/" + ConfigSectionMap("bowtie", Config)['align_cmd']
    files_to_delete = ""
    parameters = ConfigSectionMap("bowtie", Config)['parameters']

    # Align reads using Bowtie
    out_sam = align_bowtie(base_cmd, forward_paired, reverse_paired, forward_unpaired, reverse_unpaired,
                           args.output_folder, reference,
                           split_field, args.analysis_name, files_to_delete, logger, Config, args.type, parameters)

    # Post-process BAM files: Sort, index and run mark duplicates
    out_sorted_bam = post_align(out_sam, reference)

    out_sorted_bam = "%s/%s_aln_sort.bam" % (args.output_folder, args.analysis_name)
    print (out_sorted_bam)

    polished_pilon_fasta = pilon(out_sorted_bam, reference, args.output_folder, args.analysis_name, files_to_delete, logger, Config)

    plasmid_contigs = args.output_folder + "spades_plasmid_results" + "/contigs.fasta"

    (final_l500_contig, final_l500_plasmid_contig) = bioawk(polished_pilon_fasta, plasmid_contigs, args.output_folder,
                                                            args.analysis_name, logger, Config, args.assembly)

    return final_l500_contig, final_l500_plasmid_contig


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
    global files_to_delete
    files_to_delete = ""
    log_unique_time = datetime.now().strftime('%Y_%m_%d_%H_%M_%S')
    logger = generate_logger(args.output_folder, args.analysis_name, log_unique_time)
    global Config
    Config = configparser.ConfigParser()
    Config.read(config_file)
    pipeline(args, logger)
    keep_logging('End: Pipeline', 'End: Pipeline', logger, 'info')
    time_taken = datetime.now() - start_time_2
    keep_logging('Total Time taken: {}'.format(time_taken), 'Total Time taken: {}'.format(time_taken), logger, 'info')





