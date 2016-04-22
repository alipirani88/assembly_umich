__author__ = 'alipirani'


####################################################################################################################################################################################
# Usage: pipeline.py [-h] [-f1 FILE_1] [-f2 FILE_2] [-o OUTPUT_FOLDER] [-start_step START] [-end_step END] [--qa] [-A Assembler] [--err]
#
# optional arguments:
#
#    -h, --help		show this help message and exit
#
#    -v, --version         	show program's version number and exit
#
#    -f1 FILE_1            	Paired End fastq file 1
#
#    -f2 FILE_2            	Paired End fastq file 2
#
#    -O OUTPUT_FOLDER        Output Path to save the results
#
#    -start_step START_STEP  Provide the start step
#
#    -end_step END_STEP    	Provide the end step
#
#    -A ASSEMBLER          	Choose the assembler to assemble the sample reads.
#                         	Velvet Optimiser or Spades
#
#    --qa                  	Run FastQC for quality check
#
#    --err                   Run Spades using built-in BayesHammer Error Corrector
#
#    -reference                 Provide reference genome in case you select step 4 that involves Abacas contig ordering and Annotation
#
#    Optional:
#
#         	FastQC: User can run FastQC on raw sample data using the --qa flag.
###################################################################################################################################################################################

# Declaring required python modules
import argparse
from config_settings import ConfigSectionMap
from check_subroutines import *
from trimmomatic import *
from assembly import *
from qa_fastqc import qa_fastqc
from quast import quast_evaluation
from reapr import reapr
from abacas import abacas
from bioawk import bioawk
from prokka import prokka

# Print initial pipeline messages
print "\nGenome Assembly Pipeline for Illumina PE datasets.\n"

# Command line argument parsing
parser = argparse.ArgumentParser(description='Assembly pipeline for Illumina PE data')
parser.add_argument('-f1', action='store', dest="file_1", help='Paired End file 1')
parser.add_argument('-f2', action='store', dest="file_2", help='Paired End file 2', )
parser.add_argument('-o', action='store', dest="output_folder", help='Output Path ending with output directory name to save the results', required=True)
parser.add_argument('-start_step', action="store", dest="start_step", type=int, help='Provide the start step')
parser.add_argument('-end_step', action="store", dest="end_step", type=int, help='Provide the end step')
parser.add_argument('-A', action='store', dest="assembler", help='Choose the assembler to assemble the sample reads. Velvet Optimiser or Spades')
# parser.add_argument('-PE', action='store_true', default=False, dest='boolean_switch', help='Indicates Sample files are Paired End')
parser.add_argument('--qa', action='store_true', default=False, help='Run FastQC for quality check')
parser.add_argument('--err', action='store_true', default=False, help='Run BayesHammer error correction built in Spades')
parser.add_argument('-c', action='store', dest="crop", help='choose crop value to crop the reads')
parser.add_argument('-reference', action='store', dest="reference", help='Please provide a reference genome for Abacas Contig ordering')
args = parser.parse_args()
out_path = args.output_folder + "/"
filename_base = os.path.basename(args.file_1)
#first_part = filename_base[0:15]
first_part = filename_base.replace(".fastq.gz", "")
#crop_length = str(args.crop)
# END of Argument Parsing step

####################################################################SUBROUTINES###################################################################################################
# check if the output folder exists; create the output folder for saving results and call appropriate functions
make_sure_path_exists(out_path)
# check if the input file exists
file_exists(args.file_1, args.file_2)
# check java availability
java_check()
# create an extra config section namely [check_clean] at the start of pipeline to bring it in default format.
create_config_section()
reference_genome_path = args.reference

# start the pipeline with the given start and end steps:
if args.start_step and args.end_step:
    Forward_read = args.file_1
    Reverse_read = args.file_2
    # Temporarily save the path for clean reads to be used in subroutines where Trimmomatic and Assembly are not called.
    # In this case, the pipeline assumes the required clean reads(forward PE, Reverse PE, Forward unpaired, Reverse Unpaired) or the contigs/scaffolds fasta file are present in output directory
    forward_paired = out_path + ConfigSectionMap("Trimmomatic")['f_p']
    reverse_paired = out_path + ConfigSectionMap("Trimmomatic")['r_p']
    forward_unpaired = out_path + ConfigSectionMap("Trimmomatic")['f_up']
    reverse_unpaired = out_path + ConfigSectionMap("Trimmomatic")['r_up']
    #contigs = out_path + "/contigs.fa"
    #scaffolds = out_path + "/scaffolds.fa"
    if args.start_step == 1 and args.end_step == 1:
        print "\n################## Pre-processing step using Trimmomatic. ##################\n"
        (forward_paired, reverse_paired, forward_unpaired, reverse_unpaired) = clean_reads(Forward_read, Reverse_read, out_path, args.crop)
    elif args.start_step == 1 and args.end_step == 2:
        print "\n################## Pre-processing and Assembly step. ##################\n"
        (forward_paired, reverse_paired, forward_unpaired, reverse_unpaired) = clean_reads(Forward_read, Reverse_read, out_path, args.crop)
        (contigs, scaffolds) = assembly(forward_paired, reverse_paired, forward_unpaired, reverse_unpaired, args.assembler, out_path)
    elif args.start_step == 1 and args.end_step == 3:
        print "\n################## Pre-processing, Assembly and Evaluation. ##################\n"
        (forward_paired, reverse_paired, forward_unpaired, reverse_unpaired) = clean_reads(Forward_read, Reverse_read, out_path, args.crop)
        assembly(forward_paired, reverse_paired, forward_unpaired, reverse_unpaired, args.assembler, out_path)
        (contigs, scaffolds) = get_contigs(out_path, args.assembler)
        quast_evaluation(out_path, contigs, scaffolds)
    elif args.start_step == 1 and args.end_step == 4:
        print "\n################## Pre-processing, Assembly, Evaluation and Misassembly detection. ##################\n"
        if args.reference:
            reference_genome_path = args.reference
            (forward_paired, reverse_paired, forward_unpaired, reverse_unpaired) = clean_reads(Forward_read, Reverse_read, out_path, args.crop)
            assembly(forward_paired, reverse_paired, forward_unpaired, reverse_unpaired, args.assembler, out_path)
            (contigs, scaffolds) = get_contigs(out_path, args.assembler)
            quast_evaluation(out_path, contigs, scaffolds)

            # reapr(forward_paired, reverse_paired, out_path, contigs, scaffolds)

            final_l500_contig = bioawk(contigs, out_path, first_part)
            final_ordered_contigs = abacas(reference_genome_path, final_l500_contig, out_path, first_part)
	    
            #final_ordered_contigs = "/nfs/esnitkin/Ali/Project_MRSA_analysis/MRSA_assembly///6154_R1.fastq.g_contigs_ordered.fasta"
	    sed_cmd = "sed -i 's/>.*/>%s/g' %s" % (first_part, final_ordered_contigs)
	    os.system(sed_cmd)
	    final_annotation_folder = prokka(final_ordered_contigs, out_path, first_part)
            print "Final Prokka Annotation files are in: %s" % final_annotation_folder

        else:
            print "\nPlease provide a path to reference genome for Abacas\n"

    elif args.start_step == 2 and args.end_step == 3:
        print "\n################## Assembly and Evaluation step. ##################\n "
        (contigs, scaffolds) = assembly(forward_paired, reverse_paired, forward_unpaired, reverse_unpaired, args.assembler, out_path)
        quast_evaluation(out_path, contigs, scaffolds)
    elif args.start_step == 2 and args.end_step == 2:
        print "\n################## Assembly step. ##################\n"
        (contigs, scaffolds) = assembly(forward_paired, reverse_paired, forward_unpaired, reverse_unpaired, args.assembler, out_path)
    elif args.start_step == 2 and args.end_step == 4:
        if args.reference:
            print "\n################## Assembly, Evaluation and Misassembly detection step. ##################\n "
            (contigs, scaffolds) = assembly(forward_paired, reverse_paired, forward_unpaired, reverse_unpaired, args.assembler, out_path)
            quast_evaluation(out_path, contigs, scaffolds)
            #reapr(forward_paired, reverse_paired, out_path, contigs, scaffolds)
            final_l500_contig = bioawk(contigs, out_path, first_part)
            final_ordered_contigs = abacas(reference_genome_path, final_l500_contig, out_path, first_part)
	    sed_cmd = "sed -i 's/>.*/>%s/g' %s" % (first_part, final_ordered_contigs)
            os.system(sed_cmd)
            final_annotation_folder = prokka(final_ordered_contigs, out_path, first_part)
            print "Final Prokka Annotation files are in: %s" % final_annotation_folder
        else:
            print "\nPlease provide a path to reference genome for Abacas\n"
    elif args.start_step == 3 and args.end_step == 3:
        print "\n################## Evaluation step. ##################\n"
        (contigs, scaffolds) = get_contigs(out_path, args.assembler)
        quast_evaluation(out_path, contigs, scaffolds)
    elif args.start_step == 3 and args.end_step == 4:
        if args.reference:
            print "\n################## Evaluation and Misassembly detection step. ##################\n"
            (contigs, scaffolds) = get_contigs(out_path, args.assembler)
            #filename = os.path.basename(out_path)
            final_l500_contig = bioawk(contigs, out_path, first_part)
            quast_evaluation(out_path, contigs, scaffolds)
            #reapr(forward_paired, reverse_paired, out_path, contigs, scaffolds)
            final_l500_contig = bioawk(contigs, out_path, first_part)
            final_ordered_contigs = abacas(reference_genome_path, final_l500_contig, out_path, first_part)
            sed_cmd = "sed -i 's/>.*/>%s/g' %s" % (first_part, final_ordered_contigs)
            os.system(sed_cmd)
	    final_annotation_folder = prokka(final_ordered_contigs, out_path, first_part)
            print "Final Prokka Annotation files are in: %s" % final_annotation_folder
        else:
            print "\nPlease provide a path to reference genome for Abacas\n"

    elif args.start_step == 4 and args.end_step == 4:
        if args.reference:
            #print "\n################## Misassembly Detection. ##################\n"
            (contigs, scaffolds) = get_contigs(out_path, args.assembler)
            #reapr(forward_paired, reverse_paired, out_path, contigs, scaffolds)
            final_l500_contig = bioawk(contigs, out_path, first_part)
            final_ordered_contigs = abacas(reference_genome_path, final_l500_contig, out_path, first_part)
            print "This is final ordered contigs fasta file: %s" % final_ordered_contigs
            sed_cmd = "sed -i 's/>.*/>%s/g' %s" % (first_part, final_ordered_contigs)
            os.system(sed_cmd)
	    final_annotation_folder = prokka(final_ordered_contigs, out_path, first_part)
            print "Final Prokka Annotation files are in: %s" % final_annotation_folder
        else:
            print "\nPlease provide a path to reference genome for Abacas\n"

else:
    print "ERROR: Please provide start and end steps for the pipeline.\n"


#Quality Assessment using FastQC
if args.qa:
    qa_fastqc(out_path, args.file_1, args.file_2)
###################################################################################################################################################################################


