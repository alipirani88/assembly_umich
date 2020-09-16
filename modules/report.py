__author__ = 'alipirani'
"""
Generate Assembly and Ariba reports.
- This script generates different types of reports from assembly, prokka, and ariba(AMR/MLSt) results.
- Run this script after all the assembly jobs are finished.
- It will aggregate all the results that it finds in out_dir and generate a Report folder in out_dir/Report.

Input: out_dir where all the assembly results from assembly pipeline are located.

Output:
    
    1. 
    2. 
    3. 
"""

""" Declaring required python modules """
import os
import argparse
import re
import errno
import glob
from datetime import datetime
import readline
from joblib import Parallel, delayed
import multiprocessing

"""Parse command line arguments"""
parser = argparse.ArgumentParser(description='Aggregate results from assembly pipeline and generate Assembly and Ariba reports')
parser.add_argument('-out_dir', action='store', dest="out_dir", help='Provide a path where you want to save the assembly output')
args = parser.parse_args()

""" Function: Make sure the output folder exists or create at given path """
def make_sure_path_exists(out_path):
    try:
        os.makedirs(out_path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            print "Errors in output folder path! please change the output path or analysis name\n"
            exit()


def run_command(i):
    print "Running: %s" % i
    os.system(i)
    done = "Completed: %s" % i
    return done

""" Generate MultiQC assembly report """
def assembly_report(args):
    make_sure_path_exists("%s/Report/Quast" % args.out_dir)
    list_of_quast_dir = glob.glob("%s/*/quast_results" % args.out_dir)
    move_cmd_array = []
    for i in list_of_quast_dir:
        sample_name = os.path.basename(os.path.dirname(i))
        move_cmd = "cp -r %s %s/Report/Quast/%s" % (i, args.out_dir, sample_name)
        move_cmd_array.append(move_cmd)
        #os.system(move_cmd)
    results = Parallel(n_jobs=num_cores)(delayed(run_command)(i) for i in move_cmd_array)
    uniq_time_prefix = datetime.now().strftime('%Y_%m_%d_%H_%M_%S')
    print "Generating MultiQC assembly report for samples in %s" % args.out_dir
    os.chdir("%s/Report/Quast/" % args.out_dir)
    multiqc_cmd = "/nfs/esnitkin/bin_group/anaconda2/bin/multiqc -f --outdir %s/Report/Quast/ -n %s_assembly_report -i %s_assembly_report %s/Report/Quast/" % (args.out_dir, uniq_time_prefix, uniq_time_prefix, args.out_dir)
    os.system(multiqc_cmd)

    make_sure_path_exists("%s/Report/Quast_plasmid" % args.out_dir)
    list_of_quast_plasmid_dir = glob.glob("%s/*/quast_plasmid_results" % args.out_dir)
    move_cmd_array = []
    for i in list_of_quast_plasmid_dir:
        sample_name = os.path.basename(os.path.dirname(i))
        move_cmd = "cp -r %s %s/Report/Quast_plasmid/%s" % (i, args.out_dir, sample_name)
        #os.system(move_cmd)
        move_cmd_array.append(move_cmd)
    results = Parallel(n_jobs=num_cores)(delayed(run_command)(i) for i in move_cmd_array)
    uniq_time_prefix = datetime.now().strftime('%Y_%m_%d_%H_%M_%S')
    print "Generating MultiQC plasmid assembly report for samples in %s" % args.out_dir
    multiqc_cmd = "/nfs/esnitkin/bin_group/anaconda2/bin/multiqc -f --outdir %s/Report/Quast_plasmid/ -n %s_plasmid_assembly_report -i %s_plasmid_assembly_report %s/Report/Quast_plasmid/" % (args.out_dir, uniq_time_prefix, uniq_time_prefix, args.out_dir)
    os.system(multiqc_cmd)

    start_time = datetime.now().strftime('%Y-%m-%d')
    make_sure_path_exists("%s/Results/%s_assembly" % (args.out_dir, start_time))
    make_sure_path_exists("%s/Results/%s_plasmid_assembly" % (args.out_dir, start_time))
    list_of_assembly_fasta = glob.glob("%s/*/*_l500_contigs.fasta" % args.out_dir)
    for i in list_of_assembly_fasta:
        move_cmd = "cp %s %s/Results/%s_assembly/%s.fasta" % (i, args.out_dir, start_time, os.path.basename(i).replace('_l500_contigs.fasta', ''))
        print move_cmd
        os.system(move_cmd)
    list_of_plasmid_assembly_fasta = glob.glob("%s/*/*_plasmid_contigs.fasta" % args.out_dir)
    for i in list_of_plasmid_assembly_fasta:
        move_cmd = "cp %s %s/Results/%s_plasmid_assembly/%s.fasta" % (i, args.out_dir, start_time, os.path.basename(i).replace('_plasmid_contigs.fasta', ''))
        print move_cmd
        os.system(move_cmd)



def ariba_report(args):
    uniq_time_prefix = datetime.now().strftime('%Y_%m_%d_%H_%M_%S')
    print "Generating Ariba report for samples in %s" % args.out_dir
    os.chdir(args.out_dir)
    make_sure_path_exists("%s/Report/" % args.out_dir)
    os.system("ls */*_MLST/mlst_report.tsv > %s/Report/mlst_filenames" % (args.out_dir))
    os.system(
        "for i in */*_MLST; do cat $i/mlst_report.tsv | grep '^ST'; done | sort | uniq > %s/Report/header.txt" % (args.out_dir))
    os.system("for i in */*_MLST; do cat $i/mlst_report.tsv | grep -v '^ST'; done > %s/Report/mlst_temp.txt" % (args.out_dir))
    os.system("paste %s/Report/mlst_filenames %s/Report/mlst_temp.txt > %s/Report/mlst_report.txt" % (args.out_dir, args.out_dir, args.out_dir))
    os.system("cat %s/Report/header.txt %s/Report/mlst_report.txt > %s/Report/MLST_report.txt" % (
    args.out_dir, args.out_dir, args.out_dir))
    # ariba_cmd = "/nfs/esnitkin/bin_group/anaconda3/bin/ariba summary --preset minimal %s/Report/%s_minimal_AMR_report %s/*/*_AMR/report.tsv" % (args.out_dir, uniq_time_prefix, args.out_dir)
    # ariba_cmd = "/nfs/esnitkin/bin_group/anaconda3/bin/ariba summary --preset all %s/Report/%s_all_AMR_report %s/*/*_AMR/report.tsv" % (
    # args.out_dir, uniq_time_prefix, args.out_dir)
    ariba_cmd = "/nfs/esnitkin/bin_group/anaconda3/bin/python /nfs/esnitkin/bin_group/ariba/scripts/ariba summary --preset minimal %s/Report/%s_minimal_AMR_report %s/*/*_AMR/report.tsv" % (
    args.out_dir, uniq_time_prefix, args.out_dir)
    ariba_cmd = "/nfs/esnitkin/bin_group/anaconda3/bin/python /nfs/esnitkin/bin_group/ariba/scripts/ariba summary --preset all %s/Report/%s_all_AMR_report %s/*/*_AMR/report.tsv" % (
        args.out_dir, uniq_time_prefix, args.out_dir)
    print ariba_cmd
    os.system(ariba_cmd)



"""Main Function"""
if __name__ == '__main__':
    """ Generate script start time variables for logging"""
    start_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    start_time_2 = datetime.now()

    """ Generate variables for logging modules"""
    #global logger
    #log_unique_time = datetime.now().strftime('%Y_%m_%d_%H_%M_%S')
    #logger = generate_logger(args.output_folder, args.analysis_name, log_unique_time)
    num_cores = multiprocessing.cpu_count()

    """ Generate Report folder to save the results"""
    make_sure_path_exists("%s/Report" % args.out_dir)

    """ Run assembly_report module to aggregate the assembly results and generate various assembly and annotation reports"""
    assembly_report(args)

    """ Run ariba_report module to aggregate the Ariba results and generate Ariba MLST and AMR reports"""
    #ariba_report(args)

    time_taken = datetime.now() - start_time_2
    #keep_logging('Total Time taken: {}'.format(time_taken), 'Total Time taken: {}'.format(time_taken), logger, 'info')
