__author__ = 'alipirani'

import os
import argparse
import re
import errno
import glob
from datetime import datetime

parser = argparse.ArgumentParser(description='Generate Assembly and Ariba report')
parser.add_argument('-out_dir', action='store', dest="out_dir", help='Provide a path where you want to save the assembly output')
parser.add_argument('-ariba_preset', action='store', dest="preset", help='Provide ariba summary preset option for generating the report.')
args = parser.parse_args()

""" Make sure the output folder exists or create at given path """
def make_sure_path_exists(out_path):
    try:
        os.makedirs(out_path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            print "Errors in output folder path! please change the output path or analysis name\n"
            exit()


# Generate MultiQC assembly report
def assembly_report(args):
    quast_dir = "%s/"
    make_sure_path_exists("%s/Report/Quast" % args.out_dir)
    list_of_quast_dir = glob.glob("%s/*/quast_results" % args.out_dir)
    for i in list_of_quast_dir:
        sample_name = os.path.basename(os.path.dirname(i))
        move_cmd = "cp -r %s %s/Report/Quast/%s" % (i, args.out_dir, sample_name)
        os.system(move_cmd)
    uniq_time_prefix = datetime.now().strftime('%Y_%m_%d_%H_%M_%S')
    print "Generating MultiQC assembly report for samples in %s" % args.out_dir
    multiqc_cmd = "/nfs/esnitkin/bin_group/anaconda2/bin/multiqc -f --outdir %s/Report/Quast/ -n %s_assembly_report -i %s_assembly_report %s/Report/Quast/" % (args.out_dir, uniq_time_prefix, uniq_time_prefix, args.out_dir)
    os.system(multiqc_cmd)
    print "Generating Ariba report for samples in %s" % args.out_dir
    ariba_cmd = "/nfs/esnitkin/bin_group/anaconda3/bin/ariba summary --preset %s %s/Report/%s_AMR_report %s/*/*_AMR/report.tsv" % (args.preset, args.out_dir, uniq_time_prefix, args.out_dir)
    os.system(ariba_cmd)
"""
1. Create MultiQC assembly report with quast results
2. Create Ariba AMR and MLST report
"""
if __name__ == '__main__':
    start_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    start_time_2 = datetime.now()
    #args = parser().parse_args()
    #global logger
    #log_unique_time = datetime.now().strftime('%Y_%m_%d_%H_%M_%S')
    #logger = generate_logger(args.output_folder, args.analysis_name, log_unique_time)
    make_sure_path_exists("%s/Report" % args.out_dir)
    assembly_report(args)
    #keep_logging('End: Pipeline', 'End: Pipeline', logger, 'info')
    time_taken = datetime.now() - start_time_2
    #keep_logging('Total Time taken: {}'.format(time_taken), 'Total Time taken: {}'.format(time_taken), logger, 'info')
