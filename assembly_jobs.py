__author__ = 'alipirani'

import os
import argparse
import re
import glob
import errno
#import ConfigParser
import configparser
from config_settings import ConfigSectionMap
from modules.log_modules import  *
from modules.logging_subprocess import *


parser = argparse.ArgumentParser(description='Microbial Genome Assembly Pipeline.')
required = parser.add_argument_group('Required arguments')
optional = parser.add_argument_group('Optional arguments')
required.add_argument('-dir', action='store', dest="dir", help='directory of fastq files')
optional.add_argument('-filenames', action='store', dest="filenames", help='These file should contain name of forward fastq files. \
One file per line. \
These can be obtained by running: ls *_R1_001.fastq.gz > filenames \
Make Sure your Forward reads extension ends with \"_R1_001.fastq.gz\" ')
required.add_argument('-out_dir', action='store', dest="out_dir", help='Provide a path where you want to save the assembly output')
optional.add_argument('-reference', action='store', dest="reference", help='Reference Genome to be used for pipeline')
required.add_argument('-type', action='store', dest="type", help='Type of Fastq files: PE or SE')
optional.add_argument('-ariba', action='store', dest="ariba", help='Run Ariba to find AMR / MLST / both: possible options are AMR, MLST, both')
optional.add_argument('-ariba_db', action='store', dest="ariba_db", help='Path to Ariba Database')
optional.add_argument('-email', action='store', dest="email", help='email for PBS job notifications')
optional.add_argument('-resources', action='store', dest="resources", help='PBS/SLURM resources for jobs')
optional.add_argument('-assembly', action='store', dest="assembly", help='assembly: wga/plasmid/both. Default is both ')
optional.add_argument('-suffix', action='store', dest="suffix", help='suffix for fastq read files ')
optional.add_argument('-config', action='store', dest="config", help='Path to Config file, Make sure to check config settings before running pipeline', required=False)
optional.add_argument('-downsample', action='store', dest="downsample",
                          help='yes/no: Downsample Reads data to default depth of 100X or user specified depth')
optional.add_argument('-coverage_depth', action='store', dest="coverage_depth",
                          help='Downsample Reads to this user specified depth')
optional.add_argument('-scheduler', action='store', dest="scheduler",
                          help='Type of Scheduler for generating cluster jobs: PBS, SLURM, LOCAL')
args = parser.parse_args()


def get_filenames(dir, type, filenames, suffix):
    """Takes fastq file directory, optional list of filenames, type of sequence data and  generate assembly jobs.
            Args:
                filenames: optional list of fastq file names
                dir: fastq files directory
                type: type of reads. PE/SE
                suffix: fastq file name suffix if different then fastq.gz

            Output:
                List of fastq file names from input directory
    """
    if not filenames:
        if not suffix:
            suffix = "_R1_*.fastq.gz"
        try:
            list_of_files = glob.glob("%s/*%s" % (dir, suffix))
            if len(list_of_files) < 1:
                print("No fastq files with suffix %s found in reads directory %s" % (suffix, dir))
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                keep_logging('Error while listing files in reads directory.', 'Error while listing files in reads directory.', logger, 'exception')
                exit()
    else:
        list_of_files = []
        with open(filenames) as fp:
            for line in fp:
                line = line.strip()
                line = dir + "/" + line
                list_of_files.append(line)
    return list_of_files


def make_sure_path_exists(out_path):
    """Checks the directory path exists. If not, creates a new directory.
        Args:
            path: Path to Directory

        Output:
            True, if the directory exists or if not, a new directory is created.
    """
    try:
        os.makedirs(out_path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            keep_logging('Errors in output folder path! please change the output path or analysis name.', 'Errors in output folder path! please change the output path or analysis name', logger, 'exception')
            exit()


def get_scheduler_directive(scheduler, Config):
    """Generate Cluster Directive lines for a scheduler provided with args.scheduler
            Args:
                path: scheduler name, Config object

            Output:
                variables associated with scheduler
    """
    if scheduler and scheduler == "SLURM":
        script_Directive = "#SBATCH"
        job_name_flag = "--job-name="
        scheduler_directives = "#SBATCH --mail-user=%s\n#SBATCH --mail-type=%s\n#SBATCH --export=ALL\n#SBATCH --partition=%s\n#SBATCH --account=%s\n#SBATCH %s\n" \
                          % (ConfigSectionMap("slurm", Config)['email'],
                             ConfigSectionMap("slurm", Config)['notification'],
                             ConfigSectionMap("slurm", Config)['partition'],
                             ConfigSectionMap("slurm", Config)['flux_account'],
                             ConfigSectionMap("slurm", Config)['resources'])
    elif scheduler and scheduler == "PBS":
        script_Directive = "#PBS"
        job_name_flag = "-N"
        scheduler_directives = "#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l %s\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\n" \
                          % (ConfigSectionMap("scheduler", Config)['email'],
                             ConfigSectionMap("scheduler", Config)['notification'],
                             ConfigSectionMap("scheduler", Config)['resources'],
                             ConfigSectionMap("scheduler", Config)['queue'],
                             ConfigSectionMap("scheduler", Config)['flux_account'])
    else:
        script_Directive = "#SBATCH"
        job_name_flag = "--job-name="
        scheduler_directives = "#SBATCH --mail-user=%s\n#SBATCH --mail-type=%s\n#SBATCH --export=ALL\n#SBATCH --partition=%s\n#SBATCH --account=%s\n#SBATCH %s\n" \
                               % (ConfigSectionMap("slurm", Config)['email'],
                                  ConfigSectionMap("slurm", Config)['notification'],
                                  ConfigSectionMap("slurm", Config)['partition'],
                                  ConfigSectionMap("slurm", Config)['flux_account'],
                                  ConfigSectionMap("slurm", Config)['resources'])
    return scheduler_directives, script_Directive, job_name_flag



def create_new_assembly_jobs(list_of_files):
    """Takes fastq file directory, optional list of filenames, type of sequence data and  generate assembly jobs.
            Args:
                filenames: Generate assembly jobs for a list of fastq files

            Output:
    """
    jobs_temp_dir = "%s/temp_jobs" % args.out_dir
    make_sure_path_exists(jobs_temp_dir)
    keep_logging('Generating assembly cluster jobs in temporary directory %s' % jobs_temp_dir,
                 'Generating assembly cluster jobs in temporary directory %s' % jobs_temp_dir, logger, 'exception')

    scheduler_directives, script_Directive, job_name_flag = get_scheduler_directive(args.scheduler, Config)

    job_array = []
    for file in list_of_files:
        filename_base = os.path.basename(file)
        if re.search("R1_001_final.fastq.gz", filename_base) or re.search("R1.fastq.gz", filename_base) or re.search("1_combine.fastq.gz", filename_base) or re.search("1_sequence.fastq.gz", filename_base) or re.search("_forward.fastq.gz", filename_base) or re.search("R1_001.fastq.gz", filename_base) or re.search("_1.fastq.gz", filename_base) or re.search(".1.fastq.gz", filename_base) or re.search("_R1.fastq.gz", filename_base):

            # Forward reads file name and get analysis name from its name
            first_file = file
            # Get the name of reverse reads files
            if "R1_001_final.fastq.gz" in filename_base:
                second_part = filename_base.replace("R1_001_final.fastq.gz", "R2_001_final.fastq.gz")
                first_part_split = filename_base.split('R1_001_final.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S.*_", "", first_part)
            if "_R1.fastq.gz" in filename_base:
                second_part = filename_base.replace("_R1.fastq.gz", "_R2.fastq.gz")
                first_part_split = filename_base.split('_R1.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S.*_", "", first_part)
            if "R1.fastq.gz" in filename_base:
                second_part = filename_base.replace("R1.fastq.gz", "R2.fastq.gz")
                first_part_split = filename_base.split('R1.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S.*_", "", first_part)
            if "1_combine.fastq.gz" in filename_base:
                second_part = filename_base.replace("1_combine.fastq.gz", "2_combine.fastq.gz")
                first_part_split = filename_base.split('1_combine.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S.*_", "", first_part)
            if "1_sequence.fastq.gz" in filename_base:
                second_part = filename_base.replace("1_sequence.fastq.gz", "2_sequence.fastq.gz")
                first_part_split = filename_base.split('1_sequence.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S.*_", "", first_part)
            if "_forward.fastq.gz" in filename_base:
                second_part = filename_base.replace("_forward.fastq.gz", "_reverse.fastq.gz")
                first_part_split = filename_base.split('_forward.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S.*_", "", first_part)
            if "R1_001.fastq.gz" in filename_base:
                second_part = filename_base.replace("R1_001.fastq.gz", "R2_001.fastq.gz")
                first_part_split = filename_base.split('R1_001.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S.*_", "", first_part)
            if re.search("R1_001.fastq.gz", filename_base):
                second_part = filename_base.replace("_R1_001.fastq.gz", "_R2_001.fastq.gz")
                first_part_split = filename_base.split('_R1_001.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S.*_", "", first_part)

            if "_1.fastq.gz" in filename_base:
                second_part = filename_base.replace("_1.fastq.gz", "_2.fastq.gz")
                first_part_split = filename_base.split('_1.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S.*_", "", first_part)
            if ".1.fastq.gz" in filename_base:
                second_part = filename_base.replace(".1.fastq.gz", ".2.fastq.gz")
                first_part_split = filename_base.split('.1.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S.*_", "", first_part)
            if "_1.fastq.gz" in filename_base:
                second_part = filename_base.replace("_1.fastq.gz", "_2.fastq.gz")
                first_part_split = filename_base.split('_1.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S.*_", "", first_part)
            else:
                second_part = filename_base.replace("_R1_001.fastq.gz", "_R2_001.fastq.gz")
                first_part_split = filename_base.split('_R1_001.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S.*_", "", first_part)
            # Get the name of reverse reads files
            # second_part = filename_base.replace("_R1.fastq.gz", "_R2.fastq.gz")
            # first_part_split = filename_base.split('_R1.fastq.gz')
            # first_part = first_part_split[0].replace('_L001', '')
            # first_part = re.sub("_S.*_", "", first_part)
            first_file = file
            second_file = args.dir + "/" + second_part
        else:
            print ("Suffix not found in sample name.")
            exit()
        # Prepare command line optional arguments
        if args.assembly:
            assembly_para = "-assembly " + args.assembly
        else:
            assembly_para = "-assembly both"

        if args.reference:
            command = "python %s/pipeline.py -f1 %s -f2 %s -o %s/%s -start_step 1 -end_step 4 -A spades -reference %s -type PE -analysis %s -config %s %s" % (os.path.dirname(os.path.abspath(__file__)), first_file, second_file, args.out_dir, first_part, args.reference, first_part, config_file, assembly_para)
        else:
            command = "python %s/pipeline.py -f1 %s -f2 %s -o %s/%s -start_step 1 -end_step 4 -A spades -type PE -analysis %s -config %s %s" % (os.path.dirname(os.path.abspath(__file__)), first_file, second_file, args.out_dir, first_part, first_part, config_file, assembly_para)

        if args.ariba:
            ariba_para = "-ariba " + args.ariba
            command = "python %s/pipeline.py -f1 %s -f2 %s -o %s/%s -start_step 1 -end_step 4 -A spades -reference %s -type PE -analysis %s -config %s %s %s" % (os.path.dirname(os.path.abspath(__file__)), first_file, second_file, args.out_dir, first_part, args.reference, first_part, config_file, assembly_para, ariba_para)

        if args.downsample == "yes":
            if args.coverage_depth:
                depth = args.coverage_depth
            else:
                depth = 100
            command = command + " -downsample yes -coverage_depth %s" % depth

        if args.scheduler == "SLURM":
            job_name = jobs_temp_dir + "/" + first_part + ".sbat"
            submit_cmd = "sbatch"
        else:
            job_name = jobs_temp_dir + "/" + first_part + ".pbs"
            submit_cmd = "qsub"

        # Generate cluster jobs
        with open(job_name, 'w') as out:
            job_title = "%s %s%s" % (script_Directive, job_name_flag, first_part)
            out.write("#!/bin/sh" + '\n')
            out.write(job_title + '\n')
            out.write(scheduler_directives + '\n')
            out.write("cd %s/temp_jobs" % args.out_dir + '\n')
            out.write(command + '\n')

        if job_name not in job_array:
            job_array.append(job_name)
            print("{} {}".format(submit_cmd, job_name))
            #os.system("{} {}".format(submit_cmd, job_name))



print("Generating individual Assembly jobs for samples in %s" % args.dir)
list_of_files = get_filenames(args.dir, args.type, args.filenames, args.suffix)
global config_file
global Config
global logger
global log_unique_time
if args.config:
    config_file = args.config
else:
    config_file = os.path.dirname(os.path.abspath(__file__)) + "/config"
print ("%s" % config_file)
log_unique_time = datetime.now().strftime('%Y_%m_%d_%H_%M_%S')
Config = configparser.ConfigParser()
Config.read(config_file)
logger = generate_logger(args.dir, "assembly", log_unique_time)
print (list_of_files)
create_new_assembly_jobs(list_of_files)
