__author__ = 'alipirani'

import os
import argparse
import re
import glob
parser = argparse.ArgumentParser(description='Generate PBS scripts for Different Pipelines for multiple samples: varcall/assembly/new_assembly/rna/ptr')
parser.add_argument('-dir', action='store', dest="dir", help='directory of fastq files')
parser.add_argument('-filenames', action='store', dest="filenames", help='These file should contain name of forward fastq files. \
One file per line. \
These can be obtained by running: ls *_R1_001.fastq.gz > filenames \
Make Sure your Forward reads extension ends with \"_R1_001.fastq.gz\" ')
parser.add_argument('-out_dir', action='store', dest="out_dir", help='Provide a path where you want to save the assembly output')
parser.add_argument('-pipeline', action='store', dest="pipeline", help='Generating Jobs for which pipeline? varcall/assembly/new_assembly/new_assembly_SE/rna/ptr/ariba')
parser.add_argument('-reference', action='store', dest="reference", help='Reference Genome to be used for pipeline')
parser.add_argument('-type', action='store', dest="type", help='Type of Fastq files: PE or SE')
parser.add_argument('-ariba_db', action='store', dest="ariba_db", help='Path to Ariba Database')
parser.add_argument('-depthofcoverage', action='store', dest="depthofcoverage", help='Varcall pipeline until GATK Depth of Coverage')
parser.add_argument('-steps', action='store', dest="steps", help='Variant Calling Steps in sequential order.\n'
                                                                     '1.   All : This will run all the steps starting from cleaning the reads to variant calling;\n'
                                                                     '2.   clean,align,post-align,varcall,filter,stats : This will also run all steps starting from cleaning to variant calling. \nYou can also run part of the pipeline by giving "align,post-align,varcall,filter,stats" which will skip the cleaning part.\nThe order is required to be sequential. Also, while skipping any of the step make sure you have results already present in your output folder.\n'
                                                                     '3.   coverage_depth_stats: Run Only Depth of Coverage Stats module after cleaning and read mapping steps')
parser.add_argument('-ariba', action='store', dest="ariba", help='Run Ariba to find AMR / MLST / both: possible options are AMR, MLST, both')
parser.add_argument('-email', action='store', dest="email", help='email for PBS job notifications')
parser.add_argument('-resources', action='store', dest="resources", help='PBS resources for jobs')
parser.add_argument('-assembly', action='store', dest="assembly", help='assembly: wga/plasmid/both. Default is both ')
parser.add_argument('-suffix', action='store', dest="suffix", help='suffix for fastq read files ')
parser.add_argument('-config', action='store', dest="config", help='Path to Config file, Make sure to check config settings before running pipeline', required=False)
parser.add_argument('-downsample', action='store', dest="downsample",
                          help='yes/no: Downsample Reads data to default depth of 100X or user specified depth')
parser.add_argument('-coverage_depth', action='store', dest="coverage_depth",
                          help='Downsample Reads to this user specified depth')
args = parser.parse_args()

def get_filenames(dir, type, filenames, suffix):
    if not filenames:
        if not suffix:
            suffix = ".fastq.gz"
        try:
            list_of_files = glob.glob("%s/*%s" % (dir, suffix))
            if len(list_of_files) < 1:
                print "No fastq files with suffix %s found in reads directory %s" % (suffix, dir)
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

# Function to create assembly pbs scripts for each fastq files
def create_new_assembly_jobs(list_of_files):
    ##Change the email address; nodes and processor requirements accordingly
    Pbs_model_lines = "#PBS -M %s\n" \
        "#PBS -m a\n" \
        "#PBS -V\n" \
        "#PBS -l %s\n" \
        "#PBS -q flux\n" \
        "#PBS -A esnitkin_flux\n" \
        "#PBS -l qos=flux\n" % (args.email, args.resources)
    os.system("cd %s" % args.out_dir)
    job_array = []
    for file in list_of_files:
        filename_base = os.path.basename(file)
	#print filename_base
        #if "R1_001_final.fastq.gz" in filename_base or "R1.fastq.gz" in filename_base or "1_combine.fastq.gz" in filename_base or "1_sequence.fastq.gz" in filename_base or "_forward.fastq.gz" in filename_base or "R1_001.fastq.gz" in filename_base or "_1.fastq.gz" in filename_base or ".1.fastq.gz" in filename_base or "_R1.fastq.gz" in filename_base or "_R2.fastq.gz" not in filename_base:
        if re.search("R1_001_final.fastq.gz", filename_base) or re.search("R1.fastq.gz", filename_base) or re.search("1_combine.fastq.gz", filename_base) or re.search("1_sequence.fastq.gz", filename_base) or re.search("_forward.fastq.gz", filename_base) or re.search("R1_001.fastq.gz", filename_base) or re.search("_1.fastq.gz", filename_base) or re.search(".1.fastq.gz", filename_base) or re.search("_R1.fastq.gz", filename_base) or re.search("_R1_001.fastq.gz", filename_base):
            # Forward reads file name and get analysis name from its name
            first_file = file
            #print first_file
            # Get the name of reverse reads files
            if "R1_001_final.fastq.gz" in filename_base:
                second_part = filename_base.replace("R1_001_final.fastq.gz", "R2_001_final.fastq.gz")
                first_part_split = filename_base.split('R1_001_final.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S.*_", "", first_part)
            elif "_R1.fastq.gz" in filename_base:
                second_part = filename_base.replace("_R1.fastq.gz", "_R2.fastq.gz")
                first_part_split = filename_base.split('_R1.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S.*_", "", first_part)
            elif "R1.fastq.gz" in filename_base:
                second_part = filename_base.replace("R1.fastq.gz", "R2.fastq.gz")
                first_part_split = filename_base.split('R1.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S.*_", "", first_part)
            elif "1_combine.fastq.gz" in filename_base:
                second_part = filename_base.replace("1_combine.fastq.gz", "2_combine.fastq.gz")
                first_part_split = filename_base.split('1_combine.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S.*_", "", first_part)
            elif "1_sequence.fastq.gz" in filename_base:
                second_part = filename_base.replace("1_sequence.fastq.gz", "2_sequence.fastq.gz")
                first_part_split = filename_base.split('1_sequence.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S.*_", "", first_part)
            elif "_forward.fastq.gz" in filename_base:
                second_part = filename_base.replace("_forward.fastq.gz", "_reverse.fastq.gz")
                first_part_split = filename_base.split('_forward.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S.*_", "", first_part)
            elif "_R1_001.fastq.gz" in filename_base:
		print "here"
                second_part = filename_base.replace("_R1_001.fastq.gz", "_R2_001.fastq.gz")
                first_part_split = filename_base.split('_R1_001.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S.*_", "", first_part)
            elif "_1.fastq.gz" in filename_base:
                second_part = filename_base.replace("_1.fastq.gz", "_2.fastq.gz")
                first_part_split = filename_base.split('_1.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S.*_", "", first_part)
            elif ".1.fastq.gz" in filename_base:
                second_part = filename_base.replace(".1.fastq.gz", ".2.fastq.gz")
                first_part_split = filename_base.split('.1.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S.*_", "", first_part)
            elif "_1.fastq.gz" in filename_base:
                second_part = filename_base.replace("_1.fastq.gz", "_2.fastq.gz")
                first_part_split = filename_base.split('_1.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S.*_", "", first_part)
            # Get the name of reverse reads files
            #second_part = filename_base.replace("_R1_", "_R2_")
            #second_file = args.dir + "/" + second_part
            second_part = filename_base.replace("_R1.fastq.gz", "_R2.fastq.gz")
            first_part_split = filename_base.split('_R1.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
            first_part = re.sub("_S.*_", "", first_part)
            first_file = file
            second_file = args.dir + "/" + second_part

        # Change these directory paths to where your pipeline code is located. Also make sure the path to config file is correct.
        cd_command = "cd /nfs/esnitkin/bin_group/pipeline/Github/assembly_umich_dev/"
        if args.config:
            config = args.config
        else:
            config = "/nfs/esnitkin/bin_group/pipeline/Github/assembly_umich_dev/config"
        job_name = "%s" % args.out_dir + first_part + ".pbs"
        if args.assembly:
            assembly_para = "-assembly " + args.assembly
        else:
            assembly_para = "-assembly both"
        if args.reference:
            command = "module load perl-modules\n/nfs/esnitkin/bin_group/anaconda2/bin//python pipeline.py -f1 %s -f2 %s -o %s/%s -start_step 1 -end_step 4 -A spades -reference %s -type PE -analysis %s -config %s %s" % (first_file, second_file, args.out_dir, first_part, args.reference, first_part, config, assembly_para)
        else:
            command = "module load perl-modules\n/nfs/esnitkin/bin_group/anaconda2/bin//python pipeline.py -f1 %s -f2 %s -o %s/%s -start_step 1 -end_step 4 -A spades -type PE -analysis %s -config %s %s" % (first_file, second_file, args.out_dir, first_part, first_part, config, assembly_para)
        if args.ariba:
            ariba_para = "-ariba " + args.ariba
            command = "module load perl-modules\n/nfs/esnitkin/bin_group/anaconda2/bin//python pipeline.py -f1 %s -f2 %s -o %s/%s -start_step 1 -end_step 4 -A spades -reference %s -type PE -analysis %s -config %s %s %s" % (first_file, second_file, args.out_dir, first_part, args.reference, first_part, config, assembly_para, ariba_para)

        # # Adding Downsampling support 2019-06-20
        if args.downsample == "yes":
            if args.coverage_depth:
                depth = args.coverage_depth
            else:
                depth = 100
            command = command + " -downsample yes -coverage_depth %s" % depth
        with open(job_name, 'w') as out:
            job_title = "#PBS -N %s" % first_part
            out.write(job_title+'\n')
            out.write(Pbs_model_lines+'\n')
            out.write(cd_command+'\n')
            out.write(command+'\n')
        if job_name not in job_array:
            job_array.append(job_name)
    for jobs in job_array:
        print "qsub %s" % jobs
        #os.system("qsub %s" % jobs)
# Main Function; Check which Pipeline function to use based on command-line argument: pipeline
if args.pipeline == "varcall":
    print "Generating Variant Calling PBS scripts for samples in args.filenames"
    create_varcall_jobs()
elif args.pipeline == "old_assembly":
    print "Generating Assembly PBS scripts for samples in args.filenames"
    create_assembly_jobs()
elif args.pipeline == "assembly":
    print "Generating individual Assembly jobs for samples in %s" % args.dir
    list_of_files = get_filenames(args.dir, args.type, args.filenames, args.suffix)
    create_new_assembly_jobs(list_of_files)
elif args.pipeline == "a5_assembly":
    print "Generating a5 Assembly PBS scripts for samples in args.filenames"
    create_a5_jobs()
elif args.pipeline == "rna":
    print "Generating RNA Seq PBS scripts for samples in args.filenames"
    create_RNA_seq_jobs()
elif args.pipeline == "ptr":
    print "Generating PTR pipeline PBS scripts for samples in args.filenames"
    create_PTR_jobs()
elif args.pipeline == "new_assembly_SE":
    print "Generating new Assembly PBS scripts for samples in args.filenames"
    create_new_assembly_SE_jobs()
elif args.pipeline == "ariba":
    print "Generating new Ariba PBS scripts for samples in args.filenames"
    create_ariba_jobs()
elif args.pipeline == "srst2_plasmid":
    print "Generating new SRST2 Plasmid PBS scripts for samples in args.filenames"
    create_srst2_jobs()
elif args.pipeline == "recombination":
    print "Generating recombination PBS scripts for samples in args.filenames"
    create_recombination_jobs()
else:
    print "Please provide args.pipeline command-line argument"
