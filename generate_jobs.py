__author__ = 'alipirani'

import os
import argparse
import re

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
parser.add_argument('-reference_key', action='store', dest="reference_key", help='Reference Genome to be used for each sample. \
Provide a tab-seperated file with fastq filename in first column and reference genome name in corresponding second column')
parser.add_argument('-depthofcoverage', action='store', dest="depthofcoverage", help='Varcall pipeline until GATK Depth of Coverage')
parser.add_argument('-steps', action='store', dest="steps", help='Variant Calling Steps in sequential order.\n'
                                                                     '1.   All : This will run all the steps starting from cleaning the reads to variant calling;\n'
                                                                     '2.   clean,align,post-align,varcall,filter,stats : This will also run all steps starting from cleaning to variant calling. \nYou can also run part of the pipeline by giving "align,post-align,varcall,filter,stats" which will skip the cleaning part.\nThe order is required to be sequential. Also, while skipping any of the step make sure you have results already present in your output folder.\n'
                                                                     '3.   coverage_depth_stats: Run Only Depth of Coverage Stats module after cleaning and read mapping steps')

parser.add_argument('-email', action='store', dest="email", help='email for PBS job notifications')
parser.add_argument('-resources', action='store', dest="resources", help='PBS resources for jobs')
parser.add_argument('-assembly', action='store', dest="assembly", help='assembly: wga/plasmid/both. Default is both ')

args = parser.parse_args()

# read filenames file
filenames_array = []
with open(args.filenames) as fp:
    for line in fp:
        line = line.strip()
        line = args.dir + "/" + line
        filenames_array.append(line)



# Function to create (old)assembly pipeline pbs scripts for each fastq files
def create_assembly_jobs():
    ##Change the email address; nodes and processor requirements accordingly
    Pbs_model_lines = "#PBS -M %s\n" \
        "#PBS -m a\n" \
        "#PBS -V\n" \
        "#PBS -l %s\n" \
        "#PBS -q fluxod\n" \
        "#PBS -A esnitkin_fluxod\n" \
        "#PBS -l qos=flux\n" % (args.email, args.resources)

    for file in filenames_array:
        # Forward reads file name and get analysis name from its name
        filename_base = os.path.basename(file)
        first_part_split = filename_base.split('.')
        first_part = first_part_split[0]
        first_file = file

        second_part = filename_base.replace("_R1_", "_R2_")
        #second_part = filename_base.replace("forward", "reverse")
        #second_part = filename_base.replace("1_sequence", "2_sequence")
        second_file = args.dir + "/" + second_part

        # Change these directory path to where your pipeline code is located
        cd_command = "cd /nfs/esnitkin/bin_group/pipeline/assembly_umich/"
        command = "/nfs/esnitkin/bin_group/anaconda2/bin//python pipeline.py -f1 %s -f2 %s -o %s/%s -start_step 1 -end_step 4 -A spades -reference %s" % (first_file, second_file, args.out_dir, first_part, args.reference)
        job_name = "./" + first_part + ".pbs"

        with open(job_name, 'w') as out:
            job_title = "#PBS -N %s" % first_part
            out.write(job_title+'\n')
            out.write(Pbs_model_lines+'\n')
            out.write(cd_command+'\n')
            out.write(command+'\n')

# Function to create assembly pbs scripts for each fastq files
def create_new_assembly_jobs():
    ##Change the email address; nodes and processor requirements accordingly
    Pbs_model_lines = "#PBS -M %s\n" \
        "#PBS -m a\n" \
        "#PBS -V\n" \
        "#PBS -l %s\n" \
        "#PBS -q fluxod\n" \
        "#PBS -A esnitkin_fluxod\n" \
        "#PBS -l qos=flux\n" % (args.email, args.resources)

    for file in filenames_array:
        first_file = file
        # Forward file name and get analysis name from its name
        filename_base = os.path.basename(file)
        if "R1_001_final.fastq.gz" in filename_base:
            second_part = filename_base.replace("R1_001_final.fastq.gz", "R2_001_final.fastq.gz")
            first_part_split = filename_base.split('R1_001_final.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
        elif "R1.fastq.gz" in filename_base:
            second_part = filename_base.replace("R1.fastq.gz", "R2.fastq.gz")
            first_part_split = filename_base.split('R1.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
        elif "1_combine.fastq.gz" in filename_base:
            second_part = filename_base.replace("1_combine.fastq.gz", "2_combine.fastq.gz")
            first_part_split = filename_base.split('1_combine.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
        elif "1_sequence.fastq.gz" in filename_base:
            second_part = filename_base.replace("1_sequence.fastq.gz", "2_sequence.fastq.gz")
            first_part_split = filename_base.split('1_sequence.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
        elif "_forward.fastq.gz" in filename_base:
            second_part = filename_base.replace("_forward.fastq.gz", "_reverse.fastq.gz")
            first_part_split = filename_base.split('_forward.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
        elif "R1_001.fastq.gz" in filename_base:
            second_part = filename_base.replace("R1_001.fastq.gz", "R2_001.fastq.gz")
            first_part_split = filename_base.split('R1_001.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
            #first_part = first_part.replace('_S*_', '')
            first_part = re.sub(r'_S.*_', '', first_part)
        elif "_1.fastq.gz" in filename_base:
            second_part = filename_base.replace("_1.fastq.gz", "_2.fastq.gz")
            first_part_split = filename_base.split('_1.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
        elif "_R1.fastq.gz" in filename_base:
            second_part = filename_base.replace("_R1.fastq.gz", "_R2.fastq.gz")
            first_part_split = filename_base.split('_R1.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
        else:
            print "Using Standard second file naming convention"
            second_part = filename_base.replace("_R1_", "_R2_")
            first_part_split = filename_base.split('_R1.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')

        # Get the name of reverse reads files
        #second_part = filename_base.replace("_R1_", "_R2_")
        second_file = args.dir + "/" + second_part


        # Change these directory paths to where your pipeline code is located. Also make sure the path to config file is correct.
        cd_command = "cd /nfs/esnitkin/bin_group/pipeline/Github/assembly_umich/"
        config = "/nfs/esnitkin/bin_group/pipeline/Github/assembly_umich/config"
        job_name = "./" + first_part + ".pbs"
        if args.assembly:
            assembly_para = "-assembly " + args.assembly
        else:
            assembly_para = "-assembly both"
        if args.reference:
            command = "module load perl-modules\n/nfs/esnitkin/bin_group/anaconda2/bin//python pipeline.py -f1 %s -f2 %s -o %s/%s -start_step 1 -end_step 4 -A spades -reference %s -type PE -analysis %s -config %s %s" % (first_file, second_file, args.out_dir, first_part, args.reference, first_part, config, assembly_para)
        else:
            command = "module load perl-modules\n/nfs/esnitkin/bin_group/anaconda2/bin//python pipeline.py -f1 %s -f2 %s -o %s/%s -start_step 1 -end_step 4 -A spades -type PE -analysis %s -config %s %s" % (first_file, second_file, args.out_dir, first_part, first_part, config, assembly_para)
        with open(job_name, 'w') as out:
                job_title = "#PBS -N %s" % first_part
                out.write(job_title+'\n')
                out.write(Pbs_model_lines+'\n')
                out.write(cd_command+'\n')
                out.write(command+'\n')

# Function to create assembly pbs scripts for each fastq files
def create_new_assembly_SE_jobs():
    ##Change the email address; nodes and processor requirements accordingly
    Pbs_model_lines = "#PBS -M %s\n" \
        "#PBS -m a\n" \
        "#PBS -V\n" \
        "#PBS -l %s\n" \
        "#PBS -q fluxod\n" \
        "#PBS -A esnitkin_fluxod\n" \
        "#PBS -l qos=flux\n" % (args.email, args.resources)

    for file in filenames_array:
        first_file = file
        # Forward file name and get analysis name from its name
        filename_base = os.path.basename(file)
        if "R1_001_final.fastq.gz" in filename_base:
            second_part = filename_base.replace("R1_001_final.fastq.gz", "R2_001_final.fastq.gz")
            first_part_split = filename_base.split('R1_001_final.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
        elif "R1.fastq.gz" in filename_base:
            second_part = filename_base.replace("R1.fastq.gz", "R2.fastq.gz")
            first_part_split = filename_base.split('R1.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
        elif "1_combine.fastq.gz" in filename_base:
            second_part = filename_base.replace("1_combine.fastq.gz", "2_combine.fastq.gz")
            first_part_split = filename_base.split('1_combine.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
        elif "1_sequence.fastq.gz" in filename_base:
            second_part = filename_base.replace("1_sequence.fastq.gz", "2_sequence.fastq.gz")
            first_part_split = filename_base.split('1_sequence.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
        elif "_forward.fastq.gz" in filename_base:
            second_part = filename_base.replace("_forward.fastq.gz", "_reverse.fastq.gz")
            first_part_split = filename_base.split('_forward.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
        elif "R1_001.fastq.gz" in filename_base:
            second_part = filename_base.replace("R1_001.fastq.gz", "R2_001.fastq.gz")
            first_part_split = filename_base.split('R1_001.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
            #first_part = first_part.replace('_S*_', '')
            first_part = re.sub(r'_S.*_', '', first_part)
        elif "_1.fastq.gz" in filename_base:
            second_part = filename_base.replace("_1.fastq.gz", "_2.fastq.gz")
            first_part_split = filename_base.split('_1.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
        elif "_R1.fastq.gz" in filename_base:
            second_part = filename_base.replace("_R1.fastq.gz", "_R2.fastq.gz")
            first_part_split = filename_base.split('_R1.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
        elif "_unclassified.fastq.gz" in filename_base:
            second_part = filename_base.replace("_unclassified.fastq.gz", "_unclassified.fastq.gz")
            first_part_split = filename_base.split('_unclassified.fastq.gz')
            first_part = first_part_split[0].replace('_L00', '')
        else:
            print "Using Standard second file naming convention"
            second_part = filename_base.replace("_R1_", "_R2_")
            first_part_split = filename_base.split('_R1.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')

        # Get the name of reverse reads files
        #second_part = filename_base.replace("_R1_", "_R2_")
        second_file = args.dir + "/" + second_part


        # Change these directory paths to where your pipeline code is located. Also make sure the path to config file is correct.
        cd_command = "cd /nfs/esnitkin/bin_group/pipeline/Github/assembly_umich/"
        config = "/nfs/esnitkin/bin_group/pipeline/Github/assembly_umich/config"
        job_name = "./" + first_part + ".pbs"

        if args.reference:
            command = "module load perl-modules\n/nfs/esnitkin/bin_group/anaconda2/bin//python pipeline.py -f1 %s -o %s/%s -start_step 1 -end_step 4 -A spades -reference %s -type SE -analysis %s -config %s" % (first_file, args.out_dir, first_part, args.reference, first_part, config)
        else:
            command = "module load perl-modules\n/nfs/esnitkin/bin_group/anaconda2/bin//python pipeline.py -f1 %s -o %s/%s -start_step 1 -end_step 4 -A spades -type SE -analysis %s -config %s" % (first_file, args.out_dir, first_part, first_part, config)
        with open(job_name, 'w') as out:
                job_title = "#PBS -N %s" % first_part
                out.write(job_title+'\n')
                out.write(Pbs_model_lines+'\n')
                out.write(cd_command+'\n')
                out.write(command+'\n')



# Function to create assembly pbs scripts for each fastq files
def create_ariba_jobs():
    ##Change the email address; nodes and processor requirements accordingly
    Pbs_model_lines = "#PBS -M %s\n" \
        "#PBS -m a\n" \
        "#PBS -V\n" \
        "#PBS -l nodes=1:ppn=4,pmem=4000mb,walltime=60:00:00\n" \
        "#PBS -q fluxod\n" \
        "#PBS -A esnitkin_fluxod\n" \
        "#PBS -l qos=flux\n" % (args.email, args.resources)

    for file in filenames_array:
        first_file = file
        # Forward file name and get analysis name from its name
        filename_base = os.path.basename(file)
        if "R1_001_final.fastq.gz" in filename_base:
            second_part = filename_base.replace("R1_001_final.fastq.gz", "R2_001_final.fastq.gz")
            first_part_split = filename_base.split('R1_001_final.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
        elif "R1.fastq.gz" in filename_base:
            second_part = filename_base.replace("R1.fastq.gz", "R2.fastq.gz")
            first_part_split = filename_base.split('R1.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
        elif "1_combine.fastq.gz" in filename_base:
            second_part = filename_base.replace("1_combine.fastq.gz", "2_combine.fastq.gz")
            first_part_split = filename_base.split('1_combine.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
        elif "1_sequence.fastq.gz" in filename_base:
            second_part = filename_base.replace("1_sequence.fastq.gz", "2_sequence.fastq.gz")
            first_part_split = filename_base.split('1_sequence.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
        elif "_forward.fastq.gz" in filename_base:
            second_part = filename_base.replace("_forward.fastq.gz", "_reverse.fastq.gz")
            first_part_split = filename_base.split('_forward.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
        elif "R1_001.fastq.gz" in filename_base:
            second_part = filename_base.replace("R1_001.fastq.gz", "R2_001.fastq.gz")
            first_part_split = filename_base.split('R1_001.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
            #first_part = first_part.replace('_S*_', '')
            first_part = re.sub(r'_S.*_', '', first_part)
        elif "_1.fastq.gz" in filename_base:
            second_part = filename_base.replace("_1.fastq.gz", "_2.fastq.gz")
            first_part_split = filename_base.split('_1.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
        elif "_R1.fastq.gz" in filename_base:
            second_part = filename_base.replace("_R1.fastq.gz", "_R2.fastq.gz")
            first_part_split = filename_base.split('_R1.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')
        else:
            print "Using Standard second file naming convention"
            second_part = filename_base.replace("_R1_", "_R2_")
            first_part_split = filename_base.split('_R1.fastq.gz')
            first_part = first_part_split[0].replace('_L001', '')

        # Get the name of reverse reads files
        #second_part = filename_base.replace("_R1_", "_R2_")
        second_file = args.dir + "/" + second_part


        # Change these directory paths to where your pipeline code is located. Also make sure the path to config file is correct.
        #mkdir_cmd = "mkdir %s/%s" % (args.out_dir, first_part)
        cd_command = "cd %s" % (args.out_dir)
        job_name = "./" + first_part + ".pbs"

        command = "module load cd-hit\n/nfs/esnitkin/bin_group/anaconda3/bin/ariba run %s %s %s %s" % (args.ariba_db, first_file, second_file, first_part)
        with open(job_name, 'w') as out:
                job_title = "#PBS -N %s" % first_part
                out.write(job_title+'\n')
                out.write(Pbs_model_lines+'\n')
                #out.write(mkdir_cmd+'\n')
                out.write(cd_command+'\n')
                out.write(command+'\n')




# Main Function; Check which Pipeline function to use based on command-line argument: pipeline
if args.pipeline == "varcall":
    print "Generating Variant Calling PBS scripts for samples in args.filenames"
    create_varcall_jobs()
elif args.pipeline == "assembly":
    print "Generating Assembly PBS scripts for samples in args.filenames"
    create_assembly_jobs()
elif args.pipeline == "new_assembly":
    print "Generating new Assembly PBS scripts for samples in args.filenames"
    create_new_assembly_jobs()
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
else:
    print "Please provide args.pipeline command-line argument"
