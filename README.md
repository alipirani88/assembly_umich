# Microbial Genome Assembly Pipeline

## Synopsis

This pipeline runs various assembly and post-assembly steps on Illumina PE/SE reads. The different steps of the pipeline are - downsample reads using Seqtk/Mash (optional), clean reads with Trimmomatic, assemble clean reads with SPAdes assembler, perform post-assembly correction with PILON, perform AMR and MLST typing using ARIBA (Deprecated).

## Contents:

- [Installation](#installation)
- [Input](#input)
- [Steps](#steps)
- [Quick Start](#quick-start)
- [Output Files](#output-files)
- [Customizing Config file](#customizing-config-file)
- [Log](#log)

## Installation

The pipeline can be set up in two easy steps:

> 1. Clone the github directory onto your system.

```
git clone https://github.com/alipirani88/assemblage.git

```

> 2. Use assemblage/assemblage.yml and assemblage/assemblage_report.yml files to create conda environment.

Create two new environments - varcall and varcall_gubbins
```
conda env create -f assemblage/assemblage.yml -n assemblage
conda env create -f assemblage/assemblage_report.yml -n assemblage_report
```

Check installation

```
conda activate assemblage

python assemblage/assemblage.py -h
```

## Input

Input is a directory(-readsdir) containing SE/PE reads and a config file where all the configuration settings for the pipeline are set. This config file settings will be used universally on all samples available in readsdir. An example [config](https://github.com/alipirani88/assembly_umich/blob/master/config) file with default parameters are included in the pipeline folder. You can customize this config file and provide it with the -config argument.

Detailed information in section [Customizing Config file](#customizing-config-file)

Note: Apart from standard Miseq/Hiseq fastq naming extensions (R1_001_final.fastq.gz), other acceptable fastq extensions are: R1.fastq.gz/_R1.fastq.gz, 1_combine.fastq.gz, 1_sequence.fastq.gz, _forward.fastq.gz, _1.fastq.gz/.1.fastq.gz. 

<!---
To generate assembly jobs, you need a filename with fastq read sample names. The script only recognises one filename per line. To generate this filenames input, run the below command. Replace path-to- and PATH-to-save with the path to input reads directory and path to save filenames respectively

```

ls /path-to-/test_readsdir/*_R1_*.fastq.gz | awk -F'/' '{print $(NF)}' > /PATH-to-save/filenames

```

## Steps

Image here
-->


## Quick Start


### Generate and run assembly jobs for a set of PE reads. 


```

python assemblage/assemblage.py -dir /Path-t-/Reads-dir/ -out_dir /Path-to/output-dir/ -type PE -email username@umich.edu -resources "--nodes=1 --ntasks=1 --cpus-per-task=1 --mem=5g --time=50:00:00" -scheduler SLURM -coverage_depth 150 -downsample yes -config assemblage/config

```

 The above command will generate and run assembly jobs for a set of PE reads residing in Reads-dir. The results will be saved in output directory output-dir. 

- Optional - The config file contains options for some frequently used reference genome to use with ABACAS for reordering. To know which reference genomes are included in config file, look up the [config]() file or check the help menu of the pipeline. 

- The assembly will be placed in an individual folder generated for each sample in output directory. 
- A log file for each sample will be generated and can be found in each sample folder inside the output directory. A single log file of this step will be generated in main output directory. For more information on log file prefix and convention, please refer [log](#log) section below.

### Gather and Generate a Multiqc report for the assembly results.

```

conda activate assemblage_report

python assemblage/report.py -out_dir /Path-to/output-dir/


```

## Output

- The report script will gather the assembly fasta file in Results/YYYY-MM-DD_assembly and YYYY-MM-DD_plasmid_assembly folders. It will run and generate QUAST/multiqc results in /Report/Quast/

- Results for each sample can be found in its own individual folder. Each sample folder will contain the assembly fasta file with a suffix \_l500_contigs.fasta and \_l500_plasmid_contigs.fasta. 

- Prokka results can be found in \_prokka directory.

<!---
A script is provided with the pipeline, assembly_jobs.py that will take this filenames and other arguments to generate assembly jobs. To generate assembly jobs for flux, run the below command:

```
/nfs/esnitkin/bin_group/anaconda2/bin/python /nfs/esnitkin/bin_group/scripts/generate_jobs.py -dir /path-to/test_readsdir/ -filenames filenames -out_dir /path-to-output-dir/ -pipeline new_assembly -type PE -email username@umich.edu -resources nodes=1:ppn=4,mem=47000mb,walltime=24:00:00

/nfs/esnitkin/bin_group/anaconda2/bin/python /nfs/esnitkin/bin_group/pipeline/Github/assembly_umich/generate_jobs.py -dir /scratch/esnitkin_fluxod/apirani/varcall_testing/reads_dir/ -out_dir /scratch/esnitkin_fluxod/apirani/varcall_testing/assembly_demo/ -pipeline new_assembly -type PE -email apirani@umich.edu -resources nodes=1:ppn=4,mem=47000mb,walltime=24:00:00

```

Note: Spades assembler requires higher memory cluster and the above resources would be sufficient to run the analysis.

After running the above command, you will find \*.pbs script for each of the sample. You can submit these jobs using a for loop. Before running the loop, make sure the PBS parameters are mentioned correctly.

```
for i in *.pbs; do qsub $i; done
```

or if you want to run it locally:

```
for i in *.pbs; do bash $i; done
```
-->


## Customizing Config file:

By default, the pipeline uses config file that comes with the pipeline. Make sure to edit this config file or copy it to your local system, edit it and provide path of this edited config file with -config argument.

```

cp assembly_umich/config /Path-to-local/config_edit

```

The pipeline implements customisable variant calling configurations using config file. Config file can be customised to use your choice of tools and custom parameters.



If you wish to run the jobs on cluster, make sure you change scheduler parameters in scheduler section shown below: for more information, visit [flux](http://arc-ts.umich.edu/systems-and-services/flux/) homepage.

```

[scheduler]
resources: nodes=1:ppn=4,pmem=4000mb,walltime=24:00:00
email: username@umich.edu
queue: XXX
flux_account: XXX
notification: a

```

Every tool has its own *_bin option where you can set the folder name in which the tool resides. For example, in the below Trimmomatic section example, the Trimmomatic tool resides in /Trimmomatic/ folder that is set with trimmomatic_bin option which in itself resides in /nfs/esnitkin/bin_group/assembly_umich/ folder that was set in binbase option above.

```
[Trimmomatic]
trimmomatic_bin: /Trimmomatic/
adaptor_filepath: adapters/TruSeq3-Nextera_PE_combined.fa
seed_mismatches: 2
palindrome_clipthreshold: 30
simple_clipthreshold: 10
minadapterlength: 8
keep_both_reads: true
window_size: 4
window_size_quality: 20
minlength: 40
headcrop_length: 0
colon: :
targetlength: 125
crop_length: 40
f_p: forward_paired.fq.gz
f_up: forward_unpaired.fq.gz
r_p: reverse_paired.fq.gz
r_up: reverse_unpaired.fq.gz
```

Parameters for each tools can be customised under the 'tool_parameter' attribute of each tool in config file.


For example, to change the minadapterlength parameter of Trimmomatic from 8 to 10, replace minadapterlength of 8 with suppose 10 and restart the pipeline.

## Log:

The pipeline generates a log file following the naming convention: yyyy_mm_dd_hrs_mins_secs_analysisname.log.txt and tracks each event/command. The log file sections follow standard [Python logging conventions](https://docs.python.org/2/howto/logging.html): 

***INFO*** to print STDOUT messages; 

***DEBUG*** to print commands ran by pipeline, 

***ERROR*** to print STDERR messages and 

***EXCEPTION*** to print an exception that occured while the pipeline was running.



















<!---
## Steps:

The different steps of the pipeline are cleaning of reads using [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic), assembling clean reads using [Spades](http://bioinf.spbau.ru/spades)/[Velvet](https://www.ebi.ac.uk/~zerbino/velvet/)(not functional yet) assembler(, assembly evaluation using [Quast](http://bioinf.spbau.ru/quast), contig reordering in case reference genome provided using [ABACAS](http://abacas.sourceforge.net/) and finally annotation using [PROKKA](http://www.vicbioinformatics.com/software.prokka.shtml).

- Step 1: Pre-processing using Trimmomatic
- Step 2: Assembly using Spades/Velvet (Spades assembly steps also involves assembling the plasmids seperately)
- Step 3: Assembly evaluation using QUAST
- Step 4: Contig reordering using ABACAS and Annotation using Prokka


## Usage for single local run:

```

usage: pipeline.py [-h] [-f1 FILE_1] [-f2 FILE_2] -config CONFIG -analysis
                   ANALYSIS_NAME -o OUTPUT_FOLDER [-start_step START_STEP]
                   [-end_step END_STEP] [-A ASSEMBLER] [-type TYPE] [-c CROP]
                   [-reference REFERENCE] [-ariba ARIBA] [-assembly ASSEMBLY]

Assembly pipeline for Illumina SE/PE data

optional arguments:
  -h, --help            show this help message and exit

Required arguments:
  -f1 FILE_1            Paired End file 1
  -config CONFIG        Path to Config file
  -analysis ANALYSIS_NAME
                        Unique analysis name to save the results
  -o OUTPUT_FOLDER      Output Path ending with output directory name to save
                        the results
  -start_step START_STEP
                        Provide the start step. Only 1 works for now.
  -end_step END_STEP    Provide the end step. 2/3/4 All three works
  -A ASSEMBLER          Choose the assembler to assemble the sample reads.
                        Velvet Optimiser or Spades
  -type TYPE            Type of analysis: SE or PE

Optional arguments:
  -f2 FILE_2            Paired End file 2
  -c CROP               choose crop value to crop the reads
  -reference REFERENCE  Provide a reference genome for Abacas Contig ordering
  -ariba ARIBA          Run ariba AMR or MLST on clean reads. expected values:
                        AMR/MLST/BOTH
  -assembly ASSEMBLY    Select one of the following assembly options:
                        "wga":Only Spades Whole Genome Assembly or "plasmid":
                        Only Plasmid Assembly or "Both": Perform both wga and
                        plasmid assembly. Default:wga Options:
                        wga/plasmid/both

```      	

The script can be invoked at any step provided it is supplied with valid -start_step and -end-step flags. 

For example; to run only Trimmomatic on the reads:
***

```
python pipeline.py -f1 PATHtoFile1 -f2 PATHtoFile2 -o path_to_outfolder/ -start_step 1 -end_step 1 -analysis analysis_name -config path_to_config_file -type PE -A spades
```

**Note:**
***

- Before running the pipeline, Make sure the bin directory path and other tool directory path in config file are correct. More Details in section [How to set up config file?](https://github.com/alipirani88/assembly_umich/blob/master/README.md#How to set up config file) below.
- Also edit the reference genome path required for ABACAS reordering. The header name for Reference fasta file should be provided with reference parameter.
- Input file format: Either fastq or gzipped fastq
- Output Directory: Pipeline creates output folder for saving the results. -o option expects path where the output directory is required to be created.
- Assembler: -A option expects the name of assembler. Either Spades or Velvet. (Velvet is not tested)
 
**How to set up config file?**
***

The config file used for this pipeline is a YAML type file containing specific details such as path to your bin directory, tools and reference fasta file, parameters used for each tools and other system details.

This config file makes it easy to control various parameters used in the entire pipeline for different tools in a single file. An example config file is included in the project.

The path to your bin directory where all the tools required for this pipeline are installed should be specified under the section [bin_path] and variable 'binbase'.

```
[bin_path]
binbase = /home/apirani/bin/assembly_umich/bin/
```
Change the '_bin' variable under each tool section accordingly to the folder name of each tool. e.g: If trimmomatic was installed in a bin directory specified under bin_path section by the name 'Trimmomatic', then the variable 'trimmomatic_bin' should be changed to '/Trimmomatic/' 

```
[Trimmomatic]
trimmomatic_bin = /Trimmomatic/
```

In a similar fashion, the reference genome can be specified in the following way:

```
[KPNIH1]
ref_name = KPNIH1.fasta
ref_path = /home/apirani/bin/assembly_umich/bin/reference/KPNIH1/
```

Here, the main header section [KPNIH1] represents the title for the reference genome. This title is required with the parameter -reference while running the pipeline for contig reordering step. 'ref_name' is the reference fasta filename saved under the path/dir 'ref_path'

-->

