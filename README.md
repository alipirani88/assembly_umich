# De Novo Microbial Genome Assembly Pipeline

This pipeline takes Illumina PE FastQ reads as input for various steps of pre-processing, assembly, evaluation, assembly improvement and annotation steps.

The different steps of the pipeline involves quality control and cleaning of reads using [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) & [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic), assembling the clean reads using [Spades](http://bioinf.spbau.ru/spades)/[Velvet](https://www.ebi.ac.uk/~zerbino/velvet/) assembler, assembly evaluation using [Quast](http://bioinf.spbau.ru/quast), contig reordering in case reference genome provided using [ABACAS](http://abacas.sourceforge.net/) and finally annotation using [PROKKA](http://www.vicbioinformatics.com/software.prokka.shtml).

**Steps:**
    
- Step 1: Pre-processing using Trimmomatic
- Step 2: Assembly using Spades/Velvet (Spades assembly steps also involves assembling the plasmids seperately)
- Step 3: Assembly evaluation using QUAST
- Step 4: Contig reordering using ABACAS and Annotation using Prokka

(Assembly Improvement using REAPR removed)

**Usage:** 

```
python pipeline.py [-h] [-f1 FILE_1] [-f2 FILE_2] [-o OUTPUT_FOLDER] [-start_step START] [-end_step END] [--qa] [-A Assembler] [-reference Reference_File_Name]
```

**Optional Arguments:**

        -h, --help		show this help message and exit
        
        -v, --version         	show program's version number and exit
        
        -f1 FILE_1            	Paired End fastq file 1
        
        -f2 FILE_2            	Paired End fastq file 2
        
        -o OUTPUT_FOLDER        Output Path to save the results
        
        -start_step START_STEP  Provide the start step
        
        -end_step END_STEP    	Provide the end step
        
        -A ASSEMBLER          	Choose an assembler to assemble the sample reads.
                                Velvet Optimiser or Spades(prefered)
                                
        --qa                  	Run FastQC for quality check
        
        -reference              Provide reference genome in case you select step 4 that involves Abacas contig reordering and Annotation

   Optional:
            
            FastQC: FastQC can be run on raw sample data using --qa flag.
        	

***The script pipeline.py can be invoked at any step provided it is supplied with valid -start_step and -end-step flags. 
For e.g: To run only Trimmomatic on the reads, the valid options are:***

```
python pipeline.py -f1 PATHtoFile1 -f2 PATHtoFile2 -o path_to_outfolder/Output_Foldername/ -start_step 1 -end_step 1
```

**Note:**

- Before running the pipeline, Make sure the bin directory path and other tool directory path in config file are correct. More Details in section How to set up config file?  
- Also edit the reference genome path required for ABACAS reordering. The header name for Reference fasta file should be provided with reference parameter.
- Input file format: Either fastq or gzipped fastq
- Output Directory: Pipeline creates output folder for saving the results. -o option expects path followed by output directory name. The name provided at the end of directory will be used to create a new directory for saving results.
- Assembler: -A option expects the name of assembler. Either Spades or Velvet. (Velvet is not tested)
        



## How to set up config file?

The config file used for this pipeline is a YAML type file containing specific details such as path to your bin directory, tools and reference fasta file, parameters used for each tools and other system details.

This config file makes it easy to control various parameters used in the entire pipeline for different tools in a single file. An example config file is included in the project.

The path to your bin directory where all the tools required for this pipeline are installed should be specified under the section [bin_path] and variable 'binbase'.

Change the '_bin' variable under each tool section accordingly to the folder name of each tool. e.g: If trimmomatic was installed in a bin directory specified under bin_path section by the name 'Trimmomatic', then the variable 'trimmomatic_bin' should be changed to '/Trimmomatic/' 

In a similar fashion, the reference genome can be specified in the following way:

```
[KPNIH1]
ref_name = KPNIH1.fasta
ref_path = /home/apirani/bin/assembly_umich/bin/reference/KPNIH1/
```

Here, the main header section [KPNIH1] represents the title for the reference genome. This title is required with the parameter -reference while running the pipeline for contig reordering step. 'ref_name' is the reference fasta filename saved under the path/dir 'ref_path'
