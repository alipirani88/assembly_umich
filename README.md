# Assembly Pipeline for illumina PE reads

***The pipeline takes Illumina PE FastQ reads as input for various steps of pre-processing, assembly, evaluation and annotation.***

***Steps:***

-1. Pre-processing using Trimmomatic
-2. Assembly using Spades/Velvet
-3. Assembly evaluation using QUAST
-4. Abacas contig reordering and Prokka Annotation

usage: pipeline.py [-h] [-f1 FILE_1] [-f2 FILE_2] -config CONFIG -analysis
                   ANALYSIS_NAME -o OUTPUT_FOLDER [-start_step START_STEP]
                   [-end_step END_STEP] [-A ASSEMBLER] [-type TYPE] [-c CROP]
                   [-reference REFERENCE]
```
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
```      	

The script can be invoked at any step provided it is supplied with valid -start_step and -end-step flags. 

For e.g: To run only Trimmomatic on the reads, the valid options are:

```
python pipeline.py -f1 PATHtoFile1 -f2 PATHtoFile2 -o path_to_outfolder/ -start_step 1 -end_step 1 -analysis analysis_name -config path_to_config_file -type PE -A spades
```

***Note:***

Before running pipeline, Make sure the bin directory path and other tool directory path in config file are correct.

Input file format: Either fastq or gzipped fastq

Output Directory: Pipeline creates output folder for saving the results. -o option expects path where you want to create the results output file. A directory by the name provided with -analysis argument will be created in path provided with -o argument.

Assembler: -A option expects the name of assembler. Either Spades or Velvet.


         
    
