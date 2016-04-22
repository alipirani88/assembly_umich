# Bacterial Genome Assembly Pipeline

The pipeline takes Illumina PE FastQ reads as input for various steps of pre-processing, assembly, evaluation and assembly improvement steps.

**Steps:**
    
- Pre-processing using Trimmomatic
- Assembly using Spades/Velvet
- Assembly evaluation using QUAST
- Contig reordering using ABACAS
- Annotation using Prokka

Usage: 

```
python pipeline.py [-h] [-f1 FILE_1] [-f2 FILE_2] [-o OUTPUT_FOLDER] [-start_step START] [-end_step END] [--qa] [-A Assembler] [--err]
```

Optional Arguments:

        -h, --help		show this help message and exit
        
        -v, --version         	show program's version number and exit
        
        -f1 FILE_1            	Paired End fastq file 1
        
        -f2 FILE_2            	Paired End fastq file 2
        
        -O OUTPUT_FOLDER        Output Path to save the results
        
        -start_step START_STEP  Provide the start step
        
        -end_step END_STEP    	Provide the end step
        
        -A ASSEMBLER          	Choose the assembler to assemble the sample reads.
                                Velvet Optimiser or Spades
                                
        --qa                  	Run FastQC for quality check
        
        --err                   Run Spades using built-in BayesHammer Error Corrector
  
  
   Optional:
            
            FastQC: FastQC can be run on raw sample data using --qa flag.
        	

The script can be invoked at any step provided it is supplied with valid -start_step and -end-step flags. 
For e.g: To run only Trimmomatic on the reads, the valid options are:

    python pipeline.py -f1 PATHtoFile1 -f2 PATHtoFile2 -O path_to_outfolder/Output_Foldername/ -start_step 1 -end_step 1

    Note:

        Before running pipeline, Make sure the bin directory path and other tool directory path in config file are correct. Also edit the reference genome path required for ABACAS reordering.
        
        Input file format: Either fastq or gzipped fastq
        
        Output Directory: Pipeline creates output folder for saving the results. -o option expects path followed by output directory name. The name provided at the end of directory will be used to create a new directory for saving results.
        
        Assembler: -A option expects the name of assembler. Either Spades or Velvet. (Velvet is not tested)
        
         
    
