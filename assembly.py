__author__ = 'alipirani'


from check_subroutines import *
from config_settings import ConfigSectionMap

# Prepare which assembler to use based on user-provided argument.
def assembly(forward_paired, reverse_paired, forward_unpaired, reverse_unpaired, assembler, out_path):
    # Make sure the output folder exists or create at given path
    make_sure_path_exists(out_path)
    # check the assembler argument provided and call relevant subroutines
    if assembler:
            if assembler == "velvet":
                print "Assembling the input files using Velvet.\n"
                velvetoptimiser(forward_paired, reverse_paired, forward_unpaired, reverse_unpaired, out_path)
            elif assembler == "spades":
                print "Assembling the input files using SPADES assembler.\n"
                spades_assembly(forward_paired, reverse_paired, forward_unpaired, reverse_unpaired, out_path)
    else:
        print "Please provide the assembler argument for this step.\n"
        usage()
        exit()

# Assembly using Spades
def spades_assembly(forward_paired, reverse_paired, forward_unpaired, reverse_unpaired, out_path):
    print "\n################## Running SPADES on input files ##################\n"
    # check if the clean reads from Trimmomatic exists in the output folder.
    # Set the paired and unpaired string constants based on their availability
    (paired, unpaired) = check_cleanreads(forward_paired, reverse_paired, forward_unpaired, reverse_unpaired)
    if paired == "0" and unpaired == "0":
        # Clean Paired and unpaired reads doesn't exist. Take raw Input PE files for assembly
        message = "No clean Paired and unpaired reads. Considering forward_paired and reverse_paired as raw Fastq files for assembly.\n"
        print message
        cmdstring = ConfigSectionMap("bin_path")['binbase'] + ConfigSectionMap("spades")['spades_bin'] + ConfigSectionMap("spades")['base_cmd'] + " --pe1-1 " + forward_paired + " --pe1-2 " + reverse_paired + " -o " + out_path + "spades_results " + ConfigSectionMap("spades")['spades_parameters']
        print "Running: %s \n" % cmdstring
        os.system(cmdstring)
        print "Spades assembly results can be found in " + out_path + "spades_results"
        contigs = out_path + "spades_results" + "/contigs.fasta"
        scaffolds = out_path + "spades_results" + "/scaffolds.fasta"
        # Copy final contigs/scaffolds file to output directory
        cp_cmdstring = "cp %s %s %s" % (contigs, scaffolds, out_path)
        os.system(cp_cmdstring)
        print "\n################## End: SPADES ASSEMBLY ##################\n"
        return contigs, scaffolds
    elif paired == "1" and unpaired == "0":
        # Only clean Paired PE files exists. Take these files for assembly input.
        message = "Taking only paired reads for assembly.\n"
        print message
        cmdstring = ConfigSectionMap("bin_path")['binbase'] + ConfigSectionMap("spades")['spades_bin'] + ConfigSectionMap("spades")['base_cmd'] + " --pe1-1 " + forward_paired + " --pe1-2 " + reverse_paired + " -o " + out_path + "spades_results/ " + ConfigSectionMap("spades")['spades_parameters']
        print "Running: %s \n" % cmdstring
        os.system(cmdstring)
        print "Spades assembly results can be found in " + out_path + "spades_results"
        contigs = out_path + "spades_results" + "/contigs.fasta"
        scaffolds = out_path + "spades_results" + "/scaffolds.fasta"
        # Copy final contigs/scaffolds file to output directory
        cp_cmdstring = "cp %s %s %s" % (contigs, scaffolds, out_path)
        os.system(cp_cmdstring)
        print "\n################## End: SPADES ASSEMBLY ##################\n"
        return contigs, scaffolds
    elif paired == "0" and unpaired == "1":
        # Only clean unpaired PE files exists. Pending...
        cmdstring = "This can be single reads......"
        print "Running: %s \n" % cmdstring
        os.system(cmdstring)
        print "Spades assembly results can be found in " + out_path + "spades_results"
        contigs = out_path + "spades_results" + "/contigs.fasta"
        scaffolds = out_path + "spades_results" + "/scaffolds.fasta"
        # Copy final contigs/scaffolds file to output directory
        cp_cmdstring = "cp %s %s %s" % (contigs, scaffolds, out_path)
        os.system(cp_cmdstring)
        print "\n################## End: SPADES ASSEMBLY ##################\n"
        return contigs, scaffolds
    else:
        # Clean paired and unpaired files exists. Take all these files as input.
        cmdstring = ConfigSectionMap("bin_path")['binbase'] + ConfigSectionMap("spades")['spades_bin'] + ConfigSectionMap("spades")['base_cmd'] + " " + ConfigSectionMap("spades")['spades_parameters'] + " --pe1-1 " + forward_paired + " --pe1-2 " + reverse_paired + " --pe1-s " + forward_unpaired + " --pe1-s " + reverse_unpaired + " -o " + out_path + "spades_results"
        print "Running Spades with all input file parameters:\n"
        print "Running: %s \n" % cmdstring
        os.system(cmdstring)
        print "Spades assembly results can be found in " + out_path + "spades_results"
        contigs = out_path + "spades_results" + "/contigs.fasta"
        scaffolds = out_path + "spades_results" + "/scaffolds.fasta"
        # Copy final contigs/scaffolds file to output directory
        cp_cmdstring = "cp %s %s %s" % (contigs, scaffolds, out_path)
        os.system(cp_cmdstring)
	directory = os.path.dirname(out_path)
        filename = os.path.basename(directory)
        contig_l500 = "/home2/apirani/bin/bioawk-master/bioawk -c fastx '{ if(length($seq) > 500) { print \">\"$name; print $seq }}' %s > %s/%s_l500_contigs.fasta" % (contigs, out_path, filename)
        scaffold_l500 = "/home2/apirani/bin/bioawk-master/bioawk -c fastx '{ if(length($seq) > 500) { print \">\"$name; print $seq }}' %s > %s/%s_l500_scaffolds.fasta" % (scaffolds, out_path, filename)
        os.system(contig_l500)
        os.system(scaffold_l500)
        print "\n################## End: SPADES ASSEMBLY ##################\n"
        #return contigs, scaffolds




# Assembly using VelvetOptimiser
def velvetoptimiser(forward_paired, reverse_paired, forward_unpaired, reverse_unpaired, out_path):
    print "\n################## Running VELVET on input files ##################\n"
    velvet_dir = out_path + "velvet_results"
    Vforward_paired = "-shortPaired -fastq.gz " + forward_paired
    Vforward_unpaired = " -short -fastq.gz " + forward_unpaired
    Vreverse_paired = " -shortPaired2 -fastq.gz " + reverse_paired
    Vreverse_unpaired = " -short2 -fastq.gz " + reverse_unpaired
    (paired, unpaired) = check_cleanreads(forward_paired, reverse_paired, forward_unpaired, reverse_unpaired)
    contigs = out_path + "/contigs.fa"
    scaffolds = ""
    if paired == 0 and unpaired == 0:
        # Clean Paired and unpaired reads doesn't exist. Take raw Input PE files for assembly
        message = "No clean Paired and unpaired reads. Considering forward_paired and reverse_paired as raw Fastq files for assembly.\n"
        print message
        cmdstring = ConfigSectionMap("bin_path")['binbase'] + "VelvetOptimiser/VelvetOptimiser.pl -s 71 -e 121 -x 20" + " --d " + velvet_dir + " -f '" + Vforward_paired + " " + Vreverse_paired + "\'"
        print "Running: %s \n" % cmdstring
        os.system(cmdstring)
        # Copy final contigs/scaffolds file to output directory
        cp_cmdstring = "cp " + velvet_dir + "/contigs.fa " + out_path
        os.system(cp_cmdstring)
        print "\n################## END: VELVET ASSEMBLY ##################\n"
        return contigs, scaffolds

    elif paired == 1 and unpaired == 0:
        # Only clean Paired PE files exists. Take these files for assembly input.
        message = "Taking only paired reads for assembly.\n"
        print message
        cmdstring = ConfigSectionMap("bin_path")['binbase'] + "VelvetOptimiser/VelvetOptimiser.pl -s 71 -e 121 -x 20" + " --d " + velvet_dir + " -f '" + Vforward_paired + Vreverse_paired + "\'"
        print "Running: %s \n" % cmdstring
        os.system(cmdstring)
        # Copy final contigs/scaffolds file to output directory
        cp_cmdstring = "cp " + velvet_dir + "/contigs.fa " + out_path
        os.system(cp_cmdstring)
        print "\n################## END: VELVET ASSEMBLY ##################\n"
        return contigs, scaffolds

    elif paired == 0 and unpaired == 1:
        # Only clean unpaired PE files exists. Pending...
        cmdstring = "This can be single reads......"
        print "Running: %s \n" % cmdstring
        #return contigs, scaffolds
        #os.system(cmdstring)
    else:
        # Clean paired and unpaired files exists. Take all these files as input.
        os.chdir(out_path)
        cmdstring = ConfigSectionMap("bin_path")['binbase'] + "VelvetOptimiser/VelvetOptimiser.pl -s 71 -e 121 -x 20" + " --d " + velvet_dir + " -f '" + Vforward_paired + " " + Vreverse_paired + " " + Vforward_unpaired + " " + Vreverse_unpaired + "\'"
        print "Running with all input file parameters.\n"
        print "Running: %s \n" % cmdstring
        os.system(cmdstring)
        # Copy final contigs/scaffolds file to output directory
        cp_cmdstring = "cp " + velvet_dir + "/contigs.fa " + out_path
        os.system(cp_cmdstring)
        print "\n################## END: VELVET ASSEMBLY ##################\n"
        return contigs, scaffolds


