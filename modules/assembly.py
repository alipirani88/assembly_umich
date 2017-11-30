__author__ = 'alipirani'
from config_settings import ConfigSectionMap
from modules.log_modules import keep_logging
from modules.logging_subprocess import *
import os

# Prepare which assembler to use based on user-provided argument.
def assembly(forward_paired, reverse_paired, forward_unpaired, reverse_unpaired, assembler, out_path, logger, Config, do_assembly):
    # check the assembler argument provided and call relevant subroutines
    if assembler:
        if assembler == "velvet":
            # Have to insert the Velvet Functionality later.
            print ""
            #velvetoptimiser(forward_paired, reverse_paired, forward_unpaired, reverse_unpaired, out_path)
        elif assembler == "spades":
            (contigs, scaffolds, plasmid_contigs, plasmid_scaffolds) = spades_assembly(forward_paired, reverse_paired, forward_unpaired, reverse_unpaired, out_path, logger, Config, do_assembly)
            return contigs, scaffolds, plasmid_contigs, plasmid_scaffolds
    else:
        keep_logging('Please provide the assembler argument for this step.', 'Please provide the assembler argument for this step.', logger, 'exception')
        exit()

# Assembly using Spades
def spades_assembly(forward_paired, reverse_paired, forward_unpaired, reverse_unpaired, out_path, logger, Config, do_assembly):
    # check if the clean reads from Trimmomatic exists in the output folder.
    # Set the paired and unpaired string constants based on their availability
    (paired, unpaired) = check_cleanreads(forward_paired, reverse_paired, forward_unpaired, reverse_unpaired)
    # Pending Changes
    if paired == "0" and unpaired == "0":
        # Clean Paired and unpaired reads doesn't exist. Take raw Input PE files for assembly
        message = "No clean Paired and unpaired reads. Considering forward_paired and reverse_paired as raw Fastq files for assembly.\n"
        print message
        cmdstring = ConfigSectionMap("bin_path")['binbase'] + ConfigSectionMap("spades")['spades_bin'] + ConfigSectionMap("spades")['base_cmd'] + " --pe1-1 " + forward_paired + " --pe1-2 " + reverse_paired + " -o " + out_path + "spades_results " + ConfigSectionMap("spades")['spades_parameters']
        plasmid_cmdstring = ConfigSectionMap("bin_path")['binbase'] + ConfigSectionMap("spades")['spades_bin'] + ConfigSectionMap("spades")['base_cmd'] + " --pe1-1 " + forward_paired + " --pe1-2 " + reverse_paired + " -o " + out_path + "spades_plasmid_results " + ConfigSectionMap("spades")['plasmid_spades_parameters']
        print "Running: %s \n" % cmdstring
        print "Running: %s \n" % plasmid_cmdstring
        os.system(cmdstring)
        os.system(plasmid_cmdstring = ConfigSectionMap("bin_path")['binbase'] + ConfigSectionMap("spades")['spades_bin'] + ConfigSectionMap("spades")['base_cmd'] + " --pe1-1 " + forward_paired + " --pe1-2 " + reverse_paired + " -o " + out_path + "spades_plasmid_results " + ConfigSectionMap("spades")['plasmid_spades_parameters'])
        print "Spades assembly results can be found in " + out_path + "spades_results"
        print "plasmid Spades assembly results can be found in " + out_path + "spades_plasmid_results"
        contigs = out_path + "spades_results" + "/contigs.fasta"
        scaffolds = out_path + "spades_results" + "/scaffolds.fasta"
        plasmid_contigs = out_path + "spades_plasmid_results" + "/contigs.fasta"
        plasmid_scaffolds = out_path + "spades_plasmid_results" + "/contigs.fasta"
        # Copy final contigs/scaffolds file to output directory
        cp_cmdstring = "cp %s %s %s" % (contigs, scaffolds, out_path)
        os.system(cp_cmdstring)
        print "\n################## End: SPADES ASSEMBLY ##################\n"
        return contigs, scaffolds
    # Pending Changes
    elif paired == "1" and unpaired == "0":
        # Only clean Paired PE files exists. Take these files for assembly input.
        message = "Taking only paired reads for assembly.\n"
        print message
        if reverse_paired == "None" and reverse_unpaired == "None":
            cmdstring = ConfigSectionMap("bin_path", Config)['binbase'] + ConfigSectionMap("spades", Config)['spades_bin'] + ConfigSectionMap("spades", Config)['base_cmd'] + " --s1 " + forward_paired + " -o " + out_path + "spades_results/ " + ConfigSectionMap("spades", Config)['spades_parameters']
            plasmid_cmdstring = ConfigSectionMap("bin_path", Config)['binbase'] + ConfigSectionMap("spades", Config)['spades_bin'] + ConfigSectionMap("spades", Config)['base_cmd'] + " --s1 " + forward_paired + " -o " + out_path + "spades_plasmid_results/ " + ConfigSectionMap("spades", Config)['plasmid_spades_parameters']
            print "Running: %s \n" % cmdstring
            print "Running: %s \n" % plasmid_cmdstring
            os.system(cmdstring)
            os.system(plasmid_cmdstring)
            print "Spades assembly results can be found in " + out_path + "spades_results"
            print "plasmid Spades assembly results can be found in " + out_path + "spades_plasmid_results"
            contigs = out_path + "spades_results" + "/contigs.fasta"
            scaffolds = out_path + "spades_results" + "/scaffolds.fasta"
            plasmid_contigs = out_path + "spades_plasmid_results" + "/contigs.fasta"
            plasmid_scaffolds = out_path + "spades_plasmid_results" + "/contigs.fasta"
            # Copy final contigs/scaffolds file to output directory
            cp_cmdstring = "cp %s %s %s" % (contigs, scaffolds, out_path)
            os.system(cp_cmdstring)
            print "\n################## End: SPADES ASSEMBLY ##################\n"
        else:
            ##pending changes
            cmdstring = ConfigSectionMap("bin_path")['binbase'] + ConfigSectionMap("spades")['spades_bin'] + ConfigSectionMap("spades")['base_cmd'] + " --pe1-1 " + forward_paired + " --pe1-2 " + reverse_paired + " -o " + out_path + "spades_results/ " + ConfigSectionMap("spades")['spades_parameters']
            plasmid_cmdstring = ConfigSectionMap("bin_path")['binbase'] + ConfigSectionMap("spades")['spades_bin'] + ConfigSectionMap("spades")['base_cmd'] + " --pe1-1 " + forward_paired + " --pe1-2 " + reverse_paired + " -o " + out_path + "spades_plasmid_results/ " + ConfigSectionMap("spades")['plasmid_spades_parameters']
            print "Running: %s \n" % cmdstring
            print "Running: %s \n" % plasmid_cmdstring
            os.system(cmdstring)
            os.system(plasmid_cmdstring)
            print "Spades assembly results can be found in " + out_path + "spades_results"
            print "plasmid Spades assembly results can be found in " + out_path + "spades_plasmid_results"
            contigs = out_path + "spades_results" + "/contigs.fasta"
            scaffolds = out_path + "spades_results" + "/scaffolds.fasta"
            plasmid_contigs = out_path + "spades_plasmid_results" + "/contigs.fasta"
            plasmid_scaffolds = out_path + "spades_plasmid_results" + "/contigs.fasta"
            # Copy final contigs/scaffolds file to output directory
            cp_cmdstring = "cp %s %s %s" % (contigs, scaffolds, out_path)
            os.system(cp_cmdstring)
            print "\n################## End: SPADES ASSEMBLY ##################\n"
        return contigs, scaffolds, plasmid_contigs, plasmid_scaffolds
    # Pending Changes
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
        return contigs, scaffolds, plasmid_contigs, plasmid_scaffolds
    else:
        # Clean paired and unpaired files exists. Take all these files as input.
        cmdstring = ConfigSectionMap("bin_path", Config)['binbase'] + ConfigSectionMap("spades", Config)['spades_bin'] + ConfigSectionMap("spades", Config)['base_cmd'] + " " + ConfigSectionMap("spades", Config)['spades_parameters'] + " --pe1-1 " + forward_paired + " --pe1-2 " + reverse_paired + " --pe1-s " + forward_unpaired + " --pe1-s " + reverse_unpaired + " -o " + out_path + "spades_results"
        plasmid_cmdstring = ConfigSectionMap("bin_path", Config)['binbase'] + ConfigSectionMap("spades", Config)['spades_bin'] + ConfigSectionMap("spades", Config)['base_cmd'] + " " + ConfigSectionMap("spades", Config)['plasmid_spades_parameters'] + " --pe1-1 " + forward_paired + " --pe1-2 " + reverse_paired + " --pe1-s " + forward_unpaired + " --pe1-s " + reverse_unpaired + " -o " + out_path + "spades_plasmid_results"


        if do_assembly == "both":
            keep_logging('Running Spades and plasmid Spades assembly', 'Running Spades and plasmid Spades assembly', logger, 'debug')
            try:
                keep_logging(cmdstring, cmdstring, logger, 'debug')
                call(cmdstring, logger)
                #Check if they are empty
                contigs = out_path + "spades_results" + "/contigs.fasta"
                scaffolds = out_path + "spades_results" + "/scaffolds.fasta"
                # Copy final contigs/scaffolds file to output directory
                cp_cmdstring = "cp %s %s %s" % (contigs, scaffolds, out_path)
                os.system(cp_cmdstring)
                print ""
                keep_logging('Spades assembly results can be found in {}spades_results'.format(out_path), 'Spades assembly results can be found in {}spades_results'.format(out_path), logger, 'info')
            except sp.CalledProcessError:
                keep_logging('Error in Spades Assembly step. Exiting. Please check spades.log file in spades_results folder', 'Error in Spades Assembly step. Exiting. Please check spades.log file in spades_results folder', logger, 'exception')
                sys.exit(1)


            try:
                keep_logging(plasmid_cmdstring, plasmid_cmdstring, logger, 'debug')
                call(plasmid_cmdstring, logger)
                print ""
                plasmid_contigs = out_path + "spades_plasmid_results" + "/contigs.fasta"
                plasmid_scaffolds = out_path + "spades_plasmid_results" + "/contigs.fasta"
                keep_logging('Spades plasmid assembly results can be found in {}spades_plasmid_results'.format(out_path), 'Spades plasmid assembly results can be found in {}spades_plasmid_results'.format(out_path), logger, 'info')
            except sp.CalledProcessError:
                keep_logging('Error in Spades Plasmid Assembly step. Exiting. Please check spades.log file in spades_results folder', 'Error in Spades Plasmid Assembly step. Exiting. Please check spades.log file in spades_results folder', logger, 'exception')
                sys.exit(1)

        if do_assembly == "wga":
            keep_logging('Running Spades assembly', 'Running Spades assembly', logger, 'debug')
            try:
                keep_logging(cmdstring, cmdstring, logger, 'debug')
                call(cmdstring, logger)
                #Check if they are empty
                contigs = out_path + "spades_results" + "/contigs.fasta"
                scaffolds = out_path + "spades_results" + "/scaffolds.fasta"
                # Copy final contigs/scaffolds file to output directory
                cp_cmdstring = "cp %s %s %s" % (contigs, scaffolds, out_path)
                os.system(cp_cmdstring)
                print ""
                keep_logging('Spades assembly results can be found in {}spades_results'.format(out_path), 'Spades assembly results can be found in {}spades_results'.format(out_path), logger, 'info')
            except sp.CalledProcessError:
                keep_logging('Error in Spades Assembly step. Exiting. Please check spades.log file in spades_results folder', 'Error in Spades Assembly step. Exiting. Please check spades.log file in spades_results folder', logger, 'exception')
                sys.exit(1)

            plasmid_contigs = out_path + "spades_plasmid_results" + "/contigs.fasta"
            plasmid_scaffolds = out_path + "spades_plasmid_results" + "/contigs.fasta"


        if do_assembly == "plasmid":
            keep_logging('Running plasmid Spades assembly', 'Running plasmid Spades assembly', logger, 'debug')
            try:
                keep_logging(plasmid_cmdstring, plasmid_cmdstring, logger, 'debug')
                call(plasmid_cmdstring, logger)
                print ""
                plasmid_contigs = out_path + "spades_plasmid_results" + "/contigs.fasta"
                plasmid_scaffolds = out_path + "spades_plasmid_results" + "/contigs.fasta"
                keep_logging('Spades plasmid assembly results can be found in {}spades_plasmid_results'.format(out_path), 'Spades plasmid assembly results can be found in {}spades_plasmid_results'.format(out_path), logger, 'info')
            except sp.CalledProcessError:
                keep_logging('Error in Spades Plasmid Assembly step. Exiting. Please check spades.log file in spades_results folder', 'Error in Spades Plasmid Assembly step. Exiting. Please check spades.log file in spades_results folder', logger, 'exception')
                sys.exit(1)

            contigs = out_path + "spades_results" + "/contigs.fasta"
            scaffolds = out_path + "spades_results" + "/scaffolds.fasta"

        return contigs, scaffolds, plasmid_contigs, plasmid_scaffolds

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


# check if the paired and unpaired clean reads exists in output directory. Set the paired and unpaired constants accordingly.
def check_cleanreads(forward_paired, reverse_paired, forward_unpaired, reverse_unpaired):
    # Initialize the constants
    unpaired = "1"
    paired = "1"
    if not os.path.isfile(forward_paired):
        print "The pre-processed paired reads file does not exists. This step requires clean reads namely: forward_paired.fq.gz, reverse_paired.fq.gz.\n"
        paired = "0"
    if not os.path.isfile(forward_unpaired):
        print "The pre-processed unpaired reads file does not exists. This step requires (optional) unpaired clean reads named: forward_unpaired.fq.gz, reverse_unpaired.fq.gz.\n"
        unpaired = "0"
    return paired, unpaired



