__author__ = 'alipirani'

from config_settings import ConfigSectionMap
from modules.check_subroutines import *

def reapr(forward_paired, reverse_paired, out_path, contigs, scaffoldfile):
    print "\n################## Misassembly detection using Reapr. ##################\n"
    try:
        cmdstring = ConfigSectionMap("bin_path")['binbase'] + ConfigSectionMap("reapr")['reapr_bin'] + ConfigSectionMap("reapr")['base_cmd'] + " facheck " + contigs
        print "Running: %s" % cmdstring
        os.system(cmdstring)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            print "Errors in Final Assembly!\n"
            exit()
    reapr_dir = out_path + "/reapr_results"
    # Create a new directory to save reapr results
    make_sure_path_exists(reapr_dir)
    # decompress the input files for reapr input
    decom_cmdstring1 = "gunzip -c " + forward_paired + " > " + reapr_dir + "/forward.fastq"
    decom_cmdstring2 = "gunzip -c " + reverse_paired + " > " + reapr_dir + "/reverse.fastq"
    os.system(decom_cmdstring1)
    os.system(decom_cmdstring2)
    # Map reads to contigs/scaffolds file.
    map_cmdstring = ConfigSectionMap("bin_path")['binbase'] + ConfigSectionMap("reapr")['reapr_bin'] + ConfigSectionMap("reapr")['base_cmd'] + " smaltmap " + contigs + " " +  reapr_dir + "/forward.fastq" + " " + reapr_dir + "/reverse.fastq" + " " + reapr_dir + "/mapped.bam"
    print "********Running Smalt map:\n"
    print "Running: %s" % map_cmdstring
    os.system(map_cmdstring)
    # Run Perfectfrommap  to generate perfect and uniquely mapping read coverage, for input into the REAPR pieline.
    print "Running PerfectfromBAM:\n"
    perfect_cmdstring = ConfigSectionMap("bin_path")['binbase'] + ConfigSectionMap("reapr")['reapr_bin'] + ConfigSectionMap("reapr")['base_cmd'] + " perfectfrombam " + reapr_dir + "/mapped.bam " + ConfigSectionMap("reapr")['perfectfrombam_parameters']
    print "Running: %s" % perfect_cmdstring
    os.system(perfect_cmdstring)
    # Run reapr pipeline to detect misassembly
    print "********Running reapr pipeline:\n"
    reapr_cmdstring = ConfigSectionMap("bin_path")['binbase'] + ConfigSectionMap("reapr")['reapr_bin'] + ConfigSectionMap("reapr")['base_cmd'] + " pipeline " + contigs + " " + reapr_dir + "/mapped.bam " + reapr_dir + "/final_reapr perfect"
    print "Running: %s" % reapr_cmdstring
    os.system(reapr_cmdstring)
    # remove decompressed fastq file; no longer needed for reapr
    rm_string = "rm " + reapr_dir + "/*.fastq"
    os.system(rm_string)
    print "The final Reapr results are in the output directory:\n"
    print reapr_dir + "/final_reapr"
    print "\n################## END: Misassembly detection using Reapr. ##################\n"
