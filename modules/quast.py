__author__ = 'alipirani'
from config_settings import ConfigSectionMap
import os
from modules.log_modules import *
from modules.logging_subprocess import *


#Assembled contigs evaluation using QUAST
def quast_evaluation(out_path, final_l500_contig, final_l500_plasmid_contig, logger, Config):
    if not os.path.isfile(final_l500_contig):
        print ("Can't find scaffolds file in output folder. \nSpades:Please check the availability of scaffolds.fa file in the given output folder or change the name of contigs file to contigs.fasta\nVelvet:Please check the availability of contigs.fa file in the given output folder or change the name of contigs file to contigs.fa\n")
        #exit()
    if not os.path.isfile(final_l500_plasmid_contig):
        print ("Can't find scaffolds file in output folder. \nSpades:Please check the availability of scaffolds.fa file in the given output folder or change the name of contigs file to contigs.fasta\nVelvet:Please check the availability of contigs.fa file in the given output folder or change the name of contigs file to contigs.fa\n")
        #exit()

    quast_out = out_path + "/quast_results"
    quast_plasmid_out = out_path + "/quast_plasmid_results"

    #cmdstring = ConfigSectionMap("bin_path", Config)['binbase'] + ConfigSectionMap("quast", Config)['quast_bin'] + ConfigSectionMap("quast", Config)['base_cmd'] + " " + contigfile + " " + scaffoldfile + " " + plasmid_contigs_file + " " + plasmid_scaffolds_file + " " + final_l500_contig + " " + final_l500_plasmid_contig + " -o " + quast_out + " -t 0,1000,5000,10000,25000,50000"
    cmd_assembly = ConfigSectionMap("quast", Config)['base_cmd'] + " " + final_l500_contig + " -o " + quast_out + " -t 0,1000,5000,10000,25000,50000"
    plasmid_cmd_assembly = ConfigSectionMap("bin_path", Config)['binbase'] + ConfigSectionMap("quast", Config)['quast_bin'] + ConfigSectionMap("quast", Config)['base_cmd'] + " " + final_l500_plasmid_contig + " -o " + quast_plasmid_out + " -t 0,1000,5000,10000,25000,50000"
    try:
        keep_logging(cmd_assembly, cmd_assembly, logger, 'debug')
        call(cmd_assembly, logger)
        print ("")
    except sp.CalledProcessError:
        keep_logging('Error in Quast Evaluation. Check log file.', 'Error in Quast Evaluation. Check log file.', logger, 'exception')
        sys.exit(1)
    try:
        keep_logging(plasmid_cmd_assembly, plasmid_cmd_assembly, logger, 'debug')
        call(plasmid_cmd_assembly, logger)
        print ("")
    except sp.CalledProcessError:
        keep_logging('Error in Quast Evaluation. Check log file.', 'Error in Quast Evaluation. Check log file.', logger, 'exception')
        sys.exit(1)