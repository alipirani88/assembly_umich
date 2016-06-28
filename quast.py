__author__ = 'alipirani'
from config_settings import ConfigSectionMap
import os

#Assembled contigs evaluation using QUAST
def quast_evaluation(out_path, contigfile, scaffoldfile, plasmid_contigs_file, plasmid_scaffolds_file, final_l500_contig, final_l500_plasmid_contig):
    print "\n################## Running QUAST on Assemblies. ##################\n"
    if not os.path.isfile(contigfile):
        print "Can't find contigs file in output folder. \nSpades:Please check the availability of contigs.fa file in the given output folder or change the name of contigs file to contigs.fasta\nVelvet:Please check the availability of contigs.fa file in the given output folder or change the name of contigs file to contigs.fa\n"
        #exit()
    if not os.path.isfile(scaffoldfile):
        print "Can't find scaffolds file in output folder. \nSpades:Please check the availability of scaffolds.fa file in the given output folder or change the name of contigs file to contigs.fasta\nVelvet:Please check the availability of contigs.fa file in the given output folder or change the name of contigs file to contigs.fa\n"
        #exit()
    if not os.path.isfile(plasmid_contigs_file):
        print "Can't find contigs file in output folder. \nSpades:Please check the availability of contigs.fa file in the given output folder or change the name of contigs file to contigs.fasta\nVelvet:Please check the availability of contigs.fa file in the given output folder or change the name of contigs file to contigs.fa\n"
        #exit()
    if not os.path.isfile(plasmid_scaffolds_file):
        print "Can't find scaffolds file in output folder. \nSpades:Please check the availability of scaffolds.fa file in the given output folder or change the name of contigs file to contigs.fasta\nVelvet:Please check the availability of contigs.fa file in the given output folder or change the name of contigs file to contigs.fa\n"
        #exit()
    if not os.path.isfile(final_l500_contig):
        print "Can't find scaffolds file in output folder. \nSpades:Please check the availability of scaffolds.fa file in the given output folder or change the name of contigs file to contigs.fasta\nVelvet:Please check the availability of contigs.fa file in the given output folder or change the name of contigs file to contigs.fa\n"
        #exit()
    if not os.path.isfile(final_l500_plasmid_contig):
        print "Can't find scaffolds file in output folder. \nSpades:Please check the availability of scaffolds.fa file in the given output folder or change the name of contigs file to contigs.fasta\nVelvet:Please check the availability of contigs.fa file in the given output folder or change the name of contigs file to contigs.fa\n"
        #exit()

    quast_out = out_path + "/quast"
    cmdstring = ConfigSectionMap("bin_path")['binbase'] + ConfigSectionMap("quast")['quast_bin'] + ConfigSectionMap("quast")['base_cmd'] + " " + contigfile + " " + scaffoldfile + " " + plasmid_contigs_file + " " + plasmid_scaffolds_file + " " + final_l500_contig + " " + final_l500_plasmid_contig + " -o " + quast_out
    print "Running: %s" % cmdstring
    os.system(cmdstring)
    print "QUAST evaluation completed.\n"
    print "QUAST results can be found at: %s" % quast_out
    print "\n################## END: QUAST evaluation. ##################\n"
