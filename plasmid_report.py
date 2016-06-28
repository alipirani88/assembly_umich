__author__ = 'alipirani'
import os
import argparse
import re
import subprocess
import statistics
from collections import defaultdict
from collections import OrderedDict
parser = argparse.ArgumentParser(description='Generate Plasmid Assembly Report')
parser.add_argument('-filename', action='store', dest="filename", help='filename of assembly in fasta format')
args = parser.parse_args()


sample_name = args.filename
command = "/home/apirani//bin/bioawk-master/bioawk -c fastx '{ print $name, length($seq) }' < %s" % sample_name
output = subprocess.check_output(command, shell=True)
#print output

uniq_plasmid_names = []
for line in output.splitlines():
    linesplit = line.split('\t')
    plasmid_name = linesplit[0]
    plasmid_name_split = plasmid_name.split('_')
    plasmid_name_uniq = str(plasmid_name_split[6]) + "_" + str(plasmid_name_split[7])
    if plasmid_name_uniq not in uniq_plasmid_names:
        uniq_plasmid_names.append(plasmid_name_uniq)
plasmid_number = len(uniq_plasmid_names)
command_2 = "grep \'>\' %s | wc -l" % sample_name
no_of_contigs = subprocess.check_output(command_2, shell=True)

plasmid_length = OrderedDict()
median_coverage = OrderedDict()
longest_contig = defaultdict(list)
for uniq_plasmid in uniq_plasmid_names:
    uniq_plasmid_length = 0
    uniq_plasmid_coverage = 0
    median_array = []
    for line in output.splitlines():
        linesplit = line.split('\t')
        plasmid_name = linesplit[0]
        plasmid_name_split = plasmid_name.split('_')
        plasmid_name_uniq = str(plasmid_name_split[6]) + "_" + str(plasmid_name_split[7])

        if uniq_plasmid == plasmid_name_uniq:
            uniq_plasmid_length = uniq_plasmid_length + int(linesplit[1])
            uniq_plasmid_coverage = uniq_plasmid_coverage + float(plasmid_name_split[5])
            median_array.append(uniq_plasmid_coverage)
            longest_contig[uniq_plasmid].append(linesplit[1])

    median_coverage[uniq_plasmid] = statistics.median(median_array)
    #print statistics.median(median_array)
    command_3 = "grep \'%s\' %s | wc -l" % (uniq_plasmid, sample_name)
    no_of_uniq_plasmid_contigs = subprocess.check_output(command_3, shell=True)
    plasmid_length[uniq_plasmid] = uniq_plasmid_length
    #print uniq_plasmid + "\t" + str(uniq_plasmid_length)
    mean_coverage = float(uniq_plasmid_coverage) / int(no_of_uniq_plasmid_contigs)
    #print uniq_plasmid + "\t" + str(uniq_plasmid_coverage) + "\t" + str(statistics.median(median_array)) + "\t" + str(mean_coverage)

#print plasmid_length


plasm_comp_size = sum(plasmid_length.values())
#print plasm_comp_size
#print longest_contig
#for i in longest_contig:
    #maximum = max(longest_contig[i], key=longest_contig.get)
    #print i + "\t" + str(longest_contig[maximum])
    #print i + "\t" + longest_contig[i][0]



command_4 = "/home/apirani/bin/ncbi-blast-2.2.30+/bin/blastn -db /nfs/esnitkin/Ali/Project_Sample_2008_KPC_analysis//KPC2.fasta -query %s -outfmt 6" % sample_name
blast_output = subprocess.check_output(command_4, shell=True)

kpc_positive = OrderedDict()

for line in blast_output.splitlines():
    linesplit = line.split('\t')
    query_id = linesplit[0]
    align_length = linesplit[3]
    if align_length > 800:
        query_id_split = query_id.split('_')
        query_id_split_uniq = str(query_id_split[6]) + "_" + str(query_id_split[7])
        kpc_positive[query_id_split_uniq] = "KPC+"
#print kpc_positive



# kpc_string = ""
# for key in kpc_positive:
#     kpc_string = kpc_string + ",,,,,," + key + ":\t" + kpc_positive[key] + "\n"
#

#print kpc_string







final_string = ""
no_of_plasmids = len(uniq_plasmid_names)
total_no_of_contigs = no_of_contigs.strip()
plasmid_name_length_string = ""
for i in plasmid_length:
    plasmid_name_length_string = plasmid_name_length_string + i + ":" + str(plasmid_length[i]) + ","
plasmid_median_string = ""
# for key in median_coverage:
#     plasmid_median_string = '\n'.join(str(median_coverage[key]))
for i in median_coverage:
    plasmid_median_string = plasmid_median_string + i + ":" + str(median_coverage[i]) + ","
longest_contig_string = ""
for i in longest_contig:
    longest_contig_string = longest_contig_string + i + ":" + str(longest_contig[i][0]) + ","

length_median_longest_string = "\t\n"
for i in plasmid_length:
    if i in kpc_positive.keys():
        kpc_string = "KPC+"
    else:
        kpc_string = ""
    length_median_longest_string = length_median_longest_string + ",,,,," + i + "," + str(plasmid_length[i]) + "," + str(median_coverage[i]) + "," + str(longest_contig[i][0]) + "," + kpc_string + "\n"













header = "sample\tTotal # of contigs\t# of plasmids\tplasmid_length(bp)\tplasmid_median_cov\ttotal_contigs_length(bp)\tlongest_contigs_in_plasmid"
print header
final_string = os.path.basename(sample_name) + "\t" + str(total_no_of_contigs) + "\t" + str(no_of_plasmids) + "\t" + plasmid_name_length_string + "\t" + plasmid_median_string + "\t" + str(plasm_comp_size) + "\t" + longest_contig_string
#print final_string


print sample_name + "\t" + str(total_no_of_contigs) + "\t" + str(no_of_plasmids) + "\t" + str(plasm_comp_size) + "\t\t\t" + length_median_longest_string























