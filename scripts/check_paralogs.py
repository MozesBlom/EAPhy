from Bio import AlignIO
from Bio.Seq import *
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import *
import matplotlib
import matplotlib.pyplot as plt
from pylab import *
import os

def findHetSites(align_name, align_path):
    alignment_path = os.path.join(align_path, align_name)
    file_handle = open(alignment_path, 'r')
    alignment = AlignIO.read(file_handle, 'fasta')
    alignment.sort()
    indiv_number = 0
    for record in alignment:
    	indiv_number += 1
    alignment_length = alignment.get_alignment_length()
    alignment_length_counter = 0
    heterozygous_sites_counter = 0
    proportion_indiv_heterozygous_at_site = []
    while alignment_length_counter < alignment_length:
    	site = alignment[:, alignment_length_counter]
    	if (not '-' in site) and (len(site) == indiv_number):
    		a = site.count('A')
    		t = site.count('T')
    		c = site.count('C')
    		g = site.count('G')
    		homozygous_sites = a + t + c + g
    		if (homozygous_sites != indiv_number):
    			heterozygous_sites_counter += 1
    			proportion_indiv_heterozygous_at_site.append(1 - (float(homozygous_sites)/float(indiv_number)))
    			alignment_length_counter += 1
    		else:
    			alignment_length_counter += 1
    	else:
    		alignment_length_counter += 1
    file_handle.close()
    proportion_sites_heterozygous = float(heterozygous_sites_counter) / float(alignment_length)
    if (heterozygous_sites_counter != 0):
    	average_proportion_individuals_heterozygous_at_site = float(sum(proportion_indiv_heterozygous_at_site)) / float(heterozygous_sites_counter)
    else:
    	average_proportion_individuals_heterozygous_at_site = float(0)
    return proportion_sites_heterozygous, average_proportion_individuals_heterozygous_at_site

def calcAverageIndivHet(alignment_path):
    file_handle = open(alignment_path, 'r')
    alignment = AlignIO.read(file_handle, 'fasta')
    alignment.sort()
    het_count = []
    indiv_number = 0
    for record in alignment:
    	indiv_number += 1
    for record in alignment:
    	sequence = record.seq
    	seq_length = len(sequence)
    	a = sequence.count('A')
    	t = sequence.count('T')
    	c = sequence.count('C')
    	g = sequence.count('G')
    	unknown = sequence.count('-')
    	het_sites = seq_length - (a + t + c + g + unknown)
    	prop_het_sites = float(het_sites)/float(seq_length)
    	het_count.append(prop_het_sites)
    average_indiv_heterozygosity = float(sum(het_count))/float(indiv_number)
    return average_indiv_heterozygosity

def plotAverageHet(list_het_values, cut_off_value, outDir_align_check_paralog):
    figure(1)
    plt.hist(list_het_values, bins = 20, color = 'c')
    plt.title("Paralog filtering based on heterozygosity level, vertical blue line current cut-off")
    plt.xlabel("Average individivual heterozygosity per alignment")
    plt.ylabel("Alignment number")
    plt.axvline(cut_off_value, color='b', linestyle='dashed', linewidth=2)
    hist_out = os.path.join(outDir_align_check_paralog, 'findParalogsHist.png')
    plt.savefig(hist_out, dpi=128)
    return 'Plot with average heterozygosity per individual, for each alignment saved'



