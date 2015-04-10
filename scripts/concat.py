from Bio import AlignIO
import subprocess
import os
import sys
from Bio import Alphabet

def generatePartitionInput(concat_align_path, out_file_name, partition, branch, model_evol, model_sel, search): # Partition is in this case a list
    output_path = os.path.join(concat_align_path, "partition_finder.cfg")
    pf_file = open('%s' % (output_path), 'w')
    pf_file.write("alignment = %s;\n" % (out_file_name))
    pf_file.write("\n")
    pf_file.write("branchlengths = %s;\n" % (branch))
    pf_file.write("\n")
    pf_file.write("models = %s;\n" % (model_evol))
    pf_file.write("\n")
    pf_file.write("model_selection =  %s;\n" % (model_sel))
    pf_file.write("[data_blocks]\n")
    for record in partition:
        pf_file.write(record)
    pf_file.write("\n")
    pf_file.write("[schemes]\n")
    pf_file.write("search = %s;\n" % (search))
    pf_file.close()


def concatenateLoci(geneList, filt_align_path, min_align_length, concat_align_path, out_name, final_align_format, concat_align_format, branch, model_evol, model_sel, search):
    try:
        gene_file = open (geneList, 'rU')
    except IOError:
        print "Could not open file with final align names and/or list with individuals"
    else:
        gene_counter = 0
        gene_present = 0
        genes = []
        partition = []
        for gene in gene_file:
            genes.append(gene)
            gene_present = gene_present + 1
        if gene_present != 0:
            for gene in genes:
                if gene_counter == 0:
                    name = gene.rsplit('\n', 1)[0]
                    gene_name = name.rsplit(('_final_align'), 1)[0]
                    alignPath = os.path.join(filt_align_path, name)
                    alignment = AlignIO.read(alignPath, final_align_format)
                    alignment.sort()
                    alignLength = alignment.get_alignment_length()
                    if alignLength > min_align_length:
                        concat_align = alignment
                        partition.append("%s_pos1 = 1-%s\\3;\n" % (gene_name, alignLength))
                        partition.append("%s_pos2 = 2-%s\\3;\n" % (gene_name, alignLength))
                        partition.append("%s_pos3 = 3-%s\\3;\n" % (gene_name, alignLength))
                        gene_counter = gene_counter + 1
                    else:
                        pass
                else:
                    name = gene.rsplit('\n', 1)[0]
                    gene_name = name.rsplit('_final_align', 1)[0]
                    alignPath = os.path.join(filt_align_path, name)
                    alignment = AlignIO.read(alignPath, final_align_format)
                    alignment.sort()
                    alignLength = alignment.get_alignment_length()
                    if alignLength > min_align_length:
                        concatLength_old = concat_align.get_alignment_length()
                        concat_align = concat_align + alignment
                        concatLength_new = concat_align.get_alignment_length()
                        partition.append("%s_pos1 = %s-%s\\3;\n" % (gene_name, (concatLength_old + 1) ,concatLength_new))
                        partition.append("%s_pos2 = %s-%s\\3;\n" % (gene_name, (concatLength_old + 2) ,concatLength_new))
                        partition.append("%s_pos3 = %s-%s\\3;\n" % (gene_name, (concatLength_old + 3) ,concatLength_new))
                    else:
                        pass
            out_file_name = out_name + "_concat_loci.phy"
            file_handle_phylip = os.path.join(concat_align_path, out_file_name)
            phylip_file = open(file_handle_phylip, "w")
            AlignIO.write(concat_align, phylip_file, 'phylip-sequential')
            phylip_file.close()
            phylip_alignment = open(file_handle_phylip, 'rU')
            align_old = AlignIO.parse(phylip_alignment, 'phylip-sequential', alphabet = Alphabet.generic_dna)
            file_handle_out = os.path.join(concat_align_path, (out_name + "_concat_loci." + concat_align_format))
            new_alignment = open(file_handle_out, 'w')
            AlignIO.write(align_old, new_alignment, concat_align_format)
            phylip_alignment.close()
            new_alignment.close()
            gene_file.close()
            generatePartitionInput(concat_align_path, out_file_name, partition, branch, model_evol, model_sel, search)
        else:
            file_handle_out = os.path.join(concat_align_path, (out_name + "_error_report.txt"))
            summary_file = open(file_handle_out, "w")
            summary_file.write('No genes were present in the final alignment file (perhaps due to stringency requirement missing data allowed??)')
            summary_file.close()

