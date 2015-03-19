from Bio import AlignIO
import subprocess
import os
import sys

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


def concatenateLoci(geneList, indivList, filt_align_path, min_align_length, concat_align_path, out_name,  branch, model_evol, model_sel, search):
    try:
        gene_file = open (geneList, 'rU')
        id_file = open (indivList, 'rU')
    except IOError:
        print "Could not open file with final align names and/or list with individuals"
    else:
        gene_counter = 0
        gene_present = 0
        genes = []
        indivs = []
        indivs_counter = 0
        partition = []
        for gene in gene_file:
            genes.append(gene)
            gene_present = gene_present + 1
        if gene_present != 0:
            for individual in id_file:
                name = individual.rsplit('\n', 1)[0]
                indivs.append(name)
                indivs_counter = indivs_counter + 1
            for gene in genes:
                if gene_counter == 0:
                    name = gene.rsplit('\n', 1)[0]
                    gene_name = name.rsplit('_final_align.fasta', 1)[0]
                    alignPath = os.path.join(filt_align_path, name)
                    alignment = AlignIO.read(alignPath, 'fasta')
                    alignment.sort()
                    indivPresence = 0
                    alignLength = alignment.get_alignment_length()
                    if alignLength > min_align_length:
                        for record in alignment:
                            if record.id in indivs:
                                indivPresence = indivPresence + 1
                                continue
                            else:
                                break
                        if indivPresence == indivs_counter:
                            concat_align = alignment
                            partition.append("%s_pos1 = 1-%s\\3;\n" % (gene_name, alignLength))
                            partition.append("%s_pos2 = 2-%s\\3;\n" % (gene_name, alignLength))
                            partition.append("%s_pos3 = 3-%s\\3;\n" % (gene_name, alignLength))
                            gene_counter = gene_counter + 1
                        else:
                            print "Not all designated individuals present in alignment '%s', therefore gene not appended" % (name)
                    else:
                        pass
                else:
                    name = gene.rsplit('\n', 1)[0]
                    gene_name = name.rsplit('_final_align.fasta', 1)[0]
                    alignPath = os.path.join(filt_align_path, name)
                    alignment = AlignIO.read(alignPath, 'fasta')
                    alignment.sort()
                    indivPresence = 0
                    alignLength = alignment.get_alignment_length()
                    if alignLength > min_align_length:
                        for record in alignment:
                            if record.id in indivs:
                                indivPresence = indivPresence + 1
                                continue
                            else:
                                break
                        if indivPresence == indivs_counter:
                            concatLength_old = concat_align.get_alignment_length()
                            concat_align = concat_align + alignment
                            concatLength_new = concat_align.get_alignment_length()
                            partition.append("%s_pos1 = %s-%s\\3;\n" % (gene_name, (concatLength_old + 1) ,concatLength_new))
                            partition.append("%s_pos2 = %s-%s\\3;\n" % (gene_name, (concatLength_old + 2) ,concatLength_new))
                            partition.append("%s_pos3 = %s-%s\\3;\n" % (gene_name, (concatLength_old + 3) ,concatLength_new))
                        else:
                            print "Not all designated individuals present in alignment '%s', therefore not appended" % (name)
                    else:
                        pass
            for record in concat_align:
                record.id = record.id.split("ndexing", 1)[0] + record.id.split("ndexing", 1)[1]
            file_handle_fasta = os.path.join(concat_align_path, (out_name + "_concat_loci.fasta"))
            fasta_file = open(file_handle_fasta, "w")
            align_fasta = AlignIO.write(concat_align, fasta_file, 'fasta')
            fasta_file.close()
            out_file_name = out_name + "_concat_loci.phy"
            file_handle_phylip = os.path.join(concat_align_path, out_file_name)
            fasta_file = open(file_handle_fasta, "r")
            phylip_file = open(file_handle_phylip, "w")
            align_phylip = AlignIO.convert(fasta_file, 'fasta', phylip_file, 'phylip-sequential')
            phylip_file.close()
            fasta_file.close()
            gene_file.close()
            id_file.close()
            generatePartitionInput(concat_align_path, out_file_name, partition, branch, model_evol, model_sel, search)
        else:
            file_handle_out = os.path.join(concat_align_path, (out_name + "_error_report.txt"))
            summary_file = open(file_handle_out, "w")
            summary_file.write('No genes were present in the final alignment file (perhaps due to stringency requirement missing data allowed??)')
            summary_file.close()

