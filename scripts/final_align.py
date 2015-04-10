#Required modules
from Bio import AlignIO
from Bio.Seq import *
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio import Alphabet
import os
import subprocess
import sys
import csv

def select_final_alignments(filt_alignment_list_path, filt_align_path, individuals_list_path, change_id, id_change_list_path, number_missing_indivs_allowed, datatype, final_align_format, out_folder):
    outDir_alignments = os.path.join(out_folder, "alignments")
    subprocess.call("mkdir '%s'" % (outDir_alignments), shell=True)
    indivs = []
    indiv_number = 0
    list_loci_included = os.path.join(out_folder, ('loci_included_%s_indivs_max_missing.txt') % number_missing_indivs_allowed)
    list_loci_not_included = os.path.join(out_folder, ('loci_not_included_%s_indivs_max_missing.txt') % number_missing_indivs_allowed)	
    file_handle_included = open(list_loci_included, 'w')
    file_handle_not_included = open(list_loci_not_included, 'w')
    loci_included_counter = 0
    try:
        individuals_file = open(individuals_list_path, 'rU')
        alignments_file = open (filt_alignment_list_path, 'rU')
    except IOError:
        sys.exit('Could not open list with individuals and/or alignments!!')
    else:
        for individual in individuals_file:
                name = individual.rsplit('\n', 1)[0]
                indivs.append(name)
                indiv_number += 1
        for alignment in alignments_file:
        	indivs_in_align = 0
        	name = alignment.rsplit('_gaps_analyzed.fasta', 1)[0]
        	alignPath = os.path.join(filt_align_path, (name + '_gaps_analyzed.fasta'))
        	alignment = AlignIO.read(alignPath, 'fasta')
        	alignment.sort()
        	for record in alignment:
        		if record.id in indivs:
        			indivs_in_align += 1
        		else:
        			pass
                file_handle_out = open(alignPath, 'w')
                AlignIO.write(alignment, file_handle_out, 'fasta')
                file_handle_out.close()
        	if datatype == 'ambiguous':
        	    if (indiv_number - indivs_in_align) > number_missing_indivs_allowed:
        	        file_handle_not_included.write(name + "_gaps_analyzed.fasta" + '\n')
        	    else:
                        generate_final_alignment(name, alignPath, indivs, change_id, id_change_list_path, final_align_format, outDir_alignments)
        	        file_handle_included.write(name + "_final_align." + final_align_format + '\n')
        	        loci_included_counter += 1
        	elif datatype == 'haplotype':
        	    if (indiv_number - indivs_in_align) > (number_missing_indivs_allowed * 2):
        	        file_handle_not_included.write(name + "_gaps_analyzed.fasta" + '\n')
        	    else:
                        generate_final_alignment(name, alignPath, indivs, change_id, id_change_list_path, final_align_format, outDir_alignments)
                        file_handle_included.write(name + "_final_align." + final_align_format + '\n')
                        loci_included_counter += 1
                else:
                    sys.exit("datatype not correctly specified, final alignment generation aborted")               		
        individuals_file.close()
        alignments_file.close()
        file_handle_included.close()
        file_handle_not_included.close()
        return loci_included_counter


def findSpeciesID(sequence_id_name, dict_reader_object):
    for line in dict_reader_object:
        if line['old_ID'] == sequence_id_name:
            return line['new_ID']
        else:
            pass


def generate_final_alignment(base_name, alignment_path, indivs_list, change_id, id_change_list_path, final_align_format, out_folder):
    alignment = AlignIO.read(alignment_path, 'fasta')
    alignLength = alignment.get_alignment_length()
    indivs_in_align = []
    msa_temp = []
    for record in alignment:
        indivs_in_align.append(record.id)
    for individual in indivs_list:
        if individual in indivs_in_align:
            for record in alignment:
                if record.id == individual:
                    sequence = str(record.seq)
                    sequence = sequence.replace('-', 'N')
                    if change_id == 1:
                        names = open(id_change_list_path, 'rU')
                        reader = csv.DictReader(names, delimiter = '\t')
                        msa_temp.append(SeqRecord(Seq(sequence, alphabet = Alphabet.generic_dna), id=findSpeciesID(record.id, reader)))
                        names.seek(0)
                        names.close()
                    else:
                        msa_temp.append(SeqRecord(Seq(sequence, alphabet = Alphabet.generic_dna), id=individual))
                else:
                    pass
        else:
            if change_id == 1:
                names = open(id_change_list_path, 'rU')
                reader = csv.DictReader(names, delimiter = '\t')
                missing_individual = SeqRecord(Seq(("N"* alignLength), alphabet = Alphabet.generic_dna), id=findSpeciesID(individual, reader))
                msa_temp.append(missing_individual)
                names.seek(0)
                names.close()
            else:
                missing_individual = SeqRecord(Seq("N"* alignLength), id=individual)
                msa_temp.append(missing_individual)
    final_alignment = MultipleSeqAlignment(msa_temp)       
    file_out = os.path.join(out_folder, (base_name + "_final_align." + final_align_format))
    file_out_handle = open(file_out, "w")
    AlignIO.write(final_alignment, file_out_handle, final_align_format)
    file_out_handle.close()



