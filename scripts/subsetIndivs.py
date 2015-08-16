import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def checkContigs(Seq_record_object, min_contig_ln, n_ratio):
    if len(Seq_record_object.seq) >= min_contig_ln:
        seq_string = str(Seq_record_object.seq)
        Nrel = float(seq_string.count('N'))/float(len(seq_string))
        if Nrel <= n_ratio:
            return True
        else:
            return False
    else:
        return False

def createSubsetContigs(name, min_contig_ln, n_ratio, individuals, bestcontigPath, outDir_subset):
    input_handle = open(individuals, "rU")
    list_with_indivs = []
    for indiv in input_handle:
        new_fn = indiv.rsplit('\n')[0]
        list_with_indivs.append(new_fn)
    input_handle.close()
    contig = os.path.join(bestcontigPath, (name + '.fasta'))
    fasta_file_in = open(contig, "rU")
    fasta_sequences_in = SeqIO.parse(fasta_file_in, 'fasta')
    out_file_name = os.path.join(outDir_subset, (name + '_subset.fasta'))
    fasta_file_out = open(out_file_name, 'wt')
    for record in fasta_sequences_in:
        if (record.id in list_with_indivs) and (checkContigs(record, min_contig_ln, n_ratio) == True):
            SeqIO.write(record, fasta_file_out, "fasta")
        else:
            pass
    fasta_sequences_in.close()
    fasta_file_out.close()
    return (name + '_subset.fasta')

