import os
from Bio import SeqIO

def createSubsetContigs(name, individuals, bestcontigPath, outDir_subset):
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
        if record.id in list_with_indivs:
            SeqIO.write(record, fasta_file_out, "fasta")
        else:
            pass
    fasta_sequences_in.close()
    fasta_file_out.close()
    return (name + '_subset.fasta')
