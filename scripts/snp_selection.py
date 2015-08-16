import os
from Bio import AlignIO
from Bio.Seq import *
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import *


def findSNP_leastN(alignment):
    alignment.sort()
    indiv_count = 0
    msa_temp = []
    aln_len = alignment.get_alignment_length()
    for record in alignment:
        msa_temp.append(SeqRecord(Seq('', generic_alphabet), id=record.id))
        indiv_count += 1
    final_snp_align = MultipleSeqAlignment(msa_temp)
    aln_len_counter = 0
    N_min = indiv_count
    while aln_len_counter < aln_len:
        site = alignment[:, aln_len_counter]
        N = site.count('N')
        aln_len_counter += 1
        if N < N_min:
            N_min = N
        else:
            pass
    aln_len_counter = 0
    while aln_len_counter < aln_len:
        site = alignment[:, aln_len_counter]
        N = site.count('N')
        if N == N_min:
            final_snp_align = final_snp_align + alignment[:, (aln_len_counter):(aln_len_counter+1)]
            aln_len_counter += 1
        else:
            aln_len_counter += 1
    return final_snp_align


def findSite_Biallelic(alignment_path, alignment_format, N_ratio):
    align = AlignIO.read(alignment_path, alignment_format)
    align.sort()
    align_length = align.get_alignment_length()
    align_length_counter = 0
    indiv_number = 0
    msa_temp = []
    for record in align:
        msa_temp.append(SeqRecord(Seq('', generic_alphabet), id=record.id))
        indiv_number += 1
    final_snp_align = MultipleSeqAlignment(msa_temp)
    while align_length_counter < align_length:
        site = align[:, align_length_counter]
        N = site.count('N')
        if ((float(N)/float(indiv_number)) <= N_ratio):
            snp = findSNP_Biallelic(site, indiv_number)
            if snp == True:
                final_snp_align = final_snp_align + align[:, (align_length_counter):(align_length_counter+1)]
                align_length_counter += 1
            else:
                align_length_counter += 1
        else:
            align_length_counter += 1
    return final_snp_align          


def findSNP_Biallelic(site, indiv_number):
    snp = []
    a = site.count('A')
    t = site.count('T')
    c = site.count('C')
    g = site.count('G')
    r = site.count('R')
    y = site.count('Y')
    s = site.count('S')
    w = site.count('W')
    k = site.count('K')
    m = site.count('M')
    b = site.count('B')
    d = site.count('D')
    h = site.count('H')
    v = site.count('V')
    n = site.count('N')
    if (a != (indiv_number - n)) and (a != 0):
        snp.append('yes')
    elif (t != (indiv_number - n)) and (t != 0):
        snp.append('yes')
    elif (c != (indiv_number - n)) and (c != 0):
        snp.append('yes')
    elif (g != (indiv_number - n)) and (g != 0):
        snp.append('yes')
    elif (r != (indiv_number - n)) and (r != 0):
        snp.append('yes')
    elif (y != (indiv_number - n)) and (y != 0):
        snp.append('yes')
    elif (s != (indiv_number - n)) and (s != 0):
        snp.append('yes')
    elif (w != (indiv_number - n)) and (w != 0):
        snp.append('yes')
    elif (k != (indiv_number - n)) and (k != 0):
        snp.append('yes')
    elif (m != (indiv_number - n)) and (m != 0):
        snp.append('yes')
    elif (b != (indiv_number - n)) and (b != 0):
        snp.append('yes') 
    elif (d != (indiv_number - n)) and (d != 0):
        snp.append('yes')
    elif (h != (indiv_number - n)) and (h != 0):
        snp.append('yes')
    elif (v != (indiv_number - n)) and (v != 0):
        snp.append('yes')
    else: 
        snp.append('no')
    if (snp.count('yes') == 1) and correctHeterozygote_Biallelic(site):
        return True
    else:
        return False

def correctHeterozygote_Biallelic(input_string):
    no_missing_data_string = input_string.replace("N", "")
    genotype_calls = list(set(no_missing_data_string))
    if len(genotype_calls) < 4:
        if len(genotype_calls) == 3:
            if 'R' in genotype_calls:
                if ('A' in genotype_calls) and ('G' in genotype_calls):
                    return True
                else:
                    return False
            elif 'Y' in genotype_calls:
                if ('C' in genotype_calls) and ('T' in genotype_calls):
                    return True
                else:
                    return False
            elif 'S' in genotype_calls:
                if ('G' in genotype_calls) and ('C' in genotype_calls):
                    return True
                else:
                    return False
            elif 'W' in genotype_calls:
                if ('A' in genotype_calls) and ('T' in genotype_calls):
                    return True
                else:
                    return False
            elif 'K' in genotype_calls:
                if ('G' in genotype_calls) and ('T' in genotype_calls):
                    return True
                else:
                    return False
            elif 'M' in genotype_calls:
                if ('A' in genotype_calls) and ('C' in genotype_calls):
                    return True
                else:
                    return False
            elif 'B' in genotype_calls:
                return False
            elif 'D' in genotype_calls:
                return False               
            elif 'H' in genotype_calls:
                return False               
            elif 'V' in genotype_calls:
                return False
            else:
                return False    
        else:
            pass
        if len(genotype_calls) == 2:
            if 'R' in genotype_calls:
                if ('A' in genotype_calls) or ('G' in genotype_calls):
                    return True
                else:
                    return False
            elif 'Y' in genotype_calls:
                if ('C' in genotype_calls) or ('T' in genotype_calls):
                    return True
                else:
                    return False
            elif 'S' in genotype_calls:
                if ('G' in genotype_calls) or ('C' in genotype_calls):
                    return True
                else:
                    return False
            elif 'W' in genotype_calls:
                if ('A' in genotype_calls) or ('T' in genotype_calls):
                    return True
                else:
                    return False
            elif 'K' in genotype_calls:
                if ('G' in genotype_calls) or ('T' in genotype_calls):
                    return True
                else:
                    return False
            elif 'M' in genotype_calls:
                if ('A' in genotype_calls) or ('C' in genotype_calls):
                    return True
                else:
                    return False
            elif 'B' in genotype_calls:
                return False
            elif 'D' in genotype_calls:
                return False               
            elif 'H' in genotype_calls:
                return False               
            elif 'V' in genotype_calls:
                return False
            else:
                return True
        else:
            pass
    else:
        return False

def findSite(alignment_path, alignment_format, N_ratio):
    align = AlignIO.read(alignment_path, alignment_format)
    align.sort()
    align_length = align.get_alignment_length()
    align_length_counter = 0
    indiv_number = 0
    msa_temp = []
    for record in align:
        msa_temp.append(SeqRecord(Seq('', generic_alphabet), id=record.id))
        indiv_number += 1
    final_snp_align = MultipleSeqAlignment(msa_temp)
    while align_length_counter < align_length:
        site = align[:, align_length_counter]
        N = site.count('N')
        if ((float(N)/float(indiv_number)) <= N_ratio):
            snp = findSNP(site, indiv_number)
            if snp == True:
                final_snp_align = final_snp_align + align[:, (align_length_counter):(align_length_counter+1)]
                align_length_counter += 1
            else:
                align_length_counter += 1
        else:
            align_length_counter += 1
    return final_snp_align          
    
def findSNP(site, indiv_number):
    snp = []
    a = site.count('A')
    t = site.count('T')
    c = site.count('C')
    g = site.count('G')
    r = site.count('R')
    y = site.count('Y')
    s = site.count('S')
    w = site.count('W')
    k = site.count('K')
    m = site.count('M')
    b = site.count('B')
    d = site.count('D')
    h = site.count('H')
    v = site.count('V')
    n = site.count('N')
    if (a != (indiv_number - n)) and (a != 0):
        snp.append('yes')
    elif (t != (indiv_number - n)) and (t != 0):
        snp.append('yes')
    elif (c != (indiv_number - n)) and (c != 0):
        snp.append('yes')
    elif (g != (indiv_number - n)) and (g != 0):
        snp.append('yes')
    elif (r != (indiv_number - n)) and (r != 0):
        snp.append('yes')
    elif (y != (indiv_number - n)) and (y != 0):
        snp.append('yes')
    elif (s != (indiv_number - n)) and (s != 0):
        snp.append('yes')
    elif (w != (indiv_number - n)) and (w != 0):
        snp.append('yes')
    elif (k != (indiv_number - n)) and (k != 0):
        snp.append('yes')
    elif (m != (indiv_number - n)) and (m != 0):
        snp.append('yes')
    elif (b != (indiv_number - n)) and (b != 0):
        snp.append('yes') 
    elif (d != (indiv_number - n)) and (d != 0):
        snp.append('yes')
    elif (h != (indiv_number - n)) and (h != 0):
        snp.append('yes')
    elif (v != (indiv_number - n)) and (v != 0):
        snp.append('yes')
    else: 
        snp.append('no')
    if (snp.count('yes') == 1):
        return True
    else:
        return False
