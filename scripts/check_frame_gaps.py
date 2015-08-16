from Bio import AlignIO
from Bio.Seq import *
from Bio.Align import AlignInfo
from Bio.Alphabet import *
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
import os

def aln_len_check(alignment_path, indiv_gap_window, column_gap_window):
    file_handle = open(alignment_path, 'r')
    temp = AlignIO.read(file_handle, 'fasta')
    alignment_length = temp.get_alignment_length()
    aa_len = float(alignment_length)/float(3)
    file_handle.close()
    if (aa_len >= indiv_gap_window) and (aa_len >= column_gap_window):
        return True
    else:
        return False


def removeGaps(alignment, alignment_path, outDir_align_check_temp, indiv_gap_window, indiv_gap_ratio, column_gap_window, column_gap_ratio, max_insertion_cut_off):
    name = alignment.rsplit('_basic_align.fasta', 1)[0]
    file_in_handle = open(alignment_path, 'rU')
    muscle_align = AlignIO.read(file_in_handle, 'fasta')
    msa_temp = []
    align_counter = 0
    for record in muscle_align:
        align_counter += 1
    if align_counter > 1:
        for record in muscle_align:
            sequence = str(record.seq)
            sequence_length = len(sequence)
            sliding_window = (((indiv_gap_window / 2) + 1) * 3) #This is the number of codons that you shift every jump
            replacement = "-" * indiv_gap_window * 3
            counter = 0
            length_sequence_checked = 0
            while (length_sequence_checked + sliding_window) < sequence_length:
                gap_codon_window = sequence[(counter * sliding_window) : ((counter * sliding_window) + (indiv_gap_window * 3))]
                length_sequence_checked = (counter * sliding_window) + (indiv_gap_window * 3)
                ratio_gap = float(gap_codon_window.count('-'))/float(len(gap_codon_window))
                if ratio_gap > indiv_gap_ratio:
                    sequence = sequence.replace(gap_codon_window, replacement)
                    counter += 1
                else:
                    counter += 1
            counter += 1
            sequence_length = len(sequence)
            seq_stretch_remain = sequence_length - length_sequence_checked
            sequence = sequence[0:-seq_stretch_remain]      # Remove remaining codons (less than 1 x sliding window)
            msa_temp.append(SeqRecord(Seq(sequence, generic_alphabet), id=record.id))
        temp_alignment = MultipleSeqAlignment(msa_temp)
        temp_alignment.sort
        counter = 0
        sequence_length = temp_alignment.get_alignment_length()
        length_sequence_checked = 0
        sliding_window = (((column_gap_window / 2) + 1) * 3)
        append_counter = 0
        msa_temp = []
        in_frame_align = MultipleSeqAlignment(msa_temp)
        while length_sequence_checked < sequence_length:
            gap_codon_window = temp_alignment[:, (counter * sliding_window) : ((counter * sliding_window) + (column_gap_window * 3))]
            gapped_column_count = 0
            for i in range(column_gap_window):
                codon = gap_codon_window[:, (i * 3) : ((i * 3) + 3)]
                position = []
                for record in codon:
                    sequence = str(record.seq)
                    sequence = sequence.replace("-", "N")
                    convert = translate(sequence)
                    position.append(convert)
                ratio_gap = float(position.count('X'))/float(len(position))
                if ratio_gap > column_gap_ratio:
                    gapped_column_count += 1
                else:
                    pass
            if column_gap_window == 1:
                if gapped_column_count != 1:
                    codons_to_append = temp_alignment[:, (counter * sliding_window) : ((counter * sliding_window) + (sliding_window))]
                    if append_counter == 0:
                        in_frame_align = codons_to_append
                        append_counter += 1
                    else:
                        in_frame_align = in_frame_align + codons_to_append
                        append_counter += 1
                else:
                    pass
            else:
                if gapped_column_count < ((column_gap_window / 2) + 1):
                    codons_to_append = temp_alignment[:, (counter * sliding_window) : ((counter * sliding_window) + (sliding_window))]
                    if append_counter == 0:
                        in_frame_align = codons_to_append
                        append_counter += 1
                    else:
                        in_frame_align = in_frame_align + codons_to_append
                        append_counter += 1
                else:
                    pass
            length_sequence_checked = (counter * sliding_window) + (column_gap_window * 3)
            counter += 1
        alignment_length = in_frame_align.get_alignment_length()
        length_alignment_checked = 0
        append_counter = 0
        msa_temp_2 = []
        final_align = MultipleSeqAlignment(msa_temp_2)
        while length_alignment_checked < alignment_length:
            base_column = in_frame_align[:, length_alignment_checked]
            if ((len(base_column)) - base_column.count('-')) > max_insertion_cut_off:
                if append_counter == 0:
                    final_align = in_frame_align[:, length_alignment_checked:(length_alignment_checked+1)]
                    append_counter += 1
                else:
                    final_align = final_align + in_frame_align[:, length_alignment_checked:(length_alignment_checked+1)]
                    append_counter += 1
            else:
                pass
            length_alignment_checked += 1
        final_align_length = final_align.get_alignment_length()
        final_align = final_align[:, 0:final_align_length]
        out_path_file = os.path.join(outDir_align_check_temp, (name + "_gaps_analyzed.fasta"))
        file_out_handle = open(out_path_file, "w")
        AlignIO.write(final_align, file_out_handle, "fasta")
        file_in_handle.close()
        file_out_handle.close()
        return final_align_length
    elif align_counter == 1:
        for record in muscle_align:
            sequence = str(record.seq)
            sequence_length = len(sequence)
            sliding_window = (((indiv_gap_window / 2) + 1) * 3)
            replacement = "-" * indiv_gap_window * 3
            counter = 0
            length_sequence_checked = 0
            while (length_sequence_checked + sliding_window) < sequence_length:
                gap_codon_window = sequence[(counter * sliding_window) : ((counter * sliding_window) + (indiv_gap_window * 3))]
                length_sequence_checked = (counter * sliding_window) + (indiv_gap_window * 3)
                ratio_gap = float(gap_codon_window.count('-'))/float(len(gap_codon_window))
                if ratio_gap > indiv_gap_ratio:
                    sequence = sequence.replace(gap_codon_window, replacement)
                    counter += 1
                else:
                    counter += 1
            counter -= 1
            sequence_length = len(sequence)
            seq_stretch_remain = sequence_length - length_sequence_checked
            sequence = sequence[0:-seq_stretch_remain]      # Remove remaining codons (less than 1 x sliding window)
            msa_temp.append(SeqRecord(Seq(sequence, generic_alphabet), id=record.id))
            temp_alignment = MultipleSeqAlignment(msa_temp)
        temp_alignment_length = temp_alignment.get_alignment_length()
        temp_alignment = temp_alignment[:, 0:temp_alignment_length]
        out_path_file = os.path.join(outDir_align_check_temp, (name + "_gaps_analyzed.fasta"))
        file_out_handle = open(out_path_file, "w")
        AlignIO.write(temp_alignment, file_out_handle, "fasta")
        file_in_handle.close()
        file_out_handle.close()
        return temp_alignment_length
    else:
        return 0



def identifyDivAlign(alignment_path, sliding_window, max_heterozygosity_ratio):
    file_in_handle = open(alignment_path, 'r')
    nt_alignment = AlignIO.read(file_in_handle, 'fasta')
    msa_temp = []
    for record in nt_alignment:
        sequence = str(record.seq)
        nt_sequence = sequence.replace("-", "N")
        aa_sequence = translate(nt_sequence)
        msa_temp.append(SeqRecord(Seq(aa_sequence, generic_alphabet), id=record.id))
    aa_alignment = MultipleSeqAlignment(msa_temp)
    summary_align = AlignInfo.SummaryInfo(aa_alignment)
    aa_alignment_length = aa_alignment.get_alignment_length()
    concensus_seq = str(summary_align.dumb_consensus(threshold=0.5))
    for indiv in aa_alignment:
        msa_temp = []
        record_seq = str(indiv.seq)
        msa_temp.append(SeqRecord(Seq(record_seq), id=indiv.id))
        msa_temp.append(SeqRecord(Seq(concensus_seq), id="consensus"))
        temp_alignment = MultipleSeqAlignment(msa_temp)
        alignment_length = temp_alignment.get_alignment_length()
        length_alignment_checked = 0
        sliding_window_counter = 0
        while length_alignment_checked < alignment_length:
            window = temp_alignment[:, sliding_window_counter:(sliding_window_counter + sliding_window)]
            heterozygous_counter = 0
            for i in range(sliding_window):
                codon_column = window[:, i]
                aa = list(set(codon_column))
                if len(aa) < 2:
                    pass
                else:
                    if 'X' in aa:
                        pass
                    else:
                        heterozygous_counter = heterozygous_counter + 1
            heterozygous_ratio = float(heterozygous_counter)/float(sliding_window)
            if heterozygous_ratio > max_heterozygosity_ratio:
                return False
            else:
                pass
            sliding_window_counter += 1
            length_alignment_checked = (sliding_window_counter + sliding_window)
    file_in_handle.close()
    return True

def identifyStopCodons(alignment_path, max_stop):
    file_handle = open(alignment_path, 'r')
    temp = AlignIO.read(file_handle, 'fasta')
    alignment_length = temp.get_alignment_length()
    record_count = 0
    first_frame_count = 0
    second_frame_count = 0
    third_frame_count = 0
    for record in temp:
        record_count += 1
        sequence = str(record.seq)
        sequence_length = len(sequence)
        sequence = sequence.replace("-", "N")
        codon_n = 0
        stretch_check = 0
        stop_count_first_frame = 0
        while stretch_check < sequence_length:                                                                      # Count the number of stop codons for this given sequence, in all three reading frames
            codon = sequence[((0 + codon_n)*3):((codon_n + 1)*3)]
            if (codon == 'TGA') or (codon == 'TAA') or (codon == 'TAG'):
                stop_count_first_frame += 1
                codon_n += 1
            else:
                codon_n += 1
            stretch_check = (codon_n + 1) * 3
        sequence_second_frame = sequence[1:(sequence_length - 2)]
        sequence_length_second_frame = len(sequence_second_frame)
        stop_count_second_frame = 0
        codon_n = 0
        stretch_check = 0
        while stretch_check < sequence_length_second_frame:
            codon = sequence_second_frame[((0 + codon_n)*3):((codon_n + 1)*3)]
            if (codon == 'TGA') or (codon == 'TAA') or (codon == 'TAG'):
                stop_count_second_frame += 1
                codon_n += 1
            else:
                codon_n += 1
            stretch_check = (codon_n + 1) * 3
        sequence_third_frame = sequence[2:(sequence_length - 1)]
        sequence_length_third_frame = len(sequence_third_frame)
        stop_count_third_frame = 0
        codon_n = 0
        stretch_check = 0
        while stretch_check < sequence_length_third_frame:
            codon = sequence_third_frame[((0 + codon_n)*3):((codon_n + 1)*3)]
            if (codon == 'TGA') or (codon == 'TAA') or (codon == 'TAG'):
                stop_count_third_frame += 1
                codon_n += 1
            else:
                codon_n += 1
            stretch_check = (codon_n + 1) * 3
        if (stop_count_first_frame <= stop_count_second_frame) and (stop_count_first_frame <= stop_count_third_frame):  # Check whether the number of stop codon for this individual is higher in firs, second or third frame
            if stop_count_first_frame > max_stop:
                return False
            else:
                if stop_count_first_frame == stop_count_second_frame:
                    first_frame_count += 1
                    second_frame_count += 1
                elif stop_count_first_frame == stop_count_third_frame:
                    first_frame_count += 1
                    third_frame_count += 1
                else:
                    first_frame_count += 1
        elif (stop_count_second_frame <= stop_count_third_frame):
            if stop_count_second_frame > max_stop:
                return False
            else:
                if stop_count_second_frame == stop_count_third_frame:
                    second_frame_count += 1
                    third_frame_count += 1
                else:
                    second_frame_count += 1
        else: 
            if stop_count_third_frame > max_stop:
                return False
            else:
                third_frame_count += 1
    if (first_frame_count >= second_frame_count) and (first_frame_count >= third_frame_count):  # Check whether the number of stop codons for each individuals, across all individuas, was higher in first, second, or third frame
        return True
    elif second_frame_count >= third_frame_count:
        file_handle = open(alignment_path, 'r')
        alignment_old = AlignIO.read(file_handle, 'fasta')
        file_handle.close()
        alignment_new = alignment_old[:, 1:(alignment_length - 2)]
        file_handle_out = open(alignment_path, 'w')
        AlignIO.write(alignment_new, file_handle_out, "fasta")
        file_handle_out.close()
    else:
        file_handle = open(alignment_path, 'r')
        alignment_old = AlignIO.read(file_handle, 'fasta')
        file_handle.close()
        alignment_new = alignment_old[:, 2:(alignment_length - 1)]
        file_handle_out = open(alignment_path, 'w')
        AlignIO.write(alignment_new, file_handle_out, "fasta")
        file_handle_out.close()
    return True


