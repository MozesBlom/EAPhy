#Copyright (C) 2015 Mozes Blom
#
#This collection of scripts, EAPhy, is free software: you can redistribute it and/or modify it
#under the terms of the GNU General Public License as published by the
#Free Software Foundation, either version 3 of the License, or (at your
#option) any later version.
#
#This collection of scripts, EAPhy, is distributed in the hope that it will be useful, but
#WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#General Public License for more details see <http://www.gnu.org/licenses/>.
#
#This workflow depends on several modules of the BioPython package, many thanks to all contributors!
#
#Please direct all questions to: mozes.blom@gmail.com

#########################
# Required modules
#########################
import os
import sys
import Bio
from Bio import AlignIO
import subprocess
import csv
from random import randrange
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pylab import *


#########################
# Input
#########################

# Main Directory
homePath = "/Specify/path/to/EAPhy-master/"

# Specify folder containing the best contigs for each locus. If haplotype data, both haplotypes for each indiv should be included in one locus fasta file, but with seq ID/names indicating H0 or H1 (see individuals/datatype)
bestcontigPath = os.path.join(homePath, 'examples', 'haplotype', 'test_haplo_contigs')

# Specify the list with all individuals to be included in the analysis, names should exactly match the Fasta sequence identifiers!!! If haplotype data, make sure that H1 and H0 are specified (see 'datatype')
individuals = os.path.join(homePath,'examples', 'test_indivs_haplo_list.txt')

# Specify the list with all loci to be included in the analysis, names should exactly match the Fasta files
exons = os.path.join(homePath, 'examples', 'test_genes_list.txt')

# Basename of output folder and all outputfiles 
out_name = "Test"

# Specify the folder with all dependent Python scripts
path_to_scripts = os.path.join(homePath, 'scripts')
sys.path.append(path_to_scripts)

# Specify the exact path to the aligner executable
aligner = "/Specify/path/to/muscle3.8.31"


#########################
# Which Analysis to run? 1 = yes, 0 = no
#########################

complete = 1        # Run complete analysis

subset_indivs = 0   # Create contigs with individuals as specified in 'individuals'
alignments = 0      # Generate the Muscle alignments
align_check = 0     # Conduct the amino-acid coding alignment checks, as described in REF
align_stats = 0     # Calculate alignment statistics
align_final = 0     # Generate the final alignments, this can include a subset of the 'individuals', to be specified in 'indivs included' and can be done with various degrees of missing data (individuals)
concatenation = 0   # Concatenation, only to be run when datatype = ambiguous!! Does not work if datatype = haplotype
snp_selection = 0   # SNP selection, only to be run when datatype = ambiguous!! Does not work if datatype = haplotype


#########################
# Datatype? Either 'haplotype' or 'ambiguous'
#########################

datatype = 'haplotype'  # If Datatype == 'haplotype', then this should be also indicated in the individuals file, with a single indiv having two lines followed by '_h0 and _h1': i.e. indiv_A_h0 and indiv_A_h1. The contig files should still be per gene and contain two entries for each individual, for each haplotype


#########################
# Parameter settings
#########################

## Alignment filtering
max_stop = 1                    # Max number of stop codons allowed for each sequence in each alignment (checked simultaneously for 1st, 2nd and 3rd frame, if required filt_align put in correct frame), otherwise flagged in 'alignments_to_check.txt' output file
indiv_gap_window = 7            # Jump-sliding-window length (in amino acids) used to filter out individual gaps
indiv_gap_ratio = 0.5           # Gap ratio (gaps/window length) used to filter out individual gaps
column_gap_window = 3           # Jump-sliding-window length (in amino acids) used to filter out codon columns with missing data, USE UNEVEN NUMBER SINCE ACCEPT/REJECT RATIO DEPENDS ON WHETHER MORE THAN HALF COLUMNS GAPS PRESENT/ABSENT + IF YOU SET 1, NOTE THAT CODONS WITH HETEROZYGOUS SITES LEADING TO A NON-SYNON. SUBSTITUTION DEPENDENT ON ALLELE ARE REMOVED! NEEDS FIXING
column_gap_ratio = 0.01         # Gap ratio (gaps/column) used to filter out codon columns with missing data
max_insertion_cut_off = 1       # The max-number of indiv. sequenced for each nt-column for it to be considered an insertion and thus removed
sliding_window_div_check = 7    # Normal sliding-window length (in amino acids) used to compare individual sequences to consensus sequence
max_heterozygosity_ratio = 0.5  # Max ratio of amino acids that differ from consensus sequence, otherwise flagged in 'alignments_to_check.txt' output file

average_heterozygosity_cut_off = 1  # Cut-off ratio used to flag potential paralogs (0 = all flagged, 1 = non-flagged), alignments with highest average heterozygosity levels, flagged in 'potential_paralogs.txt' output file

## Number of missing indivs allowed
indivs_included = individuals   # "path/to/list/"   # NOTE normally this should be the same as the individual list as specified above, however I kept this input file here so if you really want you can subselect individuals from the above list and you don't have to run the alignments again
number_missing_indivs_allowed = [0, 2, 4] #(List of values between 0 (all indivs incl.) and total number of indivs, generates datasets with various amounts of missing data (individuals in this case))

## Minimum alignment length
min_align_length = 150          # Minimum length of alignments, to be included in the concatenated alignment

## Set parameters settings for input file Partition Finder
branch = "linked"       ## BRANCHLENGTHS: linked | unlinked ##
model_evol = "raxml"    ## MODELS OF EVOLUTION for PartitionFinder: all | raxml | mrbayes | beast | <list> ##
model_sel = "BIC"       ## MODEL SELECCTION: AIC | AICc | BIC #
search = "rcluster"     ## SCHEMES, search: all | greedy | rcluster | hcluster | user ##


## SNP selection
biallelic = 1           # Biallelic enforced? 1 = yes, 0 = no. If yes, then SNP positions can only be biallelic, sites that contain genotypes that would result in more than 2 alleles are not considered. If no, ALL polymorphic sites are considered.
reps = 3                # How often should the random SNP selection be repeated?
unique = 1              # Sample only one SNP per gene? 1 = yes, 0 = no. If yes, then MAKE SURE THAT YOUR EXONS ARE NAMED IN FOLLOWING FORMAT (and accordingly specified in the exon list): gene_exonNumber_species.fasta. So EXACTLY 2 underscores and starting with the gene name. It will then make sure that only one snp per gene is sampled (i.e. a SNAPP requirement)

#########################
# Automated Analysis
#########################

# First check whether the Data type is specified correctly:
if datatype == 'ambiguous':
    print 'Diplotype sequence data with heterozygous sites coded according to IUPAC'
elif datatype == 'haplotype':
    print 'Phased haplotype data'
else:
    sys.exit("The datatype is incorrectly specified, please check manual and double check the config script")

# Create a general results folder in case it doesn't exist yet
resultsDir = os.path.join(homePath, 'results')
if (complete == 1) or (subset_indivs == 1):
    if os.path.isdir(resultsDir):
        print 'A general results folder already exists, continue with run specific analysis'
    else:
        subprocess.call("mkdir '%s'" % (resultsDir), shell=True)
else:
    pass

# Create a run specific output folder and datatype folder (the actual output folder)
outDir_run = os.path.join(resultsDir, out_name)
outDir = os.path.join(outDir_run, datatype)
if (complete == 1) or (subset_indivs == 1):
    if os.path.isdir(outDir_run):
        if os.path.isdir(outDir):
            sys.exit("The run specific output folder %s already exists, please specify a new output name" % (outDir))
        else:
            subprocess.call("mkdir %s" % (outDir), shell=True)
    else:
        subprocess.call("mkdir '%s'" % (outDir_run), shell=True)
        subprocess.call("mkdir '%s'" % (outDir), shell=True)
else:
    pass


## Run analysis

# get subset of individuals for alignment
outDir_subset = os.path.join(outDir, "1.best_contigs")
exon_subset_list = []
if (complete == 1) or (subset_indivs == 1):
    import subsetIndivs
    if os.path.isdir(outDir_subset):
        pass
    else:
        subprocess.call("mkdir '%s'" % (outDir_subset), shell=True)
    target_exons_file = open(exons, 'rU')
    for exon in target_exons_file:
        name = exon.rsplit('.fasta', 1)[0]
        exon_subset_list.append(subsetIndivs.createSubsetContigs(name, individuals, bestcontigPath, outDir_subset))
    target_exons_file.close()
    print 'Best contigs subsetted to include only individuals as specified in %s' % (individuals)
else:
    pass


# run alignments
outDir_align = os.path.join(outDir, "2.alignments")
outDir_align_basic = os.path.join(outDir_align, "1.basic_alignments")
if (complete == 1) or (alignments == 1):
    import runAlignment
    if os.path.isdir(outDir_align):
        pass
    else:
        subprocess.call("mkdir '%s'" % (outDir_align), shell=True)
    if os.path.isdir(outDir_align_basic):
        pass
    else:
        subprocess.call("mkdir '%s'" % (outDir_align_basic), shell=True)
        for subset_contig in exon_subset_list:
            runAlignment.alignments(subset_contig, outDir_subset, outDir_align_basic, aligner)
else:
    pass



# filter Alignments
outDir_align_check = os.path.join(outDir_align, "2.checked_alignments")
outDir_align_check_temp = os.path.join(outDir_align_check, "alignments_temp")
outDir_align_check_passed = os.path.join(outDir_align_check, "1.alignments_passed")
outDir_align_check_passed_alignments = os.path.join(outDir_align_check_passed, "alignments")
outDir_align_check_inspect = os.path.join(outDir_align_check, "2.alignments_to_check")
outDir_align_check_inspect_alignments = os.path.join(outDir_align_check_inspect, "alignments")
if datatype == 'haplotype':
    outDir_align_check_filtered = os.path.join(outDir_align_check, "3.alignments_filtered")
else:
    outDir_align_check_filtered = os.path.join(outDir_align_check, "4.alignments_filtered")
if (complete == 1) or (align_check == 1):
    import check_frame_gaps
    import check_paralogs
    if os.path.isdir(outDir_align_check):
        sys.exit("A checked alignments folder already exists, to avoid overwriting of existing data, pipeline exit")
    else:
        subprocess.call("mkdir '%s'" % (outDir_align_check), shell=True)
        subprocess.call("mkdir '%s'" % (outDir_align_check_temp), shell=True)
        subprocess.call("mkdir '%s'" % (outDir_align_check_passed), shell=True)
        subprocess.call("mkdir '%s'" % (outDir_align_check_passed_alignments), shell=True)
        subprocess.call("mkdir '%s'" % (outDir_align_check_inspect), shell=True)
        subprocess.call("mkdir '%s'" % (outDir_align_check_inspect_alignments), shell=True)
        subprocess.call("mkdir '%s'" % (outDir_align_check_filtered), shell=True)
    list_passed_align = os.path.join(outDir_align_check_passed, 'alignments_passed.txt')
    list_inspect_align = os.path.join(outDir_align_check_inspect, 'alignments_require_inspection.txt')
    list_filtered_align = os.path.join(outDir_align_check_filtered, 'alignments_filtered.txt')
    file_handle_correct = open(list_passed_align, 'w')
    file_handle_wrong = open(list_inspect_align, 'w')
    file_handle_filtered = open(list_filtered_align, 'w')
    empty_align = []
    for alignment in os.listdir(outDir_align_basic):
        if alignment.endswith("_basic_align.fasta"):
            alignment_path = os.path.join(outDir_align_basic, alignment)
            align_length_zero = check_frame_gaps.removeGaps(alignment, alignment_path, outDir_align_check_temp, indiv_gap_window, indiv_gap_ratio, column_gap_window, column_gap_ratio, max_insertion_cut_off)
            if align_length_zero == True: #If the alignment is empty we don't want to include it and therefore needs to be flagged!!
                empty_align_base_name = alignment.rsplit('_basic_align.fasta', 1)[0]
                empty_align.append(empty_align_base_name)
            else:
                pass                    
        else:
            pass
    loci_correct = []
    for new_alignment in os.listdir(outDir_align_check_temp):
        if new_alignment.endswith("_gaps_analyzed.fasta"):
            new_align_base_name = new_alignment.rsplit('_gaps_analyzed.fasta', 1)[0]
            new_alignment_path = os.path.join(outDir_align_check_temp, new_alignment) 
            if new_align_base_name in empty_align:
                file_handle_wrong.write(new_alignment + '\n')
                subprocess.call("mv '%s' '%s'" % (new_alignment_path, outDir_align_check_inspect_alignments), shell=True)
            else:
                if (check_frame_gaps.identifyStopCodons(new_alignment_path, max_stop) == True) and (check_frame_gaps.identifyDivAlign(new_alignment_path, sliding_window_div_check, max_heterozygosity_ratio) == True):
                    loci_correct.append(new_alignment)
                else:
                    file_handle_wrong.write(new_alignment + '\n')
                    subprocess.call("mv '%s' '%s'" % (new_alignment_path, outDir_align_check_inspect_alignments), shell=True)
        else:
            pass
    if datatype == 'ambiguous':
        outDir_align_check_paralog = os.path.join(outDir_align_check, "3.alignments_paralogy")
        subprocess.call("mkdir '%s'" % (outDir_align_check_paralog), shell=True)
        list_paralog_align = os.path.join(outDir_align_check_paralog, 'alignments_potential_paralogs.txt')
        file_handle_paralog = open(list_paralog_align, 'w')
        average_heterozygosity = []
        filtered_align_counter = 0
        for filtered_align in loci_correct:
            filtered_align_counter += 1
            new_alignment_path = os.path.join(outDir_align_check_temp, filtered_align)
            average_heterozygosity.append(check_paralogs.calcAverageIndivHet(new_alignment_path))
        average_heterozygosity.sort()
        cut_off_value = int(average_heterozygosity_cut_off * filtered_align_counter)
        cut_off = average_heterozygosity[cut_off_value - 1]
        check_paralogs.plotAverageHet(average_heterozygosity, cut_off, outDir_align_check_paralog)
        for filtered_align in loci_correct:
            new_alignment_path = os.path.join(outDir_align_check_temp, filtered_align)
            if (check_paralogs.calcAverageIndivHet(new_alignment_path) <= cut_off) == True:
                file_handle_correct.write(filtered_align + '\n')
                file_handle_filtered.write(filtered_align + '\n')
                subprocess.call("mv '%s' '%s'" % (new_alignment_path, outDir_align_check_passed_alignments), shell=True)
            else:
                file_handle_paralog.write(new_alignment + '\n')
                subprocess.call("mv '%s' '%s'" % (new_alignment_path, outDir_align_check_inspect_alignments), shell=True)
        file_handle_paralog.close()
    else:
        for filtered_align in loci_correct:
            new_alignment_path = os.path.join(outDir_align_check_temp, filtered_align)
            file_handle_correct.write(filtered_align + '\n')
            file_handle_filtered.write(filtered_align + '\n')
            subprocess.call("mv '%s' '%s'" % (new_alignment_path, outDir_align_check_passed_alignments), shell=True)
    file_handle_wrong.close()
    file_handle_correct.close()
    file_handle_filtered.close()
    subprocess.call("rm -r '%s'" % (outDir_align_check_temp), shell=True)
    subprocess.call("cp -r '%s' '%s'" % (outDir_align_check_passed_alignments, outDir_align_check_filtered), shell=True)
else:
    pass


# Calculate alignment statistics
if datatype == 'ambiguous':
    outDir_align_check_stats = os.path.join(outDir_align_check, "5.filtered_align_stats")
else:
    outDir_align_check_stats = os.path.join(outDir_align_check, "4.filtered_align_stats")
if (complete == 1) or (align_stats == 1):
    import filt_align_stats
    if os.path.isdir(outDir_align_check_stats):
        pass
    else:
        subprocess.call("mkdir '%s'" % (outDir_align_check_stats), shell=True)
    align_stats_out = os.path.join(outDir_align_check_stats, (out_name + "_contig_stats.csv"))
    file_handle = open(align_stats_out, 'wt')
    writer = csv.writer(file_handle)
    writer.writerow(('Locus', 'Number of samples', 'Average contig length'))
    list_filtered_loci = os.path.join(outDir_align_check_filtered, 'alignments_filtered.txt')
    filtered_loci_file = open(list_filtered_loci, 'rU')
    for gene in filtered_loci_file:
        name = gene.rsplit('\n', 1)[0]
        if name.endswith("_gaps_analyzed.fasta"):
            base_name = name.rsplit('_gaps_analyzed.fasta', 1)[0]
            align_path = os.path.join(outDir_align_check_filtered, 'alignments', name)
            x, y, z = filt_align_stats.getContigStats(base_name, align_path)
            writer.writerow((x, y, z))
        else:
            pass
    file_handle.close()
else:
    pass


# Select loci to be used; account for % missing individuals
outDir_final_align = os.path.join(outDir_align, "3.final_alignments")
if (complete == 1) or (align_final == 1):
    import final_align
    if os.path.isdir(outDir_final_align):
        pass
    else:
        subprocess.call("mkdir '%s'" % (outDir_final_align), shell=True)
    list_loci_in_frame = os.path.join(outDir_align_check_filtered, 'alignments_filtered.txt')
    outDir_align_check_filtered_alignments = os.path.join(outDir_align_check_filtered, 'alignments')
    values_used = []
    loci_number_recovered = []
    for number in number_missing_indivs_allowed:
        outDir_final_align_number = os.path.join(outDir_final_align, ('%s_missing_indivs_allowed') % number)
        subprocess.call("mkdir '%s'" % (outDir_final_align_number), shell=True)
        loci_number_recovered.append(final_align.select_final_alignments(list_loci_in_frame, outDir_align_check_filtered_alignments, indivs_included, number, datatype, outDir_final_align_number))
        values_used.append(number)        
    if len(number_missing_indivs_allowed) > 1:
        figure(2)
        plt.plot(values_used, loci_number_recovered, 'ro')
        plt.axis([0, (max(number_missing_indivs_allowed) + 1), 0, (max(loci_number_recovered) + 50)])
        plt.title('Effect of missing data on loci number recovered')
        plt.xlabel('Number of missing individuals per alignment allowed')
        plt.ylabel('Number of loci recovered')
        plot_out = os.path.join(outDir_final_align, 'NumberOfLociRecoveredMissingData.png')
        plt.savefig(plot_out, dpi=128)
    else:
        pass
else:
    pass



if datatype == 'ambiguous':

# Concatenate alignments
    outDir_concat = os.path.join(outDir, "3.concatenation")
    if (complete == 1) or (concatenation == 1):
        import concat
        if os.path.isdir(outDir_concat):
            pass
        else:
            subprocess.call("mkdir '%s'" % (outDir_concat), shell=True)
        for number in number_missing_indivs_allowed:
            outDir_final_align_number = os.path.join(outDir_final_align, ('%s_missing_indivs_allowed') % number)
            list_final_loci = os.path.join(outDir_final_align_number, ('loci_included_%s_indivs_max_missing.txt') % number)
            outDir_final_align_number_alignments = os.path.join(outDir_final_align_number, 'alignments')
            outDir_concat_number = os.path.join(outDir_concat, ('%s_missing_indivs_allowed') % number)
            subprocess.call("mkdir '%s'" % (outDir_concat_number), shell=True)
            concat.concatenateLoci(list_final_loci, indivs_included, outDir_final_align_number_alignments, min_align_length, outDir_concat_number, out_name,  branch, model_evol, model_sel, search)
    else:
        pass


    # SNP selection
    outDir_snp = os.path.join(outDir, "4.SNPselection")
    if (complete == 1) or (snp_selection == 1):
        import snp_selection
        if os.path.isdir(outDir_snp):
            pass
        else:
            subprocess.call("mkdir '%s'" % (outDir_snp), shell=True)
        for number in number_missing_indivs_allowed:
            outDir_final_align_number = os.path.join(outDir_final_align, ('%s_missing_indivs_allowed') % number)
            list_final_loci = os.path.join(outDir_final_align_number, ('loci_included_%s_indivs_max_missing.txt') % number)
            outDir_final_align_number_alignments = os.path.join(outDir_final_align_number, 'alignments')
            outDir_snp_number = os.path.join(outDir_snp, ('%s_missing_indivs_allowed') % number)
            subprocess.call("mkdir '%s'" % (outDir_snp_number), shell=True)
            try:
                loci_file = open(list_final_loci, 'rU')
            except IOError:
                print "Could not open file with final loci"
            else:
                loci_counter = 0
                for line in loci_file:
                    loci_counter += 1
                loci_file.seek(0)
                if loci_counter > 1:
                    for rep in range(reps):
                        snp_align_counter = 0
                        loci_included = []
                        for alignment in loci_file:
                            align_name = alignment.rsplit('\n', 1)[0]
                            if unique == 1:
                                align_base_name = alignment.rsplit('_', 4)[0]
                            else:
                                align_base_name = alignment.rsplit('_final_align.fasta', 1)[0]
                            if not align_base_name in loci_included:
                                alignment_path = os.path.join(outDir_final_align_number, 'alignments', align_name)
                                if biallelic == 1:
                                    output_site_align = snp_selection.findSite_Biallelic(alignment_path)
                                    output_site_align.sort()
                                elif biallelic == 0:
                                    output_site_align = snp_selection.findSite(alignment_path)
                                    output_site_align.sort()
                                else:
                                    sys.exit("Biallelic yes/no not correctly specified, SNP selection failed")
                                if (output_site_align.get_alignment_length()) > 0:
                                    if snp_align_counter == 0:
                                        upper_limit = output_site_align.get_alignment_length()
                                        random_snp = randrange(0, (upper_limit), 1)
                                        snp_align_one = output_site_align[:, (random_snp):(random_snp + 1)]
                                        snp_align_all = output_site_align
                                        snp_align_counter += 1
                                        loci_included.append(align_base_name)
                                    else:
                                        upper_limit = output_site_align.get_alignment_length()
                                        random_snp = randrange(0, (upper_limit), 1)
                                        snp_align_one = snp_align_one + output_site_align[:, (random_snp):(random_snp + 1)]
                                        snp_align_all = snp_align_all + output_site_align
                                        snp_align_counter += 1
                                        loci_included.append(align_base_name)
                                else:
                                    pass
                            else:
                                alignment_path = os.path.join(outDir_final_align_number,'alignments', align_name)
                                output_site_align = snp_selection.findSite(alignment_path)
                                output_site_align.sort()
                                if (output_site_align.get_alignment_length()) > 0:
                                    snp_align_all = snp_align_all + output_site_align
                                else:
                                    pass                       
                        out_file_name = os.path.join(outDir_snp_number, ('snp_align_one_snp_gene_r_%s.fasta' % rep))
                        file_out_handle = open(out_file_name, 'w')
                        AlignIO.write(snp_align_one, file_out_handle, 'fasta')
                        file_out_handle.close()
                        out_file_name = os.path.join(outDir_snp_number, ('included_loci_one_snp_gene_r_%s.txt' % rep))
                        file_out_handle = open(out_file_name, 'w')
                        for line in loci_included:
                            file_out_handle.write(line + '\n')
                        file_out_handle.close()
                        out_file_name = os.path.join(outDir_snp_number, ('snp_align_all_snps_gene_r_%s.fasta' % rep))
                        file_out_handle = open(out_file_name, 'w')
                        AlignIO.write(snp_align_all, file_out_handle, 'fasta')
                        file_out_handle.close()
                        loci_file.seek(0)
                else:
                    file_handle_out = os.path.join(outDir_snp_number, ('%s_missing_indivs_allowed_error_report.txt') % number)
                    summary_file = open(file_handle_out, "w")
                    summary_file.write('No genes were present in the final alignment file (perhaps due to stringency requirement missing data allowed??)')
                    summary_file.close()
    else:
        pass
else:
    print 'Haplotype data cannot be concatenated or used for SNP selection, if aim is concatenation or snp selection, please specify and use ambiguous coded diplotype sequence data'
