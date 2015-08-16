import os
import subprocess

def alignments(subset_contig, outDir_subset, outDir_align_basic, aligner, subset):
	input_file_path = os.path.join(outDir_subset, subset_contig)
	if subset == True:
		base_name = subset_contig.rsplit('_subset.fasta', 1)[0]
	else:
		base_name = subset_contig.rsplit('.fasta', 1)[0]
	output_file_path = os.path.join(outDir_align_basic, (base_name + "_basic_align.fasta"))
	subprocess.call("%s -in %s -out %s" % (aligner, input_file_path, output_file_path), shell = True)
