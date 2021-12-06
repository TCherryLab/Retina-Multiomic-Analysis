#!/usr/bin/env python
import subprocess
import os
import glob

##parameters
delim = '\t'
##programs
bedtools = 'bedtools'

##methods
def get_list_of_genes(infile):
	genes = []
	with open(infile, "r") as in_fh:
		for line in in_fh:
			line = line.rstrip().split(delim)
			gene = line[0]
			genes.append(gene)
	return(genes)

def make_dict_from_info_file(infile, file_type):
	peak_dict = {}
	with open(infile, "r") as in_fh:
		lc = 0
		for line in in_fh:
			lc += 1
			if lc > 1:
				line = line.rstrip().split(delim)
				idx_num = line[0]
				if file_type == 'rna':
					peak_info = line[6]
				elif file_type == 'peak':
					peak_info = line[1] + ":" + line[2] + "-" + line[3]
				peak_dict[idx_num] = peak_info
	return(peak_dict)

def format_p2g_file(infile, outfile, p2g_dict, atac_dict):
	with open(infile, "r") as in_fh, open(outfile, "w") as out_fh:
		lc = 0
		for line in in_fh:
			lc += 1
			line = line.rstrip().split(delim) 
			if lc == 1:
				header = line + ['peak', 'gene']
				out_fh.write(delim.join(header) + '\n')
			else:
				idxATAC = line[0]
				idxRNA = line[1]
				atac_info = atac_dict[idxATAC]
				rna_info = p2g_dict[idxRNA]
				line_out = line + [atac_info, rna_info]
				out_fh.write(delim.join(line_out) + '\n')

def filter_p2g_file_get_peaks(infile, outfile, ret_genes):
	peak_list = []
	with open(infile, "r") as in_fh, open(outfile, "w") as out_fh:
		lc = 0
		for line in in_fh:
			lc += 1
			line = line.rstrip().split(delim) 
			if lc == 1:
				out_fh.write(delim.join(line) + '\n')
			else:
				gene = line[7]
				peak = line[6]
				if gene in ret_genes:
					out_fh.write(delim.join(line) + '\n')
					peak_list.append(peak)
	return(peak_list)

def filter_heatmap_data(infile, outfile, ret_peaks):
	temp_file = outfile.rsplit('.', 1)[0] + '.temp.txt'
	temp2_file = outfile.rsplit('.', 1)[0] + '.bed'
	##filter for ret peaks
	with open(infile, "r") as in_fh, open(temp_file, "w") as out_fh:
		lc = 0
		for line in in_fh:
			lc += 1
			line = line.replace('"', "").rstrip().split(',')
			if lc == 1:
				out_fh.write(delim.join(['chr', 'start', 'end'] + line) + '\n')
			else:
				peak = line[0]
				peak_split = peak.replace('-', ":").split(':')
				if peak in ret_peaks:
					out_fh.write(delim.join(peak_split + line) + '\n')
	##sort file
	with open(temp2_file, "w") as t2_fh:
		sort_file = subprocess.Popen(['sort', '-k1,1', '-k2,2n', temp_file], stdout=t2_fh)
		sort_file.wait()
	##remove bed type cols used to sort
	with open(temp2_file, "r") as in_fh, open(outfile, "w") as out_fh:
		for line in in_fh:
			line = line.split(delim)
			out_fh.write(delim.join(line[3:]))

def make_file_for_heatmap(file_prefix, data_all, data_markers, peak2gene, rna_index, atac_index, gene_list):
	##get peak/gene info for p2g data
	rna_dict = make_dict_from_info_file(rna_index, 'rna')
	atac_dict = make_dict_from_info_file(atac_index, 'peak')
	##reformat p2g file, and then filter/get peaks related to ret genes
	formatted_p2g = file_prefix + '.peak2gene.all.xls' 
	filtered_p2g = file_prefix + '.peak2gene.ret_genes.xls' 
	format_p2g_file(peak2gene, formatted_p2g, rna_dict, atac_dict)
	ret_gene_peaks = filter_p2g_file_get_peaks(formatted_p2g, filtered_p2g, gene_list)
	##filter the heatmap data
	filtered_all_heatmap_file = file_prefix + '.ret_genes.heatmap_all.txt' 
	filtered_marker_heatmap_file = file_prefix + '.ret_genes.heatmap_markers.txt' 
	filter_heatmap_data(data_all, filtered_all_heatmap_file, ret_gene_peaks)
	filter_heatmap_data(data_markers, filtered_marker_heatmap_file, ret_gene_peaks)

def compare_groups_pairwise(bed1, bed2, out_file, outfile_cols, sample_prefixes):
	##get header for bed
	with open(bed1, 'r') as b1, open(bed2, 'r') as b2:
		b1_head = b1.readline().rstrip().split(delim)
		b1_head = [sample_prefixes[0] + i for i in b1_head]
		b2_head = b2.readline().rstrip().split(delim)
		b2_head = [sample_prefixes[1] + i for i in b2_head]
	header = b1_head + b2_head
	print(header)
	temp_bed = out_file.rsplit('.', 1)[0] + '.bed'
	with open(temp_bed, "w") as out_fh: 
		bt_intersect = subprocess.Popen([bedtools, 'intersect', '-a', bed1, '-b', bed2, '-f', '0.50', '-wa', '-wb'], stdout=out_fh)
		bt_intersect.wait()
	with open(temp_bed, "r") as in_fh, open(out_file, "w") as out_fh:
		out_fh.write(delim.join(header[3:outfile_cols[0]] + header[outfile_cols[1]:]) + '\n')
		for line in in_fh:
			line = line.split(delim)
			out_fh.write(delim.join(line[3:outfile_cols[0]] + line[outfile_cols[1]:]))	


##run methods

##get list of retina genes from file
gene_file = 'retina_genes.txt'
gene_list = get_list_of_genes(gene_file)

##human data
heatmap_data_all = 'heatmap_cell_class_all_peaks.macs2_q0.000001.matrix.csv'
heatmap_data_markers = 'heatmap_cell_class_marker_peaks.macs2_q0.000001.matrix.csv'
human_peak2gene = 'human_all.min_cell20.p2g.m2_q0.000001.cor4.txt'
human_rna_index = 'human_all.min_cell20.m2_q0.000001.rnaseq_info.txt'
human_atac_index = 'human_all.min_cell20.peak_info.m2_q0.000001.txt'
human_prefix = 'human'
make_file_for_heatmap(human_prefix, heatmap_data_all, heatmap_data_markers, human_peak2gene, human_rna_index, human_atac_index, gene_list)


##organoid data
heatmap_data_all =  'heatmap_cell_class_all_peaks.macs2_q0.000001.matrix.csv'
heatmap_data_markers = 'heatmap_cell_class_marker_peaks.macs2_q0.000001.matrix.csv'
organoid_peak2gene = 'organoid_all.min_cell20.p2g.m2_q0.000001.cor4.txt'
organoid_rna_index = 'organoid_all.min_cell20.m2_q0.000001.rnaseq_info.txt'
organoid_atac_index = 'organoid_all.min_cell20.peak_info.m2_q0.000001.txt'
organoid_prefix = 'organoid'
make_file_for_heatmap(organoid_prefix, heatmap_data_all, heatmap_data_markers, organoid_peak2gene, organoid_rna_index, organoid_atac_index, gene_list)


##mouse data
heatmap_data_all = 'heatmap_cell_class_all_peaks.macs2_q0.000001.matrix.csv'
heatmap_data_markers = 'heatmap_cell_class_marker_peaks.macs2_q0.000001.matrix.csv'
mouse_peak2gene = 'mouse_all.min_cell20.p2g.m2_q0.000001.cor4.txt'
mouse_rna_index = 'mouse_all.min_cell20.m2_q0.000001.rnaseq_info.txt'
mouse_atac_index = 'mouse_all.min_cell20.peak_info.m2_q0.000001.txt'
mouse_prefix = 'mouse'
gene_file = 'mouse_retina_genes.txt'
gene_list = get_list_of_genes(gene_file)
make_file_for_heatmap(mouse_prefix, heatmap_data_all, heatmap_data_markers, mouse_peak2gene, mouse_rna_index, mouse_atac_index, gene_list)

	 
##compare human and organoid data 
human_heatmap_all = 'human.ret_genes.heatmap_all.bed'
human_heatmap_marker = 'human.ret_genes.heatmap_markers.bed'
organoid_heatmap_all = 'organoid.ret_genes.heatmap_all.bed'
organoid_heatmap_marker = 'organoid.ret_genes.heatmap_markers.bed'
human_org_heatmap_all = 'human_organoid.ret_genes.heatmap_all.int.txt'
human_org_heatmap_marker = 'human_organoid.ret_genes.heatmap_markers.int.txt'
human_org_heatmap_col_numbers = [22,26]
human_org_sample_prefix = ['hu_', 'org_']
compare_groups_pairwise(human_heatmap_all, organoid_heatmap_all, human_org_heatmap_all, human_org_heatmap_col_numbers, human_org_sample_prefix)
compare_groups_pairwise(human_heatmap_marker, organoid_heatmap_marker, human_org_heatmap_marker, human_org_heatmap_col_numbers, human_org_sample_prefix)

##compare human and mouse data 
human_heatmap_all = 'human.ret_genes.heatmap_all.bed'
human_heatmap_marker = 'human.ret_genes.heatmap_markers.bed'
mouse_heatmap_all = 'mouse.ret_genes.heatmap_all.hg38.bed'
mouse_heatmap_marker = 'mouse.ret_genes.heatmap_markers.hg38.bed'
human_mouse_heatmap_all = 'human_mouse.ret_genes.heatmap_all.int.txt'
human_mouse_heatmap_marker = 'human_mouse.ret_genes.heatmap_markers.int.txt'
human_mouse_heatmap_col_numbers = [22,26]
human_mouse_sample_prefix = ['hu_', 'mm_']
compare_groups_pairwise(human_heatmap_all, mouse_heatmap_all, human_mouse_heatmap_all, human_mouse_heatmap_col_numbers, human_mouse_sample_prefix)
compare_groups_pairwise(human_heatmap_marker, mouse_heatmap_marker, human_mouse_heatmap_marker, human_mouse_heatmap_col_numbers, human_mouse_sample_prefix)

