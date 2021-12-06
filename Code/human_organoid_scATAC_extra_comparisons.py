#!/usr/bin/env python
import subprocess
import os
import glob


##parameters
delim = '\t'

##methods

def get_correlation_info_homer_pairwise(cc_dict, d_values, size_values, bed_suffixes, tag_dir_suffix):
	for bed_suffix in bed_suffixes:
		for cell_class in cc_dict:
			samples = cc_dict[cell_class]
			##combine bed files
			merge_prefix = cell_class + '_merged.'
			annotate_prefix = cell_class + '_annotated.'
			annotate_suffix = bed_suffix.rsplit('.',1)[0] + '.txt'
			corr_prefix = cell_class + '_correlation.'
			# '''
			comb_beds = []
			for s in samples:
				bed = s + bed_suffix
				comb_beds.append(bed)
			for d in d_values:
				merge_peaks(merge_prefix, d, comb_beds, bed_suffix)
			##annotate files
			tag_dir_suffix = '.tag_dir'
			# correlation_prefix = 'correlation.'
			for d in d_values:
				merged_bed = merge_prefix + d + 'd' + bed_suffix
				##annotate those peaks, and then make a correlation file
				annotate_peaks(samples, merged_bed, d, tag_dir_suffix, annotate_prefix, size_values, annotate_suffix)
			# '''
			##format annotation files
			for size_req in size_values:
				for d_value in d_values:
					infile = annotate_prefix + d_value + 'd.' + size_req + 'size' + annotate_suffix
					out_file = corr_prefix + infile.split('.', 1)[1]
					# print(infile, out_file)
					with open(out_file, 'w') as out_fh, open(infile, 'r') as in_fh:
						lc = 0
						for line in in_fh:
							lc += 1
							line = line.strip('\n').split(delim)
							if lc == 1:
								sample_info = line[19:]
								sample_info = [s.split('.')[0] + '.' + s.split('.')[1] for s in sample_info]
								out_fh.write(delim.join(['peak_id'] + sample_info) + '\n')
								# print(sample_info)
							else:
								chrom = line[1]
								if '_' not in chrom and chrom != 'chrM' and chrom != 'chrY':
									human_count = float(line[19])
									org_count = float(line[20])
									# if min([human_count, org_count]) > 1:
									line_out = [line[0]] + line[19:]
									out_fh.write(delim.join(line_out) + '\n')


def combine_tag_dirs(samples_to_combine, combined_name):
	outdir = combined_name + '.tag_dir'
	tag_dir_cmds = []
	for sample in samples_to_combine:
		tag_dir_cmd = [ '-d', sample + '.tag_dir']
		tag_dir_cmds.extend(tag_dir_cmd)
	print(tag_dir_cmds)
	mk_tag_dir = subprocess.Popen(['makeTagDirectory', outdir] + tag_dir_cmds)
	mk_tag_dir.wait()


def get_diff_peaks_diff_sig(out_prefix1, out_prefix2, d_values, q_values, p_values, fc_values, sizes_req, human_name, org_name):
	human_tag_dir = human_name  + '.tag_dir'
	org_tag_dir = org_name  + '.tag_dir'
	for d_value in d_values:
		for q_value in q_values:
			for p_value in p_values:
				for fc_value in fc_values:
					for size_req in sizes_req:
						merged_file = 'merged.' + d_value + 'd.macs2.bampe_q' + q_value + '_keepdups_summits.txt'
						out_file1 = out_prefix1 + 'size' + size_req + '.p_' + p_value + '.fc_' + fc_value + '.'  + merged_file.split('.', 1)[1]
						out_file2 = out_prefix2 + 'size' + size_req + '.p_' + p_value + '.fc_' + fc_value + '.'  + merged_file.split('.', 1)[1]
						with open(out_file1, 'w') as out_fh:
							# mk_tag_dir = subprocess.Popen(['mergePeaks', '-d', d_value, '-matrix', out_prefix + d_value] + peak_files, stdout=out_fh)
							run_diff_peaks1 = subprocess.Popen(['getDifferentialPeaks', merged_file, human_tag_dir, org_tag_dir, '-F', fc_value, '-P', p_value, '-size', size_req], stdout=out_fh)
							run_diff_peaks1.wait()
						with open(out_file2, 'w') as out_fh:
							# mk_tag_dir = subprocess.Popen(['mergePeaks', '-d', d_value, '-matrix', out_prefix + d_value] + peak_files, stdout=out_fh)
							run_diff_peaks2 = subprocess.Popen(['getDifferentialPeaks', merged_file, org_tag_dir, human_tag_dir, '-F', fc_value, '-P', p_value, '-size', size_req], stdout=out_fh)
							run_diff_peaks2.wait()


def get_diff_peaks_replicates(out_prefix1, out_prefix2, human_names, org_names, all_peaks):
	human_tag_dirs = [i  + '.tag_dir' for i in human_names]
	org_tag_dirs = [i  + '.tag_dir' for i in org_names]
	out_file1 = out_prefix1 + 'balanced.txt'
	out_file2 = out_prefix2 + 'balanced.txt'
	if all_peaks == 'no':
		with open(out_file1, 'w') as out_fh:
			run_diff_peaks1 = subprocess.Popen(['getDifferentialPeaksReplicates.pl', '-t'] + human_tag_dirs + ['-b'] + org_tag_dirs + ['-balanced', '-genome', 'hg38'], stdout=out_fh)
			run_diff_peaks1.wait()
		with open(out_file2, 'w') as out_fh:
			run_diff_peaks2 = subprocess.Popen(['getDifferentialPeaksReplicates.pl', '-t'] + org_tag_dirs + ['-b'] + human_tag_dirs + ['-balanced', '-genome', 'hg38'], stdout=out_fh)
			run_diff_peaks2.wait()
	elif all_peaks == 'yes':
		with open(out_file1, 'w') as out_fh:
			run_diff_peaks1 = subprocess.Popen(['getDifferentialPeaksReplicates.pl', '-t'] + human_tag_dirs + ['-b'] + org_tag_dirs + ['-balanced', '-genome', 'hg38', '-all'], stdout=out_fh)
			run_diff_peaks1.wait()
		with open(out_file2, 'w') as out_fh:
			run_diff_peaks2 = subprocess.Popen(['getDifferentialPeaksReplicates.pl', '-t'] + org_tag_dirs + ['-b'] + human_tag_dirs + ['-balanced', '-genome', 'hg38', '-all'], stdout=out_fh)
			run_diff_peaks2.wait()

def find_peaks_cell_classes(names, peak_suffix):
	for name in names:
		tag_dir = name + '.tag_dir'
		peak_file = name + peak_suffix
		mk_tag_dir = subprocess.Popen(['findPeaks', tag_dir, '-o', peak_file])
		mk_tag_dir.wait()



def make_beds_from_diff_peaks_file(homer_dp_files):
	for homer_dp_file in homer_dp_files:
		out_bed = homer_dp_file.rsplit('.', 1)[0] + '.bed'
		with open(out_bed, 'w') as out_fh, open(homer_dp_file, 'r') as in_fh:
			lc = 0
			for line in in_fh:
				lc += 1
				line = line.strip('\n').split(delim)
				if lc > 1:
					chrom = line[1]
					if '.' not in chrom and '_' not in chrom:
						out_fh.write(delim.join(line[1:4]) + '\n')

def make_beds_from_find_peaks_files(samples, infile_suffix):
	for sample in samples:
		infile = sample + infile_suffix
		out_bed = infile.rsplit('.', 1)[0] + '.bed'
		with open(out_bed, 'w') as out_fh, open(infile, 'r') as in_fh:
			for line in in_fh:
				if line[0] != '#':
					line = line.strip('\n').split(delim)
					chrom = line[1]
					if '.' not in chrom and '_' not in chrom:
						out_fh.write(delim.join(line[1:4]) + '\n')

def bt_int_homer_beds_vs_cc_beds(homer_result_beds, cc_samples, cc_bed_suffix):
	cc_beds = [i + cc_bed_suffix for i in cc_samples]
	for homer_result_bed in homer_result_beds:
		out_bti = homer_result_bed.rsplit('.', 1)[0] + '.bt_int.mature_cc.txt'
		##bedtools intersect 
		with open(out_bti, "w") as out_fh: 
			hom_bt_intersect = subprocess.Popen([bedtools, 'intersect', '-a', homer_result_bed, '-b'] + cc_beds + ['-C', '-filenames'], stdout=out_fh)
			hom_bt_intersect.wait()

def filter_and_make_beds_from_diff_peaks_file(homer_dp_files, cut_offs_list):
	for homer_dp_file in homer_dp_files:
		for cut_offs in cut_offs_list:
			fc_wanted = cut_offs[0]
			p_wanted = cut_offs[1]
			out_bed = homer_dp_file.rsplit('.', 1)[0] +'.lfc' + str(fc_wanted) + '_p' + str(p_wanted) + '.bed'
			with open(out_bed, 'w') as out_fh, open(homer_dp_file, 'r') as in_fh:
				lc = 0
				for line in in_fh:
					lc += 1
					line = line.strip('\n').split(delim)
					if lc > 1:
						chrom = line[1]
						logfc = float(line[24])
						adjp = float(line[26])
						if '.' not in chrom and '_' not in chrom:
							if logfc >= fc_wanted and adjp <= p_wanted:
								out_fh.write(delim.join(line[1:4]) + '\n')

def format_int_files_for_upset(infiles):
	for infile in infiles:
		outfile = infile.rsplit('.', 1)[0] + '.upset.txt'
		with open(outfile, 'w') as out_fh, open(infile, 'r') as in_fh:
			out_fh.write(delim.join(['peak', 'cell_class', 'count']) + '\n')
			for line in in_fh:
				line = line.strip('\n').split(delim)
				peak = '_'.join(line[:3])
				cc = line[3].rsplit('.',3)[0]
				count = line[4]
				if int(count) > 1:
					count = '1'
				out_fh.write(delim.join([peak, cc, count]) + '\n')

def combine_peak_info(int_files, ann_file):
	##make dict from annfile
	ann_dict = {}
	with open(ann_file, 'r') as in_fh:
		lc = 0
		for line in in_fh:
			line = line.strip('\n').split(delim)
			lc += 1
			if lc == 1:
				ann_header = line
			else:
				peak = '_'.join(line[1:4])
				ann_dict[peak] = line
	for infile in int_files:
		outfile = infile.rsplit('.', 5)[0] + '.peak_info.txt'
		with open(outfile, 'w') as out_fh, open(infile, 'r') as in_fh:
			lc = 0
			for line in in_fh:
				line = line.strip('\n').split(',')
				lc += 1
				if lc == 1:
					int_header = [i.strip('"') for i in line[2:]]
					print(int_header)
					out_fh.write(delim.join(ann_header + int_header) + '\n')
				else:
					int_peak = line[1].strip('"')
					line_out = ann_dict[int_peak] + line[2:]
					out_fh.write(delim.join(line_out) + '\n')

def filter_peak_file_prs(infiles, cc_wanted):
	for infile in infiles:
		outfile = infile.rsplit('.', 2)[0] + '.' + cc_wanted + '.homer_peak'
		with open(outfile, 'w') as out_fh, open(infile, 'r') as in_fh:
			lc = 0
			for line in in_fh:
				line = line.strip('\n').split(delim)
				lc += 1
				if lc >1:
					# print(line[27:])
					counts = [int(i) for i in line[27:39]]
					if cc_wanted == 'org_rods_cones':
						if counts[6] != 0 or counts[8] != 0:
							counts.pop(8)
							counts.pop(6)
							# print(line[:5], counts)
							if sum(counts) == 0:
								out_fh.write(delim.join(line[:5]) + '\n')
					elif cc_wanted == 'human_rods_cones':
						if counts[2] != 0 or counts[5] != 0:
							counts.pop(5)
							counts.pop(2)
							# print(line[:5], counts)
							if sum(counts) == 0:
								out_fh.write(delim.join(line[:5]) + '\n')
					elif cc_wanted == 'human_rods':
						if counts[5] != 0:
							counts.pop(5)
							# print(line[:5], counts)
							if sum(counts) == 0:
								out_fh.write(delim.join(line[:5]) + '\n')

def findMotifsGenome_genomic_background(peak_files, size_values):
	for peak_file in peak_files:
		for size_req in size_values:
			out_dir = peak_file.rsplit('.', 1)[0] + '.size_' + size_req + '.find_motifs'
			#findMotifsGenome.pl <peak/BED file> <genome> <output directory> -size # [options]
			find_motifs = subprocess.Popen(['findMotifsGenome.pl', peak_file, 'hg38', out_dir, '-size', size_req, '-preparsedDir', './temp4'])
			find_motifs.wait()

def findMotifsGenome_human_org_comp(peak_pair_dict, size_values):
	for peak_pair in peak_pair_dict:
		for size_req in size_values:
			peak_file = peak_pair_dict[peak_pair][0]
			bg_peak_file = peak_pair_dict[peak_pair][1]
			out_dir = 'diff_peaks_rep.bg.human_vs_org_rods_cones.size_' + size_req + '.' + peak_file.split('.')[4] + '.' + peak_file.split('.')[5] + '.find_motifs'
			#findMotifsGenome.pl <peak/BED file> <genome> <output directory> -size # [options]
			find_motifs = subprocess.Popen(['findMotifsGenome.pl', peak_file, 'hg38', out_dir, '-size', size_req, '-bg', bg_peak_file, '-preparsedDir', './temp4'])
			find_motifs.wait()

def findMotifsGenome_org_human_comp(peak_pair_dict, size_values):
	for peak_pair in peak_pair_dict:
		for size_req in size_values:
			peak_file = peak_pair_dict[peak_pair][1]
			bg_peak_file = peak_pair_dict[peak_pair][0]
			out_dir = 'diff_peaks_rep.bg.org_vs_human_rods_cones.size_' + size_req + '.' + peak_file.split('.')[4] + '.' + peak_file.split('.')[5] + '.find_motifs'
			#findMotifsGenome.pl <peak/BED file> <genome> <output directory> -size # [options]
			find_motifs = subprocess.Popen(['findMotifsGenome.pl', peak_file, 'hg38', out_dir, '-size', size_req, '-bg', bg_peak_file, '-preparsedDir', './temp4'])
			find_motifs.wait()


##run methods

##step1. pairwise comparisons
cell_class_dict = {'cones':['human.Mature_Cones', 'organoid.Cones'], 'rods':['human.Mature_Rods', 'organoid.Rods'],
		'bipolars':['human.Mature_Bipolars', 'organoid.Bipolar_Cells'], 'mullers':['human.Mature_Mullers', 'organoid.Muller_Glia'],
		'horizontals':['human.Mature_Horizontals', 'organoid.Amacrine_Horizontal_Cells'], 'amacrines':['human.Mature_Amacrines', 'organoid.Amacrine_Horizontal_Cells'],
		'early_progenitor':['human.Early_Progenitors', 'organoid.Early_RPCs'], 'late_progenitor':['human.Late_Progenitors', 'organoid.Late_RPCs']}
d_values_wanted = ['100', '500']
size_values_wanted = ['500', '2000'] ##does none work for homer
bed_suffixes = ['.macs2_q0.01_summits.bed','.macs2_q0.000001_summits.bed']
tag_dir_suffix = '.tag_dir'
get_correlation_info_homer_pairwise(cell_class_dict, d_values_wanted, size_values_wanted, bed_suffixes, tag_dir_suffix)
##rpt for none size
d_values_wanted = ['100']
size_values_wanted = ['none']
get_correlation_info_homer_pairwise(cell_class_dict, d_values_wanted, size_values_wanted, bed_suffixes, tag_dir_suffix)


##step 2. compare human/organoid data using homer getDifferentialPeaks
human_ccs = ['human.AC_HC_GC_Precursors', 'human.Developing_Amacrines', 'human.Developing_Bipolars', 'human.Developing_Cones', 
	'human.Developing_Ganglions', 'human.Developing_Horizontals', 'human.Developing_Rods', 'human.Early_Progenitors', 
	'human.Ganglion_Precursors', 'human.Late_Progenitors', 'human.Mature_Amacrines', 'human.Mature_Bipolars', 
	'human.Mature_Cones', 'human.Mature_Horizontals', 'human.Mature_Mullers', 'human.Mature_Rods', 'human.Photoreceptor_Bipolar_Precursors']
org_ccs =['organoid.AC_HC_Precursors', 'organoid.Cones', 'organoid.Early_RPCs', 'organoid.Late_RPCs', 'organoid.PR_BC_Precursors', 
	'organoid.RGCs', 'organoid.Rods', 'organoid.Amacrine_Horizontal_Cells', 'organoid.Bipolar_Cells', 'organoid.Muller_Glia', 'organoid.Developing_RGCs']
human_combined = 'human_combined'
org_combined = 'organoid_combined'
q_values_req = ['0.01', '0.000001']
d_values_req = ['100', '500']
size_values_wanted = ['500', '2000']
p_values_req = ['0.01', '0.001']
fc_values_req = ['3', '2']
diff_peak_prefix_human_org = 'diff_peaks.human_organoid.'
diff_peak_prefix_org_human = 'diff_peaks.organoid_human.'
##combine individul tag dir for human/org
combine_tag_dirs(human_ccs, human_combined)
combine_tag_dirs(human_ccs, org_combined)
##run homer getDifferentialPeaks
get_diff_peaks_diff_sig(diff_peak_prefix_human_org, diff_peak_prefix_org_human, d_values_req, q_values_req, p_values_req, fc_values_req, size_values_wanted, human_combined, org_combined)


##step 3. compare adult human/organoid data using homer getDifferentialPeaksReplicates
human_samples = ['hu5_bulk', 'hu7_bulk', 'hu8_bulk']
org_samples = ['28-1_bulk', '28-2_bulk']
diff_peak_rep_bg_prefix_human_org = 'diff_peaks_rep.bg.adult_human_28wk_organoid.'
diff_peak_rep_bg_prefix_org_human = 'diff_peaks_rep.bg.28wk_organoid_adult_human.'
diff_peak_rep_bg_prefix_human_org_all = 'diff_peaks_rep.bg.adult_human_28wk_organoid.all_peaks.'
diff_peak_rep_bg_prefix_org_human_all = 'diff_peaks_rep.bg.28wk_organoid_adult_human.all_peaks.'
get_diff_peaks_replicates(diff_peak_rep_bg_prefix_human_org, diff_peak_rep_bg_prefix_org_human, human_samples, org_samples, 'no')
get_diff_peaks_replicates(diff_peak_rep_bg_prefix_human_org_all, diff_peak_rep_bg_prefix_org_human_all, human_samples, org_samples, 'yes')


##step 4. compare results from 10 against cell classes
mature_ccs = ['human.Mature_Amacrines', 'human.Mature_Bipolars', 'human.Mature_Cones', 'human.Mature_Horizontals', 
	'human.Mature_Mullers', 'human.Mature_Rods', 'organoid.Cones', 'organoid.RGCs', 'organoid.Rods', 
	'organoid.Amacrine_Horizontal_Cells', 'organoid.Bipolar_Cells', 'organoid.Muller_Glia']
homer_peak_suffix = '.homer_peaks.default.txt'
homer_peak_bed_suffix = '.homer_peaks.default.bed'
diff_peak_rep_files = ['diff_peaks_rep.28wk_organoid_adult_human.all_peaks.balanced.txt', 'diff_peaks_rep.28wk_organoid_adult_human.balanced.txt', 
		'diff_peaks_rep.adult_human_28wk_organoid.all_peaks.balanced.txt', 'diff_peaks_rep.adult_human_28wk_organoid.balanced.txt', 
		'diff_peaks_rep.bg.28wk_organoid_adult_human.all_peaks.balanced.txt', 'diff_peaks_rep.bg.28wk_organoid_adult_human.balanced.txt', 
		'diff_peaks_rep.bg.adult_human_28wk_organoid.all_peaks.balanced.txt', 'diff_peaks_rep.bg.adult_human_28wk_organoid.balanced.txt']
diff_peak_rep_beds = [i.rsplit('.', 1)[0] + '.bed' for i in diff_peak_rep_files]
##get homer peaks for cell classes
find_peaks_cell_classes(mature_ccs, homer_peak_suffix)
##make bed files from homer getDifferentialPeaksReplicates results and cell class findpeaks
make_beds_from_diff_peaks_file(diff_peak_rep_files)
make_beds_from_find_peaks_files(mature_ccs, homer_peak_suffix)
##filter background files by fc and adj pvalue and make bed
files_to_filter = ['diff_peaks_rep.bg.28wk_organoid_adult_human.balanced.txt', 'diff_peaks_rep.bg.adult_human_28wk_organoid.balanced.txt']
#logfc/p value pairs, first is default
params_to_filter = [[1.0, 0.05], [2.0, 0.001], [3.0, 0.001]]
filter_and_make_beds_from_diff_peaks_file(files_to_filter, params_to_filter)
##compare cell class peaks with getDifferentialPeaksReplicates results
diff_peak_rep_beds_cutoffs = ['diff_peaks_rep.bg.28wk_organoid_adult_human.balanced.lfc1_p0.05.bed', 'diff_peaks_rep.bg.28wk_organoid_adult_human.balanced.lfc2_p0.001.bed', 
		'diff_peaks_rep.bg.28wk_organoid_adult_human.balanced.lfc3_p0.001.bed', 'diff_peaks_rep.bg.adult_human_28wk_organoid.balanced.lfc1_p0.05.bed', 
		'diff_peaks_rep.bg.adult_human_28wk_organoid.balanced.lfc2_p0.001.bed', 'diff_peaks_rep.bg.adult_human_28wk_organoid.balanced.lfc3_p0.001.bed']
bt_int_homer_beds_vs_cc_beds(diff_peak_rep_beds_cutoffs, mature_ccs, homer_peak_bed_suffix)
##reformat for graphing
files_to_reformat = ['diff_peaks_rep.bg.28wk_organoid_adult_human.balanced.lfc1_p0.05.bt_int.mature_cc.txt', 'diff_peaks_rep.bg.28wk_organoid_adult_human.balanced.lfc2_p0.001.bt_int.mature_cc.txt', 'diff_peaks_rep.bg.28wk_organoid_adult_human.balanced.lfc3_p0.001.bt_int.mature_cc.txt', 'diff_peaks_rep.bg.adult_human_28wk_organoid.balanced.lfc1_p0.05.bt_int.mature_cc.txt', 'diff_peaks_rep.bg.adult_human_28wk_organoid.balanced.lfc2_p0.001.bt_int.mature_cc.txt', 'diff_peaks_rep.bg.adult_human_28wk_organoid.balanced.lfc3_p0.001.bt_int.mature_cc.txt', 'diff_peaks_rep.bg.28wk_organoid_adult_human.all_peaks.balanced.bt_int.mature_cc.txt', 'diff_peaks_rep.bg.adult_human_28wk_organoid.all_peaks.balanced.bt_int.mature_cc.txt']
format_int_files_for_upset(files_to_reformat)


##step 5. continuation of step11, combined peaks files, filter and then find motifs in diff peaks
intersected_files_org_human = ['diff_peaks_rep.bg.28wk_organoid_adult_human.balanced.lfc1_p0.05.bt_int.mature_cc.upset.wide.txt', 
		'diff_peaks_rep.bg.28wk_organoid_adult_human.balanced.lfc2_p0.001.bt_int.mature_cc.upset.wide.txt', 
		'diff_peaks_rep.bg.28wk_organoid_adult_human.balanced.lfc3_p0.001.bt_int.mature_cc.upset.wide.txt']
intersected_files_human_org = ['diff_peaks_rep.bg.adult_human_28wk_organoid.balanced.lfc1_p0.05.bt_int.mature_cc.upset.wide.txt', 
		'diff_peaks_rep.bg.adult_human_28wk_organoid.balanced.lfc2_p0.001.bt_int.mature_cc.upset.wide.txt', 
		'diff_peaks_rep.bg.adult_human_28wk_organoid.balanced.lfc3_p0.001.bt_int.mature_cc.upset.wide.txt']
homer_ann_org_human = 'diff_peaks_rep.bg.28wk_organoid_adult_human.balanced.txt'
homer_ann_human_org = 'diff_peaks_rep.bg.adult_human_28wk_organoid.balanced.txt'
peak_info_files_org_human = [i.rsplit('.', 5)[0] + '.peak_info.txt' for i in intersected_files_org_human]
peak_info_files_human_org = [i.rsplit('.', 5)[0] + '.peak_info.txt' for i in intersected_files_human_org]
homer_peaks_human_org = ['diff_peaks_rep.bg.adult_human_28wk_organoid.balanced.lfc1_p0.05.human_rods_cones.homer_peak', 
	'diff_peaks_rep.bg.adult_human_28wk_organoid.balanced.lfc1_p0.05.human_rods.homer_peak', 
	'diff_peaks_rep.bg.adult_human_28wk_organoid.balanced.lfc2_p0.001.human_rods_cones.homer_peak', 
	'diff_peaks_rep.bg.adult_human_28wk_organoid.balanced.lfc2_p0.001.human_rods.homer_peak', 
	'diff_peaks_rep.bg.adult_human_28wk_organoid.balanced.lfc3_p0.001.human_rods_cones.homer_peak', 
	'diff_peaks_rep.bg.adult_human_28wk_organoid.balanced.lfc3_p0.001.human_rods.homer_peak']
findmotifs_sizes = ['200', 'given']
homer_peaks_human_org_dict = {'lfc1_p0.05': ['diff_peaks_rep.bg.adult_human_28wk_organoid.balanced.lfc1_p0.05.human_rods_cones.homer_peak', 'diff_peaks_rep.bg.28wk_organoid_adult_human.balanced.lfc1_p0.05.org_rods_cones.homer_peak'],
		'lfc2_p0.001': ['diff_peaks_rep.bg.adult_human_28wk_organoid.balanced.lfc2_p0.001.human_rods_cones.homer_peak', 'diff_peaks_rep.bg.28wk_organoid_adult_human.balanced.lfc2_p0.001.org_rods_cones.homer_peak'],
		'lfc3_p0.001': ['diff_peaks_rep.bg.adult_human_28wk_organoid.balanced.lfc3_p0.001.human_rods_cones.homer_peak', 'diff_peaks_rep.bg.28wk_organoid_adult_human.balanced.lfc3_p0.001.org_rods_cones.homer_peak']}
##combine the annotated peak info from homer with the bt int data (from the filtered peaks)
combine_peak_info(intersected_files_org_human, homer_ann_org_human)
combine_peak_info(intersected_files_human_org, homer_ann_human_org)
##filter peak files for those specific to photorecepters i.e. 3 ways human rods, human rods/cones, org rods/cones
filter_peak_file_prs(peak_info_files_org_human, 'org_rods_cones')
filter_peak_file_prs(peak_info_files_human_org, 'human_rods')
filter_peak_file_prs(peak_info_files_human_org, 'human_rods_cones')
##run findMotifsGenome a few ways
##compare against genomic background for the human data
findMotifsGenome_genomic_background(homer_peaks_human_org, findmotifs_sizes)
##comapare human rods/cones vs org rods/cones
findMotifsGenome_human_org_comp(homer_peaks_human_org_dict, findmotifs_sizes)
##comapare org rods/cones vs human rods/cones
findMotifsGenome_org_human_comp(homer_peaks_human_org_dict, findmotifs_sizes)



