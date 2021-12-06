#!/usr/bin/env python
import subprocess

##parameters
delim = '\t'

##methods
def merge_peaks(out_prefix, d_value, peak_files, out_suffix):
	out_file = out_prefix + d_value + 'd' + out_suffix
	print(len(peak_files), out_file)
	with open(out_file, 'w') as out_fh:
		# mk_tag_dir = subprocess.Popen(['mergePeaks', '-d', d_value, '-matrix', out_prefix + d_value] + peak_files, stdout=out_fh)
		run_merge_peaks = subprocess.Popen(['mergePeaks', '-d', d_value] + peak_files, stdout=out_fh)
		run_merge_peaks.wait()

def annotate_peaks(samples, peak_file, d_value, tag_suffix, out_prefix, sizes_req, peak_suffix):
	tag_dirs = []
	##get list of tag_dirs
	for sample in samples:
		tag_dir = sample + tag_suffix
		tag_dirs.append(tag_dir)
	print(len(tag_dirs), peak_file)
	##annotate
	for size_req in sizes_req:
		out_file_using_size = out_prefix + d_value + 'd.' + size_req + 'size' + peak_suffix
		with open(out_file_using_size, 'w') as out_fh:
			mk_tag_dir = subprocess.Popen(['annotatePeaks.pl', peak_file, 'hg38', '-size', size_req, '-d'] + tag_dirs, stdout=out_fh)
			mk_tag_dir.wait()

def get_correlation_info_homer(cc_dict, d_values, size_values, bed_suffix, tag_dir_suffix):
	for cell_class in cc_dict:
		samples = cc_dict[cell_class]
		##combine bed files
		merge_prefix = cell_class + '_merged.'
		annotate_prefix = cell_class + '_annotated.'
		annotate_suffix = bed_suffix.rsplit('.',1)[0] + '.txt'
		corr_prefix = cell_class + '_correlation.'
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



##run methods
cell_class_dict = {'cones':['human.Mature_Cones', 'organoid.Cones'], 'rods':['human.Mature_Rods', 'organoid.Rods'],
		'bipolars':['human.Mature_Bipolars', 'organoid.Bipolar_Cells'], 'mullers':['human.Mature_Mullers', 'organoid.Muller_Glia'],
		'horizontals':['human.Mature_Horizontals', 'organoid.Horizontal_Cells'], 'amacrines':['human.Mature_Amacrines', 'organoid.Amacrine_Cells'],
		'early_progenitor':['human.Early_Progenitors', 'organoid.Early_Progenitors'], 'late_progenitor':['human.Late_Progenitors', 'organoid.Late_Progenitors']}

d_values_wanted = ['100', '500']
size_values_wanted = ['500', '2000'] ##does none work for homer
bed_suffix = '.macs2_q0.01_summits.bed'
tag_dir_suffix = '.tag_dir'
get_correlation_info_homer(cell_class_dict, d_values_wanted, size_values_wanted, bed_suffix, tag_dir_suffix)

