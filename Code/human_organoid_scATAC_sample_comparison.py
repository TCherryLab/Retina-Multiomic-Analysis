#!/usr/bin/env python
import subprocess
import os
import glob

##parameters
delim = '\t'
thread_number = '20'


##methods
def call_macs2_on_original_bams(sample_list):
	print(sample_list)
	##run macs2
	for sample in sample_list:
		bam = sample + '.possorted.bam'
		name = sample + '_bulk'
		##outfile names
		out_bampe_keepdups_2 =  name + '.macs2.bampe_p1e-5_keepdups'
		out_bampe_keepdups_3 =  name + '.macs2.bampe_p1e-10_keepdups'
		out_bampe_keepdups_4 =  name + '.macs2.bampe_q0.01_keepdups'
		##run macs2 from bam files
		run_macs2_2 = subprocess.Popen(['macs2', 'callpeak', '-t', bam, '-f', 'BAMPE', '-n', out_bampe_keepdups_2, '-g', 'hs', '-p', '1e-5', '--keep-dup', 'all'])
		run_macs2_2.wait()		
		run_macs2_2 = subprocess.Popen(['macs2', 'callpeak', '-t', bam, '-f', 'BAMPE', '-n', out_bampe_keepdups_3, '-g', 'hs', '-p', '1e-10', '--keep-dup', 'all', '--tempdir', '.'])
		run_macs2_2.wait()
		run_macs2_2 = subprocess.Popen(['macs2', 'callpeak', '-t', bam, '-f', 'BAMPE', '-n', out_bampe_keepdups_4, '-g', 'hs', '-q', '0.01', '--keep-dup', 'all', '--tempdir', '.'])
		run_macs2_2.wait()

def make_tag_dirs(name, bam):
	outdir = name + '.tag_dir'
	mk_tag_dir = subprocess.Popen(['makeTagDirectory', outdir, bam])
	mk_tag_dir.wait()

def merge_peaks_original_bams(out_prefix, d_value, peak_files, out_suffix):
	out_file = out_prefix + d_value + 'd' + out_suffix
	print(len(peak_files), out_file)
	with open(out_file, 'w') as out_fh:
		# mk_tag_dir = subprocess.Popen(['mergePeaks', '-d', d_value, '-matrix', out_prefix + d_value] + peak_files, stdout=out_fh)
		mk_tag_dir = subprocess.Popen(['mergePeaks', '-d', d_value] + peak_files, stdout=out_fh)
		mk_tag_dir.wait()

def annotate_peaks_original_bams(samples, peak_prefix, peak_suffix, d_value, tag_suffix, out_prefix, sizes_req):
	peak_file = peak_prefix + d_value + 'd' + peak_suffix
	tag_dirs = []
	##get list of tag_dirs
	for sample in samples:
		tag_dir = sample + '_bulk' + tag_suffix
		tag_dirs.append(tag_dir)
	out_file_no_size = out_prefix + d_value + 'd' + peak_suffix
	print(len(tag_dirs), peak_file, out_file_no_size)
	with open(out_file_no_size, 'w') as out_ns_fh:
		mk_tag_dir = subprocess.Popen(['annotatePeaks.pl', peak_file, 'hg38', '-d'] + tag_dirs, stdout=out_ns_fh)
		mk_tag_dir.wait()
	for size_req in sizes_req:
		out_file_using_size = out_prefix + d_value + 'd.' + size_req + 'size' + peak_suffix
		with open(out_file_using_size, 'w') as out_fh:
			mk_tag_dir = subprocess.Popen(['annotatePeaks.pl', peak_file, 'hg38', '-size', size_req, '-d'] + tag_dirs, stdout=out_fh)
			mk_tag_dir.wait()

def compute_correlation_original_bams(annotate_prefix, annotate_suffix, d_value, sizes_req, out_prefix):
	##get list of infile
	infiles = []
	no_size_file = annotate_prefix + d_value + 'd' + annotate_suffix
	infiles.append(no_size_file)
	for size_req in sizes_req:
		out_file_using_size = annotate_prefix + d_value + 'd.' + size_req + 'size' + annotate_suffix
		infiles.append(out_file_using_size)
	print(len(infiles), infiles)
	##format file by file
	for infile in infiles:
		out_file = out_prefix  + infile.split('.', 1)[1]
		print(infile, out_file)
		with open(out_file, 'w') as out_fh, open(infile, 'r') as in_fh:
			lc = 0
			for line in in_fh:
				lc += 1
				line = line.strip('\n').split(delim)
				if lc == 1:
					sample_info = line[19:]
					sample_info = [s.split('.')[0] for s in sample_info]
					out_fh.write(delim.join(['peak_id'] + sample_info + ['\n']))
					print(sample_info)
				else:
					line_out = [line[0]] + line[19:]
					out_fh.write(delim.join(line_out + ['\n']))

def homer_correlation_on_original_bams(samples, d_values, size_values, bed_suffixes, tag_dirs_needed):
	##make tag dirs
	if tag_dirs_needed == 'yes':
		for sample in samples:
			bam = sample + '.possorted.bam'
			name = sample + '_bulk'
			##make tag dirs so can analyze
			make_tag_dirs(name, bam)
	##for different parameters merge peaks and then annoate
	tag_dir_suffix = '.tag_dir'
	merge_prefix = 'merged.'
	annotate_prefix = 'annotated.'
	correlation_prefix = 'correlation.'
	for d in d_values:
		for bed_suffix in bed_suffixes:
			##merge peak files for annotating the peaks
			peak_beds = [s + '_bulk' + bed_suffix  for s in samples]
			merge_file_suffix = bed_suffix.rsplit('.', 1)[0] + '.txt'
			merge_peaks_original_bams(merge_prefix, d, peak_beds, merge_file_suffix)
			##annotate those peaks, and then make a correlation file
			annotate_peaks_original_bams(samples, merge_prefix, merge_file_suffix, d, tag_dir_suffix, annotate_prefix, size_values)
			##make file to use in r for heatmaps etc
			compute_correlation_original_bams(annotate_prefix, merge_file_suffix, d, size_values, correlation_prefix)

##run methods

##correlation on original bams 
##sample names
combined_sample_list = ['d53', 'd59', 'd74', 'd78', 'd113', 'd132', 'hu5', 'hu7', 'hu8', 
		'IPSC_c4', 'IPSC_c5_1', '5wk', '5wk_c5_1', '20wk', '20wk_c5_1', '28-1', '28-2',
		'12wk1', '12wk2', '12wk3']
human_sample_list = ['d53', 'd59', 'd74', 'd78', 'd113', 'd132', 'hu5', 'hu7', 'hu8']
org_sample_list = ['IPSC_c4', 'IPSC_c5_1', '5wk', '5wk_c5_1', '20wk', '20wk_c5_1', 
		'28-1', '28-2', '12wk1', '12wk2', '12wk3']
##homer values
d_values = ['100', '200', '500']
size_values = ['500', '2000']
size_values_plus = ['500', '2000', 'none']
peak_bed_suffices = ['.macs2.bampe_p1e-10_keepdups_summits.bed', '.macs2.bampe_q0.01_keepdups_summits.bed']

##call peaks using macs2
call_macs2_on_original_bams(combined_sample_list)

##use homer to look for correlation
homer_correlation_on_original_bams(combined_sample_list, d_values, size_values, peak_bed_suffices, 'yes')





