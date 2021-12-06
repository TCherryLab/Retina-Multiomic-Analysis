#!/usr/bin/env python
import subprocess
import os
import glob


##parameters
delim = '\t'
thread_number = '20'

##programs
sinto = 'sinto'


##methods
def run_sinto_method(sample_name, bam_file, cells_file, w_dir):
	directory = './' + sample_name
	if not os.path.exists(sample_name):
		os.makedirs(sample_name)
	os.chdir(directory)
	run_sinto = subprocess.Popen([sinto, 'filterbarcodes', '-b', w_dir + bam_file,  '-c', w_dir + cells_file, '-p', thread_number])
	run_sinto.wait()
	os.chdir(w_dir)

def get_cell_info_ipscs(infile, outfile_suffix):
	##make dict of all cell barcode by sample
	cell_info_dict = {}
	with open(infile, "r") as in_fh:
		lc = 0
		for line in in_fh:
			lc += 1
			if lc > 1:
				line = line.replace('"', '').rstrip().split(',')
				sample = line[1].split('#')[0]
				cell_barcode = line[1].split('#')[1]
				# print(sample, cell_barcode)
				if sample in cell_info_dict:
					cell_info_dict[sample].append(cell_barcode)
				else:
					cell_info_dict[sample] = [cell_barcode]
	##print out just for ipscs with ipsc as cell type
	out_files = []
	for s in cell_info_dict:
		print(s, len(cell_info_dict[s]))
		if s.startswith('i'):
			outfile = s + outfile_suffix
			out_files.append(outfile)
			with open(outfile, "w") as out_fh:
				for bc in cell_info_dict[s]:
					out_fh.write(delim.join([bc, 'IPSCs']) + '\n')
	return(out_files)

def get_cell_info_org_human(infile, outfile_suffix):
	##make dict of all cell barcode by sample
	cell_info_dict = {}
	with open(infile, "r") as in_fh:
		lc = 0
		for line in in_fh:
			lc += 1
			if lc > 1:
				line = line.replace('"', '').rstrip().split(',')
				# print(line)
				sample = line[1].split('#')[0]
				cell_barcode = line[1].split('#')[1]
				cell_class = line[2].replace(' ', '_').replace('/', '_')
				# print(sample, cell_barcode)
				if sample in cell_info_dict:
					cell_info_dict[sample].append(cell_barcode + ',' +  cell_class)
				else:
					cell_info_dict[sample] = [cell_barcode + ',' +  cell_class]
	##print out just for ipscs with ipsc as cell type
	out_files = []
	for s in cell_info_dict:
		print(s, len(cell_info_dict[s]))
		outfile = s + outfile_suffix
		out_files.append(outfile)
		with open(outfile, "w") as out_fh:
			for bc in cell_info_dict[s]:
				cell_info = bc.split(',')
				if len(cell_info) != 2:
					print(cell_info, 'looks odd')
				out_fh.write(delim.join(cell_info) + '\n')
	return(out_files)

def run_macs2_on_all(samples, qs_req):
	##run macs2
	for sample in samples:
		for q_req in qs_req:
			bams = glob.glob(sample + '/' + '*bam')
			for bam in bams:
				outname = sample + '.' + bam.split('/')[1].split('.')[0] + '.macs2_q' + q_req
				print(sample, bam, outname)
				##run macs2
				run_macs2 = subprocess.Popen(['macs2', 'callpeak', '-t', bam, '-f', 'BAMPE', '-n', outname, '-g', 'hs', '-q', '0.01', '--keep-dup', 'all'])
				run_macs2.wait()


def run_macs2_cell_class_collapsed(infile, qs_req):
	for q_req in qs_req:
		macs2_dict = {}
		with open(infile, "r") as in_fh:
			lc = 0
			for line in in_fh:
				lc += 1
				if lc > 1:
					line = line.rstrip().split(delim)
					name = '.'.join([line[0], line[2]]) + '.macs2_q' + q_req
					bam = line[3]
					if name in macs2_dict:
						macs2_dict[name].append(bam)
					else:
						macs2_dict[name] = [bam]
		for outname in macs2_dict:
			# bams = ' '.join(macs2_dict[outname])
			bams = macs2_dict[outname]
			print(outname, bams)
			##run macs2
			# run_macs2 = subprocess.Popen(['macs2', 'callpeak', '-t', bams, '-f', 'BAMPE', '-n', outname, '-g', 'hs', '-q', '0.01', '--keep-dup', 'all'])
			run_macs2 = subprocess.Popen(['macs2', 'callpeak', '-t'] + bams + ['-f', 'BAMPE', '-n', outname, '-g', 'hs', '-q', '0.01', '--keep-dup', 'all'])
			run_macs2.wait()

def make_tag_dirs(name, bam):
	outdir = name + '.tag_dir'
	mk_tag_dir = subprocess.Popen(['makeTagDirectory', outdir, bam])
	mk_tag_dir.wait()

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

	for size_req in sizes_req:
		if size_req == 'none':
			out_file_no_size = out_prefix + d_value + 'd.' + size_req + 'size' + peak_suffix
			with open(out_file_no_size, 'w') as out_ns_fh:
				mk_tag_dir = subprocess.Popen(['annotatePeaks.pl', peak_file, 'hg38', '-d'] + tag_dirs, stdout=out_ns_fh)
				mk_tag_dir.wait()
		else:
			out_file_using_size = out_prefix + d_value + 'd.' + size_req + 'size' + peak_suffix
			with open(out_file_using_size, 'w') as out_fh:
				mk_tag_dir = subprocess.Popen(['annotatePeaks.pl', peak_file, 'hg38', '-size', size_req, '-d'] + tag_dirs, stdout=out_fh)
				mk_tag_dir.wait()

def get_correlation_info_homer(infile, d_values, size_values, bed_suffix):
	##read in info file
	bam_dict = {}
	ind_beds, comb_beds = [], []
	with open(infile, "r") as in_fh:
		lc = 0
		for line in in_fh:
			lc += 1
			if lc > 1:
				line = line.rstrip().split(delim)
				name = '.'.join([line[1], line[2]])
				combined_bed = '.'.join([line[0], line[2]]) + bed_suffix
				ind_beds.append(name + bed_suffix)
				if combined_bed not in comb_beds:
					comb_beds.append(combined_bed)
				bam = line[3]
				if name in bam_dict:
					print(name, bam, 'name seen mutiple times')
				else:
					bam_dict[name] = bam
	##make tag dirs
	for td_name in bam_dict:
		bam = bam_dict[td_name]
		print(td_name, bam)
		##make tag dirs so can analyze
		make_tag_dirs(td_name, bam)
	##combine bed files
	print(len(ind_beds), len(comb_beds))
	for d in d_values:
		ind_merge_prefix = 'merged_ind_samples.'
		comb_merge_prefix = 'merged_combined.'
		merge_peaks(ind_merge_prefix, d, ind_beds, bed_suffix)
		merge_peaks(comb_merge_prefix, d, comb_beds, bed_suffix)
	##annotate files
	samples = bam_dict.keys()
	tag_dir_suffix = '.tag_dir'
	ind_annotate_prefix = 'annotated.ind_sample_beds.'
	comb_annotate_prefix = 'annotated.combined_beds.'
	# correlation_prefix = 'correlation.'
	for d in d_values:
		ind_merge_bed = 'merged_ind_samples.' + d + 'd' + bed_suffix
		comb_merge_bed = 'merged_combined.' + d + 'd' + bed_suffix
		outfile_suffix = bed_suffix.rsplit('.',1)[0] + '.txt'
		##annotate those peaks, and then make a correlation file
		annotate_peaks(samples, ind_merge_bed , d, tag_dir_suffix, ind_annotate_prefix, size_values, outfile_suffix)
		annotate_peaks(samples, comb_merge_bed , d, tag_dir_suffix, comb_annotate_prefix, size_values, outfile_suffix)

def collapse_samples_by_cell_class(infile):
	coll_dict = {}
	with open(infile, "r") as in_fh:
		lc = 0
		for line in in_fh:
			lc += 1
			if lc > 1:
				line = line.rstrip().split(delim)
				cell_class = line[2]
				sample = line[1]
				bam_file = sample + '/' + cell_class + '.bam'
				type_cc = line[0] + '.' + cell_class
				if type_cc in coll_dict:
					coll_dict[type_cc].append(bam_file)
				else:
					coll_dict[type_cc] = [bam_file]
	return(coll_dict)

def merge_bams(collapsed_dict, w_dir):
	for sample in collapsed_dict:
		merged_bam = sample + '.bam'
		# input_bams = ['I=' + b for b in collapsed_dict[sample]]
		input_bams = collapsed_dict[sample]
		print(merged_bam, input_bams)
		##can use picard or samtools to do this
		# picard_md = subprocess.Popen(['picard', 'MergeSamFiles', 'CREATE_INDEX=true', 'VALIDATION_STRINGENCY=SILENT'] + input_bams + ['O=' + merged_bam, 'TMP_DIR=' + w_dir])
		# picard_md.wait()
		st_merge = subprocess.Popen(['samtools', 'merge', '-O', 'bam', '-@', '16', merged_bam] + input_bams)
		st_merge.wait()
		st_index = subprocess.Popen(['samtools', 'index', merged_bam])
		st_index.wait()



##run methods

##step1. make cell info files for sinto
ipsc_raw_info = 'org_all_w_ipscs.092420.csv'
org_raw_info = 'org_all.cell_class_info.092420.csv'
human_raw_info = 'human_all.cell_class_info.092520.csv'
cell_info_file_suffix = '.cell_types.txt'
##get files x3
ipsc_cell_info_files = get_cell_info_ipscs(ipsc_raw_info, cell_info_file_suffix)
print(ipsc_cell_info_files)
org_cell_info_files = get_cell_info_org_human(org_raw_info, cell_info_file_suffix)
print(org_cell_info_files)
human_cell_info_files = get_cell_info_org_human(human_raw_info, cell_info_file_suffix)
print(human_cell_info_files)
all_cell_info_files = ipsc_cell_info_files + org_cell_info_files + human_cell_info_files

##step2. run sinto
bam_dict = {'20wk1': '20wk.possorted.bam', '20wk2': '20wk_c5_1.possorted.bam', '28wk1': '28-1.possorted.bam', 
		'28wk2': '28-2.possorted.bam', '5wk1': '5wk.possorted.bam', '5wk2': '5wk_c5_1.possorted.bam', 
		'12wk1': '12wk1.possorted.bam', '12wk2': '12wk2.possorted.bam', '12wk3': '12wk3.possorted.bam',
		'd113': 'd113.possorted.bam', 'd132': 'd132.possorted.bam', 'd53': 'd53.possorted.bam', 
		'd59': 'd59.possorted.bam', 'd74': 'd74.possorted.bam', 'd78': 'd78.possorted.bam', 
		'Hu5': 'hu5.possorted.bam', 'Hu7': 'hu7.possorted.bam', 'Hu8': 'hu8.possorted.bam', 
		'ipsc1': 'IPSC_c4.possorted.bam', 'ipsc2': 'IPSC_c5_1.possorted.bam'}
for cell_info_file in all_cell_info_files:
	sample = cell_info_file.split('.')[0]
	bam = bam_dict[sample]
	# print(sample, bam, cell_info_file)
	run_sinto_method(sample, bam, cell_info_file, working_dir)


##step3. merge bams to get cll class specific bams
##cell class info file made from 'cell_types.txt' files and remove any cell class with less than 20 cells per sample
sample_cell_class_info = 'human_org_samples_cell_classes_gt20.txt'
##get info
coll_cell_class_dict = collapse_samples_by_cell_class(sample_cell_class_info)
print(coll_cell_class_dict)
##merge bams
merge_bams(coll_cell_class_dict, working_dir)

##step4. run macs2 many ways
samples = bam_dict.keys()
bam_file_combo_info = 'bam_file_and_cell_class_info.txt'
q_values = ['0.01', '0.000001']
##on all samples and cell classes individually
run_macs2_on_all(samples. q_values)
##on all cell classes
run_macs2_cell_class_collapsed(bam_file_combo_info, q_values)

##step5. run homer i.e. make tag dirs, combine beds and annotate
d_values_wanted = ['100', '500']
size_values_wanted = ['500', '2000'] ##does none work for homer
bed_suffix = '.macs2_q0.01_summits.bed'
get_correlation_info_homer(bam_file_combo_info, d_values_wanted, size_values_wanted, bed_suffix)









