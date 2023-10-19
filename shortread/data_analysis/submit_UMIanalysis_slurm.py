#!/usr/bin/python
from grp import getgrnam, getgrgid
import argparse
import os
import re
import subprocess

dirname =  os.path.dirname(os.path.abspath(__file__))

parser = argparse.ArgumentParser(description='Submit UMI analysis jobs to the UGent HPC (slurm)')

# Read arguments
parser.add_argument('--version', action='version')
parser.add_argument('-m', nargs=1, choices=['fastq'], required=True, help='Specify the format of the input files (bam currently not supported)')
parser.add_argument('-t', nargs=1, choices=['se','pe'], required=True, help='Single end (se) or paired end (pe) data (SE currently not supported)')
parser.add_argument('-c', nargs=1, choices=['yes', 'no'], required=True, help='Perform consensus calling with calib? For SE data, always choose "no"')
parser.add_argument('-a', nargs=1, choices=['RNAexome', 'picoV3'], required=True, help='Analysis type')
parser.add_argument('-b', nargs=1, required=True, help='The base directory where the sample subdirectories are located', metavar='base_dir')
parser.add_argument('-g', nargs=1, choices=['GRCh38', 'GRCm39'], required=True, help='Genome build (GRCh38, GRCm39)')
parser.add_argument('-l', nargs=1, required=True, help='A text file with UMI sequence on each line. Use the keyword \'none\' for picoV3')
parser.add_argument('-n', nargs=1, choices=['50', '75', '100','150'], required=True, help='Length of reads (50, 75, 100 or 150 nts)')
# parser.add_argument('--stranded', nargs=1, choices=['yes','no'], required=True, help='(reverse) stranded sequencing or not (unstranded)?')
parser.add_argument('-o', nargs=1, required=True, help='The directory where output should be created', metavar='output_dir')
parser.add_argument('-s', nargs=1, required=True, help='String match for sample directories e.g. RNA0', metavar='string_match')
parser.add_argument('-u', nargs=1, required=True, help='The submitters email address (used for error reporting)', metavar="user_email")

# Parse arguments
args = parser.parse_args()
mode = args.m[0]
data_type = args.t[0]
ccalling = args.c[0]
analysis_type = args.a[0]
base_dir = args.b[0]
genome_build = args.g[0]
umi_file = args.l[0]
read_length = int(args.n[0])
string_match = args.s[0]
email = args.u[0]
# stranded = args.stranded[0]

# # HUMAN files
ref_fasta = "$VSC_DATA_VO/resources/Ensembl_genomes/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.chrIS_spikes_45S.fa"
gtf = "$VSC_DATA_VO/resources/Ensembl_transcriptomes/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.109.chrIS_spikes_45S.gtf"
bed_regions = "$VSC_DATA_VO/resources/Ensembl_bedregions/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.109.chrIS_spikes_45S_exons_sorted_merged.bed"
star_index = "$VSC_DATA_VO/resources/STAR_index/Homo_sapiens/GRCh38.dna.primary_assembly.chrIS_spikes_45S.star_index"
if read_length == 50 :
    kal_index = "$VSC_DATA_VO/resources/Kallisto_index/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.109.chrIS_spikes_45S.kallisto_index_25"
else :
    kal_index = "$VSC_DATA_VO/resources/Kallisto_index/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.109.chrIS_spikes_45S.kallisto_index"

# # MOUSE files
# ref_fasta = "$VSC_DATA_VO/resources/Ensembl_genomes/Mus_musculus/GRCm39/Mus_musculus.GRCm39.dna.primary_assembly.chrIS_spikes_45S.fa"
# gtf = "$VSC_DATA_VO/resources/Ensembl_transcriptomes/Mus_musculus/GRCm39/Mus_musculus.GRCm39.109.chrIS_spikes_45S.gtf"
# bed_regions = "$VSC_DATA_VO/resources/Ensembl_bedregions/Mus_musculus/GRCm39/Mus_musculus.GRCm39.109.chrIS_spikes_45S_exons_sorted_merged.bed"
# star_index = "$VSC_DATA_VO/resources/STAR_index/Mus_musculus/GRCm39.dna.primary_assembly.chrIS_spikes_45S.star_index"
# kal_index = "$VSC_DATA_VO/resources/Kallisto_index/Mus_musculus/GRCm39/Mus_musculus.GRCm39.109.chrIS_spikes_45S.kallisto_index"



project_dir = dirname
empty_umi = "UMI_empty.txt"

module_cutadapt = "cutadapt/3.5-GCCcore-11.2.0"
module_calib = "Calib/0.3.4-GCC-9.3.0"
module_python = "Python/3.9.6-GCCcore-11.2.0"
module_umitools = "UMI-tools/1.0.1-foss-2019b-Python-3.7.4"
module_star = "STAR/2.7.6a-GCC-10.2.0"
module_samtools = "SAMtools/1.15-GCC-11.2.0"
module_pysam = "Pysam/0.17.0-GCC-11.2.0"
module_scipy = "SciPy-bundle/2021.10-foss-2021b"
module_biopython = "Biopython/1.79-foss-2021b"
module_java = "Java/13.0.2"
module_picard = "picard/2.25.1-Java-11"
module_fastqc = "FastQC/0.11.9-Java-11"
module_htseq = "HTSeq/0.11.3-foss-2021b"
module_kallisto = "kallisto/0.46.1-iimpi-2020b"
module_bx = "bx-python/0.8.9-foss-2020a-Python-3.8.2"
module_rseqc = "RSeQC/4.0.0-foss-2020a-Python-3.8.2"

print(
	"""
	Submitting jobs to the HPC Infrastructure...
		mode = %s
		data_type = %s
		consensus calling with calib = %s
		analysis type = %s
		base directory of input files = %s ; searching for directories starting with %s
		genome build = %s
		UMI file = %s
		errors will be sent to %s
	""" % (mode, data_type, ccalling, analysis_type, base_dir, string_match, genome_build, umi_file, email)
)

print(
	"""
	Resources that will be used:
		reference fasta file: %s
		star index: %s
		gtf used for HTseq: %s
		bed file for rseqc infer strandedness: %s
		pipeline directory: %s
	""" % (ref_fasta, star_index, gtf, bed_regions, project_dir)
)

tasks = 8
mem = 40


def sbatch(job_name, command, index, mem = mem, tasks = tasks, workdir = '.', time='11:00:00', dep=''):
	printlines = [
		"#!/bin/bash",
		"",
		"#SBATCH -J {}".format(job_name),
		"#SBATCH -D {}".format(workdir),
		"#SBATCH --mem={}G".format(mem),
		"#SBATCH --ntasks-per-node={}".format(tasks),
		"#SBATCH -t {}".format(time),
		"#SBATCH --mail-user={}".format(email),
		"#SBATCH --mail-type=FAIL",
		"#SBATCH -o {1}/logs/{2:02d}_{0}.out".format(job_name, workdir, index),
		"#SBATCH -e {1}/logs/{2:02d}_{0}.err".format(job_name, workdir, index),
	]
	if dep != '':
		if type(dep) is list:
			depstring =''
			for element in dep:
				depstring = depstring + 'afterok:' + element + ','
			dep = depstring.rstrip(',')
		else :
			dep = 'afterok:'+dep
		deplines = [
			"#SBATCH --dependency={}".format(dep),
			"#SBATCH --kill-on-invalid-dep=yes"
		]
		printlines.extend(deplines)
	printlines.extend([""])
	printlines.extend(command.split("; "))
	printjobfile('{0}/scripts/{1:02d}_{2}.sh'.format(output_dir, index, job_name),printlines)

	sbatch_command = 'sbatch {0}/scripts/{1:02d}_{2}.sh'.format(output_dir, index, job_name)
	sbatch_response = subprocess.getoutput(sbatch_command)
	print(sbatch_response)
	job_id = sbatch_response.split(' ')[3].strip()
	# job_id = job_name
	return job_id

def printjobfile(filename, printlines):
	with open(filename, 'w') as the_file:
		for line in printlines:
			if line.startswith('#'):
				the_file.write(line+"\n")
			else:
				the_file.write(line.replace('/scratch/gent/vo/000/gvo00027','$VSC_SCRATCH_VO')+"\n")


def combinecopy(sampleID):
	if data_type == 'se':
		command = 'test -s {0}/{2}.fastq.gz && cp {0}/{2}.fastq.gz {1}/{2}.fastq.gz || cat {0}/*R1*.fastq.gz > {1}/{2}.fastq.gz; '.format(input_dir, output_dir, sampleID)
	elif data_type == 'pe':
		command = 'test -s {0}/{2}*_1.fastq.gz && cp {0}/{2}*_1.fastq.gz {1}/{2}_1.fastq.gz || cat {0}/*R1*.fastq.gz > {1}/{2}_1.fastq.gz; '.format(input_dir, output_dir, sampleID)
		command = command + 'test -s {0}/{2}*_2.fastq.gz && cp {0}/{2}*_2.fastq.gz {1}/{2}_2.fastq.gz || cat {0}/*R2*.fastq.gz > {1}/{2}_2.fastq.gz; '.format(input_dir, output_dir, sampleID)
	job_id = sbatch('combinecopy_'+sampleID,command,index,workdir=output_dir,tasks=4,mem=16,time='1:00:00')
	return job_id

def adaptortrimming(sampleID, dep=''):
	command = 'module purge; module load {0}; '.format(module_cutadapt)
	command = command + 'mkdir -p {0}/{1:02d}_trimadaptout; '.format(output_dir, index)
	if data_type == 'se':
		command = command + 'cutadapt -q 30 -l {3} -a AGATCGGAAGAGC -m 20 -o {0}/{1:02d}_trimadaptout/{2}_trim.fastq.gz {0}/{1:02d}_trimadaptout/{2}.fastq.gz; '.format(output_dir, index, sampleID, read_length)
	elif data_type == 'pe':
		if analysis_type == 'RNAexome':
			command = command + 'cutadapt -q 30 -l {3} --pair-filter=any -a AGATCGGAAGAGC -A AGATCGGAAGAGC -m 20 -o {0}/{1:02d}_trimadaptout/{2}_1_trim.fastq.gz -p {0}/{1:02d}_trimadaptout/{2}_2_trim.fastq.gz {0}/{2}_1.fastq.gz {0}/{2}_2.fastq.gz; '.format(output_dir, index, sampleID, read_length)
			command = command + 'mv {0}/{1:02d}_trimadaptout/{2}_1_trim.fastq.gz {0}/{2}_1_trim.fastq.gz; '.format(output_dir, index, sampleID)
			command = command + 'mv {0}/{1:02d}_trimadaptout/{2}_2_trim.fastq.gz {0}/{2}_2_trim.fastq.gz; '.format(output_dir, index, sampleID)
		elif analysis_type == 'picoV3':
			command = command + 'cutadapt -q 30 -l {3} --pair-filter=any -a AGATCGGAAGAGC -A AGATCGGAAGAGC -m 20 -o {0}/{1:02d}_trimadaptout/{2}_1_trim.fastq.gz -p {0}/{1:02d}_trimadaptout/{2}_2_trim.fastq.gz {0}/{2}_1.fastq.gz {0}/{2}_2.fastq.gz; '.format(output_dir, index, sampleID, read_length)
			command = command + 'mv {0}/{1:02d}_trimadaptout/{2}_1_trim.fastq.gz {0}/{2}_1_trim.fastq.gz; '.format(output_dir, index, sampleID)
			command = command + 'mv {0}/{1:02d}_trimadaptout/{2}_2_trim.fastq.gz {0}/{2}_2_trim.fastq.gz; '.format(output_dir, index, sampleID)
	job_id = sbatch('adaptortrimming_'+sampleID,command,index,workdir=output_dir,dep=dep,time='6:00:00')
	return job_id

def consensuscalling(sampleID, dep=''):
	cc_par = '-l 6 -e 0 -k 8 -m 7 -t 4' if analysis_type == 'RNAexome' else '-l1 0 -l2 8 -e 0 -k 4 -m 5 -t 4'
	command = 'module purge; module load {0}; '.format(module_calib)
	command = command + 'mkdir -p {0}/{1:02d}_ccallingout; '.format(output_dir, index)
	command = command + 'gunzip -c {0}/{1}_1_trim.fastq.gz > {0}/{2:02d}_ccallingout/{1}_1_trim.fastq; '.format(output_dir, sampleID, index)
	command = command + 'gunzip -c {0}/{1}_2_trim.fastq.gz > {0}/{2:02d}_ccallingout/{1}_2_trim.fastq; '.format(output_dir, sampleID, index)
	command = command + 'bash {0}/polyNG_filter.sh {1}/{2:02d}_ccallingout/{3}_1_trim.fastq {1}/{2:02d}_ccallingout/{3}_2_trim.fastq {1}/{2:02d}_ccallingout/{3}_1_filter.fastq {1}/{2:02d}_ccallingout/{3}_2_filter.fastq ;'.format(project_dir, output_dir, index, sampleID)
	command = command + 'calib -f {0}/{1:02d}_ccallingout/{2}_1_filter.fastq -r {0}/{1:02d}_ccallingout/{2}_2_filter.fastq {3} -o {0}/{1:02d}_ccallingout/{2}.; '.format(output_dir, index, sampleID, cc_par)
	command = command + 'calib_cons -c {0}/{1:02d}_ccallingout/{2}.cluster -q {0}/{1:02d}_ccallingout/{2}_1_filter.fastq {0}/{1:02d}_ccallingout/{2}_2_filter.fastq -o {0}/{1:02d}_ccallingout/{2}_1_cons {0}/{1:02d}_ccallingout/{2}_2_cons; '.format(output_dir, index, sampleID)
	command = command + 'cat {0}/{1:02d}_ccallingout/{2}_1_cons.fastq  | cut  -f1 > {0}/{1:02d}_ccallingout/{2}_1_cons2.fastq; '.format(output_dir, index, sampleID)
	command = command + 'cat {0}/{1:02d}_ccallingout/{2}_2_cons.fastq  | cut  -f1 > {0}/{1:02d}_ccallingout/{2}_2_cons2.fastq; '.format(output_dir, index, sampleID)
	command = command + 'rm {0}/{1:02d}_ccallingout/{2}_1_cons.fastq; '.format(output_dir, index, sampleID)
	command = command + 'rm {0}/{1:02d}_ccallingout/{2}_2_cons.fastq; '.format(output_dir, index, sampleID)
	command = command + 'rm {0}/{1:02d}_ccallingout/{2}_1_filter.fastq; '.format(output_dir, index, sampleID)
	command = command + 'rm {0}/{1:02d}_ccallingout/{2}_2_filter.fastq; '.format(output_dir, index, sampleID)
	command = command + 'mv {0}/{1:02d}_ccallingout/{2}_1_cons2.fastq {0}/{1:02d}_ccallingout/{2}_1_cons.fastq; '.format(output_dir, index, sampleID)
	command = command + 'mv {0}/{1:02d}_ccallingout/{2}_2_cons2.fastq {0}/{1:02d}_ccallingout/{2}_2_cons.fastq; '.format(output_dir, index, sampleID)
	command = command + 'gzip {0}/{1:02d}_ccallingout/*cons*fastq; '.format(output_dir, index)
	command = command + 'mv {0}/{1:02d}_ccallingout/{2}_1_cons.fastq.gz {0}/{2}_1_cons.fastq.gz; '.format(output_dir, index, sampleID)
	command = command + 'mv {0}/{1:02d}_ccallingout/{2}_2_cons.fastq.gz {0}/{2}_2_cons.fastq.gz; '.format(output_dir, index, sampleID)
	command = command + 'rm {0}/{1:02d}_ccallingout/{2}*.fastq; '.format(output_dir, index, sampleID)
	job_id = sbatch('consensuscalling_'+sampleID,command,index,workdir=output_dir,dep=dep,tasks=10,mem=60)
	return job_id

#def extract(sampleID, dep=''):
#	command = 'module purge; module load {0}; '.format(module_python)
#	command = command + 'mkdir -p {0}/{1:02d}_extractout; '.format(output_dir, index)
#	if data_type == 'se':
#		command = command + 'python {0}/UMI_trim_RNA_exome.py -i {1}/{2}.fastq.gz -o {1}/{3:02d}_extractout/{2}_extract.fastq.gz -u {4}; '.format(project_dir, output_dir, sampleID, index, umi_file)
#	elif data_type == 'pe':
#		command = command + 'python {0}/UMI_trim_RNA_exome.py -i {1}/{2}_1.fastq.gz -o {1}/{3:02d}_extractout/{2}_1_extract.fastq.gz -u {4}; '.format(project_dir, output_dir, sampleID, index, umi_file)
#		command = command + 'python {0}/UMI_trim_RNA_exome.py -i {1}/{2}_2.fastq.gz -o {1}/{3:02d}_extractout/{2}_2_extract.fastq.gz -u {4}; '.format(project_dir, output_dir, sampleID, index, umi_file)
#	job_id = sbatch('extract_'+sampleID,command,index,workdir=output_dir,dep=dep,tasks=3,mem=60)
#	return job_id

def UMIextract(sampleID, dep=''):
	fname_ext = '_cons' if ccalling == 'yes' else '_trim'
	command = 'module purge; module load {0}; '.format(module_umitools)
	command = command + 'mkdir -p {0}/{1:02d}_UMIextractout; '.format(output_dir, index)
	if data_type == 'se':
		command = command + 'umi_tools extract --extract-method=string -I {0}/{2}.fastq.gz --bc-pattern=NNNNNN -L {0}/{1:02d}_UMIextractout/UMIextract.log -S {0}/{1:02d}_UMIextractout/{2}_extract.fastq.gz; '.format(output_dir, index, sampleID, index-1)
	elif data_type == 'pe':
		if analysis_type == 'RNAexome':
			command = command + 'umi_tools extract --extract-method=string -I {0}/{2}_1{3}.fastq.gz --bc-pattern=NNNNNN --bc-pattern2=NNNNNN -L {0}/{1:02d}_UMIextractout/UMIextract.log --read2-in={0}/{2}_2{3}.fastq.gz -S {0}/{1:02d}_UMIextractout/{2}_1_extract.fastq.gz --read2-out={0}/{1:02d}_UMIextractout/{2}_2_extract.fastq.gz; '.format(output_dir, index, sampleID, fname_ext)
		elif analysis_type == 'picoV3':
			command = command + 'umi_tools extract --extract-method=string -I {0}/{2}_2{3}.fastq.gz --bc-pattern=NNNNNNNN -L {0}/{1:02d}_UMIextractout/UMIextract.log --read2-in={0}/{2}_1{3}.fastq.gz -S {0}/{1:02d}_UMIextractout/{2}_2_extract.fastq.gz --read2-out={0}/{1:02d}_UMIextractout/{2}_1_extract.fastq.gz; '.format(output_dir, index, sampleID, fname_ext)
	job_id = sbatch('UMIextract_'+sampleID,command,index,workdir=output_dir,dep=dep)
	return job_id

def removeUMI(sampleID, dep=''):
	command = 'module purge; module load {0}; '.format(module_cutadapt)
	command = command + 'mkdir -p {0}/{1:02d}_cutadaptout; '.format(output_dir, index)
	if data_type == 'se':
		print("No UMI present in SE data")
	elif data_type == 'pe':
		command = command + 'cutadapt -u 6 -l {4} -o {0}/{1:02d}_cutadaptout/{2}_1_umitrim.fastq.gz {0}/{3:02d}_UMIextractout/{2}_1_extract.fastq.gz; '.format(output_dir, index, sampleID, index-1, read_length-20)
		if analysis_type == 'RNAexome':
			command = command + 'cutadapt -u 6 -l {4} -o {0}/{1:02d}_cutadaptout/{2}_2_umitrim.fastq.gz {0}/{3:02d}_UMIextractout/{2}_2_extract.fastq.gz; '.format(output_dir, index, sampleID, index-1, read_length-20)
		elif analysis_type == 'picoV3':
			command = command + 'cutadapt -u 6 -l {4} -o {0}/{1:02d}_cutadaptout/{2}_2_umitrim.fastq.gz {0}/{3:02d}_UMIextractout/{2}_2_extract.fastq.gz; '.format(output_dir, index, sampleID, index-1, read_length-20)
	job_id = sbatch('cutadapt_'+sampleID,command,index,workdir=output_dir,dep=dep)
	return job_id

def star(sampleID, dep=''):
	command = 'module purge; module load {0}; '.format(module_star)
	command = command + 'mkdir -p {0}/{1:02d}_{2}_srout; '.format(output_dir, index, sampleID)
	if data_type == 'se':
		command = command + 'STAR --chimSegmentMin 10 --runThreadN 10 --genomeDir {0} --readFilesCommand gunzip -c --outSAMtype BAM SortedByCoordinate --outFileNamePrefix {1}/{2:02d}_{4}_srout/ --readFilesIn {1}/{3:02d}_UMIextractout/{4}_extract.fastq.gz; '.format(star_index, output_dir, index, index-3, sampleID)
	elif data_type == 'pe':
		command = command + 'STAR --chimSegmentMin 10 --runThreadN 10 --genomeDir {0} --readFilesCommand gunzip -c --outSAMtype BAM SortedByCoordinate --outFileNamePrefix {1}/{2:02d}_{4}_srout/ --readFilesIn {1}/{3:02d}_cutadaptout/{4}_1_umitrim.fastq.gz {1}/{3:02d}_cutadaptout/{4}_2_umitrim.fastq.gz; '.format(star_index, output_dir, index, index-2, sampleID)
	command = command + 'mv {0}/{1:02d}_{2}_srout/Aligned.sortedByCoord.out.bam {0}/{1:02d}_{2}_srout/{2}_Aligned.sortedByCoord.out.bam; '.format(output_dir, index, sampleID)
	job_id = sbatch('star_'+sampleID,command,index,workdir=output_dir,dep=dep,tasks=10,mem=80,time='71:59:00')
	return job_id

def denoise(sampleID, dep=''):
	command = 'module purge; module load {0}; module load {1}; module load {2}; module load {3}; '.format(module_samtools, module_pysam, module_scipy, module_biopython)
	command = command + 'samtools index {0}/{1:02d}_{2}_srout/{2}_Aligned.sortedByCoord.out.bam; '.format(output_dir, index-1, sampleID)
	if analysis_type == 'RNAexome':
		command = command + 'python {0}/denoiseUMI.py -i {1}/{4:02d}_{2}_srout/{2}_Aligned.sortedByCoord.out.bam -o {1} -m fixed -u {0}/{3} -p -c 10 -s 1000000; '.format(project_dir, output_dir, sampleID, umi_file, index-1)
	elif analysis_type == 'picoV3':
		command = command + 'python {0}/denoiseUMI.py -i {1}/{3:02d}_{2}_srout/{2}_Aligned.sortedByCoord.out.bam -o {1} -m random -p -c 10 -s 1000000; '.format(project_dir, output_dir, sampleID, index-1)
	job_id = sbatch('denoise_'+sampleID,command,index,workdir=output_dir,dep=dep,tasks=12,mem=80,time='23:59:00')
	return job_id

def determineErrorRate(sampleID, dep=''):
	command = 'module purge; module load {0}; '.format(module_java)
	command = command + 'mkdir -p {0}/{1:02d}_errorraterout; '.format(output_dir, index)
	command = command + 'java -jar {0}/fgbio-2.1.0.jar ErrorRateByReadPosition -i {1}/{3}_Aligned_noNoise.bam -r {4} -l {0}/ERCC_chrIS.interval_list -o {1}/{5:02d}_errorraterout/{3}_'.format(project_dir, output_dir, index-2, sampleID, ref_fasta, index)
	job_id = sbatch('errorRate_'+sampleID,command,index,workdir=output_dir,dep=dep,tasks=1,mem=60)
	return job_id

def bamindex(sampleID, dep=''):
	command = 'module purge; module load {0}; '.format(module_samtools)
	command = command + 'samtools index {0}/{2}_Aligned_noNoise.bam; '.format(output_dir, index-3, sampleID)
	job_id = sbatch('bamindex_'+sampleID,command,index,workdir=output_dir,dep=dep,time='6:00:00')
	return job_id

def addUMtag(sampleID, dep=''):
	command = 'module purge; module load {0}; '.format(module_pysam)
	command = command + 'mkdir -p {0}/{1:02d}_tagout; '.format(output_dir, index, sampleID)
	command = command + 'python {0}/addUMtag.py -i {1}/{2}_Aligned_noNoise.bam -o {1}/{3:02d}_tagout/{2}_Aligned_noNoise_tagged.bam; '.format(project_dir, output_dir, sampleID, index)
	command = command + 'module purge; module load SAMtools/1.15-GCC-11.2.0; '
	command = command + 'samtools index {0}/{1:02d}_tagout/{2}_Aligned_noNoise_tagged.bam; '.format(output_dir, index, sampleID)
	job_id = sbatch('addUMtag_'+sampleID,command,index,workdir=output_dir,dep=dep)
	return job_id

def dedup(sampleID, dep=''):
	paired = '--paired ' if data_type == 'pe' else ''
	command = 'module purge; module load {0}; '.format(module_umitools)
	command = command + 'umi_tools dedup -I {0}/{1:02d}_tagout/{2}_Aligned_noNoise_tagged.bam --buffer-whole-contig --spliced-is-unique --output-stats=deduplicated {3}-S {0}/{2}_Aligned_noNoise_tagged_deduplicated.bam; '.format(output_dir, index-1, sampleID, paired)
	job_id = sbatch('dedup_'+sampleID,command,index,workdir=output_dir,dep=dep,time='24:00:00',mem=80)
	return job_id

def picard(sampleID, dep=''):
	command = 'module purge; module load {0}; '.format(module_picard)
	command = command + 'mkdir -p {0}/{1:02d}_picardout; '.format(output_dir, index)
	command = command + 'java -jar ${{EBROOTPICARD}}/picard.jar MarkDuplicates I={0}/{1:02d}_tagout/{3}_Aligned_noNoise_tagged.bam O={0}/{3}_Aligned_noNoise_tagged_picard_dedup.bam M={0}/{2:02d}_picardout/picard_dedup_metrics_barcode.txt VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true; '.format(output_dir, index-2, index, sampleID)
	job_id = sbatch('picard_'+sampleID,command,index,workdir=output_dir,dep=dep)
	return job_id

def fastqc(sampleID, dep=''):
	command = 'module purge; module load {0}; '.format(module_fastqc)
	command = command + 'fastqc {0}/*.fastq.gz; '.format(output_dir)
	command = command + 'fastqc {0}/{1:02d}_UMIextractout/*.fastq.gz; '.format(output_dir, index-10)
	if data_type == 'pe':
		command = command + 'fastqc {0}/{1:02d}_cutadaptout/*_umitrim.fastq.gz; '.format(output_dir, index-9)
	job_id = sbatch('fastqc_'+sampleID,command,index,workdir=output_dir,dep=dep,time='6:00:00')
	return job_id

#  op umtagged, picard_dedup, mapped_sorted_dedup, noNoise
def htseq(sampleID, dep=''):
	command = 'module purge; module load {0}; '.format(module_samtools)
	command = command + 'mkdir -p {0}/{1:02d}_htseqout; '.format(output_dir, index)
	command = command + 'samtools sort -o {0}/{1:02d}_tagout/{2}_nameSorted_noNoise_tagged.bam -n {0}/{1:02d}_tagout/{2}_Aligned_noNoise_tagged.bam; '.format(output_dir, index-4, sampleID)
	command = command + 'samtools sort -o {0}/{1:02d}_htseqout/{2}_nameSorted_noNoise_tagged_picard_dedup.bam -n {0}/{2}_Aligned_noNoise_tagged_picard_dedup.bam; '.format(output_dir, index, sampleID)
	command = command + 'samtools sort -o {0}/{1:02d}_htseqout/{2}_nameSorted_noNoise_tagged_deduplicated.bam -n {0}/{2}_Aligned_noNoise_tagged_deduplicated.bam; '.format(output_dir, index, sampleID)
	# command = command + 'samtools sort -o {0}/{1:02d}_htseqout/{2}_nameSorted_noNoise.bam -n {0}/{2}_Aligned_noNoise.bam; '.format(output_dir, index, sampleID)
	command = command + 'module purge; module load {0}; '.format(module_htseq)
	command = command + 'htseq-count --format bam --order name --nonunique none --stranded reverse {0}/{1:02d}_tagout/{2}_nameSorted_noNoise_tagged.bam {3} > {0}/{4:02d}_htseqout/{2}_calib_htseq_counts.txt; '.format(output_dir, index-4, sampleID, gtf, index)
	command = command + 'htseq-count --format bam --order name --nonunique none --stranded reverse {0}/{1:02d}_htseqout/{2}_nameSorted_noNoise_tagged_picard_dedup.bam {3} > {0}/{1:02d}_htseqout/{2}_PicardMarkDuplicates_htseq_counts.txt; '.format(output_dir, index, sampleID, gtf)
	command = command + 'htseq-count --format bam --order name --nonunique none --stranded reverse {0}/{1:02d}_htseqout/{2}_nameSorted_noNoise_tagged_deduplicated.bam {3} > {0}/{1:02d}_htseqout/{2}_UMItoolsDedup_htseq_counts.txt; '.format(output_dir, index, sampleID, gtf)
	#command = command + 'htseq-count --format bam --order name --nonunique none --stranded reverse {0}/{1:02d}_htseqout/{2}_nameSorted_noNoise.bam {3} > {0}/{1:02d}_htseqout/{2}_htseq_counts_noNoise.txt; '.format(output_dir, index, sampleID, gtf)
	job_id = sbatch('htseq_'+sampleID,command,index,workdir=output_dir,dep=dep)
	return job_id

def kallisto(sampleID, dep=''):
	command = 'module purge; module load {0}; '.format(module_kallisto)
	command = command + 'mkdir -p {0}/{1:02d}_klout; '.format(output_dir, index)
	if data_type == 'se':
		command = command + 'kallisto quant -t 10 -i {0} --single -o {1}/{2:02d}_klout {1}/{3:02d}_UMIextractout/{4}_extract.fastq.gz; '.format(kal_index, output_dir, index, index-2, sampleID)
	elif data_type == 'pe':
		command = command + 'kallisto quant -t 10 -i {0} --rf-stranded -o {1}/{2:02d}_klout {1}/{3:02d}_cutadaptout/{4}_1_umitrim.fastq.gz {1}/{3:02d}_cutadaptout/{4}_2_umitrim.fastq.gz; '.format(kal_index, output_dir, index, index-1, sampleID)
	job_id = sbatch('kallisto_'+sampleID,command,index,workdir=output_dir,dep=dep)
	return job_id

# op srout
def rseqc(sampleID, dep=''):
	command = 'module purge; module load {0}; module load {1}; '.format(module_bx, module_rseqc)
	command = command + 'mkdir -p {0}/{1:02d}_reseqcout; '.format(output_dir, index)
	command = command + 'infer_experiment.py -r {0} -i {1}/{2}_Aligned_noNoise.bam > {1}/{3:02d}_reseqcout/{2}_RSeQC_output_all.txt; '.format(bed_regions, output_dir, sampleID, index)
	command = command + 'out=`cat {0}/{1:02d}_reseqcout/{2}_RSeQC_output_all.txt | grep "1+-" | cut -d":" -f2`; '.format(output_dir, index, sampleID)
	command = command + 'echo {2} $out > {0}/{1:02d}_reseqcout/{2}_RSeQC_output.txt; '.format(output_dir, index, sampleID)
	job_id = sbatch('reseqc_'+sampleID,command,index,workdir=output_dir,dep=dep)
	return job_id

#  op srout, picard_dedup, mapped_sorted_dedup
def idxstat(sampleID, dep=''):
	command = 'module purge; module load {0}; '.format(module_samtools)
	command = command + 'mkdir -p {0}/{1:02d}_idxstatout; '.format(output_dir, index)
	command = command + 'samtools idxstats {0}/*_tagout/{1}_Aligned_noNoise_tagged.bam > {0}/{2:02d}_idxstatout/{1}_idxstat_total.txt; '.format(output_dir, sampleID, index)
	command = command + 'samtools index {0}/{1}_Aligned_noNoise_tagged_picard_dedup.bam; '.format(output_dir, sampleID)
	command = command + 'samtools idxstats {0}/{1}_Aligned_noNoise_tagged_picard_dedup.bam > {0}/{2:02d}_idxstatout/{1}_idxstat_dedup.txt; '.format(output_dir, sampleID, index)
	command = command + 'samtools idxstats {0}/{1}_Aligned_noNoise_tagged_deduplicated.bam > {0}/{2:02d}_idxstatout/{1}_idxstat_mappedSorted.txt; '.format(output_dir, sampleID, index)
	job_id = sbatch('idxstat_'+sampleID,command,index,workdir=output_dir,dep=dep)
	return job_id

def cleanup(sampleID, dep=''):
	command = 'module purge; '
	command = command + 'rm -R {0}/*out/; '.format(output_dir)
	job_id = sbatch('cleanup_'+sampleID,command,index,workdir=output_dir,dep=dep,tasks=1,mem=5)
	return job_id

def open2group(dep=''):
	command = 'module purge; '
	command = command + 'chgrp -R gvandesompele_lab {0}; '.format(output_dir)
	job_id = sbatch('open2group_'+samplename,command,index,workdir=output_dir,dep=dep,tasks=1,mem=5)
	return job_id

for samplename in os.listdir(base_dir):
	if os.path.isdir(os.path.join(base_dir,samplename)):
		if re.search(string_match,samplename):
			output_dir = args.o[0]
			input_dir = base_dir+"/"+samplename
			output_dir = output_dir+"/"+samplename

			os.makedirs(output_dir, exist_ok=True)
			os.makedirs(output_dir+"/scripts", exist_ok=True)
			os.makedirs(output_dir+"/tmp", exist_ok=True)
			os.makedirs(output_dir+"/logs", exist_ok=True)

			index = 0
			if mode == 'fastq':
				combine_jobid = combinecopy(samplename)
				index += 1
				cutadapt_jobid = adaptortrimming(samplename, combine_jobid)
				index += 1
				if ccalling == 'yes':
					ccalling_jobid = consensuscalling(samplename,cutadapt_jobid)
					index += 1
					cutadapt_jobid = ccalling_jobid
				extract_jobid = UMIextract(samplename, cutadapt_jobid)
				index += 1
				if data_type == 'pe':
					removeUMI_jobid = removeUMI(samplename, extract_jobid)
				else :
					removeUMI_jobid = extract_jobid
				index += 1
				kallisto_jobid = kallisto(samplename, removeUMI_jobid)
				index += 1
				# index = 5
				star_jobid = star(samplename, removeUMI_jobid)
				index += 1
				denoise_jobid = denoise(samplename, star_jobid)
				index += 1
				errorrate_jobid = determineErrorRate(samplename, denoise_jobid)
				index += 1
				index_jobid = bamindex(samplename, denoise_jobid)
				index += 1
				addUMtag_jobid = addUMtag(samplename, index_jobid)
				index += 1
				dedup_jobid = dedup(samplename, addUMtag_jobid)
				index += 1
				picard_jobid = picard(samplename, addUMtag_jobid)
				index += 1
				fastqc_jobid = fastqc(samplename, picard_jobid)
				index += 1
				# index = 13
				htseq_jobid = htseq(samplename, [addUMtag_jobid, dedup_jobid, denoise_jobid, picard_jobid])
				index += 1
				reseqc_job_id = rseqc(samplename, addUMtag_jobid)
				index += 1
				idxstat_jobid = idxstat(samplename, [addUMtag_jobid, dedup_jobid, picard_jobid])
				index += 1
			# elif mode == 'BAM':
			# 	extract_jobid = extract(samplename)
			# 	index += 1
			# 	dedup_jobid = dedup(samplename, extract_jobid)
			# 	index += 1

			# cleanup_jobid = cleanup(samplename, dedup_jobid)
			# index += 1
			open2group_jobid = open2group(idxstat_jobid)
