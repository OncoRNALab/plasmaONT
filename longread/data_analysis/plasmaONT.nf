#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.junctions_chm13 = "/scratch/gent/vo/000/gvo00027/vsc42458/JAV2201_plasmapaper/code/chm13.spikes.test.bed"
params.database_file_hg38 = "$projectDir/bin/database_file_notails"
params.junctionsGTF_chm13 = "/scratch/gent/vo/000/gvo00027/vsc42458/JAV2201_plasmapaper/code/chm13.spikes.test.gtf"
params.genome_chm13 = "/scratch/gent/vo/000/gvo00027/vsc42458/JAV2201_plasmapaper/code/chm13.spikes.fa.gz"
params.junctions_hg38 = "/scratch/gent/vo/000/gvo00027/vsc42458/JAV2201_plasmapaper/code/hg38.spikes.notail.bed"
params.database_file_s = "/data/gent/vo/000/gvo00027/RNA_seq_pipeline/hg38/database_genome"
params.junctionsGTF_hg38 = "/scratch/gent/vo/000/gvo00027/vsc42458/JAV2201_plasmapaper/code/hg38.spikes.notail.gtf"
params.transcriptome_hg38 = "/data/gent/vo/000/gvo00027/RNA_seq_pipeline/hg38/Homo_sapiens.GRCh38.cdna.all.fa"
params.genome_hg38 = "/scratch/gent/vo/000/gvo00027/vsc42458/JAV2201_plasmapaper/code/hg38.spikes.notails.fa.gz"
params.genome_hg38_2 = "/scratch/gent/vo/000/gvo00027/vsc42458/fusion_experiment/code/sequences_fusions.fa.gz"
params.transcriptome = "/data/gent/vo/000/gvo00027/RNA_seq_pipeline/ensembl_transcriptomes/Homo_sapiens.GRCh38.91_withspikes_rDNA.fa"
params.transcriptome_gtf = "/data/gent/vo/000/gvo00027/RNA_seq_pipeline/ensembl_transcriptomes/Homo_sapiens.GRCh38.91_withspikes.gtf"
params.database = "/data/gent/vo/000/gvo00027/vsc42458/polyAfilter/database_file"
params.db_isoquant_hg38 = "/data/gent/vo/000/gvo00027/RNA_seq_pipeline/hg38/hg38.ensGene.db"
params.primers = "/data/gent/vo/000/gvo00027/vsc42458/JAV2001_PlasmaONT/code/primers.fas"
params.sample = ""
params.reads = "/scratch/gent/vo/000/gvo00027/vsc42458/JAV2201_plasmapaper/basecalled_reads/rEV/rEV/*"
params.outdir = "/scratch/gent/vo/000/gvo00027/vsc42458/JAV2201_plasmapaper/output/rEV"
params.SMART = "/scratch/gent/vo/000/gvo00027/vsc42458/JAV2201_plasmapaper/code/plasmaONT/bin/SMART_primers.fas"

log.info """\
		 O N T - P I P E L I N E
		 --------------------------------------
         transcriptome: ${params.transcriptome}
         genome: ${params.genome_hg38}
         GTF: ${params.junctionsGTF}
         database: ${params.database_file}
         """
         .stripIndent()



process COMBINE {

	tag "$sample_id"
	publishDir "$params.outdir/${sample_id}/01_combine_out", mode:'copy'

	input:
	tuple val(sample_id), path(reads)

	output:
	tuple val(sample_id), path("${sample_id}_combined.fastq.gz")

	script:
	"""
	cat $reads > ${sample_id}_combined.fastq.gz
	sample_id=$sample_id
	"""
}

process PORECHOP_ABI {

	tag "$sample_id"
    publishDir "$params.outdir/${sample_id}/03_porechopabi_out", mode:'copy'

	input:
	tuple val(sample_id), path(reads)

	output:
	tuple val(sample_id), file("${sample_id}_trimmed.fastq.gz"), emit: trimmed
	tuple val(sample_id), file("${sample_id}_porechopabilog.txt"), emit: porechop_log

	script:
	"""
	porechop_abi --ab_initio -i $reads -o ${sample_id}_trimmed.fastq.gz > ${sample_id}_porechopabilog.txt
	"""
}


process PYCHOPPER {

	tag "$sample_id"
        publishDir "$params.outdir/${sample_id}/04_pychopper_out", mode:'copy'

	input:
	tuple val(sample_id), path(combined_reads), path(primer_file)

	output:
	tuple val(sample_id), path("${sample_id}_full_length.fastq.gz"), emit: fl_reads
	path "${sample_id}_stats_output.txt", emit: pychopper_log
	path "${sample_id}_primers.fas", emit: primers

	script:
	"""
	gunzip -c $combined_reads > unzipped_reads
	cat $primer_file | grep -A 1 "^Consensus_" | sed 's/Consensus.*start.*/>VNP/g' | sed 's/Consensus.*end.*/>SSP/g' | sed '/--/d' > ${sample_id}_primers.fas
	pychopper -m edlib -b $params.SMART -B 100000 -c $projectDir/bin/primer_config.txt -t $task.cpus -S ${sample_id}_stats_output.txt unzipped_reads full_length_reads
	gzip -c full_length_reads > ${sample_id}_full_length.fastq.gz
	"""
}

process MINIMAP_GENOME {

	tag "$sample_id"
	publishDir "$params.outdir/${sample_id}/05_minimap_out", mode:'copy'

	input:
	tuple val(sample_id), path(splitted_reads), path(full_length_reads)

	output:
	tuple val(sample_id), file("${sample_id}_genome_sorted.bam"), file("${sample_id}_genome_sorted_fl.bam")

	script:
	"""
	gunzip -c $params.genome_hg38 > tmp_unzip_genome
	minimap2 -ax splice --cs=long --junc-bed $params.junctions_hg38 tmp_unzip_genome $splitted_reads > ${sample_id}_genome.sam
	minimap2 -ax splice --cs=long --junc-bed $params.junctions_hg38 tmp_unzip_genome $full_length_reads > ${sample_id}_genome_fl.sam
	rm tmp_unzip_genome
	samtools view -bo genome_bam ${sample_id}_genome.sam
	samtools sort genome_bam > ${sample_id}_genome_sorted.bam
	samtools view -bo genome_fl_bam ${sample_id}_genome_fl.sam
	samtools sort genome_fl_bam > ${sample_id}_genome_sorted_fl.bam
	"""

}

process POLYAFILTER {

	tag "$sample_id"
    publishDir "$params.outdir/${sample_id}/06_polyAfilter_out", mode:'copy'

    input:
		tuple val(sample_id), path(genome_reads), path(genome_reads_fl)

    output:
        tuple val(sample_id), file("${sample_id}_genome_filtered.bam"), file("${sample_id}_genome_filtered_fl.bam")

	script:
	"""
	gunzip -c $params.genome_hg38 > tmp_unzip_genome
	samtools index $genome_reads
	#python $projectDir/bin/polyAfilter/polyAfilter.py createDB $params.junctionsGTF_hg38 $params.database_file_hg38
	python $projectDir/bin/polyAfilter/polyAfilter.py createTRANS $params.database_file_hg38 $genome_reads trans_file
	python $projectDir/bin/polyAfilter/polyAfilter.py BAMfilter -d $params.database_file_hg38 -o out_bam 300 10 $genome_reads tmp_unzip_genome trans_file
	samtools sort out_bam -o ${sample_id}_genome_filtered.bam
	samtools index $genome_reads_fl
	python $projectDir/bin/polyAfilter/polyAfilter.py createTRANS $params.database_file_hg38 $genome_reads_fl trans_file
	python $projectDir/bin/polyAfilter/polyAfilter.py BAMfilter -d $params.database_file_hg38 -o out_bam_fl 300 10 $genome_reads_fl tmp_unzip_genome trans_file
	samtools sort out_bam_fl -o ${sample_id}_genome_filtered_fl.bam
	rm tmp_unzip_genome
	"""
}

process ISOQUANT {
	tag "$sample_id"
    publishDir "$params.outdir/${sample_id}/07_isoquant_out", mode:'copy'

	input:
	tuple val(sample_id), path(filtered_reads), path(filtered_reads_fl)

	output:
	path "*"
	
	script:
	"""
	samtools index $filtered_reads
	samtools index $filtered_reads_fl
	isoquant.py --reference $params.genome_hg38 --bam $filtered_reads_fl --genedb $params.junctionsGTF_hg38 --data_type nanopore -o . --force
	isoquant.py --reference $params.genome_hg38 --bam $filtered_reads --genedb $params.junctionsGTF_hg38 --data_type nanopore -o . --force
	"""
}

process NANOPLOT {
	tag "$sample_id"
    publishDir "$params.outdir/${sample_id}/08_nanoplot_out", mode:'copy'

	input:
	tuple val(sample_id), path(mapped_all), path(mapped_fl), path(polya_all), path(polya_fl)

	output:
	path "*"
	
	script:
	"""
	NanoPlot --bam $mapped_all -o ./all_mapped_reads
	NanoPlot --bam $polya_all -o ./all_polya_reads
	"""
}

process JAFFAL {
	errorStrategy 'finish'
	tag "$sample_id"
    publishDir "$params.outdir/${sample_id}/09_jaffal_out", mode:'copy'

	input:
	tuple val(sample_id), path(reads)
	
	output:
	path "*"
	
	script:
	"""
	bpipe run /JAFFA-version-2.3/JAFFAL.groovy $reads
	"""

}


workflow {
    Channel
	.fromFilePairs( params.reads , size:-1 , checkIfExists:true) { file -> file.name.replaceAll(".fastq.gz", '') }
	.set{reads_ch}

	combined_reads_ch = COMBINE( reads_ch )
    trimmed_reads_ch = PORECHOP_ABI( combined_reads_ch )
    fl_reads_ch = combined_reads_ch \
	| join( trimmed_reads_ch.porechop_log ) \
	| PYCHOPPER
	mapped_reads_ch = trimmed_reads_ch.trimmed \
	| combine( fl_reads_ch.fl_reads, by: 0 ) \
	| MINIMAP_GENOME
	polya_reads_ch = POLYAFILTER( mapped_reads_ch )
	ISOQUANT( polya_reads_ch )
	mapped_reads_ch \
	| combine( polya_reads_ch, by : 0 ) \
	| NANOPLOT
}

    
