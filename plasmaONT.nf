#!/usr/bin/env nextflow

params.junctions_chm13 = "/data/gent/vo/000/gvo00027/RNA_seq_pipeline/chm13v2/junctions_chm13.bed"
params.database_file_chm13 = "$projectDir/bin/database_file"
params.junctionsGTF_chm13 = "/data/gent/vo/000/gvo00027/RNA_seq_pipeline/chm13v2/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf"
params.genome_chm13 = "/data/gent/vo/000/gvo00027/RNA_seq_pipeline/chm13v2/chm13v2.0.fa.gz"
params.junctions_hg38 = "/data/gent/vo/000/gvo00027/RNA_seq_pipeline/hg38/hg38.ensGene.bed"
params.database_file_hg38 = "/data/gent/vo/000/gvo00027/RNA_seq_pipeline/hg38/database_genome"
params.junctionsGTF_hg38 = "/data/gent/vo/000/gvo00027/RNA_seq_pipeline/hg38/hg38.ensGene.gtf"
params.transcriptome_hg38 = "/data/gent/vo/000/gvo00027/RNA_seq_pipeline/hg38/Homo_sapiens.GRCh38.cdna.all.fa"
params.genome_hg38 = "/data/gent/vo/000/gvo00027/RNA_seq_pipeline/hg38/hg38.fa.gz"
params.transcriptome = "/data/gent/vo/000/gvo00027/RNA_seq_pipeline/ensembl_transcriptomes/Homo_sapiens.GRCh38.91_withspikes_rDNA.fa"
params.transcriptome_gtf = "/data/gent/vo/000/gvo00027/RNA_seq_pipeline/ensembl_transcriptomes/Homo_sapiens.GRCh38.91_withspikes.gtf"
params.database = "/data/gent/vo/000/gvo00027/vsc42458/polyAfilter/database_file"
params.db_isoquant_hg38 = "/data/gent/vo/000/gvo00027/RNA_seq_pipeline/hg38/hg38.ensGene.db"
params.primers = "/data/gent/vo/000/gvo00027/vsc42458/JAV2001_PlasmaONT/code/primers.fas"
params.sample = ""
params.reads = "/data/gent/vo/000/gvo00027/vsc42458/JAV2201_plasmapaper/data/large_plasma/*/fastq_guppy/pass/*.fastq.gz"
params.outdir = "/data/gent/vo/000/gvo00027/vsc42458/JAV2201_plasmapaper/output/large_plasma"

log.info """\
		 O N T - P I P E L I N E
		 --------------------------------------
         transcriptome: ${params.transcriptome}
         genome: ${params.genome_hg38}
         GTF: ${params.junctionsGTF}
         database: ${params.database_file}
         """
         .stripIndent()

Channel
	.fromFilePairs( params.reads , size:-1 , checkIfExists:true) { file -> file.name.replaceAll("_[0-9]*_0.fastq.gz", '') }
	.set{reads_ch}

process COMBINE {

	tag "$sample_id"
	publishDir "$params.outdir/${sample_id}/01_combine_out", mode:'copy'

	input:
	tuple sample_id, path(reads) from reads_ch

	output:
	tuple sample_id, "${sample_id}_combined.fastq.gz" into (combined_reads_pychopper_ch, combined_reads_porechop_ch)

	script:
	"""
	cat $reads > ${sample_id}_combined.fastq.gz
	"""
}



process PORECHOP {

	tag "$sample_id"
        publishDir "$params.outdir/${sample_id}/02b_porechop_out", mode:'copy'

	input:
	tuple sample_id, path(combined_reads) from combined_reads_porechop_ch

	output:
	tuple sample_id, "${sample_id}_trimmed.fastq.gz" into trimmed_reads_ch
	path "${sample_id}_porechoplog.txt" into (porechop_log_a_ch, porechop_log_b_ch)

	script:
	"""
	porechop_abi --ab_initio -i $combined_reads -o ${sample_id}_trimmed.fastq.gz > ${sample_id}_porechoplog.txt
	"""
}


process PYCHOPPER {

	tag "$sample_id"
        publishDir "$params.outdir/${sample_id}/02a_pychopper_out", mode:'copy'

	input:
	tuple sample_id, path(combined_reads) from combined_reads_pychopper_ch
	path(primer_file) from porechop_log_b_ch

	output:
	tuple sample_id, "${sample_id}_full_length.fastq.gz" into full_length_reads_ch
	path "${sample_id}_stats_output.txt" into pychopper_stats_ch
	path "${sample_id}_primers.fas"

	script:
	"""
	gunzip -c $combined_reads > unzipped_reads
	cat $primer_file | grep -A 1 "^Consensus_" | sed 's/Consensus.*start.*/>VNP/g' | sed 's/Consensus.*end.*/>SSP/g' | sed '/--/d' > ${sample_id}_primers.fas
	pychopper -m edlib -b ${sample_id}_primers.fas -c $projectDir/bin/primer_config.txt -t $task.cpus -S ${sample_id}_stats_output.txt unzipped_reads full_length_reads
	gzip -c full_length_reads > ${sample_id}_full_length.fastq.gz
	"""
}

process MINIMAP_GENOME {

	tag "$sample_id"
	publishDir "$params.outdir/${sample_id}/03_minimap_out", mode:'copy'

	input:
	tuple sample_id, path(splitted_reads) from trimmed_reads_ch
	tuple sample_id, path(full_length_reads) from full_length_reads_ch

	output:
	tuple sample_id, "${sample_id}_genome_sorted.bam" into sorted_bam_ch
	tuple sample_id, "${sample_id}_genome_sorted.bam" into bam_for_nanoplot_ch
	tuple sample_id, "${sample_id}_genome_sorted_fl.bam" into sorted_bam_fl_ch
	tuple sample_id, "${sample_id}_genome_sorted_fl.bam" into bam_for_nanoplot_fl_ch

	script:
	"""
	gunzip -c $params.genome_hg38 > tmp_unzip_genome
	minimap2 -ax splice --junc-bed $params.junctions_hg38 tmp_unzip_genome $splitted_reads > ${sample_id}_genome.sam
	minimap2 -ax splice --junc-bed $params.junctions_hg38 tmp_unzip_genome $full_length_reads > ${sample_id}_genome_fl.sam
	rm tmp_unzip_genome
	samtools view -bo genome_bam ${sample_id}_genome.sam
	samtools sort genome_bam > ${sample_id}_genome_sorted.bam
	samtools view -bo genome_fl_bam ${sample_id}_genome_fl.sam
	samtools sort genome_fl_bam > ${sample_id}_genome_sorted_fl.bam
	"""

}


process NANOPLOT {

	tag "$sample_id"
	publishDir "$params.outdir/${sample_id}/04_nanoplot_out", mode:'copy'

	input:
	tuple sample_id, path(mapped_reads) from bam_for_nanoplot_ch

	output:
	path "*"
	path "${sample_id}_NanoStats.txt" into nanostats_ch

	script:
	"""
	samtools view -F 4 $mapped_reads > test
	if [ -s test ]
	then
		samtools index $mapped_reads
		NanoPlot -t $task.cpus --color yellow --bam $mapped_reads --downsample 10000 -o . -p ${sample_id}_
	else
		touch ${sample_id}_NanoStats.txt
	fi
	"""

}

process POLYAFILTER {

	tag "$sample_id"
    publishDir "$params.outdir/${sample_id}/05_polyAfilter_out", mode:'copy'

    input:
		tuple sample_id, path(genome_reads) from sorted_bam_ch
		tuple sample_id, path(genome_reads_fl) from sorted_bam_fl_ch

    output:
        tuple sample_id, "${sample_id}_genome_filtered.bam" into genome_filtered_ch
		tuple sample_id, "${sample_id}_genome_filtered_fl.bam" into (genome_filtered_fl_a_ch, genome_filtered_fl_b_ch)

	script:
	"""
	gunzip -c $params.genome_hg38 > tmp_unzip_genome
	samtools index $genome_reads
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

process STRINGTIE2 {
	
	tag "$sample_id"
    publishDir "$params.outdir/${sample_id}/06_stringtie_out", mode:'copy'

	input:
	tuple sample_id, path(filtered_reads) from genome_filtered_fl_a_ch
	
	output:
	tuple sample_id, "${sample_id}_stringtie.gtf"


	script:
	"""
	samtools sort $filtered_reads -o sorted_bam
	stringtie -G $params.junctionsGTF_hg38 -o ${sample_id}_stringtie.gtf -p $task.cpus -L sorted_bam 
	"""
}

process ISOQUANT {
	tag "$sample_id"
    publishDir "$params.outdir/${sample_id}/07_isoquant_out", mode:'copy'

	input:
	tuple sample_id, path(filtered_reads) from genome_filtered_fl_b_ch

	output:
	path "*"
	
	script:
	"""
	samtools index $filtered_reads
	isoquant.py --reference $params.genome_hg38 --genedb $params.junctionsGTF_hg38 --bam $filtered_reads --data_type nanopore -o . --force
	"""
}

process HTSEQ {

	tag "$sample_id"
    publishDir "$params.outdir/${sample_id}/08_htseq_out", mode:'copy'

    input:
		tuple sample_id, path(filtered_reads) from genome_filtered_ch	
		
    output:
        path "${sample_id}_counts.txt" into counts_ch
		path "${sample_id}_annotated.sam" 

	script:
	"""
	samtools sort -n $filtered_reads -o sorted_bam
	samtools view -h sorted_bam | awk 'length(\$10) > 1 || \$1 ~ /^@/' > sorted_sam
	htseq-count --format sam --order name --mode union --stranded no --nonunique none -o ${sample_id}_annotated.sam sorted_sam $params.junctionsGTF_hg38 > ${sample_id}_counts.txt
"""
}


process MULTIQC {
    publishDir "$params.outdir/multiqc_out", mode:'copy'

	input:

	path counts_input from counts_ch.collect().ifEmpty([])
	path porechop from porechop_log_a_ch.collect().ifEmpty([])
	path nanoplot from nanostats_ch.collect().ifEmpty([])
	path pych_input from pychopper_stats_ch.collect().ifEmpty([])


	output:
	path "multiqc_report.html"

	script:
	"""
	multiqc $counts_input $porechop $nanoplot $pych_input
	"""

}
