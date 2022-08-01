jamestown = [
	"JR0847",
	"JR118232",
	"JR118231",
	"JR118230",
	"JR118236",
	"JR68100"]

hatch = [
	"AL3773",
	"ALXXX"]

all_samples = jamestown + hatch
assemblies = ["AAEX03", "ROS_Cfam_1"]

# Tool paths
bwa_path = "bwa"
samtools_path = "samtools"
fastqc_path = "fastqc"
cutadapt_path = "cutadapt"
gatk_path = "/home/ariane/Programs/gatk-4.2.6.1/gatk"
multiqc_path = "multiqc"
picard_path = "/home/ariane/Programs/picard/build/libs/picard.jar"

rule all:
	input:
		expand(
			"refs/{assembly}.fasta.fai",
			assembly=assemblies),
		"multiqc_results/multiqc_report.html",
		"multiqc_trimmed_results/multiqc_report.html",
		expand(
			"bams/{sample}.{assembly}.sorted.bam.bai",
			sample=all_samples,
			assembly=assemblies),
		expand(
			"bams/{sample}.{assembly}.sorted.mkdup.bam.bai",
			sample=all_samples,
			assembly=assemblies),
		expand(
			"stats/{sample}.{assembly}.sorted.mkdup.bam.stats",
			sample=all_samples,
			assembly=assemblies),
		expand(
			"genotyped_vcfs/{assembly}.gatk.called.rawvariants.vcf.gz",
			assembly=assemblies)

rule index_reference:
	input:
		ref = "reference/{assembly}.fasta"
	output:
		fai = "reference/{assembly}.fasta.fai",
		amb = "reference/{assembly}.fasta.amb",
		dict = "reference/{assembly}.dict"
	params:
		samtools = samtools_path,
		bwa = bwa_path
	run:
		# faidx
		shell(
			"{params.samtools} faidx {input.ref}")
		# .dict
		shell(
			"{params.samtools} dict -o {output.dict} {input.ref}")
		# bwa
		shell(
			"{params.bwa} index {input.ref}")

rule fastqc_analysis:
	input:
		"fastq/{sample}_full.{read}.fastq.gz"
	output:
		"fastqc_results/{sample}_full.{read}_fastqc.html"
	params:
		fastqc = fastqc_path
	shell:
		"{params.fastqc} -o fastqc_results {input}"

rule multiqc_analysis:
	input:
		expand(
			"fastqc_results/{sample}_full.{read}_fastqc.html",
			sample=all_samples,
			read=["R1", "R2"])
	output:
		"multiqc_results/multiqc_report.html"
	params:
		multiqc = multiqc_path
	shell:
		"export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 && "
		"{params.multiqc} --interactive -f "
		"-o multiqc_results fastqc_results"

rule trim_adapters_paired:
	input:
		fq1 = "fastq/{sample}_full.R1.fastq.gz",
		fq2 = "fastq/{sample}_full.R2.fastq.gz"
	output:
		out_fq1 = "trimmed_fastqs/{sample}_trimmed_read1.fastq.gz",
		out_fq2 = "trimmed_fastqs/{sample}_trimmed_read2.fastq.gz"
	params:
		cutadapt = cutadapt_path
	shell:
		"{params.cutadapt} -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o trimmed_fastqs {input.fq1} {input.fq2}"

rule fastqc_analysis_trimmed:
	input:
		fq1 = "trimmed_fastqs/{sample}_trimmed_read1.fastq.gz",
		fq2 = "trimmed_fastqs/{sample}_trimmed_read2.fastq.gz"
	output:
		html1 = "fastqc_trimmed_results/{sample}_trimmed_read1_fastqc.html",
		html2 = "fastqc_trimmed_results/{sample}_trimmed_read2_fastqc.html"
	params:
		fastqc = fastqc_path
	shell:
		"{params.fastqc} -o fastqc_trimmed_results {input.fq1} {input.fq2}"

rule multiqc_analysis_trimmed:
	input:
		expand(
			"fastqc_trimmed_results/{sample}_trimmed_{read}_fastqc.html",
			sample=all_samples, read=["read1", "read2"])
	output:
		"multiqc_trimmed_results/multiqc_report.html"
	params:
		multiqc = multiqc_path
	shell:
		"export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 && "
		"{params.multiqc} --interactive -f -o multiqc_trimmed_results fastqc_trimmed_results"

rule bwa_map_sort:
	input:
		fq1 = "trimmed_fastqs/{sample}_trimmed_read1.fastq.gz",
		fq2 = "trimmed_fastqs/{sample}_trimmed_read2.fastq.gz",
		ref = "reference/{assembly}.fasta",
		fai = "reference/{assembly}.fasta.fai"
	output:
		"bams/{sample}.{assembly}.sorted.bam"
	params:
		id = "{sample}",
		sm = "{sample}",
		lb = "{sample}",
		pu = "{sample}",
		pl = "Illumina",
		bwa = bwa_path,
		samtools = samtools_path
	shell:
		" {params.bwa} mem -l 1000 -n 0.01 -R "
	 	"'@RG\\tID:{params.id}\\tSM:{params.sm}\\tLB:{params.lb}\\tPU:{params.pu}\\tPL:{params.pl}' "
		"{input.ref} {input.fq1} {input.fq2}"
		"| {params.samtools} fixmate -O bam - - | {params.samtools} sort "
		"-O bam -o {output}"

rule index_bam:
	input:
		"bams/{sample}.{assembly}.sorted.bam"
	output:
		"bams/{sample}.{assembly}.sorted.bam.bai"
	params:
		samtools = samtools_path
	shell:
		"{params.samtools} index {input}"

rule picard_mkdups:
	input:
		bam = "bams/{sample}.{assembly}.sorted.bam",
		bai = "bams/{sample}.{assembly}.sorted.bam.bai"
	output:
		bam = "bams/{sample}.{assembly}.sorted.mkdup.bam",
		metrics = "stats/{sample}.{assembly}.picard_mkdup_metrics.txt"
	params:
		picard = picard_path
	shell:
		"{params.picard} -Xmx1g MarkDuplicates I={input.bam} O={output.bam} "
		"M={output.metrics}"

rule index_mkdup_bam:
	input:
		"bams/{sample}.{assembly}.sorted.mkdup.bam"
	output:
		"bams/{sample}.{assembly}.sorted.mkdup.bam.bai"
	params:
		samtools = samtools_path
	shell:
		"{params.samtools} index {input}"

rule bam_stats:
	input:
		bam = "bams/{sample}.{assembly}.sorted.mkdup.bam",
		bai = "bams/{sample}.{assembly}.sorted.mkdup.bam.bai"
	output:
		"stats/{sample}.{assembly}.sorted.mkdup.bam.stats"
	params:
		samtools = samtools_path
	shell:
		"{params.samtools} stats {input.bam} | grep ^SN | cut -f 2- > {output}"

rule gatk_gvcf:
	input:
		ref = "reference/{assembly}.fasta",
		bam = "bams/{sample}.{assembly}.sorted.mkdup.bam",
		bai = "bams/{sample}.{assembly}.sorted.mkdup.bam.bai"
	output:
		"gvcfs/{sample}.{assembly}.g.vcf.gz"
	params:
		gatk = gatk_path
	shell:
		"""{params.gatk} --java-options "-Xmx1g" """
		"""HaplotypeCaller -R {input.ref} -I {input.bam} """
		"""-ERC GVCF -O {output}"""

rule gatk_combinegvcfs:
	input:
		ref = "reference/{assembly}.fasta",
		gvcfs = lambda wildcards: expand(
			"gvcfs/{sample}.{genome}.g.vcf.gz",
			sample=all_samples,
			genome=wildcards.assembly)
	output:
		"combined_gvcfs/{assembly}.gatk.combinegvcf.g.vcf.gz"
	params:
		gatk = gatk_path
	run:
		variant_files = []
		for i in input.gvcfs:
			variant_files.append("--variant " + i)
		variant_files = " ".join(variant_files)
		print(
			"""{params.gatk} --java-options "-Xmx1g" """
			"""CombineGVCFs -R {input.ref} {variant_files} -O {output}"""
		)
		shell(
			"""{params.gatk} --java-options "-Xmx1g" """
			"""CombineGVCFs -R {input.ref} {variant_files} -O {output}""")

rule gatk_genotypegvcf:
	input:
		ref = "reference/{assembly}.fasta",
		gvcf = "combined_gvcfs/{assembly}.gatk.combinegvcf.g.vcf.gz"
	output:
		"genotyped_vcfs/{assembly}.gatk.called.rawvariants.vcf.gz"
	params:
		gatk = gatk_path
	shell:
		"""{params.gatk} --java-options "-Xmx1g" """
		"""GenotypeGVCFs -R {input.ref} -V {input.gvcf} -O {output}"""
