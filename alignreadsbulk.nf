
// Runs FastQC
process CheckReadsQuality {
	
	container 'rnaseq:latest'
	publishDir "${params.fastqc_dir}", mode: "move"

	input:
	path read

	output:
	file "${read.simpleName}_fastqc.html"
	file "${read.simpleName}_fastqc.zip"

	"""
	fastqc ${read}
	"""

}

// Generates STAR index, if absent
process GenerateStarIndex {

	cpus params.cpus

	publishDir "${projectDir}", mode: 'copy'

	input:
	val index
	path fasta
	path gtf
	val overhang

	output:
	env index_toggle

	script:
	if ( index == "empty" ) {
		"""
		mkdir "${params.genome_index_dir}" && \
		STAR \
		--runThreadN ${params.cpus} \
		--runMode genomeGenerate \
		–-limitBAMsortRAM ${params.ram} \
		--limitGenomeGenerateRAM ${params.ram} \
		--genomeDir ${params.genome_index_dir}/ \
		--genomeFastaFiles ${fasta} \
		--sjdbGTFfile ${gtf} \
		--sjdbOverhang "${overhang}" \
		--genomeSAindexNbases 12
		"""
		println "STAR index generated"
		index_toggle = true
	} else {
		println "STAR index found"
		index_toggle = true
	}

}

// Runs STAR for bulk single-read RNASeq samples
process AlignReadsBulkSingle() {
	
	cpus params.cpus

	publishDir "${params.reads_mapped_dir}/${read.simpleName}", mode: "move"

	input:
	each index_toggle // Mock variable to make the process wait for GenerateStarIndex. Defined with each, otherwise only one read will be processed
	path read
	path index

	output:
	path "Aligned.sortedByCoord.out.bam"
	path "Unmapped.out.mate1"
	path "ReadsPerGene.out.tab"
    path "Log.final.out"
    path "Log.out"
    path "Log.progress.out"
    path "SJ.out.tab"

    script:
    println "Processing sample ${read.simpleName}"
	"""
	if [[ "\$(echo ${read.getExtension()})" == "gz" ]];
    then
        gzipped="--readFilesCommand zcat"
    else
        gzipped=""
    fi

	STAR \
	--runThreadN ${params.cpus} \
	--runMode alignReads \
	--genomeDir ${index} \
	--readFilesIn ${read} \
	--outSAMtype BAM SortedByCoordinate \
	--outReadsUnmapped Fastx \
	--quantMode GeneCounts \
	\${gzipped}
	"""

}

// Runs STAR for bulk paired-end RNASeq samples
process AlignReadsBulkPaired() {
	
	cpus params.cpus

	container 'rnaseq:latest'
	publishDir "${params.reads_mapped_dir}/${pair_id}", mode: "move"

	input:
	each index_toggle // Mock variable to make the process wait for GenerateStarIndex. Defined with each, otherwise only one read pair will be processed
	tuple val(pair_id), path(read1), path(read2)
	path index

	output:
	path "Aligned.sortedByCoord.out.bam"
	path "Unmapped.out.mate1"
	path "Unmapped.out.mate2"
	path "ReadsPerGene.out.tab"
    path "Log.final.out"
    path "Log.out"
    path "Log.progress.out"
    path "SJ.out.tab"

    script:
    println "Processing sample ${read1.simpleName} and ${read2.simpleName}"
	"""
	if [[ "\$(echo ${read1.getExtension()})" == "gz" ]];
    then
        gzipped="--readFilesCommand zcat"
    else
        gzipped=""
    fi

	STAR \
	--runThreadN ${params.cpus} \
	--runMode alignReads \
	--genomeDir ${index} \
	--readFilesIn ${read1} ${read2} \
	--outSAMtype BAM SortedByCoordinate \
	--outReadsUnmapped Fastx \
	--quantMode GeneCounts \
	\${gzipped}
	"""

}

workflow {

	// Paramters for genome index build
	genome_fasta = Channel.fromPath("${projectDir}/${params.resources_dir}/*.{fa,fasta}")
	genome_gtf = Channel.fromPath("${projectDir}/${params.resources_dir}/*.gtf")
	genome_index_ch = Channel.fromPath("${projectDir}/${params.genome_index_dir}/*Parameters.txt").ifEmpty("empty")
	sjdboverhang = params.read_length - 1

	// Creates star index if it doesn't exist
	GenerateStarIndex(genome_index_ch, genome_fasta, genome_gtf, sjdboverhang)

	// Running FastQC
	reads_for_fastq = Channel.fromPath("${params.reads_raw_dir}/*.{fastq,fq,fastq.gz,gz}")
	CheckReadsQuality(reads_for_fastq)

	// Running STAR alignment
	star_index_dir = Channel.value("${projectDir}/${params.genome_index_dir}")
	if (params.library == "single-read") { // Bulk single-read samples
		
		reads_for_star = Channel.fromPath("${params.reads_raw_dir}/*.{fastq,fq,fastq.gz,gz}")
		AlignReadsBulkSingle(GenerateStarIndex.out, reads_for_star, star_index_dir)

	} else if (params.library == "paired-end") { // Bulk paired-end samples
		
		reads_for_star = Channel.fromFilePairs("${params.reads_raw_dir}/${params.reads_pattern}", size: -1, flat: true)
		AlignReadsBulkPaired(GenerateStarIndex.out, reads_for_star, star_index_dir)

	} else {
		
		println "ERROR: cannot determine library type"

	}

}