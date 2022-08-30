
// Creates necessary Docker images, if absent
process CheckDockerImages() {

    """
    echo "Checking Docker images" && \
    cd ${projectDir}/${params.dockerfiles_dir}
    /bin/bash check_docker_images.sh
    cd ../
    """

}

// Runs FastQC
process CheckReadsQuality() {
	
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

// Determines if a STAR genome index exists or not already
process CheckGenomeIndex() {

    input:
    val index_dir
    val index_file

    output:
    env create_index

    script:
    if ( index_dir == 'empty' )
        """
        create_index=true
        """
    else if ( index_file == 'empty')
        """
        create_index=true
        """
    else
        """
        create_index=false
        """
	
}

// Generates STAR index, if absent
process GenerateStarIndex() {

    cpus params.cpus

    container 'rnaseq:latest'
    publishDir "${params.genome_index_dir}", mode: "copy"

    input:
    path resources
    env create_index

    output:
    path "*"

    shell:
    '''
    if [ "$(echo !{params.type})" == "single-cell" ];
    then
        match_cell_ranger="--genomeSAsparseD 3"
    else
        match_cell_ranger=""
    fi

    if [ "${create_index}" = true ];
    then
        STAR \
        --runThreadN !{params.cpus} \
        --runMode genomeGenerate \
        –-limitBAMsortRAM !{params.ram} \
        --limitGenomeGenerateRAM !{params.ram} \
        --genomeDir . \
        --genomeFastaFiles !{resources}/*.fa* \
        --sjdbGTFfile !{resources}/*.gtf \
        --sjdbOverhang $((!{params.read_length} - 1)) \
        --genomeSAindexNbases 12 \
        ${match_cell_ranger}
    fi
    '''

}

// Runs STAR for bulk single-read RNASeq samples
process AlignReadsBulkSingle() {
	
    cpus params.cpus

    container 'rnaseq:latest'
    publishDir "${params.reads_mapped_dir}/${read.simpleName}", mode: "move"

    input:
    file read
    path index
    path index_dir

    output:
    path "Aligned.sortedByCoord.out.bam"
    path "Unmapped.out.mate1"
    path "ReadsPerGene.out.tab"
    path "Log.final.out"
    path "Log.out"
    path "Log.progress.out"
    path "SJ.out.tab"

    shell:
    '''
    if [[ "$(echo !{read.getExtension()})" == "gz" ]];
    then
        gzipped="--readFilesCommand zcat"
    else
        gzipped=""
    fi

    STAR \
    --runThreadN !{params.cpus} \
    --runMode alignReads \
    --genomeDir !{params.genome_index_dir} \
    --readFilesIn !{read} \
    --outSAMtype BAM SortedByCoordinate \
    --outReadsUnmapped Fastx \
    --quantMode GeneCounts \
    ${gzipped}
    '''

}

// Runs STAR for bulk paired-end RNASeq samples
process AlignReadsBulkPaired() {
	
    cpus params.cpus

    container 'rnaseq:latest'
    publishDir "${params.reads_mapped_dir}/${pair_id}", mode: "move"

    input:
    tuple val(pair_id), path(read1), path(read2)
    path index
    path index_dir

    output:
    path "Aligned.sortedByCoord.out.bam"
    path "Unmapped.out.mate1"
    path "Unmapped.out.mate2"
    path "ReadsPerGene.out.tab"
    path "Log.final.out"
    path "Log.out"
    path "Log.progress.out"
    path "SJ.out.tab"

    shell:
    '''
    if [[ "$(echo !{read1.getExtension()})" == "gz" ]];
    then
        gzipped="--readFilesCommand zcat"
    else
        gzipped=""
    fi

    STAR \
    --runThreadN !{params.cpus} \
    --runMode alignReads \
    --genomeDir !{params.genome_index_dir} \
    --readFilesIn !{read1} !{read2} \
    --outSAMtype BAM SortedByCoordinate \
    --outReadsUnmapped Fastx \
    --quantMode GeneCounts \
    ${gzipped}
    '''

}

// Runs STAR for single-cell paired-end RNASeq samples
process AlignReadsSingleSingle() {
	
    cpus params.cpus

    container 'rnaseq:latest'
    publishDir "${params.reads_mapped_dir}/${read.simpleName}", mode: "move"

    input:
    path read
    path cellid
    path index
    path index_dir

    output:
    path "Aligned.sortedByCoord.out.bam"
    path "Unmapped.out.mate1"
    path "ReadsPerGene.out.tab"
    path "Log.final.out"
    path "Log.out"
    path "Log.progress.out"
    path "SJ.out.tab"
    path "Features.stats"
    path "Summary.csv"
    path "UMIperCellSorted.txt"
    path "features.tsv"
    path "barcodes.tsv"
    path "matrix.mtx"

    shell:
    '''
    if [[ "$(echo !{read.getExtension()})" == "gz" ]];
    then
        gzipped="--readFilesCommand zcat"
    else
        gzipped=""
    fi

    STAR \
    --runThreadN !{params.cpus} \
    --runMode alignReads \
    --genomeDir !{params.genome_index_dir} \
    --readFilesIn !{read} !{cellid} \
    --outSAMtype BAM SortedByCoordinate \
    --outReadsUnmapped Fastx \
    --quantMode GeneCounts \
    ${gzipped} \
    !{params.starsolo}
    '''

}

// Runs STAR for single-cell paired-end RNASeq samples
process AlignReadsSinglePaired() {
	
    cpus params.cpus

    container 'rnaseq:latest'
    publishDir "${params.reads_mapped_dir}/${read1.simpleName}", mode: "move"

    input:
    path read1
    path read2
    path cellid
    path index
    path index_dir

    output:
    path "Aligned.sortedByCoord.out.bam"
    path "Unmapped.out.mate1"
    path "ReadsPerGene.out.tab"
    path "Log.final.out"
    path "Log.out"
    path "Log.progress.out"
    path "SJ.out.tab"
    path "Features.stats"
    path "Summary.csv"
    path "UMIperCellSorted.txt"
    path "features.tsv"
    path "barcodes.tsv"
    path "matrix.mtx"

    shell:
    '''
    if [[ "$(echo !{read1.getExtension()})" == "gz" ]];
    then
        gzipped="--readFilesCommand zcat"
    else
        gzipped=""
    fi

    STAR \
    --runThreadN !{params.cpus} \
    --runMode alignReads \
    --genomeDir !{index_dir} \
    --readFilesIn !{read1} !{read2} !{cellid} \
    --outSAMtype BAM SortedByCoordinate \
    --outReadsUnmapped Fastx \
    --quantMode GeneCounts \
    ${gzipped} \
    !{params.starsolo}
    '''

}

workflow {
	
    // If necessary Docker images do not exist, they'll be created here
    CheckDockerImages()

    // Read quality is checked with FastQC
    if (params.type == "bulk") {
        rawreads = Channel.fromPath("${params.reads_raw_dir}/*.{fastq,fq,fastq.gz,gz}")
        CheckReadsQuality(rawreads)
    } else {
        rawreads1 = Channel.fromPath("${params.reads_raw_dir}/${params.mate1_pattern}")
        rawreads2 = Channel.fromPath("${params.reads_raw_dir}/${params.mate2_pattern}")
        CheckReadsQuality(rawreads1.mix(rawreads2))
    }
	
    // Generate STAR index if absent
    index_dir_ch = Channel.fromPath("${projectDir}/${params.genome_index_dir}").ifEmpty("empty")
    index_files_ch = Channel.fromPath("${projectDir}/${params.genome_index_dir}/*").ifEmpty("empty")
    CheckGenomeIndex("${projectDir}", "${projectDir}/${params.genome_index_dir}")
    GenerateStarIndex("${projectDir}/${params.resources_dir}", CheckGenomeIndex.out)
	
    // STAR alignment
    if (params.type == "bulk") {
        // Bulk RNASeq samples
        if (params.library == "single-read") {
            // Bulk single-read samples
            rawreads = Channel.fromPath("${params.reads_raw_dir}/*.{fastq,fq,fastq.gz,gz}")
            AlignReadsBulkSingle(rawreads, GenerateStarIndex.out, "${projectDir}/${params.genome_index_dir}")
        } else if (params.library == "paired-end") {
            // Bulk paired-end samples
            rawreads = Channel.fromFilePairs("${params.reads_raw_dir}/${params.paired_pattern}", size: 2, flat: true)
            AlignReadsBulkPaired(rawreads, GenerateStarIndex.out, "${projectDir}/${params.genome_index_dir}")
        } else {
            println "ERROR: cannot determine library type"
        }
    } else if (params.type == "single-cell") {
        // Single-cell RNASeq samples
        if (params.library == "single-read") {
            // Single-cell single-read samples
            rawreads = Channel.fromPath("${params.reads_raw_dir}/${params.mate1_pattern}")
            cellids = Channel.fromPath("${params.reads_raw_dir}/${params.cell_id_pattern}")
            AlignReadsSingleSingle(rawreads, cellids, GenerateStarIndex.out, "${projectDir}/${params.genome_index_dir}")
        } else if (params.library == "paired-end") {
            // Single-cell paired-end samples
            rawreads1 = Channel.fromPath("${params.reads_raw_dir}/${params.mate1_pattern}")
            rawreads2 = Channel.fromPath("${params.reads_raw_dir}/${params.mate2_pattern}")
            cellids = Channel.fromPath("${params.reads_raw_dir}/${params.cell_id_pattern}")
            AlignReadsSinglePaired(rawreads1, rawreads2, cellids, GenerateStarIndex.out, "${projectDir}/${params.genome_index_dir}")
        } else {
            println "ERROR: cannot determine library type"
        }
    } else {
        println "ERROR: cannot determine sequencing mode"
    }
	
}
