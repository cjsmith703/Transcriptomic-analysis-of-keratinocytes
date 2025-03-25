//Bulk RNA-seq pipeline using HISAT2 for alignment
//Author: Dr Chris J Smith

//Conda environment yaml file can be found in conda_envs/nextflow_rnaseq.yaml
/*
 * Pipeline parameters
 */

// Execution environment setup
params.scratchDir = "/path/to/data/star_counts/"
params.outdir = "/path/to/data/output"
params.star = "/path/to/reference/star/data"

// Primary input
params.inputDir = "${params.scratchDir}/input"

// Accessory files
params.gtf = "/path/to/data/Homo_sapiens.GRCh38.101.gtf.gz"

// Create channel from input files
Channel
    .fromFilePairs("${params.inputDir}/*-{1,2}.fastq")
    .map { filename, files -> 
         def sampleName = filename.toString().split('-')[0]  // Extract everything before the first dash
        tuple(
            files[0],                    // fastq1
            files[1],                    // fastq2
            sampleName,                  // e.g., "KO_A"
        )
    }
    .set { fastq_ch }

/*
 * FASTQC
 */

process fastqc {
    beforeScript "module load miniforge"

    input:
       tuple path(fastq1), 
             path(fastq2), 
             val(sampleName)

    // Define the output as tuple to keep sample name
    output:    
        tuple val(sampleName),
              path("${sampleName}-{1,2}_fastqc.html")

    // Copy to output directory
    publishDir "${params.outdir}/fastqc", mode: 'copy'
    
    // Set the command to be executed
    script:
    """
    mamba activate nextflow_rnaseq

    fastqc ${fastq1} ${fastq2}

    """
}

/*
 * Trim Fastq files
 */

process TrimGalore {
    beforeScript "module load miniforge"        
    
    input:
       tuple path(fastq1), 
             path(fastq2), 
             val(sampleName)

    // Define the output as tuple to keep sample name
    output:
        tuple val(sampleName),
              path("${sampleName}-{1,2}_val_{1,2}.fq")
 
    // Set the command to be executed
    script:
    """
    mamba activate nextflow_rnaseq

    trim_galore --paired ${fastq1} ${fastq2}

    """
}


/*
 * Align to reference genome
 */

process star {
    
    cpus 4
    memory '100 GB'
    
    beforeScript "module load miniforge"     

    input:
        tuple val(sampleName),
              path(trimmed_reads)  // This will receive both files as a list

    output:
        path("${sampleName}ReadsPerGene.out.tab")

    //copy the output to the output directory
    publishDir "${params.outdir}/gene_counts", mode: 'copy'
 
    script:
    """
    mamba activate nextflow_rnaseq

    STAR --runThreadN 4 \
    --genomeDir ${params.star} \
    --readFilesIn ${trimmed_reads[0]} ${trimmed_reads[1]} \
    --outFileNamePrefix ${sampleName} \
    --quantMode GeneCounts 
    """
}

process mergecounts {
    beforeScript "module load miniforge"      

    input:
        path(counts_files)  // Collect all count files

    output:
        path("star_counts_matrix.csv")

    //copy the output to the output directory
    publishDir "${params.outdir}/gene_counts", mode: 'copy'
   
      script:
    """
    mamba activate nextflow_rnaseq
    
    cat <<EOF > merge_counts.R

    library(tidyverse)
    counts_files <- list.files(pattern = "*ReadsPerGene.out.tab", full.names = TRUE)
    counts_list <- lapply(counts_files, function(file) {
    df <- read.delim(file, skip = 4)
    sample_name <- gsub("ReadsPerGene.out.tab", "", basename(file))
    df <- df[, c(1,2)]
    colnames(df)[2] <- sample_name
    return(df)
    })
    merged_counts <- Reduce(function(x, y) merge(x, y, all=TRUE), counts_list)
    
    colnames(merged_counts)[1] <- "Gene_id"

    write.csv(merged_counts, "star_counts_matrix.csv", row.names=FALSE)

    EOF

    Rscript merge_counts.R
    """
}


workflow {

  
    //Perform fastqc on the input fastq files
    fastqc(fastq_ch)

    //Trim the fastq files
    TrimGalore(fastq_ch)

    //Align the trimmed fastq files to the reference genome
    star(TrimGalore.out)
    
    //Merge the counts files
    mergecounts(star.out.collect())
    
    
    }


