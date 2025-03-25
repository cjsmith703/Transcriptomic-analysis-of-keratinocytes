//Bulk RNA-seq pipeline using HISAT2 for alignment
//Author: Dr Chris J Smith

//Conda environment yaml file can be found in conda_envs/nextflow_rnaseq.yaml

/*
 * Pipeline parameters
 */

// Execution environment setup
params.scratchDir = "/path/to/data/"
params.outdir = "/path/to/data/output"
params.reference = "/path/to/reference/data/hisat2_index_human/homosapiens.GRCh38"
params.gtf = "/path/to/reference/data/Homo_sapiens.GRCh38.101.gtf.gz"

// Primary input
params.inputDir = "${params.scratchDir}/input"

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

process HISAT2 {
    beforeScript "module load miniforge"        
    
    input:
        tuple val(sampleName),
              path(trimmed_reads)  // This will receive both files as a list

    output:
        tuple val(sampleName),
              path("${sampleName}.sam")
 
    script:
    """
    mamba activate nextflow_rnaseq

    hisat2 -x ${params.reference} \
    -1 ${trimmed_reads[0]} \
    -2 ${trimmed_reads[1]} \
    -S ${sampleName}.sam
    """
}

// /*
//  * Convert SAM to BAM
//  */
process SamToBam {
    beforeScript "module load samtools"
    
    //input as tuple to keep sample name
    input:
        tuple val(sampleName),
              path(aligned_sam)

    output:
        tuple val(sampleName),
              path("${sampleName}.aligned.bam"),
              path("${sampleName}.aligned.bai")

    //copy the output to the output directory
    publishDir "${params.outdir}/bam", mode: 'copy'
     
    script:
    """
    samtools view -Sb \
    ${sampleName}.sam | \
    samtools sort -o ${sampleName}.aligned.bam;
    """
}

/*
 * Quality assess BAM files
 */

process qualimap {
    beforeScript "module load miniforge"        
    
    input:
        tuple val(sampleName),
              path(aligned_bam) 

    output:
        tuple val(sampleName),
              path("${sampleName}/qualimapReport.html"),
              path("${sampleName}/genome_results.txt"),
              path("${sampleName}/raw_data_qualimapReport/*")
 
    //copy the output to the output directory
    publishDir "${params.outdir}/qualimap", mode: 'copy'
    
    script:
    """
    mamba activate nextflow_rnaseq

    qualimap bamqc \
        -bam ${aligned_bam} \
        -outdir ${sampleName} \
        --java-mem-size=8G
    """
}

/*
 * Read quantification
 */

process featureCounts {
    beforeScript "module load miniforge"        
    
    input:
        tuple val(sampleName),
              path(aligned_bam) 

    output:
        path("${sampleName}_counts.txt")

    //copy the output to the output directory
    publishDir "${params.outdir}/featureCounts", mode: 'copy'
   
    script:
    """
    mamba activate nextflow_rnaseq

    featureCounts -p --countReadPairs -a ${params.gtf} -g gene_id -o ${sampleName}_counts.txt ${aligned_bam}
    """
}

process mergecounts {
    beforeScript "module load miniforge"        
    
    input:
        path(counts_txt)  // Collect all count files

    output:
        path("raw_counts_matrix.csv")

    //copy the output to the output directory
    publishDir "${params.outdir}/featureCounts", mode: 'copy'
   
      script:
     """
    mamba activate nextflow_rnaseq
    
    cat <<EOF > merge_counts.R

    library(tidyverse)
    counts_files <- list.files(pattern = "*_counts.txt", full.names = TRUE)
    counts_list <- lapply(counts_files, function(file) {
    df <- read.delim(file, skip = 1)
    sample_name <- gsub("_counts.txt", "", basename(file))
    df <- df[, c(1,7)]
    colnames(df)[2] <- sample_name
    return(df)
    })
    merged_counts <- Reduce(function(x, y) merge(x, y, by="Geneid", all=TRUE), counts_list)
    write.csv(merged_counts, "raw_counts_matrix.csv", row.names=FALSE)

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
    HISAT2(TrimGalore.out)
    
    //Convert the SAM files to BAM
    SamToBam(HISAT2.out)

    //Quality assess the BAM files
    qualimap(SamToBam.out)

    //Read quantification
    featureCounts(SamToBam.out)

    //Merge the counts files
    mergecounts(featureCounts.out.collect())
    
    
    }


