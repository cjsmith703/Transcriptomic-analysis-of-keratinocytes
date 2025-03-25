# Transcriptomic-analysis-of-keratinocytes

This repository contains scripts and data used in bulk RNAseq analysis of SGPL1 KO nTERT keratinocytes. 

Nextflow scripts for processing of FASTQ files through to featurecounts are included, using STAR or HISAT2. The YAML file to create a conda environment with the required dependencies is also available. 

Scripts to perform DESEQ2 analysis and create heatmaps, dotplots and other graphs included in the paper are also available. DESEQ2 data has been included and  as well as TPM counts. Further data and information is available at the [GEO repository](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE207499).

The full open access paper can be found [here](https://www.sciencedirect.com/science/article/pii/S002222752300024X)

Smith CJ, Williams JL, Hall C, Casas J, Caley MP, O'Toole EA, Prasad R, Metherell LA. Ichthyosis linked to sphingosine 1-phosphate lyase insufficiency is due to aberrant sphingolipid and calcium regulation. J Lipid Res. 2023 Apr;64(4):100351. doi: 10.1016/j.jlr.2023.100351. Epub 2023 Mar 2. PMID: 36868360; PMCID: PMC10123262.Smith CJ, Williams JL, Hall C, Casas J, Caley MP, O'Toole EA, Prasad R, Metherell LA. Ichthyosis linked to sphingosine 1-phosphate lyase insufficiency is due to aberrant sphingolipid and calcium regulation. J Lipid Res. 2023 Apr;64(4):100351. doi: 10.1016/j.jlr.2023.100351.PMID: 36868360; PMCID: PMC10123262.
