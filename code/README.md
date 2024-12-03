
# Individual Analysis Code 



## eQTL enrichment steps 

<b> 1. download_eqtl.sh </b> 

Downloads, and unzips susie credible intervals from eQTL catalogue v7. 

Outputs:
- data/eqtl/*.tsv

<b> 2. prep_eqtl_data.R </b> 

Converts susie credible interval datasets into GenomicRanges objects. 

Inputs:
- data/eqtl/*.tsv

Outputs: 
- output/eqtl/*.rds

<b> 3. prep_scent_data.R </b>

Converts output from scent analysis into GenomicRanges objects. 

Inputs:
- data/scent/*.txt

Outputs: 
- output/scent/*.rds

<b> 4. overlap_eqtl.R </b> 

Uses the regioneR bioconductor R package to measure the overlap between eQTLs (eVariants), and scent-identified enhancers. 

Inputs: 
- output/scent/*.rds
- output/eqtl/*.rds

Output: 
- output/enrichment/recall_eqtl.tsv (each line refers to a single recall measurement)

<br> 
<br>

# GWAs enrichment steps 

<b> 1. prep_gwas_data.R </b>

Converts relevant association data from the gwas catalog to GenomicRanges Objects.

Output:
- output/gwas/*.rds

<b> 2. overlap_gwas.R </b> 

Uses the regioneR bioconductor R package to measure the overlap between GWAs hits / variants, and scent-identified enhancers. 

Inputs: 
- output/scent/*.rds
- output/gwas/*.rds

Output: 
- output/enrichment/recall_gwas.tsv (each line refers to a single recall measurement)


