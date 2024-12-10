# Overlap peaks with eqtl : 

############# scent gr load + filter ##################

sig_scent = 0.1

scent_files_substrings = c("500kb_fibroblast_allcvar",
                           "500kb_Tcell_allcvar")

############## testing enrichment significance functions ###################

library(regioneR, quietly =  T)

# number of bootstraps
ntimes = 100

# make overlap function only counting overlaps once
numOverlaps_once = function(A, B, ...){
  
  A <- toGRanges(A)
  B <- toGRanges(B)
  resampled <- numOverlaps(A,B, count.once = T)
  return(resampled)
}

# making permutation value for hg38 ... 
randomizeRegions_hg38 = function(A, ...) {
  randomizeRegions(A, genome = "hg38", ...)
}


############### gwas enrichment for loop ############

# what study names to do enrichment for

gwas_study_accessions = 
  c(# T1D
    "GCST90454277", # Zeng - european 
    "GCST90018705", # Sakaue - east asian
    "GCST009916", # Onengut-Gumusc African
    "GCST008377", # Zhu east asian
    "GCST90000529", #Inshaw European
    # RA
    "GCST90132222", # all Ishigaki 
    "GCST90132223", # european Ishigaki 
    "GCST90132224", # east asian Ishigaki 
    "GCST011389", #Kwon  east asian 
    "GCST90018690", # Sakaue  east asian 
    "GCST90446530" # Guo
  )

t1d_associations_gr = readRDS("output/gwas/gwas-catalog-associations-T1D.rds")
ra_associations_gr = readRDS("output/gwas/gwas-catalog-associations-RA.rds")

associations_gr = c(t1d_associations_gr,
                    ra_associations_gr)

for(scent_file_name in scent_files_substrings){
  
  scent_gr = regioneR::toGRanges(readRDS(paste0("output/SCENT/", 
                                                scent_file_name,
                                                ".rds"
  )
  ),
  genome = "hg38"
  )
  
  scent_gr_filt =  scent_gr[mcols(scent_gr)$boot_basic_p < sig_scent]

for(study in gwas_study_accessions){

gwas_gr = associations_gr[mcols(associations_gr)$STUDY_ACCESSION == study]

authors = c(mcols(gwas_gr)$FIRST_AUTHOR)[1]
cohort = c(mcols(gwas_gr)$INITIAL_SAMPLE_SIZE)[1]
trait = c(mcols(gwas_gr)$'DISEASE/TRAIT')[1]

# Let's now measure scent peak recall of fine-mapped eqtls
overlap = regioneR::numOverlaps(gwas_gr, 
                                scent_gr_filt, 
                                count.once = T
                                )

# total positives
n_gwas = length(gwas_gr)

# total 'labelled' positive
n_sig_scent = length(scent_gr_filt)

# recall  = tp / p 
recall = overlap / n_gwas

sig_testing = regioneR::permTest(
         A=gwas_gr,
         B=scent_gr_filt,
         ntimes=ntimes,
         randomize.function=randomizeRegions_hg38, 
         evaluate.function=numOverlaps_once
        )

recall_pval = sig_testing$numOverlaps_once$pval
recall_zscore = sig_testing$numOverlaps_once$zscore

scent_tissue = stringr::str_remove(scent_file_name,
                                   pattern = "500kb_") |>
  stringr::str_remove(pattern = "_allcvar")


# save summary results from recall: 
recall_summary_res  = data.frame(
  study = study,
  trait = trait,
  authors = authors,
  cohort = cohort,
  n_recall = overlap,
  recall = recall,
  # recall_zscore = recall_zscore,
  # recall_pval = recall_pval,
  n_gwas = n_gwas,
  n_sig_scent = n_sig_scent,
  sig_scent = sig_scent,
  scent_tissue = scent_tissue,
  ntimes = ntimes
)

data.table::fwrite(recall_summary_res, 
                   file = "output/enrichment/recall_gwas.tsv", 
                   append = TRUE, 
                   col.names = FALSE
)



message("\n For study: ", study,
        "\n trait: ", trait,
        "\n Cohort: ", cohort, 
        "\n bootstrap, p-value < ", sig_scent,
        "\n Recall: ", round(recall, digits = 2),
        "\n Recalled #n: ", overlap, 
        "\n N gwas associations: ", n_gwas, 
        "\n N scent: ", n_sig_scent,
        "\n"
        )

}
  
}
