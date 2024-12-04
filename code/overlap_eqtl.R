# Overlap peaks with eqtl : 

############# scent gr files ##################

sig_scent = 0.1

scent_files_substrings = c("500kb_fibroblast_allcvar",
                           "500kb_Tcell_allcvar")

################ eqtl gr files ###############

# what study names to do enrichment for
eqtl_file_substring = c("commonmind_brain",
                        "rosmap_brain",
                        "geuvadis_lcl",
                        "twinsuk_lcl"
)

# pip value threshold to filter our casual/fine-mapped eQTLs by
sig_eqtl = 0.2

############## permutation numbers for enrichment testing ###################

library(regioneR)

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

############### eqtl enrichment for loop ############


for(scent_file_name in scent_files_substrings){
  
  scent_gr = regioneR::toGRanges(readRDS(paste0("output/SCENT/", 
                                                scent_file_name,
                                                ".rds"
  )
  ),
  genome = "hg38"
  )
  
  scent_gr_filt = scent_gr[mcols(scent_gr)$boot_basic_p < sig_scent]

for(eqtl_file_name in eqtl_file_substring){
  

eqtl_gr = regioneR::toGRanges(readRDS(paste0("output/eqtl/", 
                                             eqtl_file_name,
                                            "_credible_sets.rds"
                                            )
                                     ),
                              genome = "hg38"
                              )

# filter eqtl for the most strongly / likely casual variants
eqtl_gr_filt = eqtl_gr[mcols(eqtl_gr)$pip > sig_eqtl]

# Let's now measure scent peak recall of fine-mapped eqtls
overlap = regioneR::numOverlaps(eqtl_gr_filt, 
                                scent_gr_filt, 
                                count.once = T
                                )

# total positives
n_eqtl = length(eqtl_gr_filt)

# total 'labelled' positive
n_sig_scent = length(scent_gr_filt)

# recall  = tp / p 
recall = overlap / n_eqtl

sig_testing = regioneR::permTest(
         A=eqtl_gr_filt,
         B=scent_gr_filt,
         ntimes=ntimes,
         randomize.function=randomizeRegions_hg38,
         evaluate.function=numOverlaps_once
        )


print(sig_testing)

recall_pval = sig_testing$numOverlaps_once$pval
recall_zscore = sig_testing$numOverlaps_once$zscore


scent_tissue = stringr::str_remove(scent_file_name,
                                    pattern = "500kb_") |>
                stringr::str_remove(pattern = "_allcvar")

# save summary results from recall: 
recall_summary_res  = data.frame(
  study = eqtl_file_name,
  n_recall = overlap,
  recall = recall,
  recall_zscore = recall_zscore,
  recall_pval = recall_pval,
  n_eqtl = n_eqtl,
  n_sig_scent = n_sig_scent,
  sig_scent = sig_scent,
  sig_eqtl = sig_eqtl,
  scent_tissue = scent_tissue, 
  ntimes = ntimes
)

data.table::fwrite(recall_summary_res, 
                   file = "output/enrichment/recall_eqtl.tsv", 
                   append = TRUE, 
                   col.names = FALSE
)





message("\n For study: ", eqtl_file_name,
        "\n pip < ", sig_eqtl, " bootstrap p-value < ", sig_scent,
        "\n Recall: ", round(recall, digits = 2),
        "\n Recalled #n: ", overlap, 
        "\n N eQTL: ", n_eqtl, 
        "\n N scent: ", n_sig_scent,
        "\n"
        )

}
  
}
