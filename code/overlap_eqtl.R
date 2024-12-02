# Overlap peaks with eqtl : 

############# scent gr load + filter ##################

sig_scent = 0.1

scent_gr = regioneR::toGRanges(readRDS(paste0("output/SCENT/", 
                                              "500kb_fibroblast_allcvar_part",
                                             ".rds"
                                             )
                                       ),
                               genome = "hg38"
                               )


scent_gr_filt = scent_gr[mcols(scent_gr)$boot_basic_p < sig_scent]

############## testing enrichment significance functions ###################

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
randomizeRegions_new = function(A, ...) {
  randomizeRegions(A, genome = "hg38", ...)
}

# set up test - create function to perform bootstrapping
# different from regioneR::resampleRegions as resampleRegions
# sets replace = F
bootstrapRegions <- function(A, 
                             universe, 
                             per.chromosome=FALSE, ...) { 
  
  if(!hasArg(A)) stop("A is missing")
  if(!hasArg(universe)) stop("universe is missing")
  if(!is.logical(per.chromosome)) stop("per.chromosome must be logical")
  
  
  A <- toGRanges(A)
  universe <- toGRanges(universe)
  
  
  if(per.chromosome){
    chrResample <- function(chr) {
      Achr <- A[seqnames(A) == chr]
      universe.chr <- universe[seqnames(universe) == chr]
      resample.chr <- universe.chr[sample(1:length(universe.chr), length(Achr))]
      return(resample.chr)
    }
    
    chr.resampled <- lapply(as.list(seqlevels(A)), chrResample, replace = T)
    resampled <- do.call(c, chr.resampled)
    
  }else{
    resampled <- universe[sample(1:length(universe), length(A), replace = T)]
  }
  
  return(resampled)
  
}

############### eqtl enrichment for loop ############

# what study names to do enrichment for
eqtl_file_substring = c("commonmind_brain",
                        "rosmap_brain",
                        "geuvadis_lcl",
                        "twinsuk_lcl"
                        )

# pip value threshold to filter our casual/fine-mapped eQTLs by
sig_eqtl = 0.2

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


# bootstrapRegions = function(A, ...){
#   
#   A <- toGRanges(A)
#   resampled <- sample(A, size=length(A), replace = TRUE)
#   return(resampled)
#   
# }

############ To be fixed: #####################

# 1. should have universe as all ATAC Peaks ... 
# 2. use bootstrap regions ... 

# sig_testing = regioneR::permTest(A=eqtl_gr_filt, 
#          B=scent_gr_filt, 
#          ntimes=20, 
#          randomize.function=randomizeRegions_new, #resampleRegions, 
#         # universe=all, all ATAC 
#          evaluate.function=numOverlaps_once)


sig_testing = regioneR::permTest(A=eqtl_gr_filt,
         B=scent_gr_filt,
         ntimes=ntimes,
         randomize.function=randomizeRegions_new, #resampleRegions,
        # universe=all, all ATAC
         evaluate.function=numOverlaps_once
        )

recall_pval = sig_testing$numOverlaps_once$pval
recall_zscore = sig_testing$numOverlaps_once$zscore


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
