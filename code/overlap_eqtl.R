# overlap peaks: 

# scent_file = data.table::fread("data/SCENT/500kb_fibroblast_allcvar_part.txt")
# 
# scent_file = scent_file |> 
#              dplyr::filter(boot_basic_p < 0.1)
# 
# scent_file = scent_file |> 
#   tidyr::separate(col = "peak", 
#                   sep  = ":",
#                   into = c("chrom", 
#                            "pos")
#   ) |> 
#   tidyr::separate(col = "pos", 
#                   sep  = "-",
#                   into = c("start", 
#                            "end")
#   )
# 
# scent_file = scent_file |>
#              dplyr::relocate(end) |>
#              dplyr::relocate(start) |>
#              dplyr::relocate(chrom)
# 
# scent_gr = regioneR::toGRanges(scent_file, 
#                                genome = "hg38")

#scent_gr_filt = scent_gr[mcols(scent_gr)$]

# load scent gr: 

sig_scent = 0.05

scent_gr = regioneR::toGRanges(readRDS(paste0("output/SCENT/", 
                                              "500kb_fibroblast_allcvar_part",
                                             ".rds")
),
genome = "hg38")


scent_gr_filt = scent_gr[mcols(scent_gr)$boot_basic_p < sig_scent]

eqtl_file_substring = c("commonmind_brain",
                        "rosmap_brain",
                        "geuvadis_lcl",
                        "twinsuk_lcl")

# loading eqtl gr: 


for(eqtl_file_name in eqtl_file_substring){

eqtl_gr = regioneR::toGRanges(readRDS(paste0("output/eqtl/", 
                                             eqtl_file_name,
                                    "_credible_sets.rds")
                            ),
                    genome = "hg38")
# filter eqtl gr 

eqtl_gr_filt = eqtl_gr[mcols(eqtl_gr)$pip > 0.9]

#regioneR::numOverlaps(scent_gr, eqtl_gr, count.once = T) 
# eqtl recall ... 
overlap = regioneR::numOverlaps(eqtl_gr_filt, 
                      scent_gr_filt, 
                      count.once = T) 
n_eqtl = length(eqtl_gr_filt)

n_sig_scent = length(scent_gr_filt)

message("\n Recalled: ", overlap, 
        "\n N eQTL: ", n_eqtl, 
        "\n N scent: ", n_sig_scent,
        "\n"
        )

}
# regioneR::overlapRegions(eqtl_gr,
#                          scent_gr)
# 
# regioneR::overlapRegions(scent_gr,
#                          eqtl_gr)
