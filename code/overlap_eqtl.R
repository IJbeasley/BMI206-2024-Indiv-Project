# convert scent table: 

# overlap peaks: 

scent_file = data.table::fread("data/SCENT/500kb_fibroblast_allcvar_part.txt")

scent_file = scent_file |> 
             dplyr::filter(boot_basic_p < 0.1)

scent_file = scent_file |> 
  tidyr::separate(col = "peak", 
                  sep  = ":",
                  into = c("chrom", 
                           "pos")
  ) |> 
  tidyr::separate(col = "pos", 
                  sep  = "-",
                  into = c("start", 
                           "end")
  )

scent_file = scent_file |>
             dplyr::relocate(end) |>
             dplyr::relocate(start) |>
             dplyr::relocate(chrom)

scent_gr = regioneR::toGRanges(scent_file, 
                               genome = "hg38")

#scent_gr_filt = scent_gr[mcols(scent_gr)$]

eqtl_gr = regioneR::toGRanges(readRDS(paste0("output/eqtl/", 
                                    "commonmind_brain",
                                    "credible_sets.rds")
                            ),
                    genome = "hg38")

eqtl_gr_filt = eqtl_gr[mcols(eqtl_gr)$common_pip > 0.2]

regioneR::numOverlaps(scent_gr, eqtl_gr, count.once = T) 
# eqtl recall ... 
regioneR::numOverlaps(eqtl_gr_filt, scent_gr, count.once = T) 
length(eqtl_gr_filt)

regioneR::overlapRegions(eqtl_gr,
                         scent_gr)

regioneR::overlapRegions(scent_gr,
                         eqtl_gr)
