
# combine eqtl ... 

eqtl_file_substring = c("commonmind_brain",
                        "rosmap_brain",
                        "geuvadis_lcl",
                        "twinsuk_lcl")

for(eqtl_file_name in eqtl_file_substring){
  eqtl_dt = dtplyr::lazy_dt(
                data.table::fread(
                                paste0("data/eqtl/",
                                      eqtl_file_name,
 #                                     commonmind_brain
                                      "_credible_sets.tsv"
                                      )
                                  )
                             )
  
  
  eqtl_dt = eqtl_dt |>
    tidyr::separate(col = "variant", 
                    sep  = "_",
                    into = c("chrom", 
                             "pos", 
                             "ref", 
                             "alt")
    )
  
  # eqtl_dt = eqtl_dt |> 
  #   dplyr::rename_with(.cols = !c("molecular_trait_id", 
  #                                 "gene_id",
  #                                 "chrom",
  #                                 "rsid",
  #                                 "pos"),
  #                      .fn = ~paste0("common_", .x))
  
  eqtl_dt = eqtl_dt  |> 
    dplyr::mutate(end = pos) |>
    dplyr::rename(start = pos)
  
  eqtl_dt = eqtl_dt |>
    dplyr::relocate(end) |>
    dplyr::relocate(start) |>
    dplyr::relocate(chrom)
  
  eqtl_dt  = as.data.frame(eqtl_dt)
  
  eqtl_gr = regioneR::toGRanges(eqtl_dt,
                                genome = "hg38")
  
  
  saveRDS(eqtl_gr, paste0("output/eqtl/", 
                          eqtl_file_name,
                          "_credible_sets.rds"
                          )
          )
  
  rm(list = ls())
  gc()
}

