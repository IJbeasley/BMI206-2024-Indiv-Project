
# combine eqtl ... 


############# sorting common mind eqtl #############
{
common_mind = dtplyr::lazy_dt(
              data.table::fread("data/eqtl/commonmind_brain_credible_sets.tsv")
              )


common_mind = common_mind |>
              tidyr::separate(col = "variant", 
                               sep  = "_",
                                          into = c("chrom", 
                                                    "pos", 
                                                    "ref", 
                                                    "alt")
                                          )

common_mind = common_mind |> 
  dplyr::rename_with(.cols = !c("molecular_trait_id", 
                                "gene_id",
                                "chrom",
                                "rsid",
                                "pos"),
                     .fn = ~paste0("common_", .x))

common_mind = common_mind |> 
              dplyr::mutate(end = pos) |>
              dplyr::rename(start = pos)

common_mind = common_mind |>
              dplyr::relocate(end) |>
              dplyr::relocate(start) |>
              dplyr::relocate(chrom)

common_mind = as.data.frame(common_mind)

common_mind_gr = regioneR::toGRanges(common_mind,
                                     genome = "hg38")


saveRDS(common_mind_gr, paste0("output/eqtl/", 
                               "commonmind_brain",
                               "credible_sets.rds"
                               ))
}


############## sorting rosmap eqtl ##################

{
  
  rosmap = dtplyr::lazy_dt(
    data.table::fread("data/eqtl/rosmap_brain_credible_sets.tsv")
  )  
  
rosmap = rosmap |>
         tidyr::separate(col = "variant", 
                         sep  = "_",
                         into = c("chrom", 
                                  "pos", 
                                  "ref", 
                                  "alt")
                         )

rosmap = rosmap |> 
  dplyr::rename_with(.cols = !c("molecular_trait_id", 
                                "gene_id",
                                "chrom",
                                "rsid",
                                "pos"),
                     .fn = ~paste0("rosmap_", .x))

}

combined = dplyr::full_join(common_mind, 
                            rosmap)

combined = combined |> 
           dplyr::mutate(
                  comparable_vars = 
                  dplyr::if_else(
                                 common_ref == rosmap_ref & common_alt == rosmap_alt |
                                 common_ref == rosmap_alt & common_alt == rosmap_ref,
                                 1,
                                 0
                                 )
                         )

combined = as.data.frame(combined)

nrow(combined)

rosmap = as.data.frame(rosmap)

nrow(rosmap)

common_mind = as.data.frame(common_mind)

nrow(common_mind)
