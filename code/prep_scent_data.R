
# Prepare scent files into appropriate Genomic Ranges 
# convert scent table

#scent_file = data.table::fread("data/SCENT/500kb_fibroblast_allcvar_part.txt")

scent_files_substrings = c("500kb_fibroblast_allcvar_part",
                           "500kb_Tcell_nocvar_parts")

for(scent_file_name in scent_files_substrings){
  
scent_file = dtplyr::lazy_dt(
                     data.table::fread(paste0("data/SCENT/",
                                              scent_file_name,
                                             ".txt"
                                             )
                                       )
                               )

# scent_file = scent_file |> 
#   dplyr::filter(boot_basic_p < 0.1)

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

scent_file = as.data.frame(scent_file)

scent_gr = regioneR::toGRanges(scent_file, 
                               genome = "hg38")

saveRDS(scent_gr,
        paste0("output/SCENT/",
                scent_file_name,
                ".rds"))

rm(list = ls())
gc()

}
