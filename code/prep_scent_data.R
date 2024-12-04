
# Prepare scent files into appropriate form
# and convert scent table to Genomic Ranges 

# what  scent  output to run this for:
scent_files_substrings = c("500kb_fibroblast_allcvar",
                           "500kb_Tcell_allcvar")

for(scent_file_name in scent_files_substrings){
  
scent_file = dtplyr::lazy_dt(
                     data.table::fread(paste0("data/SCENT/",
                                              scent_file_name,
                                             ".txt"
                                             )
                                       )
                               )


# make columns match bed format
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

# then make the bed columns at the start of the data.frame
# so they are correctly read in as chromosome, start and end of each range
scent_file = scent_file |>
  dplyr::relocate(end) |>
  dplyr::relocate(start) |>
  dplyr::relocate(chrom)

scent_file = as.data.frame(scent_file)

# now make as genomicRanges
scent_gr = regioneR::toGRanges(scent_file, 
                               genome = "hg38")

saveRDS(scent_gr,
        paste0("output/SCENT/",
                scent_file_name,
                ".rds"))

rm(list = ls())
gc()

}
