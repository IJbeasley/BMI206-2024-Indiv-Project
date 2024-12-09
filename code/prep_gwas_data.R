gwas_associations_files = c("data/gwas/gwas-catalog-associations-T1D.tsv",
                            "data/gwas/gwas-catalog-associations-RA.tsv"
                            )

study_accessions = 
c(# T1D
#  "GCST90432067",
  "GCST90454277", # Zeng - european 
  "GCST90018705", # Sakaue - east asian
  "GCST009916", # Onengut-Gumusc African
  "GCST005536",
  "GCST008377", # Zhu east asian
  "GCST90000529", #Inshaw European
  # RA
  "GCST90132222", # all Ishigaki 
  "GCST90132223", # european Ishigaki 
  "GCST90132224", # east asian Ishigaki 

  "GCST011389", #Kwon  east asian 
 # "GCST90258644",
  "GCST90018690", # Sakaue  east asian 
 "GCST90446530" # Guo
  )

for(associations_file in gwas_associations_files){

  message("\n \n Getting associations for: ", 
           associations_file)
  
  associations = data.table::fread(associations_file)


# t1d_associations  = t1d_associations |>
#                     dplyr::filter(MAPPED_TRAIT == "type 1 diabetes mellitus" )
# 
# t1d_associations =  t1d_associations |>  
#                      dplyr::filter(!!rlang::sym("DISEASE/TRAIT") == "Type 1 diabetes")

  # filter by study accession 
associations = associations |> 
                   dplyr::filter(!!rlang::sym("STUDY ACCESSION") %in% 
                                   # c("GCST90432067",
                                   #   "GCST90018705",
                                   #   "GCST009916",
                                   #   "GCST005536",
                                   #   "GCST008377",
                                   #   "GCST90000529")
                                   study_accessions
                                 )

# t1d_associations = t1d_associations |> 
#   dplyr::rename(STUDY_ACCESSION = `STUDY ACCESSION`)

# remove annoying spaces in column names
associations = associations |>
                   dplyr::rename_with(.fn = ~gsub(" ", "_", .x))

# remove snps with odd ids 
associations = associations |> 
              dplyr::filter(!grepl("esv", SNPS))

if(class(associations$CHR_POS) != "integer"){

# fix weird example of chrX:78464616
associations  <- associations |> 
  dplyr::mutate(CHR_POS = dplyr::if_else(
    SNPS == "chrX:78464616",
    "78464616",
    CHR_POS
  )
  ) |> 
  dplyr::mutate(CHR_ID = 
                  dplyr::if_else(
                    SNPS == "chrX:78464616",
                    "X",
                    CHR_ID
                               )
  )

# fix weird examples for rsIDs 60681359 & 386399572
associations  <- associations |>
  dplyr::mutate(
    CHR_POS = dplyr::case_when(
      SNP_ID_CURRENT == "60681359" ~ "69880215",
      SNP_ID_CURRENT == "386399572" ~ "26115952",
      TRUE ~ CHR_POS # Preserve original value for other cases
                              )
    ) |>
      dplyr::mutate(
        CHR_ID = dplyr::case_when(
          SNP_ID_CURRENT == "60681359" ~ "18",
          SNP_ID_CURRENT == "386399572" ~ "4",
          TRUE ~ CHR_ID # Preserve original value for other cases
        )
      )

associations$CHR_POS = as.integer(associations$CHR_POS)

}
#   dplyr::rows_update(assocations,
#                      tibble(CHR_POS = 78464616,
#                             SNP_ID_CURRENT = "chrX:78464616"
#                             ), 
#                     by = "SNP_ID_CURRENT"
#                     )
# #chrX:78464616
# 
# assocations = assocations |>
# dplyr::rows_update(tibble(CHR_POS = c(69880215,
#                                      26115952),
#                           SNP_ID_CURRENT = c(60681359, 
#                                       386399572
#                           )
#                           ), 
#                    by = "SNP_ID_CURRENT")
# t1d_associations = t1d_associations |> 
#                    dplyr::rename(!!rlang::sym("STUDY ACCESSION") == STUDY_ACCESSION)

associations = associations |> 
               dplyr::distinct()

message("\n Number of associations per selected study")
print(
associations |>
dplyr::group_by(STUDY_ACCESSION) |>
dplyr::summarise(n_associations = n())
)

# format into genomic ranges object
# by renaming and reordering required columns

associations = associations |> 
  dplyr::mutate(end = CHR_POS) |>
  dplyr::rename(start = CHR_POS,
                chrom = CHR_ID
                )

associations = associations |>
dplyr::relocate(end) |>
  dplyr::relocate(start) |>
  dplyr::relocate(chrom)


associations_gr = regioneR::toGRanges(associations,
                                      genome = "hg38"
)

#print(length(associations_gr))

# now save genomic ranges object: as output
output = stringr::str_replace(associations_file,
                              pattern = ".tsv",
                              replacement = ".rds"
)

output = stringr::str_replace(output,
                              pattern = "data",
                              replacement = "output")


saveRDS(associations_gr, 
        output
        )

message("\n Genomic Ranges Object saved as: ",
        output)

} 
