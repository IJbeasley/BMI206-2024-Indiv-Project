
library(ggplot2)

############### eQTLs ########################################

# Prepare data for enrichment plot 
# Relative enrichment of eQTLs
{
eqtl_enrich = data.table::fread("output/enrichment/recall_eqtl.tsv")

names(eqtl_enrich) = 
  c("study",
      "n_recall",
      "recall",
      "recall_zscore",
      "recall_pval",
      "n_eqtl",
     "n_sig_scent",
      "sig_scent",
      "sig_eqtl",
      "ntimes"
  )

eqtl_name_fix = function(eqtl_study) {
  dplyr::case_when(eqtl_study == "commonmind_brain" ~ "CommonMind",
                   eqtl_study == "rosmap_brain" ~ "ROSMAP",
                   eqtl_study == "geuvadis_lcl" ~ "GEUVADIS",
                   eqtl_study == "twinsuk_lcl" ~ "TwinsUK",
                   TRUE ~ eqtl_study)
}

eqtl_tissue_get = function(eqtl_study){
  dplyr::case_when(grepl(tolower(eqtl_study), pattern = "lcl") ~ "LCL",
                   grepl(tolower(eqtl_study), pattern = "brain") ~"Brain",
                   TRUE ~ eqtl_study)
  
}

eqtl_enrich = eqtl_enrich |>
              dplyr::mutate(eqtl_tissue = eqtl_tissue_get(study)) |>
              dplyr::mutate(study = eqtl_name_fix(study)) 

}

###################### Plot enrichment of eQTLs #####################

eqtl_enrich |>
  ggplot(aes(y = n_recall, 
             x = study, 
             group = eqtl_tissue,
             col = as.factor(sig_scent))) + 
  geom_point() + 
  labs(y = "# Recalled \n (Fine-Mapped eQTLs overlapping SCENT enhancers)",
       x = "Study") + 
  facet_wrap(~eqtl_tissue, scales = "free_x") + 
  scale_colour_brewer(name = "SCENT \nsignifcance \nthreshold", palette = "Dark2") + 
  theme_bw()


eqtl_enrich |>
  ggplot(aes(y = recall_zscore, 
             x = study, 
             group = eqtl_tissue,
             col = as.factor(sig_scent))) + 
  geom_point() + 
  labs(y = "Z-score Recalled \n (Fine-Mapped eQTLs overlapping SCENT enhancers)",
       x = "Study") + 
  facet_wrap(~eqtl_tissue, scales = "free_x") + 
  scale_colour_brewer(name = "SCENT \nsignifcance \nthreshold", palette = "Dark2") + 
  theme_bw()


###################### Prepare Data for gwas ########################

gwas_enrich = data.table::fread("output/enrichment/recall_gwas.tsv")

names(gwas_enrich) = 
  c(
    "study",
    "trait",
    "authors",
    "cohort",
    "n_recall",
    "recall",
    # "recall_zscore",
    # "recall_pval",
    "n_gwas",
    "n_sig_scent",
    "sig_scent",
    "ntimes"
  )


gwas_enrich |>
  ggplot(aes(y = n_recall, 
             x = study, 
             #group = eqtl_tissue,
             col = as.factor(sig_scent))) + 
  geom_point() + 
  labs(y = "# Recalled \n (Bonferroni-significant overlapping SCENT enhancers)",
       x = "Study") + 
  facet_wrap(~trait, scales = "free_x") + 
  scale_colour_brewer(name = "SCENT \nsignifcance \nthreshold", palette = "Dark2") + 
  theme_bw()
