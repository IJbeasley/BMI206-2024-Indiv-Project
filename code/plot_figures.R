
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
    "scent_tissue",
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


eqtl_enrich = eqtl_enrich |> 
              dplyr::mutate(precision = n_recall / n_sig_scent) |> 
              dplyr::mutate(f1 = 2 / (
                                      (1/precision) + (1/recall)
                                      )
                            )


eqtl_enrich = eqtl_enrich |>
              dplyr::filter(sig_scent < 1)

}

###################### Plot enrichment of eQTLs #####################

{
  
  pdf(file = "figures/eqtl_enrich_f1.pdf",
      width=3 *50 / 24.5,
      height=3* 30 / 24.5)
  
  
  
  plot = eqtl_enrich |>
    ggplot(aes(y = f1, 
               x = study, 
               group = scent_tissue,
               col = as.factor(scent_tissue))) + 
    geom_point(size = 5,
               position = position_dodge(0.75)) +
    ylim(0, 0.015) + 
    labs(y = "F1",
         x = "Study") + 
    facet_wrap(~eqtl_tissue, scales = "free_x") + 
    scale_colour_brewer(name =  "SCENT \npeak tissue", palette = "Dark2") + 
    theme_bw(base_size = 12) + 
    theme(       strip.text = element_text(size = rel(1.25)), 
                 legend.text = element_text(size = rel(1.15)),
                 legend.title = element_text(size = rel(1.17))
    )
  
  print(plot)
  
  
  dev.off()
  
}

############# Number Recalled ########################
{
  
  pdf(file = "figures/eqtl_enrich_n_recall.pdf",
      width=3 *50 / 24.5,
      height=3* 30 / 24.5)
  
  

plot = eqtl_enrich |>
  ggplot(aes(y = n_recall, 
             x = study, 
             group = scent_tissue,
             col = as.factor(scent_tissue))) + 
  geom_point(size = 5,
             position = position_dodge(0.75)) +
  labs(y = "# Recalled \n (Fine-mapped eQTLs)",
       x = "Study") + 
  facet_wrap(~eqtl_tissue, scales = "free_x") + 
  scale_colour_brewer(name =  "SCENT \npeak tissue", palette = "Dark2") + 
  theme_bw(base_size = 12) + 
  theme(       strip.text = element_text(size = rel(1.25)), 
               legend.text = element_text(size = rel(1.15)),
               legend.title = element_text(size = rel(1.17))
  )

print(plot)


dev.off()
  
}

################ Z-score recalled #######################

{
  
  pdf(file = "figures/eqtl_enrich_zscore_recall.pdf",
      width=3 *50 / 24.5,
      height=3* 30 / 24.5)
  
  
plot = eqtl_enrich |>
  ggplot(aes(y = recall_zscore, 
             x = study, 
             group = eqtl_tissue,
             col = as.factor(scent_tissue))) + 
  geom_point(size = 5,
             position = position_dodge(0.75)) +
  labs(y = "Z-score for # Recalled \n (Fine-mapped eQTLs)",
       x = "Study") + 
  facet_wrap(~eqtl_tissue, scales = "free_x") + 
  scale_colour_brewer(name = "SCENT \npeak tissue", palette = "Dark2") + 
  theme_bw(base_size = 12) + 
  theme(       strip.text = element_text(size = rel(1.25)), 
    legend.text = element_text(size = rel(1.15)),
                legend.title = element_text(size = rel(1.17))
                )

print(plot)

dev.off()
  
} 

################# Proportion Recalled ####################

{
  
  pdf(file = "figures/eqtl_enrich_prop_recall.pdf",
      width=3 *50 / 24.5,
      height=3* 30 / 24.5)
  
plot = eqtl_enrich |>
  ggplot(aes(y = recall, 
             x = study, 
             group = scent_tissue,
             col = as.factor(scent_tissue))) + 
  ylim(0, 0.05) + 
  geom_point(size = 5,
             position = position_dodge(0.75)) +
  labs(y = "Recall \n(Fine-mapped eQTLs)",
       x = "Study") +
  facet_wrap(~eqtl_tissue, scales = "free_x") + 
  scale_colour_brewer(name = "SCENT \npeak tissue", palette = "Dark2") + 
  theme_bw(base_size = 12) + 
  theme(       strip.text = element_text(size = rel(1.25)), 
               legend.text = element_text(size = rel(1.15)),
               legend.title = element_text(size = rel(1.17))
  )


print(plot) 


dev.off()

}

###################### Prepare Data for gwas ########################

{
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
    "scent_tissue",
    "ntimes"
  )


gwas_enrich = gwas_enrich |> 
              dplyr::mutate(trait = dplyr::if_else(trait == "Insulin-dependent diabetes mellitus",
                                                   "Type 1 diabetes",
                                                   trait)
              )

}

################# Plot gwas enrichment #####################

{
  
  pdf(file = "figures/gwas_enrich_n_recall.pdf",
      width=5 *50 / 24.5,
      height=3.5* 30 / 24.5)
  

plot  = gwas_enrich |>
  ggplot(aes(y = n_recall, 
             x = study, 
             group = scent_tissue,
             col = as.factor(scent_tissue))) + 
  geom_point(size = 5,
             position = position_dodge(0.75)) +
  # geom_jitter(size = 5,
  #        #    position = "jitter"
  #             height = 0,
  #             width = 0.5
  #            ) + 
  labs(y = "# Recalled \n (Genome-wide significant variants)",
       x = "Study") + 
  facet_wrap(~trait, scales = "free_x") + 
  scale_colour_brewer(name = "SCENT tissue", palette = "Dark2") + 
  theme_bw(
           base_size = 12) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = rel(1.35)),
        legend.text = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.25))
        )

print(plot)


dev.off()
  
  
}



