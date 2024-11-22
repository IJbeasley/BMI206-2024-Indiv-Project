

dataset_names=("michalek_t1d"
              "sakaue_t1d"
)

urls=("ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90432001-GCST90433000/GCST90432067/harmonised/GCST90432067.h.tsv.gz"
      "ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90018001-GCST90019000/GCST90018705/harmonised/GCST90432067.h.tsv.gz")
      
for i in "${!urls[@]}"; do
    url="${urls[$i]}"
    
    wget -nc -O "data/gwas/${dataset_name}_full_ss.tsv.gz" "$url"
    
    gunzip "data/gwas/${dataset_name}_full_ss.tsv.gz" 
    
done