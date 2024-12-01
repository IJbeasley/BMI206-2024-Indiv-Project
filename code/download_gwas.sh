

dataset_names=("michalek_t1d"
               "sakaue_t1d"
               )

urls=("https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90432001-GCST90433000/GCST90432067/harmonised/GCST90432067.h.tsv.gz"
      "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90018001-GCST90019000/GCST90018705/harmonised/34594039-GCST90018705-EFO_0001359.h.tsv.gz"
      )
      
for i in "${!urls[@]}"; do
    url="${urls[$i]}"
    dataset_name="${dataset_names[$i]}"
    
    echo " "
    echo "Downloading: $dataset_name"
    echo "With url: $url"
    echo " "
    	

    
    curl -C - -s -o  "data/gwas/${dataset_name}_full_ss.tsv.gz" "$url"
    
done

for file in data/gwas/*.gz; do
    if [ -f "$file" ]; then  # Check if the file exists
        echo " "
        echo "Unzipping $file..."
        gunzip "$file"
        #gunzip "data/gwas/${dataset_name}_full_ss.tsv.gz" 
    fi
done
