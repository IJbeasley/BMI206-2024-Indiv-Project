
### bash code 

### download eqtl summary stats: 

study_ids=("QTS000008" "QTS000025" "QTS000013" "QTS000029")
dataset_ids=("QTD000075" "QTD000434" "QTD000110" "QTD000539")
dataset_names=("commonmind_brain" "rosmap_brain" "geuvadis_lcl" "twinsuk_lcl")

# Loop over each element in the arrays
for i in "${!study_ids[@]}"; do
    study_id="${study_ids[$i]}"
    dataset_id="${dataset_ids[$i]}"
    dataset_name="${dataset_names[$i]}"
    
    # credible intervals:
    wget -nc -c -q -nc -O "data/eqtl/${dataset_name}_credible_sets.tsv.gz" "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/susie/$study_id/$dataset_id/$dataset_id.credible_sets.tsv.gz"
    
    # lbf: 
   # wget -nc -c -q -nc -O "data/eqtl/${dataset_name}_lbf_variable.txt.gz" "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/susie/$study_id/$dataset_id/$dataset_id.lbf_variable.txt.gz"
  
done

for file in data/eqtl/*.gz; do
    if [ -f "$file" ]; then  # Check if the file exists
        echo "Unzipping $file..."
        gunzip "$file"
    fi
done
