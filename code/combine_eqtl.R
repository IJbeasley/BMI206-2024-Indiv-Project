
# combine eqtl ... 

rows=$(wc -l < data/eqtl/commonmind_brain_credible_sets.tsv)
columns=$(head -n 1 data/eqtl/commonmind_brain_credible_sets.tsv | awk '{print NF}')
echo "Rows: $rows, Columns: $columns"

common_mind = data.table::fread("data/eqtl/commonmind_brain_credible_sets.tsv") |>
              tidyr::separate()
rosmap = data.table::fread("data/eqtl/rosmap_brain_credible_sets.tsv")