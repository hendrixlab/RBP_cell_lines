#!/bin/bash
BASE=/nfs6/BB/Hendrix_Lab/RBPs/RBP_cell_lines

cd $BASE

hqsub -q '*' -r run_rbp_expression_no_eclip "uv run --no-project --with biopython --with numpy python rbp_list_untested_sr.py $BASE/July_2025/uniprotkb_GO_0003723_NOT_GO_0003735_NOT_2025_07_04.tsv $BASE/July_2025/idmapping_2025_07_04.tsv $BASE/rna_celline.tsv $BASE/rna_tissue_consensus.tsv $BASE/tested_rbp_list.txt"
