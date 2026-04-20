There are two main programs/approaches:

I. Cluster-based Processing of RBP data:

Core command: 

python scripts/study_rbp_celltypes.py July_2025/uniprotkb_GO_0003723_NOT_GO_0003735_NOT_2025_07_04.tsv July_2025/idmapping_2025_07_04.tsv rna_celline.tsv rna_celline_description.tsv July_2025/DAVID_table_July_2025.txt tested_rbp_list.txt RBP_vs_transcriptType_count_HepG2.tsv

This method takes in DAVID clusters and rna cell-line expression data, generates heatmaps and other figures, outputs multiple tables in tab-delimited format. Produces lists of genes, mostly around cluster groupings with the exception of full_rpb_list.txt, which is a list of all UniProt genes


(There is a similar command for with_ribo included)



II. Print All Untested RBPs with expression:

Core command:

python rbp_list_untested_sr.py July_2025/uniprotkb_GO_0003723_NOT_GO_0003735_NOT_2025_07_04.tsv July_2025/idmapping_2025_07_04.tsv rna_celline.tsv rna_tissue_consensus.tsv tested_rbp_list.txt

This produces a list of all RBPs that have not been tested by SPIDR, or ENCORE data sets. 




1. Files downloaded for proteins:

uniprotkb_GO_0003723_AND_model_organism_2023_12_15.fasta
uniprotkb_GO_0003723_AND_model_organism_2023_12_15.tsv

downloaded from: https://www.uniprot.org/uniprotkb?query=GO%3A0003723&facets=model_organism%3A9606
On Dec 15th 2023

Found: 4230 uniprot IDs

UPDATE: downloaded again on 7/19/24, got more hits using the same link:

Found: 5,905 uniprot IDs

UPDATE2: downloaded again 1/22/25, but excluded ribosomal proteins (search = "GO:0003723 NOT GO:0003735 NOT GO:0005762 NOT GO:0030688")
downloaded from: https://www.uniprot.org/uniprotkb?query=GO%3A0003723+NOT+GO%3A0003735+NOT+GO%3A0005762+NOT+GO%3A0030688&facets=model_organism%3A9606

Found 4,035 uniprot IDs

UPDATE3: downloaded again 4/21/25 excluding ribosomal proteins: (search = "GO:0003723 NOT GO:0003735 NOT GO:0005762 NOT GO:0030688")
https://www.uniprot.org/uniprotkb?query=GO%3A0003723+NOT+GO%3A0003735+NOT+GO%3A0005762+NOT+GO%3A0030688&facets=model_organism%3A9606
Downloaded TSV, selected download all, Found 4,127
Didn't update id mapping for this update

UPDATE4: downloaded yet again on 7/4/24 same search. I also added GO terms to the TSV as I did before.
Downloaded 4,156 terms


2. Human protein atlas downloaded from here: https://v22.proteinatlas.org/about/download

Specific file downloaded (from entry 23): https://v22.proteinatlas.org/download/rna_celline.tsv.zip

Transcript expression levels summarized per gene in 1055 cell lines. The tab-separated file includes Ensembl gene identifier ("Gene"), analysed sample ("Cell line"), transcripts per million ("TPM"), protein-coding transcripts per million ("pTPM") and normalized expression ("nTPM"). The data is based on The Human Protein Atlas version 22.0 and Ensembl version 103.38.


3a. Downloaded mapping from Ensembl to Uniprot from: http://feb2021.archive.ensembl.org/biomart/martview
This particular version of Ensembl was used to match the Ensembl version used in the cell atlas data

biomart returns a full list of all uniprot IDs: Ensembl_to_UniprotKB.txt

result when compared to the 4230 uniprot IDs:

Found  2700 out of 4230

in the biomart table


3b. Download table from https://www.uniprot.org/id-mapping

result is a list of 3504 ID mapping. This method is not specifically using Ensembl103 to match the cell atlas.

Result when using this mapping

Found  3207 out of 4230

I updated this with the new data downloaded on 7/19/24 and got:

4,200 IDs were mapped to 4,568 results
1,705 ID were not mapped:

However, it looks like these 4,200 IDs were mapped to only 1,925 gene IDs

UPDATE: 1/22/25 re-downloaded using the ID maping tool, found 3,342 IDs mapped
idmapping_2025_01_23

UPDATE 2: 4/21/25 redownloaded using ID mapping tool found 3,434
idmapping_2025_04_21.tsv

update 3 not made. Going with update4 below:

UPDATE 4: 7/4/25
Extracted the IDs with awk:
less uniprotkb_GO_0003723_NOT_GO_0003735_NOT_2025_07_04.tsv | awk '{print $1}' > id_list.txt

Mapped these IDs using Uniprot's tool, got 3470 IDs mapped to Ensembl
However, it said this:

3,146 IDs were mapped to 3,470 results
1,010 ID were not mapped

So there are 3,146 distinct Ensembl genes 

STEP 4: Find Expressed genes
The expression data was incorporated to see what Ensembl IDs are unique, and which ones are expressed
Found  1559 unique ensembl out of 4156 unique uniprot
total diseases/subtypes: 29
Cell Atlas: 1455 with expression data out of 1559 total: 93.0%
104 genes with no expression data

used the following command and moved log.txt into they July_2025 folder.
python scripts/print_expressed_ensembl.py July_2025/uniprotkb_GO_0003723_NOT_GO_0003735_NOT_2025_07_04.tsv July_2025/idmapping_2025_07_04.tsv rna_celline.tsv rna_celline_description.tsv > log.txt

STEP 4: DAVID Clusters
using the expressed Ensembl IDs, I clustered them using DAVID. The cluster file is plain text and stored in the directory with the uniprot genes. 

