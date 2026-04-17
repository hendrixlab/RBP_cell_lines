# Usage: run_rbp_list_untested_sr.sh

import sys
import os
import numpy as np
from Bio import SeqIO


def main():
    usage = ("Usage: " + sys.argv[0] +
             " <uniprot_table> <id_mapping> <cell_line_expression>"
             " <tissue_expression> <tested_rbp_list>")
    if len(sys.argv) != 6:
        print(usage)
        sys.exit(1)

    uniprot_table_file = sys.argv[1]
    id_mapping_file = sys.argv[2]
    cell_line_exp_file = sys.argv[3]
    tissue_exp_file = sys.argv[4]
    tested_rbp_file = sys.argv[5]

    # Step 1: Load inputs
    print("Reading UniProt table and FASTA...")
    rbp_uniprot, id_map, protein_name, sequences, reviewed_status = read_uniprot_table_file(uniprot_table_file)

    print("Reading UniProt to Ensembl mapping...")
    uniprot_to_ensembl, ensembl_to_uniprot = read_uniprot_to_ensembl_file(id_mapping_file)

    print("Reading cell line expression data...")
    cell_atlas, cellline_list = read_cell_atlas_file(cell_line_exp_file)

    print("Reading tissue expression data...")
    tissue_atlas, tissue_list = read_tissue_atlas_file(tissue_exp_file)

    print("Reading tested RBPs...")
    tested_rbps = read_tested_rbps(tested_rbp_file)

    # Step 2: Map UniProt -> Ensembl
    print("Mapping UniProt to Ensembl...")
    rbp_ensembl_set, ensembl_to_uniprot_id = map_to_ensembl(
        rbp_uniprot, uniprot_to_ensembl, reviewed_status
    )

    # Step 3: Compute expression stats
    print("Computing cell line expression stats...")
    cellline_stats = compute_expression_stats(cell_atlas, cellline_list, rbp_ensembl_set)

    print("Computing tissue expression stats...")
    tissue_stats = compute_expression_stats(tissue_atlas, tissue_list, rbp_ensembl_set)

    # Step 4: Build records, filter, write
    print("Building records and filtering tested RBPs...")
    records = build_records(
        rbp_ensembl_set, ensembl_to_uniprot_id, id_map, protein_name,
        sequences, cellline_stats, tissue_stats, tested_rbps
    )

    output_file = "rbp_list_untested.tsv"
    print(f"Writing {len(records)} RBPs to {output_file}...")
    write_compendium_tsv(output_file, records)
    print("Done.")


#####################
# INPUT SUBROUTINES #
#####################

def read_uniprot_table_file(uniprot_table_file):
    """Parse UniProt TSV and paired FASTA file."""
    rbp_uniprot = []
    protein_name = {}
    id_map = {}
    reviewed_status = {}
    with open(uniprot_table_file) as F:
        for line in F:
            if not line.startswith('Entry'):
                terms = line.strip().split('\t')
                uniprot_id = terms[0]
                reviewed = terms[1]
                this_protein_name = terms[3]
                gene_names = terms[4]
                id_map[uniprot_id] = gene_names
                protein_name[uniprot_id] = this_protein_name
                reviewed_status[uniprot_id] = reviewed
                rbp_uniprot.append(uniprot_id)
    sequences = {}
    fasta_file = os.path.splitext(uniprot_table_file)[0] + '.fasta'
    if os.path.exists(fasta_file):
        with open(fasta_file) as FASTA:
            for record in SeqIO.parse(FASTA, 'fasta'):
                terms = record.id.split('|')
                uniprot_id = terms[1]
                sequences[uniprot_id] = str(record.seq)
    else:
        print(fasta_file, 'could not be found!')
        sys.exit(1)
    return rbp_uniprot, id_map, protein_name, sequences, reviewed_status


def read_uniprot_to_ensembl_file(uniprot_to_ensembl_file):
    """Parse ID mapping file (UniProt accession -> Ensembl gene ID)."""
    uniprot_to_ensembl = {}
    ensembl_to_uniprot = {}
    with open(uniprot_to_ensembl_file) as F:
        for line in F:
            if not line.startswith('From') and not line.startswith('Gene stable ID'):
                terms = line.strip().split('\t')
                if len(terms) == 2:
                    uniprot_acc, ensemble_id_plus = terms
                    ensemble_id = ensemble_id_plus.split('.')[0]
                    uniprot_to_ensembl[uniprot_acc] = ensemble_id
                    ensembl_to_uniprot[ensemble_id] = uniprot_acc
    return uniprot_to_ensembl, ensembl_to_uniprot


def read_cell_atlas_file(cell_atlas_file):
    """Parse HPA cell line expression TSV.
    Columns: Gene, Gene name, Cell line, TPM, pTPM, nTPM
    """
    cell_atlas = {}
    celltype_set = set()
    with open(cell_atlas_file) as F:
        for line in F:
            if not line.startswith('Gene'):
                terms = line.strip().split('\t')
                ensemble_id = terms[0]
                celltype = terms[2]
                TPM = float(terms[3])
                if celltype not in cell_atlas:
                    cell_atlas[celltype] = {}
                cell_atlas[celltype][ensemble_id] = TPM
                celltype_set.add(celltype)
    celltype_list = list(celltype_set)
    return cell_atlas, celltype_list


def read_tissue_atlas_file(tissue_atlas_file):
    """Parse HPA tissue consensus expression TSV.
    Columns: Gene, Gene name, Tissue, nTPM
    """
    tissue_atlas = {}
    tissue_set = set()
    with open(tissue_atlas_file) as F:
        for line in F:
            if not line.startswith('Gene'):
                terms = line.strip().split('\t')
                ensemble_id = terms[0]
                tissue = terms[2]
                nTPM = float(terms[3])
                if tissue not in tissue_atlas:
                    tissue_atlas[tissue] = {}
                tissue_atlas[tissue][ensemble_id] = nTPM
                tissue_set.add(tissue)
    tissue_list = list(tissue_set)
    return tissue_atlas, tissue_list


def read_tested_rbps(tested_rbp_file):
    """Load tested RBP gene symbols (one per line)."""
    tested_rbps = set()
    for line in open(tested_rbp_file):
        rbp_id = line.strip()
        if rbp_id:
            tested_rbps.add(rbp_id)
    return tested_rbps


######################
# PROCESSING LOGIC   #
######################

def map_to_ensembl(rbp_uniprot, uniprot_to_ensembl, reviewed_status):
    """Map UniProt IDs to Ensembl, deduplicating by Ensembl ID.
    Prefer 'reviewed' (Swiss-Prot) over 'unreviewed' (TrEMBL).
    """
    ensembl_to_uniprot_id = {}
    mapped_count = 0
    for uniprot_id in rbp_uniprot:
        if uniprot_id in uniprot_to_ensembl:
            mapped_count += 1
            ensembl_id = uniprot_to_ensembl[uniprot_id]
            if ensembl_id not in ensembl_to_uniprot_id:
                ensembl_to_uniprot_id[ensembl_id] = uniprot_id
            else:
                existing = ensembl_to_uniprot_id[ensembl_id]
                if (reviewed_status.get(uniprot_id) == 'reviewed' and
                        reviewed_status.get(existing) != 'reviewed'):
                    ensembl_to_uniprot_id[ensembl_id] = uniprot_id

    rbp_ensembl_set = set(ensembl_to_uniprot_id.keys())
    print(f"Mapped {mapped_count} UniProt IDs to {len(rbp_ensembl_set)} unique Ensembl genes")
    return rbp_ensembl_set, ensembl_to_uniprot_id


def compute_expression_stats(atlas_dict, context_list, rbp_ensembl_set):
    """Compute per-RBP expression statistics across contexts."""
    stats = {}
    for ensembl_id in rbp_ensembl_set:
        values = []
        max_context = None
        max_value = -1.0
        for context in context_list:
            if context in atlas_dict and ensembl_id in atlas_dict[context]:
                val = atlas_dict[context][ensembl_id]
                values.append(val)
                if val > max_value:
                    max_value = val
                    max_context = context
        if len(values) > 0:
            avg_value = np.mean(values)
            avg_log_value = np.mean([np.log10(v + 1.0) for v in values])
            stats[ensembl_id] = (avg_value, avg_log_value, max_context, max_value)
    return stats


def build_records(rbp_ensembl_set, ensembl_to_uniprot_id, id_map, protein_name,
                  sequences, cellline_stats, tissue_stats, tested_rbps):
    """Build output records, filtering out tested RBPs and those with no expression."""
    records = []
    tested_count = 0
    no_expression_count = 0

    for ensembl_id in rbp_ensembl_set:
        uniprot_id = ensembl_to_uniprot_id[ensembl_id]

        gene_names_str = id_map.get(uniprot_id, '')
        gene_names = set(gene_names_str.split())
        if gene_names.intersection(tested_rbps):
            tested_count += 1
            continue

        has_cellline = ensembl_id in cellline_stats
        has_tissue = ensembl_id in tissue_stats
        if not has_cellline and not has_tissue:
            no_expression_count += 1
            continue

        short_name = gene_names_str.split()[0] if gene_names_str else 'NA'

        record = {
            'ensembl_id': ensembl_id,
            'uniprot_id': uniprot_id,
            'protein_name': short_name,
            'protein_description': protein_name.get(uniprot_id, 'NA'),
        }

        if has_cellline:
            avg_val, avg_log, max_ctx, max_val = cellline_stats[ensembl_id]
            record['avg_cellline_tpm'] = f"{avg_val:.2f}"
            record['avg_cellline_log_tpm'] = f"{avg_log:.4f}"
            record['max_cellline'] = max_ctx
            record['max_cellline_tpm'] = f"{max_val:.2f}"
        else:
            record['avg_cellline_tpm'] = 'NA'
            record['avg_cellline_log_tpm'] = 'NA'
            record['max_cellline'] = 'NA'
            record['max_cellline_tpm'] = 'NA'

        if has_tissue:
            avg_val, avg_log, max_ctx, max_val = tissue_stats[ensembl_id]
            record['avg_tissue_ntpm'] = f"{avg_val:.2f}"
            record['avg_tissue_log_ntpm'] = f"{avg_log:.4f}"
            record['max_tissue'] = max_ctx
            record['max_tissue_ntpm'] = f"{max_val:.2f}"
        else:
            record['avg_tissue_ntpm'] = 'NA'
            record['avg_tissue_log_ntpm'] = 'NA'
            record['max_tissue'] = 'NA'
            record['max_tissue_ntpm'] = 'NA'

        record['amino_acid_sequence'] = sequences.get(uniprot_id, 'NA')

        records.append(record)

    def sort_key(r):
        try:
            return -float(r['avg_cellline_tpm'])
        except ValueError:
            return 0.0
    records.sort(key=sort_key)

    print(f"Filtered out {tested_count} tested RBPs")
    print(f"Excluded {no_expression_count} RBPs with no expression data")
    return records


def write_compendium_tsv(output_file, records):
    """Write output TSV with header."""
    columns = [
        'ensembl_id', 'uniprot_id', 'protein_name', 'protein_description',
        'avg_cellline_tpm', 'avg_cellline_log_tpm', 'max_cellline', 'max_cellline_tpm',
        'avg_tissue_ntpm', 'avg_tissue_log_ntpm', 'max_tissue', 'max_tissue_ntpm',
        'amino_acid_sequence'
    ]
    with open(output_file, 'w') as OUT:
        OUT.write('\t'.join(columns) + '\n')
        for record in records:
            OUT.write('\t'.join(record[col] for col in columns) + '\n')


if __name__ == '__main__':
    main()
