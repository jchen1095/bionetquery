import csv
from collections import defaultdict

data_file = "./data/cellphonedb_data/gene_input.csv"

uniprot_counts = defaultdict(list)

with open(data_file, 'r', newline='') as file:
    reader = csv.DictReader(file)
    for row in reader:
        uniprot = row['uniprot']
        gene_name = row['gene_name']
        hgnc_symbol = row['hgnc_symbol']
        ensembl_id = row['ensembl']

        uniprot_counts[uniprot].append((gene_name, hgnc_symbol, ensembl_id))

final_rows = []

for uniprot, rows in uniprot_counts.items():
    if len(rows) > 1:
        first_gene_name, first_hgnc_symbol, _ = rows[0]
        if all(gene_name == first_gene_name and hgnc_symbol == first_hgnc_symbol for gene_name, hgnc_symbol, _ in rows):
            # Keep only the first row and fetch Ensembl ID using pyensembl
            gene_name = first_gene_name
            # ensembl_id = gene_name_to_ensg_id(gene_name)
            if ensembl_id:
                final_rows.append({
                    'gene_name': gene_name,
                    'uniprot': uniprot,
                    'hgnc_symbol': first_hgnc_symbol,
                    # 'ensembl': ensembl_id
                })
        else:
            # If gene_names or hgnc_symbols are different, keep all rows as is
            for gene_name, hgnc_symbol, ensembl_id in rows:
                final_rows.append({
                    'gene_name': gene_name,
                    'uniprot': uniprot,
                    'hgnc_symbol': hgnc_symbol,
                    'ensembl': ensembl_id
                })
    else:
        # If only one row exists for the uniprot, add it directly
        gene_name, hgnc_symbol, ensembl_id = rows[0]
        final_rows.append({
            'gene_name': gene_name,
            'uniprot': uniprot,
            'hgnc_symbol': hgnc_symbol,
            'ensembl': ensembl_id
        })

# Write the final rows to a new CSV file or process them as needed
output_file = "./data/cellphonedb_data/processed_gene_input.csv"
fieldnames = ['gene_name', 'uniprot', 'hgnc_symbol', 'ensembl']

with open(output_file, 'w', newline='') as file:
    writer = csv.DictWriter(file, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(final_rows)

print(f"Processed data written to {output_file}")
