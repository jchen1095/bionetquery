import csv
from collections import defaultdict

data_file = "./data/cellphonedb_data/gene_input.csv"
uniprot_counts = defaultdict(lambda: {'count': 0, 'gene_names': set()})

with open(data_file, 'r', newline='') as file:
    reader = csv.DictReader(file)
    for row in reader:
        uniprot = row['uniprot']
        gene_name = row['gene_name']
        uniprot_counts[uniprot]['count'] += 1
        uniprot_counts[uniprot]['gene_names'].add(gene_name)

duplicate_uniprots = []
for uniprot, data in uniprot_counts.items():
    if data['count'] > 1 and len(data['gene_names']) > 1:
        duplicate_uniprots.append(uniprot)

output_file = "./data/cellphonedb_data/uniprot_mult.csv"
with open(output_file, 'w', newline='') as csvfile:
    fieldnames = ['UniProt Symbol']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    
    writer.writeheader()
    for uniprot in duplicate_uniprots:
        writer.writerow({'UniProt Symbol': uniprot})  

print(f"UniProt symbols with different gene_names and appearing more than once have been written to '{output_file}'.")
