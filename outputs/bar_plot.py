import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import ast
import numpy as np

gt_genes = 'sorted_symbols.txt'
model_genes = 'processed_genes.txt'

def read_genes_from_file(file_path):
    with open(file_path, 'r') as file:
        return [line.strip() for line in file]

def read_txt(filepath): 
    print("Loading TXT data...")
    # read txt file
    with open(filepath, 'r') as file:
        lines = file.readlines()

    #convert each line from a string to a list
    list_of_lists = [ast.literal_eval(line.strip()) for line in lines if line.strip()]
    return list_of_lists

#reading data
ranked_gt = list(read_genes_from_file(gt_genes))
ranked_model_lists = read_txt(model_genes)

#Finding the set of genes from all 50 lists that intersect with the GT. In other words, looking through each list and appending an intersected gene.
common_genes = set(ranked_gt)
final_set = set()

for model_list in ranked_model_lists:
    curr_common_genes = common_genes.intersection(model_list)
    final_set.update(curr_common_genes)

#the genes in GT that are not represented in ANY of the 50 lists
missing_genes = common_genes.difference(final_set)

print(f"Number of genes in ground truth not found in model predictions: {len(missing_genes)}")
if missing_genes:
    print("Genes not found in model predictions:", missing_genes)

#Finding the avg number of genes per list (aroudn 300)
average_genes_per_list = sum(len(lst) for lst in ranked_model_lists) / len(ranked_model_lists)
print(f"Average number of genes per list: {average_genes_per_list:.2f}")

#Create dictionary for ranking ground truth by whatever metric is there
ground_truth_ranks = {gene: rank for rank, gene in enumerate(ranked_gt) if gene in common_genes}

#dictionary of gene rankings throughout all the list per gene
model_ranks = {gene: [] for gene in common_genes}
for model_list in ranked_model_lists:
    for rank, gene in enumerate(model_list):
        if gene in common_genes:
            model_ranks[gene].append(rank)

#computing average rank per gene
average_model_ranks = {gene: sum(ranks) / len(ranks) if ranks else np.nan for gene, ranks in model_ranks.items()}
min_model_ranks = {gene: min(ranks) if ranks else np.nan for gene, ranks in model_ranks.items()}
max_model_ranks = {gene: max(ranks) if ranks else np.nan for gene, ranks in model_ranks.items()}

#NORMALIZING
def normalize_ranks(ranks, scale=len(final_set)):
    valid_ranks = {gene: rank for gene, rank in ranks.items() if not np.isnan(rank)}
    max_rank = max(valid_ranks.values())
    return {gene: (rank / max_rank) * scale for gene, rank in valid_ranks.items()}

normalized_average_model_ranks = normalize_ranks(average_model_ranks)
normalized_min_model_ranks = normalize_ranks(min_model_ranks)
normalized_max_model_ranks = normalize_ranks(max_model_ranks)

# Filter out genes with NaN ranks in model predictions
valid_genes = [gene for gene in common_genes if gene in normalized_average_model_ranks]

#sort genes by ground truth ranking
valid_genes_sorted = sorted(valid_genes, key=lambda g: ground_truth_ranks[g])

#df for plotting
data = {
    'Gene': valid_genes_sorted,
    'Ground Truth Rank': [ground_truth_ranks[gene] for gene in valid_genes_sorted],
    'Average Model Rank': [normalized_average_model_ranks[gene] for gene in valid_genes_sorted],
    'Min Model Rank': [normalized_min_model_ranks[gene] for gene in valid_genes_sorted],
    'Max Model Rank': [normalized_max_model_ranks[gene] for gene in valid_genes_sorted]
}

df = pd.DataFrame(data)

# Plot
plt.figure(figsize=(14, 8))
line_color = 'b'
dot_color = 'r'
error_color = 'k'
alpha_value = 0.5  
dot_alpha = 1  

# Set positions for overlapping bars
x_positions = np.arange(len(df))

# Plot Ground Truth Rank as a line
plt.plot(x_positions, df['Ground Truth Rank'], color=line_color, marker='o', linestyle='-', linewidth=2, markersize=5, label='Ground Truth Rank')

# Plot Average Model Rank with dots and whiskers
plt.scatter(x_positions, df['Average Model Rank'], color=dot_color, alpha=dot_alpha, label='Average Model Rank', zorder=5)
plt.errorbar(x_positions, df['Average Model Rank'], yerr=[df['Average Model Rank'] - df['Min Model Rank'], df['Max Model Rank'] - df['Average Model Rank']], fmt='none', ecolor=error_color, alpha=0.7, capsize=5)

plt.xlabel('Gene')
plt.ylabel('Rank (Normalized)')
plt.title('Comparison of Ground Truth and Model Predicted Gene Ranks')
plt.legend()
plt.xticks(x_positions, df['Gene'], rotation=90)
plt.grid(axis='y')
plt.tight_layout()
plt.show()

# Show data for specific genes
genes_of_interest = ['GPR171', 'PDCD4']
print(df[df['Gene'].isin(genes_of_interest)])
