import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import ast
import numpy as np
import altair as alt
gt_genes = 'me_sorted_symbols_av_mac.txt'
model_genes = 'c2s_genes_av_mac.txt'

#This file produces old scatterplot
def read_genes_from_file(file_path):
    with open(file_path, 'r') as file:
        return [line.strip() for line in file]

def read_txt(filepath): 
    print("Loading TXT data...")
    with open(filepath, 'r') as file:
        lines = file.readlines()

    # Convert each line from a string to a list
    list_of_lists = [ast.literal_eval(line.strip()) for line in lines if line.strip()]
    return list_of_lists


ranked_gt = list(read_genes_from_file(gt_genes))
ranked_model_lists = read_txt(model_genes)



gt_gene_set = set(ranked_gt)
model_gene_set = set()

for model_list in ranked_model_lists:
    model_gene_set.update(model_list)

common_genes = gt_gene_set.intersection(model_gene_set)
missing_genes = gt_gene_set.difference(model_gene_set)
print(f"Number of genes in ground truth not found in at least one model prediction: {len(missing_genes)}")




cell_type = 'Alveolar Macrophage'

gt_df = pd.DataFrame(data=[], columns=['Gene', 'Rank'])
gt_df['Gene'] = ranked_gt
gt_df['Rank'] = list(range(len(ranked_gt)))
gt_df['Execution'] = 'NA'
gt_df['Cell Type'] = cell_type
gt_df['Rank Type'] = 'Ground Truth'

model_df = pd.DataFrame(index=[], data=[], columns=['Gene', 'Rank', 'Execution'])
for i, model_list in enumerate(ranked_model_lists):
    model_list = [gene for gene in model_list if gene in common_genes]
    model_list_df = pd.DataFrame(data=[], columns=['Gene', 'Rank', 'Execution'])
    model_list_df['Gene'] = model_list
    model_list_df['Rank'] = list(range(len(model_list)))
    model_list_df['Execution'] = i
    model_list_df['Cell Type'] = cell_type
    model_list_df['Rank Type'] = 'Model Prediction'
    model_df = pd.concat([model_df, model_list_df])

gt_df.to_csv('av_mac_gt_df.csv', index=False)
model_df.to_csv('av_mac_model_df.csv', index=False)




"""
print(ranked_gt)
print(ranked_model_lists)



print(f"Number of genes in ground truth not found in model predictions: {len(missing_genes)}")
if missing_genes:
    print("Genes not found in model predictions:", missing_genes)


average_genes_per_list = sum(len(lst) for lst in ranked_model_lists) / len(ranked_model_lists)
print(f"Average number of genes per list: {average_genes_per_list:.2f}")

#rank dictionaries for gt
ground_truth_ranks = {gene: rank for rank, gene in enumerate(ranked_gt) if gene in common_genes}

#rank dictionaries for model predictions
model_ranks = {gene: [] for gene in common_genes}
for model_list in ranked_model_lists:
    for rank, gene in enumerate(model_list):
        if gene in common_genes:
            model_ranks[gene].append(rank)


average_model_ranks = {gene: sum(ranks) / len(ranks) if ranks else np.nan for gene, ranks in model_ranks.items()}


differences = [ground_truth_ranks[gene] - average_model_ranks[gene] for gene in common_genes if not np.isnan(average_model_ranks[gene])]


abs_differences = np.abs(differences)

mean_difference = np.mean(abs_differences)
largest_difference = np.max(abs_differences)
std_dev_differences = np.std(abs_differences)
print(f"mean difference: {mean_difference:.2f}")
print(f"largest difference: {largest_difference:.2f}")
print(f"std dev difference: {std_dev_differences:.2f}")

data = {
    'Gene': list(common_genes),
    'Ground Truth Rank': [ground_truth_ranks[gene] for gene in common_genes],
    'Average Model Rank': [average_model_ranks[gene] for gene in common_genes]
}


df = pd.DataFrame(data)
# min_rank = df['Average Model Rank'].min(skipna=True)
# max_rank = df['Average Model Rank'].max(skipna=True)
# df['Model Rank Normalized'] = 100 * (df['Average Model Rank'] - min_rank) / (max_rank - min_rank)
sorted_model_genes = sorted(df['Gene'], key=lambda g: average_model_ranks[g])[:100]
df_top_100 = df[df['Gene'].isin(sorted_model_genes)]

plt.figure(figsize=(10, 8))
sns.scatterplot(x='Ground Truth Rank', y='Average Model Rank', data=df)
plt.plot([0, len(ranked_gt)-1], [0, 100], 'r--', lw=2)  # Diagonal line


for i in range(len(df)):
    plt.annotate(df['Gene'][i], (df['Ground Truth Rank'][i]-1, df['Average Model Rank'][i] + 1), fontsize=7, alpha=1)

plt.xlabel('Ground Truth Rank')
plt.ylabel('Average Model Rank')
plt.title('Alveolar Macrophage Rank Correlation Plot')
plt.grid(True)
plt.show()
"""