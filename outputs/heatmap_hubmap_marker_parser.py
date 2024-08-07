import json
import ast
import numpy as np 
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.colors as mcolors



def filter_and_process_entries(entries):
    filtered_entries = []

    for entry in entries:
        if "additionalMetadata" in entry and "hubmap" in entry["additionalMetadata"]:
            filtered_entries.append(entry["symbol"])

    return filtered_entries

def process_txt_and_plot_heatmap(txt_path, gt_genes):
    print(len(gt_genes))
    print("Loading TXT data...")
    with open(txt_path, 'r') as file:
        lines = file.readlines()

    # Convert each line from a string to a list
    list_of_lists = [ast.literal_eval(line.strip()) for line in lines if line.strip()]

    print(f"Number of lists: {len(list_of_lists)}")

    final_set = set()
    intersections = []

    for model_list in list_of_lists:
        curr_common_genes = gt_genes.intersection(model_list)
        intersections.append(len(curr_common_genes))
        final_set.update(curr_common_genes)
    
    print(f"Number of intersected genes: {len(final_set)}")
    print(f"Average number of intersected genes: {np.mean(intersections)}")

    max_intersections_idx = np.argmax(intersections)
    print(f"List with the most intersections is List {max_intersections_idx + 1} with {intersections[max_intersections_idx]} intersections")

    filtered_gt_genes = [gene for gene in gt_genes if gene in final_set]

    #get the order of genes from the list with the most intersections
    best_list_order = {gene: idx for idx, gene in enumerate(list_of_lists[max_intersections_idx]) if gene in filtered_gt_genes}

    #sort the filtered_gt_genes according to the order in the best list
    sorted_filtered_gt_genes = sorted(filtered_gt_genes, key=lambda gene: best_list_order.get(gene, len(filtered_gt_genes)))

    #move the list with the most intersections to the first position
    list_of_lists.insert(0, list_of_lists.pop(max_intersections_idx))

    #df for the heatmap
    heatmap_data = np.full((len(list_of_lists), len(sorted_filtered_gt_genes)), np.nan) 

    #map gene to its index in the sorted_filtered_gt_genes list
    gene_to_index = {gene: idx for idx, gene in enumerate(sorted_filtered_gt_genes)}

    #fill the heatmap data with sorted gene order based on the best list
    for list_idx, gene_list in enumerate(list_of_lists):
        for index, gene in enumerate(gene_list):
            if gene in gene_to_index:
                sorted_idx = gene_to_index[gene]
                heatmap_data[list_idx, sorted_idx] = index + 1  # Use index + 1 for color scaling

    #replace NaNs with 0 for cells where the gene does not appear
    heatmap_data = np.nan_to_num(heatmap_data)

    
    plt.figure(figsize=(12, 8))
    
    
    colors = ["blue", "purple", "red"]
    cmap = mcolors.LinearSegmentedColormap.from_list("custom_cmap", colors)
    
    #mask non-appearances (cells with value 0)
    mask = (heatmap_data == 0)
    
    #plot
    ax = sns.heatmap(heatmap_data, cmap=cmap, annot=False, cbar_kws={'label': 'Gene Index'},
                xticklabels=sorted_filtered_gt_genes, yticklabels=[f'List {i+1}' for i in range(len(list_of_lists))],
                vmin=1, vmax=300, mask=mask)  # Apply mask for non-appearances
    
    for i, label in enumerate(ax.get_xticklabels()):
        label.set_visible(False)
    
    plt.xlabel('HuBMAP Ground Truth Alveolar Macrophage Biomarkers')
    plt.ylabel('C2S Output Instances')
    plt.title('Heatmap of Biomarker Appearances and Expression Ranking in C2S Outputs')
    plt.xticks(rotation=90)  
    plt.tight_layout()
    plt.show()

    plt.savefig('heatmap_output.png')

#groundtruth
data_file = 'av_macrophage_ground_truth.json'

dataset_tissues = ["blood", "bone marrow", "caecum", "duodenum", "ileum", "jejunal epithelium", "lamina propria",
                   "liver", "lung", "mesenteric lymph node", "omentum", "sigmoid colon", "skeletal muscle tissue",
                   "spleen", "thoracic lymph node", "thymus", "transverse colon"]

#50 lists from model
txt_filepath = "c2s_genes_av_mac.txt"

def main():
    with open(data_file, 'r') as file:
        json_data = json.load(file)

    # Get filtered and sorted entries
    processed_entries = set(filter_and_process_entries(json_data))

    print(f"Number of GT genes: {len(processed_entries)}")
    print("called")

    process_txt_and_plot_heatmap(txt_filepath, processed_entries)

main()
