import ast

gt_genes = 'av_phage.txt'
model_genes = 'processed_genes_av_mac.txt'

def read_genes_from_file(file_path):
    with open(file_path, 'r') as file:
        return [line.strip() for line in file]

def read_txt(filepath): 
    print("Loading TXT data...")
    # Read the TXT file as a text file
    with open(filepath, 'r') as file:
        lines = file.readlines()

    # Convert each line from a string to a list
    list_of_lists = [ast.literal_eval(line.strip()) for line in lines if line.strip()]
    return list_of_lists

# Read data
ranked_gt = list(read_genes_from_file(gt_genes))
ranked_model_lists = read_txt(model_genes)


# Find common genes
common_genes = set(ranked_gt)
final_set = set()

for model_list in ranked_model_lists:
    print(len(model_list))
    curr_common_genes = common_genes.intersection(model_list)
    final_set.update(curr_common_genes)

print(len(ranked_gt))

print(len(final_set))