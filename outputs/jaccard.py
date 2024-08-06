import csv
import ast

def read_symbols_from_file(filename):
    with open(filename, 'r') as f:
        symbols = [line.strip() for line in f.readlines()]
    return symbols

def load_gene_list_from_csv(csv_file):
    gene_list = []
    with open(csv_file, newline='') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            gene_list.extend(row)
    return gene_list

def jaccard_index(list1, list2):
    set1 = set(list1)
    set2 = set(list2)
    
    intersection = set1.intersection(set2)
    union = set1.union(set2)
    
    intersection_size = len(intersection)
    print(intersection_size)
    
    union_size = len(union)
    print(union_size)
    
    if union_size == 0:
        jaccard_index = 0.0  # Prevent division by zero
    else:
        jaccard_index = intersection_size / union_size
    
    return jaccard_index

# Example usage
file1_symbols = read_symbols_from_file('t_cell.csv')
# file2_symbols = read_symbols_from_file("related.csv")

# combined_symbols = file1_symbols + file2_symbols
c2s_output = load_gene_list_from_csv("processed_text.csv")

index = jaccard_index(c2s_output, file1_symbols)
print(f"Jaccard Index: {index}")
