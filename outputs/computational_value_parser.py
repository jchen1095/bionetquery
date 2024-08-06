#Gets the computational values from cellxgene and the ranks the genes based on marker score 
import json

data_file = 'av_macrophage_ground_truth.json'

dataset_tissues = ["blood", "bone marrow", "caecum", "duodenum", "ileum", "jejunal epithelium", "lamina propria",
                   "liver", "lung", "mesenteric lymph node", "omentum", "sigmoid colon", "skeletal muscle tissue",
                    "spleen","thoracic lymph node", "thymus"," transverse colon"]
def filter_and_process_entries(entries):
    filtered_entries = {}

    for entry in entries:
        if entry["additionalMetadata"]["cellxgene_computational"]:
            computational_metadata = entry["additionalMetadata"]["cellxgene_computational"]
            marker_scores = [
                metadata["pc"] 
                for metadata in computational_metadata 
                if ("tissue_ontology_term_label" in metadata["groupby_dims"] and metadata["groupby_dims"]["organism_ontology_term_label"] == "Homo sapiens" and
                    metadata["groupby_dims"]["tissue_ontology_term_label"] in dataset_tissues)
            ]
            if marker_scores:  # Check if marker_scores is not empty
                max_score = max(marker_scores)
                min_score = min(marker_scores)
                if max_score - min_score >= 1:
                    print("large difference")
                    print(entry["symbol"])
                average_marker_score = sum(marker_scores) / len(marker_scores)

                filtered_entries[entry["symbol"]] = average_marker_score

        sorted_entries = sorted(filtered_entries.items(), key=lambda x: x[1], reverse=True)

    return sorted_entries


with open(data_file, 'r') as file:
    json_data = json.load(file)

#command call
processed_entries = filter_and_process_entries(json_data)

output_file = "pc_sorted_gt_av_mac.txt"
# Write the sorted entries to a new file
with open(output_file, 'w') as file:
    file.write(str(processed_entries))

output_file = "pc_sorted_symbols.txt"
# Write only the symbols to a new file
with open(output_file, 'w') as file:
    for entry in processed_entries:
        file.write(f"{entry[0]}\n")  # Write only the symbol


print(f"Sorted entries have been written to {output_file}")
#if there are multiple instances of cellxgene computational data then

#ALL of the tissues represented in the Dominguez dataset that C2S is trained on
# blood
# bone marrow
# caecum
# duodenum
# ileum
# jejunal epithelium
# lamina propria
# liver
# lung
# mesenteric lymph node
# omentum
# sigmoid colon
# skeletal muscle tissue
# spleen
# thoracic lymph node
# thymus
# transverse colon