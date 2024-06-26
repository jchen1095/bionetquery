import csv
import os 

# Specify the folder containing your CSV files
folder_path = "/Users/JenChen/Desktop/SIBMI/bionetquery/hubmap_asct_data"

# Get all files in the folder (including hidden files)
csv_files = [os.path.join(folder_path, filename) for filename in os.listdir(folder_path) if filename.endswith('.csv')]

fieldnames = []
    
for filename in csv_files:
    with open(filename, "r", newline="") as f_in:
        reader = csv.reader(f_in)
        headers = next(reader)
        for h in headers:
          if h not in fieldnames:
            fieldnames.append(h)
    
    # Then copy the data
with open("out.csv", "w", newline="") as f_out:
    writer = csv.DictWriter(f_out, fieldnames=fieldnames)
    writer.writeheader() #this is the addition.       
    for filename in csv_files:
        with open(filename, "r", newline="") as f_in:
            reader = csv.DictReader(f_in)  # Uses the field names in this file
            for line in reader:
                writer.writerow(line)