import csv
import os 

folder_path = "./data/hubmap_asct_data"

csv_files = [os.path.join(folder_path, filename) for filename in os.listdir(folder_path) if filename.endswith('.csv')]

fieldnames = []
    
for filename in csv_files:
    with open(filename, "r", newline="") as f_in:
        reader = csv.reader(f_in)
        headers = next(reader)
        for h in headers:
          if h not in fieldnames:
            fieldnames.append(h)
    
   
with open("out.csv", "w", newline="") as f_out:
    writer = csv.DictWriter(f_out, fieldnames=fieldnames)
    writer.writeheader()       
    for filename in csv_files:
        with open(filename, "r", newline="") as f_in:
            reader = csv.DictReader(f_in) 
            for line in reader:
                writer.writerow(line)