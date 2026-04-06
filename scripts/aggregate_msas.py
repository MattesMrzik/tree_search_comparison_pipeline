import os
import csv
from utils.fasta_utils import get_fasta_length, get_gap_stats

def main():
    config = snakemake.params.sn_config
    msa_dirs = snakemake.params.msa_dirs
    
    all_rows = []
    
    for d in msa_dirs:
        row = {"msa_dir": d}
        msa_path = os.path.join(d, "msa.fasta")
        
        if os.path.exists(msa_path):
            row["msa_len"] = get_fasta_length(msa_path)
            gap_stats = get_gap_stats(msa_path)
            row["gap%"] = gap_stats["gap_percentage"]
        else:
            row["msa_len"] = "NA"
            row["gap%"] = "NA"
            
        all_rows.append(row)
    
    # Write to TSV
    fieldnames = ["msa_dir", "msa_len", "gap%"]
    with open(snakemake.output.tsv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(all_rows)

if __name__ == "__main__":
    main()
