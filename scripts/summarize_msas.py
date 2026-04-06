import os
import re
import csv
from utils.fasta_utils import get_fasta_length, get_gap_stats
from utils.summarize_utils import get_tree_params, get_msa_sim_params

def main():
    config = snakemake.params.sn_config
    msa_dirs = snakemake.params.msa_dirs
    
    all_rows = []
    all_keys = set()
    
    for d in msa_dirs:
        row = {"msa_dir": d}
        
        # Extract parameters from path
        row.update(get_tree_params(d, config["tree_match"]))
        row.update(get_msa_sim_params(d, config["msa_sim_tools"]))
        
        msa_path = os.path.join(d, "msa.fasta")
        if os.path.exists(msa_path):
            row["msa_len"] = get_fasta_length(msa_path)
            gap_stats = get_gap_stats(msa_path)
            row.update(gap_stats)
            
            # File URL
            abs_msa_path = os.path.abspath(msa_path)
            row["msa"] = f"[msa](file://{abs_msa_path})"
        else:
            row["msa_len"] = "NA"
            row["gap%"] = "NA"
            row["gap_col%"] = "NA"
            row["avg_gap_len"] = "NA"
            row["msa"] = "NA"
            
        all_rows.append(row)
        all_keys.update(row.keys())
    
    # Define column order (similar to summarize_results.py but focused on MSAs)
    column_order = [
        "seed", "species", "birth_rate", "death_rate", "sampling_fraction", "mutation_rate",
        "msa_sim_tool", "root_length", "tkf_lambda", "tkf_mu", "tkf_r", "max_ins", "ir", "ip",
        "msa_len", "gap%", "gap_col%", "avg_gap_len", "msa", "msa_dir"
    ]
    
    final_columns = [col for col in column_order if col in all_keys]
    final_columns += [col for col in sorted(list(all_keys)) if col not in final_columns]

    # Write to TSV
    with open(snakemake.output.tsv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=final_columns, delimiter='\t', extrasaction='ignore', restval='NA')
        writer.writeheader()
        writer.writerows(all_rows)

if __name__ == "__main__":
    main()
