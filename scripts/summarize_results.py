import json
import os
import re
import csv
from utils.summarize_utils import get_tree_params, get_msa_sim_params, get_inference_params

def load_msa_stats(msa_summary_path):
    """Loads MSA stats from the pre-aggregated msa_summary.tsv into a lookup dictionary."""
    lookup = {}
    if os.path.exists(msa_summary_path):
        try:
            with open(msa_summary_path, 'r') as f:
                reader = csv.DictReader(f, delimiter='\t')
                for row in reader:
                    lookup[row["msa_dir"]] = row
        except Exception:
            pass
    return lookup

def get_msa_stats_and_link(path, msa_stats_lookup):
    parts = path.split(os.sep)
    res = {}
    try:
        inf_idx = parts.index("inference")
        msa_parts = list(parts)
        # Convert the inference path to the corresponding MSA directory path by replacing 'inference' with 'msas'
        msa_parts[inf_idx] = "msas"
        # Directory containing msa.fasta (backtracking two levels from the tool-specific inference folders)
        msa_dir_path = os.sep.join(msa_parts[:-2])
        msa_file_path = os.path.join(msa_dir_path, "msa.fasta")
        
        # Get MSA stats from lookup table
        stats = msa_stats_lookup.get(msa_dir_path, {})
        res.update(stats)
        
        # File URL
        abs_msa_path = os.path.abspath(msa_file_path)
        res["msa"] = f"[msa](file://{abs_msa_path})"
    except (ValueError, IndexError):
        pass
    return res

def get_tree_visual_link(path):
    parts = path.split(os.sep)
    try:
        inf_idx = parts.index("inference")
        tree_params = parts[inf_idx + 1]
        tree_png_path = os.path.join("results", "trees", f"{tree_params}.nwk.png")
        abs_png_path = os.path.abspath(tree_png_path)
        return {"true": f"[tree](file://{abs_png_path})"}
    except (ValueError, IndexError):
        return {}

def get_distances_from_json(dir_path):
    dist_path = os.path.join(dir_path, "distances.json")
    if os.path.exists(dist_path):
        try:
            with open(dist_path, 'r') as f:
                return json.load(f)
        except Exception:
            pass
    return {}

def get_last_line_value(file_path):
    if not os.path.exists(file_path):
        return "NA"
    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()
            if lines:
                return float(lines[-1].strip())
    except (ValueError, IndexError):
        pass
    return "NA"

def main():
    config = snakemake.params.sn_config
    msa_stats_lookup = load_msa_stats(snakemake.input.msa_summary)
    
    # We'll collect all possible fieldnames to ensure consistent TSV columns
    all_rows = []
    all_keys = set()

    for d in snakemake.params.dirs:
        row = {}
        row.update(get_tree_params(d, config["tree_match"]))
        row.update(get_msa_sim_params(d, config["msa_sim_tools"]))
        row.update(get_inference_params(d, config["inference_tools"]))
        row.update(get_msa_stats_and_link(d, msa_stats_lookup))
        row.update(get_tree_visual_link(d))
        row.update(get_distances_from_json(d))
        row["runtime_seconds"] = get_last_line_value(os.path.join(d, "time.txt"))
        row["logl"] = get_last_line_value(os.path.join(d, "logl.out"))

        # Add to collection
        all_rows.append(row)
        all_keys.update(row.keys())

    column_order = [
        # Global
        "seed", 
        # Tree Generation
        "species", "birth_rate", "death_rate", "sampling_fraction", "mutation_rate", "true",
        # MSA Simulation
        "msa_sim_tool", "root_length", "tkf_lambda", "tkf_mu", "tkf_r", "max_ins", "ir", "ip", "msa_len", "gap%", "gap_col%", "avg_gap_len", "msa",
        # Inference
        "inference_tool", "model", "gap_strategy", "logl", "runtime_seconds",
        # Results/Metrics
        "rf", "kf", "start_rf", "start_kf",
    ]
    
    # Ensure any unexpected keys are still included at the end
    final_columns = [col for col in column_order if col in all_keys]
    final_columns += [col for col in sorted(list(all_keys)) if col not in final_columns]

    write_rows_to_tsv(snakemake.output.tsv_path, all_rows, final_columns)

def write_rows_to_tsv(output_path, rows, fieldnames):
    if not rows:
        return
    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t', extrasaction='ignore', restval='NA')
        writer.writeheader()
        writer.writerows(rows)

if __name__ == "__main__":
    main()
