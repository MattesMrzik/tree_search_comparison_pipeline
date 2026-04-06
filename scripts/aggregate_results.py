import json
import os
import re
import csv

def main():
    config = snakemake.params.sn_config
    # Generic regexes for shared path components
    tree_regex = config["tree_match"]
    
    # Cache for MSA stats to avoid redundant I/O
    msa_cache = {}

    # We'll collect all possible fieldnames to ensure consistent TSV columns
    all_rows = []
    all_keys = set()

    for d in snakemake.params.dirs:
        row = {}
        # 1. Extract shared tree parameters from path
        tree_match = re.search(tree_regex, d)
        if tree_match:
            row.update(tree_match.groupdict())

        # 2. Extract tool-specific simulation parameters
        for tool_name, tool_conf in config["msa_sim_tools"].items():
            sim_regex = tool_conf["match_regex"]
            match = re.search(sim_regex, d)
            if match:
                row["msa_sim_tool"] = tool_name
                row.update(match.groupdict())
                break
        
        # 3. Extract inference params from path based on tool
        for inf_tool_name, inf_conf in config["inference_tools"].items():
            inf_regex = inf_conf["match_regex"]
            match = re.search(inf_regex, d)
            if match:
                row["inference_tool"] = inf_tool_name
                row.update(match.groupdict())
                break
        
        # 4. Load results from files
        row.update(get_distances_from_json(d))
        row["runtime_seconds"] = get_last_line_value(os.path.join(d, "time.txt"))
        row["log_likelihood"] = get_last_line_value(os.path.join(d, "logl.out"))
        
        # Determine MSA path - it is relative to the results directory structure
        # d is results/inference/{tree_params}/{msa_sim_tool}/{tool_params}/{inference_tool}/{inf_params}
        # MSA is results/msas/{tree_params}/{msa_sim_tool}/{tool_params}/msa.fasta
        # Tree PNG is results/trees/{tree_params}.nwk.png
        parts = d.split(os.sep)
        # Assuming path starts with 'results/inference/...'
        # We find 'inference' and swap it for 'msas', then remove the last part (jati_path)
        try:
            inf_idx = parts.index("inference")
            msa_parts = list(parts)
            msa_parts[inf_idx] = "msas"
            # The structure is .../inference/tree_params/msa_sim_tool/tool_params/inference_tool/inf_params
            # We want .../msas/tree_params/msa_sim_tool/tool_params/msa.fasta
            # So we need to remove two levels (inference_tool and inf_params)
            msa_path = os.sep.join(msa_parts[:-2] + ["msa.fasta"])
            
            tree_params = parts[inf_idx + 1]
            tree_png_path = os.path.join("results", "trees", f"{tree_params}.nwk.png")
            # Create a file:// URL that is clickable in many terminal emulators and spreadsheet apps
            abs_png_path = os.path.abspath(tree_png_path)
            row["tree_visualization"] = f"[tree](file://{abs_png_path})"
        except (ValueError, IndexError):
            msa_path = None
            row["tree_visualization"] = "NA"

        row["alignment_length"] = "NA"
        row["gap_percentage"] = "NA"

        if msa_path and os.path.exists(msa_path):
            row["alignment_length"] = get_fasta_length(msa_path)
            # Create a file:// URL for the MSA file
            abs_msa_path = os.path.abspath(msa_path)
            row["msa"] = f"[msa](file://{abs_msa_path})"
            
            if msa_path not in msa_cache:
                msa_cache[msa_path] = get_gap_stats(msa_path)
            
            gap_stats = msa_cache[msa_path]
            row["gap_percentage"] = gap_stats["gap_percentage"]
        
        # Add to collection
        all_rows.append(row)
        all_keys.update(row.keys())

    # 6. Write to TSV with a stable header
    # Define column order: Global -> Tree -> MSA -> Inference -> Metrics -> Visuals
    column_order = [
        # Global
        "seed", 
        # Tree Generation
        "species", "birth_rate", "death_rate", "sampling_fraction", "mutation_rate",
        # MSA Simulation
        "msa_sim_tool", "tkf_lambda", "tkf_mu", "tkf_r", "max_ins", "ir", "ip", "alignment_length", "gap_percentage", "msa",
        # Inference
        "inference_tool", "model", "gap_strategy", "log_likelihood", "runtime_seconds",
        # Results/Metrics
        "robinson_foulds", "kuhner_felsenstein", "start_robinson_foulds", "start_kuhner_felsenstein",
        # Visuals
        "tree_visualization"
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

def get_distances_from_json(dir_path):
    dist_path = os.path.join(dir_path, "distances.json")
    if os.path.exists(dist_path):
        try:
            with open(dist_path, 'r') as f:
                return json.load(f)
        except Exception:
            pass
    return {}

def get_fasta_length(msa_path):
    """Returns the length of the first sequence in the FASTA."""
    if os.path.exists(msa_path):
        try:
            with open(msa_path, 'r') as f:
                seq = ""
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    if line.startswith(">"):
                        if seq:
                            return len(seq)
                    else:
                        seq += line
                return len(seq)
        except Exception:
            pass
    return "NA"

def get_gap_stats(msa_path):
    """Calculates gap count and percentage by scanning the MSA file."""
    if not os.path.exists(msa_path):
        return {"gap_count": "NA", "gap_percentage": "NA"}
    try:
        total_chars = 0
        gaps = 0
        with open(msa_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith(">"):
                    continue
                total_chars += len(line)
                gaps += line.count('-')
        
        pct = (gaps / total_chars * 100) if total_chars > 0 else 0
        return {"gap_count": gaps, "gap_percentage": round(pct, 2)}
    except Exception:
        pass
    return {"gap_count": "NA", "gap_percentage": "NA"}

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

if __name__ == "__main__":
    main()
