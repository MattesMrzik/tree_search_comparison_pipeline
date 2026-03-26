import json
import os
import re
import csv

def extract_params_from_path(path, sim_regex, inf_regex):
    """
    Parses parameters from paths using regex with named groups from config.
    """
    params = {}
    
    sim_match = re.search(sim_regex, path)
    if sim_match:
        params.update(sim_match.groupdict())
    
    inf_match = re.search(inf_regex, path)
    if inf_match:
        params.update(inf_match.groupdict())
        
    return params

def get_last_line_value(file_path):
    """Reads the last line of a file and returns it as a float if possible."""
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
    # Use Snakemake's automatic object
    # snakemake.params.dirs, snakemake.params.full_config, snakemake.output[0]
    
    config = snakemake.params.full_config
    sim_regex = config["sim_match"]
    inf_regex = config["inf_match"]

    global_params = {
        "tkf_lambda": config.get("tkf_lambda", "NA"),
        "tkf_mu": config.get("tkf_mu", "NA"),
        "tkf_r": config.get("tkf_r", "NA"),
        "max_insertion_length": config.get("max_insertion_length", "NA"),
        "max_iterations": config.get("max_iterations", "NA")
    }

    rows = []
    for d in snakemake.params.dirs:
        # 1. Start with wildcards extracted from the path using named groups
        row = extract_params_from_path(d, sim_regex, inf_regex)
        
        # 2. Add global constants from config
        row.update(global_params)
        
        # 3. Add metrics from files
        # Distances
        dist_path = os.path.join(d, "distances.json")
        if os.path.exists(dist_path):
            try:
                with open(dist_path, 'r') as f:
                    dists = json.load(f)
                    # Use a consistent order or just update
                    for key in ["robinson_foulds", "kuhner_felsenstein"]:
                        row[key] = dists.get(key, "NA")
            except Exception:
                row.update({"robinson_foulds": "NA", "kuhner_felsenstein": "NA"})
        else:
            row.update({"robinson_foulds": "NA", "kuhner_felsenstein": "NA"})

        # Time
        time_path = os.path.join(d, "time.txt")
        row["runtime_seconds"] = get_last_line_value(time_path)

        # Log-Likelihood (from logl.out)
        logl_path = os.path.join(d, "logl.out")
        row["log_likelihood"] = get_last_line_value(logl_path)
        
        rows.append(row)

    if not rows:
        return

    # Write to TSV
    # We use the keys from the first row as headers to ensure all columns are included
    fieldnames = list(rows[0].keys())
    
    with open(snakemake.output[0], 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)

if __name__ == "__main__":
    main()
