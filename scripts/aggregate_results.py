import json
import os
import re
import csv

def main():
    config = snakemake.params.sn_config
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
        row = extract_params_from_inf_path(d, sim_regex, inf_regex)
        row.update(global_params)
        row.update(get_distances_from_json(d))
        row["runtime_seconds"] = get_last_line_value(os.path.join(d, "time.txt"))
        row["log_likelihood"] = get_last_line_value(os.path.join(d, "logl.out"))
        rows.append(row)
    write_rows_to_tsv(snakemake.output.tsv_path, rows)

def write_rows_to_tsv(output_path, rows):
    """
    Writes a list of dictionaries to a TSV file.
    Uses the keys from the first row as headers.
    """
    if not rows:
        return
    fieldnames = list(rows[0].keys())
    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)

def get_distances_from_json(dir_path):
    """
    Loads all metrics from distances.json if it exists.
    """
    dist_path = os.path.join(dir_path, "distances.json")
    if os.path.exists(dist_path):
        try:
            with open(dist_path, 'r') as f:
                return json.load(f)
        except Exception:
            pass
    return {}

def extract_params_from_inf_path(path, sim_regex, inf_regex):
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
    """Reads the last line of a file and returns it as a float if possible. Might be slow for large files"""
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
