import os
import yaml
import csv

project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

def add_to_ordered_set(ordered_set, new_keys):
    for key in new_keys:
        if key not in ordered_set:
            ordered_set.append(key)

def get_last_line_value(file_path):
    if not os.path.exists(file_path):
        return "NA"
    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()
            if lines:
                return lines[-1].strip()
    except (ValueError, IndexError):
        pass
    return "NA"

def load_snakemake_config_yaml(config_path = os.path.join(project_root, "config.yaml")):
    with open(config_path, "r") as f:
        return yaml.safe_load(f)
    
def write_table(rows, column_order, output_path):
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=column_order, delimiter='\t', extrasaction='ignore', restval='NA')
        writer.writeheader()
        writer.writerows(rows)
