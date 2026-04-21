import os

project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))

from viz.utils import load_snakemake_config_yaml, get_last_line_value, add_to_ordered_set, write_table
from utils import RESULTS_INF_DIR, all_inf_dirs, distances_for_true_vs_inferred, distances_for_true_vs_start_nj_tree
from calculate_time import parse_jati_time, parse_iqtree_time
from snakemake_helpers import get_tool_params

def main():
    config = load_snakemake_config_yaml()
    
    results_inf_dir = os.path.join(project_root, RESULTS_INF_DIR)
    
    if not os.path.exists(results_inf_dir):
        print(f"Directory {results_inf_dir} does not exist.")
        return

    inf_dirs = all_inf_dirs(results_inf_dir)
    
    all_rows = []
    all_keys = set()
    tree_col_names = []
    msa_col_names = []
    inf_col_names = []
    
    for d in inf_dirs:
        row = {"inf_dir": d}
        print(f"Processing {d}...")
        
        # Extract parameters
        t_params = get_tool_params(d, config, "tree_sim")
        add_to_ordered_set(tree_col_names, t_params.keys())
        row.update(t_params)
        
        msa_params = get_tool_params(d, config, "msa_sim")
        add_to_ordered_set(msa_col_names, msa_params.keys())
        row.update(msa_params)
        
        i_params = get_tool_params(d, config, "tree_inf")
        add_to_ordered_set(inf_col_names, i_params.keys())
        row.update(i_params)
        
        row.update(distances_for_true_vs_inferred(d))
        log_path = os.path.join(d, "log.txt")
        row["runtime_seconds"] = parse_jati_time(log_path) if "jati" in d else parse_iqtree_time(log_path)
        row["logl"] = get_last_line_value(os.path.join(d, "logl.out"))
        if row.get("inference_tool") == "jati":
            distances = distances_for_true_vs_start_nj_tree(d)
            for dist_name, dist_value in distances.items():
                row[f"true_vs_start_nj_{dist_name}"] = dist_value
            
        # true tree lolg, ie model search under jati
        # also nj with just presence vs absence, ie no actual chars to compare it to the tkf indel tree inference

        all_rows.append(row)
        all_keys.update(row.keys())

    column_order = ["inf_dir"] + tree_col_names + msa_col_names + inf_col_names
    remaining_cols = sorted(list(all_keys - set(column_order)))
    column_order += remaining_cols

    write_table(all_rows, column_order, os.path.join(project_root, "results/tree_inf_summary.tsv"))

if __name__ == "__main__":
    main()
