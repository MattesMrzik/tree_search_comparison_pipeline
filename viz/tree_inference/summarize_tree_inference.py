import os
import sys

project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

from viz.tree.utils import get_tree_params
from viz.msa.utils import get_msa_sim_params
from viz.tree_inference.utils import get_tree_inference_params
from viz.utils import load_snakemake_config_yaml, get_last_line_value, add_to_ordered_set, write_table
from utils import RESULTS_INF_DIR, all_inf_dirs, distances_for_true_vs_inferred

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
        t_params = get_tree_params(d, config["tree_match"])
        add_to_ordered_set(tree_col_names, t_params.keys())
        row.update(t_params)
        
        msa_params = get_msa_sim_params(d, config["msa_sim_tools"])
        add_to_ordered_set(msa_col_names, msa_params.keys())
        row.update(msa_params)
        
        i_params = get_tree_inference_params(d, config["inference_tools"])
        add_to_ordered_set(inf_col_names, i_params.keys())
        row.update(i_params)
        
        row.update(distances_for_true_vs_inferred(d))
        row["runtime_seconds"] = get_last_line_value(os.path.join(d, "time.txt"))
        row["logl"] = get_last_line_value(os.path.join(d, "logl.out"))
        # start tree distances in case we have the jati tool
        # true tree lolg, ie model search under jati

        all_rows.append(row)
        all_keys.update(row.keys())

    column_order = ["inf_dir"] + tree_col_names + msa_col_names + inf_col_names
    remaining_cols = sorted(list(all_keys - set(column_order)))
    column_order += remaining_cols

    write_table(all_rows, column_order, os.path.join(project_root, "results/summary.tsv"))

if __name__ == "__main__":
    main()
