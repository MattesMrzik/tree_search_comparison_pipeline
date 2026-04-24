import os

from viz.tree_inference.calculate_time import parse_jati_time
from viz.utils import load_snakemake_config_yaml, add_to_ordered_set, write_table, get_last_line_value
from viz.indel_points_inference.summarize_utils import (
    RESULTS_INF_DIR,
    all_indel_inf_dirs,
    compare_indel_events,
    get_msa_dir_from_inf
)
from snakemake_helpers import get_tool_params

project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))

def main():
    config = load_snakemake_config_yaml()

    results_inf_dir = os.path.join(project_root, RESULTS_INF_DIR)

    if not os.path.exists(results_inf_dir):
        print(f"Directory {results_inf_dir} does not exist.")
        return

    inf_dirs = all_indel_inf_dirs(results_inf_dir)

    all_rows = []
    all_keys = set()
    tree_col_names = []
    msa_col_names = []
    inf_col_names = []

    for d in inf_dirs:
        row = {"inf_dir": d}
        print(f"Processing {d}...")

        tree_params = get_tool_params(d, config, "tree_sim")
        add_to_ordered_set(tree_col_names, tree_params.keys())
        row.update(tree_params)

        msa_params = get_tool_params(d, config, "msa_sim")
        add_to_ordered_set(msa_col_names, msa_params.keys())
        row.update(msa_params)

        inf_params = get_tool_params(d, config, "asr")
        add_to_ordered_set(inf_col_names, inf_params.keys())
        row.update(inf_params)

        # TODO also compare against the parsimony asr
        row.update(compare_indel_events(d))

        row["logl"] = get_last_line_value(os.path.join(d, "logl.out"))
        row["logl_true"] = get_last_line_value(os.path.join(get_msa_dir_from_inf(d), "sim_indel_logl.out"))

        log_path = os.path.join(d, "log.txt")
        row["time"] = parse_jati_time(log_path)

        all_rows.append(row)
        all_keys.update(row.keys())

    column_order = ["inf_dir"] + tree_col_names + msa_col_names + inf_col_names
    remaining_cols = sorted(list(all_keys - set(column_order)))
    column_order += remaining_cols

    write_table(
        all_rows, column_order, os.path.join(project_root, "results/indel_inf_summary.tsv")
    )

if __name__ == "__main__":
    main()
