import os

RESULTS_DIR = "results"
TREE_INF_DIR = "tree_inference"
RESULTS_INF_DIR = os.path.join(RESULTS_DIR, TREE_INF_DIR)
MSA_DIR = "msas"

project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))

from viz.tree.calculate_distances import calculate_distances

def all_inf_dirs(base_dir = os.path.join(project_root, RESULTS_INF_DIR)):
    inf_dirs = []
    for root, _, files in os.walk(base_dir):
        if "final_tree.nwk" in files:
            inf_dirs.append(root)
    return inf_dirs

def get_msa_dir_from_inf(inf_dir):
    parts = inf_dir.split(os.sep)
    inf_idx = parts.index(TREE_INF_DIR)
    msa_parts = list(parts)
    msa_parts[inf_idx] = MSA_DIR
    # the last 3 parts are tree_inf_tool, tree_inf_params and seed, we want to keep seed
    msa_dir_path = os.sep.join(msa_parts[:-3] + [msa_parts[-1]])
    return msa_dir_path

def distances_for_true_vs_inferred(d):
        msa_dir = get_msa_dir_from_inf(d)
        true_tree_path = os.path.join(msa_dir, "tree.nwk")
        inferred_tree_path = os.path.join(d, "final_tree.nwk")
        return calculate_distances(true_tree_path, inferred_tree_path)

def distances_for_true_vs_start_nj_tree(d):
        msa_dir = get_msa_dir_from_inf(d)
        true_tree_path = os.path.join(msa_dir, "tree.nwk")
        inferred_tree_path = os.path.join(d, "start_tree.nwk")
        return calculate_distances(true_tree_path, inferred_tree_path)
