import re
import os
import sys

RESULTS_INF_DIR = "results/inference"

project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

from viz.tree.calculate_distances import calculate_distances

def all_inf_dirs(base_dir = os.path.join(project_root, RESULTS_INF_DIR)):
    inf_dirs = []
    for root, _, files in os.walk(base_dir):
        if "final_tree.nwk" in files:
            inf_dirs.append(root)
    return inf_dirs

def get_tree_inference_params(path, inference_tools_config):
    for inf_tool_name, inf_conf in inference_tools_config.items():
        match = re.search(f'{inf_tool_name}/{inf_conf["match_regex"]}', path)
        if match:
            params = {"inference_tool": inf_tool_name}
            params.update(match.groupdict())
            return params
    return {}

def get_msa_dir_from_inf(inf_dir):
    parts = inf_dir.split(os.sep)
    inf_idx = parts.index("inference")
    msa_parts = list(parts)
    msa_parts[inf_idx] = "msas"
    msa_dir_path = os.sep.join(msa_parts[:-2])
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
