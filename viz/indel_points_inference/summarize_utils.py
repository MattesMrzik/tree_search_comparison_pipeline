import os

RESULTS_DIR = "results"
INDEL_INF_DIR = "asr"
RESULTS_INF_DIR = os.path.join(RESULTS_DIR, INDEL_INF_DIR)
MSA_DIR = "msas"

project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))

from viz.indel_points_inference.compare import compare_from_files

def all_indel_inf_dirs(base_dir=os.path.join(project_root, RESULTS_INF_DIR)):
    inf_dirs = []
    for root, _, files in os.walk(base_dir):
        if "masa.fasta" in files:
            inf_dirs.append(root)
    return inf_dirs

# TODO: this is basically the same as the one in tree inference utils
def get_msa_dir_from_inf(inf_dir, inf_type = INDEL_INF_DIR):
    parts = inf_dir.split(os.sep)
    inf_idx = parts.index(inf_type)
    msa_parts = list(parts)
    msa_parts[inf_idx] = MSA_DIR
    # the last 3 parts are tree_inf_tool, tree_inf_params and seed, we want to keep seed
    msa_dir_path = os.sep.join(msa_parts[:-3] + [msa_parts[-1]])
    return msa_dir_path

def compare_indel_events(d, inf_type = INDEL_INF_DIR):
    msa_dir = get_msa_dir_from_inf(d, inf_type)
    tree_dir = get_msa_dir_from_inf(d, inf_type)
    tree_path = os.path.join(tree_dir, "tree.nwk")
    true_msa_path = os.path.join(msa_dir, "masa.fasta")
    inferred_msa_path = os.path.join(d, "masa.fasta")
    return compare_from_files(tree_path, true_msa_path, inferred_msa_path)
