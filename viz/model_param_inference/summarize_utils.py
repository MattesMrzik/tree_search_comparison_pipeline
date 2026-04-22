import json
import os

RESULTS_DIR = "results"
MODEL_INF_DIR = "model_inference"
RESULTS_INF_DIR = os.path.join(RESULTS_DIR, MODEL_INF_DIR)

project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))

def all_model_inf_dirs(base_dir=os.path.join(project_root, RESULTS_INF_DIR)):
    inf_dirs = []
    for root, _, files in os.walk(base_dir):
        if "logl.out" in files:
            inf_dirs.append(root)
    return inf_dirs


def get_msa_dir_from_inf(inf_dir):
    parts = inf_dir.split(os.sep)
    inf_idx = parts.index(MODEL_INF_DIR)
    msa_parts = list(parts)
    msa_parts[inf_idx] = "msas"
    msa_dir_path = os.sep.join(msa_parts[:-3] + [msa_parts[-1]])
    return msa_dir_path


def load_params_json(path):
    if not os.path.exists(path):
        return {}
    with open(path, "r") as f:
        data = json.load(f)
    return data

