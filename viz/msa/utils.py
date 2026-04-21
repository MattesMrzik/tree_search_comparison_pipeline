import os
from typing import Dict

RESULTS_MSA_DIR = "results/msas"

project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))

def all_msa_dirs(base_dir = os.path.join(project_root, RESULTS_MSA_DIR)):
    msa_dirs = []
    for root, _, files in os.walk(base_dir): # _ is dirs 
        if "msa.fasta" in files:
            msa_dirs.append(root)
    return msa_dirs

def load_msa(fasta_path: str) -> Dict[str, str]:
    msa = {}
    current_name = None
    current_seq = []

    with open(fasta_path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_name is not None:
                    msa[current_name] = "".join(current_seq)
                current_name = line[1:]
                current_seq = []
            else:
                current_seq.append(line)

        if current_name is not None:
            msa[current_name] = "".join(current_seq)

    assert len(msa) > 0, "MSA is empty"
    seq_len = len(next(iter(msa.values())))
    assert all(len(seq) == seq_len for seq in msa.values()), "All sequences must have the same length"

    return msa
