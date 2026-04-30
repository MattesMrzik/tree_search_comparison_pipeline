import os
import pandas as pd
import dendropy
import glob

from snakemake_helpers import get_tool_params
from viz.utils import PROJECT_ROOT, load_snakemake_config_yaml

def tree_height(tree):
    """Sum of all branch lengths (tree height)."""
    return sum(edge.length for edge in tree.edges() if edge.length is not None)

def num_leaves(tree):
    """Number of taxa in the tree."""
    return len(tree.leaf_nodes())

def num_nodes(tree):
    """Total number of nodes."""
    return len(list(tree.nodes()))

def main():
    config = load_snakemake_config_yaml()
    base_dir = os.path.join(PROJECT_ROOT, "results/sim/tree")
    os.makedirs(base_dir, exist_ok=True)
    out_path = os.path.join(base_dir, "tree_summary.tsv")

    rows = []
    for root, _, _ in os.walk(base_dir):
        nwk_files = glob.glob(os.path.join(root, "seed*.nwk"))
        for nwk_file in nwk_files:
            row = {}
            rel_dir = os.path.relpath(root, PROJECT_ROOT)
            t = dendropy.Tree.get(path=nwk_file, schema="newick")
            tree_params = get_tool_params(nwk_file, config, "tree_sim")
            row.update(tree_params)

            print(f"Processing {rel_dir}...")
            row["tree_path"] = nwk_file
            row["tree_height"] = tree_height(t)
            row["num_leaves"] = num_leaves(t)
            row["num_nodes"] = num_nodes(t)
            rows.append(row)

    df = pd.DataFrame(rows)
    df.to_csv(out_path, sep="\t", index=False)

if __name__ == "__main__":
    main()
