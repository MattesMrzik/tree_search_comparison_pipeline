import dendropy
import sys
import json
import argparse

def calculate_distances(true_tree_path, final_tree_path, start_tree_path, output_path):
    tree_paths = [true_tree_path, final_tree_path]
    if start_tree_path:
        tree_paths.append(start_tree_path)
    
    trees = read_trees(tree_paths)
    true_tree = trees[0]
    final_tree = trees[1]
    
    # Calculate distances for final tree
    final_rf = dendropy.calculate.treecompare.symmetric_difference(true_tree, final_tree)
    final_kf = dendropy.calculate.treecompare.euclidean_distance(true_tree, final_tree)

    results = {
        "rf": float(final_rf),
        "kf": float(final_kf),
        "start_rf": "NA",
        "start_kf": "NA"
    }

    # Calculate distances for start tree if provided
    if start_tree_path:
        start_tree_obj = trees[2]
        start_rf = dendropy.calculate.treecompare.symmetric_difference(true_tree, start_tree_obj)
        start_kf = dendropy.calculate.treecompare.euclidean_distance(true_tree, start_tree_obj)
        
        results["start_rf"] = float(start_rf)
        results["start_kf"] = float(start_kf)
    
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=4)

def read_trees(tree_paths):
    tns = dendropy.TaxonNamespace()
    loaded_trees = []
    
    def reconcile_taxa(tree, namespace):
        for leaf in tree.leaf_node_iter():
            if leaf.label:
                taxon = namespace.require_taxon(label=leaf.label)
                leaf.taxon = taxon

    for path in tree_paths:
        try:
            tree = dendropy.Tree.get(
                path=path,
                schema="newick",
                preserve_underscores=True,
                suppress_internal_node_taxa=True,
                suppress_leaf_node_taxa=True,
                taxon_namespace=tns
            )
            reconcile_taxa(tree, tns)
            tree.is_rooted = False
            loaded_trees.append(tree)
        except Exception as e:
            print(f"Error loading tree from {path}: {e}", file=sys.stderr)
            sys.exit(1)

    return loaded_trees

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate RF and KF distances between two trees.")
    parser.add_argument("--true-tree", required=True, help="Path to the true tree (Newick)")
    parser.add_argument("--final-tree", required=True, help="Path to the final inferred tree (Newick)")
    parser.add_argument("--start-tree", help="Path to the starting tree (NJ) (Newick)")
    parser.add_argument("--output", required=True, help="Path to the output JSON file")
    
    args = parser.parse_args()
    calculate_distances(args.true_tree, args.final_tree, args.start_tree, args.output)
