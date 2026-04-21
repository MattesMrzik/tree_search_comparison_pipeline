from dendropy import Tree, TaxonNamespace
from dendropy.calculate import treecompare

def calculate_distances(tree_path_a, tree_path_b):
    try: 
        tree_paths = [tree_path_a, tree_path_b]
        
        trees = read_trees(tree_paths)
        true_tree = trees[0]
        final_tree = trees[1]
        
        # Calculate distances for final tree
        final_rf = treecompare.symmetric_difference(true_tree, final_tree)
        final_kf = treecompare.euclidean_distance(true_tree, final_tree)

        results = {
            "rf": float(final_rf),
            "kf": float(final_kf),
        }

        return results
    except Exception as e:
        print(f"Error calculating distances: {e}")
        return {}
    

def read_trees(tree_paths):
    tns = TaxonNamespace()
    loaded_trees = []
    
    def reconcile_taxa(tree, namespace):
        for leaf in tree.leaf_node_iter():
            if leaf.label:
                taxon = namespace.require_taxon(label=leaf.label)
                leaf.taxon = taxon

    for path in tree_paths:
        try:
            tree = Tree.get(
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
            print(f"Error loading tree from {path}: {e}")

    return loaded_trees


