import os

from viz.indel_points_inference.compare import compare_from_files

def test_compare_insertion_annotations():
    base_dir = os.path.join(os.path.dirname(__file__), "data")
    tree_path = os.path.join(base_dir, "tree.nwk")
    true_msa_path = os.path.join(base_dir, "true_msa.fasta")
    inferred_msa_path = os.path.join(base_dir, "inferred_msa.fasta")

    result = compare_from_files(tree_path, true_msa_path, inferred_msa_path)
    # for k, v in result.items():
    #     print(f"{k}: {v}")

    assert result["long_nit"] == 3
    assert result["long_ndt"] == 2 
    assert result["long_nie"] == 3 
    assert result["long_nde"] == 2 
    assert result["long_annotation_agreement"] == 0.0
    assert result["long_indel_ratio"] == 1 
    assert result["long_indel_agreement"] == 0

    assert result["short_nit"] == 9
    assert result["short_ndt"] == 2
    assert result["short_nie"] == 7 
    assert result["short_nde"] == 6 
    assert result["short_annotation_agreement"] == 5/11
    assert result["short_indel_ratio"] == (9/2) / (7/6)
    assert result["short_indel_agreement"] == ((2**2+4**2)/(9**2+2**2)) ** 0.5

    assert result["short_ins_step_diff_mean"] == 1/6
    assert result["short_ins_len_diff_mean"] == 2/6
    assert result["short_ins_n"] == 6
    assert result["short_ins_root_ins_at_true_step_diff_mean"] == 1
    assert result["short_ins_root_ins_at_true_len_diff_mean"] == 7
    assert result["short_ins_root_ins_at_true_n"] == 1
    assert result["short_ins_root_ins_at_inf_step_diff_mean"] == -2
    assert result["short_ins_root_ins_at_inf_len_diff_mean"] == -3 -5 
    assert result["short_ins_root_ins_at_inf_n"] == 3

    assert result["short_nie"] == result["short_ins_root_ins_at_true_n"] + result["short_ins_n"]

    print("All test_compare_insertion_annotations passed!")

def test_compare_deletion_annotations():
    base_dir = os.path.join(os.path.dirname(__file__), "data/deletions")
    tree_path = os.path.join(base_dir, "tree.nwk")
    true_msa_path = os.path.join(base_dir, "true_msa.fasta")
    inferred_msa_path = os.path.join(base_dir, "inferred_msa.fasta")

    result = compare_from_files(tree_path, true_msa_path, inferred_msa_path)
    # for k, v in result.items():
    #     print(f"{k}: {v}")

    assert result["short_del_too_high_step_diff_mean"] == -1
    assert result["short_del_too_high_len_diff_mean"] == -2.5
    assert result["short_del_too_high_n"] == 1
    assert result["short_del_too_low_step_diff_mean"] == 1
    assert result["short_del_too_low_len_diff_mean"] == 2.5
    assert result["short_del_too_low_n"] == 1
    assert result["short_del_correct"] == 3
    assert result["short_del_missing"] == 2
    assert result["short_del_too_much"] == 1

    print("All test_compare_deletion_annotations passed!")

if __name__ == "__main__":
    test_compare_insertion_annotations()
    test_compare_deletion_annotations()
