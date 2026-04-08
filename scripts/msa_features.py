import numpy as np
# import pandas as pd
# from Bio import AlignIO
# from Bio import Phylo
# from io import StringIO
#




# the newick tree form alisim has internal node labels that we can use





# def get_clade_label(node):
#     """
#     Generate a label for a node by concatenating its leaf labels.
#     If it's a leaf, just returns the leaf name.
#     """
#     if node.is_terminal():
#         return node.name
#     # Sort to ensure consistent labeling
#     leaves = sorted([leaf.name for leaf in node.get_terminals()])
#     return "_".join(leaves)
#
# def collect_indels_from_alignment(msa_path, tree_path):
#     """
#     Collects insertion and deletion events by comparing parent and child sequences
#     along the tree. Assumes ancestral sequences are present in the MSA with names
#     matching the tree's node names (or labels).
#     """
#     alignment = AlignIO.read(msa_path, "fasta")
#     seq_dict = {record.id: str(record.seq) for record in alignment}
#     
#     tree = Phylo.read(tree_path, "newick")
#     
#     events = []
#     
#     # Traverse the tree from root
#     for parent in tree.find_clades(order='level'):
#         parent_label = parent.name # Assumes alignment IDs match node names
#         if parent_label not in seq_dict:
#             # If internal nodes don't have names in the tree, we might need 
#             # to handle how simulation tools (like JATI/simulate_tkf) label them.
#             # For now, assume they are there.
#             continue
#             
#         parent_seq = seq_dict[parent_label]
#         
#         for child in parent.clades:
#             child_label = child.name
#             if child_label not in seq_dict:
#                 continue
#                 
#             child_seq = seq_dict[child_label]
#             clade_label = get_clade_label(child)
#             
#             # Compare parent and child sequences
#             # 1. Ignore positions where both are gaps
#             # 2. Find stretches: parent gap, child char (Insertion)
#             # 3. Find stretches: parent char, child gap (Deletion)
#             
#             in_insertion = False
#             insertion_start = 0
#             in_deletion = False
#             deletion_start = 0
#             
#             seq_len = len(parent_seq)
#             for i in range(seq_len):
#                 p_char = parent_seq[i]
#                 c_char = child_seq[i]
#                 
#                 # Both gaps - ignore
#                 if p_char == '-' and c_char == '-':
#                     # If we were in an event, it ends here because we only look for
#                     # parent char vs child gap or vice versa.
#                     if in_insertion:
#                         events.append({'type': 'insertion', 'length': i - insertion_start, 'node': clade_label})
#                         in_insertion = False
#                     if in_deletion:
#                         events.append({'type': 'deletion', 'length': i - deletion_start, 'node': clade_label})
#                         in_deletion = False
#                     continue
#                 
#                 # Check for Insertion (Parent gap, Child char)
#                 if p_char == '-' and c_char != '-':
#                     if not in_insertion:
#                         in_insertion = True
#                         insertion_start = i
#                     if in_deletion:
#                         events.append({'type': 'deletion', 'length': i - deletion_start, 'node': clade_label})
#                         in_deletion = False
#                 
#                 # Check for Deletion (Parent char, Child gap)
#                 elif p_char != '-' and c_char == '-':
#                     if not in_deletion:
#                         in_deletion = True
#                         deletion_start = i
#                     if in_insertion:
#                         events.append({'type': 'insertion', 'length': i - insertion_start, 'node': clade_label})
#                         in_insertion = False
#                 
#                 # Both chars (match/mismatch)
#                 else:
#                     if in_insertion:
#                         events.append({'type': 'insertion', 'length': i - insertion_start, 'node': clade_label})
#                         in_insertion = False
#                     if in_deletion:
#                         events.append({'type': 'deletion', 'length': i - deletion_start, 'node': clade_label})
#                         in_deletion = False
#             
#             # Close pending events at the end of the sequence
#             if in_insertion:
#                 events.append({'type': 'insertion', 'length': seq_len - insertion_start, 'node': clade_label})
#             if in_deletion:
#                 events.append({'type': 'deletion', 'length': seq_len - deletion_start, 'node': clade_label})
#                 
#     return pd.DataFrame(events)

def gap_concentration(df):
    df['gap_concentration'] = np.where(
        df['gap_col%'] > 0,
        df['gap%'] / df['gap_col%'],
        0
    )
    return df
