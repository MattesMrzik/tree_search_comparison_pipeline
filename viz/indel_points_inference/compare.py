import dendropy
from typing import Dict, Set

from viz.msa.utils import load_msa
from viz.indel_points_inference.utils import (
    load_tree,
    infer_indels,
    EventType,
    IndelEvent,
    IndelEvents,
)

def _compute_indel_measures(
    true_events: IndelEvents,
    inferred_events: IndelEvents,
    suffix: str,
) -> Dict[str, float]:
    row = {}
    true_set = set(true_events.events)
    inferred_set = set(inferred_events.events)
    row.update(kimIndelignProbabilisticFramework2007(suffix, true_set, inferred_set, true_events, inferred_events))

    # for now only for the short events. But here how do we weight the mistakes, 
    # since long events that are here considered as many single site events will be 
    # more heavily penalized than short events.
    if suffix == "short": 
        row.update(short_insertion_statistics(true_events, inferred_events))

    return row

# TODO: when i have these the kimIndelignProbabilisticFramework2007 could also be calculated in the notebook and not here
def kimIndelignProbabilisticFramework2007(suffix,
                                          true_set: Set[IndelEvent],
                                          inferred_set: Set[IndelEvent],
                                          true_events: IndelEvents,
                                          inferred_events: IndelEvents
                                          ) -> Dict[str, float]:
    row = {}
    matched_true = true_set & inferred_set

    nit = true_events.count_by_type(EventType.INSERTION)
    ndt = true_events.count_by_type(EventType.DELETION)
    nie = inferred_events.count_by_type(EventType.INSERTION)
    nde = inferred_events.count_by_type(EventType.DELETION)
    row[f"{suffix}_nit"] = nit
    row[f"{suffix}_ndt"] = ndt
    row[f"{suffix}_nie"] = nie
    row[f"{suffix}_nde"] = nde

    # see kimIndelignProbabilisticFramework2007, best 1
    row[f"{suffix}_annotation_agreement"] = len(matched_true) / len(true_set) if len(true_set) > 0 else 0.0

    if ndt > 0 and nie > 0 and nde > 0:
        true_ratio = nit / ndt
        est_ratio = nie / nde
        # see kimIndelignProbabilisticFramework2007, best 1, can be larger than one, but also smaller i think
        row[f"{suffix}_indel_ratio"] = true_ratio / est_ratio

    nom = ((nit - nie) ** 2 + (ndt - nde) ** 2)
    denom = nit**2 + ndt**2
    if denom > 0:
        # see kimIndelignProbabilisticFramework2007, best 0
        row[f"{suffix}_indel_agreement"] = (nom / denom) ** 0.5
    return row

# For insertions we look if they are higher or lower than the true events
# these distances are positive if the inferred event is further from the root, ie lower in the tree
def short_insertion_statistics(true_events: IndelEvents, inferred_events: IndelEvents) -> Dict[str, float]:
    suffix = "short"
    row = {}
    ins_step_diff = []
    ins_len_diff = []
    ins_step_diff_root_ins_at_true = []
    ins_len_diff_root_ins_at_true = []
    ins_step_diff_root_ins_at_inf = []
    ins_len_diff_root_ins_at_inf = []

    # if the inferred is an insertion but the true has the char also at the root, 
    # then we dont have an insertion event in the true. 
    for inf_event in inferred_events.events:
        if inf_event.event_type == EventType.INSERTION:
            true_events_at_pos = true_events.get_by_column(inf_event.start)
            if len(true_events_at_pos) == 0 or all(e.event_type != EventType.INSERTION for e in true_events_at_pos):
                # therefore we have a char at the root 
                ins_step_diff_root_ins_at_true.append(inf_event.distance_steps)
                ins_len_diff_root_ins_at_true.append(inf_event.distance_length)
            else:
                # we have a true event at this position, so we look at the difference in steps and length
                for true_event in true_events_at_pos:
                    if true_event.event_type == EventType.INSERTION:
                        ins_step_diff.append(abs(inf_event.distance_steps - true_event.distance_steps))
                        ins_len_diff.append(abs(inf_event.distance_length - true_event.distance_length))
    if len(ins_step_diff) > 0:
        print(f"insertion step differences: {ins_step_diff}")
        row[f"{suffix}_ins_step_diff_mean"] = sum(ins_step_diff) / len(ins_step_diff)
        row[f"{suffix}_ins_len_diff_mean"] = sum(ins_len_diff) / len(ins_len_diff)
        row[f"{suffix}_ins_n"] = len(ins_step_diff)
    if len(ins_step_diff_root_ins_at_true) > 0:
        row[f"{suffix}_ins_step_diff_root_ins_at_true_mean"] = sum(ins_step_diff_root_ins_at_true) / len(ins_step_diff_root_ins_at_true)
        row[f"{suffix}_ins_len_diff_root_ins_at_true_mean"] = sum(ins_len_diff_root_ins_at_true) / len(ins_len_diff_root_ins_at_true)
        row[f"{suffix}_ins_root_ins_at_true_n"] = len(ins_step_diff_root_ins_at_true)

    # there can also be an insertion in the true that is not in the inferred.
    for true_event in true_events.events:
        if true_event.event_type == EventType.INSERTION:
            inferred_events_at_pos = inferred_events.get_by_column(true_event.start)
            if len(inferred_events_at_pos) == 0 or all(e.event_type != EventType.INSERTION for e in inferred_events_at_pos):
                # therefore we have a char at the root 
                ins_step_diff_root_ins_at_inf.append(-true_event.distance_steps)
                ins_len_diff_root_ins_at_inf.append(-true_event.distance_length)
    if len(ins_step_diff_root_ins_at_inf) > 0:
        row[f"{suffix}_ins_step_diff_root_ins_at_inf_mean"] = sum(ins_step_diff_root_ins_at_inf) / len(ins_step_diff_root_ins_at_inf)
        row[f"{suffix}_ins_len_diff_root_ins_at_inf_mean"] = sum(ins_len_diff_root_ins_at_inf) / len(ins_len_diff_root_ins_at_inf)
        row[f"{suffix}_ins_root_ins_at_inf_n"] = len(ins_step_diff_root_ins_at_inf)
    return row


def compare_indel_annotations(
    tree: dendropy.Tree,
    true_msa: Dict[str, str],
    inferred_msa: Dict[str, str],
) -> Dict[str, float]:
    true_events = infer_indels(true_msa, tree)
    inferred_events = infer_indels(inferred_msa, tree)
    long_measures = _compute_indel_measures(true_events, inferred_events, "long")

    true_short = true_events.split_to_single_site()
    inferred_short = inferred_events.split_to_single_site()
    print("true short events:")
    for e in true_short.events:
        print(e)
    print("inferred short events:")
    for e in inferred_short.events:
        print(e)
    short_measures = _compute_indel_measures(true_short, inferred_short, "short")

    return {**long_measures, **short_measures}


def compare_from_files(
    tree_path: str,
    true_msa_path: str,
    inferred_msa_path: str,
) -> Dict[str, float]:
    tree = load_tree(tree_path)
    true_msa = load_msa(true_msa_path)
    inferred_msa = load_msa(inferred_msa_path)
    return compare_indel_annotations(tree, true_msa, inferred_msa)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("tree", help="Path to tree Newick file")
    parser.add_argument("true_msa", help="Path to true MSA FASTA")
    parser.add_argument("inferred_msa", help="Path to inferred MSA FASTA")
    args = parser.parse_args()

    result = compare_from_files(args.tree, args.true_msa, args.inferred_msa)
    print(result)
