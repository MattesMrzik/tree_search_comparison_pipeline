import dendropy
import statistics
from typing import Dict, List, Set

from viz.msa.utils import load_msa
from viz.indel_points_inference.utils import (
    load_tree,
    infer_indels,
    EventType,
    IndelEvent,
    IndelEvents,
)

def _compute_indel_measures(
    tree: dendropy.Tree,
    true_events: IndelEvents,
    inferred_events: IndelEvents,
    suffix: str,
) -> Dict[str, float]:
    row = {}
    true_set = set(true_events.events)
    inferred_set = set(inferred_events.events)
    row.update(kimIndelignProbabilisticFramework2007(suffix, true_set, inferred_set, true_events, inferred_events))

    if suffix == "short":
        row.update(short_insertion_statistics(true_events, inferred_events))
        row.update(short_deletion_statistics(tree, true_events, inferred_events))

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

def _compute_diff_stats(step_diffs: List[float], len_diffs: List[float], prefix: str, row: dict) -> None:
    if len(step_diffs) == 0:
        return

    row[f"{prefix}_step_diff_mean"] = statistics.mean(step_diffs)
    row[f"{prefix}_step_diff_median"] = statistics.median(step_diffs)
    if len(step_diffs) > 1:
        row[f"{prefix}_step_diff_stdev"] = statistics.stdev(step_diffs)
    row[f"{prefix}_step_diff_min"] = min(step_diffs)
    row[f"{prefix}_step_diff_max"] = max(step_diffs)
    quartiles = statistics.quantiles(step_diffs, n=4)
    row[f"{prefix}_step_diff_q25"] = quartiles[0]
    row[f"{prefix}_step_diff_q75"] = quartiles[2]

    positive_step = [d for d in step_diffs if d > 0]
    negative_step = [d for d in step_diffs if d < 0]
    if positive_step:
        row[f"{prefix}_step_diff_positive_mean"] = statistics.mean(positive_step)
    if negative_step:
        row[f"{prefix}_step_diff_negative_mean"] = statistics.mean(negative_step)

    row[f"{prefix}_len_diff_mean"] = statistics.mean(len_diffs)
    row[f"{prefix}_len_diff_median"] = statistics.median(len_diffs)
    if len(len_diffs) > 1:
        row[f"{prefix}_len_diff_stdev"] = statistics.stdev(len_diffs)
    row[f"{prefix}_len_diff_min"] = min(len_diffs)
    row[f"{prefix}_len_diff_max"] = max(len_diffs)
    quartiles = statistics.quantiles(len_diffs, n=4)
    row[f"{prefix}_len_diff_q25"] = quartiles[0]
    row[f"{prefix}_len_diff_q75"] = quartiles[2]

    positive_len = [d for d in len_diffs if d > 0]
    negative_len = [d for d in len_diffs if d < 0]
    if positive_len:
        row[f"{prefix}_len_diff_positive_mean"] = statistics.mean(positive_len)
    if negative_len:
        row[f"{prefix}_len_diff_negative_mean"] = statistics.mean(negative_len)

    row[f"{prefix}_n"] = len(step_diffs)


# For insertions we look if they are higher or lower than the true events
# these distances are positive if the inferred event is further from the root, ie lower in the tree
def short_insertion_statistics(true_events: IndelEvents, inferred_events: IndelEvents) -> Dict[str, float]:
    suffix = "short"
    row = {}
    ins_step_diff = []
    ins_len_diff = []
    ins_step_diff_root_ins_at_true = [] # the len (ie n) of this the number of too many insertions
    ins_len_diff_root_ins_at_true = []
    ins_step_diff_root_ins_at_inf = [] # the len (ie n) of this the number of missing insertions
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
                        ins_step_diff.append(inf_event.distance_steps - true_event.distance_steps)
                        ins_len_diff.append(inf_event.distance_length - true_event.distance_length)
                        print("the diff is ", inf_event.distance_steps - true_event.distance_steps)
    # there can also be an insertion in the true that is not in the inferred.
    for true_event in true_events.events:
        if true_event.event_type == EventType.INSERTION:
            inferred_events_at_pos = inferred_events.get_by_column(true_event.start)
            if len(inferred_events_at_pos) == 0 or all(e.event_type != EventType.INSERTION for e in inferred_events_at_pos):
                # therefore we have a char at the root
                ins_step_diff_root_ins_at_inf.append(-true_event.distance_steps)
                ins_len_diff_root_ins_at_inf.append(-true_event.distance_length)

    _compute_diff_stats(ins_step_diff, ins_len_diff, f"{suffix}_ins", row)
    _compute_diff_stats(
        ins_step_diff_root_ins_at_true,
        ins_len_diff_root_ins_at_true,
        f"{suffix}_ins_root_ins_at_true",
        row,
    )
    _compute_diff_stats(
        ins_step_diff_root_ins_at_inf,
        ins_len_diff_root_ins_at_inf,
        f"{suffix}_ins_root_ins_at_inf",
        row,
    )
    return row


def _get_descendant_nodes(node: dendropy.Node) -> Set[str]:
    descendants = set()
    if node.label is not None:
        descendants.add(node.label)
    for child in node.child_node_iter():
        descendants.update(_get_descendant_nodes(child))
    return descendants

def short_deletion_statistics(
    tree: dendropy.Tree,
    true_events: IndelEvents,
    inferred_events: IndelEvents,
) -> Dict[str, float]:
    suffix = "short"
    row = {}

    del_correct = 0
    del_too_high_step = []
    del_too_high_len = []
    del_too_low_step = []
    del_too_low_len = []
    missing_dels = 0 # can happen if the insertion was inferred too low
    too_much_dels = 0 # can happen if the insertion was inferred too high

    columns = set()
    for event in true_events.events:
        if event.event_type == EventType.DELETION:
            columns.add(event.start)
    for event in inferred_events.events:
        if event.event_type == EventType.DELETION:
            columns.add(event.start)

    for col in columns:
        true_dels = set([e for e in true_events.get_by_column(col) if e.event_type == EventType.DELETION])
        inf_dels = set([e for e in inferred_events.get_by_column(col) if e.event_type == EventType.DELETION])

        for true_del in true_dels:
            descendant_events = inferred_events.get_events_below_node_for_column(tree, true_del.node, col)
            if len(descendant_events) == 0:
                missing_dels += 1
                continue
            assert all(e.event_type == EventType.DELETION for e in descendant_events)
            if len(descendant_events) == 1 and descendant_events[0].node == true_del.node:
                del_correct += 1
                inf_dels.remove(descendant_events[0])
                continue
            else:
                # TODO: do this without creating the list of differences 
                del_too_low_step.append(statistics.mean(e.distance_steps - true_del.distance_steps for e in descendant_events))
                del_too_low_len.append(statistics.mean(e.distance_length - true_del.distance_length for e in descendant_events))
                for e in descendant_events:
                    inf_dels.remove(e)

        for inf_del in inf_dels:
            descendant_events = true_events.get_events_below_node_for_column(tree, inf_del.node, col)
            if len(descendant_events) == 0:
                too_much_dels += 1
                continue
            assert all(e.event_type == EventType.DELETION for e in descendant_events)
            missing_dels -= len(descendant_events)
            del_too_high_step.append(statistics.mean(inf_del.distance_steps - e.distance_steps for e in descendant_events))
            del_too_high_len.append(statistics.mean(inf_del.distance_length - e.distance_length for e in descendant_events))


    _compute_diff_stats(del_too_high_step, del_too_high_len, f"{suffix}_del_too_high", row)
    _compute_diff_stats(del_too_low_step, del_too_low_len, f"{suffix}_del_too_low", row)
    row[f"{suffix}_del_correct"] = del_correct
    row[f"{suffix}_del_missing"] = missing_dels
    row[f"{suffix}_del_too_much"] = too_much_dels

    return row


def compare_indel_annotations(
    tree: dendropy.Tree,
    true_msa: Dict[str, str],
    inferred_msa: Dict[str, str],
) -> Dict[str, float]:
    true_events = infer_indels(true_msa, tree)
    inferred_events = infer_indels(inferred_msa, tree)
    long_measures = _compute_indel_measures(tree, true_events, inferred_events, "long")

    true_short = true_events.split_to_single_site()
    inferred_short = inferred_events.split_to_single_site()
    short_measures = _compute_indel_measures(tree, true_short, inferred_short, "short")

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
