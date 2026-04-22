import copy
from dataclasses import dataclass, replace
import dendropy
from enum import Enum
from typing import Dict, List, Set, Optional
from intervaltree import IntervalTree

class EventType(Enum):
    INSERTION = "insertion"
    DELETION = "deletion"

# This annotation is immutable to ensure that events can be safely used in sets and as dictionary keys.
@dataclass(frozen=True)
class IndelEvent:
    node: str
    start: int
    end: int
    event_type: EventType
    distance_steps: float = 0.0
    distance_length: float = 0.0

    def overlaps_column(self, col: int) -> bool:
        return self.start <= col < self.end

class IndelEvents:
    def __init__(self, events: Optional[List[IndelEvent]] = None):
        self.events: List[IndelEvent] = events or []
        self._by_node: Dict[str, List[IndelEvent]] = {}
        self._interval_tree_for_events = IntervalTree()
        self._build_indices()

    def _build_indices(self):
        self._by_node.clear()
        self._interval_tree_for_events = IntervalTree()
        for event in self.events:
            if event.node not in self._by_node:
                self._by_node[event.node] = []
            self._by_node[event.node].append(event)
            self._interval_tree_for_events.addi(event.start, event.end, event)

    def add(self, event: IndelEvent):
        self.events.append(event)
        if event.node not in self._by_node:
            self._by_node[event.node] = []
        self._by_node[event.node].append(event)
        self._interval_tree_for_events.addi(event.start, event.end, event)

    def get_by_node(self, node: str) -> List[IndelEvent]:
        return self._by_node.get(node, [])

    def get_by_column(self, col: int) -> List[IndelEvent]:
        intervals = self._interval_tree_for_events.overlap(col, col + 1)
        return [interval.data for interval in intervals]

    def get_columns_in_region(self, start: int, end: int) -> Set[int]:
        intervals = self._interval_tree_for_events.overlap(start, end)
        columns = set()
        for interval in intervals:
            ev = interval.data
            columns.update(range(max(ev.start, start), min(ev.end, end)))
        return columns

    def is_dollo(self) -> bool:
        covered_columns = set()
        for event in self.events:
            if event.event_type == EventType.INSERTION:
                for col in range(event.start, event.end):
                    if col in covered_columns:
                        return False
                    covered_columns.add(col)
        return True

    def count_by_type(self, event_type: EventType) -> int:
        return sum(1 for e in self.events if e.event_type == event_type)

    def split_to_single_site(self) -> "IndelEvents":
        new_events = IndelEvents()
        for event in self.events:
            for col in range(event.start, event.end):
                single_site_event = copy.deepcopy(event)
                single_site_event = replace(single_site_event, start=col, end=col + 1)
                new_events.add(single_site_event)
        return new_events

    # including the node itself
    def get_events_below_node_for_column(self, tree: dendropy.Tree, node_label: str, col: int) -> List[IndelEvent]:
        node = tree.find_node_with_label(node_label)
        if node is None:
            print(f"Warning: Node with label {node_label} not found in tree.")
            return []
        descendant_labels = {desc.label for desc in node.preorder_iter() if desc.label}
        events = []
        events_in_col = self.get_by_column(col)
        for event in events_in_col:
            if event.node in descendant_labels:
                events.append(event)
        return events

# TODO: if run-time becomes an issue, we can optimize this by caching distances to root for all nodes in the tree.
def _compute_distance_to_root(node: dendropy.Node) -> tuple:
    steps = 0
    length = 0.0
    current = node
    # we are checking is label since dendropy might create an extra dummy root node
    while current.parent_node.label is not None:
        steps += 1
        if current.edge_length is not None:
            length += current.edge_length
        current = current.parent_node
    return (steps, length)


def infer_indels(msa: Dict[str, str], tree: dendropy.Tree) -> IndelEvents:
    events = IndelEvents()
    node_to_seq = {node.label: msa[node.label] for node in tree.preorder_node_iter() if node.label}

    msa_len = len(next(iter(msa.values())))
    for node in tree.preorder_node_iter():
        if node.label is None:
            continue
        if node.parent_node is None:
            continue

        parent_label = node.parent_node.label
        if parent_label is None or parent_label not in node_to_seq:
            continue

        child_seq = node_to_seq[node.label]
        parent_seq = node_to_seq[parent_label]

        distance_steps, distance_length = _compute_distance_to_root(node)

        i = 0
        while i < msa_len:
            if parent_seq[i] == "-" and child_seq[i] != "-":
                start = i
                while i < msa_len and parent_seq[i] == "-" and child_seq[i] != "-":
                    i += 1
                end = i
                events.add(IndelEvent(
                    node=node.label,
                    start=start,
                    end=end,
                    event_type=EventType.INSERTION,
                    distance_steps=distance_steps,
                    distance_length=distance_length
                    ))
            elif parent_seq[i] != "-" and child_seq[i] == "-":
                start = i
                while i < msa_len and parent_seq[i] != "-" and child_seq[i] == "-":
                    i += 1
                end = i
                events.add(IndelEvent(
                    node=node.label,
                    start=start,
                    end=end,
                    event_type=EventType.DELETION,
                    distance_steps=distance_steps,
                    distance_length=distance_length
                    ))
            else:
                i += 1
    assert events.is_dollo(), "Inferred events violate Dollo's law, which should not happen with a correct implementation."
    return events

def load_tree(newick_path: str) -> dendropy.Tree:
    return dendropy.Tree.get(
        path=newick_path,
        schema="newick",
        preserve_underscores=True,
        suppress_internal_node_taxa=True,
        suppress_leaf_node_taxa=True,
    )
