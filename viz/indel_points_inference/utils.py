import dendropy
from dataclasses import dataclass
from enum import Enum
from typing import Dict, List, Set, Optional
from intervaltree import IntervalTree


class EventType(Enum):
    INSERTION = "insertion"
    DELETION = "deletion"


@dataclass
class IndelEvent:
    node: str
    start: int
    end: int
    event_type: EventType

    def overlaps_column(self, col: int) -> bool:
        return self.start <= col < self.end


class IndelEvents:
    def __init__(self, events: Optional[List[IndelEvent]] = None):
        self.events: List[IndelEvent] = events or []
        self._by_node: Dict[str, List[IndelEvent]] = {}
        self._tree = IntervalTree()
        self._build_indices()

    def _build_indices(self):
        self._by_node.clear()
        self._tree = IntervalTree()
        for event in self.events:
            if event.node not in self._by_node:
                self._by_node[event.node] = []
            self._by_node[event.node].append(event)
            self._tree.addi(event.start, event.end, event)

    def add(self, event: IndelEvent):
        self.events.append(event)
        if event.node not in self._by_node:
            self._by_node[event.node] = []
        self._by_node[event.node].append(event)
        self._tree.addi(event.start, event.end, event)

    def get_by_node(self, node: str) -> List[IndelEvent]:
        return self._by_node.get(node, [])

    def get_by_column(self, col: int) -> List[IndelEvent]:
        intervals = self._tree.overlap(col, col + 1)
        return [interval.data for interval in intervals]

    def get_columns_in_region(self, start: int, end: int) -> Set[int]:
        intervals = self._tree.overlap(start, end)
        columns = set()
        for interval in intervals:
            ev = interval.data
            columns.update(range(max(ev.start, start), min(ev.end, end)))
        return columns


def infer_indels(msa: Dict[str, str], tree: dendropy.Tree) -> IndelEvents:
    events = IndelEvents()
    node_to_seq = {node.label: msa[node.label] for node in tree.preorder_node_iter() if node.label}

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

        min_len = min(len(child_seq), len(parent_seq))
        i = 0
        while i < min_len:
            if parent_seq[i] == "-" and child_seq[i] != "-":
                start = i
                while i < min_len and parent_seq[i] == "-" and child_seq[i] != "-":
                    i += 1
                end = i
                events.add(IndelEvent(node=node.label, start=start, end=end, event_type=EventType.INSERTION))
            elif parent_seq[i] != "-" and child_seq[i] == "-":
                start = i
                while i < min_len and parent_seq[i] != "-" and child_seq[i] == "-":
                    i += 1
                end = i
                events.add(IndelEvent(node=node.label, start=start, end=end, event_type=EventType.DELETION))
            else:
                i += 1

        if len(parent_seq) > min_len and len(child_seq) == min_len:
            start = min_len
            end = len(parent_seq)
            if any(c != "-" for c in parent_seq[start:]):
                events.add(IndelEvent(node=node.label, start=start, end=end, event_type=EventType.DELETION))

        if len(child_seq) > min_len and len(parent_seq) == min_len:
            start = min_len
            end = len(child_seq)
            if any(c != "-" for c in child_seq[start:]):
                events.add(IndelEvent(node=node.label, start=start, end=end, event_type=EventType.INSERTION))

    return events


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

    return msa


def load_tree(newick_path: str) -> dendropy.Tree:
    return dendropy.Tree.get(
        path=newick_path,
        schema="newick",
        preserve_underscores=True,
        suppress_internal_node_taxa=True,
        suppress_leaf_node_taxa=True,
    )
