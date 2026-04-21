import os

from viz.msa.utils import load_msa
from viz.indel_points_inference.utils import load_tree, infer_indels, EventType

def test_infer_indels():
    base_dir = os.path.join(os.path.dirname(__file__), "data")
    msa_path = os.path.join(base_dir, "msa.fasta")
    tree_path = os.path.join(base_dir, "tree.nwk")

    msa = load_msa(msa_path)
    tree = load_tree(tree_path)

    events = infer_indels(msa, tree)

    print(f"Total events: {len(events.events)}")
    for e in events.events:
        print(e)

    assert len(events.events) == 4, f"Expected 4 events, got {len(events.events)}"

    i3_events = events.get_by_node("I3")
    assert len(i3_events) == 1 
    assert i3_events[0].event_type == EventType.INSERTION
    assert i3_events[0].start == 1
    assert i3_events[0].end == 3 
    assert i3_events[0].distance_steps == 1
    assert i3_events[0].distance_length == 5

    b2_events = events.get_by_node("B2")
    assert len(b2_events) == 2
    assert b2_events[0].event_type == EventType.DELETION
    assert b2_events[0].start == 2
    assert b2_events[0].end == 3
    assert b2_events[1].event_type == EventType.INSERTION
    assert b2_events[1].start == 7 
    assert b2_events[1].end == 10
    assert b2_events[1].distance_steps == 2
    assert b2_events[1].distance_length == 3 + 5

    a1_events = events.get_by_node("A1")
    assert len(a1_events) == 1
    assert a1_events[0].event_type == EventType.INSERTION
    assert a1_events[0].start ==3 
    assert a1_events[0].end == 7 
    assert a1_events[0].distance_steps == 2
    assert a1_events[0].distance_length == 2 + 5

    col2_events = events.get_by_column(2)
    assert len(col2_events) == 2 
    nodes = {e.node for e in col2_events}
    assert "I3" in nodes
    assert "B2" in nodes

    print("All tests passed!")

if __name__ == "__main__":
    test_infer_indels()
