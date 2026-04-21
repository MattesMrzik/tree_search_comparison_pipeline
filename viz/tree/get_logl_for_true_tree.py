import os

from viz.utils import get_last_line_value

def get_true_tree_logl(inf_dir, row):
    if row.get("inference_tool") == "true_tree":
        return row.get("logl", "NA")
    model = row.get("model")
    gap = row.get("gap")
    if not (model and gap):
        return "NA"
    parts = inf_dir.split(os.sep)
    try:
        inf_idx = parts.index("inference")
        true_tree_parts = parts[:inf_idx+4]
        true_tree_parts.append("true_tree")
        true_tree_parts.append(f"{model}_{gap}")
        true_tree_dir = os.sep.join(true_tree_parts)
        return get_last_line_value(os.path.join(true_tree_dir, "logl.out"))
    except (ValueError, IndexError):
        return "NA"
