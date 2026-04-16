from itertools import product
from viz.utils import load_snakemake_config_yaml
from snakemake.io import expand

def infer_wildcard_constraints(cfg_dict):
    constraints = {"seed": r"\d+"}
    for tool_name, tool_cfg in cfg_dict.items():
        for k, v in tool_cfg.items():
            if k == "path_snippet":
                continue
            if isinstance(v, list) and len(v) > 0 and isinstance(v[0], dict):
                for inner_k, inner_v in v[0].items():
                    constraint = infer_constraint(inner_v)
                    if inner_k in constraints and constraints[inner_k] != constraint:
                        raise ValueError(f"Conflicting constraints for wildcard '{inner_k}': '{constraints[inner_k]}' vs '{constraint}'")
                    constraints[inner_k] = constraint
            elif isinstance(v, list) and len(v) > 0 and isinstance(v[0], list):
                raise ValueError(f"Use list of dicts instead of list of lists for {k}")
            else:
                val = v[0] if isinstance(v, list) else v
                constraint = infer_constraint(val)
                if k in constraints and constraints[k] != constraint:
                    raise ValueError(f"Conflicting constraints for wildcard '{k}': '{constraints[k]}' vs '{constraint}'")
                constraints[k] = constraint
    return constraints

def infer_constraint(val):
    if isinstance(val, bool):
        return r"[^_/]+"
    elif isinstance(val, int):
        return r"\d+"
    elif isinstance(val, float):
        return r"[\d.]+"
    else:
        return r"[^_/]+"

def expand_tool_combos(cfg_dict):
    for tool, cfg in cfg_dict.items():
        param_lists  = {}
        paired_lists = {}

        for k, v in cfg.items():
            if k == "path_snippet":
                continue
            if isinstance(v, list) and len(v) > 0 and isinstance(v[0], dict):
                paired_lists[k] = v
            else:
                param_lists[k] = v if isinstance(v, list) else [v]

        for ind_combo in product(*param_lists.values()):
            ind_dict = dict(zip(param_lists.keys(), ind_combo))
            if paired_lists:
                for paired_combo in product(*[v for v in paired_lists.values()]):
                    paired_dict = {}
                    for d in paired_combo:
                        paired_dict.update(d)
                    combo_dict = {**ind_dict, **paired_dict}
                    snippet = cfg["path_snippet"].format(**combo_dict)
                    yield tool, snippet
            else:
                snippet = cfg["path_snippet"].format(**ind_dict)
                yield tool, snippet

config = load_snakemake_config_yaml()

# To match the wildcards in the paths to the correct stage
STAGES = {
    "tree_sim": {"tool_wc": "tree_sim_tool", "params_wc": "tree_params"},
    "msa_sim": {"tool_wc": "msa_sim_tool", "params_wc": "msa_params"},
    "model_param_inf": {"tool_wc": "inference_tool", "params_wc": "inf_params"},
    "tree_inf": {"tool_wc": "inference_tool", "params_wc": "inf_params"},
}

def make_targets(cfg, *stages, primary, suffix=""):
    all_stages = list(stages) + [primary]
    path_template = cfg[primary]["dir"] + suffix

    combos = [list(expand_tool_combos(cfg[stage]["tools"])) for stage in all_stages]
    wc_pairs = [(STAGES[stage]["tool_wc"], STAGES[stage]["params_wc"]) for stage in all_stages]

    targets = []
    for combo in product(*combos):
        kwargs = {wc: val
                  for (tool_wc, params_wc), (tool, params) in zip(wc_pairs, combo)
                  for wc, val in [(tool_wc, tool), (params_wc, params)]}
        targets += expand(path_template, seed=cfg["seeds"], **kwargs)
    return targets

def get_tree_path(tool_name):
    tool_conf = config["tree_sim"]["tools"][tool_name]
    snippet = tool_conf.get("path_snippet", tool_name)
    return config["tree_sim"]["dir"].replace("{tree_sim_tool}", tool_name).replace("{tree_params}", snippet)

def get_msa_output(tool_name):
    tool_conf = config["msa_sim"]["tools"][tool_name]
    snippet = tool_conf.get("path_snippet", tool_name)
    return config["msa_sim"]["dir"].replace("{msa_sim_tool}", tool_name).replace("{msa_params}", snippet)

def get_inf_output(tool_type, tool_name):
    tool_conf = config[tool_type]["tools"][tool_name]
    snippet = tool_conf.get("path_snippet", tool_name)
    return config[tool_type]["dir"].replace("{inference_tool}", tool_name).replace("{inf_params}", snippet)
