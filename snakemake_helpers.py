from itertools import product
from viz.utils import load_snakemake_config_yaml
from snakemake.io import expand

def infer_wildcard_constraints(cfg_dict):
    constraints = {"seed": r"\d+"}
    for tool_cfg in cfg_dict.values():
        for k, v in tool_cfg.items():
            if k == "path_snippet":
                continue
            if isinstance(v, list) and len(v) > 0 and isinstance(v[0], dict):
                # list of dicts -> flattened into separate wildcards, same as expand_tool_combos
                for inner_k, inner_v in v[0].items():
                    constraint = infer_constraint(inner_v)
                    if inner_k in constraints and constraints[inner_k] != constraint:
                        raise ValueError(f"Conflicting constraints for wildcard '{inner_k}': '{constraints[inner_k]}' vs '{constraint}'")
                    constraints[inner_k] = constraint
            elif isinstance(v, list) and len(v) > 0 and isinstance(v[0], list):
                # list of lists -> skip, should be restructured as list of dicts
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
        return r"[^_/]+"  # string: can't contain _ or /

def expand_tool_combos(cfg_dict):
    for tool, cfg in cfg_dict.items():
        param_lists  = {}  # independent params -> full product
        paired_lists = {}  # paired params -> zipped

        for k, v in cfg.items():
            if k == "path_snippet":
                continue
            if isinstance(v, list) and len(v) > 0 and isinstance(v[0], dict):
                # paired -> keep together as list of dicts
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

STAGES = {
    "tree": ("tree_sim_tool",  "tree_params", config["tree_sim_tools"]),
    "msa":  ("msa_sim_tool",   "msa_params",  config["msa_sim_tools"]),
    "tinf": ("inference_tool", "inf_params",  config["tree_inf_tools"]),
    "minf": ("inference_tool", "inf_params",  config["model_parameters_inf_tools"]),
}


def make_targets(path_template, *stage_keys):
    combos = [list(expand_tool_combos(STAGES[k][2])) for k in stage_keys]
    wc_pairs = [(STAGES[k][0], STAGES[k][1]) for k in stage_keys]
    targets = []
    for combo in product(*combos):
        kwargs = {wc: val
                  for (tool_wc, params_wc), (tool, params) in zip(wc_pairs, combo)
                  for wc, val in [(tool_wc, tool), (params_wc, params)]}
        targets += expand(path_template, seed=config["seeds"], **kwargs)
    return targets

def get_tree_path(tool_name):
    tool_conf = config["tree_sim_tools"][tool_name]
    return config["tree_path"].replace("{tree_sim_tool}", tool_name).replace("{tree_params}", tool_conf["path_snippet"])

# e.g. results/msas/{tree_sim_tool}/{tree_params}/tkf/l{lambda}_mu{mu}_r{r}_max{max_insertion}_len{root_length}_{model}/seed{seed}/msa.fasta
def get_msa_output(tool_name):
    tool_conf = config["msa_sim_tools"][tool_name]
    return config["msa_dir_path"].replace("{msa_sim_tool}", tool_name).replace("{msa_params}", tool_conf["path_snippet"])

def get_inf_output(tool_name):
    # Look up in both model parameters and tree inference tools
    if tool_name in config["model_parameters_inf_tools"]:
        tool_conf = config["model_parameters_inf_tools"][tool_name]
        return config["model_param_inf_dir"].replace("{inference_tool}", tool_name).replace("{inf_params}", tool_conf["path_snippet"])
    elif tool_name in config["tree_inf_tools"]:
        tool_conf = config["tree_inf_tools"][tool_name]
        return config["tree_inf_dir"].replace("{inference_tool}", tool_name).replace("{inf_params}", tool_conf["path_snippet"])
    else:
        raise ValueError(f"Unknown inference tool: {tool_name}")
