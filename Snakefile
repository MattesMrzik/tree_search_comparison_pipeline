configfile: "config.yaml"

import os

# Path to common tools
EVOLVER = config["evolver_path"]
PYTHON = config["python_bin"]
JATI = config["jati_binary_path"]

# Tree Generation Parameters
SPECIES = config["species"]
BIRTH_DEATH_PAIRS = config["birth_death_rates"]
SEEDS = config["seeds"]
SAMPLING = config["sampling_fraction"]
MUTATION = config["mutation_rate"]

# JATI Inference Parameters
JATI_MODELS = config["jati_models"]
JATI_PARA_LIST = config["jati_params"]
GAP_STRATEGIES = config["gap_handling_strategies"]
MAX_ITERATIONS = config["max_iterations"]

# Tool Specific Configurations
MSA_SIM_TOOLS = config["msa_sim_tools"]

# Path Templates from config
# Resolve the tree_params and jati_path_snippet once so they're available to all paths
TREE_PARAMS_PATH_SNIPPET = config["tree_params_path_snippet"]
JATI_PATH_SNIPPET = config["jati_path_snippet"]

TREE_PATH = config["tree_path"].replace("{tree_params_path_snippet}", TREE_PARAMS_PATH_SNIPPET)
MSA_PATH = config["msa_dir_path"].replace("{tree_params_path_snippet}", TREE_PARAMS_PATH_SNIPPET)
INF_PATH = config["inf_dir_path"].replace("{tree_params_path_snippet}", TREE_PARAMS_PATH_SNIPPET).replace("{jati_path_snippet}", JATI_PATH_SNIPPET)

# Helper functions for directory expansion
def get_all_dirs(template):
    dirs = []
    for tool_name, tool_conf in MSA_SIM_TOOLS.items():
        # Get all parameter combinations for this tool
        if tool_name == "tkf":
            param_sets = [tool_conf["params_path_snipped"].format(**{
                "lambda": tool_conf["lambda"], "mu": tool_conf["mu"], 
                "r": tool_conf["r"], "max_ins": tool_conf["max_ins"]
            })]
        elif tool_name == "iqtree":
            param_sets = [tool_conf["params_path_snipped"].format(ir=p[0], ip=p[1]) 
                         for p in tool_conf["indel_params"]]
        
        # Expand template with tree wildcards and tool-specific params
        exp_dict = {
            "msa_sim_tool": [tool_name],
            "tool_params": param_sets,
            "s": SPECIES, "b": [p[0] for p in BIRTH_DEATH_PAIRS], 
            "d": [p[1] for p in BIRTH_DEATH_PAIRS], 
            "f": SAMPLING, "m": MUTATION, "seed": SEEDS
        }
        
        # Add JATI wildcards if present in template
        if "{model}" in template:
            exp_dict.update({"model": JATI_MODELS, "gap": GAP_STRATEGIES})
            
        dirs.extend(expand(template, **exp_dict))
    return dirs

rule all:
    input:
        "results/summary.tsv",
        [f"{d}/distances.json" for d in get_all_dirs(INF_PATH)],
        [f"{d}/tree_plot.png" for d in get_all_dirs(MSA_PATH)]

rule generate_tree:
    output:
        TREE_PATH
    params:
        evolver = EVOLVER
    shadow: "minimal"
    shell:
        """
        printf "2\\n{wildcards.s}\\n1 {wildcards.seed} 1\\n{wildcards.b} {wildcards.d} {wildcards.f} {wildcards.m}\\n0\\n" | {params.evolver} > /dev/null 2>&1
        tail -n 1 evolver.out > {output}
        """

# Helper functions for rule-specific path generation
def get_msa_output(tool_name):
    """Returns the tool-specific MSA directory template."""
    return MSA_PATH.replace("{msa_sim_tool}", tool_name).replace(
        "{tool_params}", MSA_SIM_TOOLS[tool_name]["params_path_snipped"]
    )

rule simulate_tkf_alignment:
    input:
        tree = TREE_PATH
    output:
        msa = get_msa_output("tkf") + "/msa.fasta",
        tree_copy = get_msa_output("tkf") + "/tree.nwk"
    params:
        bin = MSA_SIM_TOOLS["tkf"]["binary_path"],
        out_dir = lambda wildcards, output: os.path.dirname(output.msa)
    shell:
        """
        {params.bin} \
            --tree-file {input.tree} \
            --lambda {wildcards.lambda} \
            --mu {wildcards.mu} \
            --r {wildcards.r} \
            --max-insertion-length {wildcards.max_ins} \
            --seed {wildcards.seed} \
            --output-dir {params.out_dir}
        cp {input.tree} {output.tree_copy}
        """

rule simulate_iqtree_alignment:
    input:
        tree = TREE_PATH
    output:
        msa = get_msa_output("iqtree") + "/msa.fasta",
        tree_copy = get_msa_output("iqtree") + "/tree.nwk"
    params:
        bin = MSA_SIM_TOOLS["iqtree"]["binary_path"],
        model = MSA_SIM_TOOLS["iqtree"]["model"],
        out_dir = lambda wildcards, output: os.path.dirname(output.msa)
    shell:
        """
        mkdir -p {params.out_dir}
        {params.bin} \
            --alisim {params.out_dir}/msa \
            -m {params.model} \
            -t {input.tree} \
            --indel {wildcards.ir},{wildcards.ip} \
            --seed {wildcards.seed} \
            --out-format fasta
        mv {params.out_dir}/msa.fa {output.msa}
        cp {input.tree} {output.tree_copy}
        """

rule visualize_msa_tree:
    input:
        tree = MSA_PATH + "/tree.nwk"
    output:
        plot = MSA_PATH + "/tree_plot.png"
    params:
        script = "scripts/visualize_trees.py",
        py_bin = PYTHON
    shell:
        "{params.py_bin} {params.script} --tree-file {input.tree} --output-file {output.plot}"

rule jati_inference:
    input:
        msa = MSA_PATH + "/msa.fasta"
    output:
        start_tree = INF_PATH + "/jati_run_out/jati_run_start_tree.newick",
        final_tree = INF_PATH + "/jati_run_out/jati_run_tree.newick",
        logl = INF_PATH + "/jati_run_out/jati_run_logl.out",
        log = INF_PATH + "/jati_run_out/jati_run.log"
    params:
        bin = JATI,
        paras = " ".join(map(str, JATI_PARA_LIST)),
        log_level = "warn",
        max_iterations = MAX_ITERATIONS,
        out_base = INF_PATH
    shell:
        """
        mkdir -p {params.out_base}
        {params.bin} \
            --out-folder {params.out_base} \
            --seq-file {input.msa} \
            --model {wildcards.model} \
            --params {params.paras} \
            --gap-handling {wildcards.gap} \
            --seed {wildcards.seed} \
            -l {params.log_level} \
            --max-iterations {params.max_iterations} \
            --no-timestamp
        """

rule jati_cleanup:
    input:
        start_tree = INF_PATH + "/jati_run_out/jati_run_start_tree.newick",
        final_tree = INF_PATH + "/jati_run_out/jati_run_tree.newick",
        logl = INF_PATH + "/jati_run_out/jati_run_logl.out",
        log = INF_PATH + "/jati_run_out/jati_run.log"
    output:
        start_tree = INF_PATH + "/start_tree.newick",
        final_tree = INF_PATH + "/final_tree.newick",
        logl = INF_PATH + "/logl.out",
        log = INF_PATH + "/log.txt"
    shell:
        """
        mv {input.start_tree} {output.start_tree}
        mv {input.final_tree} {output.final_tree}
        mv {input.logl} {output.logl}
        mv {input.log} {output.log}
        """

rule calculate_distances:
    input:
        true_tree = MSA_PATH + "/tree.nwk",
        final_tree = INF_PATH + "/final_tree.newick"
    output:
        dist_file = INF_PATH + "/distances.json"
    params:
        script = "scripts/calculate_distances.py",
        py_bin = PYTHON
    shell:
        "{params.py_bin} {params.script} --tree1 {input.true_tree} --tree2 {input.final_tree} --output {output.dist_file}"

rule calculate_time:
    input:
        log = INF_PATH + "/log.txt"
    output:
        time_file = INF_PATH + "/time.txt"
    params:
        script = "scripts/calculate_time.py",
        py_bin = PYTHON
    shell:
        "{params.py_bin} {params.script} {input.log} {output.time_file}"

rule aggregate_summary:
    input:
        [f"{d}/distances.json" for d in get_all_dirs(INF_PATH)],
        [f"{d}/time.txt" for d in get_all_dirs(INF_PATH)],
        [f"{d}/logl.out" for d in get_all_dirs(INF_PATH)]
    output:
        tsv_path = "results/summary.tsv"
    params:
        dirs = get_all_dirs(INF_PATH),
        sn_config = config
    script:
        "scripts/aggregate_results.py"
