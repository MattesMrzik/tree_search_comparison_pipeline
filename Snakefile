configfile: "config.yaml"

import os
from scripts.utils.snakemake_helpers import get_msa_length

# Path to common tools
EVOLVER = config["evolver_path"]
PYTHON = config["python_bin"]

# Tree Generation Parameters
SPECIES = config["species"]
BIRTH_DEATH_PAIRS = config["birth_death_rates"]
SEEDS = config["seeds"]
SAMPLING = config["sampling_fraction"]
MUTATION = config["mutation_rate"]

# Inference Tools & Parameters
INF_TOOLS = config["inference_tools"]
MSA_SIM_TOOLS = config["msa_sim_tools"]

# Path Templates from config
# Resolve the tree_params once so they're available to all paths
TREE_PARAMS_PATH_SNIPPET = config["tree_params_path_snippet"]

TREE_PATH = config["tree_path"].replace("{tree_params_path_snippet}", TREE_PARAMS_PATH_SNIPPET)
MSA_PATH = config["msa_dir_path"].replace("{tree_params_path_snippet}", TREE_PARAMS_PATH_SNIPPET)
INF_PATH = config["inf_dir_path"].replace("{tree_params_path_snippet}", TREE_PARAMS_PATH_SNIPPET)


# TODO if msas get simulated then if they are finished, then jati is run already, even though msa simulation is not finished, and not all msas lens are considered during the priority calculation
# we could first call the simulation rule and then if that is finished call the inference rule

# Helper functions for priority and path generation
def get_all_msa_dirs():
    dirs = []
    for tool_name, tool_conf in MSA_SIM_TOOLS.items():
        if tool_name == "tkf":
            param_sets = [tool_conf["params_path_snipped"].format(**{
                "lambda": tool_conf["lambda"], "mu": tool_conf["mu"], 
                "r": tool_conf["r"], "max_ins": tool_conf["max_ins"],
                "root_length": rl
            }) for rl in tool_conf["root_lengths"]]
        elif tool_name == "alisim":
            param_sets = [tool_conf["params_path_snipped"].format(ir=p[0], ip=p[1], root_length=rl) 
                         for p in tool_conf["indel_params"] for rl in tool_conf["root_lengths"]]
        
        exp_dict = {
            "msa_sim_tool": [tool_name],
            "tool_params": param_sets,
            "s": SPECIES, "b": [p[0] for p in BIRTH_DEATH_PAIRS], 
            "d": [p[1] for p in BIRTH_DEATH_PAIRS], 
            "f": SAMPLING, "m": MUTATION, "seed": SEEDS
        }
        dirs.extend(expand(MSA_PATH, **exp_dict))
    return dirs

def get_all_dirs(template):
    dirs = []
    for tool_name, tool_conf in MSA_SIM_TOOLS.items():
        # Get all parameter combinations for this tool
        if tool_name == "tkf":
            param_sets = [tool_conf["params_path_snipped"].format(**{
                "lambda": tool_conf["lambda"], "mu": tool_conf["mu"], 
                "r": tool_conf["r"], "max_ins": tool_conf["max_ins"],
                "root_length": rl
            }) for rl in tool_conf["root_lengths"]]
        elif tool_name == "alisim":
            param_sets = [tool_conf["params_path_snipped"].format(ir=p[0], ip=p[1], root_length=rl) 
                         for p in tool_conf["indel_params"] for rl in tool_conf["root_lengths"]]
        
        for inf_tool_name, inf_conf in INF_TOOLS.items():
            # Build inference parameters
            if inf_tool_name == "jati":
                inf_params = expand(inf_conf["path_snippet"], 
                                    model=inf_conf["models"], 
                                    gap=inf_conf["gap_strategies"])
            elif inf_tool_name == "iqtree":
                inf_params = expand(inf_conf["path_snippet"], 
                                    model=inf_conf["models"])

            # Expand template with tree wildcards, tool-specific params, and inference params
            exp_dict = {
                "msa_sim_tool": [tool_name],
                "tool_params": param_sets,
                "inference_tool": [inf_tool_name],
                "inf_params": inf_params,
                "s": SPECIES, "b": [p[0] for p in BIRTH_DEATH_PAIRS], 
                "d": [p[1] for p in BIRTH_DEATH_PAIRS], 
                "f": SAMPLING, "m": MUTATION, "seed": SEEDS
            }
            
            dirs.extend(expand(template, **exp_dict))
    return dirs

rule all:
    input:
        "results/summary.tsv",
        "results/msa_summary.tsv",
        "/Users/mrzi/Seafile/phd_obsidian/notes/pipeline_out/summary_table.md",
        [f"{d}/distances.json" for d in get_all_dirs(INF_PATH)],

rule simulate_alignments:
    input:
        [f"{d}/msa.fasta" for d in get_all_dirs(MSA_PATH)]

rule generate_trees:
    input:
        expand(TREE_PATH, 
               s=SPECIES, 
               b=[p[0] for p in BIRTH_DEATH_PAIRS], 
               d=[p[1] for p in BIRTH_DEATH_PAIRS], 
               f=SAMPLING, m=MUTATION, seed=SEEDS)

rule visualize_trees:
    input:
        expand(f"{TREE_PATH}.png", 
               s=SPECIES, 
               b=[p[0] for p in BIRTH_DEATH_PAIRS], 
               d=[p[1] for p in BIRTH_DEATH_PAIRS], 
               f=SAMPLING, m=MUTATION, seed=SEEDS)

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

rule visualize_tree:
    input:
        tree = TREE_PATH
    output:
        plot = f"{TREE_PATH}.png"
    params:
        script = "scripts/visualize_trees.py",
        py_bin = PYTHON
    shell:
        "{params.py_bin} {params.script} --tree-file {input.tree} --output-file {output.plot}"

# Helper functions for rule-specific path generation
def get_msa_output(tool_name):
    """Returns the tool-specific MSA directory template."""
    tool_conf = MSA_SIM_TOOLS[tool_name]
    return MSA_PATH.replace("{msa_sim_tool}", tool_name).replace(
        "{tool_params}", tool_conf["params_path_snipped"]
    )

def get_inf_output(tool_name):
    """Returns the tool-specific inference directory template."""
    tool_conf = INF_TOOLS[tool_name]
    return INF_PATH.replace("{inference_tool}", tool_name).replace(
        "{inf_params}", tool_conf["path_snippet"]
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
            --root-length {wildcards.root_length} \
            --seed {wildcards.seed} \
            --output-dir {params.out_dir}
        cp {input.tree} {output.tree_copy}
        """

rule simulate_alisim_alignment:
    input:
        tree = TREE_PATH
    output:
        msa = get_msa_output("alisim") + "/msa.fasta",
        tree_copy = get_msa_output("alisim") + "/tree.nwk"
    params:
        bin = MSA_SIM_TOOLS["alisim"]["binary_path"],
        model = MSA_SIM_TOOLS["alisim"]["model"],
        out_dir = lambda wildcards, output: os.path.dirname(output.msa)
    shell:
        """
        mkdir -p {params.out_dir}
        {params.bin} \
            --alisim {params.out_dir}/msa \
            -m {params.model} \
            -t {input.tree} \
            --indel {wildcards.ir},{wildcards.ip} \
            --length {wildcards.root_length} \
            --seed {wildcards.seed} \
            --out-format fasta \
            --no-unaligned 
        mv {params.out_dir}/msa.fa {output.msa}
        cp {input.tree} {output.tree_copy}
        """

rule jati_inference:
    input:
        msa = MSA_PATH + "/msa.fasta"
    output:
        start_tree = get_inf_output("jati") + "/jati_run_out/jati_run_start_tree.newick",
        final_tree = get_inf_output("jati") + "/jati_run_out/jati_run_tree.newick",
        logl = get_inf_output("jati") + "/jati_run_out/jati_run_logl.out",
        log = get_inf_output("jati") + "/jati_run_out/jati_run.log"
    params:
        bin = INF_TOOLS["jati"]["binary_path"],
        paras = " ".join(map(str, INF_TOOLS["jati"]["params"])),
        log_level = "warn",
        max_iterations = INF_TOOLS["jati"]["max_iterations"],
        out_base = get_inf_output("jati")
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
        start_tree = get_inf_output("jati") + "/jati_run_out/jati_run_start_tree.newick",
        final_tree = get_inf_output("jati") + "/jati_run_out/jati_run_tree.newick",
        logl = get_inf_output("jati") + "/jati_run_out/jati_run_logl.out",
        log = get_inf_output("jati") + "/jati_run_out/jati_run.log"
    output:
        start_tree = get_inf_output("jati") + "/start_tree.newick",
        final_tree = get_inf_output("jati") + "/final_tree.newick",
        logl = get_inf_output("jati") + "/logl.out",
        log = get_inf_output("jati") + "/log.txt"
    shell:
        """
        mv {input.start_tree} {output.start_tree}
        mv {input.final_tree} {output.final_tree}
        mv {input.logl} {output.logl}
        mv {input.log} {output.log}
        """

rule iqtree_inference:
    input:
        msa = MSA_PATH + "/msa.fasta"
    output:
        final_tree = get_inf_output("iqtree") + "/final_tree.newick",
        log = get_inf_output("iqtree") + "/log.txt",
        logl = get_inf_output("iqtree") + "/logl.out"
    params:
        bin = INF_TOOLS["iqtree"]["binary_path"],
        out_dir = lambda wildcards, output: os.path.dirname(output.final_tree),
    shell:
        """
        mkdir -p {params.out_dir}
        {params.bin} -s {input.msa} -m {wildcards.model} --prefix {params.out_dir}/iqtree -nt 1 --seed {wildcards.seed}
        mv {params.out_dir}/iqtree.treefile {output.final_tree}
        mv {params.out_dir}/iqtree.log {output.log}
        grep "Optimal log-likelihood:" {output.log} | sed 's/.*: //' > {output.logl}
        """

rule calculate_distances:
    input:
        true_tree = MSA_PATH + "/tree.nwk",
        final_tree = INF_PATH + "/final_tree.newick",
    output:
        dist_file = INF_PATH + "/distances.json"
    params:
        script = "scripts/calculate_distances.py",
        py_bin = PYTHON,
        start_tree = lambda wildcards, input: os.path.join(os.path.dirname(input.final_tree), "start_tree.newick")
    shell:
        """
        CMD="{params.py_bin} {params.script} --true-tree {input.true_tree} --final-tree {input.final_tree} --output {output.dist_file}"
        if [ -f "{params.start_tree}" ]; then
            CMD="$CMD --start-tree {params.start_tree}"
        fi
        $CMD
        """

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

rule aggregate_msas:
    input:
        [f"{d}/msa.fasta" for d in get_all_msa_dirs()]
    output:
        tsv_path = "results/msa_summary.tsv"
    params:
        msa_dirs = get_all_msa_dirs(),
        sn_config = config
    script:
        "scripts/aggregate_msas.py"

rule aggregate_summary:
    input:
        msa_summary = "results/msa_summary.tsv",
        dist_files = [f"{d}/distances.json" for d in get_all_dirs(INF_PATH)],
        time_files = [f"{d}/time.txt" for d in get_all_dirs(INF_PATH)],
        logl_files = [f"{d}/logl.out" for d in get_all_dirs(INF_PATH)]
    output:
        tsv_path = "results/summary.tsv"
    params:
        dirs = get_all_dirs(INF_PATH),
        sn_config = config
    script:
        "scripts/aggregate_results.py"

rule summary_to_obsidian:
    input:
        tsv = "results/summary.tsv"
    output:
        md = "/Users/mrzi/Seafile/phd_obsidian/notes/pipeline_out/summary_table.md"
    params:
        script = "scripts/tsv_to_md.py",
        py_bin = PYTHON
    shell:
        "{params.py_bin} {params.script} {input.tsv} {output.md}"
