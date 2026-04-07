configfile: "config.yaml"

import os

# Select environment based on config (set by profile)
ENV = config.get("env", "local")
PATHS = config["environments"][ENV]

# Tool Paths from environment-specific config
EVOLVER = PATHS["evolver"]
PYTHON = PATHS["python_bin"]
PYTHON_MODULE = PATHS.get("python_module", "")
SIMULATE_TKF = PATHS["simulate_tkf"]
IQTREE3 = PATHS["iqtree3"]
MODEL_SEARCH_PHYLO = PATHS["model_search_phylo"]
JATI = PATHS["jati"]

# Helper to load python module on HPC
LOAD_PYTHON = f"module load {PYTHON_MODULE}; " if ENV == "hpc" and PYTHON_MODULE else ""

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

def get_all_inference_dirs(template):
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
        
        for inf_tool_name, inf_conf in INF_TOOLS.items():
            if inf_tool_name == "jati":
                inf_params = []
                for gap, move in inf_conf["gap_move"]:
                    inf_params.extend(expand(inf_conf["path_snippet"], 
                                             model=inf_conf["models"], 
                                             gap=[gap],
                                             move_strategy=[move]))
            elif inf_tool_name == "true_tree":
                inf_params = expand(inf_conf["path_snippet"], 
                                   model=inf_conf["models"],
                                   gap=inf_conf["gap_strategies"])
            elif inf_tool_name == "iqtree":
                inf_params = expand(inf_conf["path_snippet"], 
                                    model=inf_conf["models"])

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
        [f"{d}/distances.json" for d in get_all_inference_dirs(INF_PATH)],

rule simulate_alignments:
    input:
        [f"{d}/msa.fasta" for d in get_all_inference_dirs(MSA_PATH)]

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
    threads: 1
    resources:
        mem_mb=1024
    shadow: "minimal"
    shell:
        """
        {LOAD_PYTHON}
        printf "2\\n{wildcards.s}\\n1 {wildcards.seed} 1\\n{wildcards.b} {wildcards.d} {wildcards.f} {wildcards.m}\\n0\\n" | {EVOLVER} > /dev/null 2>&1
        tail -n 1 evolver.out > {output}
        """

rule visualize_tree:
    input:
        tree = TREE_PATH
    output:
        plot = f"{TREE_PATH}.png"
    threads: 1
    resources:
        mem_mb=2048
    shell:
        """
        {LOAD_PYTHON}
        {PYTHON} scripts/visualize_trees.py --tree-file {input.tree} --output-file {output.plot}
        """

# Helper functions for rule-specific path generation
def get_msa_output(tool_name):
    tool_conf = MSA_SIM_TOOLS[tool_name]
    return MSA_PATH.replace("{msa_sim_tool}", tool_name).replace(
        "{tool_params}", tool_conf["params_path_snipped"]
    )

def get_inf_output(tool_name):
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
    threads: 1
    resources:
        mem_mb=4096
    shell:
        """
        {SIMULATE_TKF} \
            --tree-file {input.tree} \
            --lambda {wildcards.lambda} \
            --mu {wildcards.mu} \
            --r {wildcards.r} \
            --max-insertion-length {wildcards.max_ins} \
            --root-length {wildcards.root_length} \
            --seed {wildcards.seed} \
            --output-dir $(dirname {output.msa})
        cp {input.tree} {output.tree_copy}
        """

rule simulate_alisim_alignment:
    input:
        tree = TREE_PATH
    output:
        msa = get_msa_output("alisim") + "/msa.fasta",
        tree_copy = get_msa_output("alisim") + "/tree.nwk"
    threads: 1
    resources:
        mem_mb=1024
    shell:
        """
        mkdir -p $(dirname {output.msa})
        {IQTREE3} \
            --alisim $(dirname {output.msa})/msa \
            -m {MSA_SIM_TOOLS[alisim][model]} \
            -t {input.tree} \
            --indel {wildcards.ir},{wildcards.ip} \
            --length {wildcards.root_length} \
            --seed {wildcards.seed} \
            --out-format fasta \
            --no-unaligned 
        mv $(dirname {output.msa})/msa.fa {output.msa}
        cp {input.tree} {output.tree_copy}
        """

rule true_tree_inference:
    input:
        msa = MSA_PATH + "/msa.fasta",
        tree = MSA_PATH + "/tree.nwk"
    output:
        final_tree = get_inf_output("true_tree") + "/final_tree.newick",
        logl = get_inf_output("true_tree") + "/logl.out",
        log = get_inf_output("true_tree") + "/log.txt"
    threads: 1
    resources:
        mem_mb=4096
    params:
        epsilon = INF_TOOLS["true_tree"]["epsilon"],
        paras = " ".join(map(str, INF_TOOLS["true_tree"]["params"])),
        out_base = get_inf_output("true_tree")
    shell:
        """
        mkdir -p {params.out_base}
        {MODEL_SEARCH_PHYLO} \
            --out-folder {params.out_base} \
            --seq-file {input.msa} \
            --tree-file {input.tree} \
            --model {wildcards.model} \
            --params {params.paras} \
            --gap-handling {wildcards.gap} \
            --epsilon {params.epsilon} \
            --seed {wildcards.seed} \
            -l warn --no-timestamp
        mv {params.out_base}/jati_run_out/jati_run_tree.newick {output.final_tree}
        mv {params.out_base}/jati_run_out/jati_run_logl.out {output.logl}
        mv {params.out_base}/jati_run_out/jati_run.log {output.log}
        """

rule jati_inference:
    input:
        msa = MSA_PATH + "/msa.fasta"
    output:
        start_tree = get_inf_output("jati") + "/start_tree.newick",
        final_tree = get_inf_output("jati") + "/final_tree.newick",
        logl = get_inf_output("jati") + "/logl.out",
        log = get_inf_output("jati") + "/log.txt"
    threads: 1
    resources:
        mem_mb=4096
    params:
        paras = " ".join(map(str, INF_TOOLS["jati"]["params"])),
        max_iterations = INF_TOOLS["jati"]["max_iterations"],
        out_base = get_inf_output("jati"),
        force_nni = lambda wildcards: "--force-nni" if wildcards.move_strategy == "NNI" else ""
    shell:
        """
        mkdir -p {params.out_base}
        {JATI} \
            --out-folder {params.out_base} \
            --seq-file {input.msa} \
            --model {wildcards.model} \
            --params {params.paras} \
            --gap-handling {wildcards.gap} \
            {params.force_nni} \
            --seed {wildcards.seed} \
            -l warn \
            --max-iterations {params.max_iterations} \
            --no-timestamp
        mv {params.out_base}/jati_run_out/jati_run_start_tree.newick {output.start_tree}
        mv {params.out_base}/jati_run_out/jati_run_tree.newick {output.final_tree}
        mv {params.out_base}/jati_run_out/jati_run_logl.out {output.logl}
        mv {params.out_base}/jati_run_out/jati_run.log {output.log}
        """

rule iqtree_inference:
    input:
        msa = MSA_PATH + "/msa.fasta"
    output:
        final_tree = get_inf_output("iqtree") + "/final_tree.newick",
        log = get_inf_output("iqtree") + "/log.txt",
        logl = get_inf_output("iqtree") + "/logl.out"
    threads: 1
    resources:
        mem_mb=4096
    shell:
        """
        mkdir -p $(dirname {output.final_tree})
        {IQTREE3} -s {input.msa} -m {wildcards.model} --prefix $(dirname {output.final_tree})/iqtree -nt 1 --seed {wildcards.seed}
        mv $(dirname {output.final_tree})/iqtree.treefile {output.final_tree}
        mv $(dirname {output.final_tree})/iqtree.log {output.log}
        grep "Optimal log-likelihood:" {output.log} | sed 's/.*: //' > {output.logl}
        """

rule calculate_distances:
    input:
        true_tree = MSA_PATH + "/tree.nwk",
        final_tree = INF_PATH + "/final_tree.newick"
    output:
        dist_file = INF_PATH + "/distances.json"
    threads: 1
    resources:
        mem_mb=1024
    shell:
        """
        {LOAD_PYTHON}
        START_TREE=$(dirname {input.final_tree})/start_tree.newick
        CMD="{PYTHON} scripts/calculate_distances.py --true-tree {input.true_tree} --final-tree {input.final_tree} --output {output.dist_file}"
        if [ -f "$START_TREE" ]; then CMD="$CMD --start-tree $START_TREE"; fi
        $CMD
        """

rule calculate_time:
    input:
        log = INF_PATH + "/log.txt"
    output:
        time_file = INF_PATH + "/time.txt"
    threads: 1
    resources:
        mem_mb=512
    shell:
        """
        {LOAD_PYTHON}
        {PYTHON} scripts/calculate_time.py {input.log} {output.time_file}
        """

rule summarize_msas:
    input:
        [f"{d}/msa.fasta" for d in get_all_msa_dirs()]
    output:
        tsv_path = "results/msa_summary.tsv"
    threads: 1
    resources:
        mem_mb=2048
    params:
        msa_dirs = get_all_msa_dirs(),
        sn_config = config
    script:
        "scripts/summarize_msas.py"

rule summarize_results:
    input:
        msa_summary = "results/msa_summary.tsv",
        dist_files = [f"{d}/distances.json" for d in get_all_inference_dirs(INF_PATH)],
        time_files = [f"{d}/time.txt" for d in get_all_inference_dirs(INF_PATH)],
        logl_files = [f"{d}/logl.out" for d in get_all_inference_dirs(INF_PATH)]
    output:
        tsv_path = "results/summary.tsv"
    threads: 1
    resources:
        mem_mb=2048
    params:
        dirs = get_all_inference_dirs(INF_PATH),
        sn_config = config
    script:
        "scripts/summarize_results.py"

rule obsidian:
    input:
        tsv = "results/summary.tsv"
    output:
        md = "/Users/mrzi/Seafile/phd_obsidian/notes/pipeline_out/summary_table.md"
    threads: 1
    resources:
        mem_mb=512
    shell:
        """
        {LOAD_PYTHON}
        {PYTHON} scripts/tsv_to_md.py {input.tsv} {output.md}
        """
