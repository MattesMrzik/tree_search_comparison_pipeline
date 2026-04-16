configfile: "config.yaml"

# from viz.utils import load_snakemake_config_yaml
# cfg = load_snakemake_config_yaml()
# print(cfg)

import os
from snakemake_helpers import infer_wildcard_constraints, make_targets, get_tree_path, get_msa_output, get_inf_output

SEEDS = config["seeds"]
# Select environment based on config (set by profile)
ENV = config.get("env", "local")
PATHS = config["environments"][ENV]

# Tool Paths from environment-specific config
EVOLVER = PATHS["evolver"]
ROOTER = PATHS["root_tree"]
PYTHON = PATHS["python_bin"]
PYTHON_MODULE = PATHS.get("python_module", "")
SIMULATE_TKF = PATHS["simulate_tkf"]
IQTREE3 = PATHS["iqtree3"]
MODEL_SEARCH_PHYLO = PATHS["model_search_phylo"]
JATI = PATHS["jati"]

CONDA_ENV = PATHS.get("conda_env", "")
LOAD_PYTHON = f"source conda activate {CONDA_ENV}; " if ENV == "hpc" and CONDA_ENV else ""

TREE_PATH = config["tree_sim"]["dir"]
MSA_PATH = config["msa_sim"]["dir"]
MODEL_PARAMS_INF_TOOLS = config["model_param_inf"]["tools"]
TREE_INF_TOOLS = config["tree_inf"]["tools"]
MSA_SIM_TOOLS = config["msa_sim"]["tools"]

# apply globally
wildcard_constraints:
    **infer_wildcard_constraints({
        **config["tree_sim"]["tools"],
        **config["msa_sim"]["tools"],
        **config["tree_inf"]["tools"],
        **config["model_param_inf"]["tools"]
    })

rule all_trees:
    input:
        make_targets(config["tree_sim"]["dir"], "tree")

rule all_msas:
    input:
        make_targets(config["msa_sim"]["dir"] + "/msa.fasta", "tree", "msa"),
        make_targets(config["msa_sim"]["dir"] + "/masa.fasta", "tree", "msa")

rule all_model_infs:
    input:
        make_targets(config["model_param_inf"]["dir"] + "/logl.out", "tree", "msa", "minf"),
        [f for f in make_targets(config["msa_sim"]["dir"] + "/sim_logl.out", "tree", "msa") if "/tkf/" in f]

rule all_tree_infs:
    input:
        make_targets(config["tree_inf"]["dir"] + "/final_tree.nwk", "tree", "msa", "tinf")

rule all_tree_pngs:
    input:
        [p.replace(".nwk", ".png") for p in rules.all_trees.input]

rule all:
    input:
        rules.all_trees.input,
        rules.all_msas.input,
        rules.all_model_infs.input,
        rules.all_tree_infs.input,

#######################################################################
#######################################################################
# TREE
#######################################################################
#######################################################################

rule evolver_tree:
    output:
        tree = get_tree_path("evolver"),
        tree_raw = get_tree_path("evolver") + ".raw",
        tree_wo= get_tree_path("evolver") + ".wo"
    threads: 1
    resources:
        mem_mb=1024
    shadow: "minimal"
    shell:
        """
        printf "2\\n{wildcards.species}\\n1 {wildcards.seed} 1\\n{wildcards.birth} {wildcards.death} {wildcards.sampling_fraction} {wildcards.mutation_rate}\\n0\\n" | {EVOLVER} > /dev/null 2>&1
        mkdir -p $(dirname {output.tree})
        tail -n 1 evolver.out > {output.tree}.raw
        {ROOTER} --i {output.tree_raw} --ow {output.tree} --owo {output.tree_wo}
        """

rule iqtree_tree:
    output:
        tree = get_tree_path("iqtree"),
        tree_raw = get_tree_path("iqtree") + ".raw",
        tree_wo= get_tree_path("iqtree") + ".wo"
    threads: 1
    resources:
        mem_mb=1024
    shell:
        """
        {IQTREE3} -r {wildcards.species} {output.tree}.raw -nt 1 --seed {wildcards.seed} -redo
        {ROOTER} --i {output.tree_raw} --ow {output.tree} --owo {output.tree_wo}
        """

rule tree_png:
    input:
        tree = TREE_PATH
    output:
        plot = TREE_PATH.replace(".nwk", ".png")
    threads: 1
    resources:
        mem_mb=2048
    shell:
        """
        {LOAD_PYTHON}
        {PYTHON} viz/tree/visualize_trees.py --tree-file {input.tree} --output-file {output.plot}
        """
#######################################################################
#######################################################################
# MSA
#######################################################################
#######################################################################

rule simulate_tkf_alignment:
    input:
        tree = TREE_PATH,
        tree_wo = TREE_PATH + ".wo" # it complained about the (); around the
    output:
        msa = get_msa_output("tkf") + "/msa.fasta",
        masa = get_msa_output("tkf") + "/masa.fasta",
        tree_w_internal = get_msa_output("tkf") + "/tree.nwk",
        tree_w_internal_wo = get_msa_output("tkf") + "/tree.nwk.wo"
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
            --max-insertion-length {wildcards.max_insertion} \
            --root-length {wildcards.root_length} \
            --seed {wildcards.seed} \
            --output-dir $(dirname {output.msa})
        cp {input.tree_wo} {output.tree_w_internal_wo}
        """

rule tkf_sim_alignment_logl: 
    input:
        msa = get_msa_output("tkf") + "/masa.fasta",
        tree = get_msa_output("tkf") + "/tree.nwk"
    output:
        logl = get_msa_output("tkf") + "/sim_logl.out",
        log = get_msa_output("tkf") + "/sim_logl_log.txt"
    shell:
        """
        mkdir -p $(dirname {output.logl})
        {MODEL_SEARCH_PHYLO} \
            --out-folder $(dirname {output.logl}) \
            --seq-file {input.msa} \
            --tree-file {input.tree} \
            --model {wildcards.model}\
            --params {wildcards.lambda} {wildcards.mu} {wildcards.r} \
            --gap-handling TKF92 \
            --epsilon inf \
            --freq-opt fixed \
            --seed {wildcards.seed} \
            -l warn 
        mv $(dirname {output.logl})/*/*.out {output.logl}
        mv $(dirname {output.logl})/*/*.log {output.log}
        rm -rf $(dirname {output.logl})/*/
        """

rule simulate_alisim_alignment:
    input:
        tree = TREE_PATH,
        tree_wo = TREE_PATH + ".wo" # it complained about the (); around the whole newick tree so we use this instead
    output:
        msa = get_msa_output("alisim") + "/msa.fasta",
        tree_w_internal = get_msa_output("alisim") + "/tree.nwk",
        tree_w_internal_wo = get_msa_output("alisim") + "/tree.nwk.wo"
    threads: 1
    resources:
        mem_mb=1024
    shell:
        """
        mkdir -p $(dirname {output.msa})
        {IQTREE3} \
            --alisim $(dirname {output.msa})/msa \
            -m {MSA_SIM_TOOLS[alisim][model]} \
            -t {input.tree_wo} \
            --indel {wildcards.ir},{wildcards.ip} \
            --length {wildcards.root_length} \
            --seed {wildcards.seed} \
            --out-format fasta \
            --no-unaligned 
        mv $(dirname {output.msa})/msa.fa {output.msa}
        cp {input.tree} {output.tree_w_internal}
        cp {input.tree_wo} {output.tree_w_internal_wo}
        """

rule simulate_alisim_ancestral_alignment:
    input:
        tree = TREE_PATH,
        tree_wo = TREE_PATH + ".wo" # it complained about the (); around the whole newick tree so we use this instead
    output:
        msa = get_msa_output("alisim") + "/masa.fasta",
    threads: 1
    resources:
        mem_mb=1024
    shell:
        """
        mkdir -p $(dirname {output.msa})
        {IQTREE3} \
            --alisim $(dirname {output.msa})/masa \
            -m {MSA_SIM_TOOLS[alisim][model]} \
            -t {input.tree_wo} \
            --indel {wildcards.ir},{wildcards.ip} \
            --length {wildcards.root_length} \
            --seed {wildcards.seed} \
            --out-format fasta \
            --no-unaligned \
            --write-all
        mv $(dirname {output.msa})/masa.fa {output.msa}
        """

#######################################################################
#######################################################################
# MODEL PARAMETER INFERENCE
#######################################################################
#######################################################################

rule jati_model_param_search:
    input:
        masa = MSA_PATH + "/masa.fasta", # using this for TKF
        msa = MSA_PATH + "/msa.fasta", # using this for non-TKF
        tree = MSA_PATH + "/tree.nwk"
    output:
        logl = get_inf_output("model_param_inf", "jati_model_param_search") + "/logl.out",
        log = get_inf_output("model_param_inf", "jati_model_param_search") + "/log.txt",
        params = get_inf_output("model_param_inf", "jati_model_param_search") + "/params.json"
    threads: 1
    resources:
        mem_mb=4096
    params:
        epsilon = MODEL_PARAMS_INF_TOOLS["jati_model_param_search"]["epsilon"],
        seq_file = lambda wc, input: input.masa if wc.gap == "TKF92" else input.msa,
        paras = " ".join(map(str, MODEL_PARAMS_INF_TOOLS["jati_model_param_search"]["params"])),
        out_base = get_inf_output("model_param_inf", "jati_model_param_search")
    shell:
        """
        mkdir -p {params.out_base}
        rm -rf {params.out_base}
        {MODEL_SEARCH_PHYLO} \
            --out-folder {params.out_base} \
            --seq-file {params.seq_file} \
            --tree-file {input.tree} \
            --model {wildcards.model} \
            --params {params.paras} \
            --gap-handling {wildcards.gap} \
            --epsilon {params.epsilon} \
            --seed {wildcards.seed} \
            -l warn 
        mv {params.out_base}/*/*.out {output.logl}
        mv {params.out_base}/*/*.log {output.log}
        mv {params.out_base}/*/*.json {output.params}
        rm -rf {params.out_base}/*/
        """

rule iqtree_model_param_search:
    input:
        msa = MSA_PATH + "/msa.fasta",
        tree = MSA_PATH + "/tree.nwk.wo"
    output:
        logl = get_inf_output("model_param_inf", "iqtree_model_param_search") + "/logl.out",
        log = get_inf_output("model_param_inf", "iqtree_model_param_search") + "/log.txt"
    threads: 1
    resources:
        mem_mb=4096
    shell:
        """
        mkdir -p $(dirname {output.logl})
        {IQTREE3} -s {input.msa} -m {wildcards.model} -t {input.tree} --prefix $(dirname {output.logl})/iqtree_model_search -nt 1 --seed {wildcards.seed} -redo -blfix
        mv $(dirname {output.logl})/*.log {output.log}
        grep "Optimal log-likelihood:" {output.log} | sed 's/.*: //' > {output.logl}
        rm $(dirname {output.logl})/iqtree_model_search*
        """

#######################################################################
#######################################################################
# TREE INFERENCE
#######################################################################
#######################################################################

# TODO also ouput params and msa

rule jati_inference:
    input:
        msa = MSA_PATH + "/msa.fasta"
    output:
        start_tree = get_inf_output("tree_inf", "jati") + "/start_tree.nwk",
        final_tree = get_inf_output("tree_inf", "jati") + "/final_tree.nwk",
        logl = get_inf_output("tree_inf", "jati") + "/logl.out",
        log = get_inf_output("tree_inf", "jati") + "/log.txt"
    threads: 1
    resources:
        mem_mb=4096
    params:
        paras = " ".join(map(str, TREE_INF_TOOLS["jati"]["params"])),
        max_iterations = TREE_INF_TOOLS["jati"]["max_iterations"],
        out_base = get_inf_output("tree_inf", "jati"),
        force_nni = lambda wildcards: "--force-nni" if wildcards.move == "NNI" else ""
    shell:
        """
        mkdir -p {params.out_base}
        rm -rf {params.out_base}
        {JATI} \
            --out-folder {params.out_base} \
            --seq-file {input.msa} \
            --model {wildcards.model} \
            --params {params.paras} \
            --gap-handling {wildcards.gap} \
            {params.force_nni} \
            --seed {wildcards.seed} \
            -l warn \
            --max-iterations {params.max_iterations}
        mv {params.out_base}/*/*_start_tree.newick {output.start_tree}
        mv {params.out_base}/*/*_tree.newick {output.final_tree}
        mv {params.out_base}/*/*_logl.out {output.logl}
        mv {params.out_base}/*/*.log {output.log}
        rm -rf {params.out_base}/*/
        """

rule iqtree_inference:
    input:
        msa = MSA_PATH + "/msa.fasta"
    output:
        final_tree = get_inf_output("tree_inf", "iqtree") + "/final_tree.nwk",
        log = get_inf_output("tree_inf", "iqtree") + "/log.txt",
        logl = get_inf_output("tree_inf", "iqtree") + "/logl.out"
    threads: 1
    resources:
        mem_mb=4096
    shell:
        """
        mkdir -p $(dirname {output.final_tree})
        {IQTREE3} -s {input.msa} -m {wildcards.model} --prefix $(dirname {output.final_tree})/iqtree -nt 1 --seed {wildcards.seed} -redo
        mv $(dirname {output.final_tree})/iqtree.treefile {output.final_tree}
        mv $(dirname {output.final_tree})/iqtree.log {output.log}
        grep "Optimal log-likelihood:" {output.log} | sed 's/.*: //' > {output.logl}
        """
