configfile: "config.yaml"

# from viz.utils import load_snakemake_config_yaml
# cfg = load_snakemake_config_yaml()
# print(cfg)

import os
from snakemake_helpers import infer_wildcard_constraints, make_targets, get_tree_path, get_msa_output, get_inf_output, get_inf_output_with_msa_params

SEEDS = config["seeds"]
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
INDEL_ASR = PATHS["indel_asr"]

CONDA_ENV = PATHS.get("conda_env", "")
LOAD_PYTHON = f"source conda activate {CONDA_ENV}; " if ENV == "hpc" and CONDA_ENV else ""

TREE_PATH = config["tree_sim"]["dir"]
MSA_PATH = config["msa_sim"]["dir"]

wildcard_constraints:
    **infer_wildcard_constraints(config)

rule all_trees:
    input:
        make_targets(config, primary="tree_sim")

rule all_msas:
    input:
        make_targets(config, "tree_sim", primary="msa_sim", suffix="/msa.fasta"),
        make_targets(config, "tree_sim", primary="msa_sim", suffix="/masa.fasta"),
        [f for f in make_targets(config, "tree_sim", primary="msa_sim", suffix="/sim_logl.out") if "/tkf/" in f],
        [f for f in make_targets(config, "tree_sim", primary="msa_sim", suffix="/sim_indel_logl.out") if "/tkf/" in f]

rule all_param_infs:
    input:
        make_targets(config, "tree_sim", "msa_sim", primary="param_inf", suffix="/logl.out")

rule all_tree_infs:
    input:
        make_targets(config, "tree_sim", "msa_sim", primary="tree_inf", suffix="/final_tree.nwk")

rule all_tree_pngs:
    input:
        [p.replace(".nwk", ".png") for p in rules.all_trees.input]

indel = make_targets(config, "tree_sim", "msa_sim", primary="indel_inf")
indel_and_params = make_targets(config, "tree_sim", "msa_sim", primary="indel_and_param_inf")
rule all_tkf_msas_indel_infs:
    input:
        [t + "/masa.fasta" for t in indel if "/tkf/" in t],

rule all_tkf_msas_indel_and_param_infs:
    input:
        [t + "/masa.fasta" for t in indel_and_params if "/tkf/" in t],
        [t + "/params.json" for t in indel_and_params if "/tkf/" in t]

rule all_sanity:
    input:
        rules.all_msas.input,
        rules.all_param_infs.input,
        rules.all_tkf_msas_indel_infs.input,
        rules.all_tkf_msas_indel_and_param_infs.input

rule all_indel_and_param_infs:
    input:
        [t + "/masa.fasta" for t in indel_and_params],
        [t + "/params.json" for t in indel_and_params]

#######################################################################
#######################################################################
# TREE
#######################################################################
#######################################################################

rule evolver_tree:
    priority: lambda wildcards: -int(wildcards.seed)
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
    priority: lambda wildcards: -int(wildcards.seed)
    output:
        tree = get_tree_path("iqtree"),
        tree_raw = get_tree_path("iqtree") + ".raw",
        tree_wo= get_tree_path("iqtree") + ".wo"
    threads: 1
    resources:
        mem_mb=1024
    shell:
        """
        {IQTREE3} -r {wildcards.species} --rlen {wildcards.min} {wildcards.mean} {wildcards.max} {output.tree}.raw -nt 1 --seed {wildcards.seed} -redo
        {ROOTER} --i {output.tree_raw} --ow {output.tree} --owo {output.tree_wo}
        """

rule tree_png:
    priority: lambda wildcards: -int(wildcards.seed)
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
        {PYTHON} viz/sim/tree/visualize_trees.py --tree-file {input.tree} --output-file {output.plot}
        """

#######################################################################
#######################################################################
# MSA
#######################################################################
#######################################################################

rule simulate_tkf_alignment:
    priority: lambda wildcards: -int(wildcards.seed)
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
            --root-length {wildcards.tkf_root_length} \
            --seed {wildcards.seed} \
            --output-dir $(dirname {output.msa})
        cp {input.tree_wo} {output.tree_w_internal_wo}
        """

rule tkf_sim_alignment_logl: 
    priority: lambda wildcards: -int(wildcards.seed)
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
            --out-folder $(dirname {output.logl})/with_subst \
            --seq-file {input.msa} \
            --tree-file {input.tree} \
            --model {wildcards.model}\
            --params {wildcards.lambda} {wildcards.mu} {wildcards.r} \
            --gap-handling TKF92 \
            --epsilon inf \
            --freq-opt fixed \
            --seed {wildcards.seed} \
            -l warn 
        mv $(dirname {output.logl})/with_subst/*/*.out {output.logl}
        mv $(dirname {output.logl})/with_subst/*/*.log {output.log}
        rm -rf $(dirname {output.logl})/with_subst/
        """

rule tkf_sim_alignment_indel_logl: 
    priority: lambda wildcards: -int(wildcards.seed)
    input:
        msa = get_msa_output("tkf") + "/masa.fasta",
        tree = get_msa_output("tkf") + "/tree.nwk"
    output:
        logl = get_msa_output("tkf") + "/sim_indel_logl.out",
        log = get_msa_output("tkf") + "/sim_indel_logl_log.txt"
    shell:
        """
        mkdir -p $(dirname {output.logl})
        {MODEL_SEARCH_PHYLO} \
            --out-folder $(dirname {output.logl})/only_indels \
            --seq-file {input.msa} \
            --tree-file {input.tree} \
            --model none \
            --params {wildcards.lambda} {wildcards.mu} {wildcards.r} \
            --gap-handling TKF92 \
            --epsilon inf \
            --freq-opt fixed \
            --seed {wildcards.seed} \
            -l warn 
        mv $(dirname {output.logl})/only_indels/*/*.out {output.logl}
        mv $(dirname {output.logl})/only_indels/*/*.log {output.log}
        rm -rf $(dirname {output.logl})/only_indels/
        """

rule simulate_alisim_alignment:
    priority: lambda wildcards: -int(wildcards.seed)
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
            -m {wildcards.model} \
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
    priority: lambda wildcards: -int(wildcards.seed)
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
            -m {wildcards.model} \
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
    priority: lambda wildcards: -int(wildcards.seed)
    input:
        masa = MSA_PATH + "/masa.fasta", # using this for TKF
        msa = MSA_PATH + "/msa.fasta", # using this for non-TKF
        tree = MSA_PATH + "/tree.nwk"
    output:
        logl = get_inf_output("param_inf", "jati_model_param_search") + "/logl.out",
        log = get_inf_output("param_inf", "jati_model_param_search") + "/log.txt",
        params = get_inf_output("param_inf", "jati_model_param_search") + "/params.json"
    threads: 1
    resources:
        mem_mb=4096
    params:
        seq_file = lambda wc, input: input.masa if wc.gap == "TKF92" else input.msa,
        out_base = get_inf_output("param_inf", "jati_model_param_search")
    shell:
        """
        mkdir -p {params.out_base}
        rm -rf {params.out_base}
        {MODEL_SEARCH_PHYLO} \
            --out-folder {params.out_base} \
            --seq-file {params.seq_file} \
            --tree-file {input.tree} \
            --model {wildcards.model} \
            --gap-handling {wildcards.gap} \
            --epsilon {wildcards.epsilon} \
            --seed {wildcards.seed} \
            -l warn 
        mv {params.out_base}/*/*.out {output.logl}
        mv {params.out_base}/*/*.log {output.log}
        mv {params.out_base}/*/*.json {output.params}
        rm -rf {params.out_base}/*/
        """

rule iqtree_model_param_search:
    priority: lambda wildcards: -int(wildcards.seed)
    input:
        msa = MSA_PATH + "/msa.fasta",
        tree = MSA_PATH + "/tree.nwk.wo"
    output:
        logl = get_inf_output("param_inf", "iqtree_model_param_search") + "/logl.out",
        log = get_inf_output("param_inf", "iqtree_model_param_search") + "/log.txt"
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
# ASR 
#######################################################################
#######################################################################

rule tkf_indel_asr:
    priority: lambda wildcards: -int(wildcards.seed)
    input:
        msa = get_msa_output("tkf") + "/msa.fasta",
        tree = get_msa_output("tkf") + "/tree.nwk"
    output:
        asr = get_inf_output_with_msa_params("indel_inf", "tkf_reestimate", "tkf") + "/masa.fasta"
    threads: 1
    resources:
        mem_mb=4096
    shell:
        """
        rm -rf $(dirname {output.asr})
        mkdir -p $(dirname {output.asr})
        # the --params are from the simulation
        {INDEL_ASR} \
            --seq-file {input.msa} \
            --tree-file {input.tree} \
            --out-folder $(dirname {output.asr}) \
            --algorithm-type tkf92 \
            --params {wildcards.lambda} {wildcards.mu} {wildcards.r} \
            --seed {wildcards.seed} \
            --max-iterations {wildcards.max_iterations} \
            -l warn \
            --epsilon {wildcards.epsilon} 
        """

rule parsimony_indel_asr:
    priority: lambda wildcards: -int(wildcards.seed)
    input:
        msa = get_msa_output("tkf") + "/msa.fasta",
        tree = get_msa_output("tkf") + "/tree.nwk"
        # we only use tkf alignments here since we want to calculate the logl of the parsimony ASR under the true params
    output:
        asr = get_inf_output_with_msa_params("indel_inf", "dollo_parsimony", "tkf") + "/masa.fasta"
    threads: 1
    resources:
        mem_mb=4096
    shell:
        """
        rm -rf $(dirname {output.asr})
        mkdir -p $(dirname {output.asr})
        # the --params are used to calculate the lolg of the parsimony ASR under the true params, i.e., the ones used for the simulation
        {INDEL_ASR} \
            --seq-file {input.msa} \
            --tree-file {input.tree} \
            --out-folder $(dirname {output.asr}) \
            --params {wildcards.lambda} {wildcards.mu} {wildcards.r} \
            --algorithm-type parsimony \
            -l warn \
            --seed {wildcards.seed} 
        """

rule tkf_indel_asr_and_params:
    priority: lambda wildcards: -int(wildcards.seed)
    input:
        msa = MSA_PATH + "/msa.fasta",
        tree = MSA_PATH + "/tree.nwk"
    output:
        asr = get_inf_output("indel_and_param_inf", "jati_asr_and_params") + "/masa.fasta",
        params = get_inf_output("indel_and_param_inf", "jati_asr_and_params") + "/params.json"
    threads: 1
    resources:
        mem_mb=4096
    shell:
        """
        rm -rf $(dirname {output.asr})
        mkdir -p $(dirname {output.asr})
        {INDEL_ASR} \
            --seq-file {input.msa} \
            --tree-file {input.tree} \
            --out-folder $(dirname {output.asr}) \
            --algorithm-type tkf92 \
            --params {wildcards.s_l} {wildcards.s_m} {wildcards.s_r} \
            --seed {wildcards.seed} \
            -l warn \
            --max-iterations {wildcards.max_iterations} \
            --epsilon {wildcards.epsilon} \
            -o # optimize params
        """

#######################################################################
#######################################################################
# TREE INFERENCE
#######################################################################
#######################################################################

# TODO also ouput params and msa

rule jati_inference:
    priority: lambda wildcards: -int(wildcards.seed)
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
        out_base = get_inf_output("tree_inf", "jati"),
        force_nni = lambda wildcards: "--force-nni" if wildcards.move == "NNI" else ""
    shell:
        """
        rm -rf {params.out_base}
        mkdir -p {params.out_base}
        {JATI} \
            --out-folder {params.out_base} \
            --seq-file {input.msa} \
            --model {wildcards.model} \
            --gap-handling {wildcards.gap} \
            {params.force_nni} \
            --seed {wildcards.seed} \
            -l warn \
            --max-iterations {wildcards.max_iterations} \
            --epsilon {wildcards.epsilon}
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
