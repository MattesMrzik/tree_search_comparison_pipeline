configfile: "config.yaml"

import os
import glob
from datetime import datetime

# Path to tools
EVOLVER = config["evolver_path"]
SIMULATE_TKF = config["simulate_tkf_path"]
PYTHON = config["python_bin"]
JATI = config["jati_path"]

# PAML evolver Parameters
SPECIES = config["species"]
BIRTH_DEATH_PAIRS = config["birth_death_rates"]
SEEDS = config["seeds"]
SAMPLING = config["sampling_fraction"]
MUTATION = config["mutation_rate"]

# JATI Parameters
JATI_MODELS = config["jati_models"]
JATI_PARAS = config["jati_paras"]
GAP_STRATEGIES = config["gap_handling_strategies"]
MAX_ITERATIONS = config["max_iterations"]

# TKF92 Parameters
LAMBDA = config["tkf_lambda"]
MU = config["tkf_mu"]
R = config["tkf_r"]
MAX_INS = config["max_insertion_length"]

# Path Templates from config
SIM_DIR = config["sim_dir"]
INF_DIR = config["inf_dir"]

# Target Generator
def get_all_inference_dirs():
    return expand(INF_DIR, 
                  s=SPECIES, 
                  b=[p[0] for p in BIRTH_DEATH_PAIRS], 
                  d=[p[1] for p in BIRTH_DEATH_PAIRS], 
                  f=SAMPLING, m=MUTATION, seed=SEEDS, 
                  model=JATI_MODELS, gap=GAP_STRATEGIES)

rule all:
    input:
        "results/summary.tsv",
        expand(f"{SIM_DIR}/tree_plot.png", 
               s=SPECIES, 
               b=[p[0] for p in BIRTH_DEATH_PAIRS], 
               d=[p[1] for p in BIRTH_DEATH_PAIRS], 
               f=SAMPLING, m=MUTATION, seed=SEEDS)

rule generate_tree:
    output:
        "results/trees/s{s}_b{b}_d{d}_f{f}_m{m}_seed{seed}.nwk"
    params:
        evolver = EVOLVER
    shadow: "minimal"
    shell:
        """
        # Run evolver. Because of shadow: "minimal", evolver.out 
        # is created in an isolated directory for this specific job.
        printf "2\\n{wildcards.s}\\n1 {wildcards.seed} 1\\n{wildcards.b} {wildcards.d} {wildcards.f} {wildcards.m}\\n0\\n" | {params.evolver} > /dev/null 2>&1
        
        # Extract the tree from the local, isolated evolver.out
        tail -n 1 evolver.out > {output}
        """

rule simulate_tkf:
    input:
        tree = "results/trees/s{s}_b{b}_d{d}_f{f}_m{m}_seed{seed}.nwk"
    output:
        msa = f"{SIM_DIR}/msa.fasta",
        masa = f"{SIM_DIR}/masa.fasta",
        info = f"{SIM_DIR}/info.txt",
        tree_copy = f"{SIM_DIR}/tree.nwk"
    params:
        bin = SIMULATE_TKF,
        lambda_val = LAMBDA,
        mu_val = MU,
        r_val = R,
        max_ins = MAX_INS,
        out_dir = f"{SIM_DIR}"
    shell:
        """
        {params.bin} \
            --tree-file {input.tree} \
            --lambda {params.lambda_val} \
            --mu {params.mu_val} \
            --r {params.r_val} \
            --max-insertion-length {params.max_ins} \
            --seed {wildcards.seed} \
            --output-dir {params.out_dir}
        
        # Ensure the tree file is named tree.nwk in the output directory
        cp {input.tree} {output.tree_copy}
        """

rule visualize_msa_tree:
    input:
        tree = f"{SIM_DIR}/tree.nwk"
    output:
        plot = f"{SIM_DIR}/tree_plot.png"
    params:
        script = "scripts/visualize_trees.py",
        py_bin = PYTHON
    shell:
        """
        {params.py_bin} {params.script} --tree-file {input.tree} --output-file {output.plot}
        """

# We use this checkpoint since jati's output directory is not deterministic and we need to wait for it to be created before we can move files around
checkpoint jati_inference:
    input:
        msa = f"{SIM_DIR}/msa.fasta"
    output:
        outdir = directory(f"{INF_DIR}/out")
    params:
        bin = JATI,
        paras = " ".join(map(str, JATI_PARAS)),
        log_level = "warn",
        max_iterations = MAX_ITERATIONS
    shell:
        """
        mkdir -p {output.outdir}
        {params.bin} \
            --out-folder {output.outdir} \
            --seq-file {input.msa} \
            --model {wildcards.model} \
            --params {params.paras} \
            --gap-handling {wildcards.gap} \
            --seed {wildcards.seed} \
            -l {params.log_level} \
            --max-iterations {params.max_iterations} 
        """

def get_jati_output(wildcards):
    outdir = checkpoints.jati_inference.get(**wildcards).output.outdir
    dirs = glob.glob(f"{outdir}/*_out/")
    if not dirs:
        raise ValueError(f"JATI failed to produce an output directory in {outdir}")
    return dirs[0]

rule jati_cleanup:
    input:
        dir = get_jati_output
    output:
        start_tree = f"{INF_DIR}/start_tree.newick",
        final_tree = f"{INF_DIR}/final_tree.newick",
        logl = f"{INF_DIR}/logl.out",
        log = f"{INF_DIR}/log.txt"
    params:
        target_dir = f"{INF_DIR}"
    shell:
        """
        mv {input.dir}/* {params.target_dir}/

        # striping the time stamp prefix
        for f in {params.target_dir}/[0-9]*_*; do
            base="${{f##*/}}"
            mv "$f" "{params.target_dir}/${{base#*_}}"
        done
        
        mv {params.target_dir}/tree.newick {output.final_tree}
        mv {params.target_dir}/*.log {output.log}
        """

rule calculate_distances:
    input:
        true_tree = f"{SIM_DIR}/tree.nwk",
        final_tree = f"{INF_DIR}/final_tree.newick"
    output:
        dist_file = f"{INF_DIR}/distances.json"
    params:
        script = "scripts/calculate_distances.py",
        py_bin = PYTHON
    shell:
        "{params.py_bin} {params.script} --tree1 {input.true_tree} --tree2 {input.final_tree} --output {output.dist_file}"

rule calculate_time:
    input:
        log = f"{INF_DIR}/log.txt"
    output:
        time_file = f"{INF_DIR}/time.txt"
    params:
        script = "scripts/calculate_time.py",
        py_bin = PYTHON
    shell:
        "{params.py_bin} {params.script} {input.log} {output.time_file}"

rule aggregate_summary:
    input:
        [f"{d}/distances.json" for d in get_all_inference_dirs()],
        [f"{d}/time.txt" for d in get_all_inference_dirs()],
        [f"{d}/logl.out" for d in get_all_inference_dirs()]
    output:
        tsv_path = "results/summary.tsv"
    params:
        dirs = get_all_inference_dirs(),
        sn_config = config
    script:
        "scripts/aggregate_results.py"
