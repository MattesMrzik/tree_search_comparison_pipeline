configfile: "config.yaml"

import os

# Path to tools
EVOLVER = config["evolver_path"]
SIMULATE_TKF = config["simulate_tkf_path"]
PYTHON = config["python_bin"]

# Extracting parameters
SPECIES = config["species"]
BIRTH_DEATH_PAIRS = config["birth_death_rates"]
SEEDS = config["seeds"]
SAMPLING = config["sampling_fraction"]
MUTATION = config["mutation_rate"]

# TKF92 Parameters
LAMBDA = config["tkf_lambda"]
MU = config["tkf_mu"]
R = config["tkf_r"]
MAX_INS = config["max_insertion_length"]

rule all:
    input:
        [
            f"results/msas/s{s}_b{pair[0]}_d{pair[1]}_f{f}_m{m}_seed{seed}/msa.fasta"
            for s in SPECIES
            for pair in BIRTH_DEATH_PAIRS
            for f in SAMPLING
            for m in MUTATION
            for seed in SEEDS
        ],
        [
            f"results/msas/s{s}_b{pair[0]}_d{pair[1]}_f{f}_m{m}_seed{seed}/tree_plot.png"
            for s in SPECIES
            for pair in BIRTH_DEATH_PAIRS
            for f in SAMPLING
            for m in MUTATION
            for seed in SEEDS
        ]

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
        msa = "results/msas/s{s}_b{b}_d{d}_f{f}_m{m}_seed{seed}/msa.fasta",
        masa = "results/msas/s{s}_b{b}_d{d}_f{f}_m{m}_seed{seed}/masa.fasta",
        info = "results/msas/s{s}_b{b}_d{d}_f{f}_m{m}_seed{seed}/info.txt",
        tree_copy = "results/msas/s{s}_b{b}_d{d}_f{f}_m{m}_seed{seed}/tree.nwk"
    params:
        bin = SIMULATE_TKF,
        lambda_val = LAMBDA,
        mu_val = MU,
        r_val = R,
        max_ins = MAX_INS,
        out_dir = "results/msas/s{s}_b{b}_d{d}_f{f}_m{m}_seed{seed}"
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
        tree = "results/msas/s{s}_b{b}_d{d}_f{f}_m{m}_seed{seed}/tree.nwk"
    output:
        plot = "results/msas/s{s}_b{b}_d{d}_f{f}_m{m}_seed{seed}/tree_plot.png"
    params:
        script = "scripts/visualize_trees.py",
        py_bin = PYTHON
    shell:
        """
        {params.py_bin} {params.script} --tree-file {input.tree} --output-file {output.plot}
        """
