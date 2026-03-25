configfile: "config.yaml"

import os

# Path to evolver
EVOLVER = config["evolver_path"]

# Extracting parameters
SPECIES = config["species"]
BIRTH_DEATH_PAIRS = config["birth_death_rates"]
REPS = config["replicates"]
SAMPLING = config["sampling_fraction"]
MUTATION = config["mutation_rate"]

rule all:
    input:
        [
            f"results/trees/s{s}_b{pair[0]}_d{pair[1]}_f{f}_m{m}_r{r}.nwk"
            for s in SPECIES
            for pair in BIRTH_DEATH_PAIRS
            for f in SAMPLING
            for m in MUTATION
            for r in REPS
        ]

rule generate_tree:
    output:
        "results/trees/s{s}_b{b}_d{d}_f{f}_m{m}_r{r}.nwk"
    params:
        evolver = EVOLVER
    shadow: "minimal"
    shell:
        """
        # Run evolver. Because of shadow: "minimal", evolver.out 
        # is created in an isolated directory for this specific job.
        printf "2\\n{wildcards.s}\\n1 {wildcards.r} 1\\n{wildcards.b} {wildcards.d} {wildcards.f} {wildcards.m}\\n0\\n" | {params.evolver} > /dev/null 2>&1
        
        # Extract the tree from the local, isolated evolver.out
        tail -n 1 evolver.out > {output}
        """
