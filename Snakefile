configfile: "config.yaml"

# Path to evolver
EVOLVER = config["evolver_path"]

# Extracting parameters for expansion
SPECIES = config["species"]
BIRTH = config[ "birth_death_rates"][0]
DEATH = config["birth_death_rates"][1]
REPS = config["replicates"]
SAMPLING = config["sampling_fraction"]
MUTATION = config["mutation_rate"]

rule all:
    input:
        expand("results/tree_s{s}_b{b}_d{d}_f{f}_m{m}_r{r}.tre",
               s=SPECIES,
               b=BIRTH,
               d=DEATH,
               f=SAMPLING,
               m=MUTATION,
               r=REPS)

rule generate_tree:
    output:
        "results/tree_s{s}_b{b}_d{d}_f{f}_m{m}_r{r}.tre"
    params:
        evolver = EVOLVER
    shell:
        """
        # Option 2: Rooted trees
        # Next line: number of species
        # Next line: 1 tree, random seed (replicate), want branch lengths (1)
        # Next line: birth, death, sampling, mutation
        printf "2\\n{wildcards.s}\\n1 {wildcards.r} 1\\n{wildcards.b} {wildcards.d} {wildcards.f} {wildcards.m}\\n0\\n" | {params.evolver}
        
        # Extract the tree from the end of evolver.out
        tail -n 1 evolver.out > {output}
        rm evolver.out
        """
