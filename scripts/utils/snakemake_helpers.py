import os

def get_msa_length(wildcards, msa_path_template):
    """
    Calculates alignment length by reading the second line of the FASTA.
    
    Args:
        wildcards: The Snakemake wildcards object for the current job.
        msa_path_template: The path template (e.g. results/msas/...) containing wildcards.
        
    Returns:
        int: Length of the first sequence in the MSA, or 0 if file not found.
    """
    # Fill the template with fixed values from the current job's wildcards
    # Since wildcards contains fixed values for a specific job, .format() returns a single string
    msa_path = msa_path_template.format(**wildcards)
    
    if os.path.exists(msa_path):
        try:
            with open(msa_path, 'r') as f:
                f.readline()  # Skip header
                seq = f.readline().strip()
                return len(seq)
        except (IOError, IndexError):
            return 0
    return 0
