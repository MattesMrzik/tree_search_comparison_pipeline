import os

def get_fasta_length(msa_path):
    """Returns the length of the first sequence in the FASTA."""
    if os.path.exists(msa_path):
        try:
            with open(msa_path, 'r') as f:
                seq = ""
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    if line.startswith(">"):
                        if seq:
                            return len(seq)
                    else:
                        seq += line
                return len(seq)
        except Exception:
            pass
    return "NA"

def get_gap_stats(msa_path):
    """Calculates gap count and percentage by scanning the MSA file."""
    if not os.path.exists(msa_path):
        return {"gap_count": "NA", "gap_percentage": "NA"}
    try:
        total_chars = 0
        gaps = 0
        with open(msa_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith(">"):
                    continue
                total_chars += len(line)
                gaps += line.count('-')
        
        pct = (gaps / total_chars * 100) if total_chars > 0 else 0
        return {"gap_count": gaps, "gap_percentage": round(pct, 2)}
    except Exception:
        pass
    return {"gap_count": "NA", "gap_percentage": "NA"}
