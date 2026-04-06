import os
import re

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
    """Calculates comprehensive gap statistics for an MSA."""
    if not os.path.exists(msa_path):
        return {
            "gap%": "NA", 
            "gap_col%": "NA", 
            "avg_gap_len": "NA"
        }
    try:
        sequences = []
        with open(msa_path, 'r') as f:
            current_seq = []
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if current_seq:
                        sequences.append("".join(current_seq))
                    current_seq = []
                else:
                    current_seq.append(line)
            if current_seq:
                sequences.append("".join(current_seq))

        if not sequences:
            return {"gap%": "NA", "gap_col%": "NA", "avg_gap_len": "NA"}

        msa_len = len(sequences[0])
        num_seqs = len(sequences)
        total_chars = msa_len * num_seqs
        
        total_gaps = 0
        gap_columns = 0
        gap_events_count = 0
        
        # Column-wise check for gaps
        for j in range(msa_len):
            has_gap = False
            for i in range(num_seqs):
                if sequences[i][j] == '-':
                    has_gap = True
                    total_gaps += 1
            if has_gap:
                gap_columns += 1
        
        # Row-wise check for gap events (to calculate average gap length)
        for seq in sequences:
            # Find all continuous blocks of '-'
            gap_events = re.findall(r'-+', seq)
            gap_events_count += len(gap_events)

        gap_pct = (total_gaps / total_chars * 100) if total_chars > 0 else 0
        gap_col_pct = (gap_columns / msa_len * 100) if msa_len > 0 else 0
        avg_gap_len = (total_gaps / gap_events_count) if gap_events_count > 0 else 0

        return {
            "gap%": round(gap_pct, 2),
            "gap_col%": round(gap_col_pct, 2),
            "avg_gap_len": round(avg_gap_len, 2)
        }
    except Exception:
        pass
    return {"gap%": "NA", "gap_col%": "NA", "avg_gap_len": "NA"}
