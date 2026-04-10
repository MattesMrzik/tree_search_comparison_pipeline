import re
from datetime import datetime, timedelta

def parse_jati_time(log_content):
    lines = log_content.strip().splitlines()
    if not lines:
        return None
    
    # TODO: it could speed up if we reed in the last bits that probably contain the end time
    first_line = lines[0]
    last_line = lines[-1]
    
    fmt = "%H:%M:%S"
    try:
        start_time_str = first_line.split(' ')[0]
        end_time_str = last_line.split(' ')[0]
        
        start_dt = datetime.strptime(start_time_str, fmt)
        end_dt = datetime.strptime(end_time_str, fmt)
        
        # Handle case where end time is past midnight
        if end_dt < start_dt:
            end_dt += timedelta(days=1)
            
        return (end_dt - start_dt).total_seconds()
    except (ValueError, IndexError):
        return None

def parse_iqtree_time(log_path):
    with open(log_path, 'r') as f:
        for line in f:
            if "CPU time used for tree search:" in line:
                match = re.search(r"CPU time used for tree search:\s+([\d.]+)\s+sec", line)
                if match:
                    return float(match.group(1))
    return None

