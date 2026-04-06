#!/usr/bin/env python3
import sys
import re
from datetime import datetime, timedelta

def parse_jati_time(log_content):
    lines = log_content.strip().splitlines()
    if not lines:
        return None
    
    first_line = lines[0]
    last_line = lines[-1]
    
    fmt = "%H:%M:%S"
    try:
        start_time_str = first_line.split(' ')[0]
        end_time_str = last_line.split(' ')[0]
        
        start_dt = datetime.strptime(start_time_str, fmt)
        end_dt = datetime.strptime(end_time_str, fmt)
        
        if end_dt < start_dt:
            end_dt += timedelta(days=1)
            
        return (end_dt - start_dt).total_seconds()
    except (ValueError, IndexError):
        return None

def parse_iqtree_time(log_content):
    # Pattern: "CPU time used for tree search: 0.709 sec (0h:0m:0s)"
    match = re.search(r"CPU time used for tree search:\s+([\d.]+)\s+sec", log_content)
    if match:
        return float(match.group(1))
    
    # Fallback for total Wall-clock time if search time not found
    match = re.search(r"Total Wall-clock time:\s+([\d.]+)\s+sec", log_content)
    if match:
        return float(match.group(1))
    
    return None

def main():
    if len(sys.argv) != 3:
        print("Usage: calculate_time.py <input_log> <output_time>")
        sys.exit(1)

    input_log = sys.argv[1]
    output_time = sys.argv[2]

    with open(input_log, 'r') as f:
        content = f.read()

    # Try IQ-TREE first as it has a very specific signature
    duration = parse_iqtree_time(content)
    
    # If not IQ-TREE, try JATI (timestamp-based)
    if duration is None:
        duration = parse_jati_time(content)

    if duration is not None:
        with open(output_time, 'w') as out:
            out.write(f"{duration}")
    else:
        # If both fail, we might want to know why or just write NA
        with open(output_time, 'w') as out:
            out.write("NA")

if __name__ == "__main__":
    main()
