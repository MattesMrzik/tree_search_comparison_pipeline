#!/usr/bin/env python3
import sys
from datetime import datetime

def main():
    if len(sys.argv) != 3:
        print("Usage: calculate_time.py <input_log> <output_time>")
        sys.exit(1)

    input_log = sys.argv[1]
    output_time = sys.argv[2]

    with open(input_log, 'rb') as f:
        # Read the first line
        first_line = f.readline().decode().strip()
        if not first_line:
            raise ValueError(f"Log file {input_log} is empty")
        
        # Seek to near the end to find the last line
        f.seek(0, 2)
        pos = f.tell()
        
        # Read in chunks backwards from the end
        # A JATI log line is usually < 200 chars, so 512 is plenty
        chunk_size = 512
        if pos > chunk_size:
            f.seek(pos - chunk_size)
            buffer = f.read(chunk_size)
        else:
            f.seek(0)
            buffer = f.read(pos)
        
        lines = buffer.decode().strip().splitlines()
        last_line = lines[-1] if lines else first_line

    fmt = "%H:%M:%S"
    try:
        start_time_str = first_line.split(' ')[0]
        end_time_str = last_line.split(' ')[0]
        
        start_dt = datetime.strptime(start_time_str, fmt)
        end_dt = datetime.strptime(end_time_str, fmt)
        
        # If the end time is earlier than the start time, it likely crossed midnight
        if end_dt < start_dt:
            from datetime import timedelta
            end_dt += timedelta(days=1)
            
        duration = end_dt - start_dt
        
        with open(output_time, 'w') as out:
            out.write(f"{duration.total_seconds()}")
    except (ValueError, IndexError) as e:
        raise ValueError(f"Could not parse timestamps from log lines:\nFirst: {first_line}\nLast: {last_line}") from e

if __name__ == "__main__":
    main()
