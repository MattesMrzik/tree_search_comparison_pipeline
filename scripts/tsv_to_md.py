import pandas as pd
import sys
import os

def main():
    if len(sys.argv) != 3:
        print("Usage: tsv_to_md.py <input_tsv> <output_md>")
        sys.exit(1)

    input_path = sys.argv[1]
    output_path = sys.argv[2]

    if not os.path.exists(input_path):
        print(f"Error: {input_path} not found")
        sys.exit(1)

    try:
        df = pd.read_csv(input_path, sep='\t')
        # Replace NaN with 'NA' for cleaner markdown
        df = df.fillna('NA')
        
        # Convert to markdown
        md_table = df.to_markdown(index=False)
        
        # Ensure directory exists
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        
        with open(output_path, 'w') as f:
            f.write("---\n")
            f.write("cssclasses:\n")
            f.write("  - page-100\n")
            f.write("---\n\n")
            f.write("# Pipeline Execution Summary\n\n")
            f.write(md_table)
            
    except Exception as e:
        print(f"Error converting TSV to Markdown: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
