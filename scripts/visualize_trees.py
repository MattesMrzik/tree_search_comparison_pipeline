import toytree
import toyplot
import toyplot.png
import argparse

def main():
    parser = argparse.ArgumentParser(description="Visualize a Newick tree.")
    parser.add_argument("--tree-file", required=True, help="Path to a Newick tree file")
    parser.add_argument("--output-file", required=True, help="Path to save the output PNG")
    args = parser.parse_args()

    try:
        tree = toytree.tree(args.tree_file)
        
        # Correctly import and use toyplot.Canvas
        canvas = toyplot.Canvas(width=800, height=800, style={"background-color": "white"})
        axes = canvas.cartesian(show=False)
        
        tree.draw(
            axes=axes,
            node_labels=False,
            node_sizes=5,
            tip_labels_colors="black"
        )
        
        # Save with a high scale factor for even better print-ready resolution
        toyplot.png.render(canvas, args.output_file, scale=4.0)
        print(f"Generated visualization: {args.output_file}")
    except Exception as e:
        print(f"Error drawing tree {args.tree_file}: {e}")

if __name__ == "__main__":
    main()
