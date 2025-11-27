# Graph Preprocessing Pipeline

This document describes the graph preprocessing pipeline for the Dense-PCE project. This pipeline converts `.grh` files into the efficient binary formats required by the Dense-PCE system.

## Overview

The preprocessing pipeline transforms graph files through the following stages:

```
.grh → .edges → .clean → b_adj.bin + b_degree.bin
```

1.  **.grh to .edges**: Extracts edges from the adjacency list.
2.  **.edges to .clean**: Uses `BBkC` to remove duplicates, self-loops, and normalize node IDs.
3.  **.clean to Binary**: Uses `edgelist2binary` to create optimized binary adjacency (`b_adj.bin`) and degree (`b_degree.bin`) files.

## Usage

### Prerequisites

The `graph_pipeline.py` script automatically handles the checking and building of necessary tools (`BBkC` and `edgelist2binary`). You just need Python 3 installed.

### Single File Processing

To process a specific graph file:

```bash
python graph_pipeline.py testGraphs/fpce_graph/fpce_graph.grh
```

This will generate the following files in the same directory as the input:
- `fpce_graph.edges`
- `fpce_graph.clean`
- `b_adj.bin`
- `b_degree.bin`

### Custom Output Directory

You can specify a different directory for the output files:

```bash
python graph_pipeline.py testGraphs/fpce_graph/fpce_graph.grh ./processed_graphs/
```

### Batch Processing

To process all `.grh` files in a directory:

```bash
python graph_pipeline.py --batch synth-graphs-1000/
```

## File Formats

### .grh Format (Input)
Adjacency list format. Each line starts with a node ID (implied or explicit) followed by its neighbors.
```text
1 2 3    # Node 0 connected to 1, 2, 3
0 5      # Node 1 connected to 0, 5
...
```

### Handling Edge List Files
If your source data is in edge list format (pairs of connected nodes), you must convert it to `.grh` format before running the pipeline. Use the `transgrh.pl` script:

```bash
# Standard conversion (undirected)
perl transgrh.pl < input.txt > output.grh

# For bipartite graphs (adds offset to second column)
perl transgrh.pl B < input.txt > output.grh
```

### Binary Formats (Output)
- **`b_adj.bin`**: Compressed adjacency list representation.
- **`b_degree.bin`**: Node degree information for fast access.

These binary files allow Dense-PCE to load large graphs (millions of nodes) in seconds.

## Error Handling

The pipeline logs its progress to `graph_preprocessing.log`. If a step fails (e.g., due to missing tools or permissions), check this log file for details.

## Troubleshooting

1.  **Tools not building**: Ensure you have `cmake` and a C++ compiler (`g++` or `clang++`) installed.
2.  **WSL Issues**: If running on Windows, the script attempts to use WSL. Ensure WSL is installed and functional.
