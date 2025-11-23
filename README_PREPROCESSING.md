# Graph Preprocessing Pipeline

This document describes the graph preprocessing pipeline for the Dense-PCE project, which converts `.grh` files through multiple stages to create all necessary graph format files.

## Overview

The preprocessing pipeline converts graph files through the following stages:

```
.grh → .edges → .clean → b_adj.bin + b_degree.bin
```

## Pipeline Stages

### Stage 1: .grh to .edges
- **Input**: Graph in adjacency list format (`.grh`)
- **Output**: Edge list format (`.edges`)
- **Process**: Parses adjacency lists and extracts unique edges

### Stage 2: .edges to .clean
- **Input**: Edge list file (`.edges`)
- **Output**: Cleaned edge list (`.clean`)
- **Tool**: BBkC preprocessor
- **Process**: Removes duplicate edges, self-loops, and normalizes node IDs

### Stage 3: .clean to Binary
- **Input**: Clean edge list (`.clean`)
- **Output**: Binary adjacency (`b_adj.bin`) and degree (`b_degree.bin`) files
- **Tool**: Cohesive_subgraph_book/edgelist2binary
- **Process**: Converts to optimized binary format for fast graph operations

## Usage

### Prerequisites

Ensure the following tools are built and available:
```bash
# Build BBkC
cd EBBkC/src && mkdir -p build && cd build
cmake .. && make
cp BBkC ../../../

# Ensure edgelist2binary is available
cd EBBkC/Cohesive_subgraph_book/datasets/
make  # Build the conversion tools
```

### Single File Processing

```bash
# Process a single graph file
python graph_preprocessing_pipeline.py testGraphs/fpce_graph/fpce_graph.grh

# Process with custom output directory
python graph_preprocessing_pipeline.py graph.grh ./processed_graphs/graph_output/
```

### Batch Processing

```bash
# Process all .grh files in a directory
python graph_preprocessing_pipeline.py --batch synth-graphs-1000/

# Process batch with custom output base directory
python graph_preprocessing_pipeline.py --batch synth-graphs-1000/ ./batch_processed/
```

### Advanced Options

```bash
# Specify custom project root
python graph_preprocessing_pipeline.py --project-root /path/to/Dense-PCE-main/ input.grh

# View help and examples
python graph_preprocessing_pipeline.py --help
```

## Output Structure

For each input file `graph_name.grh`, the pipeline creates:

```
graph_name_processed/
├── graph_name.grh      # Original graph file (copied)
├── graph_name.edges    # Edge list format
├── graph_name.clean    # Cleaned edge list
├── b_adj.bin          # Binary adjacency structure
└── b_degree.bin       # Binary degree information
```

## File Formats

### .grh Format (Input)
```
1 2 3 4 6 7    # Node 0 connected to nodes 1,2,3,4,6,7
0              # Node 1 connected to node 0
0 6 7          # Node 2 connected to nodes 0,6,7
...
```

### .edges Format
```
0 1
0 2
0 3
...
```

### .clean Format
```
0 1
0 2
0 3
...

```
(Note: BBkC adds a trailing empty line)

### Binary Formats
- `b_adj.bin`: Compressed adjacency list representation
- `b_degree.bin`: Node degree information for fast access

## Error Handling

The pipeline includes comprehensive error handling and logging:

- **Tool Verification**: Checks that BBkC and edgelist2binary are available
- **File Validation**: Verifies input files exist and are readable
- **Process Monitoring**: Logs each stage with success/failure status
- **Output Validation**: Confirms expected output files are created

## Logging

All operations are logged to both console and `graph_preprocessing.log`:

```
2024-01-XX 10:30:15 - INFO - Converting graph.grh to edges format
2024-01-XX 10:30:15 - INFO - Successfully converted to graph.edges
2024-01-XX 10:30:16 - INFO - Running command: ./BBkC p graph.edges
2024-01-XX 10:30:17 - INFO - Successfully created clean file: graph.clean
2024-01-XX 10:30:18 - INFO - Successfully created binary files: b_adj.bin and b_degree.bin
```

## Integration with Algorithms

The generated files can be used with different algorithms:

### Dense-PCE (Original)
```bash
./dense-pce graph_name_processed/graph_name.grh --theta 0.8 --minimum 3
```

### Dense-PCE Modified (with Turan Filtering)
```bash
./dense-pce-mod graph_name_processed/graph_name.grh --theta 0.8 --minimum 3
```

### EBBkC (Direct Binary)
```bash
./BBkC e graph_name_processed/ 4 2
```

## Troubleshooting

### Common Issues

1. **BBkC not found**
   - Ensure BBkC is built: `cd EBBkC/src/build && make`
   - Copy to project root: `cp BBkC ../../../`

2. **edgelist2binary not found**
   - Build tools: `cd EBBkC/Cohesive_subgraph_book/datasets && make`

3. **Permission errors**
   - Make tools executable: `chmod +x BBkC edgelist2binary`

4. **Invalid .grh format**
   - Check for non-numeric values or malformed lines
   - Ensure node IDs are 0-indexed integers

### Debug Mode

Enable verbose logging by modifying the script:
```python
logging.basicConfig(level=logging.DEBUG)
```

## Performance Notes

- **Large graphs**: Processing time scales with graph size and density
- **Memory usage**: Binary conversion may require significant memory for large graphs
- **Disk space**: Output files may be larger than input for sparse graphs
- **Parallel processing**: Consider batch processing large graph collections

## Examples

### Example 1: Small Test Graph
```bash
python graph_preprocessing_pipeline.py testGraphs/fpce_graph/fpce_graph.grh
```

Output:
- `fpce_graph_processed/fpce_graph.grh` (13 lines)
- `fpce_graph_processed/fpce_graph.edges` (14 edges)
- `fpce_graph_processed/fpce_graph.clean` (14 edges + empty line)
- `fpce_graph_processed/b_adj.bin` (binary adjacency)
- `fpce_graph_processed/b_degree.bin` (binary degrees)

### Example 2: Synthetic Graph Collection
```bash
python graph_preprocessing_pipeline.py --batch synth-graphs-1000/
```

Processes all scale-free graphs in the collection, creating individual processed directories for each.
