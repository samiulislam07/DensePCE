# Dense-PCE: Efficient Pseudo-Clique Enumeration

Dense-PCE is a high-performance C++ system for finding maximal $(l, \theta)$-pseudo-cliques in large graphs. It integrates **EBBkC** (Edge-Based Branch and Bound k-Clique) as an in-memory library for fast seed generation using Turán's theorem, combined with advanced pruning techniques (order-bound and edge-bound).

## 📂 Directory Overview

Here is a quick overview of the repository structure:

- **`build_integrated.sh`**: The main build script. Compiles EBBkC as a static library and links it with Dense-PCE.
- **`graph_pipeline.py`**: A Python script for preprocessing `.grh` graph files into the formats required by the system (`.edges`, `.clean`, binary files).
- **`dense-pce-r1-m1.cpp`**: The primary source code for the Dense-PCE algorithm with integrated EBBkC support.
- **`dense-pce-r1-gated.cpp`**: A variant of the source code used for ablation studies (supports toggling specific optimizations).
- **`EBBkC/`**: Contains the EBBkC library source code, used here for efficient clique seeding.
- **`testGraph/`**: Directory containing example graph datasets.
- **`build_integrated/`**: (Created after build) Contains the compiled executable and libraries.

## 🚀 Quick Start

### Prerequisites

- **OS**: Linux or Windows (via WSL2).
- **Compiler**: GCC 7+ or Clang 5+ (must support C++17).
- **Tools**: CMake (3.6+), Python 3.x.
- **Libraries**: OpenMP, Intel TBB (Threading Building Blocks).

### 1. Build the System

The easiest way to build is using the provided script, which handles dependencies and linking automatically.

```bash
# For Linux / macOS / WSL
bash build_integrated.sh
```

This will create the executable `dense-pce-r1-m1-integrated` inside the `build_integrated/` directory.

### 2. Preprocess Your Graph

Dense-PCE requires graph data in a specific binary format for maximum efficiency. Use the included pipeline script to convert your raw graph files.

**Input Format (`.grh`)**:
Adjacency list format where each line represents a node and its neighbors.
```text
1 2 3    # Node 0 connects to 1, 2, 3
0 5      # Node 1 connects to 0, 5
...
```

**Run the Pipeline**:
```bash
# Syntax: python graph_pipeline.py <path_to_graph_file>
python graph_pipeline.py testGraph/fpce_graph/fpce_graph.grh
```

This process:
1.  Converts `.grh` to `.edges` (edge list).
2.  Cleans the graph (removes self-loops, duplicates) -> `.clean`.
3.  Generates binary files (`b_adj.bin`, `b_degree.bin`) in the same directory.

**Output**:
You will see a set of files generated in the graph's directory. The system needs the `.grh` file AND the binary files (`b_adj.bin`, `b_degree.bin`) to be present in the same folder for the fastest loading.

### 3. Run Dense-PCE

Run the compiled executable with your graph and parameters.

**Syntax**:
```bash
./build_integrated/dense-pce-r1-m1-integrated <graph_file> --minimum <min_size> --theta <density> [options]
```

- `<graph_file>`: Path to the `.grh` file (system automatically looks for binaries in the same folder).
- `--minimum <l>`: Minimum size of pseudo-cliques to find (e.g., 10).
- `--theta <θ>`: Minimum density threshold (0.0 to 1.0, e.g., 0.9).

**Example**:
```bash
./build_integrated/dense-pce-r1-m1-integrated testGraph/fpce_graph/fpce_graph.grh --minimum 3 --theta 0.6
```

## ⚙️ Modes & Advanced Usage

Dense-PCE supports different operating modes to toggle optimizations (Order Bound, Edge Bound, Turán Seeding). This is useful for performance analysis (ablation studies).

Use the `--mode` flag or individual switches:

| Mode | Description |
|------|-------------|
| **Mode 4** (Default) | Full optimizations: Order Bound + Edge Bound + Turán Seeding. |
| **Mode 3** | Order Bound + Edge Bound (No Turán). |
| **Mode 2** | Order Bound only. |
| **Mode 1** | No pruning (Baseline enumeration). |

**Example**:
```bash
# Run in Mode 2 (Order bound only)
./build_integrated/dense-pce-r1-m1-integrated testGraph/fpce_graph/fpce_graph.grh --minimum 3 --theta 0.6 --mode 2
```

You can also fine-tune specific components:
- `--no-turan`: Disable Turán seeding.
- `--no-order`: Disable order-bound pruning.
- `--no-edge`: Disable edge-bound pruning.

## 📊 Output Interpretation

The program outputs the count of maximal pseudo-cliques found for each size, along with performance metrics.

```text
Maximal pseudo-clique counts:
Size 3: 7
Size 4: 2
...
Total Maximal Pseudo-Cliques: 9

Total Iterations: 145
Edge-bound prunes saved: 42
```

- **Total Iterations**: Number of recursive calls explored.
- **Edge-bound prunes saved**: How many times the edge-bound optimization successfully cut off a search branch.

## 🛠️ Troubleshooting

- **"EBBkC binaries not found"**: Ensure you ran `graph_pipeline.py` and that `b_adj.bin` and `b_degree.bin` are in the same folder as your graph file.
- **Build Fails (OpenMP/TBB)**: Make sure you have the development libraries installed.
  - Ubuntu: `sudo apt install libomp-dev libtbb-dev`
- **WSL Paths**: If running on Windows via WSL, ensure you are executing commands within the WSL terminal and accessing files via the Linux file system for best performance.
