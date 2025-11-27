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

> **Note:** If your graph is in edge list format (e.g., `u v` pairs), use the provided `transgrh.pl` script to convert it to `.grh` first:
> ```bash
> perl transgrh.pl < input.edges > output.grh
> ```

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
./build_integrated/dense-pce-r1-m1-integrated <graph_file_path_> --minimum <min_size> --theta <density> [options]
```

- `<graph_file_path>`: Path to the `.grh` file (system automatically looks for binaries in the same folder).
- `--minimum <l>`: Minimum size of pseudo-cliques to find (e.g., 10).
- `--theta <θ>`: Minimum density threshold (0.0 to 1.0, e.g., 0.9).

**Example**:
```bash
./build_integrated/dense-pce-r1-m1-integrated testGraph/fpce_graph/fpce_graph.grh --minimum 3 --theta 0.6
```

## ⚙️ Ablation Studies & Advanced Modes

The system supports a "Gated" operation mode for ablation studies, allowing you to toggle specific optimizations (Order Bound, Edge Bound, Turán Seeding) to analyze their impact.

### Build for Ablation
While `dense-pce-r1-m1-integrated` supports these flags, you can also explicitly build the gated variant:
```bash
bash build_integrated.sh dense-pce-r1-gated.cpp
# Creates: ./build_integrated/dense-pce-r1-gated-integrated
```

### Operating Modes
Use the `--mode <N>` switch to select a preset configuration:

| Mode | Name | Description |
|------|------|-------------|
| **Mode 1** | **Baseline** | No pruning, no seeding. Pure node-by-node enumeration. |
| **Mode 2** | **Order Only** | Enables Order-Bound pruning. No Edge-Bound, no Turán. |
| **Mode 3** | **Order + Edge** | Enables Order and Edge-Bound pruning. No Turán seeding. |
| **Mode 4** | **Full FPCE** | (Default) Order + Edge + Turán seeding. |

### Fine-Grained Control
You can override specific components of a mode using flags:
- **Order Bound**: `--order` / `--no-order`
- **Edge Bound**: `--edge` / `--no-edge`
- **Turán Seeding**: `--turan` / `--no-turan`

**Examples**:
```bash
# Mode 2 (Order only)
./build_integrated/dense-pce-r1-gated-integrated testGraph/fpce_graph/fpce_graph.grh --minimum 3 --theta 0.6 --mode 2


### Batch Processing with Helper Script

For systematic experiments, use the `run.sh` script. It runs one graph per subdirectory, supports parallel execution, and records per‑graph logs.

**Usage**:
```bash
bash run.sh -d <graph_root_dir> -l <min_size> -t <theta> -m <mode> [-j N] [-x <exec_path>]
```

**Parameters**:
- `-d DIR`: Graph root directory (scans subdirectories for `.grh` files).
- `-l SIZE`: Minimum size $\ell$ (default 10).
- `-t THETA`: Density threshold $\theta$ (default 0.9).
- `-m MODE`: Ablation mode (1–4, default 4).
- `-j N`: Number of concurrent jobs (default 1).
- `-x PATH`: Executable path (default: `./build_integrated/dense-pce-r1-gated-integrated`).

**Examples**:
```bash
# Run Mode 4 (Full FPCE) on all graphs in testGraph/
bash run.sh -d testGraph -l 3 -t 0.6 -m 4

# Run Mode 2 (Order-only) with 4 parallel jobs
bash run.sh -d testGraph -l 5 -t 0.8 -m 2 -j 4
```

**Logs**:
Results are saved in `logs/run_<timestamp>/`, with one log file per graph.

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
