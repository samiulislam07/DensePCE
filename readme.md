# Dense-PCE + EBBkC Integrated System

A high-performance system for finding dense pseudo-cliques in graphs using the Dense-PCE algorithm with EBBkC (Edge-Based Branch and Bound k-Clique) as an in-memory library for seed generation.

## 🚀 Quick Start

### Prerequisites

- **Linux/WSL**: The system is designed for Linux environments. On Windows, use WSL (Windows Subsystem for Linux).
- **C++17 Compiler**: GCC 7+ or Clang 5+
- **CMake**: Version 3.6 or higher
- **OpenMP**: For parallel processing
- **Intel TBB**: Threading Building Blocks library

### Build the Integrated System

```bash
# Clone and navigate to the repository
cd Dense_PCE/Dense-PCE-main

# Build the integrated system (EBBkC library + Dense-PCE)
# Choose one of the following methods:

# Method 1: Bash script (Linux/macOS/WSL) - Recommended
bash build_integrated.sh

# Method 2: PowerShell script (Windows)
.\build_integrated.ps1

# Method 3: Batch script (Windows Command Prompt)
build_integrated.bat

# Method 4: WSL (if on Windows)
wsl bash build_integrated.sh
```

### Run the System

```bash
# Basic usage
./build_integrated/dense-pce-mod-edge-order-integrated <graph-file> --minimum <min-size> --theta <density>

# Example: Find pseudo-cliques of size ≥10 with density ≥0.9
./build_integrated/dense-pce-mod-edge-order-integrated testGraphs/gplus/gplus.grh --minimum 10 --theta 0.9
```

## 📁 System Architecture

### Components

1. **EBBkC Library** (`EBBkC/src/`):
   - **Purpose**: Generates R-clique seeds for pseudo-clique enumeration
   - **API**: `EBBkC_t::list_k_clique_mem()` - returns cliques in memory
   - **Build**: Creates `libebbkc_core.a` static library

2. **Dense-PCE** (`dense-pce-mod-edge-order.cpp`):
   - **Purpose**: Enumerates maximal pseudo-cliques using Turán's theorem
   - **Integration**: Directly calls EBBkC library (no file I/O)
   - **Features**: Edge-bound pruning, core-based optimization

3. **Graph Pipeline** (`graph_pipeline.py`):
   - **Purpose**: Preprocesses graphs for the system
   - **Output**: Generates `.edges`, `.clean`, `b_adj.bin`, `b_degree.bin` files

### Data Flow

```
Graph (.grh) → Graph Pipeline → Binary Files → EBBkC Library → R-cliques → Dense-PCE → Maximal Pseudo-cliques
```

## 🔧 Detailed Build Instructions

### Option 1: Automated Build (Recommended)

```bash
# Run the automated build script
bash build_integrated.sh
```

This script:
1. **Cleans up old build artifacts** to avoid path conflicts
2. **Builds EBBkC as a static library** (`libebbkc_core.a`)
3. **Links all dependencies** (`libcommon-utils.a`, `libgraph-pre-processing.a`)
4. **Compiles the integrated executable** with proper flags
5. **Handles location changes** gracefully (no more CMake cache issues)

**Cross-platform support:**
- **Linux/macOS**: `bash build_integrated.sh`
- **Windows PowerShell**: `.\build_integrated.ps1`
- **Windows Command Prompt**: `build_integrated.bat`

### Option 2: Manual Build

```bash
# Step 1: Build EBBkC library
cd EBBkC/src
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)

# Step 2: Compile integrated system
cd ../../../
g++ -std=c++17 -O3 -march=native -fopenmp \
    -I./EBBkC/src \
    -I./EBBkC/src/truss/dependencies/sparsepp \
    -I./EBBkC/src/truss/dependencies/libpopcnt \
    -I./EBBkC/src/truss/dependencies \
    -I./EBBkC/src/truss/util \
    -I./EBBkC/src/truss/decompose \
    dense-pce-mod-edge-order.cpp \
    EBBkC/src/build/libebbkc_core.a \
    EBBkC/src/build/libcommon-utils.a \
    EBBkC/src/build/libgraph-pre-processing.a \
    -lgomp -ltbb -lpthread \
    -o dense-pce-mod-edge-order-integrated
```

## 📊 Usage Examples

### Basic Usage

```bash
# Find pseudo-cliques with minimum size 10 and density threshold 0.9
./build_integrated/dense-pce-mod-edge-order-integrated testGraphs/gplus/gplus.grh --minimum 10 --theta 0.9
```

### Advanced Usage

```bash
# High-density pseudo-cliques (θ = 0.95)
./build_integrated/dense-pce-mod-edge-order-integrated testGraphs/tech/tech.grh --minimum 15 --theta 0.95

# Large pseudo-cliques with moderate density
./build_integrated/dense-pce-mod-edge-order-integrated testGraphs/email-EU/email-EU.grh --minimum 20 --theta 0.8
```

### Batch Processing

```bash
# Process multiple graphs
for graph in testGraphs/*/*.grh; do
    echo "Processing $graph..."
    ./build_integrated/dense-pce-mod-edge-order-integrated "$graph" --minimum 10 --theta 0.9
done
```

## 📈 Performance Features

### Optimizations

1. **In-Memory Integration**: No file I/O between EBBkC and Dense-PCE
2. **Edge-Bound Pruning**: Eliminates unnecessary recursive calls
3. **Core-Based Filtering**: Uses graph degeneracy for early termination
4. **Parallel Processing**: OpenMP support for multi-threading

### Performance Metrics

The system reports:
- **Pseudo-clique counts** by size
- **Total iterations** performed
- **Edge-bound pruning efficiency** (calls saved)
- **Build times** for main and sub-problems

### Example Output

```
Computed R: 6
Directory path: testGraphs/gplus
Detected BBkC binaries in: testGraphs/gplus, loading graph from binaries...
degeneracy: 12
NOT Pruning by Order Bound: (l,θ)-pseudo-cliques exist as l (10) <= μ (26)
|V| = 725, |E| = 3956
Truss number = 6
total nodes: 23628
Maximal pseudo-clique counts:
Size 10: 8

Total Iterations: 28055
#calls of PCE_iter saved by EDGE bound = 22823
```

## 🔍 Graph Preprocessing

### Using the Graph Pipeline

The graph pipeline has been significantly improved to handle path issues when the project is moved to different devices or locations.

```bash
# Preprocess a graph for the integrated system
python graph_pipeline.py testGraphs/gplus/gplus.grh

# Process with custom output directory
python graph_pipeline.py testGraphs/gplus/gplus.grh ./processed_graphs/

# Process all graphs in a directory
python graph_pipeline.py --batch synth-graphs-1000/

# Process batch with custom output base directory
python graph_pipeline.py --batch synth-graphs-1000/ ./batch_processed/
```

**This creates:**
- `gplus.edges` (edge list)
- `gplus.clean` (cleaned edges)
- `b_adj.bin`, `b_degree.bin` (binary format)

### Testing Path Handling

Before running the main pipeline, test the path handling:

```bash
python test_path_handling.py
```

This verifies:
- WSL availability and functionality
- Path conversion accuracy
- Tool accessibility
- Current directory setup

### Manual Preprocessing

```bash
# Convert .grh to .edges
python3 -c "
import sys
sys.path.append('.')
from graph_pipeline import GraphPipeline
pipeline = GraphPipeline()
pipeline.grhtoedges('testGraphs/gplus/gplus.grh')
"

# Clean edges using BBkC
./EBBkC/src/build/BBkC p testGraphs/gplus/gplus.edges

# Convert to binary format
./EBBkC/Cohesive_subgraph_book/datasets/edgelist2binary testGraphs/gplus gplus.clean
```

## 🛠️ Troubleshooting

### Common Issues

1. **CMake Cache Path Conflicts**:
   ```bash
   # Error: CMakeCache.txt directory is different than where it was created
   # Solution: The build scripts now automatically clean up old cache files
   bash build_integrated.sh  # This will clean and rebuild
   ```

2. **WSL Path Conversion Errors**:
   ```bash
   # Error: WSL command failed with return code 1
   # Solution: Test path handling first
   python test_path_handling.py
   ```

3. **OpenMP Library Not Found**:
   ```bash
   # Install OpenMP development package
   sudo apt-get install libomp-dev  # Ubuntu/Debian
   ```

4. **TBB Library Not Found**:
   ```bash
   # Install Intel TBB
   sudo apt-get install libtbb-dev
   ```

5. **CMake Version Too Old**:
   ```bash
   # Update CMake
   sudo apt-get install cmake
   ```

6. **WSL Not Available**:
   - Install WSL2 on Windows
   - Or use a Linux virtual machine

### Build Verification

```bash
# Check if build was successful
ls -la build_integrated/
# Should show: dense-pce-mod-edge-order-integrated and libebbkc_core.a

# Test the executable
./build_integrated/dense-pce-mod-edge-order-integrated --help

# Test path handling (if using graph pipeline)
python test_path_handling.py
```

## 📚 Algorithm Details

### Dense-PCE Algorithm

- **Input**: Graph G, minimum size l, density threshold θ
- **Output**: All maximal (l,θ)-pseudo-cliques
- **Method**: Turán's theorem with edge-oriented branching

### EBBkC Integration

- **Purpose**: Generate R-clique seeds efficiently
- **Method**: Edge-based branch and bound with truss decomposition
- **Integration**: In-memory API eliminates process overhead

### Performance Characteristics

- **Time Complexity**: O(3^(n/3)) in worst case
- **Space Complexity**: O(n²) for adjacency representation
- **Practical Performance**: Significantly faster than file-based approaches

## 🔧 Recent Improvements

### Build System Enhancements

- **Cross-platform support**: Bash, PowerShell, and Batch scripts
- **Automatic cleanup**: Removes old CMake cache files to prevent path conflicts
- **Location independence**: Works regardless of where the project is moved
- **Better error handling**: Comprehensive error reporting and recovery

### Graph Pipeline Improvements

- **Robust WSL path handling**: Properly handles paths with spaces, parentheses, and special characters
- **Enhanced error recovery**: Timeout protection and detailed logging
- **Working directory management**: Commands execute in correct directory context
- **Batch processing support**: Process multiple graphs efficiently

## 🤝 Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Test thoroughly
5. Submit a pull request

## 📄 License

This project is licensed under the MIT License - see the LICENSE file for details.

## 📞 Support

For issues and questions:
1. Check the troubleshooting section
2. Run the test scripts: `python test_path_handling.py`
3. Review the build logs and `graph_preprocessing.log`
4. Open an issue on GitHub

## 📖 Additional Documentation

- **Build Instructions**: See `BUILD_INSTRUCTIONS.md` for detailed build procedures
- **Graph Pipeline Guide**: See `GRAPH_PIPELINE_README.md` for comprehensive pipeline documentation
- **Path Handling Tests**: Use `test_path_handling.py` to verify system setup

---

**Note**: This integrated system provides significant performance improvements over the original file-based approach by eliminating I/O overhead and process spawning costs. The recent improvements ensure robust operation across different platforms and locations.