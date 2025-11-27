#!/bin/bash

# run.sh
# Helper script to run Dense-PCE batch experiments

# Defaults
MIN_SIZE=10
THETA=0.9
MODE=4
JOBS=1
EXEC_PATH="./build_integrated/dense-pce-r1-m1-integrated"
GRAPH_DIR=""

# Parse arguments
while getopts "d:l:t:m:j:x:h" opt; do
  case $opt in
    d) GRAPH_DIR="$OPTARG" ;;
    l) MIN_SIZE="$OPTARG" ;;
    t) THETA="$OPTARG" ;;
    m) MODE="$OPTARG" ;;
    j) JOBS="$OPTARG" ;;
    x) EXEC_PATH="$OPTARG" ;;
    h) 
       echo "Usage: $0 -d <graph_dir> -l <min_size> -t <theta> -m <mode> [-j <jobs>] [-x <exec_path>]"
       exit 0
       ;;
    \?) echo "Invalid option -$OPTARG" >&2; exit 1 ;;
  esac
done

if [ -z "$GRAPH_DIR" ]; then
    echo "Error: Graph directory (-d) is required."
    exit 1
fi

if [ ! -f "$EXEC_PATH" ]; then
    echo "Error: Executable not found at $EXEC_PATH"
    echo "Please build it first or specify correct path with -x"
    exit 1
fi

# Timestamp for logs
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_DIR="logs/run_${TIMESTAMP}"
mkdir -p "$LOG_DIR"

echo "=== Starting Batch Run ==="
echo "Directory: $GRAPH_DIR"
echo "Params: l=$MIN_SIZE, theta=$THETA, mode=$MODE"
echo "Logs: $LOG_DIR"
echo "Concurrency: $JOBS"

# Function to process a single graph
process_graph() {
    local grh_file="$1"
    local graph_name=$(basename "$(dirname "$grh_file")")
    if [ "$graph_name" = "." ]; then
        graph_name=$(basename "$grh_file" .grh)
    fi
    
    local log_file="${LOG_DIR}/${graph_name}.log"
    
    echo "Processing $graph_name..."
    "$EXEC_PATH" "$grh_file" --minimum "$MIN_SIZE" --theta "$THETA" --mode "$MODE" > "$log_file" 2>&1
}

export -f process_graph
export MIN_SIZE THETA MODE EXEC_PATH LOG_DIR

# Find all .grh files and run in parallel
# Using find to locate .grh files. Assuming structure: root/graph_name/graph.grh or just root/graph.grh
if command -v parallel >/dev/null 2>&1; then
    find "$GRAPH_DIR" -name "*.grh" | parallel -j "$JOBS" process_graph {}
else
    # Fallback for when GNU parallel is not installed
    # Simple sequential loop if jobs=1, or background & wait if jobs > 1 (simplified)
    
    if [ "$JOBS" -eq 1 ]; then
        find "$GRAPH_DIR" -name "*.grh" | while read -r file; do
            process_graph "$file"
        done
    else
        # Simple parallelization with xargs if available
         find "$GRAPH_DIR" -name "*.grh" | xargs -P "$JOBS" -I {} bash -c 'process_graph "$@"' _ {}
    fi
fi

echo "=== Batch Run Complete ==="
echo "Logs available in $LOG_DIR"
