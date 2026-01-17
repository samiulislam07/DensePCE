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

    # --- Performance capture (GNU time -v) ---
    # Prefer /usr/bin/time (GNU time on Ubuntu/WSL). Fallback to `time` if needed.
    local time_bin="/usr/bin/time"
    if [ ! -x "$time_bin" ]; then
        time_bin="$(command -v time || true)"
    fi

    {
        echo "========================================"
        echo "Graph: $grh_file"
        echo "Params: --minimum $MIN_SIZE --theta $THETA --mode $MODE"
        echo "Executable: $EXEC_PATH"
        echo "Start: $(date -u +"%Y-%m-%dT%H:%M:%SZ")"
        echo "========================================"
        echo ""
    } >> "$log_file"

    local tmp_file
    tmp_file="$(mktemp 2>/dev/null || mktemp -t densepce_time)"

    # Run and capture program output + GNU time stats (if available)
    if [ -n "$time_bin" ] && "$time_bin" -v -o "$tmp_file" true >/dev/null 2>&1; then
        # GNU time supports -v and -o
        "$time_bin" -v -o "$tmp_file" \
            "$EXEC_PATH" "$grh_file" --minimum "$MIN_SIZE" --theta "$THETA" --mode "$MODE" \
            >> "$log_file" 2>&1
    else
        # Fallback: no GNU time. Run without detailed stats.
        echo "[WARN] GNU time (-v) not available; running without CPU/memory metrics." >> "$log_file"
        "$EXEC_PATH" "$grh_file" --minimum "$MIN_SIZE" --theta "$THETA" --mode "$MODE" \
            >> "$log_file" 2>&1
        : > "$tmp_file"
    fi

    # Extract metrics (GNU time -v format)
    local user_time sys_time elapsed_time max_memory_kb total_cpu peak_mem_display
    user_time="$(sed -n 's/.*User time (seconds): //p' "$tmp_file" | tail -n 1)"
    sys_time="$(sed -n 's/.*System time (seconds): //p' "$tmp_file" | tail -n 1)"
    elapsed_time="$(sed -n 's/.*Elapsed (wall clock) time (h:mm:ss or m:ss): //p' "$tmp_file" | tail -n 1)"
    max_memory_kb="$(sed -n 's/.*Maximum resident set size (kbytes): //p' "$tmp_file" | tail -n 1)"

    total_cpu="N/A"
    if [[ "$user_time" =~ ^[0-9.]+$ && "$sys_time" =~ ^[0-9.]+$ ]]; then
        total_cpu="$(awk -v u="$user_time" -v s="$sys_time" 'BEGIN{printf "%.6f", (u+s)}')"
    fi

    peak_mem_display="$max_memory_kb"
    if [[ "$max_memory_kb" =~ ^[0-9]+$ ]]; then
        peak_mem_display="$(awk -v kb="$max_memory_kb" 'BEGIN{printf "%.2f MB", (kb/1024.0)}')"
    fi

    {
        echo ""
        echo "=== Performance Summary (GNU time -v) ==="
        echo "CPU Time (user+sys): ${total_cpu}s"
        echo "User Time: ${user_time}s"
        echo "System Time: ${sys_time}s"
        echo "Real Time: ${elapsed_time}"
        echo "Peak Memory: ${peak_mem_display}"
        echo "End: $(date -u +"%Y-%m-%dT%H:%M:%SZ")"
        echo ""
    } >> "$log_file"

    rm -f "$tmp_file"
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
