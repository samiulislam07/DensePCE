#!/bin/bash

# Build script for the integrated Dense-PCE + EBBkC system
# This script ONLY compiles the main source file, assuming EBBkC libraries are already built.

set -e  # Exit on any error

# --- Prerequisites ---
# This script assumes you have already run the original, full build script at least once.
# It requires the following files to exist from that initial build:
#   - ./build_integrated/libebbkc_core.a
#   - ./EBBkC/src/build/libcommon-utils.a
#   - ./EBBkC/src/build/libgraph-pre-processing.a

# Check for command line argument
if [ $# -eq 0 ]; then
    echo "Usage: $0 <source_file.cpp>"
    echo "Example: $0 dense-pce-ab.cpp"
    exit 1
fi

SOURCE_FILE="$1"

# Check if source file exists
if [ ! -f "$SOURCE_FILE" ]; then
    echo "Error: Source file '$SOURCE_FILE' not found"
    exit 1
fi

# Extract base name for output executable
BASE_NAME=$(basename "$SOURCE_FILE" .cpp)
OUTPUT_NAME="${BASE_NAME}-integrated"

echo "=== Compiling Integrated Executable (assumes EBBkC is already built) ==="
echo "Source file: $SOURCE_FILE"
echo "Output executable: $OUTPUT_NAME"

# Create the output directory if it doesn't exist
mkdir -p build_integrated

echo "=== Compiling $SOURCE_FILE with pre-built EBBkC libraries... ==="

# Compile the specified source file, linking against the pre-built EBBkC libraries.
# The '2>/dev/null ||' part attempts the first command and, if it fails, tries the
# second one with a different OpenMP library flag (-lgomp).
g++ -std=c++17 -O3 -march=native -fopenmp \
    -I./EBBkC/src \
    -I./EBBkC/src/truss/dependencies/sparsepp \
    -I./EBBkC/src/truss/dependencies/libpopcnt \
    -I./EBBkC/src/truss/dependencies \
    -I./EBBkC/src/truss/util \
    -I./EBBkC/src/truss/decompose \
    "$SOURCE_FILE" \
    ./build_integrated/libebbkc_core.a \
    ./EBBkC/src/build/libcommon-utils.a \
    ./EBBkC/src/build/libgraph-pre-processing.a \
    -ltbb -lpthread \
    -o ./build_integrated/"$OUTPUT_NAME" 2>/dev/null || \
g++ -std=c++17 -O3 -march=native -fopenmp \
    -I./EBBkC/src \
    -I./EBBkC/src/truss/dependencies/sparsepp \
    -I./EBBkC/src/truss/dependencies/libpopcnt \
    -I./EBBkC/src/truss/dependencies \
    -I./EBBkC/src/truss/util \
    -I./EBBkC/src/truss/decompose \
    "$SOURCE_FILE" \
    ./build_integrated/libebbkc_core.a \
    ./EBBkC/src/build/libcommon-utils.a \
    ./EBBkC/src/build/libgraph-pre-processing.a \
    -lgomp -ltbb -lpthread \
    -o ./build_integrated/"$OUTPUT_NAME"

echo "=== Build completed successfully! ==="
echo "Executable: ./build_integrated/$OUTPUT_NAME"
echo ""
echo "Usage example:"
echo "./build_integrated/$OUTPUT_NAME testGraph/gplus/gplus.grh --minimum 10 --theta 0.9"