#!/bin/bash

# Build script for the integrated Dense-PCE + EBBkC system
# This script compiles EBBkC as a library and links it with a specified Dense-PCE source file

set -e  # Exit on any error

# Check for command line argument, default to dense-pce-r1-m1.cpp
if [ $# -eq 0 ]; then
    SOURCE_FILE="dense-pce-r1-m1.cpp"
    echo "No source file specified, defaulting to: $SOURCE_FILE"
else
    SOURCE_FILE="$1"
fi

# Check if source file exists
if [ ! -f "$SOURCE_FILE" ]; then
    echo "Error: Source file '$SOURCE_FILE' not found"
    exit 1
fi

# Extract base name for output executable
BASE_NAME=$(basename "$SOURCE_FILE" .cpp)
OUTPUT_NAME="${BASE_NAME}-integrated"

echo "=== Building Integrated Dense-PCE + EBBkC System ==="
echo "Source file: $SOURCE_FILE"
echo "Output executable: $OUTPUT_NAME"

# Create build directory if it doesn't exist
mkdir -p build_integrated

# Check if EBBkC is already compiled
EBBKC_LIB="build_integrated/libebbkc_core.a"
EBBKC_COMMON_UTILS="EBBkC/src/build/libcommon-utils.a"
EBBKC_GRAPH_PREPROC="EBBkC/src/build/libgraph-pre-processing.a"

if [ -f "$EBBKC_LIB" ] && [ -f "$EBBKC_COMMON_UTILS" ] && [ -f "$EBBKC_GRAPH_PREPROC" ]; then
    echo "=== EBBkC library already compiled, skipping EBBkC build ==="
    cd build_integrated
else
    echo "=== Step 1: Building EBBkC as a library ==="
    
    # Clean up old build artifacts to avoid path conflicts
    echo "=== Cleaning up old build artifacts ==="
    # Only remove the EBBkC build directory, keep build_integrated with existing executables
    rm -rf EBBkC/src/build
    
    cd build_integrated
    
    # Build EBBkC core library
    cd ../EBBkC/src
    mkdir -p build
    cd build
    
    # Clean any existing CMake cache to avoid path conflicts
    rm -f CMakeCache.txt
    rm -rf CMakeFiles
    
    # Configure and build EBBkC
    echo "Configuring EBBkC with CMake..."
    cmake .. -DCMAKE_BUILD_TYPE=Release
    echo "Building EBBkC..."
    make -j$(nproc)
    
    # Copy the library to main directory
    cp libebbkc_core.a ../../../build_integrated/
    cd ../../../build_integrated
fi

echo "=== Step 2: Compiling $SOURCE_FILE with EBBkC library ==="

# Compile the specified source file with EBBkC library
# Try different OpenMP library names for different systems
echo "Compiling integrated executable..."
g++ -std=c++17 -O3 -march=native -fopenmp \
    -I../EBBkC/src \
    -I../EBBkC/src/truss/dependencies/sparsepp \
    -I../EBBkC/src/truss/dependencies/libpopcnt \
    -I../EBBkC/src/truss/dependencies \
    -I../EBBkC/src/truss/util \
    -I../EBBkC/src/truss/decompose \
    ../$SOURCE_FILE \
    libebbkc_core.a \
    ../EBBkC/src/build/libcommon-utils.a \
    ../EBBkC/src/build/libgraph-pre-processing.a \
    -ltbb -lpthread \
    -o $OUTPUT_NAME 2>/dev/null || \
g++ -std=c++17 -O3 -march=native -fopenmp \
    -I../EBBkC/src \
    -I../EBBkC/src/truss/dependencies/sparsepp \
    -I../EBBkC/src/truss/dependencies/libpopcnt \
    -I../EBBkC/src/truss/dependencies \
    -I../EBBkC/src/truss/util \
    -I../EBBkC/src/truss/decompose \
    ../$SOURCE_FILE \
    libebbkc_core.a \
    ../EBBkC/src/build/libcommon-utils.a \
    ../EBBkC/src/build/libgraph-pre-processing.a \
    -lgomp -ltbb -lpthread \
    -o $OUTPUT_NAME

echo "=== Build completed successfully! ==="
echo "Executable: ./build_integrated/$OUTPUT_NAME"
echo ""
echo "Usage example:"
echo "./build_integrated/$OUTPUT_NAME testGraph/gplus/gplus.grh --minimum 10 --theta 0.9"
