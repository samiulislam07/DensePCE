#!/usr/bin/env python3
"""
Windows-Compatible Graph Preprocessing Pipeline for Dense-PCE Project

This script provides a cross-platform preprocessing pipeline that handles the
BBkC preprocessing step internally for Windows compatibility.

Usage:
    python graph_pipeline.py <input_grh_file> [output_directory]
    python graph_pipeline.py --batch <input_directory> [output_base_directory]

Author: AI Assistant
Date: 2024
"""

import os
import sys
import subprocess
import shutil
import argparse
from pathlib import Path
import logging
import platform
import tempfile

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('graph_preprocessing.log'),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

class WindowsGraphPreprocessor:
    def __init__(self, project_root=None):
        """Initialize the graph preprocessor with project paths."""
        if project_root is None:
            # Assume script is in Dense-PCE-main directory
            self.project_root = Path(__file__).parent.absolute()
        else:
            self.project_root = Path(project_root).absolute()
        
        # Define tool paths
        self.bbkc_path = self.project_root / "BBkC"
        self.edgelist2binary_path = self.project_root / "EBBkC" / "Cohesive_subgraph_book" / "datasets" / "edgelist2binary"
        self.wsl_available = shutil.which("wsl") is not None
        
        # Check platform and tools
        self.is_windows = platform.system() == "Windows"
        self._verify_tools()
    
    def _verify_tools(self):
        """Verify that required tools are available or can be emulated."""
        if self.is_windows:
            if not self.wsl_available:
                raise EnvironmentError("WSL not detected on Windows. BBkC/edgelist2binary require WSL on Windows.")
            if not self.bbkc_path.exists():
                raise FileNotFoundError(f"BBkC executable not found at {self.bbkc_path}")
        else:
            # Check if BBkC exists and is executable
            bbkc_needs_rebuild = False
            
            if not self.bbkc_path.exists():
                logger.info("BBkC executable not found. Will attempt to build it.")
                bbkc_needs_rebuild = True
            elif not os.access(self.bbkc_path, os.X_OK):
                logger.warning(f"BBkC does not have execute permissions. Attempting to fix...")
                try:
                    os.chmod(self.bbkc_path, 0o755)
                    logger.info("Successfully added execute permissions to BBkC")
                except Exception as e:
                    logger.warning(f"Could not add execute permissions: {e}. Will rebuild.")
                    bbkc_needs_rebuild = True
            else:
                # Test if BBkC actually runs (catches library version issues)
                try:
                    result = subprocess.run(
                        [str(self.bbkc_path), "--help"],
                        capture_output=True,
                        text=True,
                        timeout=5
                    )
                    # Check for library version errors in stderr
                    if result.returncode != 0 or result.stderr:
                        error_msg = result.stderr.lower() if result.stderr else ""
                        # Check for common library version error patterns
                        if any(pattern in error_msg for pattern in [
                            "glibc", "glibcxx", "version `", "not found",
                            "library", "so.", "required by"
                        ]):
                            logger.warning("BBkC has library version compatibility issues detected in error output")
                            logger.info("Will attempt to rebuild BBkC on this system...")
                            bbkc_needs_rebuild = True
                    # If we get here and no library errors, it should work
                except (OSError, subprocess.TimeoutExpired) as e:
                    logger.warning(f"BBkC exists but cannot be executed (possibly library version mismatch): {e}")
                    logger.info("Will attempt to rebuild BBkC on this system...")
                    bbkc_needs_rebuild = True
            
            # Build BBkC if needed
            if bbkc_needs_rebuild:
                if not self._build_bbkc():
                    raise FileNotFoundError(f"Failed to build BBkC. Please build it manually: cd EBBkC/src && mkdir -p build && cd build && cmake .. && make && cp BBkC ../../../")
            
            # Check edgelist2binary
            if not self.edgelist2binary_path.exists():
                raise FileNotFoundError(f"edgelist2binary tool not found at {self.edgelist2binary_path}")
            
            # Also check execute permissions for edgelist2binary on non-Windows
            if not os.access(self.edgelist2binary_path, os.X_OK):
                logger.warning(f"edgelist2binary does not have execute permissions. Attempting to fix...")
                try:
                    os.chmod(self.edgelist2binary_path, 0o755)
                    logger.info("Successfully added execute permissions to edgelist2binary")
                except Exception as e:
                    logger.warning(f"Could not add execute permissions: {e}. Please run: chmod +x {self.edgelist2binary_path}")
        
        logger.info("Tool verification completed")

    def _build_bbkc(self):
        """Build BBkC from source on the current system."""
        logger.info("=== Building BBkC from source ===")
        
        ebbkc_src = self.project_root / "EBBkC" / "src"
        build_dir = ebbkc_src / "build"
        
        if not ebbkc_src.exists():
            logger.error(f"EBBkC source directory not found at {ebbkc_src}")
            return False
        
        try:
            # Create build directory
            build_dir.mkdir(parents=True, exist_ok=True)
            
            # Clean old CMake cache
            cmake_cache = build_dir / "CMakeCache.txt"
            if cmake_cache.exists():
                cmake_cache.unlink()
            
            # Configure with CMake
            logger.info("Configuring BBkC with CMake...")
            cmake_result = subprocess.run(
                ["cmake", ".."],
                cwd=build_dir,
                capture_output=True,
                text=True,
                timeout=60
            )
            
            if cmake_result.returncode != 0:
                logger.error(f"CMake configuration failed: {cmake_result.stderr}")
                return False
            
            # Build with make
            logger.info("Building BBkC...")
            make_result = subprocess.run(
                ["make", "-j"],
                cwd=build_dir,
                capture_output=True,
                text=True,
                timeout=300
            )
            
            if make_result.returncode != 0:
                logger.error(f"Build failed: {make_result.stderr}")
                return False
            
            # Find the built BBkC executable
            built_bbkc = build_dir / "BBkC"
            if not built_bbkc.exists():
                logger.error("BBkC executable was not created in build directory")
                return False
            
            # Copy to project root
            logger.info(f"Copying BBkC to {self.bbkc_path}...")
            shutil.copy2(built_bbkc, self.bbkc_path)
            os.chmod(self.bbkc_path, 0o755)
            
            logger.info("Successfully built and installed BBkC")
            return True
            
        except subprocess.TimeoutExpired:
            logger.error("Build process timed out")
            return False
        except Exception as e:
            logger.error(f"Error building BBkC: {e}")
            return False

    def _build_edgelist2binary(self):
        """Build edgelist2binary from source on the current system."""
        logger.info("=== Building edgelist2binary from source ===")
        
        edgelist2binary_dir = self.project_root / "EBBkC" / "Cohesive_subgraph_book" / "datasets"
        
        if not edgelist2binary_dir.exists():
            logger.error(f"edgelist2binary source directory not found at {edgelist2binary_dir}")
            return False
        
        if not (edgelist2binary_dir / "makefile").exists() and not (edgelist2binary_dir / "Makefile").exists():
            logger.error(f"Makefile not found in {edgelist2binary_dir}")
            return False
        
        try:
            # Create .obj directory if it doesn't exist
            obj_dir = edgelist2binary_dir / ".obj"
            obj_dir.mkdir(exist_ok=True)
            
            # Build with make
            logger.info("Building edgelist2binary...")
            make_result = subprocess.run(
                ["make"],
                cwd=edgelist2binary_dir,
                capture_output=True,
                text=True,
                timeout=300
            )
            
            if make_result.returncode != 0:
                logger.error(f"Build failed: {make_result.stderr or make_result.stdout}")
                return False
            
            # Check if the executable was created
            built_tool = edgelist2binary_dir / "edgelist2binary"
            if not built_tool.exists():
                logger.error("edgelist2binary executable was not created")
                return False
            
            # Ensure it's executable
            os.chmod(built_tool, 0o755)
            
            # Check if the path matches our expected location
            if built_tool != self.edgelist2binary_path:
                logger.warning(f"Built executable at {built_tool} doesn't match expected path {self.edgelist2binary_path}")
                # Path should match, but just in case
                if not self.edgelist2binary_path.exists():
                    logger.warning("Expected path doesn't exist, but build succeeded")
            
            logger.info("Successfully built edgelist2binary")
            return True
            
        except subprocess.TimeoutExpired:
            logger.error("Build process timed out")
            return False
        except Exception as e:
            logger.error(f"Error building edgelist2binary: {e}")
            return False


    def _to_wsl_path(self, path: Path) -> str:
        """Convert a Windows Path to WSL /mnt/<drive>/... form with proper escaping."""
        p = str(path.resolve())
        if len(p) > 1 and p[1] == ':' and p[0].isalpha():
            drive = p[0].lower()
            rest = p[2:].replace('\\', '/')
            # Don't escape here - we'll use single quotes around the path
            return f"/mnt/{drive}{rest}"
        return p.replace('\\', '/')
    
    def _run_wsl_command(self, command, working_dir=None, capture_output=True):
        """Run a command in WSL with proper working directory and error handling."""
        try:
            # Create a bash script that sets the working directory and runs the command
            if working_dir:
                wsl_working_dir = self._to_wsl_path(Path(working_dir))
                # Use single quotes to properly handle special characters in paths
                bash_script = f"cd '{wsl_working_dir}' && {command}"
            else:
                bash_script = command
            
            logger.info(f"Running WSL command: {bash_script}")
            
            result = subprocess.run(
                ["wsl", "bash", "-c", bash_script],
                capture_output=capture_output,
                text=True,
                timeout=300  # 5 minute timeout
            )
            
            if result.returncode != 0:
                logger.warning(f"WSL command failed with return code {result.returncode}")
                if capture_output:
                    logger.warning(f"STDOUT: {result.stdout}")
                    logger.warning(f"STDERR: {result.stderr}")
            
            return result
            
        except subprocess.TimeoutExpired:
            logger.error("WSL command timed out")
            return None
        except Exception as e:
            logger.error(f"Error running WSL command: {e}")
            return None
    
    def grh_to_edges(self, grh_file, edges_file):
        """Convert .grh file to .edges format."""
        logger.info(f"Converting {grh_file} to edges format")
        
        try:
            with open(grh_file, 'r') as f_in, open(edges_file, 'w') as f_out:
                for node_id, raw in enumerate(f_in):
                    line = raw.strip()
                    if not line or line == "[EOF]":
                        continue
                    # Robust tokenization: skip any non-integer tokens
                    tokens = line.split()
                    neighbors = []
                    for tok in tokens:
                        try:
                            neighbors.append(int(tok))
                        except ValueError:
                            # skip garbage token
                            continue
                    # Write edges (only forward edges to avoid duplicates)
                    for neighbor in neighbors:
                        if neighbor > node_id:
                            f_out.write(f"{node_id} {neighbor}\n")
            
            logger.info(f"Successfully converted to {edges_file}")
            return True
            
        except Exception as e:
            logger.error(f"Error converting grh to edges: {e}")
            return False
    
    def edges_to_clean(self, edges_file, output_dir):
        """Clean edges using BBkC preprocessing with improved path handling."""
        clean_file = Path(output_dir) / f"{Path(edges_file).stem}.clean"
        
        logger.info(f"Cleaning edges file {edges_file} using BBkC")
        try:
            if self.is_windows:
                # Use WSL to run BBkC with proper working directory
                linux_bbkc = self._to_wsl_path(self.bbkc_path)
                linux_edges = self._to_wsl_path(Path(edges_file))
                # Use single quotes around paths to handle special characters
                wsl_cmd = f"'{linux_bbkc}' p '{linux_edges}'"
                
                result = self._run_wsl_command(wsl_cmd, working_dir=output_dir)
                if result is None:
                    return False
            else:
                # BBkC creates output in the same directory as the input file
                # Use absolute path and don't set cwd
                edges_abs = Path(edges_file).resolve()
                
                # Ensure the input file exists
                if not edges_abs.exists():
                    logger.error(f"Input edges file does not exist: {edges_abs}")
                    return False
                
                cmd = [str(self.bbkc_path), "p", str(edges_abs)]
                logger.info(f"Running command: {' '.join(cmd)}")
                result = subprocess.run(cmd, capture_output=True, text=True)
                
                # Log output for debugging
                if result.stdout:
                    logger.debug(f"BBkC stdout: {result.stdout}")
                if result.stderr:
                    logger.debug(f"BBkC stderr: {result.stderr}")
                logger.debug(f"BBkC return code: {result.returncode}")

            if result.returncode != 0:
                logger.error(f"BBkC preprocessing failed: {result.stderr or result.stdout}")
                return False

            if clean_file.exists():
                logger.info(f"Successfully created clean file: {clean_file}")
                return clean_file
            else:
                logger.error("Clean file was not created by BBkC")
                return False

        except Exception as e:
            logger.error(f"Error running BBkC: {e}")
            return False
    
    def clean_to_binary(self, clean_file, output_dir):
        """Convert .clean file to binary format using edgelist2binary with improved path handling."""
        clean_path = Path(clean_file)
        out_dir = Path(output_dir)
        b_adj_file = out_dir / "b_adj.bin"
        b_degree_file = out_dir / "b_degree.bin"

        # 1) Native edgelist2binary on non-Windows
        if self.edgelist2binary_path.exists() and not self.is_windows:
            logger.info(f"Converting {clean_file} to binary format using native edgelist2binary")
            try:
                # edgelist2binary expects: <tool> <output_dir> <clean_file>
                # Since cwd is set to out_dir, use '.' for output dir and just the filename
                cmd = [str(self.edgelist2binary_path), ".", clean_path.name]
                logger.info(f"Running command: {' '.join(cmd)} (cwd={out_dir})")
                result = subprocess.run(cmd, capture_output=True, text=True, cwd=out_dir)
                
                # Check if output files were created first (some tools succeed even with non-zero return codes)
                if b_adj_file.exists() and b_degree_file.exists():
                    logger.info("Successfully created binary files using native edgelist2binary")
                    return True
                
                # If files weren't created, check for errors
                if result.returncode != 0:
                    error_msg = result.stderr or result.stdout or ""
                    # Check if this is a library version issue
                    if any(pattern in error_msg.lower() for pattern in [
                        "glibc", "glibcxx", "version `", "not found",
                        "library", "so.", "required by"
                    ]):
                        logger.error("edgelist2binary has library version compatibility issues. Attempting to rebuild...")
                        if self._build_edgelist2binary():
                            # Retry after rebuilding
                            logger.info("Retrying edgelist2binary conversion after rebuild...")
                            result = subprocess.run(cmd, capture_output=True, text=True, cwd=out_dir)
                            # Check file existence after retry
                            if b_adj_file.exists() and b_degree_file.exists():
                                logger.info("Successfully created binary files using native edgelist2binary after rebuild")
                                return True
                            if result.returncode != 0:
                                logger.error(f"edgelist2binary failed after rebuild: {result.stderr or result.stdout}")
                        else:
                            logger.error("Failed to rebuild edgelist2binary. Please build it manually.")
                    else:
                        logger.warning(f"edgelist2binary failed: {error_msg}")
                else:
                    # Return code is 0 but files don't exist - this is unexpected
                    logger.error(f"edgelist2binary completed successfully but binary files were not created")
                    logger.debug(f"edgelist2binary stdout: {result.stdout}")
                    logger.debug(f"edgelist2binary stderr: {result.stderr}")
            except Exception as e:
                logger.warning(f"edgelist2binary error: {e}")

        # 2) Windows: try WSL with improved path handling
        if self.is_windows and self.wsl_available:
            try:
                linux_tool = self._to_wsl_path(self.edgelist2binary_path)
                wsl_cmd = f"'{linux_tool}' . '{clean_path.name}'"
                logger.info(f"Converting via WSL: {wsl_cmd}")
                
                result = self._run_wsl_command(wsl_cmd, working_dir=out_dir)
                if result and result.returncode == 0:
                    if b_adj_file.exists() and b_degree_file.exists():
                        logger.info("Successfully created binary files via WSL edgelist2binary")
                        return True
                    else:
                        logger.warning("WSL edgelist2binary completed but binary files not found")
                else:
                    logger.warning("WSL edgelist2binary failed")
            except Exception as e:
                logger.warning(f"WSL conversion error: {e}")

        logger.error("edgelist2binary not available or failed. Cannot produce binaries.")
        return False
    
    def process_single_graph(self, grh_file, output_dir=None):
        """Process a single graph file through the complete pipeline."""
        grh_path = Path(grh_file)
        
        if not grh_path.exists():
            logger.error(f"Input file does not exist: {grh_file}")
            return False
        
        # Use the same directory as the input file
        if output_dir is None:
            output_dir = grh_path.parent
        else:
            output_dir = Path(output_dir)
        
        output_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f"Processing {grh_file} in directory {output_dir}")
        
        # Step 1: Convert .grh to .edges (in same directory)
        edges_file = output_dir / f"{grh_path.stem}.edges"
        if not self.grh_to_edges(grh_path, edges_file):
            return False
        
        # Step 2: Convert .edges to .clean
        clean_file = self.edges_to_clean(edges_file, output_dir)
        if not clean_file:
            return False
        
        # Step 3: Convert .clean to binary format
        if not self.clean_to_binary(clean_file, output_dir):
            return False
        
        logger.info(f"Successfully processed {grh_file}")
        self._print_output_summary(output_dir)
        return True
    
    def process_batch(self, input_dir, output_base_dir=None):
        """Process all .grh files in a directory."""
        input_path = Path(input_dir)
        
        if not input_path.exists():
            logger.error(f"Input directory does not exist: {input_dir}")
            return False
        
        # For batch processing, use the same directories as the input files
        # unless a specific output base directory is provided
        if output_base_dir is not None:
            output_base_dir = Path(output_base_dir)
        
        # Find all .grh files
        grh_files = list(input_path.glob("*.grh"))
        
        if not grh_files:
            logger.warning(f"No .grh files found in {input_dir}")
            return True
        
        logger.info(f"Found {len(grh_files)} .grh files to process")
        
        success_count = 0
        for grh_file in grh_files:
            try:
                # If output_base_dir is specified, create subdirectory for each graph
                if output_base_dir is not None:
                    output_dir = output_base_dir / grh_file.parent.name
                    output_dir.mkdir(parents=True, exist_ok=True)
                else:
                    output_dir = None  # Will use same directory as input file
                
                if self.process_single_graph(grh_file, output_dir):
                    success_count += 1
                else:
                    logger.error(f"Failed to process {grh_file}")
            except Exception as e:
                logger.error(f"Exception processing {grh_file}: {e}")
        
        logger.info(f"Batch processing complete: {success_count}/{len(grh_files)} files processed successfully")
        return success_count == len(grh_files)
    
    def _print_output_summary(self, output_dir):
        """Print summary of output files."""
        output_path = Path(output_dir)
        files = list(output_path.glob("*"))
        
        logger.info(f"Output directory: {output_dir}")
        logger.info("Generated files:")
        for file in sorted(files):
            size = file.stat().st_size if file.is_file() else "DIR"
            logger.info(f"  {file.name} ({size} bytes)")

def main():
    parser = argparse.ArgumentParser(
        description="Windows-Compatible Graph Preprocessing Pipeline for Dense-PCE Project",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process single graph file (creates .edges, .clean, .bin files in same directory)
  python graph_pipeline.py protein/Protein.grh

  # Process single graph with custom output directory
  python graph_pipeline.py graph.grh ./processed_graphs/

  # Process all .grh files in a directory (creates files in each file's directory)
  python graph_pipeline.py --batch synth-graphs-1000/

  # Process batch with custom output base directory
  python graph_pipeline.py --batch synth-graphs-1000/ ./batch_processed/
        """
    )
    
    parser.add_argument('input', help='Input .grh file or directory (with --batch)')
    parser.add_argument('output', nargs='?', help='Output directory (optional)')
    parser.add_argument('--batch', action='store_true', help='Process all .grh files in input directory')
    parser.add_argument('--project-root', help='Project root directory (auto-detected if not specified)')
    
    args = parser.parse_args()
    
    try:
        preprocessor = WindowsGraphPreprocessor(args.project_root)
        
        if args.batch:
            success = preprocessor.process_batch(args.input, args.output)
        else:
            success = preprocessor.process_single_graph(args.input, args.output)
        
        if success:
            logger.info("Pipeline completed successfully!")
            return 0
        else:
            logger.error("Pipeline failed!")
            return 1
            
    except Exception as e:
        logger.error(f"Fatal error: {e}")
        return 1

if __name__ == "__main__":
    sys.exit(main())
