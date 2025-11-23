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
            if not self.bbkc_path.exists():
                raise FileNotFoundError(f"BBkC executable not found at {self.bbkc_path}")
        
        # Check edgelist2binary - if missing we will try WSL or fallback to internal
        if not self.edgelist2binary_path.exists():
            raise FileNotFoundError(f"edgelist2binary tool not found at {self.edgelist2binary_path}")
        
        logger.info("Tool verification completed")

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
                cmd = [str(self.bbkc_path), "p", str(edges_file)]
                logger.info(f"Running command: {' '.join(cmd)}")
                result = subprocess.run(cmd, capture_output=True, text=True, cwd=output_dir)

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
                cmd = [str(self.edgelist2binary_path), str(out_dir), clean_path.name]
                logger.info(f"Running command: {' '.join(cmd)} (cwd={out_dir})")
                result = subprocess.run(cmd, capture_output=True, text=True, cwd=out_dir)
                if result.returncode != 0:
                    logger.warning(f"edgelist2binary failed: {result.stderr or result.stdout}")
                else:
                    if b_adj_file.exists() and b_degree_file.exists():
                        logger.info("Successfully created binary files using native edgelist2binary")
                        return True
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
