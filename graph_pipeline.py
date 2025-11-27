#!/usr/bin/env python3
"""
Cross-Platform Graph Preprocessing Pipeline for Dense-PCE Project

This script provides a robust preprocessing pipeline that works on Linux, macOS,
and Windows (via WSL or native if tools are available).

It automatically handles:
1. Tool detection and building (BBkC, edgelist2binary)
2. Platform-specific execution (WSL on Windows, native on Unix)
3. Path conversions and file handling

Usage:
    python graph_pipeline.py <input_grh_file> [output_directory]
    python graph_pipeline.py --batch <input_directory> [output_base_directory]

Author: AI Assistant
Date: 2025
"""

import os
import sys
import subprocess
import shutil
import argparse
from pathlib import Path
import logging
import platform
import stat

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

class GraphPipeline:
    def __init__(self, project_root=None):
        """Initialize the graph preprocessor with project paths."""
        if project_root is None:
            # Assume script is in Dense-PCE-main directory
            self.project_root = Path(__file__).parent.absolute()
        else:
            self.project_root = Path(project_root).absolute()
        
        self.system = platform.system()
        self.is_windows = self.system == "Windows"
        
        # Define tool paths
        # On Windows, these might be run via WSL, so we track the native paths here
        self.bbkc_path = self.project_root / "BBkC"
        
        # The edgelist2binary tool is inside EBBkC
        self.edgelist2binary_dir = self.project_root / "EBBkC" / "Cohesive_subgraph_book" / "datasets"
        self.edgelist2binary_src = self.edgelist2binary_dir / "main.cpp"
        self.edgelist2binary_path = self.edgelist2binary_dir / "edgelist2binary"
        
        # Check for WSL availability on Windows
        self.wsl_available = False
        if self.is_windows:
            self.wsl_available = shutil.which("wsl") is not None
            if not self.wsl_available:
                logger.warning("WSL not found on Windows. Will attempt to run tools natively.")
        
        self._verify_and_setup_tools()
    
    def _verify_and_setup_tools(self):
        """Verify that required tools are available, building them if necessary."""
        logger.info(f"Initializing pipeline on {self.system}...")
        
        # 1. Check/Build BBkC
        if not self._check_tool_runnable(self.bbkc_path, "--help"):
            logger.info("BBkC not ready or not runnable. Attempting to build...")
            if not self._build_bbkc():
                # On Windows without WSL, we can't easily build, so this is critical
                if self.is_windows and not self.wsl_available:
                    raise EnvironmentError("BBkC not found and cannot build natively on Windows without WSL.")
                elif not self.is_windows:
                    raise RuntimeError("Failed to build BBkC. Check logs for details.")
        
        # 2. Check/Build edgelist2binary
        # Note: edgelist2binary takes args, so running with no args usually exits with 1 or prints usage
        # We'll assume if it exists and we can run it (even if it returns error code), it's binary compatible
        if not self._check_tool_runnable(self.edgelist2binary_path, args=[], expect_success=False):
            logger.info("edgelist2binary not ready. Attempting to build...")
            if not self._build_edgelist2binary():
                 if self.is_windows and not self.wsl_available:
                    raise EnvironmentError("edgelist2binary not found and cannot build natively on Windows without WSL.")
                 elif not self.is_windows:
                    raise RuntimeError("Failed to build edgelist2binary.")

        logger.info("All tools verified successfully.")

    def _check_tool_runnable(self, tool_path: Path, args=None, expect_success=True):
        """
        Check if a tool exists and can be executed on the CURRENT platform.
        For Windows, this checks if it runs via WSL if configured, or natively.
        """
        if not tool_path.exists():
            return False
            
        if args is None:
            args = []
            
        if isinstance(args, str):
            args = [args]

        try:
            # Ensure executable permissions on Unix-like systems
            if not self.is_windows:
                st = os.stat(tool_path)
                os.chmod(tool_path, st.st_mode | stat.S_IEXEC)

            if self.is_windows and self.wsl_available:
                # On Windows with WSL, we try running via WSL first
                # This verifies the binary is a Linux binary (which is what we expect with WSL)
                # If the user has a native Windows binary, this might fail, and we should fallback?
                # Strategy: Try WSL. If fail, try native.
                
                # Convert path for WSL
                wsl_path = self._to_wsl_path(tool_path)
                cmd = ["wsl", str(wsl_path)] + args
                
                result = subprocess.run(
                    cmd, 
                    capture_output=True, 
                    timeout=5
                )
                
                # If exec format error (running Windows exe in WSL), return False so we rebuild/recheck
                if result.returncode == 126 or "Exec format error" in result.stderr.decode():
                    return False
                    
                if expect_success and result.returncode != 0:
                    return False
                
                return True
                
            else:
                # Native execution (Linux, Mac, or Windows native)
                cmd = [str(tool_path)] + args
                result = subprocess.run(
                    cmd,
                    capture_output=True,
                    timeout=5
                )
                
                if expect_success and result.returncode != 0:
                    return False
                    
                return True
                
        except (OSError, subprocess.SubprocessError) as e:
            logger.debug(f"Tool check failed for {tool_path}: {e}")
            return False

    def _build_bbkc(self):
        """Build BBkC from source."""
        if self.is_windows and self.wsl_available:
            # Build inside WSL
            logger.info("Building BBkC via WSL...")
            return self._build_bbkc_wsl()
        elif self.is_windows:
            logger.error("Native Windows build for BBkC not fully supported automatically. Please build manually.")
            return False
        else:
            # Native Linux/Mac Build
            return self._build_bbkc_native()

    def _build_bbkc_native(self):
        """Build BBkC on Linux/macOS."""
        ebbkc_src = self.project_root / "EBBkC" / "src"
        build_dir = ebbkc_src / "build"
        
        if not ebbkc_src.exists():
            logger.error(f"Source not found: {ebbkc_src}")
            return False
            
        # Check for CMake and Make
        if not shutil.which("cmake") or not shutil.which("make"):
            logger.error("cmake and make are required to build BBkC.")
            return False

        try:
            build_dir.mkdir(parents=True, exist_ok=True)
            
            # Clean cache
            if (build_dir / "CMakeCache.txt").exists():
                (build_dir / "CMakeCache.txt").unlink()
            
            logger.info("Running cmake...")
            subprocess.run(["cmake", ".."], cwd=build_dir, check=True, capture_output=True)
            
            logger.info("Running make...")
            subprocess.run(["make", "-j"], cwd=build_dir, check=True, capture_output=True)
            
            built_bin = build_dir / "BBkC"
            if built_bin.exists():
                shutil.copy2(built_bin, self.bbkc_path)
                return True
            return False
            
        except subprocess.CalledProcessError as e:
            logger.error(f"Build failed: {e}")
            # Try to print stderr for debugging
            if hasattr(e, 'stderr') and e.stderr:
                logger.error(e.stderr.decode())
            return False

    def _build_bbkc_wsl(self):
        """Build BBkC using WSL commands."""
        # Paths in WSL format
        src_dir = self._to_wsl_path(self.project_root / "EBBkC" / "src")
        build_dir = f"{src_dir}/build"
        
        cmds = [
            f"mkdir -p '{build_dir}'",
            f"cd '{build_dir}'",
            "cmake ..",
            "make -j",
            f"cp BBkC '{self._to_wsl_path(self.project_root)}/BBkC'"
        ]
        
        full_cmd = " && ".join(cmds)
        result = self._run_wsl_cmd_raw(full_cmd)
        
        return result.returncode == 0 and self.bbkc_path.exists()

    def _build_edgelist2binary(self):
        """Build edgelist2binary tool."""
        if self.is_windows and self.wsl_available:
            return self._build_edgelist2binary_wsl()
        else:
            return self._build_edgelist2binary_native()
            
    def _build_edgelist2binary_native(self):
        logger.info("Building edgelist2binary natively...")
        compiler = shutil.which("g++") or shutil.which("clang++") or shutil.which("c++")
        
        if not compiler:
            logger.error("No C++ compiler found (g++, clang++, or c++).")
            return False
            
        cmd = [
            compiler,
            "-O3",
            "-std=c++11",
            str(self.edgelist2binary_src),
            "-o",
            str(self.edgelist2binary_path)
        ]
        
        try:
            subprocess.run(cmd, check=True, cwd=self.edgelist2binary_dir)
            return self.edgelist2binary_path.exists()
        except subprocess.CalledProcessError as e:
            logger.error(f"edgelist2binary build failed: {e}")
            return False

    def _build_edgelist2binary_wsl(self):
        logger.info("Building edgelist2binary via WSL...")
        wsl_src = self._to_wsl_path(self.edgelist2binary_src)
        wsl_out = self._to_wsl_path(self.edgelist2binary_path)
        
        # Assuming g++ is installed in WSL
        cmd = f"g++ -O3 -std=c++11 '{wsl_src}' -o '{wsl_out}'"
        
        result = self._run_wsl_cmd_raw(cmd)
        return result.returncode == 0

    # --- Helper Methods ---

    def _to_wsl_path(self, path: Path) -> str:
        """Convert path to WSL format (/mnt/c/...)."""
        p = str(path.resolve())
        if len(p) > 1 and p[1] == ':' and p[0].isalpha():
            drive = p[0].lower()
            rest = p[2:].replace('\\', '/')
            return f"/mnt/{drive}{rest}"
        return p.replace('\\', '/')

    def _run_wsl_cmd_raw(self, bash_command):
        """Run a raw bash command string in WSL."""
        return subprocess.run(
            ["wsl", "bash", "-c", bash_command],
            capture_output=True,
            text=True
        )

    def _run_command(self, cmd_args, cwd=None):
        """
        Run a command transparently handling WSL vs Native.
        cmd_args[0] should be the executable path (Path object or string).
        """
        executable = Path(cmd_args[0])
        args = cmd_args[1:]
        
        if self.is_windows and self.wsl_available:
            # Use WSL
            wsl_exec = self._to_wsl_path(executable)
            wsl_cwd = self._to_wsl_path(cwd) if cwd else None
            
            # Convert args: if Path object, convert to WSL path
            wsl_args = []
            for arg in args:
                if isinstance(arg, Path):
                    wsl_args.append(self._to_wsl_path(arg))
                else:
                    wsl_args.append(str(arg))
            
            # Construct argument string safely
            str_args = " ".join([f"'{arg}'" for arg in wsl_args])
            
            full_cmd = f"'{wsl_exec}' {str_args}"
            if wsl_cwd:
                full_cmd = f"cd '{wsl_cwd}' && {full_cmd}"
            
            logger.debug(f"Executing WSL command: {full_cmd}")
            return self._run_wsl_cmd_raw(full_cmd)
        else:
            # Native
            # Convert executable path to string
            native_cmd = [str(executable)] + [str(a) for a in args]
            logger.debug(f"Executing native command: {native_cmd}")
            return subprocess.run(
                native_cmd,
                cwd=cwd,
                capture_output=True,
                text=True
            )

    # --- Processing Logic ---

    def grh_to_edges(self, grh_file, edges_file):
        """Convert .grh file to .edges format (Python-based, platform independent)."""
        logger.info(f"Converting {grh_file.name} to edges format")
        try:
            with open(grh_file, 'r') as f_in, open(edges_file, 'w') as f_out:
                for node_id, raw in enumerate(f_in):
                    line = raw.strip()
                    if not line or line == "[EOF]":
                        continue
                    tokens = line.split()
                    for tok in tokens:
                        try:
                            neighbor = int(tok)
                            # Forward edges only
                            if neighbor > node_id:
                                f_out.write(f"{node_id} {neighbor}\n")
                        except ValueError:
                            continue
            return True
        except Exception as e:
            logger.error(f"Error converting grh to edges: {e}")
            return False

    def edges_to_clean(self, edges_file, output_dir):
        """Run BBkC preprocessing."""
        logger.info(f"Cleaning edges file {edges_file.name}")
        
        # Command: BBkC p <file>
        # Note: We need to pass absolute paths to ensure tools find them
        result = self._run_command(
            [self.bbkc_path, "p", edges_file.absolute()],
            cwd=output_dir
        )
        
        if result.returncode != 0:
            logger.error(f"BBkC failed. Exit code: {result.returncode}")
            logger.error(f"STDOUT: {result.stdout}")
            logger.error(f"STDERR: {result.stderr}")
            return None
            
        clean_file = output_dir / f"{edges_file.stem}.clean"
        if clean_file.exists():
            return clean_file
        else:
            logger.error("BBkC did not produce .clean file")
            return None

    def clean_to_binary(self, clean_file, output_dir):
        """Run edgelist2binary."""
        logger.info(f"Creating binaries from {clean_file.name}")
        
        # Command: edgelist2binary <output_dir> <clean_filename>
        # IMPORTANT: edgelist2binary typically expects the filename to be in the CWD or just the name
        # Based on previous code: cmd = [tool, out_dir, clean_file.name] with cwd=out_dir
        
        result = self._run_command(
            [self.edgelist2binary_path, output_dir, clean_file.name],
            cwd=output_dir
        )
        
        if result.returncode != 0:
            logger.error(f"edgelist2binary failed. Exit code: {result.returncode}")
            logger.error(f"STDOUT: {result.stdout}")
            logger.error(f"STDERR: {result.stderr}")
            return False
            
        return (output_dir / "b_adj.bin").exists()

    def process_single_graph(self, grh_file, output_dir=None):
        grh_path = Path(grh_file).absolute()
        if not grh_path.exists():
            logger.error(f"File not found: {grh_path}")
            return False
            
        if output_dir:
            out_path = Path(output_dir).absolute()
        else:
            out_path = grh_path.parent
            
        out_path.mkdir(parents=True, exist_ok=True)
        
        # 1. GRH -> Edges
        edges_file = out_path / f"{grh_path.stem}.edges"
        if not self.grh_to_edges(grh_path, edges_file):
            return False
            
        # 2. Edges -> Clean
        clean_file = self.edges_to_clean(edges_file, out_path)
        if not clean_file:
            return False
            
        # 3. Clean -> Binary
        if not self.clean_to_binary(clean_file, out_path):
            return False
            
        logger.info(f"Successfully processed {grh_path.name}")
        return True

    def process_batch(self, input_dir, output_base=None):
        in_dir = Path(input_dir).absolute()
        if not in_dir.exists():
            logger.error(f"Directory not found: {in_dir}")
            return False
            
        grh_files = list(in_dir.glob("*.grh"))
        if not grh_files:
            logger.warning("No .grh files found")
            return True
            
        success_count = 0
        for f in grh_files:
            out_dir = None
            if output_base:
                # Maintain directory structure relative to input, or just flat?
                # Let's just put it in a folder named after the file if using base
                out_dir = Path(output_base).absolute() / f.parent.name
                
            if self.process_single_graph(f, out_dir):
                success_count += 1
                
        logger.info(f"Batch complete: {success_count}/{len(grh_files)}")
        return success_count == len(grh_files)

def main():
    parser = argparse.ArgumentParser(description="Robust Graph Preprocessing Pipeline")
    parser.add_argument('input', help='Input file or directory')
    parser.add_argument('output', nargs='?', help='Output directory')
    parser.add_argument('--batch', action='store_true', help='Process directory')
    parser.add_argument('--project-root', help='Project root override')
    
    args = parser.parse_args()
    
    try:
        pipeline = GraphPipeline(args.project_root)
        
        if args.batch:
            success = pipeline.process_batch(args.input, args.output)
        else:
            success = pipeline.process_single_graph(args.input, args.output)
            
        sys.exit(0 if success else 1)
        
    except Exception as e:
        logger.exception(f"Fatal error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
