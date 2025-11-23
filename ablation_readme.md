# Ablation Runs: Gated and Non‑Gated PCE Variants

This note explains how to build and run the ablation‑friendly binary (`dense-pce-ab-gated`) and the non‑gated baselines on your datasets.

## Build

Use the project build script (links EBBkC automatically):

```bash
wsl bash build.sh dense-pce-ab-gated.cpp
wsl bash build.sh dense-pce-ab.cpp
wsl bash build.sh dense-pce-mod-edge-order.cpp
```

Outputs appear under `build_integrated/`:
- `dense-pce-ab-gated-integrated`
- `dense-pce-ab-integrated`
- `dense-pce-mod-edge-order-integrated`

## Inputs

Executables accept either:
- A `.grh` text graph (one line per vertex: neighbors separated by spaces; trailing `[EOF]` optional), or
- A directory containing `b_degree.bin` and `b_adj.bin` (BBkC PKT/BSR binaries). If these exist alongside the `.grh`, they are used for zero‑copy Turán seeding.

Most commands below use `testGraph/fpce_graph/fpce_graph.grh` as an example.

## Gated variant (single binary, four modes)

`dense-pce-ab-gated-integrated` supports a single `--mode` switch and fine‑grained toggles. Defaults to mode 4.

Modes:
- Mode 1: no pruning (no order‑bound, no edge‑bound, no Turán); node‑by‑node enumeration
- Mode 2: order‑bound only (no edge‑bound, no Turán)
- Mode 3: order‑bound + edge‑bound (no Turán)
- Mode 4: full FPCE path (order + edge + Turán seeding)

Fine‑grained flags (optional, override mode):
- `--order | --no-order`
- `--edge | --no-edge`
- `--turan | --no-turan`

Examples:

```bash
# Full FPCE (default == --mode 4)
./build_integrated/dense-pce-ab-gated-integrated testGraph/fpce_graph/fpce_graph.grh --minimum 3 --theta 0.6

# Explicit modes
wsl ./build_integrated/dense-pce-ab-gated-integrated ... --minimum 3 --theta 0.6 --mode 1   # PCE-like (no pruning)
wsl ./build_integrated/dense-pce-ab-gated-integrated ... --minimum 3 --theta 0.6 --mode 2   # order only
wsl ./build_integrated/dense-pce-ab-gated-integrated ... --minimum 3 --theta 0.6 --mode 3   # order+edge
wsl ./build_integrated/dense-pce-ab-gated-integrated ... --minimum 3 --theta 0.6 --mode 4   # order+edge+Turán

# Overrides (examples)
./build_integrated/dense-pce-ab-gated-integrated ... --mode 1 --order        # enable only order-bound
./build_integrated/dense-pce-ab-gated-integrated ... --mode 4 --no-turan     # full minus Turán seeding
```

Notes:
- When Turán is off, the EBBkC CSR buffers are freed immediately after load to reduce memory.
- When both bounds are off, core numbers and edge‑bound histogram allocations are skipped.

## Non‑gated variants (for baseline comparisons)

- `dense-pce-ab-integrated`: AB variant aligned with mod‑edge‑order parity (ceil(theta_P), δ(P) partition), Turán on, order + edge bounds on.
- `dense-pce-mod-edge-order-integrated`: reference mod‑edge‑order implementation.

Examples:

```bash
./build_integrated/dense-pce-ab-integrated testGraph/fpce_graph/fpce_graph.grh --minimum 3 --theta 0.6
./build_integrated/dense-pce-mod-edge-order-integrated testGraph/fpce_graph/fpce_graph.grh --minimum 3 --theta 0.6
```

## Batch runs with the helper script

`run_single_param.sh` runs one graph per subdirectory and records per‑graph logs under `logs/run_<timestamp>/`.

Usage:
```bash
bash run_single_param.sh -d <graph_root_dir> -l <min_size> -t <theta> -m <mode> [-j N] [-x <exec_path>]
```
- `-d DIR` graph root directory (each subfolder should contain one `.grh` and optionally BBkC binaries)
- `-l SIZE` minimum size ℓ (default 10)
- `-t THETA` density threshold θ (default 0.9)
- `-m MODE` ablation mode (1–4, default 4)
- `-j N` run up to N subdirectories concurrently (default 1)
- `-x PATH` override the executable (defaults to `dense-pce-ab-gated-integrated` if present, otherwise AB legacy)

Examples:
```bash
# Full FPCE across all datasets in a folder
bash run_single_param.sh -d testGraph -l 3 -t 0.6 -m 4

# Order-only ablation across a folder using the gated binary
bash run_single_param.sh -d testGraph -l 5 -t 0.8 -m 2

# Run with a specific executable (legacy AB) in mode 3 semantics (ignored by legacy)
bash run_single_param.sh -d testGraph -l 3 -t 0.7 -m 3 -x ./build_integrated/dense-pce-ab-integrated
```

## Output

Each run prints per‑size counts, the total number of maximal pseudo‑cliques, total iterations, and (when applicable) the number of subtrees saved by the edge bound. Example tail:

```
Maximal pseudo-clique counts:
Size 3: 7
Size 5: 7
Total Maximal Pseudo-Cliques: 14

Total Iterations: 60
Edge-bound prunes saved: 0
```

## Tips
- If a directory contains `b_degree.bin` and `b_adj.bin`, EBBkC seeds are streamed directly from those CSR arrays (mode 4). Otherwise, enumeration falls back to the node‑by‑node driver.
- The `--mode` toggle is output‑neutral; it only controls pruning and seeding to explore iteration/memory tradeoffs.
