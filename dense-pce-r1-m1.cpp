#include <iostream>
#include <fstream>
#include <unordered_map>
#include <vector>
#include <sstream>
#include <string>
#include <algorithm>
#include <limits>
#include <stack>
#include <numeric>
#include <cmath>
#include <filesystem>
#include <cstring>
#include <cstdint>
#include <climits>
#include "EBBkC/src/edge_oriented.h"

// Globals required by EBBkC in library mode
int K = 0; // target clique size
int L = 1;   // early-termination t-plex parameter (strictest setting)
unsigned long long N = 0ULL; // number of vertices in the graph

// Ablation gating via command-line flags
// Defaults: full FPCE path (order+edge+Turan)
static bool g_enable_order_bound = true;
static bool g_enable_edge_bound  = true;
static bool g_enable_turan       = true;

// counting sort with tracklist
class Graph {
public:
    Graph() = default;
    // Temporary adjacency for building; cleared after CSR finalize
    std::vector<std::vector<int>> adj_map;
    // CSR storage
    std::vector<uint32_t> csr_offsets;   // size n+1 (total nodes + 1)
    std::vector<uint32_t> csr_neighbors; // size m (total edges)

    // CSR view for EBBkC zero-copy entry (from PKT/BSR binaries)
    std::vector<uint32_t> csr_node_off;  // size n+1
    std::vector<int>      csr_edge_dst;  // size m
    uint32_t csr_n = 0;
    uint32_t csr_m = 0;
    // Optional reverse map when CSR is core-shrunk: new-id -> original id
    std::vector<int> csr_rev_map;

    int total_nodes = 0 ;
    // --- Order-bound state ---
    int degeneracy = 0;                 // ξ_G  (max core)
    std::vector<int> core_numbers;      // per-vertex core number (optional for later bounds)

    void add_edge(int u, int v) {
        if (u >= (int)adj_map.size() || v >= (int)adj_map.size()) {
            adj_map.resize(std::max(u, v) + 1);
        }
        adj_map[u].push_back(v);
        adj_map[v].push_back(u);
    }

    void print_graph() const { //prints out each vertex and its list of neighbors
        int n = total_nodes;
        for (int i = 0; i < n; ++i) {
            std::cout << "Node " << i << " : ";
            uint32_t begin = (i < (int)csr_offsets.size()) ? csr_offsets[i] : 0u;
            uint32_t end = (i + 1 < (int)csr_offsets.size()) ? csr_offsets[i + 1] : begin;
            for (uint32_t e = begin; e < end; ++e) {
                std::cout << static_cast<int>(csr_neighbors[e]) << " ";
            }
            std::cout << std::endl;
        }
    }

    std::string graph_file_path;
    void finalize_adjacency() { // converts adj_map into CSR format
        // Sort and unique temporary adjacency
        for (auto& nbrs : adj_map) {
            std::sort(nbrs.begin(), nbrs.end()); // sort each vertex's neighbor list
            nbrs.erase(std::unique(nbrs.begin(), nbrs.end()), nbrs.end()); // remove duplicate edges
        }
        // Build CSR offset array
        int n = static_cast<int>(adj_map.size());
        total_nodes = std::max(total_nodes, n);
        csr_offsets.assign(static_cast<size_t>(total_nodes + 1), 0u);
        for (int u = 0; u < total_nodes; ++u) {
            csr_offsets[u + 1] = csr_offsets[u] + static_cast<uint32_t>(u < n ? adj_map[u].size() : 0u); // cumulative vertex count
        }
        // Populate CSR neighbor array
        uint32_t m = csr_offsets[total_nodes];
        csr_neighbors.clear();
        csr_neighbors.resize(m);
        for (int u = 0; u < total_nodes; ++u) {
            uint32_t off = csr_offsets[u];
            if (u < n) {
                const auto& nbrs = adj_map[u];
                for (size_t k = 0; k < nbrs.size(); ++k) csr_neighbors[off + static_cast<uint32_t>(k)] = static_cast<uint32_t>(nbrs[k]);
            }
        }
        // Free temporary adjacency to save memory
        adj_map.clear();
        adj_map.shrink_to_fit();
    }

    // Matula & Beck k-core in O(n+m) over CSR to get core_numbers and ξ_G
    void compute_core_numbers_from_csr() {
        const int n = total_nodes;
        core_numbers.assign(n, 0);

        if (n == 0) { degeneracy = 0; return; }

        std::vector<int> deg(n);
        int maxdeg = 0;
        for (int u = 0; u < n; ++u) {
            int du = static_cast<int>(csr_offsets[u + 1] - csr_offsets[u]);
            deg[u] = du;
            if (du > maxdeg) maxdeg = du;
        }

        std::vector<int> bin(maxdeg + 1, 0);
        for (int u = 0; u < n; ++u) bin[deg[u]]++;

        int start = 0;
        for (int d = 0; d <= maxdeg; ++d) {
            int num = bin[d];
            bin[d] = start;
            start += num;
        }

        std::vector<int> pos(n), vert(n);
        for (int u = 0; u < n; ++u) {
            pos[u] = bin[deg[u]];
            vert[pos[u]] = u;
            bin[deg[u]]++;
        }
        for (int d = maxdeg; d > 0; --d) bin[d] = bin[d - 1];
        bin[0] = 0;

        int degen = 0;
        for (int i = 0; i < n; ++i) {
            int v = vert[i];
            degen = std::max(degen, deg[v]);
            core_numbers[v] = deg[v];
            // decrease-neighbor step
            for (int64_t it = csr_offsets[v]; it < csr_offsets[v + 1]; ++it) {
                int u = csr_neighbors[it];
                if (deg[u] > deg[v]) {
                    int du = deg[u];
                    int pu = pos[u];
                    int pw = bin[du];
                    int w  = vert[pw];
                    if (u != w) {
                        pos[u] = pw;        vert[pu] = w;
                        pos[w] = pu;        vert[pw] = u;
                    }
                    bin[du]++;             deg[u]--;
                }
            }
        }
        degeneracy = degen;
    }

    // μ(θ, ξG): min of the two corollaries (exactly what we used in mod-edge-order)
    int compute_order_bound(double theta, int xi) const {
        if (theta <= 0.0) return total_nodes;
        const long double th = static_cast<long double>(theta);
        const long double x  = static_cast<long double>(xi);
        const long double eps = 1e-12L;

        // Corollary 1: |S| ≤ ⌊ 2 ξG / θ ⌋
        int b1 = static_cast<int>(std::floor((2.0L * x) / th + eps));

        // Corollary 2: if θ > ξG / (ξG + 1), |S| ≤ ⌊ 1 / (1 − ξG / ((ξG + 1) θ)) ⌋
        int b2 = std::numeric_limits<int>::max();
        long double thresh = x / (x + 1.0L);
        if (th > thresh + eps) {
            long double denom = 1.0L - x / ((x + 1.0L) * th);
            if (denom <= 0.0L) denom = eps;
            b2 = static_cast<int>(std::floor(1.0L / denom + eps));
        }
        return std::min(b1, b2);
    }

    bool has_edge(int u, int v) const {
        if (u < 0 || v < 0 || u >= total_nodes) return false;
        if (csr_offsets.empty()) return false;
        uint32_t begin = csr_offsets[static_cast<size_t>(u)];
        uint32_t end = csr_offsets[static_cast<size_t>(u + 1)];
        uint32_t vv = static_cast<uint32_t>(v); // type-casting
        return std::binary_search(csr_neighbors.begin() + begin, csr_neighbors.begin() + end, vv); // binary search is O(log degree) instead of O(degree).
    }

    void read_graph_from_file(const std::string& filename) {
        graph_file_path = filename;
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Error: Could not open file " << filename << std::endl;
            return;
        }

        int current_vertex = 0;
        std::string line;
        
        while (std::getline(file, line)) {
            if (line == "[EOF]") break;
            total_nodes++;
            if (line.empty()) {
                current_vertex++;
                continue;
            }

            std::stringstream ss(line);
            int neighbor;
            while (ss >> neighbor) {
                add_edge(current_vertex, neighbor);
            }
            current_vertex++;
        }
        // Ensure isolated trailing vertices exist as empty adjacency lists
        if (adj_map.size() < static_cast<size_t>(total_nodes)) {
            adj_map.resize(static_cast<size_t>(total_nodes));
        }
        finalize_adjacency();
    }

    // Load graph from BBkC binaries and build only our own sorted CSR.
    void read_graph_from_bbkcbinaries(const std::string& dir, const std::string& original_graph_path) {
        // file setup
        graph_file_path = original_graph_path;
        std::filesystem::path deg_path = std::filesystem::path(dir) / "b_degree.bin";
        std::filesystem::path adj_path = std::filesystem::path(dir) / "b_adj.bin";
        std::ifstream deg_file(deg_path, std::ios::binary);
        std::ifstream adj_file(adj_path, std::ios::binary);
        if (!deg_file.is_open() || !adj_file.is_open()) {
            std::cerr << "Error: Could not open BBkC binaries in directory " << dir << std::endl;
            return;
        }

        // Read degress file header
        // b_degree.bin: [ui tt][ui n][ui m][ui degree[0..n-1]]
        uint32_t tt=0, n=0, m=0; // tt (type size), n (vertices), m (edges)
        deg_file.read(reinterpret_cast<char*>(&tt), sizeof(uint32_t));
        deg_file.read(reinterpret_cast<char*>(&n), sizeof(uint32_t));
        deg_file.read(reinterpret_cast<char*>(&m), sizeof(uint32_t));
        if (!deg_file.good()) { std::cerr << "Error: Failed reading degree header\n"; return; }
        if (tt != sizeof(uint32_t)) std::cerr << "Warning: ui size header " << tt << "\n";

        // Read degree array
        std::vector<uint32_t> degrees(n);
        if (n>0) {
            deg_file.read(reinterpret_cast<char*>(degrees.data()), static_cast<std::streamsize>(n*sizeof(uint32_t)));
            if (!deg_file.good()) { std::cerr << "Error: Failed reading degree array\n"; return; }
        }

        // Build only our enumerator CSR (sorted rows for binary_search in has_edge())
        total_nodes = static_cast<int>(n);
        csr_offsets.assign(n+1, 0);
        for (uint32_t i=0;i<n;++i) csr_offsets[i+1] = csr_offsets[i] + degrees[i];
        csr_neighbors.resize(csr_offsets[n]);

        // Stream neighbors and sort each row slice
        for (uint32_t u=0; u<n; ++u) {
            const uint32_t du  = degrees[u];
            const uint32_t off = csr_offsets[u];
            const uint32_t end = off + du;
            for (uint32_t k=0; k<du; ++k) {
                uint32_t v;
                adj_file.read(reinterpret_cast<char*>(&v), sizeof(uint32_t));
                if (!adj_file.good()) { std::cerr << "Error: Unexpected EOF in b_adj.bin (u="<<u<<")\n"; return; }
                csr_neighbors[off + k] = v;
            }
            std::sort(csr_neighbors.begin() + off, csr_neighbors.begin() + end);
        }
    }
};

class PseudoCliqueEnumerator {
public:
    PseudoCliqueEnumerator(Graph& graph, float theta = 1.0, int min_size = 1, int max_size = -1, int order_ub = -1) : graph(graph), theta(theta), min_size(min_size) {        
        this->order_ub = order_ub; // μ (order bound); -1 means "disabled"
        total_nodes = graph.total_nodes;
        std::cout << "total nodes: " << total_nodes << " "   << "\n";

        this->max_size = max_size != -1 ? max_size : total_nodes;

        //pseudo_cliques.resize(max_size + 1);
        pseudo_cliques_count.resize(max_size + 1, 0); // count of pseudo-cliques of size i

        // Initialize intrusive buckets,, assign(count, value)
        headP.assign(max_size + 1, -1); // headP[k]: Points to first vertex inside P having degree k
        headNP.assign(max_size + 2, -1); // +1 for targetDeg access, headNP[k]: Points to first candidate vertex outside P (NP) connected to k vertices inside P
        countP.assign(max_size + 1, 0); // countP[k]: Stores the # of vertices inside P having degree k
        countNP.assign(max_size + 2, 0); // countNP[k]: Stores the # of vertices outside P (NP) connected to k vertices inside P
        prevP.assign(graph.total_nodes, -1); // prevP[v]: Points to the previous vertex in the P bucket for vertex v
        nextP.assign(graph.total_nodes, -1); // nextP[v]: Points to the next vertex in the P bucket for vertex v
        prevNP.assign(graph.total_nodes, -1);
        nextNP.assign(graph.total_nodes, -1);
        
        tracks.clear();
        tracks.reserve(graph.total_nodes); // reserve space for all vertices
        tracks.resize(graph.total_nodes); // initialize all tracks to 0
        for (int i = 0; i < graph.total_nodes; ++i) {
            tracks[i].inP = 0; // 0: vertex is not in P
            tracks[i].degP = -1; // -1: vertex has no degree in P
            tracks[i].degNP = 0; // 0: vertex has no degree in NP
            tracks[i].posP = -1; // -1: vertex has no position in P
            tracks[i].posNP = -1; // -1: vertex has no position in NP
            insertNP_front(i, 0); // insert vertex into NP bucket with degree 0
        }

        theta_P = 0; // theta_P: Theta * |P| * (|P| + 1) / 2 - |E_P|
        total_nodes_in_P = 0;
        total_edges_in_P = 0;

        // LHS is fixed for the (ℓ, θ) target
        lhs_required_edges = static_cast<int>(
            std::ceil(theta * (min_size * (min_size - 1) / 2.0))
        );// minimum # of edges P must have to satisfy theta

        // Prepare core histogram only if edge-bound can use it
        max_core_seen = 0;
        if (g_enable_edge_bound) {
            for (int c : graph.core_numbers) max_core_seen = std::max(max_core_seen, c); // max core number in the graph
            core_hist.assign(max_core_seen + 1, 0); // initialize core histogram with max core number + 1
            cur_min_core = INT_MAX; // initialize current minimum core to INT_MAX
        } else {
            core_hist.clear();
            cur_min_core = INT_MAX;
        }
        // Match mod-edge-order / fpce: gate EDGE bound with |P| >= R
        seed_R = static_cast<int>(
            std::ceil(1.0 / (1.0 - theta * (min_size - 1) / static_cast<double>(min_size)))
        );
    }

    // calculates the minimum connection required for any new node to join the current P
    void set_theta_P() { //θ(K) = θ * clq(|K| + 1) - |E[K]|
        theta_P = (theta * total_nodes_in_P * (total_nodes_in_P + 1) / 2.0) - total_edges_in_P;
    }

    // --- Intrusive list helpers ---
    // removes vertex v from the doubly linked list corresponding to degree deg inside P
    inline void unlinkP(int v, int deg) {
        int pv = prevP[v];
        int nv = nextP[v];
        if (pv != -1) nextP[pv] = nv; else headP[deg] = nv;
        if (nv != -1) prevP[nv] = pv;
        prevP[v] = nextP[v] = -1;
        if (deg >= 0 && deg < (int)countP.size()) countP[deg]--;
    }
    // inserts vertex v at the beginning of the doubly linked list for the degree bucket
    inline void insertP_front(int v, int deg) {
        prevP[v] = -1;
        nextP[v] = headP[deg];
        if (headP[deg] != -1) prevP[headP[deg]] = v;
        headP[deg] = v;
        if (deg >= 0 && deg < (int)countP.size()) countP[deg]++;
    }
    inline void unlinkNP(int v, int deg) {
        int pv = prevNP[v];
        int nv = nextNP[v];
        if (pv != -1) nextNP[pv] = nv; else headNP[deg] = nv;
        if (nv != -1) prevNP[nv] = pv;
        prevNP[v] = nextNP[v] = -1;
        if (deg >= 0 && deg < (int)countNP.size()) countNP[deg]--;
    }
    inline void insertNP_front(int v, int deg) {
        prevNP[v] = -1;
        nextNP[v] = headNP[deg];
        if (headNP[deg] != -1) prevNP[headNP[deg]] = v;
        headNP[deg] = v;
        if (deg >= 0 && deg < (int)countNP.size()) countNP[deg]++;
    }

    void print_all() const {
        std::cout << "P: ";
        for (int d = 0; d <= total_nodes_in_P; ++d) {
            for (int v = headP[d]; v != -1; v = nextP[v]) std::cout << v << " ";
            std::cout << "| ";
        }
        std::cout << "\nNP: ";
        for (int d = 0; d <= total_nodes_in_P + 1 && d < (int)headNP.size(); ++d) {
            for (int v = headNP[d]; v != -1; v = nextNP[v]) std::cout << v << " ";
            std::cout << "| ";
        }
        // std::cout << "\ntracks: ";
        // for (const auto& track : tracks) {
        //     std::cout << "\n\n" ;
        //     int key = track.first;
        //     auto& value = track.second;

        //     std::cout << key << ": {is_inside_P: " << value[0]
        //               << ", degree_in_P: " << value[1]
        //               << ", degree_in_NP: " << value[2]
        //               << ", position_inside_P: " << value[3]
        //               << ", position_inside_NP: " << value[4] << "} ";
        // }
        std::cout << "\ntheta_P: " << theta_P << " (" << total_nodes_in_P << ", " << total_edges_in_P << ")\n";
    }

    void add_to_inside_P(int v);
    void remove_from_inside_P(int v);
    void iter(int v);
    // Internal: insert v into P WITHOUT recursing (used by Turán seed build)
    void add_vertex_internal(int v);

    std::vector<int> get_pseudo_cliques_count(){
        return pseudo_cliques_count;
    };

    int get_iter_count(){
        return iter_count;
    };

    unsigned long long get_numcalls_saved_by_edge_bound() const {
        return numcalls_saved_by_edge_bound;
    }

    // Enumerate starting from a Turán seed (R-clique) provided by EBBkC, takes a "seed" clique (found by the external library) and "teleports" the algorithm to that state, effectively skipping the early levels of the recursion tree to prune the search space.
    void enumerate_with_turan(const std::vector<int>& turan_seed) {
        if (turan_seed.empty()) return;

        // A) Deterministic order: sort by (deg, id) to match mod-edge-order seed canonicalization
        auto sorted_seed = turan_seed;
        std::sort(sorted_seed.begin(), sorted_seed.end(),
                [&](int a, int b) {
                    // degree proxy from CSR slice lengths, sort by degree first, then by id to ensure deterministic ordering
                    const uint32_t da = graph.csr_offsets[a+1] - graph.csr_offsets[a]; // cumulative count of next node - current node
                    const uint32_t db = graph.csr_offsets[b+1] - graph.csr_offsets[b];
                    return std::make_pair(da, a) < std::make_pair(db, b);
                });

        // B) Build the whole seed into P with NO recursion
        for (int v : sorted_seed) add_vertex_internal(v);

        // C) Choose v* from the current δ(P) bucket (minimum-id in that bucket)
        int deltaP = 0;
        for (int d = 0; d <= total_nodes_in_P; ++d) {
            if (d >= 0 && d < (int)headP.size() && headP[d] != -1) { deltaP = d; break; } // starting scanning from 0, first non-empty bucket corrsponds to min-degree
        }
        int v_star = -1;
        if (deltaP >= 0 && deltaP < (int)headP.size()) {
            int best = -1;
            for (int x = headP[deltaP]; x != -1; x = nextP[x]) {
                if (best == -1 || x < best) best = x;
            }
            v_star = (best != -1 ? best : sorted_seed.front());
        } else {
            v_star = sorted_seed.front();
        }

        // D) Recurse from v*
        iter(v_star);

        // E) Roll back the seed in reverse order
        for (auto it = sorted_seed.rbegin(); it != sorted_seed.rend(); ++it) {
            remove_from_inside_P(*it);
        }
    }

private:
    struct Track {
        uint8_t  inP;     // 0/1 flag
        uint16_t degP;    // ∈ [0, μ]
        uint16_t degNP;   // ∈ [0, μ]
        int32_t  posP;    // index in inside_P[degP]
        int32_t  posNP;    // index in neighbors_and_P[degNP]
    };
    // Collect current pseudo-clique P as a sorted vector (via intrusive lists)
    std::vector<int> collect_current_pseudo_clique() const {
        std::vector<int> current_p;
        current_p.reserve(total_nodes_in_P);
        int maxDeg = static_cast<int>(headP.size()) - 1;
        if (maxDeg < 0) maxDeg = 0;
        for (int d = 0; d <= maxDeg; ++d) {
            for (int v = headP[d]; v != -1; v = nextP[v]) {
                current_p.push_back(v);
            }
        }
        std::sort(current_p.begin(), current_p.end());
        return current_p;
    }

    // No vertex outside P has degree into P >= ceil(theta_P)
    bool is_current_maximal() const {
        int min_add_deg = static_cast<int>(std::ceil(theta_P));
        if (min_add_deg < 0) min_add_deg = 0;
        // Compare histogram counts: all vertices vs vertices in P (using counts)
        int maxD = std::min(static_cast<int>(countNP.size()) - 1, total_nodes_in_P); // sets uppr limit for the loop
        for (int d = min_add_deg; d <= maxD; ++d) {
            size_t occ_num = (d >= 0 && d < static_cast<int>(countNP.size())) ? static_cast<size_t>(countNP[d]) : 0; // countNP[d]: total # of vertices in entire graph that connect to d vertices inside P
            size_t in_p_num = (d >= 0 && d < static_cast<int>(countP.size())) ? static_cast<size_t>(countP[d]) : 0; // countP[d]: total # of vertices inside P having degree d
            if (occ_num > in_p_num) return false; // some outside vertex can extend P
        }
        return true;
    }

    Graph& graph;
    float theta;
    int min_size;
    int max_size;
    int total_nodes;
    float theta_P;
    int total_nodes_in_P;
    int total_edges_in_P;
    int iter_count = 0;
    std::stack<int> children;
    // Intrusive bucketed lists for degrees 0..max_size
    std::vector<int32_t> headP, headNP;   // heads per degree (-1 if empty)
    std::vector<int32_t> prevP, nextP;    // per-vertex links in P buckets
    std::vector<int32_t> prevNP, nextNP;  // per-vertex links in NP buckets
    std::vector<int32_t> countP, countNP; // counts per degree
    std::vector<std::vector<int> > pseudo_cliques;
    std::vector<int> pseudo_cliques_count;
    std::vector<Track> tracks;
    int order_ub = -1;   // μ (order bound); -1 means "disabled"

    // --- Edge-bound state ---
    int lhs_required_edges = 0;             // ceil(theta * l*(l-1)/2)   [fixed for a run]
    std::vector<int> core_hist;             // histogram of core numbers among vertices in P
    int cur_min_core = INT_MAX;             // c(P) = min core within P
    int max_core_seen = 0;                  // histogram capacity

    // Optional accounting (like fpce.c/mod-edge-order)
    unsigned long long numcalls_saved_by_edge_bound = 0ULL;

    // Turán seed size used to gate EDGE bound like mod-edge-order
    int seed_R = 1;

    inline void report_if_new(const std::vector<int>& clique_sorted) {
        if (clique_sorted.size() < static_cast<size_t>(min_size)) return;
        if (clique_sorted.size() < pseudo_cliques_count.size())
            pseudo_cliques_count[clique_sorted.size()] += 1;
    }
    
    // Small helpers that work with intrusive lists
    inline int min_degree_in_P_fast() const {
        if (total_nodes_in_P == 0) return 0;
        for (int d = 0; d <= total_nodes_in_P; ++d) {
            if (d >= 0 && d < (int)headP.size() && headP[d] != -1) return d;
        }
        return 0;
    }
    inline int current_min_core_in_P_fast() const {
        return (cur_min_core == INT_MAX) ? 0 : cur_min_core;
    }

    bool prune_by_edge_bound() const {
        const int S = min_size;                 // ℓ
        const int P = total_nodes_in_P;         // |P|
        const int t = S - P;                    // vertices still needed
    
        // Only prune when 1 ≤ |P| < ℓ  (same gating as mod-edge-order / fpce.c)
        if (P < 1 || t <= 0) return false;
    
        const int delta   = min_degree_in_P_fast();         // δ(P)
        const int cP      = current_min_core_in_P_fast();   // c(P)
        const double sumE = static_cast<double>(total_edges_in_P);
        const int LHS     = lhs_required_edges;
    
        // Lemma 9 shape (fpce.c): g = c(P) - δ(P) - 1
        const int g = cP - delta - 1;
    
        if (g >= 0) {
            // RHS1 = |E[P]| + min(c(P), δ(P)+t)*(t + δ + 0.5) - 0.5*(δ+1)*(2δ+1)
            const int max_value = std::min(cP, delta + t);
            const double rhs1 = std::floor(
                sumE + max_value * (t + delta + 0.5)
                    - 0.5 * (delta + 1) * (2 * delta + 1)
            );
            if (LHS > rhs1) return true;   // EDGE_BOUND1
        } else {
            // RHS2 = |E[P]| + t*δ
            const double rhs2 = std::floor(sumE + t * static_cast<double>(delta));
            if (LHS > rhs2) return true;   // EDGE_BOUND2
        }
    
        // Lemma 8 fallback: RHS3 = |E[P]| + t*(δ + (t+1)/2)
        const double rhs3 = std::floor(sumE + t * (delta + (t + 1) * 0.5));
        if (LHS > rhs3) return true;       // EDGE_BOUND3
    
        return false;
    }
    
}; 

void PseudoCliqueEnumerator::add_vertex_internal(int v) {
    if (total_nodes_in_P > max_size) return;
    // move v into P
    tracks[v].inP  = 1;
    tracks[v].degP = static_cast<uint16_t>(tracks[v].degNP);
    insertP_front(v, tracks[v].degP);

    // core histogram (for edge-bound)
    int cv = (v >= 0 && v < (int)graph.core_numbers.size()) ? graph.core_numbers[v] : 0;
    if ((int)core_hist.size() <= cv) core_hist.resize(cv + 1, 0);
    core_hist[cv] += 1;
    if (cv < cur_min_core) cur_min_core = cv;
    tracks[v].posP = -1;

    // update neighbors’ P/NP degrees and |E[P]|
    uint32_t begin = (v >= 0 && v + 1 < (int)graph.csr_offsets.size())
                     ? graph.csr_offsets[(size_t)v] : 0u; // neighbor list for vertex v starts
    uint32_t end   = (v >= 0 && v + 1 < (int)graph.csr_offsets.size())
                     ? graph.csr_offsets[(size_t)v + 1] : begin; // neighbor list for vertex v ends
    for (uint32_t ei = begin; ei < end; ++ei) {
        int u = (int)graph.csr_neighbors[ei];
        // Update neighbors inside P
        if (tracks[u].inP) {
            total_edges_in_P += 1;
            int d = static_cast<int>(tracks[u].degP);
            unlinkP(u, d);
            insertP_front(u, d + 1);
            tracks[u].degP = static_cast<uint16_t>(d + 1);
        }
        // Update neighbors outside P
        int dnp = static_cast<int>(tracks[u].degNP);
        unlinkNP(u, dnp);
        insertNP_front(u, dnp + 1);
        tracks[u].degNP = static_cast<uint16_t>(dnp + 1);
    }

    total_nodes_in_P += 1;
    set_theta_P();
}

void PseudoCliqueEnumerator::add_to_inside_P(int v) {
    add_vertex_internal(v);   // state updates only
    iter(v);                  // now recurse once (normal growth path)
}

void PseudoCliqueEnumerator::remove_from_inside_P(int v) {
    // Remove v from P
    tracks[v].inP = 0;
    int degree = tracks[v].degP;
    int pos = tracks[v].posP;
    if (degree >= 0) unlinkP(v, degree);

    // Iterate Neighbors of v
    {
        uint32_t begin = (v >= 0 && v + 1 < (int)graph.csr_offsets.size()) ? graph.csr_offsets[static_cast<size_t>(v)] : 0u; // neighbor list for vertex v starts
        uint32_t end = (v >= 0 && v + 1 < (int)graph.csr_offsets.size()) ? graph.csr_offsets[static_cast<size_t>(v + 1)] : begin; // neighbor list for vertex v ends
        for (uint32_t ei = begin; ei < end; ++ei) {
            int adj_u = static_cast<int>(graph.csr_neighbors[ei]); // neighbor of v
            if (tracks[adj_u].inP) {
                total_edges_in_P--; // update |E[P]|
                int degree = static_cast<int>(tracks[adj_u].degP); // degree of neighbor
                unlinkP(adj_u, degree); // remove neighbor from P
                if (degree > 0) {
                    insertP_front(adj_u, degree - 1);
                    tracks[adj_u].degP = static_cast<uint16_t>(degree - 1);
                } else {
                    // Sanity: degree should not be negative; clamp to 0
                    tracks[adj_u].degP = 0;
                    tracks[adj_u].posP = -1;
                }
            }
            // Downgrade Neighbors Outside P
            int degree_np = static_cast<int>(tracks[adj_u].degNP);
            int pos_np = tracks[adj_u].posNP;
            unlinkNP(adj_u, degree_np);
            if (degree_np > 0) {
                insertNP_front(adj_u, degree_np - 1);
                tracks[adj_u].degNP = static_cast<uint16_t>(degree_np - 1);
            } else {
                // Sanity: clamp at 0
                tracks[adj_u].degNP = 0;
                tracks[adj_u].posNP = -1;
            }
        }
    }

    // Update Edge Bound Statistics and c(P) on removal
    int cv = (v >= 0 && v < (int)graph.core_numbers.size()) ? graph.core_numbers[v] : 0;
    if (cv >= 0 && cv < (int)core_hist.size()) {
        if (core_hist[cv] > 0) core_hist[cv] -= 1;
        if (core_hist[cv] == 0 && cv == cur_min_core) {
            // advance to next non-empty bucket or reset
            int nextc = cv;
            while (nextc < (int)core_hist.size() && core_hist[nextc] == 0) ++nextc;
            cur_min_core = (nextc < (int)core_hist.size()) ? nextc : INT_MAX;
        }
    }

    total_nodes_in_P--;
    set_theta_P();
}

void PseudoCliqueEnumerator::iter(int v) {
    iter_count++;

    // EDGE bound: apply only when 1 ≤ |P| < ℓ (parity-correct ob-eb)
    if (g_enable_edge_bound && total_nodes_in_P >= 1 && total_nodes_in_P < min_size) {
        if (prune_by_edge_bound()) {
            numcalls_saved_by_edge_bound++;
            return;
        }
    }

    // Integer add-degree threshold must be ceil(theta_P)
    const int min_add_deg = static_cast<int>(std::ceil(theta_P));

    // Always partition with respect to δ(P)
    const int deltaP = min_degree_in_P_fast(); // min-deg currently inside P

    std::vector<int> children_vec;

    // ---- Type 1: deg == 0..δ(P)-1 ----
    for (int deg = std::max(0, min_add_deg); deg < deltaP; ++deg) {
        for (int u = headNP[deg]; u != -1; u = nextNP[u]) {
            if (tracks[u].inP || tracks[u].degNP < min_add_deg) continue;
            children_vec.push_back(u);
        }
    }

    // ---- Type 2: deg == δ(P) ----
    for (int u = headNP[deltaP]; u != -1; u = nextNP[u]) {
        if (tracks[u].inP || tracks[u].degNP < min_add_deg) continue;
        if (u < v) {
            children_vec.push_back(u);
        } else if (graph.has_edge(v, u)) { // u must be connected to v AND to every other node x in the minimum-degree bucket that is smaller than u
            bool valid = true;
            // must be adjacent to all x in P[δ(P)] with x < u
            for (int x = headP[deltaP]; x != -1; x = nextP[x]) {
                if (x < u && !graph.has_edge(u, x)) { valid = false; break; }
            }
            if (valid) children_vec.push_back(u);
        }
    }

    // ---- Type 3: deg == δ(P)+1 ----
    {
        const int targetDeg = deltaP + 1;
        if (targetDeg < static_cast<int>(headNP.size())) {
            for (int u = headNP[targetDeg]; u != -1; u = nextNP[u]) {
                if (tracks[u].inP || tracks[u].degNP < min_add_deg) continue;
                if (!graph.has_edge(v, u) || !(v > u)) continue; // allow this child u if its ID is smaller than the current node v

                bool valid = true;
                // (1) u must be adjacent to all x in P[δ(P)]
                for (int x = headP[deltaP]; x != -1; x = nextP[x]) {
                    if (!graph.has_edge(u, x)) { valid = false; break; } 
                }
                // (2) u must also connect to all existing nodes in the delta_P + 1 bucket that have a smaller ID than itself
                if (valid && targetDeg < static_cast<int>(headP.size())) {
                    for (int x = headP[targetDeg]; x != -1; x = nextP[x]) {
                        if (x < u && !graph.has_edge(u, x)) { valid = false; break; } 
                    }
                }
                if (valid) children_vec.push_back(u);
            }
        }
    }

    // --- report maximal (existing block retained) ---
    if (total_nodes_in_P >= min_size && is_current_maximal()) {
        auto cur = collect_current_pseudo_clique();
        std::sort(cur.begin(), cur.end());
        report_if_new(cur);
    }

    // Local order-bound stop at μ: do not expand past μ
    if (g_enable_order_bound && order_ub > 0 && total_nodes_in_P == order_ub) return;

    // Deterministic expansion order
    std::sort(children_vec.begin(), children_vec.end());
    for (int u : children_vec) {
        add_to_inside_P(u);
        remove_from_inside_P(u);
    }
}


int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <filename> [--theta <theta>] [--minimum <min>] [--maximum <max>]" << std::endl;
        return 1;
    }

    // 1. Re‑assemble the (possibly space‑containing) filename.
    std::vector<std::string> filename_parts;
    int arg_idx = 1;
    for (; arg_idx < argc; ++arg_idx) {
        if (std::strncmp(argv[arg_idx], "--", 2) == 0) break; // first option detected
        filename_parts.emplace_back(argv[arg_idx]);
    }

    if (filename_parts.empty()) {
        std::cerr << "Error: filename not provided." << std::endl;
        return 1;
    }

    std::string filename;
    for (size_t k = 0; k < filename_parts.size(); ++k) {
        if (k) filename += ' ';
        filename += filename_parts[k];
    }

    // 2. Parse optional arguments that follow the filename.
    double theta = 1.0;
    int minimum = 1;
    int maximum = std::numeric_limits<int>::max();
    

    for (; arg_idx < argc; ++arg_idx) {
        std::string arg = argv[arg_idx];
        if (arg == "--theta" && arg_idx + 1 < argc) {
            theta = std::stod(argv[++arg_idx]);
        } else if (arg == "--minimum" && arg_idx + 1 < argc) {
            minimum = std::stoi(argv[++arg_idx]);
        } else if (arg == "--mode" && arg_idx + 1 < argc) {
            int mode = std::stoi(argv[++arg_idx]);
            if (mode == 1) { g_enable_order_bound = false; g_enable_edge_bound = false; g_enable_turan = false; }
            else if (mode == 2) { g_enable_order_bound = true; g_enable_edge_bound = false; g_enable_turan = false; }
            else if (mode == 3) { g_enable_order_bound = true; g_enable_edge_bound = true; g_enable_turan = false; }
            else { g_enable_order_bound = true; g_enable_edge_bound = true; g_enable_turan = true; }
        } else if (arg == "--no-order") {
            g_enable_order_bound = false;
        } else if (arg == "--no-edge") {
            g_enable_edge_bound = false;
        } else if (arg == "--no-turan") {
            g_enable_turan = false;
        } else if (arg == "--order") {
            g_enable_order_bound = true;
        } else if (arg == "--edge") {
            g_enable_edge_bound = true;
        } else if (arg == "--turan") {
            g_enable_turan = true;
        } else {
            std::cerr << "Unknown option: " << arg << std::endl;
            return 1;
        }
    }
    if (minimum <= 0) {
        std::cerr << "Error: --minimum must be positive." << std::endl;
        return 1;
    }
    if (theta <= 0.0 || theta > 1.0) {
        std::cerr << "Error: --theta must be in the open interval (0,1)." << std::endl;
        return 1;
    }

    int R = std::ceil(1.0 / (1.0 - theta * (minimum - 1) / static_cast<double>(minimum)));
    std::cout << "Computed R: " << R << std::endl;
    std::cout << "Gates => order:" << (g_enable_order_bound?"on":"off")
              << " edge:" << (g_enable_edge_bound?"on":"off")
              << " turan:" << (g_enable_turan?"on":"off") << std::endl;

    // 4. Derive directory path in a cross‑platform way.
    std::filesystem::path p(filename);
    std::string dir_path = p.parent_path().empty() ? "." : p.parent_path().string();
    std::cout << "Directory path: " << dir_path << std::endl;

    Graph graph;
    // Try to load BBkC binaries from the same directory as the provided .grh path
    {
        std::filesystem::path bdeg = std::filesystem::path(dir_path) / "b_degree.bin";
        std::filesystem::path badj = std::filesystem::path(dir_path) / "b_adj.bin";
        if (std::filesystem::exists(bdeg) && std::filesystem::exists(badj)) {
            std::cout << "Detected BBkC binaries in: " << dir_path << ", loading graph from binaries..." << std::endl;
            graph.read_graph_from_bbkcbinaries(dir_path, filename);
        } else {
            graph.read_graph_from_file(filename);
        }
    }
    
    // --- Core shrink mirroring FPCE roots: compute cores on our CSR,
    //     then build EBBkC CSR as the induced (R-1)-core subgraph (only if Turán enabled).
    graph.compute_core_numbers_from_csr();  // fills graph.core_numbers and graph.degeneracy
    std::cout << "degeneracy: " << graph.degeneracy << std::endl;
    if (g_enable_turan && R >= 3) {
        const int n = graph.total_nodes;
        std::vector<int> keep(n, 0), newid(n, -1), rev;
        int kept = 0;
        for (int u = 0; u < n; ++u) {
            if (u < (int)graph.core_numbers.size() && graph.core_numbers[u] >= (R - 1)) {
                keep[u] = 1;
                newid[u] = kept++;
            }
        }
        // Build induced CSR into EBBkC arrays
        std::vector<uint32_t> r_off((size_t)kept + 1u, 0u);
        std::vector<int>      r_dst;
        r_dst.reserve(graph.csr_neighbors.size());
        rev.assign(kept, -1);
        for (int u = 0; u < n; ++u) {
            if (!keep[u]) continue;
            rev[newid[u]] = u; // original vertex u gets newid[u] as its new index
            uint32_t begin = graph.csr_offsets[(size_t)u];
            uint32_t end   = graph.csr_offsets[(size_t)u + 1];
            for (uint32_t e = begin; e < end; ++e) {
                int v = (int)graph.csr_neighbors[e];
                if (v >= 0 && v < n && keep[v]) r_dst.push_back(newid[v]);
            }
            r_off[(size_t)(newid[u] + 1)] = (uint32_t)r_dst.size();
        }
        graph.csr_node_off.swap(r_off);
        graph.csr_edge_dst.swap(r_dst);
        graph.csr_n = (uint32_t)kept;
        graph.csr_m = (uint32_t)graph.csr_edge_dst.size();
        graph.csr_rev_map.swap(rev);
        std::cout << "EBBkC CSR restricted to (R-1)-core: n=" << graph.csr_n
                  << " m=" << graph.csr_m << std::endl;
    }
    
    // If maximum is not specified, use total vertices (existing code below also does this)
    if (maximum == std::numeric_limits<int>::max()) {
        maximum = graph.total_nodes;
    }

    // 1) Compute ξ_G from CSR, then μ(θ, ξ_G)
    if (g_enable_order_bound || g_enable_edge_bound) {
        graph.compute_core_numbers_from_csr();
    }
    int mu = 0;
    if (g_enable_order_bound) {
        mu = graph.compute_order_bound(theta, graph.degeneracy);
        std::cout << "Order bound μ = " << mu << " (degeneracy ξ_G = " << graph.degeneracy << ")\n";
        if (minimum > mu) {
            std::cout << "Pruning by Order Bound: no (ℓ,θ)-pseudo-clique can exist since ℓ (" 
                    << minimum << ") > μ (" << mu << ")\n";
            return 0;
        }
        if (maximum > mu) maximum = mu;
    } else {
        std::cout << "Order bound disabled" << std::endl;
    }

    PseudoCliqueEnumerator PC(graph, theta, minimum, maximum,
#if ENABLE_ORDER_BOUND
        mu
#else
        -1
#endif
    );

    // Seeding strategy: Turán via EBBkC when enabled and R >= 3; otherwise node-by-node
    if (g_enable_turan && R >= 3 && graph.csr_n > 0) {
        const int n = static_cast<int>(graph.csr_n);
        const int m = static_cast<int>(graph.csr_m);
        EBBkC_t::list_k_clique_mem_stream_from_csr(
            n,
            m,
            graph.csr_node_off.data(),
            graph.csr_edge_dst.data(),
            R, /*threads=*/2,
            [&](const std::vector<int>& r_clique){
                if (!graph.csr_rev_map.empty()) {
                    std::vector<int> mapped;
                    mapped.reserve(r_clique.size());
                    for (int id : r_clique) mapped.push_back(graph.csr_rev_map[id]);
                    PC.enumerate_with_turan(mapped);
                } else {
                    PC.enumerate_with_turan(r_clique);
                }
            }
        );
        graph.csr_node_off.clear(); graph.csr_node_off.shrink_to_fit();
        graph.csr_edge_dst.clear(); graph.csr_edge_dst.shrink_to_fit();
    } else {
        for (int node = 0; node < graph.total_nodes; ++node) {
            PC.add_to_inside_P(node);
            PC.remove_from_inside_P(node);
        }
    }

    std::vector<int> pseudo_clique_counts = PC.get_pseudo_cliques_count();

    // Print the clique sizes
    std::cout << "Maximal pseudo-clique counts:" << std::endl;
    bool found_any = false;
    long long total_found = 0;
    for (size_t sz = 0; sz < pseudo_clique_counts.size(); ++sz) {
        if (pseudo_clique_counts[sz] > 0) {
            std::cout << "Size " << sz << ": " << pseudo_clique_counts[sz] << "\n";
            found_any = true;
            total_found += pseudo_clique_counts[sz];
        }
    }
     
    if (!found_any) {
        std::cout << "(No pseudo-cliques found meeting the criteria)" << std::endl;
    }
 
    std::cout << "Total Maximal Pseudo-Cliques: " << total_found << std::endl;
    std::cout << "\nTotal Iterations: " << PC.get_iter_count() << std::endl;
    std::cout << "Edge-bound prunes saved: " << PC.get_numcalls_saved_by_edge_bound() << std::endl;

    return 0;
}