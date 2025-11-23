#ifndef DESCOL_EBBKC_H
#define DESCOL_EBBKC_H
#include "def.h"
#include <vector>
#include <functional>

using namespace std;

// —— new: store each found clique ——
// each inner vector holds the global vertex IDs of one clique
// extern unsigned long long N;
extern std::vector<std::vector<int>> cliques_vec;

class EBBkC_Graph_t {

public:
    // for both EBBkC & EBBkC+ framework
    int v_size = 0;
    int e_size = 0;
    int truss_num = 0;
    bool is_sub = false;
    HashMap_t edge2id;
    Edge_t* edges = nullptr;

    int* sub_e_size = nullptr;    // number of edges in the sub-branch
    int* sub_v_size = nullptr;    // number of nodes in the sub-branch
    int** sub_v = nullptr;        // nodes in the sub-branch
    int** sub_e = nullptr;        // edges in the sub-branch

    int** T = nullptr;            // T[e]: a set of nodes
    int* T_size = nullptr;        // T_size[e]: the size of T[e]
    int** C = nullptr;            // C[e]: a set of edges
    int* C_size = nullptr;        // C_size[e]: the size of C[e]

    vector<int> new2old;
    vector<int> rank;

    // for EBBkC framework only
    int* v_lab = nullptr;
    int* e_lab = nullptr;
    int** out_v_size = nullptr;
    int** out_e_size = nullptr;

    // for EBBkC+ framework only
    int* col = nullptr;           // col[u]: the color of vertex 'u'
    int* lab = nullptr;           // lab[u]: the level of vertex 'u'
    int** G_deg = nullptr;        // G_deg[l][u]: the number of neighbors of vertex 'u' at level 'l'
    int** G_adj = nullptr;        // G_adj[l][u]: the set of neighbors of vertex 'u' at level 'l'
    int** DAG_deg = nullptr;      // DAG_deg[l][u]: the number of out-neighbors of vertex 'u' at level 'l'
    int** DAG_adj = nullptr;      // DAG_adj[l][u]: the set of out-neighbors of vertex 'u' at level 'l'
    bool** used = nullptr;

    // for EBBkC++ (early-termination)
    int P_size = 0;
    int P_act = 0;
    int F_size = 0;
    int *P = nullptr;
    int *F = nullptr;
    int *lack_size = nullptr;
    int **lack = nullptr;
    int *lev = nullptr;
    int *loc = nullptr;

    std::vector<std::pair<int,int>> branch_ends; // size K+1 
    // In EBBkC_Graph_t declaration
    std::vector<int> current_clique;  // Track vertices in current clique
    // Streaming sink for found cliques (global vertex IDs). If set, results are streamed.
    std::function<void(const std::vector<int>&)> clique_sink;

    EBBkC_Graph_t();
    ~EBBkC_Graph_t();

//    void read_edges_from_file(const char* file_name);
    void truss_decompose(const char* w_file_name);
//    void read_ordered_edges_from_file(const char* file_name);
    void build(bool sub);

//    void EBBkC(int l, unsigned long long *cliques);
//    void EBBkC_plus(int, unsigned long long *cliques);
//    void EBBkC_plus_plus(int l, unsigned long long *cliques);

    void branch(int e, EBBkC_Graph_t* g);
    void EBBkC_plus_plus(int l, unsigned long long *cliques);

    bool can_terminate(int l, unsigned long long* cliques);
    void list_in_plex(int start, int p, int q, unsigned long long* cliques);
    void EBBkC_Comb_list(int *list, int  list_size, int  start, int  picked, int  k, unsigned long long *cliques);

    //file writing function for dense-pce integration
    void finalize_clique_output(const char* output_path);
};


class EBBkC_t {
public:
    // static double truss_order(const char* r_file_name, const char* w_file_name);
    static double list_k_clique(const char* file_name);
    static double list_k_clique_parallel(const char* file_name);
    // In-memory API: run EBBkC on directory `dir` with parameters (k,l),
    // and return the found k-cliques in `out_cliques` without any file I/O.
    static double list_k_clique_mem(const char* dir, int k, int l, std::vector<std::vector<int>>& out_cliques);
    // Streaming API: run EBBkC and invoke `sink` for each found k-clique; retains O(1) memory.
    static double list_k_clique_mem_stream(const char* dir, int k, int l, const std::function<void(const std::vector<int>&)>& sink);
    // CSR streaming API: avoid disk load by passing CSR arrays directly (borrowed, not freed)
    static double list_k_clique_mem_stream_from_csr(
        uint32_t n, uint32_t m,
        const uint32_t* node_off,
        const int* edge_dst,
        int k, int l,
        const std::function<void(const std::vector<int>&)>& sink);
};


#endif //DESCOL_EBBKC_H
