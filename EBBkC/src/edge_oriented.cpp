#include <set>
#include <algorithm>
#include <unordered_map>
#include <omp.h>
#include <cassert>
#include <chrono>
#include <vector>
#include "set_operation.h"
#include "edge_oriented.h"
#include "truss/util/graph/graph.h"
#include "truss/util/log/log.h"
#include "truss/util/timer.h"
#include "truss/decompose/parallel_all_edge_cnc.h"
#include "truss/util/reordering/reorder_utils.h"
#include "truss/decompose/iter_helper.h"


#include <chrono>
static long long total_build_main_ms = 0;
static long long total_build_sub_ms  = 0;


// non-const K and L so an embedding app can set K/L dynamically for in-memory API runs.
extern int K, L; 
extern unsigned long long N;
std::vector<std::vector<int>> cliques_vec;

EBBkC_Graph_t::EBBkC_Graph_t() = default; //default syntax instructs the compiler to generate a simple, default constructor that performs no actions. All member variables are initialized via the build() method instead.

EBBkC_Graph_t::~EBBkC_Graph_t() {
    int i, k = K, node_size = v_size, link_size = e_size;

    if (is_sub) { //If the graph instance is a subproblem (is_sub is true), the sizes are adjusted based on the truss_num and a reduced k, as subproblems are smaller.
        k = K - 2;
        node_size = truss_num;
        link_size = truss_num * (truss_num - 1) / 2;
    }

    delete [] edges;

    // deallocate the 2D arrays T (common neighbors) and C (common edges for the next level)
    if (T) {
        for (i = 0; i < link_size; i++) delete [] T[i];
        delete [] T;
    }

    if (C) {
        for (i = 0; i < link_size; i++) delete [] C[i];
        delete [] C;
    }

    delete [] T_size;

    delete [] C_size;

    // sub_v and sub_e were used to hold the vertices and edges of recursive subproblems.
    if (sub_v) {
        for (i = 0; i <= k; i++) delete [] sub_v[i];
        delete [] sub_v;
    }

    if (sub_e) {
        for (i = 0; i <= k; i++) delete [] sub_e[i];
        delete [] sub_e;
    }

    delete [] sub_v_size;

    delete [] sub_e_size;

    delete [] lab;

    // DAG_deg and G_deg, which stored vertex degrees at different levels of recursion.    
    if (DAG_deg) {
        for (i = 0; i <= k; i++) delete [] DAG_deg[i];
        delete [] DAG_deg;

    }

    if (G_deg) {
        for (i = 0; i <= k; i++) delete [] G_deg[i];
        delete [] G_deg;
    }


    delete [] col; //col array, which stored the color of each vertex
    // DAG_adj and G_adj, which stored the adjacency lists for the colored DAG and the undirected graph, respectively
    if (DAG_adj) {
        for (i = 0; i < node_size; i++) delete [] DAG_adj[i];
        delete [] DAG_adj;
    }

    if (G_adj) {
        for (i = 0; i < node_size; i++) delete [] G_adj[i];
        delete [] G_adj;
    }
    // used 2D array, likely used as a temporary marker during coloring or other vertex processing.
    if (used) {
        for (i = 0; i <= k; i++) delete [] used[i];
        delete [] used;
    }


    delete [] v_lab;

    if (out_v_size) {
        for (i = 0; i <= k; i++) delete [] out_v_size[i];
        delete [] out_v_size;
    }

    if (out_e_size) {
        for (i = 0; i <= k; i++) delete [] out_e_size[i];
        delete [] out_e_size;
    }
    // (F, P, lack_size) used specifically for the early termination technique to handle plexes.
    delete [] F;

    delete [] P;

    delete [] lack_size;
    // lack 2D array, which stored the missing connections for vertices in a detected plex.
    if (lack) {
        for (i = 0; i < node_size; i++) delete [] lack[i];
        delete [] lack;
    }
    // helper arrays (lev, loc) used for the plex-listing algorithm.
    delete [] lev;

    delete [] loc;
}


//file writing function for dense-pce integration
void EBBkC_Graph_t::finalize_clique_output(const char* output_path) {
    FILE* f = fopen(output_path, "w");
    if (!f) {
        perror("Error opening clique output file");
        exit(1);
    }
    
    for (auto &cl : cliques_vec) {
        for (int vid : cl) {
            fprintf(f, "%d ", vid);
        }
        fprintf(f, "\n");
    }
    fclose(f);
    cliques_vec.clear();
}

void EBBkC_Graph_t::build(bool sub) { //allocates and initializes all the data structures that the destructor deallocates. It takes a boolean sub to indicate if it's building a subproblem graph
    auto build_start = std::chrono::high_resolution_clock::now();


    int i, k = K, node_size = v_size, link_size = e_size;

    is_sub = sub;

    if (sub) {
        k = K - 2;
        node_size = truss_num;
        link_size = (truss_num) * (truss_num - 1) / 2;
    }

    sub_v = new int* [k + 1];

    sub_e = new int* [k + 1];

    sub_e_size = new int [k + 1];

    sub_v_size = new int [k + 1];

    for (i = 0; i < k; i++) sub_v[i] = new int [truss_num + 1];
    sub_v[k] = new int [node_size]; // sub_v[d][…] holds the list of candidate vertices at recursion level d

    for (i = 0; i < k; i++) sub_e[i] = new int [truss_num * (truss_num - 1) / 2];
    sub_e[k] = new int [link_size]; // sub_e[d][…] holds the list of candidate edges at level d

    sub_v_size[k] = 0;
    // sub_v_size[d] and sub_e_size[d] track how many entries are currently in each list
    sub_e_size[k] = 0;

    //Allocates the lab (label) array and initializes all values to k, the current recursion depth.
    lab = new int [node_size]; // lab[u] stores the current “level” or depth at which vertex u was added; initializing all to k marks them as un–visited relative to the deepest recursion
    for (i = 0; i < node_size; i++) lab[i] = k;

    DAG_deg = new int* [k + 1]; // DAG_deg[d][u] = number of out‐neighbors of u in the directed acyclic graph (DAG) used by the color‐oriented branch & bound at level d
    for (i = 0; i <= k; i++) DAG_deg[i] = new int [node_size];

    G_deg = new int* [k + 1]; // G_deg[d][u] = number of neighbors of u in the undirected subgraph at level d
    for (i = 0; i <= k; i++) G_deg[i] = new int [node_size];

    col = new int [node_size]; // col[u] holds the greedy‐color assigned to vertex u

    DAG_adj = new int* [node_size]; // DAG_adj[u][0…DAG_deg[0][u]-1] lists all out‐neighbors of u in the oriented DAG
    for (i = 0; i < node_size; i++) DAG_adj[i] = new int [truss_num + 1];

    G_adj = new int* [node_size]; // G_adj[u][0…G_deg[0][u]-1] lists all neighbors of u in the undirected subgraph
    for (i = 0; i < node_size; i++) G_adj[i] = new int  [truss_num + 1];

    used = new bool* [k + 1]; // used[d][c] is a temporary flag array at level d to mark whether color c has been seen among processed neighbors—critical in the coloring heuristic and early‐termination checks
    for (i = 0; i <= k; i++) used[i] = new bool [node_size + 1]();

    v_lab = new int [node_size]; // v_lab[u]: per‐vertex label used in the plex‐listing routines (initialized to k
    for (i = 0; i < node_size; i++) v_lab[i] = k;

    out_v_size = new int* [k + 1];
    for (i = 0; i <= k; i++) out_v_size[i] = new int [link_size];
    // out_v_size[d][…], out_e_size[d][…]: track the sizes of the “frontier”—vertices and edges that can still extend a plex—at each level
    out_e_size = new int* [k + 1];
    for (i = 0; i <= k; i++) out_e_size[i] = new int [link_size];

    F = new int [truss_num + 1];
    // P and F partition vertices into “pending” vs “full‐degree” sets
    P = new int [truss_num + 1];

    lack_size = new int [node_size];

    lack = new int* [node_size]; // lack[u][…] lists which vertices each u is missing to complete its plex
    for (i = 0; i < node_size; i++) lack[i] = new int [L + 1];

    lev = new int [node_size]();
    // lev[u]/loc[u] track positions and temporary markers during the recursive plex enumeration
    loc = new int [node_size];


    auto build_end = std::chrono::high_resolution_clock::now();
    auto build_duration = std::chrono::duration_cast<std::chrono::milliseconds>(build_end - build_start).count();
    if (!sub) {
        total_build_main_ms += build_duration;
        printf("[DEBUG] build(main) took %ld ms\n", build_duration);
    } else {
        total_build_sub_ms += build_duration;
        printf("[DEBUG] build(sub) took %ld ms\n", build_duration);
    }
}

void EBBkC_Graph_t::truss_decompose(const char *dir) {
    graph_t g;

    //load the graph from file
    Graph G(dir); //load the graph from the file specified by dir
    g.adj = G.edge_dst; //adapts the loaded graph G into a local graph_t structure g
    g.num_edges = G.node_off;
    g.n = G.nodemax;
    g.m = G.edgemax;
    //calls ReorderWrapper to re-label vertex IDs, likely using a core-ordering heuristic to improve memory access patterns
    string reorder_method("core");

    vector <int32_t> new_vid_dict;
    vector <int32_t> old_vid_dict;
    ReorderWrapper(g, dir, reorder_method, new_vid_dict, old_vid_dict);

    //setting up data structures required by the truss decomposition library
    //edge list array
    Timer get_eid_timer;

    Edge *edgeIdToEdge = (Edge *) malloc((g.m / 2) * sizeof(Edge)); //mapping from edge IDs to edge objects (edgeIdToEdge)
    assert(edgeIdToEdge != nullptr);
    log_info("Malloc Time: %.9lf s", get_eid_timer.elapsed());
    get_eid_timer.reset();

    //Populate the edge list array
    getEidAndEdgeList(&g, edgeIdToEdge);
    log_info("Init Eid Time: %.9lf s", get_eid_timer.elapsed());
    get_eid_timer.reset();

    int *EdgeSupport = (int *) malloc(g.m / 2 * sizeof(int)); //an array to store the support of each edge (EdgeSupport)
    assert(EdgeSupport != nullptr);

    auto max_omp_threads = omp_get_max_threads();
    omp_set_num_threads(max_omp_threads);
    log_info("Max Threads: %d", max_omp_threads);
#pragma omp parallel for
    for (auto i = 0; i < max_omp_threads; i++) {
        auto avg = g.m / 2 / max_omp_threads;
        auto iter_beg = avg * i;
        auto iter_end = (i == max_omp_threads - 1) ? g.m / 2 : avg * (i + 1);
        memset(EdgeSupport + iter_beg, 0, (iter_end - iter_beg) * sizeof(int));
    }
    log_info("Init EdgeSupport Time: %.9lf s", get_eid_timer.elapsed());
    get_eid_timer.reset();

    Timer global_timer;
    // The key call is truss_num = PKT_intersection(...), which executes a parallel truss decomposition algorithm. This function computes the truss value and a rank for every edge in the graph. The maximum truss value found is stored in truss_num.
    truss_num = PKT_intersection(&g, EdgeSupport, edgeIdToEdge);
    // DEBUG: Print truss_num immediately after decomposition
    // printf("PKT_intersection returned truss_num = %d\n", truss_num); 
    // log_info("Truss decomposition done. truss_num = %d", truss_num);

#pragma omp single
    {
        int u, v, w;
        int *old2new = new int[N_NODES];
        for (int i = 0; i < N_NODES; i++) old2new[i] = -1;

        e_size = 0;
        edges = new Edge_t[g.m / 2];

        T = new int *[g.m / 2];
        T_size = new int[g.m / 2];

        int kept_edges = 0;
        int skipped_edges = 0;
        
        for (int i = 0; i < g.m / 2; i++) {
            if (g.edge_truss[i] < K) { // Changed from <= to <
                skipped_edges++;
                continue;
            }
            kept_edges++; //This is a crucial filter. If an edge's truss value (g.edge_truss[i]) is less than or equal to the target clique size K, it's impossible for that edge to be part of a K-clique. The loop skips such edges, pruning them from the graph structure that the main algorithm will use.

            Edge e = edgeIdToEdge[i];
            u = e.u;
            v = e.v;

            if (old2new[u] == -1) {
                old2new[u] = (int) new2old.size();
                new2old.push_back(u);
            }

            if (old2new[v] == -1) {
                old2new[v] = (int) new2old.size();
                new2old.push_back(v);
            }

            //For a valid edge, it's added to the EBBkC_Graph_t object's own edges list. Its new ID (e_size) is stored in a hash map for quick lookups, and its rank from the decomposition is stored.
            edges[e_size] = Edge_t(old2new[u], old2new[v], false);
            edge2id.insert(edges[e_size], e_size);
            rank.push_back(g.edge_rank[i]);

            //The code retrieves the set of common neighbors for the current edge, which the library has stored in g.v_set[i]. This set is copied into the T array (T[e_size]), which will be used during the branching phase.
            int sz = g.v_set[i].size();
            T_size[e_size] = sz;

            T[e_size] = new int[truss_num + 1];
            for (int j = 0; j < sz; j++) {
                int w = g.v_set[i][j];

                if (old2new[w] == -1) {
                    old2new[w] = (int) new2old.size();
                    new2old.push_back(w);
                }

                T[e_size][j] = old2new[w];
            }

            e_size++;
        }

        // DEBUG: Print edge filtering results
        //printf("Edge filtering: Total=%ld, Kept=%d, Skipped=%d\n", g.m/2, kept_edges, skipped_edges);
        //printf("Edge | Truss | Rank\n");
        // for (int i = 0; i < e_size; i++) {
        //     printf("%d-%d | %d | %d\n", 
        //         new2old[edges[i].s], 
        //         new2old[edges[i].t],
        //         g.edge_truss[i],  // Truss number
        //         rank[i]);         // Edge rank
        // }

        C = new int *[e_size];
        C_size = new int[e_size];

        //After processing all edges, this second loop populates the C array. For each edge i, it iterates through all pairs of its common neighbors (T[i][j], T[i][k]). It forms an edge e between this pair. If this edge e exists in the graph and has a higher rank than edge i, its ID (idx) is stored in C[i]. C[i] therefore contains all edges of the next recursive level that are "supported" by edge i.
        for (int i = 0; i < e_size; i++) {

            C_size[i] = 0;
            int sz = T_size[i];
            C[i] = new int[sz * (sz - 1) / 2];

            for (int j = 0; j < T_size[i]; j++) {
                for (int k = j + 1; k < T_size[i]; k++) {

                    Edge_t e = Edge_t(T[i][j], T[i][k], false);
                    int idx = edge2id.exist(e);

                    if (idx != -1 && rank[idx] > rank[i]) {
                        C[i][C_size[i]++] = idx;
                    }
                }
            }
        }
    };

    v_size = new2old.size();

    printf("|V| = %d, |E| = %d\n", v_size, e_size);
    printf("Truss number = %d\n", truss_num - 2);

    //Free memory
    free_graph(&g);
    free(edgeIdToEdge);
    free(EdgeSupport);

}

//The branch function implements the core edge-oriented branching strategy. It takes an edge e from the main, truss-ordered graph and generates a self-contained subproblem g for the recursive algorithm to solve. This also completes the hybrid ordering by applying a local color-based ordering.
void EBBkC_Graph_t::branch(int e, EBBkC_Graph_t* g) {
    // DEBUG: Print edge being processed
    // printf("Branching on edge %d: (%d->%d) [Main Graph IDs: %d,%d]\n",
    //        e, edges[e].s, edges[e].t, 
    //        new2old[edges[e].s], new2old[edges[e].t]);
    // printf("T_size[%d] = %d\n", e, T_size[e]);
    int c, i, j, k, p, e_, u, v, w, s, t, end, dist, l = K;
    int *old2new = new int[v_size]; //a temporary mapping (old2new) from the main graph's vertex IDs to the new, smaller IDs of the subproblem, v_size = total nodes in G, so large enough to map any vertex from the main graph.

    g->v_size = 0;
    g->e_size = 0;
    g->edges = new Edge_t[C_size[e] + 1];
    g->new2old = vector<int>(T_size[e]);

    // Initialize old2new with -1
    for (int i = 0; i < v_size; i++) {
        old2new[i] = -1;
    }
    
    for (i = 0; i < T_size[e]; i++) { // iterates through T[e], common neighbors of the branching edge e, contains original vertex IDs
        u = T[e][i]; // vertex ID from the original graph, "old" ID, e.g., 500
        Edge_t e1(edges[e].s, u, false);
        Edge_t e2(edges[e].t, u, false);
        int id1 = edge2id.exist(e1);
        int id2 = edge2id.exist(e2);

        if (id1 != -1 && id2 != -1 && rank[id1] > rank[e] && rank[id2] > rank[e]){
            old2new[u] = g->v_size; // Assign the next available local ID (e.g., 0, then 1, ...)
            g->new2old[g->v_size] = this->new2old[u]; // Stores the old ID 'u' at the new index
            g->v_size++;
        }
        
    }

    //This loop populates the edges of the new subgraph g. It uses the pre-computed C[e] array, which contains the IDs of all edges connecting the vertices in the subproblem. The vertex IDs are translated to the new local IDs using the old2new map.
    for (i = 0; i < C_size[e]; i++) { // iterates through C[e]
        e_ = C[e][i]; // edge ID from the original graph
        s = edges[e_].s; // source vertex ID from the original graph
        t = edges[e_].t; // target vertex ID from the original graph

        if (old2new[s] != -1 && old2new[t] != -1) {
            g->edges[g->e_size].s = old2new[s]; // uses the old2new map to convert them to local IDs and adds the new edge to g->edges
            g->edges[g->e_size++].t = old2new[t];
        }
        
    }

    // DEBUG: Print subgraph size
    // printf("Subgraph created: vertices=%d, edges=%d\n", g->v_size, g->e_size);
    // if (g->v_size > 0) {
    //     for(int i=0; i<g->v_size; i++){
    //         printf("Sample subgraph vertex: %d -> Original ID: %d\n", i, g->new2old[i]);
    //     }
    // }

    // u0, u1 are the endpoints of the e-th edge in the parent graph
    // g->branch_ends.resize(K+1);
    // g->branch_ends[K] = { 
    //     this->new2old[edges[e].s],  // Convert to original ID
    //     this->new2old[edges[e].t]   // Convert to original ID
    // };  // at the top level
    // at deeper levels, K-2, K-4, … you similarly record in branch_ends[l]


    // Initialize current clique with endpoints
    g->current_clique.clear();
    g->current_clique.push_back(this->new2old[edges[e].s]);
    g->current_clique.push_back(this->new2old[edges[e].t]);


    delete [] old2new; // Since it's only needed to build this one subproblem, it's deleted at the end of the branch function

    // if (l == 3 || l == 4) {
    //     g->sub_v_size[l - 2] = g->v_size;
    //     g->sub_e_size[l - 2] = g->e_size;
    //     // Initialize candidate vertex set
    //     for (int idx = 0; idx < g->v_size; idx++) {
    //         g->sub_v[l-2][idx] = idx;  // Store local vertex indices
    //     }
    //     if (l == 4) {
    //         //g->sub_e_size[l - 2] = g->e_size;
    //         for (int idx = 0; idx < g->e_size; idx++) {
    //             g->sub_e[l-2][idx] = idx;  // Store local edge indices
    //         }
    //     }
    //     return;
    // }

    // pruning check, If the newly created subgraph g is already too small to possibly contain a (K-2)-clique, the function returns immediately, and no further processing is done for this branch.
    if (g->v_size < l - 2 || g->e_size < (l - 2) * (l - 3) / 2) { // if the subgraph has fewer than l-2 vertices or fewer than (l-2)*(l-3)/2 edges, it's too small to contain a (K-2)-clique, so the function returns immediately
        g->sub_v_size[l - 2] = g->v_size;
        g->sub_e_size[l - 2] = g->e_size;
        return;
    }

    //This loop initializes the color and degree arrays for the subproblem. It sets all vertices' colors to 0 and their degrees to 0.
    for (j = 0; j < g->v_size; j++) {
        g->col[j] = 0;
        g->DAG_deg[l - 2][j] = 0;
        g->G_deg[l - 2][j] = 0;
    }
    // This loop populates the G_adj adjacency list for the subproblem. It iterates through all edges in g and adds them to the adjacency list, incrementing the degree of each vertex.
    for (j = 0; j < g->e_size; j++) {
        s = g->edges[j].s;
        t = g->edges[j].t;
        g->G_adj[s][g->G_deg[l - 2][s]++] = t;
        g->G_adj[t][g->G_deg[l - 2][t]++] = s;
    }

    //This section begins the local ordering for the subproblem. It first builds an undirected adjacency list (G_adj) for g. Then, it creates a list of vertices and sorts them based on their degree within the subgraph. This is the first step of the greedy coloring heuristic.
    auto *list = new KeyVal_t[truss_num + 1];
    for (j = 0; j < g->v_size; j++) {
        list[j].key = j;
        list[j].val = g->G_deg[l - 2][j];
    }
    sort(list, list + g->v_size);

    //This loop implements the greedy graph coloring algorithm. It iterates through the vertices u in their sorted order. For each vertex, it marks the colors of its already-colored neighbors as "used," finds the smallest color c that is not used, and assigns g->col[u] = c.
    for (j = 0; j < g->v_size; j++) {
        u = list[j].key;
        for (k = 0; k < g->G_deg[l - 2][u]; k++) {
            v = g->G_adj[u][k];
            g->used[l - 2][g->col[v]] = true;
        }
        for (c = 1; g->used[l - 2][c]; c++);
        g->col[u] = c;
        for (k = 0; k < g->G_deg[l - 2][u]; k++) {
            v = g->G_adj[u][k];
            g->used[l - 2][g->col[v]] = false;
        }
    }
    delete[] list;

    g->sub_v_size[l - 2] = 0;

    for (j = 0; j < g->v_size; j++) {
        g->sub_v[l - 2][g->sub_v_size[l - 2]++] = j;
    }

    g->sub_e_size[l - 2] = 0;
    //The final step is to build the Directed Acyclic Graph (DAG) that the recursive EBBkC_plus_plus function will use. It iterates through every edge in the subproblem and directs it from the vertex with the higher color to the vertex with the lower color. This orientation is stored in the DAG_adj adjacency list. This completes the preparation of the subproblem.
    for (j = 0; j < g->e_size; j++) {
        g->sub_e[l - 2][g->sub_e_size[l - 2]++] = j;
        s = g->edges[j].s;
        t = g->edges[j].t;
        g->edges[j].s = (g->col[s] > g->col[t]) ? s : t;
        g->edges[j].t = (g->col[s] > g->col[t]) ? t : s;
        s = g->edges[j].s;
        t = g->edges[j].t;

        g->DAG_adj[s][g->DAG_deg[l - 2][s]++] = t;
    }

    return;
}

void EBBkC_Graph_t::EBBkC_plus_plus(int l, unsigned long long *cliques) {
    // DEBUG: Print current recursion level
    //printf("Entering EBBkC_plus_plus: level=%d, |V|=%d, |E|=%d\n", l, sub_v_size[l], sub_e_size[l]);
    int c, i, j, k, p, e, e_, u, v, w, s, t, end, dist;

    if (sub_v_size[l] < l || sub_e_size[l] < l * (l - 1) / 2){
        //printf("Pruning at level %d: graph too small\n", l);
        return;
    }  // pruning check to see if the subproblem is large enough to contain an l-clique?



    // ===================== FIX 2: REORDER AND SIMPLIFY CHECKS =====================
    // Check for K=3 first. It's a problem-wide shortcut.
    if (K == 3) {
        for (int idx = 0; idx < sub_v_size[l]; idx++) { // l will be 1 here
            int w = sub_v[l][idx];
            int orig_w = this->new2old[w];
            // Add to current clique and store
            std::vector<int> full_clique = current_clique;
            full_clique.push_back(orig_w);
            if (clique_sink) clique_sink(full_clique); else cliques_vec.push_back(full_clique);
            (*cliques)++;
        }
        return;
    }

    // if (K == 4) {
    //     for (int idx = 0; idx < sub_e_size[l]; idx++) {
    //         e = sub_e[l][idx];
    //         int a = edges[e].s, b = edges[e].t;
    //         int orig_a = this->new2old[a];
    //         int orig_b = this->new2old[b];
    //         // Add edge vertices to current clique
    //         std::vector<int> full_clique = current_clique;
    //         full_clique.push_back(orig_a);
    //         full_clique.push_back(orig_b);
    //         cliques_vec.push_back(full_clique);
    //         (*cliques)++;
    //     }
    //     return;
    // }

    // Now, handle the general recursive base cases.
    if (l == 2) {
        for (i = 0; i < sub_v_size[l]; i++) {
            u = sub_v[l][i];
            for (j = 0; j < DAG_deg[l][u]; j++) {
                v = DAG_adj[u][j];
                // Add both vertices to current clique
                std::vector<int> full_clique = current_clique;
                full_clique.push_back(this->new2old[u]);
                full_clique.push_back(this->new2old[v]);
                if (clique_sink) clique_sink(full_clique); else cliques_vec.push_back(full_clique);
                (*cliques)++;
            }
        }
        return;
    }

    if (l == 3) { // base case, highly optimized routine to count all triangles in the current DAG subproblem 
        for (i = 0; i < sub_v_size[l]; i++) {
            u = sub_v[l][i];
            if (col[u] < l) continue;
            // mark neighbors of u
            for (j = 0; j < DAG_deg[l][u]; j++) lab[v = DAG_adj[u][j]] = l-1;
            for (j = 0; j < DAG_deg[l][u]; j++) {
                v = DAG_adj[u][j];
                if (col[v] < l-1) continue;
                for (k = 0; k < DAG_deg[l][v]; k++) {
                    w = DAG_adj[v][k];
                    if (lab[w] == l-1) {
                        // Add all three vertices
                        std::vector<int> full_clique = current_clique;
                        full_clique.push_back(this->new2old[u]);
                        full_clique.push_back(this->new2old[v]);
                        full_clique.push_back(this->new2old[w]);
                        if (clique_sink) clique_sink(full_clique); else cliques_vec.push_back(full_clique);
                        (*cliques)++;
                    }
                }
            }
            // reset marks
            for (j = 0; j < DAG_deg[l][u]; j++) lab[DAG_adj[u][j]] = l;
        }
        return;
    }


    // if (K == 3) { // if the target clique size is 3, the number of cliques is simply the number of vertices in the subproblem
    //     // DEBUG: Print found triangles
    //     //printf("K=3 base: Found %d triangles\n", sub_v_size[l]);
    //     // every vertex w in sub_v[l] completes a triangle {u0,u1,w}
    //     auto &ends = branch_ends[K];
    //     for (int idx = 0; idx < sub_v_size[l]; idx++) {
    //         int w = sub_v[l][idx];                       
    //         int orig_w = this->new2old[w];  // Already original ID
    //         //printf("Storing K3: (%d,%d,%d)\n", ends.first, ends.second, orig_w);
    //         cliques_vec.push_back({ends.first, ends.second, orig_w});
    //         (*cliques)++;
    //     }
    //     return;
    // }

    // if (K == 4) { // if the target clique size is 4, the number of cliques is simply the number of edges in the subproblem
    //     // DEBUG: Print found 4-cliques
    //     //printf("K=4 base: Found %d 4-cliques\n", sub_e_size[l]);
    //     // each edge e in sub_e[l] completes a 4-clique {u0,u1,u2,u3}
    //     auto &ends2 = branch_ends[K];        // contains {u0,u1}
    //     for (int idx = 0; idx < sub_e_size[l]; idx++) {
    //         e = sub_e[l][idx];
    //         int a = edges[e].s, b = edges[e].t;
    //         int orig_a = this->new2old[a];  // Already original
    //         int orig_b = this->new2old[b];  // Already original
    //         //printf("Storing K4: (%d,%d,%d,%d)\n", ends2.first, ends2.second, orig_a, orig_b);
    //         // record and count
    //         //cliques_vec.push_back({ ends2.first, ends2.second, a, b });
    //         cliques_vec.push_back({ends2.first, ends2.second, orig_a, orig_b});
    //         (*cliques)++;
    //     }
    //     return;
    // }

    // if (l == 2) {
    //     auto &top_ends = branch_ends[K];  // Get top-level endpoints
    //     for (i = 0; i < sub_v_size[l]; i++) {
    //         u = sub_v[l][i];
    //         for (j = 0; j < DAG_deg[l][u]; j++) {
    //             v = DAG_adj[u][j];
    //             // Include top-level endpoints + current edge
    //             cliques_vec.push_back({
    //                 top_ends.first,
    //                 top_ends.second,
    //                 this->new2old[u],
    //                 this->new2old[v]
    //             });
    //             (*cliques)++;
    //         }
    //     }
    //     return;
    // }

    // if (l == 3) { // base case, highly optimized routine to count all triangles in the current DAG subproblem 
    //     // fully enumerate each triangle in the DAG by storing {u,v,w}
    //     auto &top_ends = branch_ends[K];  // Get top-level endpoints
    //     for (i = 0; i < sub_v_size[l]; i++) {
    //         u = sub_v[l][i];
    //         if (col[u] < l) continue;
    //         // mark neighbors of u
    //         for (j = 0; j < DAG_deg[l][u]; j++) lab[v = DAG_adj[u][j]] = l-1;
    //         for (j = 0; j < DAG_deg[l][u]; j++) {
    //             v = DAG_adj[u][j];
    //             if (col[v] < l-1) continue;
    //             for (k = 0; k < DAG_deg[l][v]; k++) {
    //                 w = DAG_adj[v][k];
    //                 if (lab[w] == l-1) {
    //                     // Include top-level endpoints + triangle
    //                     cliques_vec.push_back({
    //                         top_ends.first,
    //                         top_ends.second,
    //                         this->new2old[u],
    //                         this->new2old[v],
    //                         this->new2old[w]
    //                     });
    //                     (*cliques)++;
    //                 }
    //             }
    //         }
    //         // reset marks
    //         for (j = 0; j < DAG_deg[l][u]; j++) lab[DAG_adj[u][j]] = l;
    //     }
    //     return;
    // }

    if (can_terminate(l, cliques)) { // check if the current subgraph is a "t-plex" (a very dense graph). If so, it uses a combinatorial method to count cliques, which is much faster than further recursion.
        //printf("Early termination at level %d\n", l);
        return;
    }

    for (i = 0; i < sub_v_size[l]; i++) { // main recursive loop, iterates through each vertex u in the subproblem.
        u = sub_v[l][i];

        if (col[u] < l) continue; // If a vertex's color is less than the size of the clique l we're looking for, it's impossible for it to be part of an l-clique with its lower-colored neighbors, so the entire branch is skipped.

        // Add vertex to current clique
        current_clique.push_back(this->new2old[u]);

        sub_v_size[l - 1] = 0; // resets the counter for the number of vertices in the subproblem at the next recursive level (l-1)
        dist = 0; //  initializes a counter named dist

        for (j = 0; j < DAG_deg[l][u]; j++) { //  iterates through all the neighbors of the pivot vertex u
            v = DAG_adj[u][j]; // retrieves the current neighbor
            lab[v] = l - 1;
            sub_v[l - 1][sub_v_size[l - 1]++] = v;
            DAG_deg[l - 1][v] = 0;
            G_deg[l - 1][v] = 0;

            if (!used[l][col[v]]) {
                used[l][col[v]] = true;
                dist++;
            }
        }

        if (dist >= l - 1) {

            sub_e_size[l - 1] = 0;
            for (j = 0; j < sub_v_size[l - 1]; j++) {
                v = sub_v[l - 1][j];

                end = DAG_deg[l][v];
                for (k = 0; k < end; k++) {
                    w = DAG_adj[v][k];
                    if (lab[w] == l - 1) {
                        DAG_deg[l - 1][v]++;
                        sub_e_size[l - 1]++;

                        // just for early-termination
                        G_deg[l - 1][v]++;
                        G_deg[l - 1][w]++;

                    } else {
                        DAG_adj[v][k--] = DAG_adj[v][--end];
                        DAG_adj[v][end] = w;
                    }
                }
            }

            EBBkC_plus_plus(l - 1, cliques);
        }

        for (j = 0; j < sub_v_size[l - 1]; j++) {
            v = sub_v[l - 1][j];
            lab[v] = l;
            used[l][col[v]] = false;
        }
        // Remove vertex after recursion
        current_clique.pop_back();
    }
    // DEBUG: After recursion
    //printf("Exiting EBBkC_plus_plus at level %d\n", l);
}

void EBBkC_Graph_t::EBBkC_Comb_list(int *list, int  list_size, int  start, int  picked, int  k, unsigned long long *cliques) {
    // buffer for the current partial combination (local IDs)
    static std::vector<int> cur;
    if (picked == 0) 
        cur.clear();

    // BASE CASE: we’ve picked k vertices → record & count
    // if (picked == k) {
    //     std::vector<int> orig;
    //     orig.reserve(k + 2);  // Reserve space for clique + endpoints
        
    //     // Add top-level endpoints first
    //     auto &top_ends = this->branch_ends[K];
    //     orig.push_back(top_ends.first);
    //     orig.push_back(top_ends.second);
        
    //     // Add vertices from combinatorial list
    //     for (int t = 0; t < k; t++) {
    //         int main_idx = this->new2old[cur[t]];
    //         orig.push_back(main_idx);
    //     }
        
    //     cliques_vec.push_back(orig);
    //     (*cliques)++;
    //     return;
    // }

    if (picked == k) {
        std::vector<int> full_clique = current_clique;  // Start with current clique
        for (int t = 0; t < k; t++) {
            full_clique.push_back(this->new2old[cur[t]]);
        }
        if (clique_sink) clique_sink(full_clique); else cliques_vec.push_back(full_clique);
        (*cliques)++;
        return;
    }

    // RECURSION: try each remaining candidate
    for (int i = start; i < list_size; i++) {
        cur.push_back(list[i]); 
        EBBkC_Comb_list(list, list_size, i + 1, picked + 1, k, cliques);
        cur.pop_back();
    }
}

void EBBkC_Graph_t::list_in_plex(int start, int p, int q, unsigned long long *cliques) {
    if (F_size < q) return;

    if (p == 0) {
        // if (q > F_size - q)
        //     EBBkC_Comb_list(F, F_size, 0, 0, F_size - q, cliques);
        // else
        //     EBBkC_Comb_list(F, F_size, 0, 0, q, cliques);
        EBBkC_Comb_list(F, F_size, 0, 0, q, cliques);
        return;
    }

    int i, j, u, v, vis = 0;

    for (i = start; i < P_size && P_act >= p; i++) {
        u = P[i];

        if (lev[u]) continue;

        for (j = 0; j < lack_size[u]; j++) {
            v = lack[u][j];
            if (loc[v] >= i && lev[v] == 0) {
                lev[v] = p;
                P_act--;
            }
        }

        // FIX: Add the selected vertex 'u' to the current clique before recursing.
        // We use this->new2old[u] to get its original, global ID.
        current_clique.push_back(this->new2old[u]);

        list_in_plex(i + 1, p - 1, q, cliques);

        // FIX: Backtrack by removing the vertex after the recursive call returns.
        current_clique.pop_back();

        for (j = 0; j < lack_size[u]; j++) {
            v = lack[u][j];
            if (loc[v] >= i && lev[v] == p) {
                lev[v] = 0;
                P_act++;
            }
        }

        P_act--;
        vis++;
    }
    P_act += vis;
}

bool EBBkC_Graph_t::can_terminate(int l, unsigned long long *cliques) {
    // Only allow early termination at the top level (l == K-2)
    if (l != K - 2) {
        return false;
    }
    int i, j, k, u, v, end, p_;

    if (sub_e_size[l] < sub_v_size[l] * (sub_v_size[l] - L) / 2) return false;

    if (sub_e_size[l] == sub_v_size[l] * (sub_v_size[l] - 1) / 2) {
        // FIX: Always choose l vertices (needed for clique formation)
        EBBkC_Comb_list(sub_v[l], sub_v_size[l], 0, 0, l, cliques);
        return true;
    }

    if (L == 1) return false;

    for (i = 0; i < sub_v_size[l]; i++) {
        u = sub_v[l][i];

        if (sub_v_size[l] - G_deg[l][u] > L) {
            return false;
        }
    }

    F_size = 0;
    P_size = 0;

    for (i = 0; i < sub_v_size[l]; i++) {
        u = sub_v[l][i];

        if (G_deg[l][u] == sub_v_size[l] - 1) {
            loc[u] = -1;
            F[F_size++] = u;
            continue;
        }

        loc[u] = P_size;
        P[P_size++] = u;
    }

    int* e = new int [P_size * P_size]();

    for (i = 0; i < P_size; i++) {
        u = P[i];
        lack_size[u] = 0;

        end = DAG_deg[l][u];
        for (j = 0; j < end; j++) {
            v = DAG_adj[u][j];

            if (loc[v] != -1) {
                e[loc[u] * P_size + loc[v]] = 1;
                e[loc[v] * P_size + loc[u]] = 1;
            }
        }
    }

    for (i = 0; i < P_size * P_size; i++) {
        if (!e[i]) {
            j = i / P_size, k = i % P_size;
            u = P[j], v = P[k];
            lack[u][lack_size[u]++] = v;
        }
    }

    delete [] e;

    for(i = 0; i <= l; i++) {
        P_act = P_size;
        list_in_plex(0, i, l - i, cliques);
    }

    return true;
}


double EBBkC_t::list_k_clique(const char *file_name) {
    double runtime;
    struct rusage start, end;
    EBBkC_Graph_t G, g; // G is the original graph, g is the subproblem graph

    printf("Reading edges from %s ...\n", file_name);
    GetCurTime(&start);
    G.truss_decompose(file_name); // truss decomposition

    printf("Building necessary data structure ...\n");
    G.build(false); // build the data structure for the original graph, false argument indicates this is not a subproblem

    printf("Iterate over all cliques\n");

    g.truss_num = G.truss_num; // copy the maximum truss number from the original graph to the subproblem graph
    g.build(true); // build the data structure for the subproblem graph, true argument signifies it's a subproblem
    for (int i = 0; i < G.e_size; i++) { // iterate over all edges in the original graph
        G.branch(i, &g); // takes the common neighbors of edge i (T[i]) and the edges connecting them (C[i]) to build a new, small graph g. It also performs a local, color-based ordering on g.
        g.EBBkC_plus_plus(K - 2, &N); // asked to find (K-2)-cliques in the subproblem graph g. When combined with the original edge i (which has 2 vertices), any (K-2)-clique found in g forms a K-clique in the original graph. The count of found cliques is added to the global counter N.
    }

    GetCurTime(&end);
    runtime = GetTime(&start, &end);

    printf("Number of %u-cliques: %llu\n", K, N);

    printf("[DEBUG] Total build(main) time: %lld ms\n", total_build_main_ms);
    printf("[DEBUG] Total build(sub) time: %lld ms\n", total_build_sub_ms);
    printf("[DEBUG] Total build time: %lld ms\n", total_build_main_ms + total_build_sub_ms);
    printf("[DEBUG] Average build(sub) time: %.2f ms\n", G.e_size ? (double)total_build_sub_ms / G.e_size : 0.0);


    // Modified output path to include K in filename
    std::string input_path(file_name);
    if (!input_path.empty() && input_path.back() == '/') {
        input_path.pop_back();
    }
    std::string output_path = input_path + "/cliques_K" + std::to_string(K);
    G.finalize_clique_output(output_path.c_str());
    
    // for (auto &cl : cliques_vec) {
    //     printf("clique:");
    //     for (int vid : cl) printf(" %d", vid);
    //     printf("\n");
    // }
    
    // clear for next run
    cliques_vec.clear();

    return runtime;
}

double EBBkC_t::list_k_clique_mem(const char* dir, int k, int l, std::vector<std::vector<int>>&  out_cliques) {
    // Configure globals for this run
    K = k;
    L = l;
    N = 0ULL;
    cliques_vec.clear();

    EBBkC_Graph_t G, g;
    double t0 = omp_get_wtime();
    G.truss_decompose(dir);
    G.build(false);
    g.truss_num = G.truss_num;
    g.build(true);
    // Use default buffered mode: collect results
    g.clique_sink = nullptr;
    for (int i = 0; i < G.e_size; i++) {
        G.branch(i, &g);
        g.EBBkC_plus_plus(K - 2, &N);
    }
    double ms = (double)(omp_get_wtime() - t0) * 1e3;

    printf("Number of %u-cliques: %llu\n", K, N);
    printf("[DEBUG] Total build(main) time: %lld ms\n", total_build_main_ms);
    printf("[DEBUG] Total build(sub) time: %lld ms\n", total_build_sub_ms);
    printf("[DEBUG] Total build time: %lld ms\n", total_build_main_ms + total_build_sub_ms);
    printf("[DEBUG] Average build(sub) time: %.2f ms\n", G.e_size ? (double)total_build_sub_ms / G.e_size : 0.0);

    out_cliques.swap(cliques_vec);
    cliques_vec.clear();
    return ms;
}

double EBBkC_t::list_k_clique_mem_stream(const char* dir, int k, int l, const std::function<void(const std::vector<int>&)>& sink) {
    // Configure globals for this run
    K = k;
    L = l;
    N = 0ULL;
    cliques_vec.clear();

    EBBkC_Graph_t G, g;
    double t0 = omp_get_wtime();
    G.truss_decompose(dir);
    G.build(false);
    g.truss_num = G.truss_num;
    g.build(true);
    // Streaming mode: deliver each clique immediately via sink
    g.clique_sink = sink;
    for (int i = 0; i < G.e_size; i++) {
        G.branch(i, &g);
        g.EBBkC_plus_plus(K - 2, &N);
    }
    double ms = (double)(omp_get_wtime() - t0) * 1e3;

    // No accumulation; ensure buffer is empty
    cliques_vec.clear();
    return ms;
}

double EBBkC_t::list_k_clique_mem_stream_from_csr(
    uint32_t n, uint32_t m,
    const uint32_t* node_off,
    const int* edge_dst,
    int k, int l,
    const std::function<void(const std::vector<int>&)>& sink) {
    K = k; L = l; N = 0ULL;
    cliques_vec.clear();

    // Build graph_t from borrowed CSR
    graph_t gmain{};
    gmain.n = n;
    gmain.m = m;
    gmain.adj = const_cast<int*>(edge_dst);
    gmain.num_edges = const_cast<eid_t*>(reinterpret_cast<const eid_t*>(node_off));
    gmain.borrowed = true;

    // Construct wrapper Graph-like structure minimal for truss code
    EBBkC_Graph_t G, g;
    double t0 = omp_get_wtime();

    // Reuse existing build path by simulating post-decompose state:
    // We still need PKT_intersection; it expects a graph_t with fields set like in truss_decompose.
    graph_t gg{};
    gg.adj = gmain.adj;
    gg.num_edges = gmain.num_edges;
    gg.n = n;
    gg.m = m;

    Edge *edgeIdToEdge = (Edge *) malloc((gg.m / 2) * sizeof(Edge));
    assert(edgeIdToEdge != nullptr);
    getEidAndEdgeList(&gg, edgeIdToEdge);

    int *EdgeSupport = (int *) malloc(gg.m / 2 * sizeof(int));
    assert(EdgeSupport != nullptr);
    int max_threads = omp_get_max_threads();
#pragma omp parallel for
    for (int i = 0; i < max_threads; i++) {
        auto avg = gg.m / 2 / max_threads;
        auto iter_beg = avg * i;
        auto iter_end = (i == max_threads - 1) ? gg.m / 2 : avg * (i + 1);
        memset(EdgeSupport + iter_beg, 0, (iter_end - iter_beg) * sizeof(int));
    }

    // Run decomposition
    int truss_num = PKT_intersection(&gg, EdgeSupport, edgeIdToEdge);

    // Populate EBBkC_Graph_t structures similarly to truss_decompose()
    {
        int u, v;
        int *old2new = new int[n];
        for (uint32_t i = 0; i < n; i++) old2new[i] = -1;
        G.e_size = 0;
        G.edges = new Edge_t[gg.m / 2];
        G.T = new int *[gg.m / 2];
        G.T_size = new int[gg.m / 2];
        G.rank.clear();

        for (uint32_t i = 0; i < gg.m / 2; i++) {
            // Filter by K like in truss_decompose path
            if (gg.edge_truss && gg.edge_truss[i] < K) continue;
            Edge e = edgeIdToEdge[i];
            u = e.u; v = e.v;
            if (old2new[u] == -1) { old2new[u] = (int) G.new2old.size(); G.new2old.push_back(u); }
            if (old2new[v] == -1) { old2new[v] = (int) G.new2old.size(); G.new2old.push_back(v); }
            G.edges[G.e_size] = Edge_t(old2new[u], old2new[v], false);
            G.edge2id.insert(G.edges[G.e_size], G.e_size);
            G.rank.push_back(gg.edge_rank ? gg.edge_rank[i] : i);
            int sz = gg.v_set[i].size();
            G.T_size[G.e_size] = sz;
            G.T[G.e_size] = new int[truss_num + 1];
            for (int j = 0; j < sz; j++) {
                int w = gg.v_set[i][j];
                if (old2new[w] == -1) { old2new[w] = (int) G.new2old.size(); G.new2old.push_back(w); }
                G.T[G.e_size][j] = old2new[w];
            }
            G.e_size++;
        }

        G.v_size = (int)G.new2old.size();
        G.truss_num = truss_num;

        // Build C and C_size like in truss_decompose()
        G.C = new int *[G.e_size];
        G.C_size = new int[G.e_size];
        for (int i = 0; i < G.e_size; i++) {
            G.C_size[i] = 0;
            int sz = G.T_size[i];
            G.C[i] = new int[sz * (sz - 1) / 2];
            for (int j = 0; j < G.T_size[i]; j++) {
                for (int k = j + 1; k < G.T_size[i]; k++) {
                    Edge_t e = Edge_t(G.T[i][j], G.T[i][k], false);
                    int idx = G.edge2id.exist(e);
                    if (idx != -1 && G.rank[idx] > G.rank[i]) {
                        G.C[i][G.C_size[i]++] = idx;
                    }
                }
            }
        }

        // Build g (subproblem container)
        g.truss_num = G.truss_num;
        g.build(true);

        // Stream results by setting sink on g
        g.clique_sink = sink;
        for (int i = 0; i < G.e_size; i++) {
            G.branch(i, &g);
            g.EBBkC_plus_plus(K - 2, &N);
        }

        delete [] old2new;
    }

    // Cleanup
    free(edgeIdToEdge);
    free(EdgeSupport);
    if (gg.edge_truss) free(gg.edge_truss);
    if (gg.edge_rank) free(gg.edge_rank);
    double ms = (double)(omp_get_wtime() - t0) * 1e3;
    cliques_vec.clear();
    return ms;
}

double EBBkC_t::list_k_clique_parallel(const char *file_name) {
    double runtime = 0, runtime1 = 0;
    int n_edges, i;
    // struct rusage start, end;
    double start, end;
    EBBkC_Graph_t G, g;

    printf("Reading edges from %s ...\n", file_name);
    start = omp_get_wtime();
    G.truss_decompose(file_name);

    printf("Building necessary data structure ...\n");
    G.build(false);

    printf("Iterate over all cliques\n");

#pragma omp parallel private(g, start, end, n_edges, i) reduction(+:N)
    {
        n_edges = 0;
        g.truss_num = G.truss_num;
        g.build(true);
#pragma omp for schedule(dynamic, 1) nowait
        for (i = 0; i < G.e_size; i++) {
            G.branch(i, &g);
            g.EBBkC_plus_plus(K - 2, &N);
            n_edges++;
        }
        printf("Thread: %d, handled %d edges\n", omp_get_thread_num(), n_edges);
    }

    return (double) (omp_get_wtime() - start) * 1e3;
}