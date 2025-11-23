#include <iostream>
#include <cstdlib> //standard library functions
#include <cstdio> //C-style I/O
#include <cstring>
#include <ctime>
#include <set>
#include <algorithm>
#include <string>
#include <omp.h> //OpenMP, used to enable parallel programming
#include <cassert> //C-style assertions, used for debugging
#include "set_operation.h"
#include "def.h"
#include "edge_oriented.h"


using namespace std;
int K, L; // K=size of cliques, L=plex parameter for the early termination technique
unsigned long long N = 0; // global 64-bit unsigned integer to count total number of cliques found

int main(int argc, char** argv) {
    double runtime;
    string act(argv[1]); //retrieves the first command-line argument (e.g., "p", "e", or "ep")

    if (act == "p") {                           // pre-process
        printf("Pre-process %s\n", argv[2]); // argv[2] is the input filename
        string src_filename(argv[2]);
//extract the file extension (suffix) from the source filename and check if it is either .edges or .mtx. If not, the program exits.
        string suffix = src_filename.substr(src_filename.find_last_of('.'));
        if (suffix != ".edges" && suffix != ".mtx") exit(0);
// construct a new filename with a .clean extension.
        string prefix = src_filename.substr(0, src_filename.find_last_of('.'));
        string clean_filename = prefix + ".clean";
        clean_edges(src_filename.c_str(), clean_filename.c_str()); // removes duplicate edges and self-loops
//        string index_filename = prefix + ".index";
//        runtime = EBBkC_t::truss_order(clean_filename.c_str(), index_filename.c_str());
        printf("Pre-processed in %.2lf ms\n\n", runtime);
    }

    else if (act == "e") { //"e" for sequential execution EBBkC+ET = EBBkC algorithm with Early Termination (ET)
        string src_filename(argv[2]);
//        string suffix = src_filename.substr(src_filename.find_last_of('.'));
//        if (suffix != ".index") exit(0);

        K = atoi(argv[3]); // atoi = string to an integer
        L = atoi(argv[4]);

        runtime = EBBkC_t::list_k_clique(argv[2]); // main call to execute the sequential k-clique listing algorithm
        //printf("Number of %u-cliques: %llu\n", K, N);
        printf("EBBkC+ET (t = %d) runtime %.2lf ms\n\n", L, runtime);
    }

    else if (act == "ep") {                      // EBBkC+ET (parallel)
        string src_filename(argv[2]);
//        string suffix = src_filename.substr(src_filename.find_last_of('.'));
//        if (suffix != ".index") exit(0);

        K = atoi(argv[3]); 
        L = atoi(argv[4]);

        omp_set_num_threads(atoi(argv[5])); //sets the number of threads for the subsequent parallel region

        runtime = EBBkC_t::list_k_clique_parallel(argv[2]);
        printf("Number of %u-cliques: %llu\n", K, N);
        printf("EBBkC+ET (t = %d) runtime %.2lf ms\n\n", L, runtime);
    }

    else {
        printf("Wrong usage.\n");
        exit(0);
    }

    return 0;
}
