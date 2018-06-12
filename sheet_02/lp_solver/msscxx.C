/**
    This file is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    It is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this file.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <iostream>
#include <stdio.h>
#include <vector>
#include <math.h>
#include <string.h>
#include <iomanip>
#include <cassert>

#ifdef __cplusplus
extern "C" {
#include "lp.h"
}
#else

#include "lp.h"

#endif

using Edge = std::vector<int>;

/** Read a graph instance from a file.
    f       : filename containing the problem graph
    pncount : number of vertices in the problem graph (will be filled)
              the vertices are numbered by 0,...,pncount-1 (different from the DIMACS-format!)
    pecount : number of edges in the problem graph (will be filled)
              the edges are numbered by 0,...,pecount-1
    pelist  : edge-list, consisting of an array of length 2*pecount containing
              edges as consecutive node-pairs, 
              i.e. edge i is connecting pelist[2*i] and pelist[2*i + 1]
    pnweights: weights of the vertices, consisting of an array of length pncount.
 */
int read_dimacs(char *f, int *pncount, int *pecount, int **pelist, int **pnweights);

/** Build an lp for the instance (ncount, ecount, elist, nweights). */
int build_lp(LINGOlp **lp, int ncount, int ecount, int *elist, int *nweights);

/* Converts the int array 'pelist' of the edges returned by read_dimacs into
 * std::vector<std::vector<int>>
 * a vector of tupels of int
 * */
std::vector <Edge> convert_elist_to_vector(const int elist[], const int ecount);

/* Converts the int array 'pnweights' of the vertex positions returned by read_dimacs into
 * std::vector<int>
 * a vector of ints
 * */
std::vector<int> convert_to_vector(const int nweights[], const int ncount);

std::vector <std::vector<size_t>> neighbors(const size_t number_of_vertices,
                                            const std::vector <Edge> &edge_vector);

double average_position_of_neighbours(const std::vector<double> &vertex_vector,
                                      const size_t vertex_index,
                                      const std::vector <std::vector<size_t>> &neighbors_vector);

std::vector<double> iterate_average_position_of_neighbours(const std::vector<int> &vertex_vector,
                                                           const std::vector<bool> &fixed_vertices_bool,
                                                           const std::vector <Edge> &edge_vector);

/*Overloading the '<<' operator to print vectors*/
template<typename T>
std::ostream &operator<<(std::ostream &out, const std::vector <T> &v) {
    out << "[";
    size_t last = v.size() - 1;
    for (size_t i = 0; i < v.size(); ++i) {
        out << std::setw(3) << v[i];
        if (i != last)
            out << ", ";
    }
    out << "]";
    return out;
}

template<typename T>
T linear_length(const std::vector <T> &vertex_vector, const std::vector <Edge> &edge_vector) {
    T sum = 0;
    for (auto edge : edge_vector) {
        sum += std::abs(vertex_vector[edge[0]] - vertex_vector[edge[1]]);
    }
    return sum;
}

template<typename T>
T quadratic_length(const std::vector <T> &vertex_vector, const std::vector <Edge> &edge_vector) {
    T sum = 0;
    for (auto edge : edge_vector) {
        sum += std::pow(std::abs(vertex_vector[edge[0]] - vertex_vector[edge[1]]), 2);
    }
    return sum;
}

int read_dimacs(char *f, int *pncount, int *pecount, int **pelist,
                int **pnweights) {
    int rval = 0;
    int haveprob = 0;
    int icount = 0;

    int ncount = 0;
    int ecount = 0;
    int *elist = (int *) NULL;
    int *nweights = (int *) NULL;
    int nnweights = 0;

    int i, end0, end1;

    int n;
    char buf[256];
    char *p;

    FILE *in = (FILE *) NULL;

    *pncount = *pecount = 0;
    if (*pelist) free(*pelist);
    if (*pnweights) free(*pnweights);

    in = fopen(f, "r");
    if (!in) {
        fprintf(stderr, "Unable to open %s for input\n", f);
        rval = 1;
        goto CLEANUP;
    }

    while (fgets(buf, 254, in) != (char *) NULL) {
        p = buf;
        if (p[0] == 'p') {
            const char *delim = " \t\n";
            char *data = (char *) NULL;
            strtok(p, delim); /* get 'p' */

            if (haveprob) {
                fprintf(stderr, "ERROR in Dimacs file -- two p-lines.\n");
                rval = 1;
                goto CLEANUP;
            }
            haveprob = 1;
            data = strtok(NULL, delim); /* get type */
            if (strcmp(data, "edge") && strcmp(data, "edges") &&
                strcmp(data, "col") && strcmp(data, "graph")) {
                fprintf(stderr, "ERROR in Dimacs file -- not an edge file\n");
                rval = 1;
                goto CLEANUP;
            }
            data = strtok(NULL, delim);
            sscanf(data, "%d", &ncount);
            data = strtok(NULL, delim);
            sscanf(data, "%d", &ecount);

            printf("Number of Nodes: %d\n", ncount);
            printf("Number of Edges: %d\n", ecount);

            elist = (int *) malloc(2 * ecount * sizeof(int));
            LINGOcheck_NULL (elist, "out of memory for elist");
            nweights = (int *) malloc(ncount * sizeof(int));
            LINGOcheck_NULL (nweights, "out of memory for nweights");
            for (i = 0; i < ncount; i++) nweights[i] = 1;
        } else if (p[0] == 'e') {
            if (!haveprob) {
                fprintf(stderr, "ERROR in Dimacs file -- no problem defined \n");
                rval = 1;
                goto CLEANUP;
            }
            if (icount >= ecount) {
                fprintf(stderr, "ERROR in Dimacs file -- to many edges\n");
                rval = 1;
                goto CLEANUP;
            }
            p++;
            sscanf(p, "%d %d", &end0, &end1);
            elist[2 * icount] = end0 - 1;    /* Number nodes from 0, not 1 */
            elist[2 * icount + 1] = end1 - 1;
            icount++;
        } else if (p[0] == 'n') {
            int weight;
            if (!haveprob) {
                fprintf(stderr, "ERROR in Dimacs file -- n before p\n");
                rval = 1;
                goto CLEANUP;
            }
            p++;
            sscanf(p, "%d %d", &n, &weight);
            nweights[n - 1] = weight;
            nnweights++;
        }
    }

    /* Some dimacs col-instances are buggy => reduce # edges to icount*/
    *pncount = ncount;
    *pecount = icount;
    *pelist = elist;
    if (pnweights) {
        *pnweights = nweights;
    } else {
        if (nweights) free(nweights);
    }
    CLEANUP:
    if (rval) {
        if (elist) free(elist);
        if (nweights) free(nweights);
    }

    if (in) fclose(in);

    return rval;
}


int build_lp(LINGOlp **lp,
             int ncount, int ecount, int *elist, int *nweights) {
    int rval = 0;
    int i;

    rval = LINGOlp_init(lp, "MSSlp");
    LINGOcheck_rval (rval, "LINGOlp_init failed");

    rval = LINGOlp_objective_sense(*lp, LINGOlp_MAX);
    LINGOcheck_rval (rval, "LINGOlp_objective_sense");

    for (i = 0; i < ncount; i++) {
        double w = (double) nweights[i];
        /* As no rows were defined yet, nzcount, rind, and rval are 0.*/
        rval = LINGOlp_addcol(*lp, 0, (int *) NULL, (double *) NULL,
                              w, 0.0, 1.0, LINGOlp_CONTINUOUS, NULL);
        LINGOcheck_rval (rval, "LINGOlp_addcol failed");
    }

    for (i = 0; i < ecount; i++) {
        int v = elist[2 * i];
        int w = elist[2 * i + 1];
        int count = 2;
        int inodes[2];
        double coef[2] = {1.0, 1.0};
        double rhs = 1.0;
        inodes[0] = v;
        inodes[1] = w;

        rval = LINGOlp_addrow(*lp, count, inodes, coef, LINGOlp_LESS_EQUAL,
                              rhs, NULL);
        if (rval) LINGOlp_printerrorcode(rval);
        LINGOcheck_rval (rval, "LINGOlp_addrow failed");

        /** To see how you can incrementally add rows and optimize
        comment out the following block.
        */
        /* { */
        /*    double obj; */
        /*    rval = LINGOlp_optimize (*lp); */
        /*    LINGOcheck_rval (rval, "LINGOlp_optimize failed"); */

        /*    rval = LINGOlp_objval ( *lp, &obj); */
        /*    LINGOcheck_rval (rval, "LINGOlp_objval failed"); */
        /*    printf("Objective after adding %i constraints/rows: %f.\n",i+1,obj); */
        /* } */

    }
    CLEANUP:
    return rval;
}

std::vector <Edge> convert_elist_to_vector(const int elist[], const int ecount) {
    std::vector <Edge> edge_vector;
    for (int i = 0; i < 2 * ecount; i += 2) {
        Edge new_edge;
        new_edge.push_back(elist[i]);
        new_edge.push_back(elist[i + 1]);
        edge_vector.push_back(new_edge);
    }
    return edge_vector;
}

std::vector<int> convert_to_vector(const int nweights[], const int ncount) {
    std::vector<int> vertex_vector;
    for (int i = 0; i < ncount; i++) {
        vertex_vector.push_back(nweights[i]);
    }
    return vertex_vector;
}

// Returns a vector of length (# of vertices)
// where each entry is the adjacency list of the corresponding vertex
std::vector <std::vector<size_t>> neighbors(const size_t number_of_vertices,
                                            const std::vector <Edge> &edge_vector) {
    std::vector <std::vector<size_t>> neighbors_vector(number_of_vertices);
    for (auto edge : edge_vector) {
        neighbors_vector[edge[0]].push_back(edge[1]);
        neighbors_vector[edge[1]].push_back(edge[0]);
    }
    return neighbors_vector;
}

/*  This function assumes that each vertex has at least one neighbor
 *
 * */
double average_position_of_neighbours(const std::vector<double> &vertex_vector,
                                      const size_t vertex_index,
                                      const std::vector <std::vector<size_t>> &neighbors_vector) {
    double sum = 0;
    assert(neighbors_vector[vertex_index].size() >= 1);
    assert(vertex_index < vertex_vector.size());
    for (auto index : neighbors_vector[vertex_index]) {
        sum += vertex_vector[index];
    }
    return sum / (double) neighbors_vector[vertex_index].size();
}

std::vector<double> iterate_average_position_of_neighbours(const std::vector<int> &vertex_vector,
                                                           const std::vector<bool> &fixed_vertices_bool,
                                                           const std::vector <Edge> &edge_vector) {
    double maximum_circuit_movement = 0;
    // Can choose any initial position according to the exercise
    // Just copy the values for the fixed vector

    std::vector<double> new_vertex_positions;
    for (size_t i = 0; i < vertex_vector.size(); i++) {
        if (fixed_vertices_bool[i] == true) {
            new_vertex_positions.push_back((double) vertex_vector[i]);
        } else {
            // Set the initial position of all unfixed vertices to 1
            new_vertex_positions.push_back((double) 1.0);
        }
    }

    std::vector <std::vector<size_t>> neighbors_vector = neighbors(vertex_vector.size(), edge_vector);
//    std::cout << "neighbors_vector:" << std::endl;
//    std::cout << neighbors_vector << std::endl;

    std::vector <size_t> indices_that_are_updated;
    for (size_t i = 0; i < vertex_vector.size(); i++) {
        if (fixed_vertices_bool[i] == false) {
            indices_that_are_updated.push_back(i);
        }
    }
// std::cout << indices_that_are_updated << std::endl;

    size_t iteration = 0;
    do {
        maximum_circuit_movement = 0;
        for (auto index : indices_that_are_updated) {
            double new_position = average_position_of_neighbours(new_vertex_positions, index, neighbors_vector);
            if (std::abs(new_position - new_vertex_positions[index]) > maximum_circuit_movement) {
                maximum_circuit_movement = std::abs(new_position - new_vertex_positions[index]);
            }
            new_vertex_positions[index] = new_position;
        }
//        std::cout << "Iteration "
//                  << iteration
//                  << ", Quadratic length: "
//                  << quadratic_length(new_vertex_positions, edge_vector)
//                  << std::endl;
//        std::cout << "new_vertex_positions:" << std::endl;
//        std::cout << new_vertex_positions << std::endl;
        iteration++;
    } while (maximum_circuit_movement >= 0.1);

    return new_vertex_positions;
}


int main(int argc, char **argv) {
    // Initializations
    char *dimacs_fname = (char *) NULL;
    int rval = 0;
    int ncount, ecount;
    int *elist = (int *) NULL;
    int *nweights = (int *) NULL;
    LINGOlp *lp = (LINGOlp *) NULL;
    double *x = (double *) NULL;

    std::vector <Edge> edge_vector;
    std::vector<int> vertex_vector;
    // This vector indicates if the specified vertex at the position has to remain
    // fixed (i.e. it is preplaced)
    std::vector<bool> fixed_vertices_bool;

    std::vector<double> result_a;

    if (argc < 2) {
        printf("Usage mss <filename>\n");
        rval = 1;
        goto CLEANUP;
    }
    dimacs_fname = argv[1];

    // Read DIMACS file
    rval = read_dimacs(dimacs_fname, &ncount, &ecount, &elist, &nweights);
    LINGOcheck_rval (rval, "read_dimacs failed");

    edge_vector = convert_elist_to_vector(elist, ecount);
    vertex_vector = convert_to_vector(nweights, ncount);
    for (size_t i = 0; i < vertex_vector.size(); i++) {
        if (vertex_vector[i] == -1) {
            fixed_vertices_bool.push_back(false);
        } else {
            fixed_vertices_bool.push_back(true);
        }
    }

    std::cout << "edge_vector:" << std::endl;
    std::cout << edge_vector << std::endl;
    std::cout << "vertex_vector:" << std::endl;
    std::cout << vertex_vector << std::endl;
    std::cout << "fixed_vertices_bool:" << std::endl;
    std::cout << fixed_vertices_bool << std::endl;

    result_a = iterate_average_position_of_neighbours(vertex_vector, fixed_vertices_bool, edge_vector);
    std::cout << "-----(a)-----" << std::endl;
    std::cout << "result_a:" << std::endl;
    std::cout << result_a << std::endl;
    std::cout << "linear_length_a: " << linear_length(result_a, edge_vector) << std::endl;
    std::cout << "quadratic_length_a: " << quadratic_length(result_a, edge_vector) << std::endl;
    std::cout << "Positions g of the circuits C:" << std::endl;
    for(size_t i = 0; i < vertex_vector.size(); i++){
        if(fixed_vertices_bool[i] == false){
            std::cout << i << " " << result_a[i] << std::endl;
        }
    }

//  rval = build_lp (&lp,  ncount,  ecount, elist, nweights);
//  LINGOcheck_rval (rval, "build_lp failed");
//
//  rval = LINGOlp_optimize (lp);
//  LINGOcheck_rval (rval, "LINGOlp_optimize failed");
//
//  /** Allocate an array for storing the primal solution.*/
//  x = (double*) malloc(ncount * sizeof(double));
//  LINGOcheck_NULL(x,"Failed to allocate result vector x.");
//
//  /** Retrieve the primal solution.*/
//  rval =  LINGOlp_x (lp,x);
//  LINGOcheck_rval (rval, "LINGOlp_x failed");
//
//  printf ("Printing solution:\n");
//  for (int i = 0; i < ncount; ++i) {
//    printf ("node %d val %f.\n", i, x[i]);
//  }

    CLEANUP:
    if (elist) free(elist);
    if (nweights) free(nweights);
    if (x) free(x);
    LINGOlp_free(&lp);

    return rval;
}
