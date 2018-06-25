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
#include <algorithm>

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
int build_lp(LINGOlp **lp,
             const std::vector<int> &vertex_vector,
             const std::vector<bool> &fixed_vertices_bool,
             const std::vector <Edge> &edge_vector);

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

std::vector<double> positions_with_lp(const std::vector<int> &vertex_vector,
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
    // assert(neighbors_vector[vertex_index].size() >= 1);
    if (neighbors_vector[vertex_index].size() == 0) {
        std::cout << "Vertex does not have neighbors" << std::endl;
        return vertex_vector[vertex_index];
    }

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

int build_lp(LINGOlp **lp,
             const std::vector<int> &vertex_vector,
             const std::vector<bool> &fixed_vertices_bool,
             const std::vector <Edge> &edge_vector) {
    int rval = 0;
    size_t i;

    rval = LINGOlp_init(lp, "MSSlp");
    LINGOcheck_rval (rval, "LINGOlp_init failed");

    rval = LINGOlp_objective_sense(*lp, LINGOlp_MIN);
    LINGOcheck_rval (rval, "LINGOlp_objective_sense");

    /** Add a column/variable to the LP.
    p       : The LP
    nzcount : number of non-zero row-entries of the new column.
    rind    : array of length nzcount, specifying the (non-zero) row/constraint
              indices of the new column.
    rvals   : array of length nzcount, specifying the coefficients of the
              rows/contraints specified in rind.
    obj     : coefficient in the objective (c-vector) of the new variable.
    lb      : lower bound for the new variable.
    ub      : upper bound for the new variable.
    vartype : Variable type: one of LINGOlp_CONTINUOUS, LINGOlp_INTEGER and LINGOlp_BINARY.
              The latter two are not valid in this context.
    name    : name of the constraint, which can be set to (char*) NULL.

    int LINGOlp_addcol (LINGOlp *p, int nzcount, int *rind, double *rvals,
                       double obj, double lb, double ub, char vartype, char *name);
    */

    // Add a variable g_v for each vertex v
    // The variables in C (vertices that are to be placed)
    // are constrained to be in the range {1, ..., k=vertex_vector.size()}
    // The preplaced vertices are constrained to their fixed position
    // The cost/coefficient in the c-vector is zero for these variables
    for (i = 0; i < vertex_vector.size(); i++) {
        /* As no rows were defined yet, nzcount, rind, and rval are 0.*/
        if (fixed_vertices_bool[i] == false) {
            // Free vertex
            rval = LINGOlp_addcol(*lp, 0, (int *) NULL, (double *) NULL,
                                  0, 1.0, (double) vertex_vector.size(), LINGOlp_CONTINUOUS, NULL);
            LINGOcheck_rval (rval, "LINGOlp_addcol failed");
        } else {
            // Fixed vertex
            rval = LINGOlp_addcol(*lp, 0, (int *) NULL, (double *) NULL,
                                  0, (double) vertex_vector[i], (double) vertex_vector[i], LINGOlp_CONTINUOUS, NULL);
            LINGOcheck_rval (rval, "LINGOlp_addcol failed");
        }
    }
    // Add a variable delta_e for each edge
    // Later we will add the constraints for e=(v,w):
    // g_v - g_w <= delta_e
    // g_w - g_v <= delta_e
    // This ensures that delta_e is set to be greater or equal to
    // the absolute difference between start and end position
    //
    // The cost/coefficient in the c-vector is set to 1
    for (auto edge : edge_vector) {
        rval = LINGOlp_addcol(*lp, 0, (int *) NULL, (double *) NULL,
                              1, 0.0, (double) vertex_vector.size(), LINGOlp_CONTINUOUS, NULL);
        LINGOcheck_rval (rval, "LINGOlp_addcol failed");
    }

    /** Add a row/constraint to the LP (sparsely). Parameters:
    p       : The LP
    nzcount : number of (non-zero) column-entries in the new row.
    cind    : array of length nzcount, specifying the (non-zero) column/variable
              indices of the new row.
    cval    : array of length nzcount, specifying the coefficients of the
              columns/variables specified in cind.
    sense   : type of inequalitym either LINGOlp_EQUAL, LINGOlp_LESS_EQUAL, or
              LINGOlp_GREATER_EQUAL
    rhs     : the right hand side (b-component) of the constraint.
    name    : name of the constraint, which can be set to (char*) NULL.

    int LINGOlp_addrow (LINGOlp *p, int nzcount, int *cind, double *cval,
                        char sense, double rhs, char *name);
    */

    for (i = 0; i < edge_vector.size(); i++) {
        Edge edge = edge_vector[i];
        int nzcount = 3;
        double cval[3] = {1.0, -1.0, -1.0};
        double rhs = 0.0;
        int cind[3];
        cind[0] = edge[0];
        cind[1] = edge[1];
        cind[2] = vertex_vector.size() + i;

        // g_v - g_w - delta_e <= 0
        rval = LINGOlp_addrow(*lp, nzcount, cind, cval, LINGOlp_LESS_EQUAL,
                              rhs, NULL);
        if (rval) LINGOlp_printerrorcode(rval);
        LINGOcheck_rval (rval, "LINGOlp_addrow failed");

        cval[0] = -1.0;
        cval[1] = 1.0;
        // - g_v + g_w - delta_e <= 0
        rval = LINGOlp_addrow(*lp, nzcount, cind, cval, LINGOlp_LESS_EQUAL,
                              rhs, NULL);
        if (rval) LINGOlp_printerrorcode(rval);
        LINGOcheck_rval (rval, "LINGOlp_addrow failed");
    }

    CLEANUP:
    return rval;
}

std::vector<double> positions_with_lp(const std::vector<int> &vertex_vector,
                                      const std::vector<bool> &fixed_vertices_bool,
                                      const std::vector <Edge> &edge_vector) {
    int rval = 0;
    LINGOlp *lp = (LINGOlp *) NULL;
    double *x = (double *) NULL;
    std::vector<double> result_b;

    rval = build_lp(&lp, vertex_vector, fixed_vertices_bool, edge_vector);
    LINGOcheck_rval (rval, "build_lp failed");

    rval = LINGOlp_optimize(lp);
    LINGOcheck_rval (rval, "LINGOlp_optimize failed");

    /** Allocate an array for storing the primal solution.*/
    x = (double *) malloc((vertex_vector.size() + edge_vector.size()) * sizeof(double));
    LINGOcheck_NULL(x, "Failed to allocate result vector x.");

    /** Retrieve the primal solution.*/
    rval = LINGOlp_x(lp, x);
    LINGOcheck_rval (rval, "LINGOlp_x failed");
    // Convert the node coordinates to vector
    result_b = std::vector<double>(x, x + vertex_vector.size());

    CLEANUP:
    if (x) free(x);
    LINGOlp_free(&lp);
    return result_b;
}


bool are_all_vertices_fixed(const std::vector<bool> &fixed_vertices_bool) {
    for (auto value : fixed_vertices_bool) {
        if (value == false) {
            return false;
        }
    }
    return true;
}

template<typename T>
T median(std::vector <T> v) {
    size_t n = v.size() / 2;
    std::nth_element(v.begin(), v.begin() + n, v.end());
    return v[n];
}

std::vector<double> unassigned_positions(const std::vector<double> &vertex_vector,
                                         const std::vector<bool> &fixed_vertices_bool) {
    std::vector<double> unassigned_positions;
    for (size_t i = 0; i < vertex_vector.size(); i++) {
        if (fixed_vertices_bool[i] == false) {
            unassigned_positions.push_back(vertex_vector[i]);
        }
    }
    return unassigned_positions;
}

std::vector<int> available_positions(const std::vector<int> &vertex_vector,
                                     const std::vector<bool> &fixed_vertices_bool,
                                     const int left_border,
                                     __attribute__((unused)) const int right_border) {
    std::vector<int> available_positions;

    // Start with the entire range
    for (size_t i = 0; i < vertex_vector.size(); i++) {
        available_positions.push_back(left_border + i);
    }


    for (size_t i = 0; i < vertex_vector.size(); i++) {
        if (vertex_vector[i] == -1) continue;
        if (fixed_vertices_bool[i] == true) {
            auto it = std::find(available_positions.begin(), available_positions.end(), vertex_vector[i]);
            if (it != available_positions.end())
                available_positions.erase(it);
        }
    }

    return available_positions;
}

std::vector<int> exercise_c_with_a(const std::vector<int> &vertex_vector,
                                   const std::vector<bool> &fixed_vertices_bool,
                                   const std::vector <Edge> &edge_vector,
                                   const int left_border, const int right_border) {
    if (are_all_vertices_fixed(fixed_vertices_bool) == true) {
        //std::cout << "Returning because fixed_vertices_bool: " << fixed_vertices_bool << std::endl;
        return vertex_vector;
    } else if (vertex_vector.size() == 0) {
        return vertex_vector;
    }

    // Obtain some g which is not necessarily injective and integral
//    std::cout << "****** Calling iterate_average_position_of_neighbours with" << std::endl;
//    std::cout << "vertex_vector:" << vertex_vector << std::endl;
//    std::cout << "fixed_vertices_bool:" << fixed_vertices_bool << std::endl;
//    std::cout << "edge_vector:" << edge_vector << std::endl;
    std::vector<double> noninjective_fractional_solution = iterate_average_position_of_neighbours(vertex_vector,
                                                                                                  fixed_vertices_bool,
                                                                                                  edge_vector);

//    std::cout << "noninjective_fractional_solution: " << noninjective_fractional_solution << std::endl;
    std::vector<double> unassigned_pos = unassigned_positions(noninjective_fractional_solution, fixed_vertices_bool);
//    std::cout << "unassigned_pos: " << unassigned_pos << std::endl;

    double median_unassigned_value = median(unassigned_pos);

    size_t index_of_median_unassigned_value;
    // Find index of median unfixed vertex
    for (size_t i = 0; i < noninjective_fractional_solution.size(); i++) {
        if (fixed_vertices_bool[i] == true) continue;
        if (noninjective_fractional_solution[i] == median_unassigned_value) {
            index_of_median_unassigned_value = i;
            break;
        }
    }
//    std::cout << "median_unassigned_value: " << median_unassigned_value << std::endl;
//    std::cout << "index_of_median_unassigned_value: " << index_of_median_unassigned_value << std::endl;

    std::vector<int> available_pos = available_positions(vertex_vector, fixed_vertices_bool, left_border, right_border);
//    std::cout << "available_positions: " << available_pos << std::endl;

    int median_available_position = median(available_pos);
//    std::cout << "median_available_position: " << median_available_position << std::endl;

    // Pick median unassigned circuit and assign it to the median available position in f
    std::vector<int> vertex_vector_new(vertex_vector);
    std::vector<bool> fixed_vertices_bool_new(fixed_vertices_bool);
    vertex_vector_new[index_of_median_unassigned_value] = median_available_position;
    fixed_vertices_bool_new[index_of_median_unassigned_value] = true;

//    std::cout << "vertex_vector_new: " << vertex_vector_new << std::endl;


    // Make and solve subproblems

    // LEFT

    size_t left_interval_size = index_of_median_unassigned_value + 1;
    size_t right_interval_size = vertex_vector_new.size() - index_of_median_unassigned_value;

    std::vector<int> vertex_vector_left(left_interval_size, -1);
    std::vector<bool> fixed_vertices_bool_left(left_interval_size, false);
    std::vector <Edge> edge_vector_left;

    std::vector<bool> already_assigned_to_a_side(noninjective_fractional_solution.size(), false);

    // Fixed vertices
    size_t open_spots_on_left_side = left_interval_size;
    for(size_t i = 0; i < vertex_vector_left.size(); i++){
        if(fixed_vertices_bool_new[i] == true){
            vertex_vector_left[i] = vertex_vector_new[i];
            fixed_vertices_bool_left[i] = true;
            already_assigned_to_a_side[i] = true;
            open_spots_on_left_side--;
        }
    }

//    std::cout << "~~~~~~~~~~~~~~~~~~open_spots_on_left_side: " << open_spots_on_left_side << std::endl;
//    std::cout << "~~~~~~~~~~~~~~~~~~vertex_vector_left: " << vertex_vector_left << std::endl;
//    std::cout << "~~~~~~~~~~~~~~~~~~fixed_vertices_bool_left: " << fixed_vertices_bool_left << std::endl;

    // Fill remaining spots with unfixed vertices
    size_t count = 0;
    for (size_t i = 0; i < vertex_vector_left.size(); i++) {
        if (count > open_spots_on_left_side) break;
        if (vertex_vector_left[i] == -1){
            for(size_t j = 0; j < vertex_vector_new.size(); j++){
                if(already_assigned_to_a_side[j] == false){
                    already_assigned_to_a_side[j] = true;
                    vertex_vector_left[i] = noninjective_fractional_solution[j];
                    count++;
                    break;
                }
            }
        }
    }

//    std::cout << "~~~~~~~~~~~~~~~~~~open_spots_available_on_left_side: " << open_spots_on_left_side << std::endl;
//    std::cout << "~~vertex_vector_left: " << vertex_vector_left << std::endl;
//    std::cout << "~~fixed_vertices_bool_left: " << fixed_vertices_bool_left << std::endl;


    // edges
    for (auto edge : edge_vector) {
        if (edge[0] < (int) vertex_vector_left.size() && edge[1] < (int) vertex_vector_left.size()) {
            edge_vector_left.push_back(edge);
            continue;
        }
        if (edge[0] < (int) vertex_vector_left.size() - 1 && edge[1] >= (int) vertex_vector_left.size()) {
            edge[1] = vertex_vector_left.size() - 1;
            edge_vector_left.push_back(edge);
        }
    }

//    std::cout << "edge_vector_left: " << edge_vector_left << std::endl;


    // RIGHT

    std::vector<int> vertex_vector_right(right_interval_size, -1);
    std::vector<bool> fixed_vertices_bool_right(right_interval_size, false);
    std::vector <Edge> edge_vector_right;



    size_t open_spots_on_right_side = right_interval_size;
    for (size_t i = 0; i < vertex_vector_right.size(); i++) {
        if(fixed_vertices_bool_new[i + left_interval_size - 1] == true){
            vertex_vector_right[i] = vertex_vector_new[i + left_interval_size - 1];
            fixed_vertices_bool_right[i] = true;
            already_assigned_to_a_side[i + left_interval_size - 1] = true;
            open_spots_on_right_side--;
        }
    }

    for(size_t i = 0; i < vertex_vector_right.size(); i++){
        if(vertex_vector_right[i] == -1){
            for(size_t j = 0; j < vertex_vector_new.size(); j++){
                if(already_assigned_to_a_side[j] == false){
                    already_assigned_to_a_side[j] = true;
                    vertex_vector_right[i] = noninjective_fractional_solution[j];
                    break;
                }
            }
        }
    }

//    std::cout << "~~vertex_vector_right: " << vertex_vector_right << std::endl;
//    std::cout << "~~fixed_vertices_bool_right: " << fixed_vertices_bool_right << std::endl;

    for (auto edge : edge_vector) {
        if (edge[0] >= (int) vertex_vector_left.size() && edge[1] >= (int) vertex_vector_left.size()) {
            edge[0] = edge[0] - vertex_vector_left.size() + 1;
            edge[1] = edge[1] - vertex_vector_left.size() + 1;
            edge_vector_right.push_back(edge);
            continue;
        }
        if (edge[0] <= (int) vertex_vector_left.size() && edge[1] >= (int) vertex_vector_left.size()) {
            edge[0] = 0;
            edge[1] = edge[1] - vertex_vector_left.size() + 1;
            edge_vector_right.push_back(edge);
        }
    }

//    std::cout << "edge_vector_right: " << edge_vector_right << std::endl;


    std::vector<int> solution_left = exercise_c_with_a(vertex_vector_left, fixed_vertices_bool_left, edge_vector_left,
                                                       left_border, left_border + left_interval_size);
    std::vector<int> solution_right = exercise_c_with_a(vertex_vector_right, fixed_vertices_bool_right,
                                                        edge_vector_right, index_of_median_unassigned_value,
                                                        index_of_median_unassigned_value + right_interval_size);

//    std::cout << "solution_left: " << solution_left << std::endl;
//    std::cout << "solution_right: " << solution_right << std::endl;

    std::vector<int> solution(vertex_vector.size());
//    std::cout << "solution_initialized: " << solution << std::endl;

    for (size_t i = 0; i < solution_left.size(); i++) {
        solution[i] = solution_left[i];
    }
    // Skip first element of the right solution (this is the median which was already included in left solution)
    for (size_t i = 1; i < solution_right.size(); i++) {
        solution[solution_left.size() + i - 1] = solution_right[i];
    }

//    std::cout << "solution: " << solution << std::endl;

    return solution;
}

std::vector<int> exercise_c_with_b(const std::vector<int> &vertex_vector,
                                   const std::vector<bool> &fixed_vertices_bool,
                                   const std::vector <Edge> &edge_vector,
                                   const int left_border, const int right_border) {
    if (are_all_vertices_fixed(fixed_vertices_bool) == true) {
        //std::cout << "Returning because fixed_vertices_bool: " << fixed_vertices_bool << std::endl;
        return vertex_vector;
    } else if (vertex_vector.size() == 0) {
        return vertex_vector;
    }

    // Obtain some g which is not necessarily injective and integral
//    std::cout << "****** Calling iterate_average_position_of_neighbours with" << std::endl;
//    std::cout << "vertex_vector:" << vertex_vector << std::endl;
//    std::cout << "fixed_vertices_bool:" << fixed_vertices_bool << std::endl;
//    std::cout << "edge_vector:" << edge_vector << std::endl;
    std::vector<double> noninjective_fractional_solution = positions_with_lp(vertex_vector,
                                                                                                  fixed_vertices_bool,
                                                                                                  edge_vector);

//    std::cout << "noninjective_fractional_solution: " << noninjective_fractional_solution << std::endl;
    std::vector<double> unassigned_pos = unassigned_positions(noninjective_fractional_solution, fixed_vertices_bool);
//    std::cout << "unassigned_pos: " << unassigned_pos << std::endl;

    double median_unassigned_value = median(unassigned_pos);

    size_t index_of_median_unassigned_value;
    // Find index of median unfixed vertex
    for (size_t i = 0; i < noninjective_fractional_solution.size(); i++) {
        if (fixed_vertices_bool[i] == true) continue;
        if (noninjective_fractional_solution[i] == median_unassigned_value) {
            index_of_median_unassigned_value = i;
            break;
        }
    }
//    std::cout << "median_unassigned_value: " << median_unassigned_value << std::endl;
//    std::cout << "index_of_median_unassigned_value: " << index_of_median_unassigned_value << std::endl;

    std::vector<int> available_pos = available_positions(vertex_vector, fixed_vertices_bool, left_border, right_border);
//    std::cout << "available_positions: " << available_pos << std::endl;

    int median_available_position = median(available_pos);
//    std::cout << "median_available_position: " << median_available_position << std::endl;

    // Pick median unassigned circuit and assign it to the median available position in f
    std::vector<int> vertex_vector_new(vertex_vector);
    std::vector<bool> fixed_vertices_bool_new(fixed_vertices_bool);
    vertex_vector_new[index_of_median_unassigned_value] = median_available_position;
    fixed_vertices_bool_new[index_of_median_unassigned_value] = true;

//    std::cout << "vertex_vector_new: " << vertex_vector_new << std::endl;


    // Make and solve subproblems

    // LEFT

    size_t left_interval_size = index_of_median_unassigned_value + 1;
    size_t right_interval_size = vertex_vector_new.size() - index_of_median_unassigned_value;

    std::vector<int> vertex_vector_left(left_interval_size, -1);
    std::vector<bool> fixed_vertices_bool_left(left_interval_size, false);
    std::vector <Edge> edge_vector_left;

    std::vector<bool> already_assigned_to_a_side(noninjective_fractional_solution.size(), false);

    // Fixed vertices
    size_t open_spots_on_left_side = left_interval_size;
    for(size_t i = 0; i < vertex_vector_left.size(); i++){
        if(fixed_vertices_bool_new[i] == true){
            vertex_vector_left[i] = vertex_vector_new[i];
            fixed_vertices_bool_left[i] = true;
            already_assigned_to_a_side[i] = true;
            open_spots_on_left_side--;
        }
    }

//    std::cout << "~~~~~~~~~~~~~~~~~~open_spots_on_left_side: " << open_spots_on_left_side << std::endl;
//    std::cout << "~~~~~~~~~~~~~~~~~~vertex_vector_left: " << vertex_vector_left << std::endl;
//    std::cout << "~~~~~~~~~~~~~~~~~~fixed_vertices_bool_left: " << fixed_vertices_bool_left << std::endl;

    // Fill remaining spots with unfixed vertices
    size_t count = 0;
    for (size_t i = 0; i < vertex_vector_left.size(); i++) {
        if (count > open_spots_on_left_side) break;
        if (vertex_vector_left[i] == -1){
            for(size_t j = 0; j < vertex_vector_new.size(); j++){
                if(already_assigned_to_a_side[j] == false){
                    already_assigned_to_a_side[j] = true;
                    vertex_vector_left[i] = noninjective_fractional_solution[j];
                    count++;
                    break;
                }
            }
        }
    }

//    std::cout << "~~~~~~~~~~~~~~~~~~open_spots_available_on_left_side: " << open_spots_on_left_side << std::endl;
//    std::cout << "~~vertex_vector_left: " << vertex_vector_left << std::endl;
//    std::cout << "~~fixed_vertices_bool_left: " << fixed_vertices_bool_left << std::endl;


    // edges
    for (auto edge : edge_vector) {
        if (edge[0] < (int) vertex_vector_left.size() && edge[1] < (int) vertex_vector_left.size()) {
            edge_vector_left.push_back(edge);
            continue;
        }
        if (edge[0] < (int) vertex_vector_left.size() - 1 && edge[1] >= (int) vertex_vector_left.size()) {
            edge[1] = vertex_vector_left.size() - 1;
            edge_vector_left.push_back(edge);
        }
    }

//    std::cout << "edge_vector_left: " << edge_vector_left << std::endl;


    // RIGHT

    std::vector<int> vertex_vector_right(right_interval_size, -1);
    std::vector<bool> fixed_vertices_bool_right(right_interval_size, false);
    std::vector <Edge> edge_vector_right;



    size_t open_spots_on_right_side = right_interval_size;
    for (size_t i = 0; i < vertex_vector_right.size(); i++) {
        if(fixed_vertices_bool_new[i + left_interval_size - 1] == true){
            vertex_vector_right[i] = vertex_vector_new[i + left_interval_size - 1];
            fixed_vertices_bool_right[i] = true;
            already_assigned_to_a_side[i + left_interval_size - 1] = true;
            open_spots_on_right_side--;
        }
    }

    for(size_t i = 0; i < vertex_vector_right.size(); i++){
        if(vertex_vector_right[i] == -1){
            for(size_t j = 0; j < vertex_vector_new.size(); j++){
                if(already_assigned_to_a_side[j] == false){
                    already_assigned_to_a_side[j] = true;
                    vertex_vector_right[i] = noninjective_fractional_solution[j];
                    break;
                }
            }
        }
    }

//    std::cout << "~~vertex_vector_right: " << vertex_vector_right << std::endl;
//    std::cout << "~~fixed_vertices_bool_right: " << fixed_vertices_bool_right << std::endl;

    for (auto edge : edge_vector) {
        if (edge[0] >= (int) vertex_vector_left.size() && edge[1] >= (int) vertex_vector_left.size()) {
            edge[0] = edge[0] - vertex_vector_left.size() + 1;
            edge[1] = edge[1] - vertex_vector_left.size() + 1;
            edge_vector_right.push_back(edge);
            continue;
        }
        if (edge[0] <= (int) vertex_vector_left.size() && edge[1] >= (int) vertex_vector_left.size()) {
            edge[0] = 0;
            edge[1] = edge[1] - vertex_vector_left.size() + 1;
            edge_vector_right.push_back(edge);
        }
    }

//    std::cout << "edge_vector_right: " << edge_vector_right << std::endl;


    std::vector<int> solution_left = exercise_c_with_a(vertex_vector_left, fixed_vertices_bool_left, edge_vector_left,
                                                       left_border, left_border + left_interval_size);
    std::vector<int> solution_right = exercise_c_with_a(vertex_vector_right, fixed_vertices_bool_right,
                                                        edge_vector_right, index_of_median_unassigned_value,
                                                        index_of_median_unassigned_value + right_interval_size);

//    std::cout << "solution_left: " << solution_left << std::endl;
//    std::cout << "solution_right: " << solution_right << std::endl;

    std::vector<int> solution(vertex_vector.size());
//    std::cout << "solution_initialized: " << solution << std::endl;

    for (size_t i = 0; i < solution_left.size(); i++) {
        solution[i] = solution_left[i];
    }
    // Skip first element of the right solution (this is the median which was already included in left solution)
    for (size_t i = 1; i < solution_right.size(); i++) {
        solution[solution_left.size() + i - 1] = solution_right[i];
    }

//    std::cout << "solution: " << solution << std::endl;

    return solution;
}

int main(int argc, char **argv) {
    // Initializations
    char *dimacs_fname = (char *) NULL;
    int rval = 0;
    int ncount, ecount;
    int *elist = (int *) NULL;
    int *nweights = (int *) NULL;

    std::vector <Edge> edge_vector;
    std::vector<int> vertex_vector;

    // This vector indicates if the specified vertex at the position has to remain
    // fixed (i.e. it is preplaced)
    std::vector<bool> fixed_vertices_bool;

    std::vector<double> result_a;
    std::vector<double> result_b;
    std::vector<int> result_c_a;
    std::vector<int> result_c_b;

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

    std::cout << std::endl << "BEGIN -----(a)-----" << std::endl;

    result_a = iterate_average_position_of_neighbours(vertex_vector, fixed_vertices_bool, edge_vector);
    std::cout << "result_a:" << std::endl;
    std::cout << result_a << std::endl;
    std::cout << "linear_length_a: " << linear_length(result_a, edge_vector) << std::endl;
    std::cout << "quadratic_length_a: " << quadratic_length(result_a, edge_vector) << std::endl;
    std::cout << "Positions g of the circuits C for (a):" << std::endl;
    for (size_t i = 0; i < vertex_vector.size(); i++) {
        if (fixed_vertices_bool[i] == false) {
            std::cout << i << " " << result_a[i] << std::endl;
        }
    }

    std::cout << "END -----(a)-----" << std::endl << std::endl;

    //-------------------------------

    std::cout << std::endl << "BEGIN -----(b)-----" << std::endl;

    result_b = positions_with_lp(vertex_vector, fixed_vertices_bool, edge_vector);
    std::cout << "result_b:" << std::endl;
    std::cout << result_b << std::endl;
    std::cout << "linear_length_b: " << linear_length(result_b, edge_vector) << std::endl;
    std::cout << "quadratic_length_b: " << quadratic_length(result_b, edge_vector) << std::endl;
    std::cout << "Positions g of the circuits C for (b):" << std::endl;
    for (size_t i = 0; i < vertex_vector.size(); i++) {
        if (fixed_vertices_bool[i] == false) {
            std::cout << i << " " << result_b[i] << std::endl;
        }
    }

    std::cout << "END -----(b)-----" << std::endl << std::endl;

    //------------------------------

    // Method using (a)

    std::cout << std::endl << "BEGIN -----(c_a)-----" << std::endl;

    result_c_a = exercise_c_with_a(vertex_vector, fixed_vertices_bool, edge_vector, 1, vertex_vector.size());
    std::cout << "result_c_a:" << std::endl;
    std::cout << result_c_a << std::endl;
    std::cout << "linear_length_c_a: " << linear_length(result_c_a, edge_vector) << std::endl;
    std::cout << "quadratic_length_c_a: " << quadratic_length(result_c_a, edge_vector) << std::endl;
    std::cout << "Positions g of the circuits C for (c_a):" << std::endl;
    for (size_t i = 0; i < vertex_vector.size(); i++) {
            std::cout << i << " " << result_c_a[i] << std::endl;
    }

    std::cout << "Sorry, injectivity does not work yet, one would have to pass a list of all used positions or something similar" << std::endl;

    std::cout << "END -----(c_a)-----" << std::endl << std::endl;


    // Method using (b)

    std::cout << std::endl << "BEGIN -----(c_b)-----" << std::endl;

    result_c_b = exercise_c_with_b(vertex_vector, fixed_vertices_bool, edge_vector, 1, vertex_vector.size());
    std::cout << "result_c_b:" << std::endl;
    std::cout << result_c_b << std::endl;
    std::cout << "linear_length_c_b: " << linear_length(result_c_b, edge_vector) << std::endl;
    std::cout << "quadratic_length_c_b: " << quadratic_length(result_c_b, edge_vector) << std::endl;
    std::cout << "Positions g of the circuits C for (c_b):" << std::endl;
    for (size_t i = 0; i < vertex_vector.size(); i++) {
            std::cout << i << " " << result_c_b[i] << std::endl;
    }

    std::cout << "Sorry, injectivity does not work yet, one would have to pass a list of all used positions or something similar" << std::endl;

    std::cout << "END -----(c_b)-----" << std::endl << std::endl;


    CLEANUP:
    if (elist) free(elist);
    if (nweights) free(nweights);


    return rval;
}
