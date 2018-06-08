//
// Created by ben on 22.05.18.
//

#include <iostream>
#include <algorithm>
#include <sstream>
#include <vector>
#include <iomanip>
#include <cassert>
#include <string>
#include <limits.h>

// --------------------------------------------------
// Helper functions

/*Overloading the '<<' operator to print vectors*/
template<typename T>
std::ostream &operator<<(std::ostream &out, const std::vector<T> &v) {
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

// Extract the list of x-coordinates from a list of points
std::vector<long> get_x_coordinates(const std::vector<std::vector<long>> &list_of_points) {
    std::vector<long> x_coordinates;
    for (auto point : list_of_points) {
        x_coordinates.push_back(point[0]);
    }
    return x_coordinates;
}

std::vector<long> get_y_coordinates(const std::vector<std::vector<long>> &list_of_points) {
    std::vector<long> y_coordinates;
    for (auto point : list_of_points) {
        y_coordinates.push_back(point[1]);
    }
    return y_coordinates;
}

long l1_distance(const std::vector<long> a, const std::vector<long> b) {
    assert(a.size() == b.size());
    long dist = 0;
    for (size_t i = 0; i < a.size(); i++) {
        dist += labs(a[i] - b[i]);
    }
    return dist;
}

// --------------------------------------------------
// Functions for estimating netlength
// --------------------------------------------------

long bounding_box_length(const std::vector<std::vector<long>> &list_of_points) {
    long xmin = list_of_points[0][0];
    long xmax = list_of_points[0][0];
    long ymin = list_of_points[0][1];
    long ymax = list_of_points[0][1];

    // Find minimum and maximum
    for (auto point: list_of_points) {
        if (point[0] < xmin) {
            xmin = point[0];
        }
        if (point[0] > xmax) {
            xmax = point[0];
        }
        if (point[1] > ymax) {
            ymax = point[1];
        }
        if (point[0] < ymin) {
            ymin = point[1];
        }
    }
    long length = xmax - xmin + ymax - ymin;
    return length;
}

double clique_netlength(const std::vector<std::vector<long>> &list_of_points) {

    long clique_netlength = 0;

    std::vector<long> x_coordinates = get_x_coordinates(list_of_points);
    std::vector<long> y_coordinates = get_y_coordinates(list_of_points);
    // Sort coordinates independently
    std::sort(x_coordinates.begin(), x_coordinates.end());
    std::sort(y_coordinates.begin(), y_coordinates.end());
    // std::cout << "x_coordinates sorted:" << x_coordinates << std::endl;
    // std::cout << "y_coordinates sorted:" << y_coordinates << std::endl;

    // The first segment is counted (n-1)-times
    // the second segment is counted (n-2)-times
    // ...
    // the last segment is counted 1 time
    for (size_t i = 1; i < x_coordinates.size(); i++) {
        clique_netlength += (x_coordinates.size() - i) * (x_coordinates[i] - x_coordinates[i - 1]);
        clique_netlength += (y_coordinates.size() - i) * (y_coordinates[i] - y_coordinates[i - 1]);
    }

    return clique_netlength / (double) (list_of_points.size() - 1);
}

long star_netlength(const std::vector<std::vector<long>> &list_of_points) {

    // Can deal with x- and y-direction separately
    long star_netlength = 0;
    long starx;
    long stary;

    std::vector<long> x_coordinates = get_x_coordinates(list_of_points);
    std::vector<long> y_coordinates = get_y_coordinates(list_of_points);
    // Sort coordinates independently
    std::sort(x_coordinates.begin(), x_coordinates.end());
    std::sort(y_coordinates.begin(), y_coordinates.end());

    // Mean index is the point with the lowest total distance to
    // all other point
    size_t l = (list_of_points.size() - 1) / 2;
    starx = x_coordinates[l];
    stary = y_coordinates[l];

    // Add distances to the star point
    for (size_t i = 0; i < x_coordinates.size(); i++) {
        star_netlength += labs(starx - x_coordinates[i]);
        star_netlength += labs(stary - y_coordinates[i]);
    }

    return star_netlength;
}

size_t get_min_index(const std::vector<long> &key, const std::vector<bool> &included_in_MST){
    size_t min_index = 0;
    long minimum = LONG_MAX;

    for(size_t i = 0; i < key.size(); i++){
        // Entry only counts if the vertix is not in the tree already
        if(included_in_MST[i] == false &&
                key[i] < minimum){
            min_index = i;
            minimum = key[i];
        }
    }

    return min_index;
}

long MST_length(const std::vector<std::vector<long>> &list_of_points) {
    // We are using Prim's algorithm implemented with an adjacency matrix
    // to achive the runtime of O(n^2), where n = #V(G) is the number of vertices
    const size_t num_vertices = list_of_points.size();
    assert(num_vertices > 0);

    std::vector<bool> included_in_MST; // remembers if the vertex is already connected to the MST
    std::vector<long> key; // Weights of the edges in the cut, we are going to pick the minimum element later
    long MST_len = 0;

    // Initialize the values
    for (size_t i = 0; i < num_vertices; i++) {
        included_in_MST.push_back(false); // No vertex visited
        key.push_back(LONG_MAX);
    }

    // First vertex in list is starting vertex (thus distance zero)
    key[0] = 0;

    // Need num_vertices-1 steps to construct a spanning tree
    for (size_t k = 0; k <= num_vertices - 1; k++) {
        // Get vertex not included in tree with shortest distance
        size_t min_position = get_min_index(key, included_in_MST);
        //std::cout << "min_position: " << min_position << std::endl;

        // Mark as visited
        included_in_MST[min_position] = true;
        // Add length of new edge to length of MST
        MST_len += key[min_position];

        // Update the keys of all the adjacent (here: adjacent means all, because we are on the complete graph)
        // vertices of the new vertex in the MST which are not in the tree yet
        for(size_t i = 0; i < num_vertices; i++){
            if(included_in_MST[i] == false &&
               l1_distance(list_of_points[i], list_of_points[min_position]) < key[i]){
                key[i] = l1_distance(list_of_points[i], list_of_points[min_position]);
            }
        }

//        std::cout << "included_in_MST: " << included_in_MST << std::endl;
//        std::cout << "key: " << key << std::endl;
//        std::cout << "MST_len temp: " << MST_len << std::endl;

    }

    return MST_len;
}

//long algorithm_length(const std::vector<std::vector<long>> &list_of_points) {
//    //TODO
//    return -1;
//}

// --------------------------------------------------
// main
// --------------------------------------------------

int main(int argc, char* argv[]) {
    // To silence the compiler warning about unused variables
    (void)argc;
    (void)argv;

    // does main stuff
    // Read in coordinates from STDIN
    std::vector<std::vector<long>> list_of_points;

    while (!std::cin.eof()) {
        std::string line;
        std::getline(std::cin, line);

        if (std::cin.fail()) {
            // Error or EOF
            //std::cout << "An error occured while reading in" << std::endl;
            break;
        }

        std::istringstream stream(line);
        std::vector<long> point;
        int num_x;
        stream >> num_x;
        int num_y;
        stream >> num_y;
        point.push_back(num_x);
        point.push_back(num_y);
        list_of_points.push_back(point);
    }
    // list_of_points contains all the points that were read in

    std::cout << "Read in the following coordinates" << std::endl;
    std::cout << list_of_points << std::endl;

    // Calculate results
    long bounding_box_len = bounding_box_length(list_of_points);
    double clique_len = clique_netlength(list_of_points);
    long star_len = star_netlength(list_of_points);
    long MST_len = MST_length(list_of_points);
//    long algorithm_len = algorithm_length(list_of_points); //TODO
    long algorithm_len = -1;

    std::cout << "Bounding box: " << bounding_box_len << std::endl;
    std::cout << "Clique: " << clique_len << std::endl;
    std::cout << "Star: " << star_len << std::endl;
    std::cout << "MST: " << MST_len << std::endl;
    std::cout << "Algorithm: " << algorithm_len << std::endl;

    return 0;

}