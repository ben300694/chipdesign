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

std::vector<long> get_x_coordinates(const std::vector<std::vector<long>> &list_of_points){
    std::vector<long> x_coordinates;
    for(auto point : list_of_points){
        x_coordinates.push_back(point[0]);
    }
    return x_coordinates;
}

std::vector<long> get_y_coordinates(const std::vector<std::vector<long>> &list_of_points){
    std::vector<long> y_coordinates;
    for(auto point : list_of_points){
        y_coordinates.push_back(point[1]);
    }
    return y_coordinates;
}

// --------------------------------------------------
// Functions for estimating netlength

long bounding_box(const std::vector<std::vector <long>> listofpoints){
    long xmin=listofpoints[0][0];
    long xmax=listofpoints[0][0];
    long ymin=listofpoints[0][1];
    long ymax=listofpoints[0][1];

    // Find minimum and maximum
    for(auto point: listofpoints){
        if (point[0]<xmin){
            xmin=point[0];
        }
        if (point[0]>xmax){
            xmax=point[0];
        }
        if (point[1]>ymax){
            ymax=point[1];
        }
        if (point[0]<ymin){
            ymin=point[1];
        }
    }
    long length = xmax-xmin+ymax-ymin;
    return length;
}

double clique_netlength(const std::vector<std::vector<long>> &list_of_points){

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
    for(size_t i = 1; i < x_coordinates.size(); i++){
        clique_netlength += (x_coordinates.size() - i)*(x_coordinates[i] - x_coordinates[i-1]);
        clique_netlength += (y_coordinates.size() - i)*(y_coordinates[i] - y_coordinates[i-1]);
    }

    return clique_netlength/(double)(list_of_points.size() - 1);
}

long star_netlength(const std::vector<std::vector<long>> &list_of_points){

    long star_netlength = 0;
    long starx;
    long stary;

    std::vector<long> x_coordinates = get_x_coordinates(list_of_points);
    std::vector<long> y_coordinates = get_y_coordinates(list_of_points);
    // Sort coordinates independently
    std::sort(x_coordinates.begin(), x_coordinates.end());
    std::sort(y_coordinates.begin(), y_coordinates.end());
    size_t l=(list_of_points.size()-1)/2;
    starx = x_coordinates[l];
    stary = y_coordinates[l];

    for(size_t i = 0; i < x_coordinates.size(); i++){
        star_netlength += labs(starx - x_coordinates[i]);
        star_netlength += labs(stary - y_coordinates[i]);
    }

    return star_netlength;
}



int main(int argc, char** argv) {
    //does main stuff
    // Read in coordinates from STDIN
    std::vector<std::vector<long>> list_of_points;

    while (!std::cin.eof()) {
        std::string line;
        std::getline(std::cin, line);

        if (std::cin.fail()) {
            //error
            std::cout << "An error occured while reading in" << std::endl;
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

    std::cout << "Read in the following coordinates" << std::endl;
    std::cout << list_of_points << std::endl;

    // Calculate results
    long bounding_box_netlength = bounding_box(list_of_points);
    double clique_length = clique_netlength(list_of_points);
    long star_length = star_netlength(list_of_points);

    std::cout << "Bounding box: " << bounding_box_netlength << std::endl;
    std::cout << "Clique: " << clique_length << std::endl;
    std::cout << "Star: " << star_length << std::endl;

    return 0;

}