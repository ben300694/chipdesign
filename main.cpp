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

long clique_netlength(const std::vector<std::vector<long>> &list_of_points){
    long clique_netlength = 0;

    std::vector<long> x_coordinates = get_x_coordinates(list_of_points);
    std::vector<long> y_coordinates = get_y_coordinates(list_of_points);
    std::sort(x_coordinates.begin(), x_coordinates.end());
    std::sort(y_coordinates.begin(), y_coordinates.end());
    std::cout << "x_coordinates sorted:" << x_coordinates << std::endl;
    std::cout << "y_coordinates sorted:" << y_coordinates << std::endl;


    return clique_netlength;
}

int main(int argc, char **argv) {

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

    long bounding_box_netlength = 0; //TODO
    long clique_length = clique_netlength(list_of_points);

    return 0;

}