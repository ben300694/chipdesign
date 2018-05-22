//
// Created by ben on 22.05.18.
//

#include <iostream>
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

int main(int argc, char** argv) {

    std::vector<std::vector<long>> list_of_points;

    while (!std::cin.eof()){
        std::string line;
        std::getline(std::cin, line);

        if (std::cin.fail()){
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

   std::cout << list_of_points;

}