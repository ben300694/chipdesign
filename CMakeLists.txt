cmake_minimum_required(VERSION 3.10)
project(chipdesign)

set(CMAKE_CXX_STANDARD 11)

add_executable(chipdesign
        SteinerInstances/rand10.txt
        SteinerInstances/rand100.txt
        SteinerInstances/rand1000.txt
        SteinerInstances/rand200.txt
        SteinerInstances/rand30.txt
        SteinerInstances/rand500.txt
        SteinerInstances/vlsi10.txt
        SteinerInstances/vlsi100.txt
        SteinerInstances/vlsi1000.txt
        SteinerInstances/vlsi20.txt
        SteinerInstances/vlsi200.txt
        SteinerInstances/vlsi30.txt
        SteinerInstances/vlsi300.txt
        SteinerInstances/vlsi40.txt
        SteinerInstances/vlsi400.txt
        SteinerInstances/vlsi50.txt
        SteinerInstances/vlsi500.txt
        SteinerInstances/vlsi60.txt
        SteinerInstances/vlsi600.txt
        SteinerInstances/vlsi70.txt
        SteinerInstances/vlsi700.txt
        SteinerInstances/vlsi80.txt
        SteinerInstances/vlsi800.txt
        SteinerInstances/vlsi90.txt
        SteinerInstances/vlsi900.txt
        LICENSE
        README.md main.cpp)
