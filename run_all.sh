#!/bin/bash

for filename in ./SteinerInstances/*.txt; do
    echo $filename
    cat $filename | ./cmake-build-debug/chipdesign
    echo "-------------------------------"
done
