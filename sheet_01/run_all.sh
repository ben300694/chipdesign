#!/bin/bash

for filename in ./SteinerInstances/*.txt; do
    echo "----------------------------------------"  
    echo "Calling program on " $filename
    cat $filename | ./cmake-build-debug/chipdesign
    echo "----------------------------------------"
done
