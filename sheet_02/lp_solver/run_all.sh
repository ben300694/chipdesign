#!/bin/bash

for filename in ./instances/*.dimacs; do
    echo "----------------------------------------"  
    echo "Calling program on " $filename
    ./msscxx $filename
    echo "----------------------------------------"
done
