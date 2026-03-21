#!/bin/bash

cmake -S . -B build
cd build
make
echo ""
./bin/MasterThesis
cd ..