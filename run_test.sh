#!/bin/bash

cmake -S . -B build
make -C build
echo ""

cd build
ctest --output-on-failure

cd ..