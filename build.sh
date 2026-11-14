#!/bin/bash

set -e

PYTHON=$(which python)

echo "Using Python: $PYTHON"
echo "Python version: $($PYTHON --version)"

echo "Configuring project..."

cmake -S . -B build \
    -DCMAKE_BUILD_TYPE=Release \
    -DPython_EXECUTABLE="$PYTHON" \
    -DPython_FIND_STRATEGY=LOCATION

echo ""
echo "Building Python extension..."

cmake --build build --config Release -j$(nproc)

echo ""
echo "Build completed successfully!"
