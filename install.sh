#!/bin/bash

# Exit on error
set -e

# Ensure you are in the project root directory
PROJECT_ROOT=$(dirname "$0")
cd "$PROJECT_ROOT"

# Create and navigate to the build directory
echo "Creating build directory..."
mkdir -p build
cd build

# Compile Fortran sources into a Python extension module using f2py
echo "Compiling Fortran sources..."
f2py -c ../lib/_dlsodes.pyf --opt='-std=legacy -O3' -m dlsodes ../lib/blkdta000.f ../lib/opkda1.f ../lib/opkda2.f ../lib/dlsodes.f --quiet

# Move the compiled extension module to the Python package directory
# echo "Moving compiled module to the Python package directory..."
mv *.so ../

# Navigate back to the project root directory
cd ..

# Find the site-packages directory of the current Python environment
SITE_PACKAGES=$(python -c "import site; print(site.getsitepackages()[0])")

# Create the destination directory if it doesn't exist
DEST_DIR="$SITE_PACKAGES/pylsodes"
mkdir -p "$DEST_DIR"

# Copy the Python package to the site-packages directory
echo "Copying Python package to site-packages..."
cp -r ./__init__.py "$DEST_DIR/"
cp -r ./*.so "$DEST_DIR/"

echo "Installation complete!"
