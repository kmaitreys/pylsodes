#!/bin/bash

# Exit on error
set -e

# Ensure you are in the project root directory
PROJECT_ROOT=$(dirname "$0")
cd "$PROJECT_ROOT"

# Check if Python dependencies are installed
echo "Checking dependencies..."

PYTHON_CMD="python"

# Check if numpy is installed, if not install it
if ! $PYTHON_CMD -c "import numpy" &>/dev/null; then
    echo "Installing numpy..."
    $PYTHON_CMD -m pip install --quiet numpy==2.0
fi

# Check if meson is installed, if not install it
if ! $PYTHON_CMD -c "import mesonbuild" &>/dev/null; then
    echo "Installing meson..."
    $PYTHON_CMD -m pip install --quiet meson
fi

# Check if ninja is installed, if not install it
if ! command -v ninja &>/dev/null; then
    echo "Installing ninja..."
    $PYTHON_CMD -m pip install --quiet ninja
fi

echo "All dependencies are installed."

# Create and navigate to the build directory
echo "Creating build directory..."
mkdir -p build
cd build

# Compile Fortran sources into a Python extension module using f2py
echo "Compiling Fortran sources..."
CFLAGS="-O3" FFLAGS="-std=legacy -O3" LDFLAGS="-O3" python -m numpy.f2py -c ../lib/_dlsodes.pyf -m _dlsodes ../lib/blkdta000.f ../lib/opkda1.f ../lib/opkda2.f ../lib/dlsodes.f --quiet --backend="meson"
CFLAGS="-O3" FFLAGS="-std=legacy -O3" LDFLAGS="-O3" python -m numpy.f2py -c ../lib/_dlsode.pyf -m _dlsode ../lib/blkdta000.f ../lib/opkda1.f ../lib/opkda2.f ../lib/dlsode.f --quiet --backend="meson"

# Move the compiled extension module to the Python package directory
# echo "Moving compiled module to the Python package directory..."
# mv *.so ../

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
cp -r build/*.so "$DEST_DIR/"
cp -r ./*.pyi "$DEST_DIR/"

echo "Installation complete!"
