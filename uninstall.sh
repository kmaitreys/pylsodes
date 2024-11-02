#!/bin/bash

# Exit on error
set -e

# Find the site-packages directory of the current Python environment
SITE_PACKAGES=$(python -c "import site; print(site.getsitepackages()[0])")

# Define the destination directory
DEST_DIR="$SITE_PACKAGES/pylsodes"

# Check if the package exists in the site-packages directory
if [ -d "$DEST_DIR" ]; then
    echo "Removing pylsodes package from site-packages..."
    rm -rf "$DEST_DIR"
else
    echo "pylsodes package not found in site-packages."
fi

# Find and remove the compiled .so files in the project root directory if they exist
echo "Removing compiled .so files from the project root..."
find . -maxdepth 1 -name "*.so" -exec rm -f {} \;

echo "Uninstallation complete!"
