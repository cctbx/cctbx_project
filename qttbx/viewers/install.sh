#!/bin/bash

# Expects to run from qttbx/viewers/install.sh Run as 'sh install.sh'

# Get the path to the Conda base environment
conda_base=$(conda info --base)
echo "Conda base path: $conda_base"

# Define the path to your specific Conda environment
conda_env_path="$PWD/../../../../conda_base"
echo "Conda env path: $conda_env_path"

# Ensure conda command is available by sourcing the conda.sh script from the Conda base installation
source "$conda_base/etc/profile.d/conda.sh"

# Activate the Conda environment
conda activate "$conda_env_path"

# Check if the environment is activated
if [ $? -ne 0 ]; then
  echo "Failed to activate conda environment: $conda_env_path"
  exit 1
fi
echo "Conda environment activated. Running subsequent commands..."

# Install additional packages
mamba -c conda-forge install pyside2 nodejs qtawesome ipykernel qt-webengine qtconsole-base
# conda -c conda-forge install pyside2 nodejs qtawesome ipykernel qt-webengine qtconsole-base


# Get github repos

# Clone Phenix molstar
cd molstar
git clone https://github.com/phenix-project/phenix-molstar.git

# Compile phenix molstar to javascript
cd ../molstar/phenix-molstar
npm install
npm run build

