#!/bin/bash

# Print starting message
echo "Starting installation script..."

# Clone the repository if it does not exist
REPO_URL="git@github.com:phenix-project/phenix-molstar"
PROJECT_DIR="$(pwd)/phenix-molstar"

if [ -d "$PROJECT_DIR" ]; then
  echo "Project directory already exists. Skipping git clone."
else
  echo "Cloning the repository..."
  git clone $REPO_URL
  if [ $? -ne 0 ]; then
    echo "Error: Failed to clone the repository."
    exit 1
  fi
fi

# Ensure Node.js and npm from the active conda environment are used
if [ -z "$CONDA_PREFIX" ]; then
  echo "Error: No active conda environment found."
  exit 1
fi

echo "Using Node.js and npm from the active conda environment: $CONDA_PREFIX"

# Verify installation
echo "Verifying Node.js installation..."
node -v
if [ $? -ne 0 ]; then
  echo "Error: Node.js installation verification failed. Install nodejs with conda."
  exit 1
fi

echo "Verifying npm installation..."
npm -v
if [ $? -ne 0 ]; then
  echo "Error: npm installation verification failed. Install nodejs with conda."
  exit 1
fi

# Change to project directory
cd $PROJECT_DIR

# Install dependencies and build the project
echo "Installing project dependencies..."
npm install
if [ $? -ne 0 ]; then
  echo "Error: Failed to install project dependencies."
  exit 1
fi

echo "Building the project..."
npm run build
if [ $? -ne 0 ]; then
  echo "Error: Project build failed."
  exit 1
fi

echo "Installation script completed successfully."