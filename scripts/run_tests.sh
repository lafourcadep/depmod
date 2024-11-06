#!/usr/bin/bash

# Ensure bash v >= 4. is used.
if (( BASH_VERSINFO[0] < 4 )); then
    echo "${__name__} requires bash v4 or greater"
    echo "Current bash version: ${BASH_VERSION}"
    exit 1
fi

# Get the root path of the project
ROOT=$(cd "$(cd "$(dirname "${BASH_SOURCE[0]}")" &> /dev/null && pwd)/.." &> /dev/null && pwd)

# Check if pytest is available in the path
if ! [ -x "$(command -v pytest)" ]; then
    echo "pytest is not available"
    exit 1
fi

# Go to root folder and execute the tests
cd ${ROOT} && pytest
