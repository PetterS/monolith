#!/bin/bash
# Run this in a build directory with GCC and CMAKE_BUILD_TYPE=Debug.
# It will generate a report in a directory called "html".
set -e

# Delete old temp files.
lcov --zerocounters --directory .
rm *.info || true

# Run all tests. This will capture data.
ctest

# Generate the report.
lcov --capture --directory . --output-file coverage.info
lcov --output $(pwd)/coverage2.info --extract coverage.info "**/minimum/**"
lcov --output $(pwd)/my_coverage.info --remove coverage2.info "*_test.cpp"
genhtml my_coverage.info --output-directory html
