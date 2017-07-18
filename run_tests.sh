#!/usr/bin/env bash

# The "tr" command fixes the ^M characters in the output.
cd test
./test_suite.sh 2>&1 | tr '\r' '\n' > test.log
cat test.log | grep -q 'success rate: 100%'
