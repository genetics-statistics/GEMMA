#!/usr/bin/env bash

# download shunit2 in order to run tests (see INSTALL.md)

cd test
./test_suite.sh | tee /dev/stderr | grep -q 'success rate: 100%'
