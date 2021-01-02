#!/bin/bash

# See https://wiki.wxwidgets.org/Valgrind_Suppression_File_Howto

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

cd "${DIR}"/..

set -ex


COMMAND_NO_OP="python scripts/run_meshgrid_debug.py no-op"
COMMAND_TO_TEST="python scripts/run_meshgrid_debug.py"


valgrind --leak-check=full --show-reachable=yes --error-limit=no --gen-suppressions=all --log-file=/tmp/minimalraw.log ${COMMAND_NO_OP}
cat /tmp/minimalraw.log | /usr/bin/awk -f ./scripts/parse_valgrind_suppressions.sh > /tmp/minimal_mesh_to_depth.supp

echo "Running valgrind with suppressions"

valgrind --suppressions=/tmp/minimal_mesh_to_depth.supp ${COMMAND_TO_TEST}