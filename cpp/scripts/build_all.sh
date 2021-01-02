#!/usr/bin/env bash

# Uses environment variable CMAKE_BIN

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJ_DIR=$DIR/..


MODE="release"

while getopts ":d" opt; do
    case $opt in
        d)
            MODE="debug"
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            ;;
    esac
done



set -ex

cd $PROJ_DIR


if [[ -z "${CMAKE_BIN}" ]]; then
  SCRIPT_CMAKE_BIN="cmake"
else
  SCRIPT_CMAKE_BIN="${CMAKE_BIN}"
fi

if [ $MODE = "debug" ]; then
    mkdir -p cmake-build-debug
    ${SCRIPT_CMAKE_BIN} -H. -Bcmake-build-debug -DCMAKE_BUILD_TYPE=Debug
    make -Ccmake-build-debug -j 10
elif [ $MODE = "release" ]; then
    mkdir -p cmake-build-release
    ${SCRIPT_CMAKE_BIN} -H. -Bcmake-build-release -Dtest=OFF -DDEBUG=OFF -DCMAKE_BUILD_TYPE=Release
    make -Ccmake-build-release -j 10
fi

echo "OK"