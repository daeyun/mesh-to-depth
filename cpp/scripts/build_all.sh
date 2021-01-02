#!/usr/bin/env bash

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

if [ $MODE = "debug" ]; then
    mkdir -p cmake-build-debug
    cmake -H. -Bcmake-build-debug -DCMAKE_BUILD_TYPE=Debug
    make -Ccmake-build-debug -j 10
elif [ $MODE = "release" ]; then
    mkdir -p cmake-build-release
    cmake -H. -Bcmake-build-release -Dtest=OFF -DDEBUG=OFF -DCMAKE_BUILD_TYPE=Release
    make -Ccmake-build-release -j 10
fi

echo "OK"