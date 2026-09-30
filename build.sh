#!/bin/sh

set -e
cmake=cmake
type cmake3 >/dev/null 2>&1 && cmake=cmake3

main_dir=$(dirname `readlink -f $0`)
build_dir=$main_dir/.build

mkdir -p $build_dir
cd $build_dir
$cmake .. -DCMAKE_INSTALL_PREFIX=$main_dir \
       -DCMAKE_INSTALL_RPATH_USE_LINK_PATH="ON"
$cmake --build . -j `nproc`
$cmake --install .
