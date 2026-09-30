#!/bin/sh

set -e

main_dir=$(dirname `readlink -f $0`)
build_dir=$main_dir/.build
bin_dir=$main_dir/bin

rm -rf $build_dir $bin_dir
