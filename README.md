# K-mer Counting

## Requirments

Requires:
    - zlib
    - CMake
    - GCC

## Installation

Create build directory

`cd` into build directory and run `cmake .. -DTARGET_GROUP=release` to compile a relase build. If you need to specify your zlib library location you can do so like so `-DZLIB_LIBRARY="/usr/local/lib/zlib/lib/libz.a"`

To make a debug build run: `cmake .. -DTARGET_GROUP=debug`


To run tests build run: `cmake .. -DTARGET_GROUP=test`


Release build will be generated in `build/src/kmer-compressor`

