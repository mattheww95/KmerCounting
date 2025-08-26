# K-mer Counting

## Requirments

Requires:

    - zlib

    - CMake

    - GCC

    - make

## Installation

Create build directory

`cd` into build directory and run `cmake .. -DTARGET_GROUP=release` to compile a relase build. If you need to specify your zlib library location you can do so like so `-DZLIB_LIBRARY="/usr/local/lib/zlib/lib/libz.a"`

To make a debug build run: `cmake .. -DTARGET_GROUP=debug`


To run tests build run: `cmake .. -DTARGET_GROUP=test`

After running any of the cmake options, simply run `make` in the build directory


Release build will be generated in `build/src/kmer-compressor`



## Usage


```
Welcome to our slow k-mer counting program!

Options:
  [REQUIRED] -i|--input:
        Input file.
  [REQUIRED] -s|--size:
        Kmer size.
  [REQUIRED] -o|--output:
        Output file.
  -h|--help:
        Show help and exit.

```


## Purpose

Just trying to remain competent in C.
