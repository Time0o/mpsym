[![GitHub Releases](https://img.shields.io/github/release/Time0o/mpsym.svg)](https://github.com/Time0o/mpsym/releases)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Documentation](https://img.shields.io/badge/docs-doxygen-blue.svg)](https://Time0o.github.io/mpsym/doxygen)
[![Build Status](https://travis-ci.com/Time0o/mpsym.svg?branch=master)](https://travis-ci.com/Time0o/mpsym)
[![Coverage
status](https://codecov.io/gh/Time0o/mpsym/branch/master/graph/badge.svg)](https://codecov.io/gh/Time0o/mpsym?branch=master)

# MPSoC Symmetry Reduction Package #

This is an ongoing effort to implement a fast C++ library with a small set of
dependencies which uses methods from computational group theory in order to
detect symmetries in abstract MPSoC hardware descriptions. These symmetries
can then in principle be exploited by e.g. compiler pipelines employing
automatic task to processor mapping to very efficiently detect equivalent
mappings and thus avoid expensive re-estimation of performance characteristics.
This project also aims to provide an intuitive API specifically for this
purpose.

## Documentation ##

The documentation can be built from source as described below and is also
hosted [here](https://time0o.github.io/mpsym/).

## Installation ##

(Note that this project is in the early stages of it's development and thus
not yet ready to be used in another application, detailed build instructions
for this case will follow in the future)

### Dependencies ###

Compiling and running the tests should only require the following:
* A compiler (preferably _gcc_) with C++11 support
* _CMake_ version 3.0 or up
* A current version of the _Boost_ libraries
* Lua 5.3.*
* _Google Test_ (automatically downloaded as part of the build process)
* _Nauty Traces_ (automatically downloaded as part of the build process)

Generating code coverage (only supported for _gcc_) also requires:
* The `gcov`, `lcov` and `genhtml` command line tools

Executing profiling demos (only supported on Linux) also requires:
* `perf`

Building the documentation also requires:
* _Doxygen_
* _Python 3 + Breathe + Sphinx_ (optional)

### Compiling and Running Tests ###

Compile both source and tests with:

```
cd debug
cmake -DCMAKE_BUILD_TYPE=Debug ..
make all
```

And then run all tests using:

```
make test
```

Or a single test with:

```
./test/xxx_test
```

It is also possible to obtain debug output when running a single test
by passing `-v`, `-vv` or `-vvv` as command line options (resulting in
increasing debug output verbosity), e.g.:

```
./test/perm_test -vvv
```

In addition, the `-o` flag can be used to run single testcases within
a test, e.g.:

```
./test/perm_test -o CanConstructPerm
```

Run `./test/xxx_test --help` for more information.

### Building the Documentation ###

(Note that the documentation can also be found
[here](https://time0o.github.io/mpsym/))

Change into the `release` directory and run `cmake -DCMAKE_BUILD_TYPE=Release
..`, followed by `make doxygen` or `make sphinx` (which automatically executes
the former). Generated documentation will be placed in `./doc/doxygen` and/or
`./doc/sphinx`.
