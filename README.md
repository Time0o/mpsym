# Installation #

### Dependencies ###

* A compiler with _C++11_ support
* _cmake_ version 3.0 or up
* _google test_ (automatically downloaded as part of the build process)

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

It is also possbile to obtain debug output when running a single test
by passing `-v`, `-vv` or `-vvv` as command line options (resulting in
increasing debug output verbosity). e.g:

```
./test/perm_test -vvv
```
