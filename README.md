# Papillon Nuclear Data Library
[![CMake Workflow](https://github.com/HunterBelanger/papillon-ndl/actions/workflows/cmake.yml/badge.svg)](https://github.com/HunterBelanger/papillon-ndl/actions)
[![Documentation Status](https://readthedocs.org/projects/papillon-ndl/badge/?version=latest)](https://papillon-ndl.readthedocs.io/en/latest/?badge=latest)
[![License](https://img.shields.io/badge/License-GPLv3-brightgreen)](https://github.com/HunterBelanger/papillon-ndl/blob/develop/LICENSE)

The Papillon Nuclear Data Library (NDL) is a C++/Python library for reading,
sampling, and interacting with continuous energy nuclear data, stored in the ACE
nuclear data format.

For examples of how to use the library in both C++ and Python, take a look at
the [documentation site](https://papillon-ndl.readthedocs.io). That is where you
will also be able to fined more detailed installation instructions.

## License
PapillonNDL is provided under the terms and conditions of the GPLv3
license.

## Dependencies
The library may be built on Unix-like operating systems or on Windows. All that
is required is cmake >= 3.11, and a C++ compiler which supports the C++20
standard. For Unix-like systems, the recommended compilers are GCC >= 11.
On Windows, you should have MSVC >= 19.29. In order to build the
Python interface, Python >= 3.5 should be installed on your system, in
addition to the Python development libraries and header files.

Tests are not built by default, and should only be needed for developers. You
can turn them on by using ```-DPNDL_TESTS=ON``` with cmake.

## Install
To build PapillonNDL, navigate to the directory where you would like to keep the
source files, and then run the following commands:
```
$ git clone https://github.com/HunterBelanger/papillon-ndl.git
$ cd papillon-ndl
$ cmake -E make_directory build
$ cd build
$ cmake -DCMAKE_BUILD_TYPE=Release ..
$ cmake --build . --target install
```

If you do NOT want to build the Python bindings for PapillonNDL, then you should
add the flag ```-DPNDL_PYTHON=OFF``` to the cmake command.
