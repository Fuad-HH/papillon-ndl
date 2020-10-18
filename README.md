# Papillon Nuclear Data Library
![License](https://img.shields.io/badge/License-CeCILL%20v2.1-brightgreen)

The Papillon Nuclear Data Library (NDL) is used for reading, sampling,
and interacting with continuous energy nuclear data, stored in the ACE
nuclear data format.

## Dependencies
To build and install the library, cmake >= 3.9 is required, along with a C++
compiler which supports the C++17 standard. The recommended compiler is
Clang >= 5, though GCC >= 6 should suffice.

Building the unit tests (using ```-DPNDL_TESTS=ON``` when calling cmake)
requires that Google test already be installed on the system, and is not
provided. Tests are not built by default, and therefore this is not required
for most users.

## Install

## License
PapillonNDL is provided under the terms and conditions of the CeCILLv2.1
license. This is a French equivalent of the GPLv3 license, and is explicitly
compatible with both the GPLv2 and v3. The French version of this license is
provided in the LICENSE file, along with the equaly valid English version, which
may be found in the LICENSE-ENGLISH file. The CeCILLv2.1 is approved by the FSF.
For more information about this license may be found [here](https://cecill.info/).
