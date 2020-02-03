# libgmultigrid
A header-only library for geometric multigrid.

## Installation

The recommended way to use this library is by adding it as a submodule to the project that you would like to use it in. Navigate to the directory you would like to clone the submodule into, and then run:
```
git submodule add https://github.com/icethrush/libgmultigrid.git
```
After this, add `libgmultigrid/include` to the include path of your project, and you should be able to use it in your code.

## Dependencies

This library uses Eigen, and thus requires the Eigen directories to be in the include path of your project or machine.

This library has only been included in projects that were built using Clang. Usage with gcc/g++ is possible, but the compiler may emit a different set of warnings that may need to be suppressed (if, for instance, you are compiling with -Wall).


