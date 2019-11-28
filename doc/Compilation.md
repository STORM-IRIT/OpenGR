# Compilation  {#compilation}
OpenGR is supported on Linux, MacOS and Windows platforms, with continuous testing.
The library is now header only, and distributed with a cmake package.
Demo applications and tests are also provided

This page reviews the following aspects:
* Build from source: how to install the library, compile the tests and install the applications: [here](#library)
* CMake package: how to use the library in your own software: [here](#usage)
* a note about compilation mode and performances: [here](#perfs)

## <a name="library"></a> Compiling the library, applications and tests
For people in a hurry:
```bash
git clone https://github.com/STORM-IRIT/OpenGR.git
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=./install
make install
cd install/scripts/
./run-example.sh
```

### Dependencies
OpenGR takes care of its own dependencies, and downloads them when required. As a result, **you may need an internet connection for the first compilation**.

The libraries and the standalone application require:
* [Eigen](http://eigen.tuxfamily.org/). The library is automagically downloaded and configured by the cmake process. If you need to force using Eigen from a specific path, call CMAKE with `-DEIGEN3_INCLUDE_DIR=/path/to/eigen`.

IO and tests:
* [Boost-filesystem](http://www.boost.org/doc/libs/1_57_0/libs/filesystem/doc/index.htm) (optionnal), used to read dataset folders. Requires Boost version: 1.57 or more.

### CMake Options
Our CMake scripts can be configured using the following options:
* Compilation rules
    * `OPTION (OpenGR_COMPILE_TESTS "Enable testing" TRUE)`
    * `OPTION (OpenGR_COMPILE_APPS "Compile demo applications (including the Super4PCS standalone)" TRUE)`
    * `OPTION (IO_USE_BOOST "Use boost::filesystem for texture loading" TRUE)`
    * `OPTION (OpenGR_USE_CHEALPIX "Use Chealpix for orientation filtering (deprecated)" FALSE)` We recommend to keep this option to `FALSE`.
* Advanced functionalities
    * `+OPTION (OpenGR_USE_WEIGHTED_LCP "Use gaussian weights for point samples when computing LCP" FALSE)` Was implicitely set to `FALSE` in previous release.
* Extras
    * `OPTION (DL_DATASETS "Download demo datasets and associated run scripts" FALSE)`
* Debug options
    * `OPTION (ENABLE_TIMING "Enable computation time recording" FALSE)`

Options can be set by calling `cmake -DMY_OPTION=TRUE`, or by editing the file [CMakeList.txt](https://github.com/nmellado/Super4PCS/blob/master/CMakeLists.txt) (not recommended).

### Compilation targets
Our CMAKE scripts provide several targets:
* `all`: build the libraries and the demos, if enabled with `OpenGR_COMPILE_APPS`,
* `install`: build and install the libraries, the cmake package, and if requested the demo applications. This is **recommended** target.
It can be customized by calling cmake with `-DCMAKE_INSTALL_PREFIX=your_path`.
* `OpenGR-PCLWrapper`/`Super4PCS`: build the demos application,
* `buildtests`: compile tests. This target is generated only if `OpenGR_COMPILE_TESTS` option is set to `TRUE`,
* `test`: run the tests locally
* `doc`: run doxygen to compile the documentation
* `dl-datasets`: download the demo datasets.

### Installed directories structure:
After compilation and installation, you should obtain the following structure:
```
installation-directory/
  - assets     -> contains the demo hippo files, and the downloaded datasets if any
  - bin        -> contains the main Super4PCS application
  - include    -> header files of the Super4PCS library
  - lib        -> library binaries: super4pcs_algo and super4pcs_io, and chealpix if compiled
  - lib/cmake  -> cmake package files
  - scripts    -> demo scripts allowing to register the hippo meshes
  - doc        -> this documentation, compiled using Doxygen
```

***

### Linux and MacOS Builds
Nothing specific here, just follow the [aforementionned](#library) instructions.

***

### Windows Builds
OpenGR requires c++11 features and CMake support.
We recommend to use Microsoft Visual Studio 2019 (2017 should also work) with [CMake support](https://docs.microsoft.com/en-us/cpp/build/cmake-projects-in-visual-studio?view=vs-2019).
CMake standalone is also supported (and used for Continuous Integration on AppVeyor), but not recommended.

By default, the project can be opened and compiled straight away with no parameter setting.

If you want to compile with Boost support, the dependency directories must be properly configured in Visual Studio.
They can be set either when calling CMake or globally by configuring Visual Studio (see how-to [here](https://social.msdn.microsoft.com/Forums/vstudio/en-US/a494abb8-3561-4ebe-9eb0-6f644a679862/visual-studio-2010-professional-how-to-add-include-directory-for-all-projects?forum=vcgeneral#7b5ab5f2-f793-4b0e-a18a-679948d12bdd)).
![](img/VStudio-globalFolders.jpg)

***

## <a name="usage"></a> Use OpenGR library in your own application
**From release v.1.2**, OpenGR is now header-only, and provides a CMake package generated during the installation process.
The target namespace is `gr`.

To use it in your own application, add the following lines to your project file:
```cmake
project(OpenGR-externalAppTest)
cmake_minimum_required(VERSION 3.3)

set (CMAKE_CXX_STANDARD 11)

find_package(OpenGR REQUIRED)
find_package( Eigen3 REQUIRED )

add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} gr::algo Eigen3::Eigen)
```

In addition, CMake must be ran so that `CMAKE_PREFIX_PATH` contains `OpenGR_install_dir/lib/cmake` when compiling `OpenGR-externalAppTest`.

OpenGR files will be located in the `gr` folder. For instance, to use the Super4PCS algorithm:
```cpp
#include <gr/algorithms/match4pcsBase.h>
#include <gr/algorithms/FunctorSuper4pcs.h>

```
This functionality is tested by our continuous integration system. Checkout the [externalAppTest](https://github.com/STORM-IRIT/OpenGR/tree/master/tests/externalAppTest) for a working example.

## <a name="perfs"></a> Debug mode and performances
Note that we heavily use template mechanisms that requires to enable inlining in order to be efficient. Compiling in Debug mode without inlining may result in longer running time than expected.
