### Installing MSYS2

See [README_MSYS2](../README_MSYS2.md)

### Installing Intel Compiler (when using MKL)

Install the Intel Compiler by running Intel Parallel Studio XE 2019 installer.

### Getting Quantum ESPRESSO (QE)

Download version 7.1 of Quantum ESPRESSO (QE) from the [QE download page](https://www.quantum-espresso.org/download-page/). Extract the downloaded file and execute the following commands:

```shell
$ tar xvfz qe-7.1-ReleasePack.tar.gz
$ patch -p0 < qe-7.1.patch
$ cd qe-7.1
```

Then, run the `./configure` command to configure QE. To use MKL, edit the `make.inc` file:

```shell
DFLAGS         = -D_WIN32 -D__DFTI -D__LIBXC -D__MPI
IFLAGS         = -I. -I$(TOPDIR)/include -I$(TOPDIR)/FoX/finclude  -I/mingw64/include -I$(TOPDIR)/mkl/include
BLAS_LIBS      = -L$(TOPDIR)/mkl/lib/intel64 -lmkl_rt
LAPACK_LIBS    = -L$(TOPDIR)/mkl/lib/intel64 -lmkl_rt
FFT_LIBS       = -L$(TOPDIR)/mkl/lib/intel64 -lmkl_rt
```

Next, build QE using the following command:

```shell
$ rm external/lapack/.appveyor.yml
$ rm external/lapack/.travis.yml
$ cp -r /c/Program\ Files\ \(x86\)/IntelSWTools/compilers_and_libraries/windows/mkl/ mkl
$ make all
```

Copy the required DLL files:

```shell
$ cp /c/Program\ Files\ \(x86\)/IntelSWTools/compilers_and_libraries/windows/redist/intel64_win/mkl/*.dll bin/
$ cp /c/Program\ Files\ \(x86\)/Common\ Files/Intel/Shared\ Libraries/redist/intel64_win/compiler/*.dll bin/
```

Run the tests:

```shell
$ cd test-suite
$ export PATH=$PATH:/c/Program\ Files/Microsoft\ MPI/Bin
$ make run-tests
$ cd ..
```

Move to the PW/example directory and run the provided examples:

```shell
$ cd PW/example
$ bash run_all_example
```

### Building average.exe

To avoid bugs, build `average.exe` without MPI and MKL using the following steps:

```shell
$ tar xvfz qe-7.1-ReleasePack.tar.gz
$ patch -p0 < qe-7.1.patch
$ cd qe-7.1
$ ./configure ARCH=mingw64 CFLAGS="-Dsrandom=srand -Drandom=rand" --with-libxc --with-libxc-prefix=/mingw64 --disable-parallel
$ rm external/lapack/.appveyor.yml
$ rm external/lapack/.travis.yml
$ make pp
```

### Creating a Package

Create a package using the following steps:

```shell
$ mkdir -p dist/pseudo
$ cp -r bin dist
$ cd dist/bin
$ strip *x
$ for f in *.x; do mv $f ${f%.x}.exe; done
$ cp /mingw64/bin/libgfortran-5.dll .
$ cp /mingw64/bin/libgcc_s_seh-1.dll .
$ cp /mingw64/bin/libquadmath-0.dll .
$ cp /mingw64/bin/libwinpthread-1.dll .
$ cp /c/Windows/System32/msmpi.dll .
$ cp /c/Program\ Files\ \(x86\)/IntelSWTools/compilers_and_libraries/windows/redist/intel64_win/mkl/*.dll .
$ cp /c/Program\ Files\ \(x86\)/Common\ Files/Intel/Shared\ Libraries/redist/intel64_win/compiler/*.dll .
```

Copy the pseudopotential files in the official package to the `pseudo/` folder. Specifically, copy `Ti.pw-mt_fhi.UPF` to the `pseudo/` folder.

Copy the `MOLs` folder from [https://github.com/nisihara1/MOLs](https://github.com/nisihara1/MOLs).

