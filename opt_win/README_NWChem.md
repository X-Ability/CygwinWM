### Installing MSYS2

See [README_MSYS2](/README_MSYS2.md)

However, to make it 32-bit, use MSYS2 MinGW x86

### Installing MPICH2

Download and install MPICH2 from [this link](http://www.mpich.org/static/tarballs/1.4.1p1/mpich2-1.4.1p1-win-ia32.msi). If you have CygwinWM installed, ensure it is empty.

NWChem remembers the compile-time directory, so work with the actual installation folder:

```shell
$ mkdir -p /c/cygwin_wm/opt_win
$ cp -rp /c/Program\ Files\ \(x86\)/MPICH2/ /c/cygwin_wm/opt_win/
```

### Obtaining NWChem

```shell
$ cd /c/cygwin_wm/opt_win
```

Download NWChem 7.0.2:

```shell
$ wget https://github.com/nwchemgit/nwchem/archive/refs/tags/v7.0.2-release.tar.gz
$ tar xvfz v7.0.2-release.tar.gz
$ mv nwchem-7.0.2-release/ NWChem_7.0.2
```

### Compiling NWChem

In `NWChem_7.0.2/src/nwpw/nwpwlib/io/compressed_io.c`, comment out the line containing `fsync`.

In `NWChem_7.0.2/src/tools/get-tools-github`, comment out the line:

```shell
./msmpi_patch.sh $GA_DIR
```

If your gfortran major version has two digits, modify `NWChem_7.0.2/src/config/makefile.h`. Change the line inside the conditional block:

```shell
ifeq ($(TARGET),$(findstring $(TARGET),LINUX CYGNUS CYGWIN))
```

From:

```shell
        GNUMINOR=$(shell $(_FC) -dM -E - < /dev/null 2> /dev/null | egrep __GNUC_MINOR | cut -c24)
```

To:

```shell
        GNUMINOR=$(shell $(_FC) -dM -E - < /dev/null 2> /dev/null | egrep __GNUC_MINOR | cut -c25)
```

Edit the `nwchem_env.sh` file:

```shell
$ vim nwchem_env.sh
export NWCHEM_TOP=/c/cygwin_wm/opt_win/NWChem_7.0.2
export NWCHEM_TARGET=LINUX
export USE_MPI=yes
export MPI_LOC=/c/cygwin_wm/opt_win/MPICH2/
export MPI_INCLUDE=$MPI_LOC/include
export MPI_LIB=$MPI_LOC/lib
export LIBMPI="-lfmpich2g -lmpi"
export NWCHEM_MODULES=all
export DEPEND_CC=gcc
export USE_INTERNALBLAS=y
```

Source the environment variables:

```shell
$ source nwchem_env.sh
```

Generate the configuration:

```shell
$ cd $NWCHEM_TOP/src
$ make nwchem_config
```

Compile NWChem:

```shell
$ make FC=gfortran FFLAGS="-fallow-argument-mismatch" DEPEND_CC=gcc
```

To enable static linking, re-run the last three commands with `-static`:

```shell
$ gfortran -fallow-argument-mismatch -I.  -I/c/cygwin_wm/opt_win/NWChem_7.0.2/src/include -I/c/cygwin_wm/opt_win/NWChem_7.0.2/src/tools/install/include -DLINUX -DGFORTRAN -DGCC46 -DCHKUNDFLW -DGCC4 -DPARALLEL_DIAG -DCOMPILATION_DATE="'`date +%a_%b_%d_%H:%M:%S_%Y`'" -DCOMPILATION_DIR="'/c/cygwin_wm/opt_win/NWChem_7.0.2'" -DNWCHEM_BRANCH="'7.0.2'"  -c -o nwchem.o nwchem.F

$ gfortran -fallow-argument-mismatch -I.  -I/c/cygwin_wm/opt_win/NWChem_7.0.2/src/include -I/c/cygwin_wm/opt_win/NWChem_7.0.2/src/tools/install/include -DLINUX -DGFORTRAN -DGCC46 -DCHKUNDFLW -DGCC4 -DPARALLEL_DIAG -DCOMPILATION_DATE="'`date +%a_%b_%d_%H:%M:%S_%Y`'" -DCOMPILATION_DIR="'/c/cygwin_wm/opt_win/NWChem_7.0.2'" -DNWCHEM_BRANCH="'7.0.2'"  -c -o stubs.o stubs.F

$ gfortran -fno-tree-dominator-opts  -std=legacy -m32 -march=pentium4 -mtune=pentium4 -g -O0  -L/c/cygwin_wm/opt_win/NWChem_7.0.2/lib/LINUX -L/c/cygwin_wm/opt_win/NWChem_7.0.2/src/tools/install/lib  -o /c/cygwin_wm/opt_win/NWChem_7.0.2/bin/LINUX/nwchem nwchem.o stubs.o -lnwctask -lccsd -lmcscf -lselci -lmp2 -lmoints -lstepper -ldriver -loptim -lnwdft -lgradients -lcphf -lesp -lddscf -ldangchang -lguess -lhessian -lvib -lnwcutil -lrimp2 -lproperty -lsolvation -lnwints -lprepar -lnwmd -lnwpw -lofpw -lpaw -lpspw -lband -lnwpwlib -lcafe -lspace -lanalyze -lqhop -lpfft -ldplot -ldrdy -lvscf -lqmmm -lqmd -letrans -lpspw -ltce -lbq -lmm -lcons -lperfm -ldntmc -lccca -ldimqm -lnwcutil -lga -larmci -lpeigs -lperfm -lcons -lbq -lnwcutil -lnwclapack  -lnwcblas    -lnwclapack  -lnwcblas  -L/c/cygwin_wm/opt_win/MPICH2//lib -lfmpich2g -lmpi    -lcomex -lfmpich2g -lmpi -lpthread  -lwsock32 -static
```

### Creating a Package

```shell
$ cd /c/cygwin_wm/opt_win
$ mv NWChem_7.0.2 NWChem_7.0.2_
$ mkdir -p NWChem_7.0.2/bin
$ mkdir -p NWChem_7.0.2/data
$ cp NWChem_7.0.2_/bin/LINUX/nwchem NWChem_7.0.2/bin
$ chmod 755 NWChem_7.0.2/bin/nwchem
$ cp /c/Windows/SysWOW64/fmpich2g.dll NWChem_7.0.2/bin
$ cp /c/Windows/SysWOW64/mpich2mpi.dll NWChem_7.0.2/bin
$ cp /c/Windows/SysWOW64/mpich2nemesis.dll NWChem_7.0.2/bin
$ cp -r NWChem_7.0.2_/src/basis/libraries NWChem_7.0.2/data
$ cp -r NWChem_7.0.2_/src/data NWChem_7.0.2
$ cp -r NWChem_7.0.2_/src/nwpw/libraryps NWChem_7.0.2/data
$ vim NWChem_7.0.2/data/default.nwchemrc
nwchem_basis_library /c/cygwin_wm/opt_win/NWChem_7.0.2/data/libraries/
nwchem_nwpw_library /c/cygwin_wm/opt_win/NWChem_7.0.2/data/libraryps/
ffield amber
amber_1 /c/cygwin_wm/opt_win/NWChem_7.0.2/data/amber_s/
amber_2 /c/cygwin_wm/opt_win/NWChem_7.0.2/data/amber_q/
amber_3 /c/cygwin_wm/opt_win/NWChem_7.0.2/data/amber_x/
amber_4 /c/cygwin_wm/opt_win/NWChem_7.0.2/data/amber_u/
spce /c/cygwin_wm/opt_win/NWChem_7.0.2/data/solvents/spce.rst
charmm_s /c/cygwin_wm/opt_win/NWChem_7.0.2/data/charmm_s/
charmm_x /c/cygwin_wm/opt_win/NWChem_7.0.2/data/charmm_x/
```

