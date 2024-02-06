### Installing MSYS2

See [README_MSYS2](README_MSYS2.md)


### Obtaining LAMMPS

Download LAMMPS from the following link and extracting it. You might encounter errors during extraction. Since it's under the example directory, there might not be a problem:

```shell
$ wget https://download.lammps.org/tars/lammps-29Sep2021.tar.gz
$ tar xvfz lammps-29Sep2021.tar.gz
$ cd lammps-29Sep2021/
```

Modify the following files in `src/MAKE/` as shown below:

In `src/MAKE/Makefile.serial` and `src/MAKE/Makefile.mpi`:

```shell
LINKFLAGS =	-g -O -std=c++11 -static-libgcc -static
LIB = -lwsock32 -lpsapi
LMP_INC =	-DLAMMPS_GZIP
```

Then, proceed with the following commands:

```shell
$ cd src
$ make yes-molecule yes-misc yes-rigid yes-reaxff yes-extra-dump \
       yes-mc yes-kspace yes-phonon yes-qeq yes-reaction yes-replica \
       yes-dpd-basic yes-extra-molecule yes-manybody yes-meam \
       yes-extra-fix yes-fep yes-tally yes-colvars
$ make lib-colvars args="-m serial"
$ make serial
$ make lib-colvars args="-m mpi"
$ make mpi
$ cd ..
```

### Creating a Package

Create a package using the following steps:

```shell
$ mkdir -p dist/bin
$ cp src/lmp_serial.exe dist/bin
$ cp src/lmp_mpi.exe dist/bin
$ strip dist/bin/*.exe
$ cp -r potentials dist/
```

Copy `msmpi.dll` to the `dist/bin` directory.

