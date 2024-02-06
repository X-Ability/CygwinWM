### Installing MSYS2

See [README_MSYS2](/README_MSYS2.md)

Additionally, install `mingw-w64-x86_64-lapack` by running:

```shell
pacman -S mingw-w64-x86_64-lapack
```

### Obtaining Atomsk source

1. Download the Atomsk source code from the following link: [Atomsk Source](https://atomsk.univ-lille.fr/code/atomsk_b0.13.1.tar.gz).

2. Extract the downloaded archive:

```shell
$ wget https://atomsk.univ-lille.fr/code/atomsk_b0.13.1.tar.gz
$ tar xvfz atomsk_b0.13.1.tar.gz
$ cd atomsk_b0.13.1/src
```

### Modifying Makefile.static

Edit the `Makefile.static` file as follows:

Change the `FFLAGS` and `LAPACK` lines as shown below:

```shell
FFLAGS=-O2 $(OPENMP) -fno-backslash -I..$(SEP)$(OBJ) -J..$(SEP)$(OBJ) -static -fPIC -static-libgfortran
LAPACK=-llapack -lblas
```

### Building Atomsk

Compile Atomsk with the modified `Makefile.static`:

```shell
$ make -f Makefile.static atomsk
```

### Distributing Atomsk

Move the `atomsk.exe` binary to the `/opt_win/atomsk_b0.13.1/bin/` directory in the CygwinWM environment for distribution.
```


