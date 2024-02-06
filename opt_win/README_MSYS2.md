### Installing MSYS2

Install MSYS2 following the installation instructions on [MSYS2's website](https://www.msys2.org/).

1. Launch MSYS2 MinGW x64.
2. Open a command prompt and update the package database:

   ```shell
   $ pacman -Sy
   $ pacman -Su
   ```

3. Install the necessary packages:

   ```shell
   $ pacman -S patch
   $ pacman -S m4
   $ pacman -S python
   $ pacman -S mingw-w64-x86_64-gcc
   $ pacman -S mingw-w64-x86_64-gcc-fortran
   $ pacman -S make
   $ pacman -S git
   ```

4. Verify the installed GCC version:

   ```shell
   $ gcc --version
   gcc.exe (Rev5, Built by MSYS2 project) 12.2.0
   ```

5. Verify the installed Make version:

   ```shell
   $ make --version
   GNU Make 4.3
   Built for x86_64-pc-msys
   ```

### Installing libxc

Download the latest `libxc-5.2.3.tar.gz` from the [libxc download page](https://tddft.org/programs/libxc/download/).

1. Extract the downloaded file:

   ```shell
   $ tar xvfz libxc-5.2.3.tar.gz
   $ cd libxc-5.2.3
   ```

2. Run the configuration:

   ```shell
   $ ./configure
   ```

3. Build libxc:

   ```shell
   $ make
   ```

4. Install libxc:

   ```shell
   $ make install
   ```

### Installing MS-MPI

Install MS-MPI for MSYS2 using the following command:

```shell
$ pacman -S mingw-w64-x86_64-msmpi
```

However, for execution, you will need the Windows version of MS-MPI. Download it from [Microsoft's download page](https://www.microsoft.com/en-us/download/details.aspx?id=100593) and copy `mpiexec.exe` to `mpirun.exe` in C:\\Program\ Files\\Microsoft\ MPI\\Bin.

