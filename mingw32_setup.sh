#!/bin/bash -x

# pacman -S mingw-w64-i686-toolchain

DIR_MSMPI=/c/cygwin_wm/tmp/msmpi32

# Must run on 32-bit machine
function PrepareMSMPI() {
  rm -fr ${DIR_MSMPI}
  mkdir -p ${DIR_MSMPI}
  cd ${DIR_MSMPI}
  
  cp "$MSMPI_LIB32/msmpi.lib" .
  cp "$WINDIR/system32/msmpi.dll" .
  gendef msmpi.dll
  dlltool -d msmpi.def -D msmpi.dll -l libmsmpi.a
  cp "$MSMPI_INC/mpif.h" .
  cp "$MSMPI_INC/x86/mpifptr.h" .
  cp "$MSMPI_INC"/mpi.h .
  cp "$MSMPI_BIN"/mpiexec.exe .
  cp "$MSMPI_BIN"/smpd.exe .
  
  cd ..
}

function InstallMSMPI() {
  DIR_DEST=/c/cygwin_wm/mingw/MicrosoftMPI

  rm -rf ${DIR_DEST}/bin32

  mkdir -p                           ${DIR_DEST}/bin32
  
  cp ${DIR_MSMPI}/mpiexec.exe        ${DIR_DEST}/bin32      || exit 1
#  cp ${DIR_MSMPI}/msmpi.dll          ${DIR_DEST}/bin32      || exit 1
  cp ${DIR_MSMPI}/smpd.exe           ${DIR_DEST}/bin32      || exit 1
  cp ${DIR_MSMPI}/MicrosoftMPI-Redistributable-EULA.rtf ${DIR_DEST}/bin32        || exit 1
}

# https://lammps.sandia.gov/tars/
function InstallLammps() {
  DIR_DEST=/c/cygwin_wm/mingw/lammps-30Jul16

  rm -rf ${DIR_DEST}/bin32
  
  cd mingw32

  wget https://winmostar.com/wm/cygwin_wm/packages/lammps-30Jul16.tar.gz
  rm -fr lammps-30Jul16
  tar xf lammps-30Jul16.tar.gz

  patch -u -p1 -d lammps-30Jul16 < ../lammps-30Jul16_32bit_mingw.patch || exit 1

  cd lammps-30Jul16
  
  cd src
  make yes-misc
  make yes-rigid
  make yes-user-reaxc
#    NOTE: Do not use "|| exit 1" after "make" because this returns an error code
#      even though it is working properly.
  make serial
  make mpi
  cd ..
  
  mkdir -p                           ${DIR_DEST}/bin32
  mkdir -p                           ${DIR_DEST}/Doc
  mkdir -p                           ${DIR_DEST}/Potentials
  mkdir -p                           ${DIR_DEST}/Benchmarks

  cp src/lmp_serial.exe              ${DIR_DEST}/bin32/      || exit 1
  cp src/lmp_mpi.exe                 ${DIR_DEST}/bin32/      || exit 1
  cp LICENSE                         ${DIR_DEST}/            || exit 1
  cp potentials/*                    ${DIR_DEST}/Potentials/ || exit 1
  cp bench/in.*                      ${DIR_DEST}/Benchmarks/ || exit 1
  cp doc/Manual.pdf                  ${DIR_DEST}/Doc/        || exit 1

  cp /mingw32/bin/libstdc++-6.dll    ${DIR_DEST}/bin32/      || exit 1
  cp /mingw32/bin/libgcc_s_dw2-1.dll ${DIR_DEST}/bin32/      || exit 1

  cp ${DIR_MSMPI}/msmpi.dll          ${DIR_DEST}/bin32/      || exit 1
  
  cd ..
  
  cd ..
}

function InstallQE() {
  DIR_DEST=/c/cygwin_wm/mingw/q-e-qe-6.4.1

  rm -rf ${DIR_DEST}/bin32
  
  cd mingw32

  wget https://gitlab.com/QEF/q-e/-/archive/qe-6.4.1/q-e-qe-6.4.1.tar.gz
  rm -fr q-e-qe-6.4.1
  tar zxvf q-e-qe-6.4.1.tar.gz

  cd q-e-qe-6.4.1
  ./configure ARCH=mingw32 CC=gcc DFLAGS="-D_WIN32 -D__DFTI -D__MPI" || exit 1
  pwd
  ls ../../
  patch -p1 < ../../q-e-qe-6.4.1_32bit_mingw.patch || exit 1

  make all || exit 1

  cd bin
  for f in *.x; do mv $f ${f%x}exe; done
  cd ..
  
  mkdir -p                            ${DIR_DEST}/bin32
  mkdir -p                            ${DIR_DEST}/doc
  mkdir -p                            ${DIR_DEST}/pseudo
  
  cp bin/*                            ${DIR_DEST}/bin32/      || exit 1
  cp License                          ${DIR_DEST}/            || exit 1
  cp Doc/*.pdf                        ${DIR_DEST}/doc/        || exit 1
  cp PW/Doc/*.pdf                     ${DIR_DEST}/doc/        || exit 1
  cp PW/Doc/*.html                    ${DIR_DEST}/doc/        || exit 1

  cp /c/Windows/System32/ws2_32.dll   ${DIR_DEST}/bin32/      || exit 1
  cp /c/Windows/System32/msmpi.dll    ${DIR_DEST}/bin32/      || exit 1
  cp /mingw32/bin/libgcc_s_dw2-1.dll  ${DIR_DEST}/bin32/      || exit 1
  cp /mingw32/bin/libgfortran-5.dll   ${DIR_DEST}/bin32/      || exit 1
  cp /mingw32/bin/libquadmath-0.dll   ${DIR_DEST}/bin32/      || exit 1
  cp /mingw32/bin/libwinpthread-1.dll ${DIR_DEST}/bin32/      || exit 1
  cp /mingw32/bin/libgomp-1.dll       ${DIR_DEST}/bin32/      || exit 1
  cp /c/Program\ Files\ \(x86\)/IntelSWTools/compilers_and_libraries/windows/redist/ia32_win/mkl/*.dll               ${DIR_DEST}/bin32/ || exit 1
  cp /c/Program\ Files\ \(x86\)/IntelSWTools/compilers_and_libraries/windows/redist/ia32_win/compiler/libiomp5md.dll ${DIR_DEST}/bin32/ || exit 1

  cp ${DIR_MSMPI}/msmpi.dll          ${DIR_DEST}/bin32/       || exit 1

  cd ..
  
  cd ..
}

set -x

date

export LANG=C

cd /c/cygwin_wm/tmp

mkdir -p mingw32

#PrepareMSMPI

#InstallMSMPI
InstallLammps
InstallQE

date
