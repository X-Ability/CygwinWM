#!/bin/bash -x

DIR_MSMPI=/c/cygwin_wm/tmp/msmpi64

function PrepareMSMPI() {
  rm -fr ${DIR_MSMPI}
  mkdir -p ${DIR_MSMPI}
  cd ${DIR_MSMPI}
  
  cp "$MSMPI_LIB64/msmpi.lib" .
  cp "$WINDIR/system32/msmpi.dll" .
  gendef msmpi.dll
  dlltool -d msmpi.def -D msmpi.dll -l libmsmpi.a
  cp "$MSMPI_INC/mpif.h" .
  cp "$MSMPI_INC/x64/mpifptr.h" .
  cp "$MSMPI_INC"/mpi.h .
  cp "$MSMPI_BIN"/mpiexec.exe .
  cp "$MSMPI_BIN"/smpd.exe .
  
  cd ..
}

# https://lammps.sandia.gov/tars/
function InstallLammps() {
  DIR_DEST=/c/cygwin_wm/usr/local/lammps-30Jul16_64bit

  rm -rf ${DIR_DEST}
  
  mkdir -p mingw64
  cd mingw64

  wget https://winmostar.com/wm/cygwin_wm/packages/lammps-30Jul16.tar.gz
  rm -fr lammps-30Jul16
  tar xf lammps-30Jul16.tar.gz

  patch -u -p1 -d lammps-30Jul16 < ../lammps-30Jul16_64bit_mingw.patch

  cd lammps-30Jul16/src
  make yes-misc
  make yes-rigid
  make yes-user-reaxc
#    NOTE: Do not use "|| exit 1" after "make" because this returns an error code
#      even though it is working properly.
  make mpi
  cd ..
  
  mkdir -p                           ${DIR_DEST}/bin
  mkdir -p                           ${DIR_DEST}/Doc
  mkdir -p                           ${DIR_DEST}/Potentials
  mkdir -p                           ${DIR_DEST}/Benchmarks
  cp src/lmp_mpi.exe                 ${DIR_DEST}/bin/        || exit 1
  cp LICENSE                         ${DIR_DEST}/            || exit 1
  cp potentials/*                    ${DIR_DEST}/Potentials/ || exit 1
  cp bench/in.*                      ${DIR_DEST}/Benchmarks/ || exit 1
  cp doc/Manual.pdf                  ${DIR_DEST}/Doc/        || exit 1
  cp /mingw64/bin/libstdc++-6.dll    ${DIR_DEST}/bin/        || exit 1
  cp ${DIR_MSMPI}/mpiexec.exe        ${DIR_DEST}/bin/        || exit 1
  cp ${DIR_MSMPI}/msmpi.dll          ${DIR_DEST}/bin/        || exit 1
  cp ${DIR_MSMPI}/smpd.exe           ${DIR_DEST}/bin/        || exit 1
  cp ${DIR_MSMPI}/MicrosoftMPI-Redistributable-EULA.rtf ${DIR_DEST}/bin/        || exit 1
  
  cd ..
}

set -x

date

export LANG=C

cd /c/cygwin_wm/tmp

#PrepareMSMPI
InstallLammps

date
