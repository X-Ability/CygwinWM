#!/bin/bash -x

function InstallMPICH2() {
  rm -rf /usr/local/mpich2-1.5
# http://www.mpich.org/static/downloads/1.5/mpich2-1.5.tar.gz
  wget https://winmostar.com/wm/cygwin_wm/packages/mpich2-1.5.tar.gz
  rm -fr mpich2-1.5
  tar xfz mpich2-1.5.tar.gz
  cd mpich2-1.5
  ./configure --prefix=/usr/local/mpich2-1.5 || exit 1
  make -j 4 || exit 1
  make install || exit 1
  cd ..
}

function InstallGrace() {
  rm -rf /usr/local/grace
  wget https://winmostar.com/wm/cygwin_wm/packages/grace-5.1.25.tar.gz
  rm -fr grace-5.1.25
  tar xfz grace-5.1.25.tar.gz
  cd grace-5.1.25
  ./configure || exit 1
  make || exit 1
  make install || exit 1
  cd ..
}

function InstallOpenBabel() {
  rm -fr openbabel-2.4.1
  rm -f openbabel-2-4-1.obgmxtop.patch
  wget https://winmostar.com/wm/cygwin_wm/packages/openbabel-2.4.1.tar.gz
  rm -fr openbabel-2.4.1
  tar xfz openbabel-2.4.1.tar.gz
  patch -u -p1 -d openbabel-2.4.1 < openbabel-2-4-1.obgmxtop.patch
  cd openbabel-2.4.1
  mkdir build
  cd build
  cmake .. || exit 1
  make || exit 1
  make install || exit 1
  cd ../..
}

function InstallGromacs() {
  rm -rf /usr/local/gromacs_sse
  rm -rf /usr/local/gromacs_avx
  wget https://winmostar.com/wm/cygwin_wm/packages/gromacs-5.0.7.tar.gz
  rm -fr gromacs-5.0.7
  tar xfz gromacs-5.0.7.tar.gz
  cd gromacs-5.0.7

  OPT_COMMON="-DGMX_BUILD_OWN_FFTW=ON -DGMX_GPU=OFF -DGMX_USE_RDTSCP=OFF "
  OPT_SSE="-DGMX_SIMD=SSE2    -DCMAKE_INSTALL_PREFIX=/usr/local/gromacs_sse "
  OPT_AVX="-DGMX_SIMD=AVX_256 -DCMAKE_INSTALL_PREFIX=/usr/local/gromacs_avx "
  OPT_SINGLE="-DGMX_DOUBLE=OFF "
  OPT_DOUBLE="-DGMX_DOUBLE=ON  "

# To undefine HAVE_SCHED_AFFINITY
  sed -i.bak -e 's/^.*Affinity.*$//g' -e 's/^.*AFFINITY.*$//g' CMakeLists.txt 

  mkdir build_sse
  cd build_sse
  cmake .. $OPT_COMMON $OPT_SSE $OPT_SINGLE || exit 1
  make -j 4 || exit 1
  make install || exit 1
  cd ..

  mkdir build_sse_d
  cd build_sse_d
  cmake .. $OPT_COMMON $OPT_SSE $OPT_DOUBLE || exit 1
  make -j 4 || exit 1
  make install || exit 1
  cd ..

  mkdir build_avx
  cd build_avx
  cmake .. $OPT_COMMON $OPT_AVX $OPT_SINGLE || exit 1
  make -j 4 || exit 1
  make install || exit 1
  cd ..

  mkdir build_avx_d
  cd build_avx_d
  cmake .. $OPT_COMMON $OPT_AVX $OPT_DOUBLE || exit 1
  make -j 4 || exit 1
  make install || exit 1
  cd ..

  cd ..
}

function InstallAmberTools18() {
  rm -rf /usr/local/amber18/
  dir=`pwd`
  wget https://winmostar.com/wm/cygwin_wm/packages/AmberTools18.tar.bz2
  tar xfj AmberTools18.tar.bz2
  mv amber18 /usr/local/
  export AMBERHOME=/usr/local/amber18
  cd $AMBERHOME
  echo y | ./configure -cygwin -noX11 -nosse --skip-python gnu
  sed -i 's/LDFLAGS=/LDFLAGS=-Wl,--allow-multiple-definition/' AmberTools/src/cpptraj/config.h
  sed -i 's#${AMBER_PREFIX}/bin:${PATH}#${PATH}:${AMBER_PREFIX}/bin#g' amber.sh
  source amber.sh || exit 1
  cp /tmp/parmchk2.c AmberTools/src/antechamber/
  make install AMBERBUILDFLAGS="-Wl,--allow-multiple-definition -DWITHOUT_BACKTRACE" || exit 1
  cd $dir
}

function InstallAcpype18() {
  rm -rf /usr/local/acpype/ ./acpype_r10101
  wget https://winmostar.com/wm/cygwin_wm/packages/acpype_r10101.tgz
  tar xvfz acpype_r10101.tgz
  chmod a+x acpype.py
  mv acpype.py acpype_r10101/
  mv acpype_r10101 /usr/local/acpype
}

function InstallERmod() {
  rm -rf /usr/local/ermod
  wget https://winmostar.com/wm/cygwin_wm/packages/ermod-0.3.4.tar.gz
  rm -fr ermod-0.3.4
  tar xfz ermod-0.3.4.tar.gz
  cd ermod-0.3.4/vmdplugins
  make
  cp compile/*.so libexec/
  cd ..
  ./configure --prefix=/usr/local/ermod --disable-mpi --enable-openmp --disable-opt || exit 1
  make || exit 1
  make install || exit 1
  cd ..
}

function InstallMODYLAS() {
  rm -rf /usr/local/MODYLAS_1.0.4
  rm -fr MODYLAS_1.0.4
  tar xfz MODYLAS_1.0.4.tar_1.gz
  patch -u -p0 < MODYLAS_1.0.4.patch
  cd MODYLAS_1.0.4/source
  export FCFLAGS="-DMPIPARA -cpp -O3 -DCOMM_CUBE -DFJMPIDIR -DSYNC_COM -DONEPROC_AXIS"
  ./configure --with-kind-fortran-compiler=INTEL --prefix=/usr/local/MODYLAS_1.0.4 || exit 1
  cd src
  make || exit 1
  make install || exit 1
  cd ../../..
}

function InstallOpenMX() {
  rm -rf /usr/local/openmx3.8
  rm -rf /usr/local/fermisurfer
  wget https://winmostar.com/wm/cygwin_wm/packages/openmx3.8.tar.gz
  wget https://winmostar.com/wm/cygwin_wm/packages/patch3.8.4.tar.gz
  wget https://winmostar.com/wm/cygwin_wm/packages/openmx3.8.4.patch
  wget https://winmostar.com/wm/cygwin_wm/packages/fermisurfer_1.7.1.zip
  wget https://winmostar.com/wm/cygwin_wm/packages/openmx3.9.tar.gz
  rm -fr openmx3.8
  tar xvfz openmx3.8.tar.gz
  cd openmx3.8/source
  tar xvfz ../../patch3.8.4.tar.gz
  touch TRAN_Main_Analysis.c
  touch TRAN_Main_Analysis_NC.c
  cd ../..
  patch -u -p0 < openmx3.8.4.patch
  cd openmx3.8/source
  sed -i.bak -e '/^CC=/s/gcc/mpicc/g' makefile
  sed -i.bak -e '/^FC=/s/gfortran/mpif90/g' makefile
  sed -i.bak -e 's|/usr/lib/libmpi_mpifh.dll.a /usr/lib/libmpi.dll.a||g' makefile
  make || exit 1
  make install || exit 1
  gcc bandgnu13.c -lm -o ../work/bandgnu13.exe || exit 1
  make DosMain || exit 1
  cd ../..
  mv openmx3.8 /usr/local
  unzip fermisurfer_1.7.1.zip
  mv fermisurfer /usr/local/
  rm -fr openmx3.9
  tar xvfz openmx3.9.tar.gz
  echo '#!/bin/bash' > MX_TRAP.sh
  sed s/TEST/TEST_/g openmx3.9/work/MX_TRAP.sh >> MX_TRAP.sh
  chmod 755 MX_TRAP.sh
  mv MX_TRAP.sh /usr/local/openmx3.8/work/MX_TRAP.sh
}

function InstallPhonopy() {
  pip3.7 install pkgconfig || exit 1
  pip3.7 install h5py || exit 1
  
  wget https://winmostar.com/wm/cygwin_wm/packages/phonopy-1.12.6.53.tar.gz
  rm -fr phonopy-1.12.6.53
  tar xfz phonopy-1.12.6.53.tar.gz
  cd phonopy-1.12.6.53
  python3.7 setup.py install || exit 1
  cd ..
}

function InstallMDTraj() {
  wget https://winmostar.com/wm/cygwin_wm/packages/mdtraj-1.9.0.tar.gz
  rm -fr mdtraj-1.9.0
  tar xfz mdtraj-1.9.0.tar.gz
  cd mdtraj-1.9.0
  patch -u -p1 -d mdtraj < ../mdtraj-1.9.0_rev2.patch
  python3.7 setup.py install || exit 1
  cd ..
}

function InstallBoltzTraP() {
  rm -rf /usr/local/boltztrap-1.2.5
  wget https://winmostar.com/wm/cygwin_wm/packages/BoltzTraP.tar.bz2
  rm -fr boltztrap-1.2.5
  tar xfj BoltzTraP.tar.bz2
  cd boltztrap-1.2.5/src/
  rm BoltzTraP
  sed -i s/:log/log/g x_trans
  make FOPT="-g -funroll-loops -O3 -ffast-math -fgcse-lm -fgcse-sm -ffast-math -ftree-vectorize -fexternal-blas" || exit 1
  cd ../..
  mv boltztrap-1.2.5/ /usr/local/
}

function InstallMATCH() {
  rm -rf /usr/local/MATCH_RELEASE
  wget https://winmostar.com/wm/cygwin_wm/packages/MATCH_RELEASE.tar.gz
  tar xzf MATCH_RELEASE.tar.gz
  mv MATCH_RELEASE /usr/local/
}

function InstallConditionalERmod() {
  rm -rf /usr/local/ermod_conditional
  wget https://winmostar.com/wm/cygwin_wm/packages/ermod-0.3.5.tar.gz
  rm -fr ermod-0.3.5
  tar xfz ermod-0.3.5.tar.gz
  cp -r ermod_conditional/* ermod-0.3.5/
  cd ermod-0.3.5/vmdplugins
  make || exit 1
  cp compile/*.so libexec/
  cd ..
  ./configure --prefix=/usr/local/ermod_conditional --disable-mpi --enable-openmp --disable-opt || exit 1 || exit 1
  make || exit 1
  make install || exit 1
  cd ..
}

function InstallEnumlib() {
  rm -rf /usr/local/enumlib
  wget https://winmostar.com/wm/cygwin_wm/packages/enumlib-1.0.8.zip
  wget https://winmostar.com/wm/cygwin_wm/packages/symlib-1.1.0.zip
  rm -fr enumlib-1.0.8
  unzip enumlib-1.0.8.zip
  unzip symlib-1.1.0.zip
  mv symlib-1.1.0/* enumlib-1.0.8/symlib/
  cd enumlib-1.0.8/symlib/src
  export F90=gfortran
  make || exit 1
  cd ../../src
  make enum.x || exit 1
  make polya.x || exit 1
  make makestr.x || exit 1
  cd ../..
  mv enumlib-1.0.8 /usr/local/enumlib
  export F90=
}

function InstallPackmol() {
  rm -rf /usr/local/packmol-18.166
  wget https://winmostar.com/wm/cygwin_wm/packages/packmol-18.166.zip
  rm -fr packmol-18.166
  unzip packmol-18.166.zip
  cd packmol-18.166
  ./configure `which gfortran` || exit 1
  sed -i -e s/--fast-math/--fast-math\ -m32/g Makefile
  make || exit 1
  cd ..
  mv packmol-18.166 /usr/local/
}

function InstallMktop() {
  rm -rf /usr/local/mktop_2.2.1
  wget https://winmostar.com/wm/cygwin_wm/packages/mktop_2.2.1.tar
  rm -rf mktop_2.2.1
  mkdir mktop_2.2.1
  cd mktop_2.2.1
  tar xvf ../mktop_2.2.1.tar
  cp ../mktop.pl .
  chmod 755 mktop.pl
  cd ..
  mv mktop_2.2.1 /usr/local/
}

# If setup.py returns an error, launch C:\cygwin_wm\bin\ash.exe and execute "cd /bin" and "./rebaseall".
function InstallPymatgen() {
  wget https://winmostar.com/wm/cygwin_wm/packages/cython-0.29.13.zip
  rm -fr cython-0.29.13
  unzip cython-0.29.13.zip
  cd cython-0.29.13
  python3.7 setup.py install || exit 1
  cd ..
  
  wget https://winmostar.com/wm/cygwin_wm/packages/scipy-1.1.0.zip
  rm -fr scipy-1.1.0
  unzip scipy-1.1.0.zip
  cd scipy-1.1.0
  python3.7 setup.py install || exit 1
  cd ..
  
  wget https://winmostar.com/wm/cygwin_wm/packages/spglib-1.15.1.zip
  rm -fr spglib-1.15.1
  unzip spglib-1.15.1.zip
  cd spglib-1.15.1/python
  python3.7 setup.py install || exit 1
  cd ../..
  
  pip3.7 install wheel || exit 1

  wget https://winmostar.com/wm/cygwin_wm/packages/matplotlib-3.1.0.zip
  rm -fr matplotlib-3.1.0
  unzip matplotlib-3.1.0.zip
  cd matplotlib-3.1.0
  python3.7 setup.py install || exit 1
  cd ..
  
  wget https://winmostar.com/wm/cygwin_wm/packages/pandas-1.0.3.zip
  rm -fr pandas-1.0.3
  unzip pandas-1.0.3.zip
  cd pandas-1.0.3
  python3.7 setup.py install || exit 1
  cd ..
  
  pip3.7 install ase || exit 1

  wget https://winmostar.com/wm/cygwin_wm/packages/pymatgen-2020.4.29.zip
  rm -fr pymatgen-2020.4.29
  unzip pymatgen-2020.4.29.zip
  cd pymatgen-2020.4.29
  python3.7 setup.py install || exit 1
  cd ..
}

function InstallLammps() {
  rm -rf /usr/local/lammps-30Jul16
  wget https://winmostar.com/wm/cygwin_wm/packages/lammps-30Jul16.tar.gz
  rm -fr lammps-30Jul16
  tar xf lammps-30Jul16.tar.gz
  cd lammps-30Jul16/src
  make yes-misc
  make yes-rigid
  make yes-user-reaxc
  sed -i.bak -e '/^LMP_INC/s/$/ -DLAMMPS_XDR/g' MAKE/Makefile.mpi
  sed -i.bak -e '/^LMP_INC/s/$/ -DLAMMPS_XDR/g' MAKE/Makefile.serial
  make serial || exit 1
  make mpi || exit 1
  cd ..
  cd ..
  mv lammps-30Jul16 /usr/local/
}

set -x

date

SCRIPT_DIR=$(cd $(dirname $0); pwd)

export LANG=C

grep "export LANG=C" /etc/profile.d/lang.sh || echo "export LANG=C" >> /etc/profile.d/lang.sh

cat << EOF > /etc/profile.d/winmostar.sh
export PATH=\$PATH:/usr/local/mpich2-1.5/bin
export PATH=\$PATH:/usr/local/grace/bin
source /usr/local/gromacs_sse/bin/GMXRC
#source /usr/local/gromacs_avx/bin/GMXRC
source /usr/local/amber18/amber.sh
export PATH=\$PATH:/usr/local/acpype
export PATH=\$PATH:/usr/local/MODYLAS_1.0.4/bin
export PATH=\$PATH:/usr/local/openmx3.8/work:/usr/local/fermisurfer/bin
export OPENMX_DATA_PATH=/usr/local/openmx3.8/DFT_DATA13
export PATH=\$PATH:/usr/local/boltztrap-1.2.5/src
export PATH=\$PATH:/usr/local/boltztrap-1.2.5/util
export PerlChemistry=/usr/local/MATCH_RELEASE/PerlChemistry
export MATCH=/usr/local/MATCH_RELEASE/MATCH
export PATH=\$PATH:/usr/local/MATCH_RELEASE/MATCH/scripts
export PATH=\$PATH:/usr/local/packmol-18.166
export PATH=\$PATH:/usr/local/mktop_2.2.1
# To prevent calling native binary of LAMMPS, 
#   substitute /usr/local/lammps-30Jul16/src before $PATH.
export PATH=/usr/local/lammps-30Jul16/src:\$PATH
EOF

. /etc/profile.d/winmostar.sh

cd $SCRIPT_DIR
cp -f ChangeLog setup-x86.exe /cygdrive/c/cygwin_wm/
cp -rf src/* /tmp/

cd /tmp/

if [ ! -f /bin/python3.exe ]; then
  ln -s /bin/python3.7.exe /bin/python3.exe
fi

if [ ! -f /bin/python2.exe ]; then
  ln -s /bin/python2.7.exe /bin/python2.exe
fi

# Reqiured to run acpype
if [ ! -f /bin/python.exe ]; then
  ln -s /bin/python2.exe /bin/python.exe
fi

# required to build phonopy
# https://seesaawiki.jp/w/kou1okada/d/20200205%3A%20Cygwin%20-%20Python3%20-%20h5py
if [ ! -f /bin/hdf5.dll ]; then
  ln -s /bin/cyghdf5-103.dll /bin/hdf5.dll
fi

# os.symlink of python2.7 makes symbolic links to /mnt/ for windows paths
ln -s /cygdrive /mnt


InstallMPICH2
InstallGrace
InstallOpenBabel
InstallGromacs
InstallAmberTools18
InstallAcpype18
InstallERmod
InstallMODYLAS
InstallOpenMX
InstallBoltzTraP
InstallMATCH
InstallConditionalERmod
InstallEnumlib
InstallPackmol
InstallMktop
InstallPymatgen
InstallPhonopy
InstallMDTraj
#InstallLammps

date
