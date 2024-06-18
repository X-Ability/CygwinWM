#!/bin/bash -x

# MPICH2 used for MODYLAS, OpenMX and ERmod
# http://www.mpich.org/static/downloads/1.5/mpich2-1.5.tar.gz
function InstallMPICH2() {
  rm -rf /usr/local/mpich2-1.5
  wget https://winmostar.com/wm/cygwin_wm/packages/mpich2-1.5.tar.gz
  rm -fr mpich2-1.5
  tar xfz mpich2-1.5.tar.gz
  patch -u -p1 -d mpich2-1.5 < mpich2-1.5-fno-common.patch
  cd mpich2-1.5
  ./configure --prefix=/usr/local/mpich2-1.5 FFLAGS=-fallow-argument-mismatch || exit 1
  make -j 4 || exit 1
  make install || exit 1
  cd ..
}

# https://github.com/openbabel/openbabel/archive/openbabel-2-4-1.tar.gz
function InstallOpenBabel() {
  rm -fr openbabel-2.4.1
  #rm -f openbabel-2-4-1.obgmxtop.patch
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

# Since 5.1 and later are not compatible with gmx rdf, use 5.0.7.
# http://www.gromacs.org/Downloads
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


# http://ambermd.org/GetAmber.php#ambertools
function InstallAmberTools18() {
  rm -rf /usr/local/amber18/
  dir=`pwd`
  wget https://winmostar.com/wm/cygwin_wm/packages/AmberTools18.tar.bz2
  tar xfj AmberTools18.tar.bz2
  mv amber18 /usr/local/
  export AMBERHOME=/usr/local/amber18
  cd $AMBERHOME
  echo y | ./configure -cygwin -noX11 -nosse --skip-python gnu
  patch -u -p1 < $dir/AmberTools18.patch
  source amber.sh || exit 1
  make install AMBERBUILDFLAGS="-Wl,--allow-multiple-definition --allow-argument-mismatch  -DWITHOUT_BACKTRACE" || exit 1
  cd $dir
}

# http://svn.code.sf.net/p/ccpn/code/branches/stable/ccpn/python/acpype/
function InstallAcpype18() {
  rm -rf /usr/local/acpype/ ./acpype_r10101
  wget https://winmostar.com/wm/cygwin_wm/packages/acpype_r10101.tgz
  tar xvfz acpype_r10101.tgz
  chmod a+x acpype.py
  mv acpype.py acpype_r10101/
  mv acpype_r10101 /usr/local/acpype
}

# http://sourceforge.net/projects/ermod/files/?source=navbar
function InstallERmod() {
  rm -rf /usr/local/ermod
  wget https://winmostar.com/wm/cygwin_wm/packages/ermod-0.3.4.tar.gz
  rm -fr ermod-0.3.4
  tar xfz ermod-0.3.4.tar.gz
  cd ermod-0.3.4/vmdplugins
  make
  cp compile/*.so libexec/
  cd ..
  ./configure --prefix=/usr/local/ermod --disable-mpi --enable-openmp --disable-opt FCFLAGS=--allow-argument-mismatch || exit 1
  make || exit 1
  make install || exit 1
  cd ..
}

# You must agree to the MODYLAS license in order to use MODYLAS.
# https://www.modylas.org/download
function InstallMODYLAS() {
  rm -rf /usr/local/MODYLAS_1.0.4
  rm -fr MODYLAS_1.0.4
  tar xfz MODYLAS_1.0.4.tar_1.gz
  patch -u -p0 < MODYLAS_1.0.4.patch
  cd MODYLAS_1.0.4/source
  export FCFLAGS="-DMPIPARA -cpp -O3 -DCOMM_CUBE -DFJMPIDIR -DSYNC_COM -DONEPROC_AXIS -fallow-argument-mismatch"
  ./configure --with-kind-fortran-compiler=INTEL --prefix=/usr/local/MODYLAS_1.0.4 || exit 1
  cd src
  make || exit 1
  make install || exit 1
  cd ../../..
}

# Use OpenMX.3.8 and MX_TRANS.sh from OpenMX3.9
# http://www.openmx-square.org/download.html
# https://ja.osdn.net/projects/fermisurfer/releases/
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
  cd ../..
  patch -u -p0 < openmx3.8.4.patch
  cd openmx3.8/source
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

# OpenMX3.9
# http://www.openmx-square.org/download.html
# https://ja.osdn.net/projects/fermisurfer/releases/
# https://www.netlib.org/scalapack/
function InstallOpenMX3.9() {
  rm -rf /usr/local/openmx3.9
  rm -rf /usr/local/fermisurfer
  wget https://winmostar.com/wm/cygwin_wm/packages/openmx3.9.tar.gz
  wget https://winmostar.com/wm/cygwin_wm/packages/patch3.9.9.tar
  wget https://winmostar.com/wm/cygwin_wm/packages/openmx3.9.9.patch
  wget https://winmostar.com/wm/cygwin_wm/packages/fermisurfer_1.7.1.zip
  wget https://winmostar.com/wm/cygwin_wm/packages/scalapack-2.1.0.tgz
  wget https://winmostar.com/wm/cygwin_wm/packages/scalapack-2.1.0.patch

  rm -fr scalapack-2.1.0
  tar xvfz scalapack-2.1.0.tgz
  patch -u -p0 < scalapack-2.1.0.patch
  cd scalapack-2.1.0
  make
  cd ..

  rm -fr openmx3.9
  tar xvfz openmx3.9.tar.gz
  cd openmx3.9/source
  tar xvf ../../patch3.9.9.tar
  cd ../..
  patch -u -p0 < openmx3.9.9.patch
  cd openmx3.9/source
  cp ../../scalapack-2.1.0/libscalapack.a .
  make || exit 1
  make install || exit 1
  gcc bandgnu13.c -lm -o ../work/bandgnu13.exe || exit 1
  make DosMain || exit 1
  cd ../..
  mv openmx3.9 /usr/local
  unzip fermisurfer_1.7.1.zip
  mv fermisurfer /usr/local/
}

# https://pypi.python.org/pypi/phonopy/1.12.6.53
function InstallPhonopy() {
  pip3.9 install pkgconfig || exit 1
  pip3.9 install h5py || exit 1
  pip3.9 install numpy==1.24.3 || exit 1
  pip3.9 install matplotlib==3.6.3 || exit 1
  pip3.9 install phonopy || exit 1
}

# https://github.com/mdtraj/mdtraj/releases/tag/1.9.0
function InstallMDTraj() {
  pip3.9 install scipy==1.6.3 || exit 1
  pip3.9 install mdtraj==1.9.6 || exit 1
}

# https://www.imc.tuwien.ac.at/forschungsbereich_theoretische_chemie/forschungsgruppen/prof_dr_gkh_madsen_theoretical_materials_chemistry/boltztrap/
function InstallBoltzTraP() {
  rm -rf /usr/local/boltztrap-1.2.5
  wget https://winmostar.com/wm/cygwin_wm/packages/BoltzTraP.tar.bz2
  rm -fr boltztrap-1.2.5
  tar xfj BoltzTraP.tar.bz2
  patch -u -p0 < boltztrap-1.2.5.patch
  cd boltztrap-1.2.5/src/
  rm BoltzTraP
  # To eliminate the dependence on the machine type
  make FOPT="-g -funroll-loops -O3 -ffast-math -fgcse-lm -fgcse-sm -ffast-math -ftree-vectorize -fexternal-blas -fallow-argument-mismatch" || exit 1
  cd ../..
  mv boltztrap-1.2.5/ /usr/local/
}

# https://brooks.chem.lsa.umich.edu/index.php?page=match&subdir=articles/resources/software
function InstallMATCH() {
  rm -rf /usr/local/MATCH_RELEASE
  wget https://winmostar.com/wm/cygwin_wm/packages/MATCH_RELEASE.tar.gz
  tar xzf MATCH_RELEASE.tar.gz
  mv MATCH_RELEASE /usr/local/
}

# http://sourceforge.net/projects/ermod/files/?source=navbar
function InstallConditionalERmod() {
  rm -rf /usr/local/ermod_conditional
  wget https://winmostar.com/wm/cygwin_wm/packages/ermod-0.3.5.tar.gz
  rm -fr ermod-0.3.5
  tar xfz ermod-0.3.5.tar.gz
  cp -r ermod-0.3.5_conditional/* ermod-0.3.5/
  cd ermod-0.3.5/vmdplugins
  make || exit 1
  cp compile/*.so libexec/
  cd ..
  ./configure --prefix=/usr/local/ermod_conditional --disable-mpi --enable-openmp --disable-opt FCFLAGS=--allow-argument-mismatch || exit 1 || exit 1
  make || exit 1
  make install || exit 1
  cd ..
}

# Used for structure generation from occupancy in Crystal Builder
# https://github.com/msg-byu/enumlib/tree/v1.0.8
# https://github.com/msg-byu/symlib/tree/v1.1.0
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

# https://github.com/mcubeg/packmol/releases
function InstallPackmol() {
  rm -rf /usr/local/packmol-18.166
  wget https://winmostar.com/wm/cygwin_wm/packages/packmol-18.166.zip
  rm -fr packmol-18.166
  unzip packmol-18.166.zip
  cd packmol-18.166
  ./configure `which gfortran` || exit 1
  make || exit 1
  cd ..
  mv packmol-18.166 /usr/local/
}

# http://www.aribeiro.net.br/mktop/
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

function InstallParmEd() {
  pip3.9 install parmed || exit 1
}

function InstallASE() {
  pip3.9 install ase || exit 1
}

# http://towhee.sourceforge.net/
function InstallTowhee() {
  rm -rf /usr/local/towhee-8.2.3
  wget https://winmostar.com/wm/cygwin_wm/packages/towhee-8.2.3.tar.gz
  rm -fr towhee-8.2.3
  tar zxvf towhee-8.2.3.tar.gz
  cd towhee-8.2.3
  
  mkdir -p /usr/local/towhee-8.2.3
  ./configure --enable-safe-compare --enable-fix-GNU --prefix=/usr/local/towhee-8.2.3 || exit 1
  sed -i -e "s/define NTMAX 7/define NTMAX 10/" \
         -e "s/define NUMAX 2202/define NUMAX 1000/" \
         -e "s/define TTORMAX 930/define TTORMAX 3000/" \
         -e "s/define TIMPMAX 120/define TIMPMAX 200/" \
         -e "s/define NNTYPE 375/define NNTYPE 1000/" \
         Source/preproc.h || exit 1
  cd Source
  make -j4 || exit 1
  make install || exit 1
  cd ..
  
  cp -r ForceFields /usr/local/towhee-8.2.3/
  
  cd ..
}

# http://theory.cm.utexas.edu/henkelman/code/bader/
function InstallBader() {
  rm -rf /usr/local/bader-1.04
  wget https://winmostar.com/wm/cygwin_wm/packages/bader_v1.04.tar.gz
  rm -fr bader
  tar zxvf bader_v1.04.tar.gz
  cd bader
  
  make -f makefile.osx_gfortran || exit 1
  mkdir -p /usr/local/bader-1.04/bin
  cp bader.exe /usr/local/bader-1.04/bin/
  
  cd ..
}

set -x

date

SCRIPT_DIR=$(cd $(dirname $0); pwd)

export LANG=C

grep "export LANG=C" /etc/profile.d/lang.sh || echo "export LANG=C" >> /etc/profile.d/lang.sh

cat << EOF > /etc/profile.d/winmostar.sh
# To remove PATH from Windows
export PATH=\`echo \$PATH | awk -v RS=: '{print \$0}' | \\
  grep -v "Microsoft MPI" | \\
  grep -v "MPICH2" | \\
  grep -v "Quantum ESPRESSO" | \\
  grep -v "LAMMPS" | \\
  awk -v ORS=: '{print \$0}'\`
export PATH=\$PATH:/usr/local/mpich2-1.5/bin
source /usr/local/gromacs_sse/bin/GMXRC
#source /usr/local/gromacs_avx/bin/GMXRC
#source /usr/local/amber18/amber.sh
source /usr/local/amber18/amber.sh
export PATH=\$PATH:/usr/local/acpype
export PATH=\$PATH:/usr/local/MODYLAS_1.0.4/bin
#export PATH=\$PATH:/usr/local/openmx3.8/work:/usr/local/fermisurfer/bin
#export OPENMX_DATA_PATH=/usr/local/openmx3.8/DFT_DATA13
export PATH=\$PATH:/usr/local/openmx3.9/work:/usr/local/fermisurfer/bin
export OPENMX_DATA_PATH=/usr/local/openmx3.9/DFT_DATA19
export PATH=\$PATH:/usr/local/boltztrap-1.2.5/src
export PATH=\$PATH:/usr/local/boltztrap-1.2.5/util
export PerlChemistry=/usr/local/MATCH_RELEASE/PerlChemistry
export MATCH=/usr/local/MATCH_RELEASE/MATCH
export PATH=\$PATH:/usr/local/MATCH_RELEASE/MATCH/scripts
export PATH=\$PATH:/usr/local/packmol-18.166
export PATH=\$PATH:/usr/local/mktop_2.2.1
export PATH=\$PATH:/usr/local/towhee-8.2.3/bin
export TOWHEE_FF_PATH=/usr/local/towhee-8.2.3/ForceFields
export PATH=\$PATH:/usr/local/bader-1.04/bin
EOF

. /etc/profile.d/winmostar.sh

cd $SCRIPT_DIR
cp -f ChangeLog setup-x86_64.exe /cygdrive/c/cygwin_wm/
cp -rf src/* /tmp/

cd /tmp/

if [ ! -f /bin/python3.exe ]; then
  ln -s /bin/python3.9.exe /bin/python3.exe
fi

if [ ! -f /bin/python2.exe ]; then
  ln -s /bin/python2.7.exe /bin/python2.exe
fi

# Reqiured to run acpype
echo 1 | /usr/sbin/update-alternatives --config python


# required to build phonopy
# https://seesaawiki.jp/w/kou1okada/d/20200205%3A%20Cygwin%20-%20Python3%20-%20h5py
if [ ! -f /bin/hdf5.dll ]; then
  ln -s /bin/cyghdf5-103.dll /bin/hdf5.dll
fi

# os.symlink of python2.7 makes symbolic links to /mnt/ for windows paths
ln -s /cygdrive /mnt

InstallMPICH2
InstallOpenBabel
InstallGromacs
InstallAmberTools18
InstallAcpype18
InstallERmod
InstallMODYLAS
#InstallOpenMX
InstallOpenMX3.9
InstallBoltzTraP
InstallMATCH
InstallConditionalERmod
InstallEnumlib
InstallPackmol
InstallMktop
InstallPhonopy
InstallMDTraj
InstallParmEd
InstallASE
InstallTowhee
InstallBader

date
