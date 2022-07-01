#!/bin/bash

## Usage 
##
## $ sudo su -
## # apt -y update
## # apt -y upgrade
## # bash WSL_setup.sh >& WSL_setup.log
##


function InstallGromacs() {
  rm -rf /usr/local/gromacs gromacs-5.0.7.tar.gz gromacs-5.0.7
  wget https://winmostar.com/wm/cygwin_wm/packages/gromacs-5.0.7.tar.gz
  tar xfz gromacs-5.0.7.tar.gz
  cd gromacs-5.0.7
  mkdir build
  cd build
  cmake .. -DGMX_BUILD_OWN_FFTW=ON -DGMX_GPU=OFF -DGMX_SIMD=AVX_256 -DGMX_DOUBLE=OFF -DGMX_USE_RDTSCP=OFF -DCMAKE_INSTALL_PREFIX=/usr/local/gromacs
  make
  make install
  cd ..
  mkdir build_d
  cd build_d
  cmake .. -DGMX_BUILD_OWN_FFTW=ON -DGMX_GPU=OFF -DGMX_SIMD=AVX_256 -DGMX_DOUBLE=ON -DGMX_USE_RDTSCP=OFF -DCMAKE_INSTALL_PREFIX=/usr/local/gromacs
  make
  make install
  cd ..
}

function InstallAmberTools18() {
  rm -rf /usr/local/amber18/ AmberTools18.tar.bz2
  dir=`pwd`
  wget https://winmostar.com/wm/cygwin_wm/packages/AmberTools18.tar.bz2
  tar xfj AmberTools18.tar.bz2
  mv amber18 /usr/local/
  export AMBERHOME=/usr/local/amber18
  cd $AMBERHOME
  echo y | ./configure -noX11 --skip-python gnu
  source amber.sh
  make install
  cd $dir
}

function InstallAcpype18() {
  rm -rf /usr/local/acpype/ ./acpype
  svn checkout http://svn.code.sf.net/p/ccpn/code/branches/stable/ccpn/python/acpype/ acpype -r 10101
  wget https://winmostar.com/wm/cygwin_wm/packages/acpype.py
  mv acpype.py /usr/loca/acpype/
  mv acpype /usr/local/
}

function InstallERmod() {
  rm -rf /usr/local/ermod ermod-0.3.4.tar.gz ermod-0.3.4
  wget https://winmostar.com/wm/cygwin_wm/packages/ermod-0.3.4.tar.gz
  tar xfz ermod-0.3.4.tar.gz
  cd ermod-0.3.4/vmdplugins
  make
  cp compile/*.so libexec/
  cd ..
  ./configure --prefix=/usr/local/ermod --disable-mpi --enable-openmp
  make
  make install
  cd ..
}

function InstallMODYLAS() {
  rm -rf /usr/local/MODYLAS_1.0.4 MODYLAS_1.0.4.tar_1.gz MODYLAS_1.0.4.patch MODYLAS_1.0.4
  wget https://winmostar.com/wm/cygwin_wm/packages/MODYLAS_1.0.4.tar_1.gz
  wget https://winmostar.com/wm/cygwin_wm/packages/MODYLAS_1.0.4.patch
  tar xfz MODYLAS_1.0.4.tar_1.gz
  patch -u -p0 < MODYLAS_1.0.4.patch
  cd MODYLAS_1.0.4/source
  ./configure --with-kind-fortran-compiler=INTEL --prefix=/usr/local/MODYLAS_1.0.4
  cd src
  make
  make install
  cd ../../..
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
  cp mktop_2.2.1.pl mktop.pl
  chmod 755 mktop.pl
  cd ..
  mv mktop_2.2.1 /usr/local/
}

function InstallNWChem() {
  dir=`pwd`
  cd /tmp/
  rm -rf /usr/local/NWChem  Nwchem-6.6.revision27746-src.2015-10-20.tar.gz nwchem-6.6
  #wget http://www.nwchem-sw.org/download.php?f=Nwchem-6.6.revision27746-src.2015-10-20.tar.gz -O Nwchem-6.6.revision27746-src.2015-10-20.tar.gz
  wget https://winmostar.com/wm/cygwin_wm/packages/NWChem/Nwchem-6.6.revision27746-src.2015-10-20.tar.gz
  tar xfz Nwchem-6.6.revision27746-src.2015-10-20.tar.gz
  cd nwchem-6.6

  #for f in Tddft_mxvec20 Tools_lib64 Config_libs66 Cosmo_meminit Sym_abelian Xccvs98 Dplot_tolrho Driver_smalleig Ga_argv Raman_displ Ga_defs Zgesvd Cosmo_dftprint Txs_gcc6 Gcc6_optfix Util_gnumakefile Util_getppn Gcc6_macs_optfix Notdir_fc Xatom_vdw Hfmke; do wget http://www.nwchem-sw.org/images/${f}.patch.gz; done
  for f in Tddft_mxvec20 Tools_lib64 Config_libs66 Cosmo_meminit Sym_abelian Xccvs98 Dplot_tolrho Driver_smalleig Ga_argv Raman_displ Ga_defs Zgesvd Cosmo_dftprint Txs_gcc6 Gcc6_optfix Util_gnumakefile Util_getppn Gcc6_macs_optfix Notdir_fc Xatom_vdw Hfmke; do wget https://winmostar.com/wm/cygwin_wm/packages/NWChem/${f}.patch.gz; done
  for f in *.patch.gz; do gunzip $f; done
  for f in *.patch; do patch -p0 < $f; done

  export NWCHEM_TOP=/tmp/nwchem-6.6
  export USE_MPI=y
  export NWCHEM_TARGET=LINUX64
  export USE_PYTHONCONFIG=y
  export PYTHONVERSION=2.7
  export PYTHONHOME=/usr
  export BLASOPT="-lopenblas -lpthread -lrt"
  export BLAS_SIZE=4
  export USE_64TO32=y

  cd src

  make nwchem_config NWCHEM_MODULES="all python"
  make 64_to_32
  make

  mkdir /usr/local/NWChem
  mkdir /usr/local/NWChem/bin
  mkdir /usr/local/NWChem/data

  cp $NWCHEM_TOP/bin/${NWCHEM_TARGET}/nwchem /usr/local/NWChem/bin
  chmod 755 /usr/local/NWChem/bin/nwchem

  cp -r $NWCHEM_TOP/src/basis/libraries /usr/local/NWChem/data
  cp -r $NWCHEM_TOP/src/data /usr/local/NWChem
  cp -r $NWCHEM_TOP/src/nwpw/libraryps /usr/local/NWChem/data

  cat << EOF > /usr/local/NWChem/data/default.nwchemrc
nwchem_basis_library /usr/local/NWChem/data/libraries/
nwchem_nwpw_library /usr/local/NWChem/data/libraryps/
ffield amber
amber_1 /usr/local/NWChem/data/amber_s/
amber_2 /usr/local/NWChem/data/amber_q/
amber_3 /usr/local/NWChem/data/amber_x/
amber_4 /usr/local/NWChem/data/amber_u/
spce /usr/local/NWChem/data/solvents/spce.rst
charmm_s /usr/local/NWChem/data/charmm_s/
charmm_x /usr/local/NWChem/data/charmm_x/
EOF

  ln -s /usr/local/NWChem/data/default.nwchemrc $HOME/.nwchemrc

  cd $dir
}

function InstallLAMMPS() {
  rm -rf /usr/local/lammps-30Jul16 lammps-30Jul16.tar.gz lammps-30Jul16
  wget https://winmostar.com/wm/cygwin_wm/packages/lammps-30Jul16.tar.gz
  tar xfz lammps-30Jul16.tar.gz
  cd lammps-30Jul16/src
  make yes-misc
  make yes-rigid
  make yes-user-reaxc
  make serial
  make mpi
  cd ../..
  mv lammps-30Jul16 /usr/local/
}

function InstallBoltzTraP() {
  rm -rf /usr/local/boltztrap-1.2.5 BoltzTraP.tar.bz2 boltztrap-1.2.5
  wget https://winmostar.com/wm/cygwin_wm/packages/BoltzTraP.tar.bz2
  tar xfj BoltzTraP.tar.bz2
  cd boltztrap-1.2.5/src
  sed -i s/:log/log/g x_trans
  make
  cd ../..
  mv boltztrap-1.2.5 /usr/local
  sed -i.bak "s/<policy domain=\"coder\" rights=\"none\" pattern=\"PS\" \/>/<policy domain=\"coder\" rights=\"read|write\" pattern=\"PS\" \/>/g" /etc/ImageMagick-6/policy.xml
}

function InstallPymatgen() {
  wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
  bash miniconda.sh -b -p $HOME/miniconda3
  source $HOME/.bash_profile
  conda install -y python=3.7
  conda install -y --channel conda-forge pymatgen
}

set -x

export LANG=C

cd /tmp/

apt -y update
apt -y upgrade
apt -y install cmake g++ csh gfortran python2.7 flex zlib1g-dev \
  libbz2-dev subversion python libnetcdf-dev liblapack-dev python-numpy \
  openbabel grace imagemagick gnuplot bc dos2unix libopenmpi-dev python-dev \
  libopenblas-dev openmpi-bin tcsh libfftw3-dev autoconf libboost-dev libxml2-dev
  python3-dev python3-setuptools python3-pip

cat << EOF > /etc/profile.d/winmostar.sh
source /usr/local/gromacs/bin/GMXRC
source /usr/local/amber18/amber.sh
export PATH=\$PATH:/usr/local/acpype
export PerlChemistry=/usr/local/MATCH_RELEASE/PerlChemistry
export MATCH=/usr/local/MATCH_RELEASE/MATCH
export PATH=\$PATH:/usr/local/MATCH_RELEASE/MATCH/scripts
export PATH=\$PATH:/usr/local/MODYLAS_1.0.4/bin
export PATH=\$PATH:/usr/local/NWChem/bin
export NWCHEM_BASIS_LIBRARY=/usr/local/NWChem/data/libraries/
export PATH=\$PATH:/usr/local/espresso-5.2.1/bin
export ESPRESSO_PSEUDO=/usr/local/espresso-5.2.1/pseudo
export PATH=\$PATH:/usr/local/lammps-30Jul16/src
export LAMMPS_POTENTIALS=/usr/local/lammps-30Jul16/potentials
export PATH=\$PATH:/usr/local/boltztrap-1.2.5/util
export PATH=\$PATH:/usr/local/boltztrap-1.2.5/src
export PATH=\$PATH:/usr/local/packmol-18.166
export PATH=\$PATH:/usr/local/mktop_2.2.1
EOF

cat << EOF >> $HOME/.bash_profile
# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/home/furukawa/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/home/furukawa/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/home/furukawa/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/home/furukawa/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

EOF

InstallGromacs
InstallAmberTools18
InstallAcpype18
InstallERmod
InstallMODYLAS
InstallNWChem
InstallLAMMPS
InstallBoltzTraP
InstallPackmol
InstallMktop

# for Pymatgen
InstallPymatgen

#rm -rf /tmp/*
