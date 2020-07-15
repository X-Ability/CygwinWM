# CygwinWM
Cygwin Setup for Winmostar

## Description

Download and install Cygwin and various packages.

## Prerequisite

MODYLAS is freeware, but is not allowed to be redistributed.   
X-Ability Co. Ltd. has special permission to distribute it as a binary package.  
If you want to install MODYLAS, download MODYLAS_1.0.4.tar_1.gz manually and put it in src folder.  
If you don't need to install MODYLAS, comment out the line that calls InstallMODYLAS in cygwin_setup.sh.  

## Usage

1. If the C:\cygwin_wm folder exists, delete or rename it.
2. Double click cygwin_setup.vbs.

If the installation script exits abnormally and you want to resume from the middle of the script.
1. Set CLEAN_INSTALL in cygwin_setup.vbs to False.
2. Comment out the line that calls an unnecessary installation function in cygwin_setup.sh.
3. Double click cygwin_setup.vbs again.

## License

| Package          | License    | Remarks                          | Reference                                                                        |
| ---------------- | ---------- | -------------------------------- | -------------------------------------------------------------------------------- |
| Cygwin           | GPL        |                                  | https://cygwin.com/licensing.html                                                |
| MPICH            | Original   | BSD like license                 | https://github.com/pmodels/mpich/blob/master/COPYRIGHT                           |
| Grace            | GPL        |                                  | https://plasma-gate.weizmann.ac.il/Grace/doc/GPL.html                            |
| OpenBabel        | GPL        |                                  | https://openbabel.org/wiki/Frequently_Asked_Questions                            |
| Gromacs          | LGPL       |                                  | http://www.gromacs.org/About_Gromacs                                             |
| AmberTools       | GPL etc.   |                                  | https://ambermd.org/AmberTools.php                                               |
| Acpype           | GPL        |                                  | https://github.com/alanwilter/acpype/blob/master/LICENSE                         |
| ERmod            | GPL        |                                  | https://sourceforge.net/projects/ermod/                                          |
| ERmod(vmdplugin) | UIUC       | based on MIT/X11 and BSD license | https://sourceforge.net/p/ermod/code/ci/default/tree/vmdplugins/LICENSE          |
| MODYLAS          | Original   |                                  | https://www.modylas.org/node/18                                                  |
| OpenMX           | GPL        |                                  | http://www.openmx-square.org/whatisopenmx.html                                   |
| FermiSurfer      | MIT        |                                  | http://www.openmx-square.org/whatisopenmx.html                                   |
| BoltzTraP        | GPL        |                                  | https://www.imc.tuwien.ac.at/index.php?id=21094                                  |
| MATCH            | MIT        |                                  | https://brooks.chem.lsa.umich.edu/download/software/match/MATCH_Users_Manual.pdf |
| enumlib          | MIT        |                                  | https://github.com/msg-byu/enumlib/blob/master/LICENSE                           |
| symlib           | MIT        |                                  | https://github.com/msg-byu/symlib/blob/master/LICENSE                            |
| packmol          | MIT        |                                  | https://github.com/m3g/packmol/blob/master/LICENSE                               |
| MKTOP            | GPL        |                                  | https://github.com/aar2163/MKTOP/blob/master/mktop.pl                            |
| Pymatgen         | MIT        |                                  | https://pymatgen.org/                                                            |
| Cython           | Apache     |                                  | https://cython.org/                                                              |
| ScyPy            | BSD        |                                  | https://www.scipy.org/scipylib/license.html                                      |
| Spglib           | BSD        |                                  | https://spglib.github.io/spglib/                                                 |
| matplotlib       | Original   | BSD compatible                   | https://matplotlib.org/3.1.0/users/license.html                                  |
| pandas           | BSD        |                                  | https://github.com/pandas-dev/pandas/blob/master/LICENSE                         |
| Phonopy          | BSD        |                                  | https://phonopy.github.io/phonopy/                                               |
| MDTraj           | LGPL       |                                  | https://fermisurfer.osdn.jp/en/_build/html/copy.html                             |

UIUC : University of Illinois Open Source License

XA-created files in the src folder subject to the same license as the original packages.  

