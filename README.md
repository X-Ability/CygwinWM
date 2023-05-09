# CygwinWM
Cygwin Setup for Winmostar

## Description

Download and install Cygwin and various packages.

## Prerequisite

MODYLAS is a freeware, but is not allowed to be redistributed.   
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

If you get an error message like this.


>      1 [main] python3.7m 4356 child_info_fork::abort: unable to remap \??\C:\cygwin_wm\lib\python3.7\site-packages\Cython-0.29.13-py3.7-cygwin-3.1.6-i686.egg\Cython\Compiler\FlowControl.cpython-37m-i386-cygwin.dll (using C:\cygwin_wm\lib\python3.7\site-packages\Cython-0.29.13-py3.7-cygwin-3.1.6-i686.egg\Cython\Compiler\FlowControl.cpython-37m-i386-cygwin.dll) to same address as parent (0x1B80000) - try running rebaseall


1. Double click C:\cygwin_wm\bin\ash.exe
2. Execute rebase command like this.

```bash
$  /bin/rebase -s -v /lib/python3.7/site-packages/Cython-0.29.13-py3.7-cygwin-3.1.6-i686.egg/Cython/Compiler/*.dll
```
3. Restart cygwin_setup.vbs as described above.

## License

| Package          | License    | Remarks                          | Version   | Reference                                                                        |
| ---------------- | ---------- | -------------------------------- | ----------|--------------------------------------------------------------------------------- |
| Cygwin           | GPL        |                                  |           | https://cygwin.com/licensing.html                                                |
| MPICH            | Original   | BSD like license                 | 2-1.5     | https://github.com/pmodels/mpich/blob/master/COPYRIGHT                           |
| Grace            | GPL        |                                  | 5.1.25    | https://plasma-gate.weizmann.ac.il/Grace/doc/GPL.html                            |
| OpenBabel        | GPL        |                                  | 2.4.1     | https://openbabel.org/wiki/Frequently_Asked_Questions                            |
| Gromacs          | LGPL       |                                  | 5.0.7     | http://www.gromacs.org/About_Gromacs                                             |
| AmberTools       | GPL etc.   |                                  | 18        | https://ambermd.org/AmberTools.php                                               |
| Acpype           | GPL        |                                  | r10101    | https://github.com/alanwilter/acpype/blob/master/LICENSE                         |
| ERmod            | GPL        |                                  | 0.3.4     | https://sourceforge.net/projects/ermod/                                          |
| ERmod(vmdplugin) | UIUC       | based on MIT/X11 and BSD license |           | https://sourceforge.net/p/ermod/code/ci/default/tree/vmdplugins/LICENSE          |
| MODYLAS          | Original   |                                  | 1.0.4     | https://www.modylas.org/node/18                                                  |
| OpenMX           | GPL        |                                  | 3.8.4     | http://www.openmx-square.org/whatisopenmx.html                                   |
| FermiSurfer      | MIT        |                                  | 1.7.1     | https://fermisurfer.osdn.jp/en/_build/html/copy.html                             |
| BoltzTraP        | GPL        |                                  | 1.2.5     | https://www.imc.tuwien.ac.at/index.php?id=21094                                  |
| MATCH            | MIT        |                                  |           | https://brooks.chem.lsa.umich.edu/download/software/match/MATCH_Users_Manual.pdf |
| enumlib          | MIT        |                                  | 1.0.8     | https://github.com/msg-byu/enumlib/blob/master/LICENSE                           |
| symlib           | MIT        |                                  | 1.1.0     | https://github.com/msg-byu/symlib/blob/master/LICENSE                            |
| packmol          | MIT        |                                  | 18.166    | https://github.com/m3g/packmol/blob/master/LICENSE                               |
| MKTOP            | GPL        |                                  | 2.2.1     | https://github.com/aar2163/MKTOP/blob/master/mktop.pl                            |
| Pymatgen         | MIT        |                                  | 2020.4.29 | https://pymatgen.org/                                                            |
| Cython           | Apache     |                                  | 0.29.13   | https://cython.org/                                                              |
| ScyPy            | BSD        |                                  | 1.1.0     | https://www.scipy.org/scipylib/license.html                                      |
| Spglib           | BSD        |                                  | 1.15.1    | https://spglib.github.io/spglib/                                                 |
| matplotlib       | Original   | BSD compatible                   | 3.1.0     | https://matplotlib.org/3.1.0/users/license.html                                  |
| pandas           | BSD        |                                  | 1.0.3     | https://github.com/pandas-dev/pandas/blob/master/LICENSE                         |
| Phonopy          | BSD        |                                  | 1.12.6.53 | https://phonopy.github.io/phonopy/                                               |
| MDTraj           | LGPL       |                                  | 1.9.0     | https://github.com/mdtraj/mdtraj                                                 |
| ParmEd           | LGPL       |                                  | 3.2.0     | https://github.com/ParmEd/ParmEd                                                 |
| Towhee           | GPL        |                                  | 3.2.0     | http://towhee.sourceforge.net/                                                   |
| Bader            | GPL        |                                  | v1.04     | http://theory.cm.utexas.edu/henkelman/code/bader/                                |

UIUC : University of Illinois Open Source License

The files in src folder, which are created by X-Ability, subject to the same license as the original packages.  

