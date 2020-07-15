# CygwinWM
Cygwin Setup for Winmostar

## Description

Download and install Cygwin and various packages.

## Prerequisite

MODYLAS is not allowed to be redistributed. 
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






