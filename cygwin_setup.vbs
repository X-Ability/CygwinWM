
CLEAN_INSTALL = True

Set objShell    = WScript.CreateObject("WScript.Shell")
Set objFSO      = CreateObject("Scripting.FileSystemObject")
Set objShellApp = CreateObject("Shell.Application")

strScriptFile = WScript.ScriptFullName
strScriptDir  = objFSO.GetParentFolderName(strScriptFile)


If CLEAN_INSTALL And objFso.FolderExists("C:\cygwin_wm") Then
    WScript.Echo "C:\cygwin_wm already exists."
    WScript.Quit
End If

' Install Cygwin
' https://www.cygwin.com/
cmd = ".\setup-x86.exe " &_
  "--site http://cygwin.mirror.constant.com " &_
  "--no-shortcuts  " &_
  "--no-desktop " &_
  "--no-admin " &_
  "--quiet-mode " &_
  "--root C:\cygwin_wm " &_
  "--arch x86 " &_
  "--local-package-dir C:\cygwin_wm\cygwin-packages " &_
  "--verbose " &_
  "--prune-install " &_
  "--packages libbz2-devel " &_
  "--packages unzip " &_
  "--packages autoconf " &_
  "--packages cmake " &_
  "--packages flex " &_
  "--packages gcc-fortran " &_
  "--packages gcc-g++ " &_
  "--packages make " &_
  "--packages pkg-config " &_
  "--packages subversion " &_
  "--packages gnuplot " &_
  "--packages ImageMagick " &_
  "--packages ghostscript " &_
  "--packages m4 " &_
  "--packages libboost-devel " &_
  "--packages libfftw3-devel " &_
  "--packages libfreetype-devel " &_
  "--packages libmpfr4 " &_
  "--packages libnetcdf-devel " &_
  "--packages libpng-devel " &_
  "--packages zlib-devel " &_
  "--packages bc " &_
  "--packages liblapack-devel " &_
  "--packages openssh " &_
  "--packages python-h5py " &_
  "--packages python27-cython " &_
  "--packages python27-devel " &_
  "--packages python27-numpy " &_
  "--packages python27-pip " &_
  "--packages python27-yaml " &_
  "--packages tcsh " &_
  "--packages dos2unix " &_
  "--packages patchutils " &_
  "--packages libhdf5-devel " &_
  "--packages python37 " &_
  "--packages python37-devel " &_
  "--packages python37-numpy " &_
  "--packages python37-pip " &_
  "--packages python37-yaml " &_
  "--packages wget "

If CLEAN_INSTALL Then
    objShell.Run cmd, 1, True
End If

strShDir  = Replace(strScriptDir, ":\", "/")
strShDir  = Replace(strShDir, "\", "/")
strShFile  = "/cygdrive/" & strShDir & "/cygwin_setup.sh"
strLogFile = "/cygdrive/" & strShDir & "/cygwin_setup.log"

cmd = "C:\cygwin_wm\bin\bash.exe --login -c """ & strShFile & " 2>&1 | tee " & strLogFile & """"
objShell.Run cmd, 1, True

