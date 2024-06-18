
CLEAN_INSTALL = False
INSTALL_QE_LMP = False

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
cmd = ".\setup-x86_64.exe " &_
  "--site http://cygwin.mirror.constant.com " &_
  "--no-shortcuts  " &_
  "--no-desktop " &_
  "--no-admin " &_
  "--quiet-mode " &_
  "--root C:\cygwin_wm " &_
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
  "--packages python39 " &_
  "--packages python39-devel " &_
  "--packages python39-pip " &_
  "--packages libopenblas " &_
  "--packages liblapack-devel " &_
  "--packages wget " &_
  "--packages zip " &_
  "--packages libiconv " &_
  "--packages connect-proxy "



If CLEAN_INSTALL Then
    objShell.Run cmd, 1, True
End If

strShDir  = Replace(strScriptDir, ":\", "/")
strShDir  = Replace(strShDir, "\", "/")
strShFile  = "/cygdrive/" & strShDir & "/cygwin_setup.sh"
strLogFile = "/cygdrive/" & strShDir & "/cygwin_setup.log"

cmd = "C:\cygwin_wm\bin\bash.exe --login -c """ & strShFile & " 2>&1 | tee " & strLogFile & """"
objShell.Run cmd, 1, True

If INSTALL_QE_LMP Then
    strShFile  = "/" & strShDir & "/mingw32_setup.sh"
    strLogFile = "/" & strShDir & "/mingw32_setup.log"
    
    cmd = "C:\msys64\msys2_shell.cmd -mingw32 -c """ & strShFile & " 2>&1 | tee " & strLogFile & """"
    objShell.Run cmd, 1, True
    
    strShFile  = "/" & strShDir & "/mingw64_setup.sh"
    strLogFile = "/" & strShDir & "/mingw64_setup.log"
    
    cmd = "C:\msys64\msys2_shell.cmd -mingw64 -c """ & strShFile & " 2>&1 | tee " & strLogFile & """"
    objShell.Run cmd, 1, True
End If

