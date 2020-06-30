set workdir=%~dp0
cd %workdir%

..\NSIS\makensis.exe cygwin_wm_installer.nsi

..\CodeSigning\sign.bat cygwin_wm.exe


