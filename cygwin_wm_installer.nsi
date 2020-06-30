; NSIS Script

!include MUI2.nsh

!verbose 4

; /SOLIDをつけるとエラーになる(サイズが大きすぎるから？)
SetCompressor lzma

; インストーラーの名前
Name "Cygwin for Winmostar"

; インストーラーのファイル名
OutFile "cygwin_wm.exe"

; インストール先ディレクトリ
InstallDir "C:\cygwin_wm"

; 必要な実行権限（user/admin）
RequestExecutionLevel user

; インストール時に表示するものを指定
; 上書き禁止
!define MUI_PAGE_CUSTOMFUNCTION_LEAVE ValidateDirectory
!insertmacro MUI_PAGE_DIRECTORY
!define MUI_PAGE_CUSTOMFUNCTION_LEAVE CreateSymlinks
!insertmacro MUI_PAGE_INSTFILES
!insertmacro MUI_PAGE_FINISH

!insertmacro MUI_LANGUAGE "English"

; インストール時の挙動
Section "Install" ; この名前は重要ではない

  ; キャンセルできるように
  GetDlgItem $0 $HWNDPARENT 2
  EnableWindow $0 1

  SectionIn RO

  ; インストール先ディレクトリを指定
  SetOutPath $INSTDIR

  SetOverwrite off
  File /a /r "C:\cygwin_wm\*.*"

SectionEnd

Function ValidateDirectory
  StrCpy $0 $INSTDIR 16
  StrCmp $0 "C:\Program Files" folderwarning
  StrCpy $0 $INSTDIR 22
  StrCmp $0 "C:\Program Files (x86)" folderwarning
  Goto folderok

  folderwarning:

  MessageBox MB_OK|MB_ICONSTOP "Cygwin does not work properly with folders that require administrator privileges for writing."
  Abort

  folderok:

  ${If} ${FileExists} `$INSTDIR\*.*`
    MessageBox MB_OK|MB_ICONSTOP "Installation folder already exists. Move or delete the folder before installation"
    Abort
  ${EndIf}

FunctionEnd

Function CreateSymlinks
  ExecWait '"$INSTDIR\make_symlink.cmd"'
FunctionEnd

